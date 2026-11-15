# Paired-read merging (port of dada2::mergePairs) and the public aligner.

"""
    rc(seq) -> String

Reverse complement of a DNA string. IUPAC codes and case are kept.
"""
rc(s::AbstractString) = String(revcomp(codeunits(s)))

const _COMPLEMENT = let t = collect(UInt8, 0:255)
    for (a, b) in zip("ACGTMRWSYKVHDBN-acgtmrwsykvhdbn", "TGCAKYWSRMBDHVN-tgcakywsrmbdhvn")
        t[UInt8(a)+1] = UInt8(b)
    end
    t
end

"Reverse complement of a DNA byte sequence (IUPAC codes included; other bytes kept)."
revcomp(v::AbstractVector{UInt8}) = [_COMPLEMENT[v[i]+1] for i in reverse(eachindex(v))]

"""
    nwalign(s1, s2; match, mismatch, gap, homo_gap=nothing, band=-1,
            endsfree=true, vec=false) -> (String, String)

Needleman-Wunsch alignment of two sequences, as R's `dada2::nwalign`. The scores
default to the session options `MATCH`, `MISMATCH` and `GAP_PENALTY` (see
[`setDadaOpt`](@ref)).
"""
function nwalign(s1::AbstractString, s2::AbstractString; match::Integer=getDadaOpt(:MATCH),
                 mismatch::Integer=getDadaOpt(:MISMATCH), gap::Integer=getDadaOpt(:GAP_PENALTY),
                 homo_gap::Union{Nothing,Integer}=nothing, band::Integer=-1,
                 endsfree::Bool=true, vec::Bool=false)
    homo_gap === nothing && (homo_gap = gap)
    buf = AlignBuffer()
    a = Vector{UInt8}(codeunits(s1)); b = Vector{UInt8}(codeunits(s2))
    if vec
        homo_gap != gap && error("Homopolymer gap penalties are not implemented in the vectorized aligner.")
        nwalign_vectorized!(buf, a, b, match, mismatch, gap, endsfree ? 0 : gap, band)
    else
        (isACGT(s1) && isACGT(s2)) || error("Sequences must contain only A/C/G/T characters.")
        if endsfree && homo_gap != gap
            nwalign_endsfree_homo!(buf, encode_nt(s1), encode_nt(s2), match, mismatch, gap, homo_gap, band)
            return (_decode_al(buf.al0, buf.len), _decode_al(buf.al1, buf.len))
        end
        nwalign_endsfree!(buf, a, b, match, mismatch, gap, band; endsfree=endsfree)
    end
    return (String(buf.al0[1:buf.len]), String(buf.al1[1:buf.len]))
end
_decode_al(v, n) = String([c == GAP ? GAP : NT_CHAR[c] for c in view(v, 1:n)])

"""
    mergePairs(dadaF, derepF, dadaR, derepR; minOverlap=12, maxMismatch=0,
               returnRejects=false, propagateCol=String[], justConcatenate=false,
               trimOverhang=false, verbose=false, progress)

Merge denoised forward and reverse reads (port of R's `mergePairs`). Accepts
single results or vectors (one per sample). Returns a `DataFrame` (or a vector
of them) with columns `sequence, abundance, forward, reverse, nmatch,
nmismatch, nindel, prefer, accept`, sorted by decreasing abundance.

`progress=false` hides the progress bar (default: [`DenoiseDNA.setprogress!`](@ref)).
"""
function mergePairs(dadaF, derepF, dadaR, derepR; progress::Bool=_progress_default(), kwargs...)
    if dadaF isa DadaResult
        return _merge_one(dadaF, derepF, dadaR, derepR; kwargs...)
    end
    n = length(dadaF)
    (length(derepF) == n && length(dadaR) == n && length(derepR) == n) ||
        error("The dadaF/derepF/dadaR/derepR arguments must be the same length.")
    bar = _progress(n, "Merging"; unit=" samples", enabled=progress)
out = map(1:n) do i
    m = _merge_one(dadaF[i], derepF[i], dadaR[i], derepR[i]; kwargs...)
    _next!(bar)
    m
end
nin = sum(i -> sum(dadaF[i].abundance), 1:n; init=0)
nok = sum(m -> sum(m.abundance[m.accept]), out; init=0)
_finish!(bar; info="$(_fmt(nok)) reads merged ($(_pct(nok / max(nin, 1))) of denoised forward reads)")
return out
end

function _merge_one(dF::DadaResult, drF, dR::DadaResult, drR; minOverlap::Integer=12,
                    maxMismatch::Integer=0, returnRejects::Bool=false,
                    propagateCol::AbstractVector=String[], justConcatenate::Bool=false,
                    trimOverhang::Bool=false, verbose::Bool=false)
    drF isa AbstractString && (drF = derepFastq(drF))
    drR isa AbstractString && (drR = derepFastq(drR))
    mapF = drF.map; mapR = drR.map
    (length(mapF) == length(mapR) && maximum(mapF) == length(dF.map) && maximum(mapR) == length(dR.map)) ||
        error("Non-corresponding derep and dada objects.")
    nread = length(mapF)
    # unique (forward, reverse) pairs in order of first appearance; count all
    pairidx = Dict{Tuple{Int,Int},Int}()
    fw = Int[]; rv = Int[]; ab = Int[]
    @inbounds for k in 1:nread
        f = mapF[k] == 0 ? 0 : dF.map[mapF[k]]
        r = mapR[k] == 0 ? 0 : dR.map[mapR[k]]
        (f == 0 || r == 0) && continue
        key = (f, r)
        j = get(pairidx, key, 0)
        if j == 0
            push!(fw, f); push!(rv, r); push!(ab, 1)
            pairidx[key] = length(fw)
        else
            ab[j] += 1
        end
    end
    np = length(fw)
    seqs = Vector{String}(undef, np)
    nmatch = zeros(Int, np); nmis = zeros(Int, np); nind = zeros(Int, np)
    prefer = Vector{Union{Missing,Int}}(missing, np); accept = trues(np)
    if justConcatenate
        for k in 1:np
            seqs[k] = dF.clustering.sequence[fw[k]] * "NNNNNNNNNN" * rc(dR.clustering.sequence[rv[k]])
        end
    else
        m, mm, g = maxMismatch == 0 ? (1, -64, -64) : (1, -8, -8)
        buf = AlignBuffer()
        n0F = dF.clustering.n0; n0R = dR.clustering.n0
        for k in 1:np
            a = Vector{UInt8}(codeunits(dF.clustering.sequence[fw[k]]))
            b = Vector{UInt8}(codeunits(rc(dR.clustering.sequence[rv[k]])))
            nwalign_endsfree!(buf, a, b, m, mm, g, -1)
            al0 = buf.al0[1:buf.len]; al1 = buf.al1[1:buf.len]
            ev = eval_pair(al0, al1)
            nmatch[k] = ev.match; nmis[k] = ev.mismatch; nind[k] = ev.indel
            prefer[k] = 1 + (n0R[rv[k]] > n0F[fw[k]])
            accept[k] = ev.match >= minOverlap && (ev.mismatch + ev.indel) <= maxMismatch
            seqs[k] = String(pair_consensus(al0, al1, prefer[k], trimOverhang))
        end
    end
    for k in 1:np
        accept[k] || (seqs[k] = "")
    end
    df = DataFrame(sequence=seqs, abundance=ab, forward=fw, reverse=rv, nmatch=nmatch,
                   nmismatch=nmis, nindel=nind, prefer=prefer, accept=accept)
    for col in propagateCol
        c = Symbol(col)
        if c in propertynames(dF.clustering)
            df[!, Symbol("F.", col)] = dF.clustering[fw, c]
            df[!, Symbol("R.", col)] = dR.clustering[rv, c]
        end
    end
    df = df[sortperm(df.abundance; rev=true, alg=Base.Sort.DEFAULT_STABLE), :]
    verbose && println("$(sum(df.abundance[df.accept])) paired-reads (in $(count(df.accept)) unique pairings) successfully merged out of $(sum(df.abundance)) (in $(nrow(df)) pairings) input.")
    returnRejects || (df = df[df.accept, :])
    any(nonunique(df[:, [:sequence]])) && @info "Duplicate sequences in merged output."
    return df
end
