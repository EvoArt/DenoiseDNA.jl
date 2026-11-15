# High-level sample inference (port of dada2::dada and C++ dada_uniques).

"""
    DadaOptions(; kwargs...)

Algorithm options, named exactly as in R's `setDadaOpt`/`getDadaOpt`
(e.g. `OMEGA_A`, `BAND_SIZE`, `USE_QUALS`). Any of these can be passed as
keyword arguments to [`dada`](@ref) and [`learnErrors`](@ref).
"""
Base.@kwdef mutable struct DadaOptions
    OMEGA_A::Float64 = 1e-40
    OMEGA_P::Float64 = 1e-4
    OMEGA_C::Float64 = 1e-40
    DETECT_SINGLETONS::Bool = false
    USE_KMERS::Bool = true
    KDIST_CUTOFF::Float64 = 0.42
    MAX_CONSIST::Int = 10
    MATCH::Int = 5
    MISMATCH::Int = -4
    GAP_PENALTY::Int = -8
    BAND_SIZE::Int = 16
    VECTORIZED_ALIGNMENT::Bool = true
    MAX_CLUST::Int = 0
    MIN_FOLD::Float64 = 1.0
    MIN_HAMMING::Int = 1
    MIN_ABUNDANCE::Int = 1
    USE_QUALS::Bool = true
    HOMOPOLYMER_GAP_PENALTY::Union{Nothing,Int} = nothing
    GAPLESS::Bool = true
    GREEDY::Bool = true
    PSEUDO_PREVALENCE::Int = 2
    PSEUDO_ABUNDANCE::Float64 = Inf
end

Base.copy(o::DadaOptions) = DadaOptions((getfield(o, f) for f in fieldnames(DadaOptions))...)

function _setopts!(o::DadaOptions, kwargs)
    for (k, v) in kwargs
        hasfield(DadaOptions, k) || throw(ArgumentError("$k is not a valid DADA option"))
        setfield!(o, k, convert(fieldtype(DadaOptions, k), v))
    end
    return o
end

# Session defaults, as R's `dada_opts` environment.
const _DADA_OPTS = DadaOptions()

_options(; kwargs...) = _setopts!(copy(_DADA_OPTS), kwargs)

"""
    setDadaOpt(; kwargs...)
    setDadaOpt(opts::DadaOptions)

Change the session defaults of the [`DadaOptions`](@ref) used by
[`dada`](@ref), [`learnErrors`](@ref), [`nwalign`](@ref) and the chimera
functions, e.g. `setDadaOpt(OMEGA_A=1e-20, BAND_SIZE=32)`. Keywords passed to
those functions still take precedence.
"""
setDadaOpt(; kwargs...) = (_setopts!(_DADA_OPTS, kwargs); nothing)
setDadaOpt(o::DadaOptions) = setDadaOpt(; (f => getfield(o, f) for f in fieldnames(DadaOptions))...)

"""
    getDadaOpt() -> DadaOptions
    getDadaOpt(name::Symbol)

The current session defaults (a copy), or the value of one option.
"""
getDadaOpt() = copy(_DADA_OPTS)
getDadaOpt(name::Symbol) = getfield(_DADA_OPTS, name)

"""
    DadaResult

Output of [`dada`](@ref) for one sample.

- `sequence`, `abundance`: the denoised sequences (ASVs) and their read counts.
- `clustering`: per-cluster diagnostics (as R's `\$clustering`).
- `map`: for each input unique, the index of its denoised sequence (0 = not corrected).
- `pval`: final abundance p-value of each input unique within its cluster.
- `trans`: observed transition counts by quality (16 x nq).
- `quality`: mean quality of each denoised sequence (nclust x maxlen).
- `birth_subs`: substitutions that drove the birth of each cluster.
- `err_in`, `err_out`: error model(s) used and estimated.
"""
struct DadaResult
    sequence::Vector{String}
    abundance::Vector{Int}
    clustering::DataFrame
    map::Vector{Int}
    pval::Vector{Float64}
    trans::Matrix{Int}
    quality::Matrix{Float64}
    birth_subs::DataFrame
    err_in::Any
    err_out::Union{Nothing,Matrix{Float64}}
    opts::DadaOptions
end

const NT_CODE = let t = zeros(UInt8, 256)
    t[UInt8('A')+1] = 1; t[UInt8('C')+1] = 2; t[UInt8('G')+1] = 3; t[UInt8('T')+1] = 4; t[UInt8('N')+1] = 5
    t
end
const NT_CHAR = UInt8['A', 'C', 'G', 'T', 'N']

encode_nt(s::AbstractString) = [NT_CODE[c+1] for c in codeunits(s)]
decode_nt(v::AbstractVector{UInt8}) = String([NT_CHAR[c] for c in v])

isACGT(s::AbstractString) = all(c -> c == UInt8('A') || c == UInt8('C') || c == UInt8('G') || c == UInt8('T'), codeunits(s))

# ------------------------------------------------------------------ core

function _cluster_params(o::DadaOptions; max_clust::Int=o.MAX_CLUST)
    gap = o.GAP_PENALTY > 0 ? -o.GAP_PENALTY : o.GAP_PENALTY
    homo = o.HOMOPOLYMER_GAP_PENALTY === nothing ? gap : o.HOMOPOLYMER_GAP_PENALTY
    homo > 0 && (homo = -homo)
    vec = o.VECTORIZED_ALIGNMENT
    homo != gap && (vec = false)
    o.BAND_SIZE == 0 && (vec = false)
    return ClusterParams(; match=o.MATCH, mismatch=o.MISMATCH, gap_p=gap, homo_gap_p=homo,
        use_kmers=o.USE_KMERS, kdist_cutoff=o.KDIST_CUTOFF, band=o.BAND_SIZE,
        omegaA=o.OMEGA_A, omegaP=o.OMEGA_P, detect_singletons=o.DETECT_SINGLETONS,
        max_clust=max_clust, min_fold=o.MIN_FOLD, min_hamming=o.MIN_HAMMING,
        min_abund=o.MIN_ABUNDANCE, use_quals=true, vectorized_alignment=vec,
        gapless=o.GAPLESS, greedy=o.GREEDY)
end

"""
    dada_uniques(seqs, abundances, priors, err, quals, P, omegaC; nthreads=1)

Run the partitioning on one set of unique sequences and assemble the outputs
(port of the C++ `dada_uniques`). `quals` is `nunique x maxlen`.
"""
function dada_uniques(seqs::Vector{String}, abundances::Vector{Int}, priors::AbstractVector{Bool},
                      err::Matrix{Float64}, quals::Union{Nothing,Matrix{Float64}},
                      P::ClusterParams, omegaC::Float64; nthreads::Int=1, onbud=nothing,
                      cache::Union{Nothing,AlignCache}=nothing)
    nraw = length(seqs)
    nraw == 0 && error("Zero input sequences.")
    size(err, 1) == 16 || error("Error matrix must have 16 rows.")
    lens = length.(seqs)
    maxlen = maximum(lens); minlen = minimum(lens)
    minlen <= KMER_SIZE && error("Input sequences must all be longer than the kmer-size ($KMER_SIZE).")
    has_quals = quals !== nothing && size(quals, 2) > 0

    raws = Vector{Raw}(undef, nraw)
    for i in 1:nraw
        s = encode_nt(seqs[i])
        q = has_quals ? UInt8[round(UInt8, quals[i, p], RoundNearestTiesAway) for p in 1:lens[i]] : UInt8[]
        raws[i] = Raw(s, q, priors[i], lens[i], abundances[i], i, 0.0, -999.0, NO_COMP, false, true)
    end

    b = run_dada(raws, err, P; maxlen=maxlen, nthreads=nthreads, onbud=onbud, cache=cache)
    nc = nclust(b)

    # Final alignments of each raw to its center, and birth alignments.
    ws = Workspace()
    subs = member_subs(b, P, ws)
    birth_subs = Vector{Union{Nothing,Sub}}(nothing, nc)
    for (i, bi) in enumerate(b.bi)
        if i > 1
            from = b.bi[bi.birth_comp.i].center
            ok = sub_new!(ws.sub, ws.buf, b, from, bi.center, P, P.use_kmers, 1.0)
            birth_subs[i] = ok ? copy(ws.sub) : nothing
        end
    end

    # Final within-cluster p-values; uncorrected raws are dropped from outputs.
    pval = zeros(Float64, nraw)
    for bi in b.bi
        for index in bi.raws
            raw = b.raw[index]
            if index == bi.center
                raw.p = 1.0
            else
                raw.p = calc_pA(raw.reads, raw.comp.lambda * bi.reads, true)
                raw.p < omegaC && (raw.correct = false)
            end
            pval[index] = raw.p
        end
    end

    clustering = _clustering_df(b, subs, birth_subs, has_quals)
    trans = _transition_matrix(b, subs, has_quals, size(err, 2))
    cquals = _cluster_quals(b, subs, has_quals, maxlen)
    bsubs = _birth_subs_df(b, birth_subs, has_quals)
    map = zeros(Int, nraw)
    for (i, bi) in enumerate(b.bi), index in bi.raws
        map[index] = b.raw[index].correct ? i : 0
    end
    return (clustering=clustering, birth_subs=bsubs, trans=trans, quality=cquals, map=map, pval=pval)
end

function _clustering_df(b::B, subs, birth_subs, has_quals::Bool)
    nc = nclust(b)
    seqs = Vector{String}(undef, nc)
    abund = zeros(Int, nc); n0 = zeros(Int, nc); n1 = zeros(Int, nc); nunq = zeros(Int, nc)
    birth_from = Vector{Union{Missing,Int}}(missing, nc)
    birth_pval = Vector{Union{Missing,Float64}}(missing, nc)
    birth_fold = Vector{Union{Missing,Float64}}(missing, nc)
    birth_ham = Vector{Union{Missing,Int}}(missing, nc)
    birth_qave = Vector{Union{Missing,Float64}}(missing, nc)
    for (i, bi) in enumerate(b.bi)
        maxr = 0; maxraw = 0
        for index in bi.raws
            if b.raw[index].reads > maxr
                maxraw = index; maxr = b.raw[index].reads
            end
        end
        seqs[i] = maxraw == 0 ? "" : decode_nt(b.raw[maxraw].seq)
        for index in bi.raws
            raw = b.raw[index]
            raw.correct || continue
            abund[i] += raw.reads
            nunq[i] += 1
            s = subs[index]
            if s !== nothing
                s.nsubs == 0 && (n0[i] += raw.reads)
                s.nsubs == 1 && (n1[i] += raw.reads)
            end
        end
        if i > 1
            birth_from[i] = bi.birth_from
            birth_pval[i] = bi.birth_pval
            birth_fold[i] = bi.birth_fold
            birth_ham[i] = bi.birth_comp.hamming
            if has_quals
                q = 0.0
                s = birth_subs[i]
                if s !== nothing
                    center = b.raw[bi.center]
                    for k in 1:s.nsubs
                        q += center.qual[s.map[s.pos[k]]]
                    end
                    q = q / s.nsubs
                end
                birth_qave[i] = q
            end
        end
    end
    # post-hoc p-value of each center given the other clusters
    center_of = Dict(bi.center => i for (i, bi) in enumerate(b.bi))
    tot_e = zeros(Float64, nc)
    for (i, bi) in enumerate(b.bi)
        for c in bi.comp
            j = get(center_of, Int(c.index), 0)
            if j != 0 && i != j
                tot_e[j] += c.lambda * bi.reads
            end
        end
    end
    pv = [calc_pA(b.raw[bi.center].reads, tot_e[i], true) for (i, bi) in enumerate(b.bi)]
    return DataFrame(sequence=seqs, abundance=abund, n0=n0, n1=n1, nunq=nunq, pval=pv,
                     birth_from=birth_from, birth_pval=birth_pval, birth_fold=birth_fold,
                     birth_ham=birth_ham, birth_qave=birth_qave)
end

function _transition_matrix(b::B, subs, has_quals::Bool, ncol::Int)
    has_quals || (ncol = 1)
    T = zeros(Int, 16, ncol)
    for bi in b.bi
        center = b.raw[bi.center]
        for index in bi.raws
            raw = b.raw[index]
            raw.correct || continue
            s = subs[index]
            s === nothing && continue
            @inbounds for pos0 in 1:center.length
                pos1 = s.map[pos0]
                pos1 == GAP_GLYPH && continue
                t = 4 * (Int(center.seq[pos0]) - 1) + Int(raw.seq[pos1])
                q = has_quals ? Int(raw.qual[pos1]) + 1 : 1
                T[t, q] += raw.reads
            end
        end
    end
    return T
end

function _cluster_quals(b::B, subs, has_quals::Bool, maxlen::Int)
    nc = nclust(b)
    Q = zeros(Float64, nc, has_quals ? maxlen : 0)
    has_quals || return Q
    nreads = zeros(Int, maxlen)
    for (i, bi) in enumerate(b.bi)
        len = b.raw[bi.center].length
        fill!(nreads, 0)
        for index in bi.raws
            raw = b.raw[index]
            raw.correct || continue
            s = subs[index]
            s === nothing && continue
            for pos0 in 1:len
                pos1 = s.map[pos0]
                pos1 == GAP_GLYPH && continue
                nreads[pos0] += raw.reads
                Q[i, pos0] += raw.qual[pos1] * raw.reads
            end
        end
        for pos0 in 1:len
            Q[i, pos0] /= nreads[pos0]
        end
        Q[i, len+1:end] .= NaN
    end
    return Q
end

function _birth_subs_df(b::B, birth_subs, has_quals::Bool)
    pos = Int[]; ref = String[]; sb = String[]; qual = Float64[]; clust = Int[]
    for (i, s) in enumerate(birth_subs)
        s === nothing && continue
        center = b.raw[b.bi[i].center]
        for k in 1:s.nsubs
            push!(pos, s.pos[k]); push!(ref, string(Char(NT_CHAR[s.nt0[k]])))
            push!(sb, string(Char(NT_CHAR[s.nt1[k]])))
            push!(qual, has_quals ? Float64(center.qual[s.map[s.pos[k]]]) : NaN)
            push!(clust, i)
        end
    end
    return DataFrame(pos=pos, ref=ref, sub=sb, qual=qual, clust=clust)
end

# ------------------------------------------------------------------ R-level dada()

"""
    combineDereps(dereps) -> Derep

Pool several samples into one `Derep` (port of R's `combineDereps2`).
"""
function combineDereps(dereps::AbstractVector{Derep})
    maxlen = maximum(d -> size(d.quals, 2), dereps)
    index = Dict{String,Int}()
    seqs = String[]
    for d in dereps, s in d.sequence
        if !haskey(index, s)
            push!(seqs, s); index[s] = length(seqs)
        end
    end
    counts = zeros(Int, length(seqs))
    Q = zeros(Float64, length(seqs), maxlen)
    rmap = Int[]
    for d in dereps
        m = [index[s] for s in d.sequence]
        for (k, u) in enumerate(m)
            counts[u] += d.abundance[k]
            for p in 1:maxlen
                q = p <= size(d.quals, 2) ? d.quals[k, p] : NaN
                Q[u, p] += q * d.abundance[k]
            end
        end
        append!(rmap, [x == 0 ? 0 : m[x] for x in d.map])
    end
    Q ./= counts
    ord = sortperm(counts; rev=true, alg=Base.Sort.DEFAULT_STABLE)
    rank = invperm(ord)
    return Derep(seqs[ord], counts[ord], Q[ord, :], [u == 0 ? 0 : rank[u] for u in rmap])
end

_nthreads(multithread::Bool) = multithread ? Threads.nthreads() : 1
_nthreads(multithread::Integer) = Int(multithread)

"""
    dada(derep, err; selfConsist=false, pool=false, priors=String[], multithread=false,
         errorEstimationFunction=loessErrfun, verbose=false, progress, OPTION=value...)

Infer the exact sequence variants in each sample (port of R's `dada`).

- `derep`: a [`Derep`](@ref), a vector of them, or FASTQ file path(s).
- `err`: 16-row error matrix (e.g. from [`learnErrors`](@ref)), or `nothing`
  with `selfConsist=true`.
- `pool`: `false` (independent), `true` (pool all samples) or `:pseudo`.
- `priors`: sequences with prior evidence of existence.
- Options (`OMEGA_A`, `BAND_SIZE`, ...) are given as keywords; see [`DadaOptions`](@ref).

Returns a [`DadaResult`](@ref) (or a vector of them for multiple samples).

`progress=false` hides the progress bar (default: [`DenoiseDNA.setprogress!`](@ref)).

With several rounds (`selfConsist`, or `pool=:pseudo`) alignments are reused
between rounds; `cache_alignments=false` recomputes them (same results, slower).
"""
function dada(derep, err=nothing; errorEstimationFunction=loessErrfun, selfConsist::Bool=false,
              pool=false, priors::AbstractVector{<:AbstractString}=String[], multithread=false,
              verbose::Bool=false, progress::Bool=_progress_default(), cache_alignments::Bool=true,
              kwargs...)
    opts = _options(; kwargs...)
    single = derep isa Derep || derep isa AbstractString
    dereps = derep isa Derep ? [derep] :
             derep isa AbstractString ? [derepFastq(derep)] :
             eltype(derep) <: AbstractString ? derepFastq(collect(derep); progress=progress) : collect(Derep, derep)
    length(dereps) <= 1 && (pool = false)
    pseudo = false
    pooled_input = nothing
    if pool === true
        pooled_input = dereps
        dereps = [combineDereps(dereps)]
    elseif pool === :pseudo || pool == "pseudo"
        pool = false; pseudo = true
    elseif pool !== false
        throw(ArgumentError("Invalid pool argument."))
    end
    pset = Set{String}(priors)
    pseudo_priors = Set{String}()

    initializeErr = false
    if selfConsist && err === nothing
        initializeErr = true
    else
        err = getErrors(err)
    end
    (0 <= opts.OMEGA_A < 1) || error("OMEGA_A must be between zero and one.")
    (0 <= opts.OMEGA_P < 1) || error("OMEGA_P must be between zero and one.")
    opts.OMEGA_C >= 1 && selfConsist && error("Some error correction required when learning error rates.")
    if !opts.USE_QUALS
        errorEstimationFunction = noqualErrfun
    end
    nthreads = _nthreads(multithread)
    nconsist = initializeErr ? 0 : 1
    errs = Matrix{Float64}[]
    # alignments don't depend on the error model: reuse them across rounds
    caches = cache_alignments && (selfConsist || pseudo) ? [AlignCache() for _ in dereps] : nothing
    local results
    local cur
    while true
        nconsist > 0 && push!(errs, err)
        results = Vector{Any}(undef, length(dereps))
        P = _cluster_params(opts; max_clust=initializeErr ? 1 : opts.MAX_CLUST)
        # samples are independent: run them in parallel, each single-threaded
        sample_threads = nthreads > 1 && length(dereps) > 1
        inner = sample_threads ? 1 : nthreads
        ns = length(dereps)
        desc = initializeErr ? "Initialising" : selfConsist ? "Round $nconsist" :
               pseudo ? "Pseudo-pool $nconsist" : "Denoising"
        # one sample (or a pooled one): count its clusters as they are found
        bar = _progress(ns == 1 ? 0 : ns, desc; unit=ns == 1 ? " ASVs" : " samples", enabled=progress)
        run_one = function (i, onbud)
            d = dereps[i]
            all(isACGT, d.sequence) || error("Invalid derep: sequences must be made up only of A/C/G/T.")
            qmax = opts.USE_QUALS ? ceil(Int, maximum(x -> isnan(x) ? -Inf : x, d.quals)) : 0
            erri = initializeErr ? fill(1.0, 16, max(41, qmax + 1)) : err
            if size(erri, 2) < qmax + 1
                erri = hcat(erri, repeat(erri[:, end], 1, qmax + 1 - size(erri, 2)))
            end
            pri = [s in pset || s in pseudo_priors for s in d.sequence]
            dada_uniques(d.sequence, d.abundance, pri, erri, d.quals, P, opts.OMEGA_C; nthreads=inner, onbud=onbud,
                         cache=caches === nothing ? nothing : caches[i])
        end
        if sample_threads
            Threads.@threads :dynamic for i in eachindex(dereps)
                results[i] = run_one(i, nothing)
                _next!(bar)
            end
        else
            for i in eachindex(dereps)
                onbud = !progress ? nothing : ns == 1 ? (c -> _update!(bar, c)) :
                        (c -> _info!(bar, _note(bar.io, "sample $i: $c ASVs")))
                results[i] = run_one(i, onbud)
                ns == 1 ? _update!(bar, nrow(results[i].clustering)) : _next!(bar)
                verbose && _progress_println(bar, "Sample $i - $(sum(dereps[i].abundance)) reads in $(length(dereps[i])) unique sequences.")
            end
        end
        cur = accumulateTrans([r.trans for r in results])
        newerr = try
            errorEstimationFunction(cur)
        catch e
            (selfConsist || verbose) && @warn "Error rates could not be estimated (this is usually because of very few reads)." exception = e
            nothing
        end
        selfConsist && newerr === nothing && error("Error matrix is NULL.")
        if progress
            nasv = sum(r -> nrow(r.clustering), results)
            info = ns == 1 || initializeErr ? "" : "$(_fmt(nasv)) ASVs"
            if selfConsist && !initializeErr && newerr !== nothing && size(newerr) == size(err)
                if any(e -> e == newerr, errs)
                    info *= _sty(bar.io, "  converged"; color=_PAL.good)
                else
                    δ = maximum(abs(log10(max(a, 1e-300)) - log10(max(b, 1e-300))) for (a, b) in zip(newerr, err))
                    info *= _note(bar.io, "  max Δlog₁₀ err " * _g(δ))
                end
            end
            _finish!(bar; info=info)
        end
        err = newerr
        if initializeErr
            initializeErr = false
            err[[1, 6, 11, 16], :] .= 1.0
        end
        verbose && selfConsist && !progress && print(".")
        if !selfConsist || any(e -> e == err, errs) || nconsist >= opts.MAX_CONSIST
            (!pseudo || nconsist >= 2) && break
        end
        if pseudo && nconsist >= 1
            st = makeSequenceTable([(sequence=r.clustering.sequence, abundance=r.clustering.abundance) for r in results])
            cs = vec(sum(st.table; dims=1)); prev = vec(sum(st.table .> 0; dims=1))
            pseudo_priors = Set(st.sequences[(prev .>= opts.PSEUDO_PREVALENCE) .| (cs .>= opts.PSEUDO_ABUNDANCE)])
        end
        nconsist += 1
    end
    if selfConsist && verbose
        println()
        nconsist >= opts.MAX_CONSIST ? @warn("Self-consistency loop terminated before convergence.") :
            println("Convergence after $nconsist rounds.")
    end

    err_in = selfConsist ? errs : (isempty(errs) ? nothing : errs[1])
    out = [DadaResult(r.clustering.sequence, r.clustering.abundance, r.clustering, r.map, r.pval,
                      r.trans, r.quality, r.birth_subs, err_in, err, opts) for r in results]
    if pooled_input !== nothing
        out = [_unpool(out[1], dereps[1], d) for d in pooled_input]
    end
    return single ? out[1] : out
end

# Split a pooled result back into a per-sample result (as R does for pool=TRUE).
function _unpool(r::DadaResult, pooled::Derep, d::Derep)
    pidx = Dict(s => i for (i, s) in enumerate(pooled.sequence))
    pm = [r.map[pidx[s]] for s in d.sequence]
    keep = falses(length(r.sequence))
    for c in pm
        c > 0 && (keep[c] = true)
    end
    newbi = cumsum(keep)
    map = [c == 0 ? 0 : newbi[c] for c in pm]
    abund = zeros(Int, count(keep))
    for (k, c) in enumerate(map)
        c > 0 && (abund[c] += d.abundance[k])
    end
    cl = r.clustering[keep, :]
    cl.abundance = abund
    bs = r.birth_subs[keep[r.birth_subs.clust], :]
    bs.clust = newbi[bs.clust]
    # R keeps the pooled $pval and $trans unchanged
    return DadaResult(cl.sequence, abund, cl, map, r.pval, r.trans, r.quality[keep, :], bs,
                      r.err_in, r.err_out, r.opts)
end

"Unique sequences and abundances of a dada result (as R's `getUniques`)."
getUniques(r::DadaResult; collapse::Bool=true) = (sequence=r.sequence, abundance=r.abundance)
getUniques(d::Derep; collapse::Bool=true) = (sequence=d.sequence, abundance=d.abundance)
