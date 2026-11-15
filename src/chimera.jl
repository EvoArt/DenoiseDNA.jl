# Bimera (two-parent chimera) detection (port of dada2 chimeras.R / chimera.cpp).

# Left/right overlap of a query (al0) with a parent (al1) (port of get_lr).
function _get_lr(al0, al1, len::Int, allow_one_off::Bool, max_shift::Int)
    at(v, i) = 1 <= i <= len ? v[i] : 0x00
    pos = 1; left = 0
    while at(al0, pos) == GAP && pos <= len
        pos += 1
    end
    while at(al1, pos) == GAP && pos - 1 < max_shift
        pos += 1; left += 1
    end
    while pos <= len && al0[pos] == al1[pos]
        pos += 1; left += 1
    end
    left_oo = 0
    if allow_one_off
        left_oo = left
        pos += 1
        (pos <= len && al0[pos] != GAP) && (left_oo += 1)
        while pos <= len && al0[pos] == al1[pos]
            pos += 1; left_oo += 1
        end
    end
    pos = len; right = 0
    while at(al0, pos) == GAP && pos >= 1
        pos -= 1
    end
    while at(al1, pos) == GAP && pos - 1 > len - max_shift
        pos -= 1; right += 1
    end
    while pos >= 1 && al0[pos] == al1[pos]
        pos -= 1; right += 1
    end
    right_oo = 0
    if allow_one_off
        right_oo = right
        pos -= 1
        (pos >= 1 && al0[pos] != GAP) && (right_oo += 1)
        while pos >= 1 && al0[pos] == al1[pos]
            pos -= 1; right_oo += 1
        end
    end
    return left, right, left_oo, right_oo
end

# Hamming distance between aligned sequences ignoring end gaps.
function _ham_endsfree(al0, al1, len::Int)
    i = 1
    g0 = al0[i] == GAP; g1 = al1[i] == GAP
    while (g0 || g1) && i < len
        i += 1
        g0 = g0 && al0[i] == GAP; g1 = g1 && al1[i] == GAP
    end
    j = len
    g0 = al0[j] == GAP; g1 = al1[j] == GAP
    while (g0 || g1) && j > 1
        j -= 1
        g0 = g0 && al0[j] == GAP; g1 = g1 && al1[j] == GAP
    end
    ham = 0
    for p in i:j
        al0[p] != al1[p] && (ham += 1)
    end
    return ham
end

struct _ChimParams
    allow_one_off::Bool
    min_one_off_par_dist::Int
    match::Int
    mismatch::Int
    gap_p::Int
    max_shift::Int
end

function _align_lr!(buf::AlignBuffer, sq::Vector{UInt8}, par::Vector{UInt8}, C::_ChimParams)
    nwalign_vectorized!(buf, sq, par, C.match, C.mismatch, C.gap_p, 0, C.max_shift)
    left, right, loo, roo = _get_lr(buf.al0, buf.al1, buf.len, C.allow_one_off, C.max_shift)
    allowed = C.allow_one_off && _ham_endsfree(buf.al0, buf.al1, buf.len) >= C.min_one_off_par_dist
    return left, right, loo, roo, allowed
end

"""
Overlaps of query `sq` with each parent `bytes[k]`, `k ∈ ks`, computed with the
batched aligner (parents grouped by length). Calls `f(k, left, right, loo, roo, allowed)`.
"""
function _foreach_lr(f, vb::VecBuffer, buf::AlignBuffer, sq, bytes, ks, C::_ChimParams)
    order = sort(ks; by=k -> length(bytes[k]), alg=Base.Sort.DEFAULT_STABLE)
    W = lanes(vb)
    lo = 1
    while lo <= length(order)
        len = length(bytes[order[lo]])
        hi = lo
        while hi < length(order) && hi - lo + 1 < W && length(bytes[order[hi+1]]) == len
            hi += 1
        end
        batch = view(order, lo:hi)
        load!(vb, sq, [bytes[k] for k in batch])
        fill_dp!(vb, C.match, C.mismatch, C.gap_p, 0, C.max_shift)
        for (lane, k) in enumerate(batch)
            traceback!(buf, vb, lane)
            left, right, loo, roo = _get_lr(buf.al0, buf.al1, buf.len, C.allow_one_off, C.max_shift)
            allowed = C.allow_one_off && _ham_endsfree(buf.al0, buf.al1, buf.len) >= C.min_one_off_par_dist
            f(k, left, right, loo, roo, allowed)
        end
        lo = hi + 1
    end
end

"""
    isBimera(sq, parents; allowOneOff=false, minOneOffParentDistance=4, maxShift=16)

Is `sq` an exact bimera of two of the `parents`?
"""
function isBimera(sq::AbstractString, parents::AbstractVector{<:AbstractString}; allowOneOff::Bool=false,
                  minOneOffParentDistance::Integer=4, maxShift::Integer=16,
                  MATCH::Integer=getDadaOpt(:MATCH), MISMATCH::Integer=getDadaOpt(:MISMATCH),
                  GAP_PENALTY::Integer=getDadaOpt(:GAP_PENALTY))
    C = _ChimParams(allowOneOff, minOneOffParentDistance, MATCH, MISMATCH, GAP_PENALTY, maxShift)
    return _is_bimera(Vector{UInt8}(codeunits(sq)), [Vector{UInt8}(codeunits(p)) for p in parents], C, AlignBuffer())
end

_is_bimera(sq::Vector{UInt8}, pars::Vector{Vector{UInt8}}, C::_ChimParams, buf::AlignBuffer) =
    _is_bimera_lr(length(sq), (_align_lr!(buf, sq, par, C) for par in pars), C)

# Bimera test of a query of length L from the overlaps with its parents, in
# the parents' order (port of C_is_bimera).
function _is_bimera_lr(L, overlaps, C::_ChimParams)
    max_left = 0; max_right = 0
    oo_l = 0; oo_r = 0; oo_loo = 0; oo_roo = 0
    for (left, right, loo, roo, allowed) in overlaps
        left + right >= L && continue          # identical / pure shift / internal indel
        max_left = max(max_left, left); max_right = max(max_right, right)
        if allowed
            oo_l = max(oo_l, left); oo_r = max(oo_r, right)
            oo_loo = max(oo_loo, loo); oo_roo = max(oo_roo, roo)
        end
        max_left + max_right >= L && return true
        C.allow_one_off && (oo_l + oo_roo >= L || oo_loo + oo_r >= L) && return true
    end
    return false
end

"""
    isBimeraDenovo(unqs; minFoldParentOverAbundance=2, minParentAbundance=8,
                   allowOneOff=false, minOneOffParentDistance=4, maxShift=16,
                   multithread=false, verbose=false, progress) -> Vector{Bool}

Flag bimeras among unique sequences, using more abundant sequences as
potential parents (R's `isBimeraDenovo`).

`progress=false` hides the progress bar (default: [`DenoiseDNA.setprogress!`](@ref)).
"""
function isBimeraDenovo(unqs; minFoldParentOverAbundance::Real=2, minParentAbundance::Real=8,
                        allowOneOff::Bool=false, minOneOffParentDistance::Integer=4,
                        maxShift::Integer=16, multithread=false, verbose::Bool=false,
                        progress::Bool=_progress_default(),
                        MATCH::Integer=getDadaOpt(:MATCH), MISMATCH::Integer=getDadaOpt(:MISMATCH),
                        GAP_PENALTY::Integer=getDadaOpt(:GAP_PENALTY))
    u = getUniques(unqs)
    seqs_in = _sequences(unqs)
    allunique(seqs_in) || @info "Duplicate sequences detected."
    seqs = u.sequence; ab = u.abundance
    C = _ChimParams(allowOneOff, minOneOffParentDistance, MATCH, MISMATCH, GAP_PENALTY, maxShift)
    bytes = [Vector{UInt8}(codeunits(s)) for s in seqs]
    bims = falses(length(seqs))
    nt = _nthreads(multithread)
    bar = _progress(length(seqs), "Chimeras"; unit=" sequences", enabled=progress)
    Threads.@threads :static for t in 1:nt
        buf = AlignBuffer(); vb = VecBuffer{LANES}()
        lr = Dict{Int,Tuple{Int,Int,Int,Int,Bool}}()
        for i in t:nt:length(seqs)
            pars = [k for k in eachindex(seqs) if ab[k] > minFoldParentOverAbundance * ab[i] && ab[k] > minParentAbundance]
            if length(pars) >= 2
                # align to all parents in batches, then scan them in R's order
                empty!(lr)
                _foreach_lr(vb, buf, bytes[i], bytes, pars, C) do k, l, r, loo, roo, al
                    lr[k] = (l, r, loo, roo, al)
                end
                bims[i] = _is_bimera_lr(length(bytes[i]), (lr[k] for k in pars), C)
            end
            _next!(bar)
        end
    end
    _finish!(bar; info="$(_fmt(count(bims))) bimeras")
    bimset = Set(seqs[bims])
    out = [s in bimset for s in seqs_in]
    verbose && println("Identified $(count(out)) bimeras out of $(length(out)) input sequences.")
    return out
end

_sequences(x::SequenceTable) = x.sequences
_sequences(x::DadaResult) = x.sequence
_sequences(x::Derep) = x.sequence
_sequences(x::AbstractDataFrame) = String.(x.sequence)
_sequences(x::NamedTuple) = String.(x.sequence)

"""
    isBimeraDenovoTable(st; minSampleFraction=0.9, ignoreNNegatives=1,
                        minFoldParentOverAbundance=1.5, minParentAbundance=2,
                        allowOneOff=false, minOneOffParentDistance=4, maxShift=16,
                        multithread=false, verbose=false, progress) -> Vector{Bool}

Consensus chimera detection across the samples of a sequence table: a
sequence is a bimera if flagged in a sufficient fraction of the samples in
which it occurs (R's `isBimeraDenovoTable`).

`progress=false` hides the progress bar (default: [`DenoiseDNA.setprogress!`](@ref)).
"""
function isBimeraDenovoTable(st::SequenceTable; minSampleFraction::Real=0.9, ignoreNNegatives::Integer=1,
                             minFoldParentOverAbundance::Real=1.5, minParentAbundance::Real=2,
                             allowOneOff::Bool=false, minOneOffParentDistance::Integer=4,
                             maxShift::Integer=16, multithread=false, verbose::Bool=false,
                             progress::Bool=_progress_default(),
                             MATCH::Integer=getDadaOpt(:MATCH), MISMATCH::Integer=getDadaOpt(:MISMATCH),
                             GAP_PENALTY::Integer=getDadaOpt(:GAP_PENALTY))
    allunique(st.sequences) || error("Duplicate sequences detected in input.")
    C = _ChimParams(allowOneOff, minOneOffParentDistance, MATCH, MISMATCH, GAP_PENALTY, maxShift)
    mat = st.table
    nrow, ncol = size(mat)
    bytes = [Vector{UInt8}(codeunits(s)) for s in st.sequences]
    nflag = zeros(Int, ncol); nsam = zeros(Int, ncol)
    nt = _nthreads(multithread)
    bar = _progress(ncol, "Chimeras"; unit=" ASVs", enabled=progress)
    isparent(i, j, k) = mat[i, k] > minFoldParentOverAbundance * mat[i, j] && mat[i, k] >= minParentAbundance
    Threads.@threads :static for t in 1:nt
        buf = AlignBuffer(); vb = VecBuffer{LANES}()
        lefts = fill(-1, ncol); rights = fill(-1, ncol)
        lefts_oo = fill(-1, ncol); rights_oo = fill(-1, ncol); allowed = falses(ncol)
        cand = falses(ncol)
        for j in t:nt:ncol
            sqlen = length(bytes[j])
            fill!(lefts, -1); fill!(rights, -1)
            if allowOneOff
                fill!(lefts_oo, -1); fill!(rights_oo, -1); fill!(allowed, false)
            end
            # every potential parent in any sample, aligned once in batches
            # (R aligns the same pairs lazily, caching them the same way)
            fill!(cand, false)
            for i in 1:nrow
                mat[i, j] > 0 || continue
                for k in 1:ncol
                    isparent(i, j, k) && (cand[k] = true)
                end
            end
            ks = findall(cand)
            _foreach_lr(vb, buf, bytes[j], bytes, ks, C) do k, left, right, loo, roo, al
                al && (allowed[k] = true)
                if left + right < sqlen
                    lefts[k] = left; rights[k] = right
                    allowOneOff && (lefts_oo[k] = loo; rights_oo[k] = roo)
                else
                    lefts[k] = 0; rights[k] = 0
                    allowOneOff && (lefts_oo[k] = 0; rights_oo[k] = 0)
                end
            end
            ns = 0; nf = 0
            for i in 1:nrow
                v = mat[i, j]
                v <= 0 && continue
                ns += 1
                ml = 0; mr = 0; ol = 0; or = 0; oloo = 0; oroo = 0
                for k in 1:ncol
                    if isparent(i, j, k)
                        ml = max(ml, lefts[k]); mr = max(mr, rights[k])
                        if allowOneOff && allowed[k]
                            ol = max(ol, lefts[k]); or = max(or, rights[k])
                            oloo = max(oloo, lefts_oo[k]); oroo = max(oroo, rights_oo[k])
                        end
                    end
                end
                if ml + mr >= sqlen
                    nf += 1
                elseif allowOneOff && (ol + oroo >= sqlen || oloo + or >= sqlen)
                    nf += 1
                end
            end
            nflag[j] = nf; nsam[j] = ns
            _next!(bar)
        end
    end
    bims = [nflag[j] >= nsam[j] || (nflag[j] > 0 && nflag[j] >= (nsam[j] - ignoreNNegatives) * minSampleFraction)
            for j in 1:ncol]
    _finish!(bar; info="$(_fmt(count(bims))) bimeras · $(_pct(sum(st.table[:, bims]) / max(sum(st.table), 1))) of reads")
    verbose && println("Identified $(count(bims)) bimeras out of $(ncol) input sequences.")
    return bims
end

"""
    removeBimeraDenovo(x; method="consensus", verbose=false, progress, kwargs...)

Remove bimeras. For a [`SequenceTable`](@ref), `method` is `"consensus"`
(default), `"pooled"` or `"per-sample"`. For per-sample inputs (a dada result,
derep or `sequence`/`abundance` DataFrame) the bimeras are dropped directly.
"""
function removeBimeraDenovo(st::SequenceTable; method::AbstractString="consensus", verbose::Bool=false,
                            progress::Bool=_progress_default(), kwargs...)
    if method == "consensus"
        bim = isBimeraDenovoTable(st; verbose=verbose, progress=progress, kwargs...)
        return st[:, .!bim]
    elseif method == "pooled"
        bim = isBimeraDenovo(st; verbose=verbose, progress=progress, kwargs...)
        return st[:, .!bim]
    elseif method == "per-sample"
        tab = copy(st.table)
        bar = _progress(size(tab, 1), "Chimeras"; unit=" samples", enabled=progress)
        for i in axes(tab, 1)
            row = (sequence=st.sequences, abundance=tab[i, :])
            nz = row.abundance .> 0
            bim = falses(length(st.sequences))
            bim[nz] = isBimeraDenovo((sequence=st.sequences[nz], abundance=row.abundance[nz]); verbose=verbose, progress=false, kwargs...)
            tab[i, bim] .= 0
            _next!(bar)
        end
        _finish!(bar)
        keep = vec(sum(tab; dims=1)) .> 0
        return SequenceTable(st.samples, st.sequences[keep], tab[:, keep])
    else
        error("Valid values for method: \"pooled\", \"consensus\", or \"per-sample\"")
    end
end
function removeBimeraDenovo(df::AbstractDataFrame; verbose::Bool=false, kwargs...)
    bim = isBimeraDenovo(df; verbose=verbose, kwargs...)
    return df[.!bim, :]
end
function removeBimeraDenovo(x::Union{DadaResult,Derep}; verbose::Bool=false, kwargs...)
    u = getUniques(x)
    bim = isBimeraDenovo(x; verbose=verbose, kwargs...)
    return (sequence=u.sequence[.!bim], abundance=u.abundance[.!bim])
end

# ------------------------------------------------------------------ shifts

"""
    isShiftedPair(sq1, sq2; minOverlap=20, flagSubseqs=false) -> Bool

Whether `sq1` and `sq2` agree exactly (no mismatches or internal indels)
over an overlap of at least `minOverlap` bases, with each one sticking out
past the other. With `flagSubseqs`, one may lie wholly inside the other.
"""
function isShiftedPair(sq1::AbstractString, sq2::AbstractString; minOverlap::Integer=20,
                       flagSubseqs::Bool=false)
    a0, a1 = nwalign(sq1, sq2; band=-1)
    ev = eval_pair(Vector{UInt8}(codeunits(a0)), Vector{UInt8}(codeunits(a1)))
    return (ev.match < length(sq1) || flagSubseqs) && (ev.match < length(sq2) || flagSubseqs) &&
           ev.match >= minOverlap && ev.mismatch == 0 && ev.indel == 0
end

"""
    isShift(sq, pars; minOverlap=20, flagSubseqs=false) -> Bool

Whether `sq` is a shifted version of any of the sequences `pars` (see
[`isShiftedPair`](@ref)).
"""
isShift(sq::AbstractString, pars; minOverlap::Integer=20, flagSubseqs::Bool=false) =
    any(par -> isShiftedPair(sq, par; minOverlap, flagSubseqs), pars)

"""
    isShiftDenovo(unqs; minOverlap=20, flagSubseqs=false, verbose=false) -> Vector{Bool}

Flag sequences that are shifted versions of a more abundant sequence in
`unqs`. The result follows the input order (`getSequences(unqs)`).

As in R, `flagSubseqs` is accepted but not passed on, so subsequences are
never flagged.
"""
function isShiftDenovo(unqs; minOverlap::Integer=20, flagSubseqs::Bool=false, verbose::Bool=false)
    u = getUniques(unqs)
    shifted = Set{String}()
    for (sq, a) in zip(u.sequence, u.abundance)
        pars = u.sequence[u.abundance .> a]
        if isempty(pars)
            verbose && println("No possible parents.")
        elseif isShift(sq, pars; minOverlap)
            push!(shifted, sq)
        end
    end
    return [s in shifted for s in getSequences(unqs)]
end
