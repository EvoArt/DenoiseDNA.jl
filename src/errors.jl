# Error models (port of dada2 errorModels.R).

const TRANSITION_NAMES = ["A2A", "A2C", "A2G", "A2T", "C2A", "C2C", "C2G", "C2T",
                          "G2A", "G2C", "G2G", "G2T", "T2A", "T2C", "T2G", "T2T"]
const MAX_ERROR_RATE = 0.25
const MIN_ERROR_RATE = 1e-7

# Rows (1-based) of the 12 substitution types, in the order R builds `est`.
const _SUB_ROWS = [2, 3, 4, 5, 7, 8, 9, 10, 12, 13, 14, 15]

# Assemble the 16-row error matrix from the 12 estimated substitution rates.
function _expand_err(est::Matrix{Float64})
    nq = size(est, 2)
    err = Matrix{Float64}(undef, 16, nq)
    for (k, r) in enumerate(_SUB_ROWS)
        err[r, :] = est[k, :]
    end
    for (self, rows) in ((1, 1:3), (6, 4:6), (11, 7:9), (16, 10:12))
        for q in 1:nq
            # R: 1 - colSums(est[rows, ]); colSums accumulates in extended
            # precision, so round the exact sum once.
            s = Float64(sum(big(est[k, q]) for k in rows))
            err[self, q] = 1 - s
        end
    end
    return err
end

"""
    loessErrfun(trans) -> Matrix{Float64}

Estimate error rates from a 16 x nq transition count matrix by a loess fit of
log10 error frequency on quality score, weighted by the number of
observations (identical to R's `loessErrfun`).
"""
function loessErrfun(trans::AbstractMatrix)
    nq = size(trans, 2)
    qq = Float64.(0:nq-1)
    est = Matrix{Float64}(undef, 12, nq)
    row = 0
    for nti in 0:3, ntj in 0:3
        nti == ntj && continue
        row += 1
        errs = Float64.(trans[4nti+ntj+1, :])
        tot = vec(sum(Float64.(trans[4nti+1:4nti+4, :]); dims=1))
        rlogp = log10.((errs .+ 1) ./ tot)
        keep = isfinite.(rlogp)          # R: Inf -> NA, then na.omit
        pred = fill(NaN, nq)
        if any(keep)
            fit = loess(qq[keep], rlogp[keep]; weights=tot[keep], span=0.75, degree=2)
            pred = predict(fit, qq)
        end
        hi = findlast(!isnan, pred); lo = findfirst(!isnan, pred)
        (hi === nothing) && error("loess fit produced no predictions")
        pred[hi+1:end] .= pred[hi]
        pred[1:lo-1] .= pred[lo]
        est[row, :] = 10.0 .^ pred
    end
    clamp!(est, MIN_ERROR_RATE, MAX_ERROR_RATE)
    return _expand_err(est)
end

"""
    noqualErrfun(trans; pseudocount=1)

Error rates that ignore quality scores (used when `USE_QUALS=false`).
"""
function noqualErrfun(trans::AbstractMatrix; pseudocount=1)
    nq = size(trans, 2)
    obs = vec(sum(trans; dims=2)) .+ pseudocount
    est = Matrix{Float64}(undef, 12, nq)
    row = 0
    for nti in 0:3, ntj in 0:3
        nti == ntj && continue
        row += 1
        est[row, :] .= obs[4nti+ntj+1] / sum(obs[4nti+1:4nti+4])
    end
    return _expand_err(est)
end

"""
    makeBinnedQualErrfun(binnedQ) -> function

Error function for binned quality scores (e.g. NovaSeq): linear interpolation
of observed error frequencies between the binned values.
"""
function makeBinnedQualErrfun(binnedQ::AbstractVector{<:Integer})
    isempty(binnedQ) && error("The quality scores used in your data must be provided.")
    binned = sort(collect(Int, binnedQ))
    return function (trans::AbstractMatrix)
        nq = size(trans, 2)
        qq = 0:nq-1
        obsq = qq[vec(sum(trans; dims=1)) .> 0]
        qmax = maximum(obsq); qmin = minimum(obsq)
        qmax > maximum(binned) && error("Input data contains a higher quality score than the provided binned values.")
        qmin < minimum(binned) && error("Input data contains a lower quality score than the provided binned values.")
        qmax in binned || @warn "Maximum observed quality score is not in the provided binned values."
        qmin in binned || @warn "Minimum observed quality score is not in the provided binned values."
        est = Matrix{Float64}(undef, 12, nq)
        row = 0
        for nti in 0:3, ntj in 0:3
            nti == ntj && continue
            row += 1
            errs = Float64.(trans[4nti+ntj+1, :])
            tot = vec(sum(Float64.(trans[4nti+1:4nti+4, :]); dims=1))
            p = errs ./ tot
            pred = fill(NaN, nq)
            for i in 1:length(binned)-1
                loQ = binned[i]; hiQ = binned[i+1]
                loP = p[loQ+1]; hiP = p[hiQ+1]
                if !isnan(loP) && !isnan(hiP)
                    pred[loQ+1:hiQ+1] = range(loP, hiP; length=hiQ - loQ + 1)
                end
            end
            hi = findlast(!isnan, pred); lo = findfirst(!isnan, pred)
            pred[hi+1:end] .= pred[hi]
            pred[1:lo-1] .= pred[lo]
            est[row, :] = pred
        end
        clamp!(est, MIN_ERROR_RATE, MAX_ERROR_RATE)
        return _expand_err(est)
    end
end

"""
    PacBioErrfun(trans)

Error function for PacBio CCS reads: loess on quality scores below 93 and a
direct estimate for the special Q93 bin.
"""
function PacBioErrfun(trans::AbstractMatrix)
    nq = size(trans, 2)
    if nq >= 94
        nq == 94 || error("Max qual score of 93 not the last column as expected.")
        err = loessErrfun(trans[:, 1:93])
        t93 = trans[:, 94]
        tot = [sum(t93[1:4]), sum(t93[5:8]), sum(t93[9:12]), sum(t93[13:16])]
        err93 = (t93 .+ 1) ./ (repeat(tot; inner=4) .+ 4)
        return hcat(err, err93)
    else
        @info "The max qual score of 93 was not detected. Using standard error fitting."
        return loessErrfun(trans)
    end
end

"""
    accumulateTrans(trans) -> Matrix

Sum a vector of 16-row transition count matrices, padding the narrower ones
with zero columns.
"""
function accumulateTrans(trans::AbstractVector)
    maxcol = maximum(size.(trans, 2))
    out = zeros(Int, 16, maxcol)
    for t in trans
        out[:, 1:size(t, 2)] .+= t
    end
    return out
end

"""
    ErrorModel

Output of [`learnErrors`](@ref): `err_out` is the fitted error matrix,
`err_in` the matrices used in each self-consistency round and `trans` the
observed transition counts.
"""
struct ErrorModel
    err_out::Matrix{Float64}
    err_in::Vector{Matrix{Float64}}
    trans::Matrix{Int}
end

"""
    getErrors(obj; enforce=true) -> Matrix{Float64}

Extract the error matrix from a matrix, [`ErrorModel`](@ref) or [`DadaResult`](@ref)(s).
"""
function getErrors(obj; enforce::Bool=true)
    err = obj isa AbstractMatrix ? Matrix{Float64}(obj) :
          obj isa ErrorModel ? obj.err_out :
          obj isa DadaResult ? obj.err_out :
          obj isa AbstractVector{<:DadaResult} ? obj[1].err_out :
          obj === nothing ? nothing : error("Unrecognized error matrix input.")
    if enforce
        err === nothing && error("Error matrix is NULL.")
        size(err, 1) == 16 || error("Error matrix must have 16 rows (A2A, A2C, ...).")
        all(>=(0), err) || error("All error matrix entries must be >= 0.")
        all(<=(1), err) || error("All error matrix entries must be <=1.")
        any(==(0), err) && @warn "Zero in error matrix."
    end
    return err
end

"""
    inflateErr(err, inflation; inflateSelfTransitions=false)

Inflate error rates by `inflation`, keeping them below 1.
"""
function inflateErr(err, inflation::Real; inflateSelfTransitions::Bool=false)
    e = copy(getErrors(err))
    rows = inflateSelfTransitions ? (1:16) : _SUB_ROWS
    for r in rows
        e[r, :] .= (e[r, :] .* inflation) ./ (1 .+ (inflation - 1) .* e[r, :])
    end
    return e
end

"""
    learnErrors(files_or_dereps; nbases=1e8, errorEstimationFunction=loessErrfun,
                multithread=false, randomize=false, MAX_CONSIST=10, OMEGA_C=0,
                verbose=false, progress, OPTIONS...) -> ErrorModel

Learn the error model by alternating sample inference and error rate
estimation until convergence, using samples until `nbases` bases have been
read (port of R's `learnErrors`). The progress display shows one line per
self-consistency round, with the change in the error rates since the last one.

`progress=false` hides the progress bar (default: [`DenoiseDNA.setprogress!`](@ref)).
"""
function learnErrors(fls; nbases::Real=1e8, errorEstimationFunction=loessErrfun,
                     multithread=false, randomize::Bool=false, MAX_CONSIST::Integer=10,
                     OMEGA_C::Real=0, verbose::Bool=false, progress::Bool=_progress_default(), kwargs...)
    items = fls isa Derep || fls isa AbstractString ? [fls] : collect(fls)
    randomize && (items = items[Random.randperm(length(items))])
    drps = Derep[]
    NBASES = 0.0; NREADS = 0.0
    bar = _progress(0, "Reading"; unit=" samples", enabled=progress && any(it -> it isa AbstractString, items))
    for it in items
        d = it isa Derep ? it : derepFastq(it)
        push!(drps, d)
        NREADS += sum(d.abundance)
        NBASES += sum(Float64.(d.abundance) .* length.(d.sequence))
        _next!(bar; info=_note(bar.io, "$(_g(NBASES)) of $(_g(Float64(nbases))) bases"))
        NBASES > nbases && break
    end
    _finish!(bar; info="$(_g(NBASES)) bases in $(_fmt(Int(NREADS))) reads")
    verbose && println("$(Int(NBASES)) total bases in $(Int(NREADS)) reads from $(length(drps)) samples will be used for learning the error rates.")
    dds = dada(drps, nothing; errorEstimationFunction=errorEstimationFunction, selfConsist=true,
               multithread=multithread, verbose=verbose, progress=progress, MAX_CONSIST=MAX_CONSIST, OMEGA_C=OMEGA_C, kwargs...)
    dds isa DadaResult && (dds = [dds])
    return ErrorModel(dds[1].err_out, dds[1].err_in, accumulateTrans([d.trans for d in dds]))
end

# ------------------------------------------------------------------ bad bases

"""
    getBadBases(clust, birth_subs; omegaB=1e-20, minOccurence=4, verbose=false) -> Vector{Int}

Positions where the one-off clusters of a [`dada`](@ref) result (`clust =
r.clustering`, `birth_subs = r.birth_subs`) were born from more substitutions
than chance allows: a Poisson upper tail, Bonferroni-corrected over
positions, below `omegaB`, and at least `minOccurence` substitutions.
Requires equal-length sequences.
"""
function getBadBases(clust::AbstractDataFrame, birth_subs::AbstractDataFrame; omegaB::Real=1e-20,
                     minOccurence::Integer=4, verbose::Bool=false)
    oos = findall(isequal(1), clust.birth_ham)
    counts = Dict{Int,Int}()
    for (p, c) in zip(birth_subs.pos, birth_subs.clust)
        c in oos && (counts[p] = get(counts, p, 0) + 1)
    end
    allequal(length.(clust.sequence)) || error("Requires same length sequences.")
    seqlen = length(clust.sequence[1])
    bad = sort!([p for (p, n) in counts
                 if ppois_upper(n, length(oos) / seqlen) * seqlen < omegaB && n >= minOccurence])
    verbose && println(length(bad), " bad bases identified.")
    return bad
end

"""
    isBadBaseFP(clust, birth_subs; minFraction=0.51, omegaB=1e-10, minOccurence=4, verbose=false) -> Vector{Bool}

Flag clusters born mostly (a fraction of at least `minFraction` of their
birth substitutions) from bad bases found by [`getBadBases`](@ref).
"""
function isBadBaseFP(clust::AbstractDataFrame, birth_subs::AbstractDataFrame; minFraction::Real=0.51,
                     omegaB::Real=1e-10, minOccurence::Integer=4, verbose::Bool=false)
    bb = Set(getBadBases(clust, birth_subs; omegaB, minOccurence, verbose))
    nsub = zeros(Int, nrow(clust)); nbad = zeros(Int, nrow(clust))
    for (p, c) in zip(birth_subs.pos, birth_subs.clust)
        nsub[c] += 1
        nbad[c] += p in bb
    end
    fp = [n > 0 && b / n >= minFraction for (n, b) in zip(nsub, nbad)]
    verbose && println(sum(fp), " false positives caused by bad bases identified from ", nrow(clust), " input sequences.")
    return fp
end
