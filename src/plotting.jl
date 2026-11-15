# Data behind the diagnostic plots. The plots themselves live in the Makie
# extension (ext/DenoiseDNAMakieExt.jl), loaded by `using CairoMakie` (or any
# other Makie backend).

"""
    qualityProfile(file; n=500_000) -> NamedTuple

Per-cycle quality summary of (the first `n` reads of) a FASTQ file, as shown
by [`plotQualityProfile`](@ref): `counts[cycle, q+1]` reads with quality `q`
at each cycle, and per-cycle `mean`, quartiles `q25`/`q50`/`q75`, and `cum`,
the fraction of reads at least that long.
"""
function qualityProfile(file::AbstractString; n::Integer=500_000)
    fq = read_fastq(file)
    m = min(n, length(fq))
    quals = fq.qual[1:m]
    maxlen = maximum(length, quals; init=0)
    maxq = maximum(q -> isempty(q) ? 0 : Int(maximum(q)), quals; init=0)
    counts = zeros(Int, maxlen, maxq + 1)
    for q in quals, (cycle, s) in enumerate(q)
        counts[cycle, s+1] += 1
    end
    scores = 0:maxq
    tot = vec(sum(counts; dims=2))
    mean = [sum(counts[c, :] .* scores) / tot[c] for c in 1:maxlen]
    quant(c, p) = scores[findfirst(>=(p), cumsum(counts[c, :]) ./ tot[c])]
    return (file=basename(file), nreads=length(fq), sampled=m, scores=scores, counts=counts,
            mean=mean, q25=quant.(1:maxlen, 0.25), q50=quant.(1:maxlen, 0.5),
            q75=quant.(1:maxlen, 0.75), cum=tot ./ m)
end

"""
    errorProfile(x) -> DataFrame

Long table of error rates by transition and quality, as shown by
[`plotErrors`](@ref): observed frequencies (`count / tot`), the fitted
(`Estimated`) and input (`Input`) error rates and the `Nominal` rate implied by
the quality score. `x` is an [`ErrorModel`](@ref), a [`DadaResult`](@ref) (or
vector of them) or an error matrix.
"""
function errorProfile(x)
    trans, err_out, err_in = _error_parts(x)
    nq = trans !== nothing ? size(trans, 2) : size(err_out, 2)
    nq <= 1 && error("errorProfile needs an error model fitted with quality scores.")
    rows = NamedTuple[]
    for t in 1:16, q in 0:nq-1
        from = TRANSITION_NAMES[t][1]; to = TRANSITION_NAMES[t][3]
        self = from == to
        count = trans === nothing || self ? missing : trans[t, q+1]
        tot = trans === nothing || self ? missing : sum(trans[4*((t-1)÷4)+1:4*((t-1)÷4)+4, q+1])
        obs = count === missing ? missing : count / tot
        nominal = self ? 1 - 10.0^(-q / 10) : 10.0^(-q / 10) / 3
        push!(rows, (Transition=TRANSITION_NAMES[t], from=string(from), to=string(to), Qual=q,
                     count=count, tot=tot, Observed=obs,
                     Estimated=err_out === nothing || self || q >= size(err_out, 2) ? missing : err_out[t, q+1],
                     Input=err_in === nothing || self || q >= size(err_in, 2) ? missing : err_in[t, q+1],
                     Nominal=nominal))
    end
    return DataFrame(rows)
end

_error_parts(e::ErrorModel) = (e.trans, e.err_out, isempty(e.err_in) ? nothing : e.err_in[1])
_error_parts(r::DadaResult) = (r.trans, r.err_out, r.err_in isa AbstractVector ? first(r.err_in) : r.err_in)
_error_parts(rs::AbstractVector{DadaResult}) =
    (accumulateTrans([r.trans for r in rs]), _error_parts(rs[1])[2], _error_parts(rs[1])[3])
_error_parts(m::AbstractMatrix) = (nothing, Matrix{Float64}(m), nothing)

"""
    plotQualityProfile(files; n=500_000, aggregate=false) -> Figure

Heatmap of quality scores by cycle for each FASTQ file, with the mean
(green), median and quartiles (orange) and, for variable-length reads, the
fraction of reads reaching each cycle (red). Requires a Makie backend
(`using CairoMakie`).
"""
function plotQualityProfile end

"""
    plotErrors(x; nti="ACGT", ntj="ACGT", obs=true, err_out=true, err_in=false, nominalQ=false) -> Figure

Observed (points) and fitted (black line) error rates for each transition as
a function of quality score; `nominalQ=true` adds the rates implied by the
quality scores (red). Requires a Makie backend (`using CairoMakie`).
"""
function plotErrors end

"""
    plotComplexity(files; kmerSize=2, window=nothing, by=5, n=100_000, bins=100, aggregate=false) -> Figure

Histogram of sequence complexity (see [`seqComplexity`](@ref)) of the reads in
each file. Requires a Makie backend (`using CairoMakie`).
"""
function plotComplexity end

function __init__()
    Base.Experimental.register_error_hint(MethodError) do io, exc, _, _
        if exc.f in (plotQualityProfile, plotErrors, plotComplexity)
            print(io, "\nPlotting needs a Makie backend: run `using CairoMakie` (or GLMakie/WGLMakie) first.")
        end
    end
end
