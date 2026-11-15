module DenoiseDNAMakieExt

using DenoiseDNA
using DenoiseDNA: qualityProfile, errorProfile, read_fastq, seqComplexity
using Makie

const GREEN = colorant"#66C2A5"
const ORANGE = colorant"#FC8D62"
const TILE = cgrad([colorant"#F5F5F5", colorant"black"])

_files(x::AbstractString) = [x]
_files(x) = collect(x)
_grid(n) = (ncol = ceil(Int, sqrt(n)); (ceil(Int, n / ncol), ncol))

function _quality_panel!(ax, pos, p; label)
    cycles = 1:size(p.counts, 1)
    z = Float64.(p.counts)
    z[z .== 0] .= NaN                      # as ggplot, no tile where there are no reads
    heatmap!(ax, cycles, p.scores, z; colormap=TILE, nan_color=:transparent)
    lines!(ax, cycles, p.mean; color=GREEN)
    lines!(ax, cycles, p.q25; color=ORANGE, linewidth=0.6, linestyle=:dash)
    lines!(ax, cycles, p.q50; color=ORANGE, linewidth=0.6)
    lines!(ax, cycles, p.q75; color=ORANGE, linewidth=0.6, linestyle=:dash)
    if length(unique(p.cum)) > 1
        lines!(ax, cycles, 10 .* p.cum; color=:red, linewidth=0.6)
        # R's secondary axis: the red line's 0-10 is 0-100% of reads
        ax2 = Axis(pos; yaxisposition=:right, yticks=([0, 10], ["0%", "100%"]), yticklabelcolor=:red, ytickcolor=:red,
                   backgroundcolor=:transparent)
        hidexdecorations!(ax2); hidespines!(ax2); ax2.ygridvisible = false
        linkaxes!(ax, ax2)
    end
    text!(ax, 0, 0; text=label, color=:red, align=(:left, :bottom), fontsize=11)
    ylims!(ax, 0, nothing)
end

function DenoiseDNA.plotQualityProfile(files; n::Integer=500_000, aggregate::Bool=false,
                                       size=nothing)
    fl = _files(files)
    profs = [qualityProfile(f; n=n) for f in fl]
    if aggregate
        # (`size` is the figure-size keyword here, hence Base.size)
        maxc = maximum(p -> Base.size(p.counts, 1), profs); maxq = maximum(p -> Base.size(p.counts, 2), profs)
        counts = zeros(Int, maxc, maxq)
        for p in profs
            counts[1:Base.size(p.counts, 1), 1:Base.size(p.counts, 2)] .+= p.counts
        end
        scores = 0:maxq-1
        tot = vec(sum(counts; dims=2))
        quant(c, q) = scores[findfirst(>=(q), cumsum(counts[c, :]) ./ tot[c])]
        agg = (counts=counts, scores=scores, mean=[sum(counts[c, :] .* scores) / tot[c] for c in 1:maxc],
               q25=quant.(1:maxc, 0.25), q50=quant.(1:maxc, 0.5), q75=quant.(1:maxc, 0.75),
               cum=tot ./ sum(p.sampled for p in profs))
        fig = Figure(; size=something(size, (600, 450)))
        ax = Axis(fig[1, 1]; title="$(length(fl)) files (aggregated)", xlabel="Cycle", ylabel="Quality Score")
        _quality_panel!(ax, fig[1, 1], agg; label="Total reads: $(sum(p.nreads for p in profs))")
        return fig
    end
    nr, nc = _grid(length(profs))
    fig = Figure(; size=something(size, (380nc, 300nr)))
    for (k, p) in enumerate(profs)
        r, c = fldmod1(k, nc)
        ax = Axis(fig[r, c]; title=p.file, titlesize=12,
                  xlabel=r == nr ? "Cycle" : "", ylabel=c == 1 ? "Quality Score" : "")
        lab = p.nreads >= n ? "Reads >= $n" : "Reads: $(p.nreads)"
        _quality_panel!(ax, fig[r, c], p; label=lab)
    end
    return fig
end

function DenoiseDNA.plotErrors(x; nti="ACGT", ntj="ACGT", obs::Bool=true, err_out::Bool=true,
                               err_in::Bool=false, nominalQ::Bool=false, size=nothing)
    df = errorProfile(x)
    from = collect(string.(collect(nti))); to = collect(string.(collect(ntj)))
    all(in(("A", "C", "G", "T")), vcat(from, to)) || error("nti and ntj must be nucleotide(s): A/C/G/T.")
    fig = Figure(; size=something(size, (220length(to), 200length(from))))
    # shared scales, as R's facets (empty panels get the same axes)
    shown = Float64[]
    obs && append!(shown, filter(>(0), Float64.(collect(skipmissing(df.Observed)))))
    err_out && append!(shown, filter(>(0), Float64.(collect(skipmissing(df.Estimated)))))
    err_in && append!(shown, filter(>(0), Float64.(collect(skipmissing(df.Input)))))
    nominalQ && append!(shown, filter(>(0), Float64.(collect(skipmissing(df.Nominal)))))
    ylim = isempty(shown) ? (1e-4, 1.0) : (minimum(shown) / 1.5, maximum(shown) * 1.5)
    xlim = (-1.0, maximum(df.Qual) + 1.0)
    for (r, f) in enumerate(from), (c, t) in enumerate(to)
        sub = df[(df.from .== f) .& (df.to .== t), :]
        ax = Axis(fig[r, c]; title="$(f)2$(t)", titlesize=12, yscale=log10,
                  xlabel=r == length(from) ? "Consensus quality score" : "",
                  ylabel=c == 1 ? "Error frequency (log10)" : "")
        limits!(ax, xlim..., ylim...)
        f == t && !nominalQ && continue
        if obs
            ok = .!ismissing.(sub.Observed) .& (coalesce.(sub.Observed, 0.0) .> 0)
            if any(ok)
                cnt = Float64.(sub.count[ok])
                ms = 2 .+ 6 .* sqrt.(cnt ./ maximum(cnt))     # area ∝ count, like scale_size_area
                scatter!(ax, sub.Qual[ok], Float64.(sub.Observed[ok]); color=(:gray40, 0.8), markersize=ms)
            end
        end
        if err_in && !all(ismissing, sub.Input)
            ok = .!ismissing.(sub.Input)
            lines!(ax, sub.Qual[ok], Float64.(sub.Input[ok]); color=:black, linestyle=:dash)
        end
        if err_out && !all(ismissing, sub.Estimated)
            ok = .!ismissing.(sub.Estimated)
            lines!(ax, sub.Qual[ok], Float64.(sub.Estimated[ok]); color=:black)
        end
        if nominalQ && f != t
            lines!(ax, sub.Qual, sub.Nominal; color=:red)
        end
    end
    return fig
end

function DenoiseDNA.plotComplexity(files; kmerSize::Integer=2, window=nothing, by::Integer=5,
                                   n::Integer=100_000, bins::Integer=100, aggregate::Bool=false,
                                   size=nothing)
    fl = _files(files)
    cx = [let fq = read_fastq(f)
              seqComplexity([String(copy(s)) for s in fq.seq[1:min(n, end)]]; kmerSize=kmerSize, window=window, by=by)
          end for f in fl]
    groups = aggregate ? [("$(length(fl)) files (aggregated)", reduce(vcat, cx))] :
             [(basename(f), c) for (f, c) in zip(fl, cx)]
    nr, nc = _grid(length(groups))
    fig = Figure(; size=something(size, (380nc, 300nr)))
    for (k, (lab, c)) in enumerate(groups)
        r, col = fldmod1(k, nc)
        ax = Axis(fig[r, col]; title=lab, titlesize=12, xlabel="Effective Oligonucleotide Number", ylabel="Count",
                  xticks=0:(4^kmerSize ÷ 4):4^kmerSize)
        # ggplot's binning for `bins` over the scale limits: width range/(bins-1), centred on 0
        w = 4^kmerSize / (bins - 1)
        hist!(ax, filter(!isnan, c); bins=range(-w / 2, 4^kmerSize + w / 2; length=bins + 1), color=:gray30)
        xlims!(ax, -w / 2, 4^kmerSize + w / 2)
    end
    return fig
end

end # module
