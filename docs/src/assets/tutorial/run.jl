# Run with: julia -t 4 --project=<environment> run.jl MiSeq_SOP output
# Environment: DenoiseDNA, DataFrames, CairoMakie. No taxonomy database is downloaded.
using DenoiseDNA, DataFrames, CairoMakie
using Printf

path, output = abspath.(ARGS[1:2])
mkpath(output)
fnFs = sort(filter(endswith("_R1_001.fastq"), readdir(path; join=true)))
fnRs = replace.(fnFs, "_R1_001.fastq" => "_R2_001.fastq")
@assert !isempty(fnFs) && all(isfile, fnRs)
sample_names = first.(split.(basename.(fnFs), "_"))
@assert allunique(sample_names)
filtFs = joinpath.(output, "filtered", sample_names .* "_F_filt.fastq.gz")
filtRs = joinpath.(output, "filtered", sample_names .* "_R_filt.fastq.gz")
mkpath(joinpath(output, "filtered"))
save(joinpath(output, "quality-forward.png"), plotQualityProfile(fnFs[1:2]))
save(joinpath(output, "quality-reverse.png"), plotQualityProfile(fnRs[1:2]))
out = filterAndTrim(fnFs, filtFs, fnRs, filtRs; truncLen=(240,160),
    maxN=0, maxEE=(2,2), truncQ=2, rm_phix=true, compress=true, multithread=true)
errF = learnErrors(filtFs; multithread=true)
errR = learnErrors(filtRs; multithread=true)
save(joinpath(output, "errors-forward.png"), plotErrors(errF; nominalQ=true))
save(joinpath(output, "errors-reverse.png"), plotErrors(errR; nominalQ=true))
dadaFs = dada(filtFs, errF; multithread=true)
dadaRs = dada(filtRs, errR; multithread=true)
mergers = mergePairs(dadaFs, filtFs, dadaRs, filtRs)
seqtab = makeSequenceTable(mergers; names=sample_names)
seqtab_nochim = removeBimeraDenovo(seqtab; method="consensus", multithread=true)
track = DataFrame(sample=sample_names, input=out.reads_in, filtered=out.reads_out,
    denoisedF=[sum(d.abundance) for d in dadaFs],
    denoisedR=[sum(d.abundance) for d in dadaRs],
    merged=[sum(m.abundance) for m in mergers],
    nonchim=vec(sum(seqtab_nochim.table; dims=2)))
function tsv(file, df)
    open(joinpath(output, file), "w") do io
        println(io, join(names(df), '\t'))
        for r in eachrow(df)
            println(io, join(r, '\t'))
        end
    end
end
tsv("read-tracking.tsv", track)
tsv("asv-counts.tsv", DataFrame(seqtab_nochim))
lens = length.(seqtab.sequences)
lengths = sort(unique(lens))
fig = Figure(size=(800,380))
ax = Axis(fig[1,1]; xlabel="Merged sequence length (bp)", ylabel="ASVs", xticks=lengths)
barplot!(ax, lengths, [count(==(l), lens) for l in lengths]; color=:steelblue)
save(joinpath(output, "asv-lengths.png"), fig)
stages = [:input, :filtered, :denoisedF, :denoisedR, :merged, :nonchim]
fig = Figure(size=(1000,480))
ax = Axis(fig[1,1]; xticks=(1:length(stages), string.(stages)),
    ylabel="Reads remaining (% of input)", title="Each line is one sample")
for r in eachrow(track)
    lines!(ax, 1:length(stages), 100 .* [r[s] / r.input for s in stages]; color=(:steelblue,0.45))
end
ylims!(ax, 0, 105)
save(joinpath(output, "read-retention.png"), fig)
# The mock is a control, not a time point. Sort biological samples by numeric day.
keep = findall(s -> startswith(s, "F3D"), sample_names)
days = parse.(Int, replace.(sample_names[keep], "F3D" => ""))
order = sortperm(days); keep = keep[order]; days = days[order]
counts = seqtab_nochim.table[keep,:]
depth = vec(sum(counts; dims=2))
richness = vec(sum(counts .> 0; dims=2))
shannon = [begin p = row[row .> 0] ./ sum(row); -sum(p .* log.(p)) end for row in eachrow(counts)]
fig = Figure(size=(1000,380))
for (j, (values, label)) in enumerate(((depth,"Reads"),(richness,"Observed ASVs"),(shannon,"Shannon diversity")))
    ax = Axis(fig[1,j]; xlabel="Day", ylabel=label)
    scatterlines!(ax, days, values; color=:steelblue)
end
save(joinpath(output, "diversity.png"), fig)
tsv("diversity.tsv", DataFrame(sample=sample_names[keep], day=days, reads=depth, observed=richness, shannon=shannon))
i = findfirst(==("Mock"), sample_names)
mock = seqtab_nochim.sequences[seqtab_nochim.table[i,:] .> 0]
mock_ref = getSequences(joinpath(path,"HMP_MOCK.v35.fasta"))
matched = [any(r -> occursin(s,r) || occursin(rc(s),r), mock_ref) for s in mock]
open(joinpath(output,"summary.txt"),"w") do io
    println(io, "Julia: ", VERSION, "; DenoiseDNA: ", pkgversion(DenoiseDNA))
    println(io, "Before chimera removal: ", size(seqtab))
    println(io, "After chimera removal: ", size(seqtab_nochim))
    println(io, "Merged reads retained: ", sum(seqtab_nochim.table)/sum(seqtab.table))
    println(io, "Mock ASVs: ", length(mock), "; exact reference matches: ", sum(matched))
    show(io, MIME("text/plain"), track)
end
println("Tutorial complete: ", output)
