# Render the diagnostic plots for the MiSeq SOP run (visual check).
#   julia --project=<env with CairoMakie> plots.jl <r_out_dir> <out_dir>
using DenoiseDNA, CairoMakie
rdir, out = ARGS[1], ARGS[2]
mkpath(out)
samples = readlines(joinpath(rdir, "samples.txt"))
fF = [joinpath(rdir, "filtered", s * "_F_filt.fastq.gz") for s in samples]
raw = joinpath(dirname(rdir), "..", "MiSeq_SOP")
rawF = [joinpath(raw, s * "_R1_001.fastq") for s in samples]
save(joinpath(out, "quality_raw_F.png"), plotQualityProfile(rawF[1:4]))
save(joinpath(out, "quality_filt_agg.png"), plotQualityProfile(fF; aggregate=true))
e = learnErrors(fF[1:3])
save(joinpath(out, "errors_F.png"), plotErrors(e; nominalQ=true))
save(joinpath(out, "complexity.png"), plotComplexity(rawF[1:2]))
save(joinpath(out, "errors_F.pdf"), plotErrors(e; nominalQ=true))
println("plots written to $out")
