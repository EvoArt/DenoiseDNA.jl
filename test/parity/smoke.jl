# Quick end-to-end smoke run of the Julia pipeline on a few samples.
#   julia --project=<pkg> smoke.jl <data_dir> <out_dir> [nsamples]
using DenoiseDNA
data_dir, out = ARGS[1], ARGS[2]
n = length(ARGS) >= 3 ? parse(Int, ARGS[3]) : 2
mkpath(out)

function step(f, label)
    println(">> $label  (live heap $(round(Base.gc_live_bytes() / 2^20, digits=1)) MiB)"); flush(stdout)
    t = @elapsed v = f()
    println("<< $label  $(round(t, digits=2)) s  (live heap $(round(Base.gc_live_bytes() / 2^20, digits=1)) MiB)"); flush(stdout)
    v
end

fs = sort(filter(f -> endswith(f, "_R1_001.fastq"), readdir(data_dir; join=true)))[1:n]
rs = replace.(fs, "_R1_001" => "_R2_001")
ff = joinpath.(out, replace.(basename.(fs), ".fastq" => ".filt.fastq.gz"))
fr = joinpath.(out, replace.(basename.(rs), ".fastq" => ".filt.fastq.gz"))
fq = step(() -> read_fastq(fs[1]), "read_fastq")
println(length(fq), " reads"); flush(stdout)
ft = step(() -> filterAndTrim(fs, ff, rs, fr; truncLen=(240, 160), maxEE=(2, 2), rm_phix=false), "filter (no phix)")
println(ft)
step(() -> isPhiX(String.(fq.seq[1:200])), "isPhiX x200")
dF = step(() -> derepFastq(ff), "derep F"); dR = step(() -> derepFastq(fr), "derep R")
println(dF[1]); flush(stdout)
toy = [i in (1, 6, 11, 16) ? 0.97 : 0.01 for i in 1:16, _ in 1:41]
r1 = step(() -> dada(dF[1], toy), "dada one sample (toy err)")
println(r1); flush(stdout)
eF = step(() -> learnErrors(ff; verbose=true), "learnErrors F")
eR = step(() -> learnErrors(fr), "learnErrors R")
println(eF); flush(stdout)
ddF = step(() -> dada(dF, eF), "dada F"); ddR = step(() -> dada(dR, eR), "dada R")
println.(ddF); flush(stdout)
step(() -> dada(dF, eF; pool=true), "dada pool")
step(() -> dada(dF, eF; pool=:pseudo), "dada pseudo")
mg = step(() -> mergePairs(ddF, dF, ddR, dR), "mergePairs")
st = step(() -> makeSequenceTable(mg), "seqtab")
println(st)
for m in ("consensus", "pooled", "per-sample")
    nc = step(() -> removeBimeraDenovo(st; method=m), "chimeras $m")
    println(m, ": ", nc)
end
println("SMOKE OK")
