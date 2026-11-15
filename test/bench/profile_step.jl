# Profile one pipeline step on the MiSeq SOP data (run on Athena).
#
#   julia -t N --project=<benchenv> profile_step.jl <step> <r_out_dir> <out_dir> [nsamples]
#
# step ∈ filter, derep, learn, dada, dada_pool, merge, nochim
# Inputs are R's own outputs (filtered reads, error model), so each step can be
# profiled on its own. Writes <out_dir>/<step>_profile.txt (sorted by self and
# total time) and <out_dir>/<step>_time.txt.

using DenoiseDNA, Profile, ProfileToLLM, Printf
const D = DenoiseDNA
isdefined(DenoiseDNA, :setprogress!) && DenoiseDNA.setprogress!(false)

step, rdir, out = ARGS[1], ARGS[2], ARGS[3]
nsam = length(ARGS) >= 4 ? parse(Int, ARGS[4]) : 20
mkpath(out)
samples = readlines(joinpath(rdir, "samples.txt"))[1:nsam]
raw = joinpath(dirname(rdir), "..", "MiSeq_SOP")
fF = [joinpath(rdir, "filtered", s * "_F_filt.fastq.gz") for s in samples]
fR = [joinpath(rdir, "filtered", s * "_R_filt.fastq.gz") for s in samples]
rmat(p) = reduce(vcat, [permutedims(parse.(Float64, split(l, '\t'))) for l in readlines(p)])
errF = rmat(joinpath(rdir, "err_F.tsv")); errR = rmat(joinpath(rdir, "err_R.tsv"))
MT = Threads.nthreads() > 1 ? Threads.nthreads() : false

work = if step == "filter"
    rawF = [joinpath(raw, s * "_R1_001.fastq") for s in samples]
    rawR = [joinpath(raw, s * "_R2_001.fastq") for s in samples]
    tmp = mktempdir()
    () -> filterAndTrim(rawF, joinpath.(tmp, basename.(fF)), rawR, joinpath.(tmp, basename.(fR));
                        truncLen=(240, 160), maxEE=(2, 2), rm_phix=true, multithread=MT)
elseif step == "derep"
    () -> (derepFastq(fF), derepFastq(fR))
elseif step == "learn"
    dF = derepFastq(fF)
    () -> learnErrors(dF; multithread=MT)
elseif step == "dada"
    dF = derepFastq(fF)
    () -> dada(dF, errF; multithread=MT)
elseif step == "dada_pool"
    dF = derepFastq(fF)
    () -> dada(dF, errF; pool=true, multithread=MT)
elseif step in ("merge", "nochim")
    dF = derepFastq(fF); dR = derepFastq(fR)
    ddF = dada(dF, errF; multithread=MT); ddR = dada(dR, errR; multithread=MT)
    if step == "merge"
        () -> mergePairs(ddF, dF, ddR, dR)
    else
        st = makeSequenceTable(mergePairs(ddF, dF, ddR, dR))
        () -> removeBimeraDenovo(st; method="consensus", multithread=MT)
    end
else
    error("unknown step $step")
end

work()                                   # compile
t = minimum(@elapsed(work()) for _ in 1:2)
a = @allocated work()
open(joinpath(out, "$(step)_time.txt"), "w") do io
    @printf(io, "%s  threads=%d  samples=%d  min time %.3f s  alloc %.1f MiB\n", step, Threads.nthreads(), nsam, t, a / 2^20)
end
Profile.init(n=20_000_000, delay=0.0005)
Profile.clear()
Profile.@profile work()
open(joinpath(out, "$(step)_profile.txt"), "w") do io
    println(io, "== $step: sorted by self")
    print_profile(io, profile_table(Profile.fetch(); sortby=:self); max_rows=45)
    println(io, "\n== $step: sorted by total")
    print_profile(io, profile_table(Profile.fetch(); sortby=:total); max_rows=45)
    # inclusive time per DenoiseDNA function (summed over its lines)
    rows = profile_table(Profile.fetch(); sortby=:total)
    agg = Dict{Tuple{Symbol,String},Tuple{Float64,Float64}}()
    for r in rows
        f = String(r.file)
        occursin(r"src[/\\](dada[/\\])?[a-z_]+\.jl", f) || continue
        key = (r.func, basename(f))
        s, t = get(agg, key, (0.0, 0.0))
        agg[key] = (s + r.self_pct, max(t, r.total_pct))
    end
    println(io, "\n== $step: package functions (self% summed over lines, max total%)")
    for ((fn, file), (s, t)) in sort(collect(agg); by=x -> -x[2][2])
        t < 0.5 && continue
        println(io, rpad("$fn", 28), rpad(file, 16), lpad(round(s; digits=1), 6), lpad(round(t; digits=1), 7))
    end
end
println(read(joinpath(out, "$(step)_time.txt"), String))
