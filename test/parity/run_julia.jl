# Julia run of the MiSeq SOP pipeline, mirroring run_r.R.
#
#   julia -t N --project=<pkg> run_julia.jl <data_dir> <r_out_dir> <out_dir>
#
# Writes two result sets:
#   <out_dir>/e2e : Julia end-to-end from the raw reads
#   <out_dir>/iso : each step run on R's output from the previous step, so
#                   any difference is attributable to that step alone.

using DenoiseDNA, DataFrames, Printf
const D = DenoiseDNA

data_dir, rdir, out = ARGS[1], ARGS[2], ARGS[3]
const MT = Threads.nthreads() > 1 ? Threads.nthreads() : false
e2e = joinpath(out, "e2e"); iso = joinpath(out, "iso")
mkpath(joinpath(e2e, "filtered")); mkpath(iso)

g17(x) = x === missing || (x isa AbstractFloat && isnan(x)) ? "NA" : x isa AbstractFloat ? @sprintf("%.17g", x) : string(x)
function wtsv(df::DataFrame, path)
    open(path, "w") do io
        println(io, join(names(df), '\t'))
        for r in eachrow(df)
            println(io, join((g17(x) for x in r), '\t'))
        end
    end
end
function wmat(m, path)
    open(path, "w") do io
        for i in axes(m, 1)
            println(io, join((g17(x) for x in m[i, :]), '\t'))
        end
    end
end
rmat(path) = reduce(vcat, [permutedims([s == "NA" ? NaN : parse(Float64, s) for s in split(l, '\t')]) for l in readlines(path)])
rtsv(path) = begin
    L = readlines(path); h = split(L[1], '\t')
    rows = [split(l, '\t'; keepempty=true) for l in L[2:end]]
    DataFrame([Symbol(h[j]) => [r[j] for r in rows] for j in eachindex(h)])
end

timings = Pair{String,Float64}[]
macro timeit(label, ex)
    quote
        local t0 = time()
        local v = $(esc(ex))
        local dt = time() - t0
        push!(timings, $(esc(label)) => dt)
        @printf("[time] %-12s %8.2f s\n", $(esc(label)), dt)
        v
    end
end

samples = readlines(joinpath(rdir, "samples.txt"))
fnFs = [joinpath(data_dir, s * "_R1_001.fastq") for s in samples]
fnRs = [joinpath(data_dir, s * "_R2_001.fastq") for s in samples]
filtFs = [joinpath(e2e, "filtered", s * "_F_filt.fastq.gz") for s in samples]
filtRs = [joinpath(e2e, "filtered", s * "_R_filt.fastq.gz") for s in samples]
rfiltFs = [joinpath(rdir, "filtered", s * "_F_filt.fastq.gz") for s in samples]
rfiltRs = [joinpath(rdir, "filtered", s * "_R_filt.fastq.gz") for s in samples]

dump_derep(d, dir, tag) = begin
    wtsv(DataFrame(sequence=d.sequence, abundance=d.abundance), joinpath(dir, "derep_$tag.tsv"))
    wmat(d.quals, joinpath(dir, "derepq_$tag.tsv"))
    write(joinpath(dir, "derepmap_$tag.txt"), join(d.map, '\n') * "\n")
end
dump_err(e, dir, tag) = begin
    wmat(e.err_out, joinpath(dir, "err_$tag.tsv")); wmat(e.trans, joinpath(dir, "errtrans_$tag.tsv"))
    write(joinpath(dir, "errnconsist_$tag.txt"), string(length(e.err_in)))
end
dump_dada(dd, dir, tag) = begin
    wtsv(dd.clustering, joinpath(dir, "dada_$tag.tsv"))
    write(joinpath(dir, "dadamap_$tag.txt"), join(dd.map, '\n') * "\n")
    write(joinpath(dir, "dadapval_$tag.txt"), join(g17.(dd.pval), '\n') * "\n")
end
long(st) = DataFrame(st)

function run_downstream(ddF, drF, ddR, drR, dir; label="")
    mg = @timeit "merge$label" mergePairs(ddF, drF, ddR, drR; returnRejects=true)
    for (s, m) in zip(samples, mg)
        wtsv(m, joinpath(dir, "merge_$s.tsv"))
    end
    acc = [m[m.accept, :] for m in mg]
    st = @timeit "seqtab$label" makeSequenceTable(acc; names=samples)
    run_chimeras(st, dir; label=label)
end
function run_chimeras(st, dir; label="")
    wtsv(DataFrame(sequence=st.sequences), joinpath(dir, "seqtab_cols.tsv"))
    wtsv(long(st), joinpath(dir, "seqtab.tsv"))
    nc = @timeit "nochim$label" removeBimeraDenovo(st; method="consensus", multithread=MT)
    wtsv(DataFrame(sequence=nc.sequences), joinpath(dir, "nochim_consensus.tsv"))
    ncp = removeBimeraDenovo(st; method="pooled", multithread=MT)
    wtsv(DataFrame(sequence=ncp.sequences), joinpath(dir, "nochim_pooled.tsv"))
    ncs = removeBimeraDenovo(st; method="per-sample", multithread=MT)
    wtsv(long(ncs), joinpath(dir, "nochim_persample.tsv"))
end

# warm up (compilation) on one sample so timings measure run time
let d = derepFastq(rfiltFs[1])
    dada(d, D.getErrors(fill(0.01, 16, 41) .+ [i in (1, 6, 11, 16) ? 0.95 : 0.0 for i in 1:16, _ in 1:41]))
end

# ------------------------------------------------------------- end-to-end
ft = @timeit "filter" filterAndTrim(fnFs, filtFs, fnRs, filtRs; truncLen=(240, 160), maxN=0, maxEE=(2, 2),
                                    truncQ=2, rm_phix=true, compress=true, multithread=MT)
wtsv(DataFrame(sample=samples, reads_in=ft.reads_in, reads_out=ft.reads_out), joinpath(e2e, "filter.tsv"))
drF = @timeit "derep" derepFastq(filtFs); drR = derepFastq(filtRs)
for (s, a, b) in zip(samples, drF, drR)
    dump_derep(a, e2e, "F_$s"); dump_derep(b, e2e, "R_$s")
end
errF = @timeit "learn_F" learnErrors(filtFs; multithread=MT)
errR = @timeit "learn_R" learnErrors(filtRs; multithread=MT)
dump_err(errF, e2e, "F"); dump_err(errR, e2e, "R")
ddF = @timeit "dada_F" dada(drF, errF; multithread=MT)
ddR = @timeit "dada_R" dada(drR, errR; multithread=MT)
for (s, a, b) in zip(samples, ddF, ddR)
    dump_dada(a, e2e, "F_$s"); dump_dada(b, e2e, "R_$s")
end
ddFp = @timeit "dada_pool" dada(drF, errF; pool=true, multithread=MT)
ddFs = @timeit "dada_pseudo" dada(drF, errF; pool=:pseudo, multithread=MT)
for (s, a, b) in zip(samples, ddFp, ddFs)
    dump_dada(a, e2e, "Fpool_$s"); dump_dada(b, e2e, "Fpseudo_$s")
end
run_downstream(ddF, drF, ddR, drR, e2e)
wtsv(DataFrame(step=first.(timings), seconds=last.(timings)), joinpath(e2e, "timing.tsv"))
empty!(timings)

# ------------------------------------------------------------- isolated
idrF = @timeit "derep" derepFastq(rfiltFs); idrR = derepFastq(rfiltRs)
for (s, a, b) in zip(samples, idrF, idrR)
    dump_derep(a, iso, "F_$s"); dump_derep(b, iso, "R_$s")
end
ierrF = @timeit "learn_F" learnErrors(idrF; multithread=MT)
ierrR = @timeit "learn_R" learnErrors(idrR; multithread=MT)
dump_err(ierrF, iso, "F"); dump_err(ierrR, iso, "R")
rerrF = rmat(joinpath(rdir, "err_F.tsv")); rerrR = rmat(joinpath(rdir, "err_R.tsv"))
iddF = @timeit "dada_F" dada(idrF, rerrF; multithread=MT)
iddR = @timeit "dada_R" dada(idrR, rerrR; multithread=MT)
for (s, a, b) in zip(samples, iddF, iddR)
    dump_dada(a, iso, "F_$s"); dump_dada(b, iso, "R_$s")
end
iddFp = @timeit "dada_pool" dada(idrF, rerrF; pool=true, multithread=MT)
iddFs = @timeit "dada_pseudo" dada(idrF, rerrF; pool=:pseudo, multithread=MT)
for (s, a, b) in zip(samples, iddFp, iddFs)
    dump_dada(a, iso, "Fpool_$s"); dump_dada(b, iso, "Fpseudo_$s")
end

# merge on R's dada results
function r_dada(tag)
    cl = rtsv(joinpath(rdir, "dada_$tag.tsv"))
    seqs = String.(cl.sequence)
    clust = DataFrame(sequence=seqs, abundance=parse.(Int, cl.abundance), n0=parse.(Int, cl.n0))
    map = parse.(Int, readlines(joinpath(rdir, "dadamap_$tag.txt")))
    D.DadaResult(seqs, clust.abundance, clust, map, Float64[], zeros(Int, 16, 1), zeros(0, 0),
                 DataFrame(), nothing, nothing, D.DadaOptions())
end
rddF = [r_dada("F_$s") for s in samples]; rddR = [r_dada("R_$s") for s in samples]
mg = @timeit "merge" mergePairs(rddF, idrF, rddR, idrR; returnRejects=true)
for (s, m) in zip(samples, mg)
    wtsv(m, joinpath(iso, "merge_$s.tsv"))
end
# seqtab from R's merges, chimeras on R's seqtab
racc = [let m = rtsv(joinpath(rdir, "merge_$s.tsv"))
            m = m[m.accept .== "TRUE", :]
            DataFrame(sequence=String.(m.sequence), abundance=parse.(Int, m.abundance))
        end for s in samples]
ist = @timeit "seqtab" makeSequenceTable(racc; names=samples)
run_chimeras(ist, iso)
wtsv(DataFrame(step=first.(timings), seconds=last.(timings)), joinpath(iso, "timing.tsv"))
println("Julia pipeline complete")
