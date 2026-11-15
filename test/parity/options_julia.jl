# Julia half of the options parity check (see options_r.R). Each function
# runs on R's inputs from the previous step (isolated), so a difference is
# pinned to one function. Exits 1 on any mismatch.
#   julia --project=. test/parity/options_julia.jl <outdir of options_r.R> <MiSeq_SOP dir>
using DenoiseDNA, DataFrames, CodecZlib
const D = DenoiseDNA
DenoiseDNA.setprogress!(false)

rdir, sop = ARGS
f(name) = joinpath(rdir, name)
jdir = mkpath(joinpath(rdir, "julia"))
ext = readchomp(f("extdata.txt"))
nfail = 0
function check(name, ok, detail="")
    global nfail += !ok
    println(ok ? "PASS " : "FAIL ", name, isempty(detail) ? "" : "  -- " * detail)
end
readmat(path) = reduce(vcat, [permutedims(parse.(Float64, split(l, '\t'))) for l in readlines(path)])
relmax(a, b) = size(a) == size(b) ? maximum(abs.(a .- b) ./ max.(abs.(b), 1e-300); init=0.0) : Inf
gunzip(path) = isfile(path) ? String(transcode(GzipDecompressor, read(path))) : nothing
timing = Dict{String,Float64}()
# warm timing: run once to compile, then time a second run
macro timed_as(name, ex)
    quote
        $(esc(ex))
        t0 = time(); v = $(esc(ex)); timing[$(esc(name))] = time() - t0; v
    end
end

# ---------------------------------------------------------------- filterAndTrim
rawF = joinpath(sop, "F3D0_S188_L001_R1_001.fastq"); rawR = joinpath(sop, "F3D0_S188_L001_R2_001.fastq")
shuf = f("R2_shuffled.fastq.gz")
filt = Dict(
    "trim" => (rawR, (; trimLeft=(10, 5), trimRight=(5, 10), maxEE=(2, 3), minLen=50)),
    "lens" => (rawR, (; maxLen=(245, 240), minLen=(200, 150), minQ=(10, 5), truncQ=(11, 2), truncLen=(0, 0))),
    "minq" => (rawR, (; minLen=(200, 150), minQ=(10, 5), truncQ=(11, 2))),
    "orient" => (rawR, (; orient_fwd="TACGG", truncLen=(240, 160), maxEE=2)),
    "matchids" => (shuf, (; matchIDs=true, truncLen=(240, 160), maxEE=2)),
    "idfield" => (shuf, (; matchIDs=true, id_field=1, truncLen=(240, 160))),
    "maxn" => (rawR, (; maxN=1, truncQ=0, rm_phix=false)),
    "single" => (nothing, (; truncLen=200, maxEE=1)),
    "lowcomp" => (nothing, (; rm_lowcomplex=8, truncLen=150)))
for l in readlines(f("filt_counts.txt"))
    nm, rin, rout = split(l, '\t')
    rev, kw = filt[nm]
    oF = joinpath(jdir, "$(nm)_F.fastq.gz"); oR = joinpath(jdir, "$(nm)_R.fastq.gz")
    t = rev === nothing ? filterAndTrim(rawF, oF; kw...) : filterAndTrim(rawF, oF, rev, oR; kw...)
    same = gunzip(oF) == gunzip(f("filt/$(nm)_F.fastq.gz")) &&
           (rev === nothing || gunzip(oR) == gunzip(f("filt/$(nm)_R.fastq.gz")))
    check("filterAndTrim $nm", (string(t.reads_in[1]), string(t.reads_out[1])) == (rin, rout) && same,
          "J $(t.reads_in[1])/$(t.reads_out[1]) R $rin/$rout, files identical: $same")
end

if isfile(f("complexity_k2.txt"))
    seqs = [String(copy(s)) for s in D.read_fastq(rawF).seq]
    for (file, kw) in (("complexity_k2.txt", (;)), ("complexity_k3w50.txt", (; kmerSize=3, window=50, by=7)))
        r = parse.(Float64, readlines(f(file))); j = seqComplexity(seqs; kw...)
        check("seqComplexity $file", relmax(j, r) < 1e-12, "max rel diff $(relmax(j, r))")
    end
end

# ---------------------------------------------------------------- dada options
fnF = joinpath.(ext, ["sam1F.fastq.gz", "sam2F.fastq.gz"]); fnR = joinpath.(ext, ["sam1R.fastq.gz", "sam2R.fastq.gz"])
errF = readmat(f("errF.tsv")); errR = readmat(f("errR.tsv"))
drpF = derepFastq(fnF); drpR = derepFastq(fnR)
priors = readlines(f("priors.txt"))
function same_dada(dd, nm)
    bad = String[]
    for (i, d) in enumerate(dd)
        L = readlines(f("dada_$(nm)_$i.tsv"))[2:end]
        rs = [split(x, '\t')[1] for x in L]; ra = [parse(Int, split(x, '\t')[2]) for x in L]
        rmap = [x == "NA" ? 0 : parse(Int, x) for x in readlines(f("dadamap_$(nm)_$i.txt"))]
        d.sequence == rs || push!(bad, "sample $i ASVs ($(length(d.sequence)) vs $(length(rs)))")
        d.abundance == ra || push!(bad, "sample $i abundances")
        d.map == rmap || push!(bad, "sample $i map")
    end
    return isempty(bad), join(bad, "; ")
end
dopts = Dict("default" => (;), "omega" => (; OMEGA_A=1e-20, OMEGA_C=1e-10), "singletons" => (; DETECT_SINGLETONS=true),
             "priors" => (; priors=priors, OMEGA_P=1e-2, OMEGA_A=1e-60), "homo" => (; HOMOPOLYMER_GAP_PENALTY=-1),
             "novec" => (; VECTORIZED_ALIGNMENT=false), "band32" => (; BAND_SIZE=32), "band0" => (; BAND_SIZE=0),
             "unbanded" => (; BAND_SIZE=-1), "minfold" => (; MIN_FOLD=2, MIN_HAMMING=2, MIN_ABUNDANCE=3),
             "maxclust" => (; MAX_CLUST=5), "nokmers" => (; USE_KMERS=false), "kdist" => (; KDIST_CUTOFF=0.2),
             "nogapless" => (; GAPLESS=false), "nogreedy" => (; GREEDY=false),
             "gaps" => (; GAP_PENALTY=-4, MATCH=4, MISMATCH=-5))
for nm in sort(collect(keys(dopts)))
    dd = @timed_as "dada_$nm" dada(drpF, errF; dopts[nm]...)
    ok, why = same_dada(dd, nm)
    check("dada $nm", ok, why)
end
errNQ = readmat(f("errNQ.tsv"))
let e = @timed_as "learn_noqual" learnErrors(fnF; USE_QUALS=false)
    check("learnErrors USE_QUALS=false", relmax(getErrors(e), errNQ) < 1e-12, "max rel diff $(relmax(getErrors(e), errNQ))")
end
let (ok, why) = same_dada(dada(drpF, errNQ; USE_QUALS=false), "noqual")
    check("dada USE_QUALS=false", ok, why)
end

# ---------------------------------------------------------------- error functions
trans = Int.(readmat(f("transF.tsv")))
for (nm, fun) in (("loess", loessErrfun), ("noqual", noqualErrfun), ("binned", makeBinnedQualErrfun([2, 11, 25, 37, 40])),
                  ("inflate", t -> inflateErr(errF, 3)), ("inflate_self", t -> inflateErr(errF, 2; inflateSelfTransitions=true)))
    r = readmat(f("errfun_$nm.tsv")); j = fun(trans)
    check("errfun $nm", relmax(j, r) < 1e-12, "max rel diff $(relmax(j, r))")
end

# ---------------------------------------------------------------- PacBio
drpPB = derepFastq(joinpath(ext, "samPB.fastq.gz"))
transPB = Int.(readmat(f("transPB.tsv"))); errPB = readmat(f("errPB.tsv"))
let r = readmat(f("errfun_pacbio.tsv")), j = PacBioErrfun(transPB)
    check("errfun pacbio", relmax(j, r) < 1e-12, "max rel diff $(relmax(j, r))")
end
let e = @timed_as "learn_pb" learnErrors(drpPB; errorEstimationFunction=PacBioErrfun, BAND_SIZE=32)
    check("learnErrors PacBio (end-to-end)", e.trans == transPB && relmax(getErrors(e), errPB) < 1e-10,
          "trans identical: $(e.trans == transPB), err max rel diff $(relmax(getErrors(e), errPB))")
end
let dd = @timed_as "dada_pb" dada(drpPB, errPB; BAND_SIZE=32)
    ok, why = same_dada([dd], "pacbio")
    check("dada PacBio variable-length", ok, why)
end

# ---------------------------------------------------------------- mergePairs
base = dada(drpF, errF); ddR = dada(drpR, errR)
cell(x::Bool) = x ? "TRUE" : "FALSE"
cell(::Missing) = "NA"
cell(x) = string(x)
function same_table(df::DataFrame, path)
    L = readlines(path)
    hdr = split(L[1], '\t')
    names(df) == hdr || return false, "columns J $(names(df)) R $hdr"
    rows = [join(cell.(collect(r)), '\t') for r in eachrow(df)]
    rows == L[2:end] || return false, "rows differ ($(count(rows .!= L[2:end])) of $(length(rows)))"
    return true, ""
end
mopts = Dict("mismatch" => (; maxMismatch=2, minOverlap=20, returnRejects=true), "concat" => (; justConcatenate=true),
             "overhang" => (; trimOverhang=true, returnRejects=true), "propcol" => (; propagateCol=["n0", "birth_ham", "nunq"]))
for nm in sort(collect(keys(mopts)))
    mm = mergePairs(base, drpF, ddR, drpR; mopts[nm]...)
    res = [same_table(m, f("merge_$(nm)_$i.tsv")) for (i, m) in enumerate(mm)]
    check("mergePairs $nm", all(first, res), join(last.(res), "; "))
end

# ---------------------------------------------------------------- chimeras
rres = joinpath(dirname(rdir), "results_nt1", "r")
cols = readlines(joinpath(rres, "seqtab_cols.tsv"))[2:end]
lg = [split(l, '\t') for l in readlines(joinpath(rres, "seqtab.tsv"))[2:end]]
samples = unique(first.(lg))
tab = zeros(Int, length(samples), length(cols))
srow = Dict(s => i for (i, s) in enumerate(samples)); scol = Dict(s => j for (j, s) in enumerate(cols))
for (s, q, a) in lg
    tab[srow[s], scol[q]] = parse(Int, a)
end
st = SequenceTable(samples, cols, tab)
rl(name) = parse.(Int, readlines(f(name))) .== 1
check("isBimeraDenovo allowOneOff", isBimeraDenovo(st; allowOneOff=true) == rl("chim_oneoff.txt"))
check("isBimeraDenovo allowOneOff, distance 2, maxShift 8",
      isBimeraDenovo(st; allowOneOff=true, minOneOffParentDistance=2, maxShift=8) == rl("chim_oneoff2.txt"))
check("isBimeraDenovoTable minSampleFraction/ignoreNNegatives",
      isBimeraDenovoTable(st; minSampleFraction=0.5, ignoreNNegatives=2) == rl("chim_table.txt"))
check("isBimeraDenovoTable allowOneOff", isBimeraDenovoTable(st; allowOneOff=true, minFoldParentOverAbundance=1) ==
                                         rl("chim_table_oneoff.txt"))
let kept = removeBimeraDenovo(st; method="consensus", minSampleFraction=0.6, ignoreNNegatives=0)
    check("removeBimeraDenovo consensus options", [s in Set(kept.sequences) for s in st.sequences] == rl("chim_remove_kept.txt"))
end

# ---------------------------------------------------------------- timing
rt = Dict(String(split(l, '\t')[1]) => parse(Float64, split(l, '\t')[2]) for l in readlines(f("timing.tsv"))[2:end])
println("\nstep\tR (s)\tJulia (s)")
for k in sort(collect(keys(rt)))
    haskey(timing, k) && println(k, '\t', round(rt[k]; digits=2), '\t', round(timing[k]; digits=2))
end
println(nfail == 0 ? "ALL PASS" : "$nfail FAILED")
exit(nfail == 0 ? 0 : 1)
