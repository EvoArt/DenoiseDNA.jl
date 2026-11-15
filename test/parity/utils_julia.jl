# Julia half of the utilities parity check (see utils_r.R). Compares each
# output with R's and prints one line per check; exits 1 on any mismatch.
#   julia --project=. test/parity/utils_julia.jl <outdir of utils_r.R>
using DenoiseDNA, DataFrames
const D = DenoiseDNA
DenoiseDNA.setprogress!(false)

rdir = ARGS[1]
ext = readchomp(joinpath(rdir, "extdata.txt"))
f(name) = joinpath(rdir, name)
jdir = mkpath(joinpath(rdir, "julia"))
fnF = joinpath.(ext, ["sam1F.fastq.gz", "sam2F.fastq.gz"])
fnR = joinpath.(ext, ["sam1R.fastq.gz", "sam2R.fastq.gz"])

nfail = 0
function check(name, ok, detail="")
    global nfail += !ok
    println(ok ? "PASS " : "FAIL ", name, ok || isempty(detail) ? "" : "  -- " * detail)
end

function readtsv(path)
    L = readlines(path)
    h = split(L[1], '\t')
    rows = [split(l, '\t'; keepempty=true) for l in L[2:end]]
    return DataFrame([h[k] => [r[k] for r in rows] for k in eachindex(h)])
end
function readst(path)
    L = readlines(path)
    seqs = String.(split(L[1], '\t')[2:end])
    rows = [split(l, '\t') for l in L[2:end]]
    tab = isempty(rows) ? zeros(Int, 0, length(seqs)) :
          reduce(vcat, [permutedims(parse.(Int, r[2:end])) for r in rows])
    return SequenceTable(String.(first.(rows)), seqs, tab)
end
same(a::SequenceTable, b::SequenceTable) =
    a.samples == b.samples && a.sequences == b.sequences && a.table == b.table
relmax(a, b) = maximum(abs.(a .- b) ./ max.(abs.(b), floatmin()); init=0.0)

# getSequences
check("getSequences(fasta)", getSequences(joinpath(ext, "example_seqs.fa")) == readlines(f("getseq_fa.txt")))
check("getSequences(fastq)", getSequences(fnF[1]) == readlines(f("getseq_fq.txt")))

# sameOrder
check("sameOrder", string.(uppercase.(string.([D.sameOrder(fnF[1], fnR[1]), D.sameOrder(fnF[1], fnR[2])]))) ==
                   readlines(f("sameorder.txt")))

# checkConvergence (via our learnErrors; loess differs from R at ~1e-14)
for (fn, tag) in ((fnF, "F"), (fnR, "R"))
    err = learnErrors(fn)
    c = D.checkConvergence(err); r = parse.(Float64, readlines(f("checkconv_$tag.txt")))
    check("checkConvergence($tag)", length(c) == length(r) && (r[end] == 0) == (c[end] == 0) &&
          relmax(c[1:end-1], r[1:end-1]) < 1e-10, "J=$c R=$r")
end

# uniquesToFasta
for i in 1:2
    u = readtsv(f("uniques_$i.tsv"))
    nt = (sequence=String.(u.sequence), abundance=parse.(Int, u.abundance))
    for (suffix, kw) in (("", (;)), ("_w60", (; width=60)), ("_ids", (; ids=["ASV$k" for k in eachindex(nt.sequence)])))
        out = joinpath(jdir, "unq_$i$suffix.fa")
        uniquesToFasta(nt, out; kw...)
        check("uniquesToFasta $i$suffix", read(out) == read(f("unq_$i$suffix.fa")))
    end
end
let out = joinpath(jdir, "unq_append.fa")
    nts = [readtsv(f("uniques_$i.tsv")) for i in 1:2]
    uniquesToFasta((sequence=String.(nts[1].sequence), abundance=parse.(Int, nts[1].abundance)), out)
    uniquesToFasta((sequence=String.(nts[2].sequence), abundance=parse.(Int, nts[2].abundance)), out; mode="a")
    check("uniquesToFasta mode=a", read(out) == read(f("unq_append.fa")))
end

# getBadBases / isBadBaseFP
let nbad = 0, nlines = 0, first_bad = ""
    cl = [readtsv(f("clust_$i.tsv")) for i in 1:2]
    bs = [readtsv(f("bsubs_$i.tsv")) for i in 1:2]
    for d in cl; d.birth_ham = [x == "NA" ? missing : parse(Int, x) for x in d.birth_ham]; d.sequence = String.(d.sequence); end
    for d in bs; d.pos = parse.(Int, d.pos); d.clust = parse.(Int, d.clust); end
    for l in readlines(f("badbases.txt"))
        p = split(l, '\t'; keepempty=true)
        i = parse(Int, p[1]); omegaB = parse(Float64, p[2]); minOcc = parse(Int, p[3])
        b = getBadBases(cl[i], bs[i]; omegaB, minOccurence=minOcc)
        fp = isBadBaseFP(cl[i], bs[i]; omegaB, minOccurence=minOcc, minFraction=0.3)
        ok = join(b, ",") == p[4] && join(Int.(fp), ",") == p[5]
        nlines += 1
        if !ok
            nbad += 1
            isempty(first_bad) && (first_bad = "$(p[1:3]): J bb=$(join(b, ",")) R bb=$(p[4])")
        end
    end
    check("getBadBases/isBadBaseFP ($nlines settings; nonempty in R: $(count(l -> !isempty(split(l, '\t'; keepempty=true)[4]), readlines(f("badbases.txt")))))",
          nbad == 0, first_bad)
end

# mergeSequenceTables, collapseNoMismatch
st = readst(f("st.tsv")); stF = readst(f("stF.tsv")); stRC = readst(f("stRC.tsv")); stDup = readst(f("stDup.tsv"))
check("mergeSequenceTables plain", same(mergeSequenceTables(st, stF), readst(f("merge_plain.tsv"))))
check("mergeSequenceTables tryRC", same(mergeSequenceTables(st, stF, stRC; tryRC=true), readst(f("merge_rc.tsv"))))
check("mergeSequenceTables repeats=sum", same(mergeSequenceTables(st, stDup; repeats="sum", orderBy="nsamples"),
                                              readst(f("merge_sum.tsv"))))
check("mergeSequenceTables orderBy=nothing", same(mergeSequenceTables([stF, st]; orderBy=nothing),
                                                  readst(f("merge_noorder.tsv"))))
check("collapseNoMismatch", same(collapseNoMismatch(st), readst(f("collapse.tsv"))))

# isShiftDenovo
let u = readtsv(f("shift_input.tsv"))
    nt = (sequence=String.(u.sequence), abundance=parse.(Int, u.abundance))
    for (file, kw) in (("shift_default.txt", (;)), ("shift_200.txt", (; minOverlap=200)))
        r = parse.(Int, readlines(f(file)))
        j = Int.(isShiftDenovo(nt; kw...))
        check("isShiftDenovo $file ($(sum(r)) flagged in R)", j == r, "J flags $(findall(==(1), j)), R $(findall(==(1), r))")
    end
end

# nwextract, strdiff
let next = 0, ndiff = 0
    for l in readlines(f("nwextract_strdiff.txt"))
        s, q, m, ext_r, pos, nt0, nt1 = String.(split(l, '\t'; keepempty=true))
        next += D.nwextract(q, m) != ext_r
        d = D.strdiff(s, m)
        ndiff += (join(d.pos, ","), join(d.nt0, ","), join(d.nt1, ",")) != (pos, nt0, nt1)
    end
    check("nwextract", next == 0)
    check("strdiff", ndiff == 0)
end

println(nfail == 0 ? "ALL PASS" : "$nfail FAILED")
exit(nfail == 0 ? 0 : 1)
