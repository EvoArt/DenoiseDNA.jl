# Julia half of the removePrimers parity check (see primers_r.R).
#   julia --project=. test/parity/primers_julia.jl <outdir of primers_r.R>
using DenoiseDNA, CodecZlib
const D = DenoiseDNA

rdir = ARGS[1]
f(name) = joinpath(rdir, name)
jdir = mkpath(joinpath(rdir, "julia"))
nfail = 0
function check(name, ok, detail="")
    global nfail += !ok
    println(ok ? "PASS " : "FAIL ", name, ok || isempty(detail) ? "" : "  -- " * detail)
end

# 1. matcher
let nbad = 0, shown = 0, n = 0
    for l in readlines(f("match_cases.txt"))
        pat, subj, mm, indels, fixed, expect = split(l, '\t'; keepempty=true)
        m = D._match_pattern(D._dnabits(pat), D._dnabits(subj), parse(Int, mm), indels == "1", fixed == "1")
        got = join(("$(first(r)):$(last(r))" for r in m), ",")
        n += 1
        if got != expect
            nbad += 1
            shown < 10 && println("  case $pat $subj mm=$mm indels=$indels fixed=$fixed: J=$got R=$expect")
            shown += 1
        end
    end
    check("Biostrings matchPattern/vmatchPattern ($n cases)", nbad == 0, "$nbad differ")
end

# 2. removePrimers
fn = readchomp(f("input.txt"))
F27 = "AGRGTTYGATYMTGGCTCAG"; R1492 = "RGYTACCTTGTTACGACTT"
runs = Dict("default" => (;), "noorient" => (; orient=false), "indels" => (; allow_indels=true),
            "notrim" => (; trim_fwd=false, trim_rev=false), "fwdonly" => (;), "mm0" => (; max_mismatch=0),
            "mm5" => (; max_mismatch=5, allow_indels=true))
gunzip(path) = isfile(path) ? String(transcode(GzipDecompressor, read(path))) : nothing
for l in readlines(f("pb_counts.txt"))
    nm, rin, rout = split(l, '\t')
    out = joinpath(jdir, "pb_$nm.fastq.gz")
    rev = nm == "fwdonly" ? nothing : rc(R1492)
    t = removePrimers(fn, out, F27, rev; runs[nm]...)
    j = gunzip(out); r = gunzip(f("pb_$nm.fastq.gz"))
    detail = "J in/out $(t.reads_in[1])/$(t.reads_out[1]), R $rin/$rout"
    check("removePrimers $nm", string(t.reads_in[1]) == rin && string(t.reads_out[1]) == rout && j == r, detail)
end

println(nfail == 0 ? "ALL PASS" : "$nfail FAILED")
exit(nfail == 0 ? 0 : 1)
