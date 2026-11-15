# Why is dada(HOMOPOLYMER_GAP_PENALTY=...) slower than R? (run on Athena)
#   julia --project=<benchenv> homo_bench.jl <extdata dir> <errF.tsv> [out_dir]
# 1. per-alignment time of the scalar aligners on real reads (sam1F uniques)
# 2. warm dada time per alignment path
# 3. profile of the homopolymer dada run
using DenoiseDNA, BenchmarkTools, Printf, Profile, ProfileToLLM
using DenoiseDNA: AlignBuffer, nwalign_endsfree!, nwalign_endsfree_homo!, nwalign_vectorized!, encode_nt
DenoiseDNA.setprogress!(false)

ext, errfile = ARGS[1], ARGS[2]
out = length(ARGS) >= 3 ? ARGS[3] : nothing
err = reduce(vcat, [permutedims(parse.(Float64, split(l, '\t'))) for l in readlines(errfile)])
drp = derepFastq(joinpath.(ext, ["sam1F.fastq.gz", "sam2F.fastq.gz"]))

# 1. aligners
ref = encode_nt(drp[1].sequence[1])
qs = [encode_nt(s) for s in drp[1].sequence[2:129]]
buf = AlignBuffer()
per(t) = 1e6t / length(qs)
t_vec = @belapsed for q in $qs; nwalign_vectorized!($buf, $ref, q, 5, -4, -8, 0, 16); end
t_ef = @belapsed for q in $qs; nwalign_endsfree!($buf, $ref, q, 5, -4, -8, 16); end
t_homo = @belapsed for q in $qs; nwalign_endsfree_homo!($buf, $ref, q, 5, -4, -8, -1, 16); end
@printf("per alignment (L=%d, band 16): vectorized (scalar) %.2f us, endsfree %.2f us, endsfree_homo %.2f us\n",
        length(ref), per(t_vec), per(t_ef), per(t_homo))
a = @allocated nwalign_endsfree_homo!(buf, ref, qs[1], 5, -4, -8, -1, 16)
b = @allocated nwalign_endsfree!(buf, ref, qs[1], 5, -4, -8, 16)
println("allocated per call: homo $a bytes, endsfree $b bytes")

# 2. dada per path
for (nm, kw) in (("default", (;)), ("novec", (; VECTORIZED_ALIGNMENT=false)), ("homo", (; HOMOPOLYMER_GAP_PENALTY=-1)))
    dada(drp, err; kw...)
    t = @belapsed dada($drp, $err; $kw...) samples = 5 evals = 1
    @printf("dada %-8s %.3f s\n", nm, t)
end

# 3. profile
if out !== nothing
    mkpath(out)
    Profile.init(n=10_000_000, delay=0.0002)
    Profile.clear()
    Profile.@profile for _ in 1:5; dada(drp, err; HOMOPOLYMER_GAP_PENALTY=-1); end
    open(joinpath(out, "homo_profile.txt"), "w") do io
        print_profile(io, profile_table(Profile.fetch(); sortby=:self); max_rows=30)
    end
    println("profile written to ", joinpath(out, "homo_profile.txt"))
end
