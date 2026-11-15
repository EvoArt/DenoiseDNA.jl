# Microbenchmark of the anti-diagonal aligner: per-alignment time of the
# previous scalar kernel vs the W-lane kernel for several W.
#   julia --project=<benchenv> align_micro.jl
using DenoiseDNA, BenchmarkTools, Random, Printf
using DenoiseDNA: _grow!, _usz, GAP, VecBuffer, AlignBuffer, load!, fill_dp!, traceback!, nwalign_vectorized!, encode_nt
include("old_aligner.jl")

rng = MersenneTwister(1)
mut(s) = (c = copy(s); for _ in 1:3; c[rand(rng, eachindex(c))] = rand(rng, 0x01:0x04); end;
          p = rand(rng, 20:length(c)-20); deleteat!(c, p); insert!(c, rand(rng, 20:length(c)-20), 0x02); c)
for L in (240, 160)
    ref = rand(rng, 0x01:0x04, L)
    qs = [mut(ref) for _ in 1:64]
    ob = OldBuf()
    t_old = @belapsed for q in $qs; old_nwalign_vectorized!($ob, $ref, q, 5, -4, -8, 0, 16); end
    nb = AlignBuffer()
    t_one = @belapsed for q in $qs; nwalign_vectorized!($nb, $ref, q, 5, -4, -8, 0, 16); end
    @printf("L=%d  old scalar %.2f us/aln   new W=1 %.2f us/aln\n", L, 1e6t_old / 64, 1e6t_one / 64)
    for W in (8, 16, 32, 64)
        vb = VecBuffer{W}(); tb = AlignBuffer()
        f = () -> for k in 1:W:64
            batch = view(qs, k:min(k + W - 1, 64))
            load!(vb, ref, batch); fill_dp!(vb, 5, -4, -8, 0, 16)
            for l in eachindex(batch); traceback!(tb, vb, l); end
        end
        tdp = @belapsed begin
            for k in 1:$W:64
                load!($vb, $ref, view($qs, k:min(k + $W - 1, 64))); fill_dp!($vb, 5, -4, -8, 0, 16)
            end
        end
        t = @belapsed $f()
        @printf("L=%d  W=%-2d  %.2f us/aln (DP only %.2f)\n", L, W, 1e6t / 64, 1e6tdp / 64)
    end
end
