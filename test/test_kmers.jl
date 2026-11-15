@testset "k-mers" begin
    rng = MersenneTwister(5)
    seqs = [D.encode_nt(randseq(rng, rand(rng, 150:250))) for _ in 1:40]
    t = D.KmerTable(seqs)
    @test size(t.k16, 2) == 0                       # 8-bit counts suffice
    for _ in 1:200
        i, j = rand(rng, 1:40, 2)
        k1 = D.assign_kmer!(zeros(UInt16, D.N_KMER), seqs[i]); k2 = D.assign_kmer!(zeros(UInt16, D.N_KMER), seqs[j])
        @test D.kmer_dist(t, i, length(seqs[i]), j, length(seqs[j])) ==
              D.kmer_dist(k1, length(seqs[i]), k2, length(seqs[j]))
    end
    # a count above 254 switches to exact 16-bit counts
    long = [D.encode_nt("A"^300), D.encode_nt("A"^280 * "C"^20)]
    tl = D.KmerTable(long)
    @test size(tl.k16, 2) == 2
    @test D.kmer_dist(tl, 1, 300, 2, 300) ≈ 1 - (276 + 0) / 296
    @test D.kord_dist(D.assign_kmer_order!(zeros(UInt16, 300), long[1]), 300,
                      D.assign_kmer_order!(zeros(UInt16, 300), long[2]), 300) ≈ 1 - 276 / 296
end
