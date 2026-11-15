@testset "chimeras" begin
    rng = MersenneTwister(11)
    p1 = randseq(rng, 200)
    p2 = mutate(rng, p1, 20)
    chim = p1[1:100] * p2[101:end]
    other = mutate(rng, p1, 25)
    @test isBimera(chim, [p1, p2])
    @test !isBimera(other, [p1, p2])
    @test !isBimera(p1, [p1, p2])                         # identical parent is ignored
    # one-off model: chimera with an extra mismatch
    chim1 = chim[1:150] * (chim[151] == 'A' ? "C" : "A") * chim[152:end]
    @test !isBimera(chim1, [p1, p2])
    @test isBimera(chim1, [p1, p2]; allowOneOff=true)

    unq = (sequence=[p1, p2, chim, other], abundance=[1000, 600, 40, 30])
    @test isBimeraDenovo(unq) == [false, false, true, false]
    @test isBimeraDenovo(unq; minParentAbundance=1000) == [false, false, false, false]

    st = SequenceTable(["s1", "s2", "s3"], [p1, p2, chim, other],
                       [1000 600 40 30; 800 500 30 0; 900 700 0 10])
    @test isBimeraDenovoTable(st) == [false, false, true, false]
    nc = removeBimeraDenovo(st)
    @test nc.sequences == [p1, p2, other]
    @test removeBimeraDenovo(st; method="pooled").sequences == [p1, p2, other]
    ps = removeBimeraDenovo(st; method="per-sample")
    @test ps.sequences == [p1, p2, other]
    @test removeBimeraDenovo(st; method="consensus", multithread=2).sequences == nc.sequences
    @test_throws ErrorException removeBimeraDenovo(st; method="nope")
end
