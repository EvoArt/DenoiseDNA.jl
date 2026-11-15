@testset "dada" begin
    rng = MersenneTwister(2024)
    L = 150
    s1 = randseq(rng, L)
    s2 = mutate(rng, s1, 6)
    s3 = mutate(rng, s1, 15)
    truth = [s1, s2, s3]
    dir = mktempdir()
    reads = simulate_reads(rng, truth, [1500, 400, 120]; len=L)
    p = write_fastq(joinpath(dir, "sim.fastq"), reads)
    d = derepFastq(p)
    err = toy_err()

    r = dada(d, err)
    @test r isa DadaResult
    @test Set(r.sequence) == Set(truth)
    # abundance of each ASV close to the simulated number of reads
    ab = Dict(zip(r.sequence, r.abundance))
    @test abs(ab[s1] - 1500) < 40 && abs(ab[s2] - 400) < 20 && abs(ab[s3] - 120) < 10
    @test length(r.map) == length(d)
    @test all(0 .<= r.map .<= length(r.sequence))
    @test sum(r.abundance) == sum(d.abundance[r.map .> 0])
    @test size(r.trans) == (16, 41)
    @test sum(r.trans) == sum(d.abundance[r.map .> 0] .* L)
    @test nrow(r.clustering) == length(r.sequence)
    @test ismissing(r.clustering.birth_from[1])

    # options are passed through; a huge OMEGA_A gives at least as many clusters
    r2 = dada(d, err; OMEGA_A=1e-5)
    @test length(r2.sequence) >= length(r.sequence)
    @test_throws ArgumentError dada(d, err; NOT_AN_OPTION=1)
    # non-vectorized aligner gives the same partition here
    r3 = dada(d, err; VECTORIZED_ALIGNMENT=false)
    @test Set(r3.sequence) == Set(r.sequence)

    # priors can rescue a rare variant
    s4 = mutate(rng, s1, 2)
    reads4 = vcat(reads, simulate_reads(rng, [s4], [3]; len=L, qmean=39))
    d4 = derepFastq(write_fastq(joinpath(dir, "sim4.fastq"), reads4))
    @test (s4 in dada(d4, err; priors=[s4]).sequence) || !(s4 in d4.sequence)

    # multiple samples, pooling and multithreading give consistent results
    rs = dada([d, d], err)
    @test length(rs) == 2 && rs[1].sequence == r.sequence
    rp = dada([d, d], err; pool=true)
    @test Set(rp[1].sequence) ⊇ Set(r.sequence)        # pooling has more power
    @test sum(rp[1].abundance) <= sum(d.abundance)
    @test rp[1].sequence == rp[2].sequence
    rps = dada([d, d], err; pool=:pseudo)
    @test Set(rps[1].sequence) ⊇ Set(r.sequence)
    rt = dada([d, d], err; multithread=2)
    @test rt[2].sequence == r.sequence && rt[2].map == r.map

    # learnErrors recovers sensible rates from the simulated data
    e = learnErrors(d)
    @test e isa ErrorModel
    @test size(e.err_out, 1) == 16
    @test 1 <= length(e.err_in) <= 10
    @test e.err_out[2, 36] < e.err_out[2, 16]
    rl = dada(d, e)
    @test Set(rl.sequence) == Set(truth)
    # reusing alignments across rounds changes nothing
    e0 = learnErrors(d; cache_alignments=false)
    @test e0.err_out == e.err_out && e0.trans == e.trans && length(e0.err_in) == length(e.err_in)
    @test dada([d, d], err; pool=:pseudo, cache_alignments=false)[1].map == rps[1].map
    # one thread vs several: identical partitions
    @test dada(d, err; multithread=4).map == r.map
end
