@testset "filterAndTrim" begin
    dir = mktempdir()
    q(n, v=35) = fill(v, n)
    s = "ACGT"^15                                   # 60 bp
    reads = [
        ("keep", s, q(60)),
        ("lowq_tail", s, vcat(q(40), [2], q(19))),  # truncQ=2 truncates at 41 -> 40 bp
        ("hasN", "N" * s[2:end], q(60)),
        ("short", s[1:30], q(30)),
        ("highEE", s, q(60, 10)),                   # EE = 60 * 0.1 = 6
    ]
    fin = write_fastq(joinpath(dir, "in.fastq"), reads)

    out = joinpath(dir, "out.fastq.gz")
    r = filterAndTrim(fin, out; truncQ=2, minLen=35, maxN=0, maxEE=2, rm_phix=false)
    @test r.reads_in == [5]
    @test r.reads_out == [2]
    fq = read_fastq(out)
    @test fq.id == ["keep", "lowq_tail"]
    @test length(fq.seq[2]) == 40

    # truncLen: reads shorter than truncLen are dropped; others truncated
    r = filterAndTrim(fin, out; truncQ=2, truncLen=50, maxN=0, maxEE=Inf, rm_phix=false)
    fq = read_fastq(out)
    @test fq.id == ["keep", "highEE"]
    @test all(length.(fq.seq) .== 50)

    # trimLeft / trimRight
    r = filterAndTrim(fin, out; truncQ=0, trimLeft=5, trimRight=5, minLen=1, maxN=10, rm_phix=false)
    fq = read_fastq(out)
    @test length(fq.seq[1]) == 50 && String(fq.seq[1]) == s[6:55]

    # paired: a pair is kept only if both reads pass
    rev = write_fastq(joinpath(dir, "rev.fastq"), reverse(reads))
    fo = joinpath(dir, "f.fastq.gz"); ro = joinpath(dir, "r.fastq.gz")
    r = filterAndTrim(fin, fo, rev, ro; truncQ=2, minLen=35, maxEE=(2, Inf), rm_phix=false)
    # forward passes at positions {1, 2}; reverse (reversed list) passes at {1, 4, 5}
    fF = read_fastq(fo); fR = read_fastq(ro)
    @test r.reads_out == [1]
    @test fF.id == ["keep"] && fR.id == ["highEE"]

    # phiX removal
    phix = D._phix_genome()
    preads = [("phix", phix[1001:1150], q(150)), ("rcphix", rc(phix[2001:2150]), q(150)), ("other", randseq(MersenneTwister(1), 150), q(150))]
    pin = write_fastq(joinpath(dir, "p.fastq"), preads)
    @test isPhiX(String[x[2] for x in preads]) == [true, true, false]
    r = filterAndTrim(pin, out; rm_phix=true)
    @test read_fastq(out).id == ["other"]

    @test seqComplexity(["AAAAAAAAAA", "ACGTACGTACGTAGCTAGCATGCATCGATGC"])[1] ≈ 1.0
    @test seqComplexity(["ACGTACGTACGTAGCTAGCATGCATCGATGC"])[1] > 5
end
