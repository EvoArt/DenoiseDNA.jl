@testset "utilities" begin
    rng = MersenneTwister(7)

    @testset "getUniques collapse / getSequences / getAbund" begin
        nt = (sequence=["CCA", "AAC", "CCA"], abundance=[1, 2, 3])
        @test getUniques(nt) == (sequence=["AAC", "CCA"], abundance=[2, 4])
        @test getUniques(nt; collapse=false) == nt
        @test getSequences(nt) == ["CCA", "AAC", "CCA"]
        @test getSequences(nt; collapse=true) == ["AAC", "CCA"]
        @test getSequences(["acg", "ACG"]; collapse=true) == ["ACG", "ACG"]   # R: unique, then toupper
        @test D.getAbund(nt) == 6
    end

    @testset "uniquesToFasta / read_fasta" begin
        seqs = [randseq(rng, 50) for _ in 1:3]
        nt = (sequence=seqs, abundance=[30, 20, 10])
        mktempdir() do dir
            fa = joinpath(dir, "u.fa")
            uniquesToFasta(nt, fa; width=20)
            L = readlines(fa)
            @test L[1] == ">sq1;size=30;"
            @test length(L[2]) == 20 && length(L[4]) == 10 && length(L) == 12
            @test getSequences(fa) == seqs
            @test D.read_fasta(fa).id == ["sq1;size=30;", "sq2;size=20;", "sq3;size=10;"]
            uniquesToFasta(nt, fa; mode="a")
            @test length(getSequences(fa)) == 6
        end
    end

    @testset "setDadaOpt / getDadaOpt" begin
        saved = getDadaOpt()
        try
            setDadaOpt(OMEGA_A=1e-20, MATCH=2)
            @test getDadaOpt(:OMEGA_A) == 1e-20
            @test D._options().OMEGA_A == 1e-20
            @test D._options(OMEGA_A=1e-5).OMEGA_A == 1e-5      # keywords take precedence
            @test D.nwalign("ACGTACGT", "ACGAACGT") == D.nwalign("ACGTACGT", "ACGAACGT"; match=2)
            @test_throws ArgumentError setDadaOpt(NOT_AN_OPTION=1)
        finally
            setDadaOpt(saved)
        end
        @test getDadaOpt(:OMEGA_A) == 1e-40 && getDadaOpt(:MATCH) == 5
    end

    @testset "nwextract / strdiff" begin
        s = randseq(rng, 120)
        @test D.nwextract(s[31:80], s) == s[31:80]
        m = collect(s); m[10] = m[10] == 'A' ? 'C' : 'A'; m = String(m)
        d = D.strdiff(s, m)
        @test d.pos == [10] && d.nt0 == [string(s[10])] && d.nt1 == [string(m[10])]
    end

    @testset "shifts" begin
        s = randseq(rng, 200)
        @test isShiftedPair(s[1:180], s[21:200])
        @test !isShiftedPair(s[1:180], s[21:200]; minOverlap=170)
        @test !isShiftedPair(s, s[21:180])                  # subsequence
        @test isShiftedPair(s, s[21:180]; flagSubseqs=true)
        other = randseq(rng, 200)
        nt = (sequence=[s[1:190], other, s[5:200], other[1:150]], abundance=[100, 90, 5, 95])
        @test isShiftDenovo(nt) == [false, false, true, false]  # other[1:150] is a subsequence of other
    end

    @testset "mergeSequenceTables" begin
        a = SequenceTable(["s1", "s2"], ["AAAA", "CCCC"], [1 2; 3 4])
        b = SequenceTable(["s3"], ["CCCC", "GGGG"], [5 6])
        m = mergeSequenceTables(a, b)
        @test m.samples == ["s1", "s2", "s3"]
        @test m.sequences == ["CCCC", "GGGG", "AAAA"]
        @test m.table == [2 0 1; 4 0 3; 5 6 0]
        @test_throws ErrorException mergeSequenceTables(a, a)
        s = mergeSequenceTables(a, a; repeats="sum", orderBy=nothing)
        @test s.samples == ["s1", "s2"] && s.table == 2 .* a.table
        a2 = SequenceTable(["s1"], ["AAAC", "CCAG"], [1 2])
        c = SequenceTable(["s4"], ["GGGA", "GTTT"], [7 1])          # GTTT = rc(AAAC)
        r = mergeSequenceTables([a2, c]; tryRC=true, orderBy=nothing)
        @test r.sequences == ["AAAC", "CCAG", "GGGA"]
        @test r.table == [1 2 0; 1 0 7]
    end

    @testset "getBadBases / isBadBaseFP" begin
        clust = DataFrame(sequence=[repeat("A", 10) for _ in 1:12], birth_ham=[missing; fill(1, 10); 2])
        # nine one-off clusters born at position 4, one at position 7
        bs = DataFrame(pos=[fill(4, 9); 7; 4; 5], clust=[2:11; 12; 12])
        @test getBadBases(clust, bs; omegaB=1e-3) == [4]
        @test getBadBases(clust, bs; omegaB=1e-3, minOccurence=10) == Int[]
        fp = isBadBaseFP(clust, bs; omegaB=1e-3)
        @test fp == [false; fill(true, 9); false; false]      # cluster 12: 1 of 2 subs bad < 0.51
        @test isBadBaseFP(clust, bs; omegaB=1e-3, minFraction=0.5)[12]
    end
end

@testset "primer matching and removal" begin
    m(p, s, mm=0; indels=false, fixed=D.isACGT(p)) = D._match_pattern(D._dnabits(p), D._dnabits(s), mm, indels, fixed)
    # expected values are R/Biostrings 2.78 outputs (test/parity/probe_match.R)
    s = "AGTTGCAAAGTTGCATTTT"
    @test m("TTGCA", s, 1) == [3:7, 11:15]
    @test m("CAGTT", s, 1) == [0:4, 8:12, 14:18]
    @test m("GGAGTT", s, 2) == [-1:4, 7:12, 13:18]
    @test m("TTTTGG", s, 2) == [16:21]
    @test m("ACGTA", "AAACGTNAAA"; fixed=true) == UnitRange{Int}[]
    @test m("ACGTA", "AAACGTNAAA"; fixed=false) == [3:7]
    @test m("RCGT", "AAAMCGTAAA") == [4:7]
    @test m("AAA", "AAAAAA") == [1:3, 2:4, 3:5, 4:6]
    s4 = "CCCCACGTTGCACCCC"
    @test m("ACGTGCA", s4, 1; indels=true) == [5:12]
    @test m("ACGTTTGCA", s4, 1; indels=true) == [5:12]
    @test m("ACGTTGCA", s4, 2; indels=true) == [5:12]
    @test m("GCACCCCAA", s4, 2; indels=true) == [10:16]
    @test m("CACG", "CCACGACGCC", 1; indels=true) == [2:5, 6:8]

    rng = MersenneTwister(11)
    fwd = "AGRGTTYGATYMTGGCTCAG"; rev = rc("RGYTACCTTGTTACGACTT")
    pf = "AGAGTTTGATCCTGGCTCAG"; pr = rc("AGCTACCTTGTTACGACTT")     # concrete instances
    inserts = [randseq(rng, 300) for _ in 1:4]
    reads = [pf * inserts[1] * pr, rc(pf * inserts[2] * pr), randseq(rng, 340), pf * inserts[4]]
    mktempdir() do dir
        fin = write_fastq(joinpath(dir, "in.fastq"), [("r$k", r, fill(30, length(r))) for (k, r) in enumerate(reads)])
        fout = joinpath(dir, "out.fastq.gz")
        t = removePrimers(fin, fout, fwd, rev)
        @test t.reads_in == [4] && t.reads_out == [2]
        @test getSequences(fout) == inserts[1:2]
        t = removePrimers(fin, fout, fwd, rev; orient=false)
        @test t.reads_out == [1]
        t = removePrimers(fin, fout, fwd)
        @test t.reads_out == [3]
        @test getSequences(fout)[3] == inserts[4]
    end
end

@testset "filterAndTrim keeping no reads" begin
    rng = MersenneTwister(2)
    reads = [("r$k", randseq(rng, 100), fill(30, 100)) for k in 1:5]
    mktempdir() do dir
        fF = write_fastq(joinpath(dir, "F.fastq"), reads); fR = write_fastq(joinpath(dir, "R.fastq"), reads)
        t = filterAndTrim(fF, joinpath(dir, "oF.fastq.gz"), fR, joinpath(dir, "oR.fastq.gz"); maxLen=50)
        @test t.reads_out == [0]
        t = filterAndTrim(fF, joinpath(dir, "o.fastq.gz"); maxLen=50, rm_lowcomplex=4)
        @test t.reads_out == [0]
    end
end
