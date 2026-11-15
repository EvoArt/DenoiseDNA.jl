@testset "FASTQ I/O" begin
    dir = mktempdir()
    reads = [("r1 extra", "ACGTN", [30, 31, 2, 40, 0]), ("r2", "acgt", [10, 20, 30, 40])]
    p = write_fastq(joinpath(dir, "x.fastq"), reads)
    fq = read_fastq(p)
    @test length(fq) == 2
    @test fq.id == ["r1 extra", "r2"]
    # (String(v) would take ownership of v and empty it, so compare bytes)
    @test fq.seq[1] == codeunits("ACGTN")
    @test fq.seq[2] == codeunits("ACGT")        # lowercase is upper-cased
    @test fq.qual[1] == UInt8[30, 31, 2, 40, 0]

    # gzip round trip through the writer
    gz = joinpath(dir, "y.fastq.gz")
    w = D.FastqWriter(gz)
    for k in 1:2
        D.write_record(w, fq.id[k], fq.seq[k], fq.qual[k])
    end
    close(w)
    fq2 = read_fastq(gz)
    @test fq2.id == fq.id && fq2.seq == fq.seq && fq2.qual == fq.qual

    # CRLF line endings
    crlf = joinpath(dir, "crlf.fastq")
    write(crlf, "@a\r\nACGT\r\n+\r\nIIII\r\n")
    @test read_fastq(crlf).seq[1] == codeunits("ACGT")

    bad = joinpath(dir, "bad.fastq")
    write(bad, "@a\nACGT\n+\nIII\n")
    @test_throws ErrorException read_fastq(bad)
end
