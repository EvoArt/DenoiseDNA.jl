@testset "derepFastq" begin
    dir = mktempdir()
    reads = [("1", "CCCC", [30, 30, 30, 30]), ("2", "AAAA", [30, 31, 2, 40]),
             ("3", "CCCC", [31, 30, 30, 30]), ("4", "GGGG", [10, 10, 10, 10]),
             ("5", "AAAA", [31, 31, 3, 40]), ("6", "TTTT", [20, 20, 20, 20])]
    p = write_fastq(joinpath(dir, "d.fastq"), reads)
    d = derepFastq(p)
    # most abundant first; ties keep lexical order (A < C < G < T)
    @test d.sequence == ["AAAA", "CCCC", "GGGG", "TTTT"]
    @test d.abundance == [2, 2, 1, 1]
    @test d.map == [2, 1, 2, 3, 1, 4]
    # mean (not floor) quality
    @test d.quals[1, :] == [30.5, 31.0, 2.5, 40.0]
    @test d.quals[2, 1] == 30.5

    # variable length -> NaN padding
    p2 = write_fastq(joinpath(dir, "v.fastq"), [("a", "ACG", [1, 2, 3]), ("b", "ACGTT", [1, 2, 3, 4, 5])])
    d2 = derepFastq(p2)
    @test size(d2.quals) == (2, 5)
    @test isnan(d2.quals[findfirst(==("ACG"), d2.sequence), 4])

    # chunking preserves R semantics: uniques first seen in a later chunk come after
    d3 = derepFastq(p; n=2)
    @test d3.abundance == [2, 2, 1, 1]
    @test sum(d3.abundance) == length(reads)

    # combineDereps pools abundances and averages qualities by abundance
    c = combineDereps([d, d])
    @test c.abundance == 2 .* d.abundance
    @test c.quals ≈ d.quals
    @test length(c.map) == 2length(d.map)
end
