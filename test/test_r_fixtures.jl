# Regression tests against outputs produced by R's dada2 (see test/parity/make_fixtures.R).

@testset "R reference fixtures" begin
    @testset "vectorized aligner == dada2 C_nwvec" begin
        pairs = [split(l, '\t') for l in readlines(joinpath(FIXTURES, "align_pairs.tsv"))]
        ral = [split(l, '\t') for l in readlines(joinpath(FIXTURES, "align_r.tsv"))]
        buf = D.AlignBuffer()
        params = Dict("16" => (5, -4, -8, 0, 16), "nofree" => (1, -64, -64, -64, -1))
        n = length(pairs); nbad = 0
        for (k, r) in enumerate(ral)
            a, b = pairs[mod1(k, n)]
            D.nwalign_vectorized!(buf, Vector{UInt8}(codeunits(a)), Vector{UInt8}(codeunits(b)), params[r[1]]...)
            nbad += (String(buf.al0[1:buf.len]) != r[2] || String(buf.al1[1:buf.len]) != r[3])
        end
        @test nbad == 0
    end

    @testset "loess == R stats::loess" begin
        num(s) = s == "NA" ? NaN : parse(Float64, s)
        worst = 0.0
        for l in readlines(joinpath(FIXTURES, "loess_r.tsv"))
            f = split(l, '\t')
            q = num.(split(f[2], ",")); y = num.(split(f[3], ",")); w = num.(split(f[4], ","))
            r = num.(split(f[5], ","))
            p = D.predict(D.loess(q, y; weights=w), 0:40)
            @test isnan.(p) == isnan.(r)
            ok = .!isnan.(r)
            length(q) >= 12 && (worst = max(worst, maximum(abs.(p[ok] .- r[ok]) ./ abs.(r[ok]); init=0.0)))
        end
        @test worst < 1e-12
    end

    mock = joinpath(FIXTURES, "mock_F.fastq.gz")
    if isfile(mock)
        @testset "dada on the MiSeq SOP mock community == R" begin
            num(s) = s == "NA" ? NaN : parse(Float64, s)
            err = reduce(vcat, [permutedims(num.(split(l, '\t'))) for l in readlines(joinpath(FIXTURES, "mock_F_err.tsv"))])
            d = derepFastq(mock)
            rd = [split(l, '\t') for l in readlines(joinpath(FIXTURES, "mock_F_derep.tsv"))[2:end]]
            @test d.sequence == first.(rd)
            @test d.abundance == parse.(Int, last.(rd))
            r = dada(d, err)
            L = readlines(joinpath(FIXTURES, "mock_F_dada.tsv")); h = split(L[1], '\t')
            rows = [split(l, '\t'; keepempty=true) for l in L[2:end]]
            col(c) = [x[findfirst(==(c), h)] for x in rows]
            @test r.sequence == col("sequence")
            @test r.abundance == parse.(Int, col("abundance"))
            @test r.clustering.n0 == parse.(Int, col("n0"))
            @test r.map == parse.(Int, readlines(joinpath(FIXTURES, "mock_F_map.txt")))
            rp = num.(readlines(joinpath(FIXTURES, "mock_F_pval.txt")))
            @test maximum(abs.(rp .- r.pval) ./ max.(rp, 1e-300)) < 1e-10
        end
    end
end
