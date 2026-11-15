@testset "plot data" begin
    dir = mktempdir()
    reads = [("a", "ACGTACGTAC", fill(30, 10)), ("b", "ACGTACGTAC", vcat(fill(20, 5), fill(40, 5))),
             ("c", "ACGTAC", fill(10, 6))]
    p = qualityProfile(write_fastq(joinpath(dir, "q.fastq"), reads))
    @test p.nreads == 3 && p.sampled == 3
    @test size(p.counts) == (10, 41)
    @test vec(sum(p.counts; dims=2)) == [3, 3, 3, 3, 3, 3, 2, 2, 2, 2]
    @test p.mean[1] ≈ 20.0 && p.mean[10] ≈ 35.0
    @test p.q50[1] == 20 && p.q25[1] == 10 && p.q75[1] == 30
    @test p.cum[1] == 1.0 && p.cum[10] ≈ 2 / 3
    @test qualityProfile(joinpath(dir, "q.fastq"); n=1).sampled == 1

    trans = zeros(Int, 16, 41)
    trans[1, 31] = 990; trans[2, 31] = 10
    err = toy_err()
    e = errorProfile(ErrorModel(err, [err], trans))
    @test nrow(e) == 16 * 41
    r = e[(e.Transition .== "A2C") .& (e.Qual .== 30), :]
    @test r.Observed[1] ≈ 0.01 && r.Estimated[1] == err[2, 31]
    @test ismissing(e[(e.Transition .== "A2A") .& (e.Qual .== 30), :Observed][1])
    @test errorProfile(err).Observed |> x -> all(ismissing, x)
    @test_throws MethodError plotErrors(err)      # needs a Makie backend
end
