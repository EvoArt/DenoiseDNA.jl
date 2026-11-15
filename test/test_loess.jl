@testset "loess" begin
    # A local quadratic fit reproduces a quadratic exactly, and cubic Hermite
    # interpolation of exact values and slopes is exact for quadratics.
    x = collect(0.0:40.0)
    y = 0.5 .- 0.1 .* x .+ 0.002 .* x .^ 2
    fit = D.loess(x, y; weights=rand(MersenneTwister(1), 41) .* 100)
    @test D.predict(fit, x) ≈ y atol = 1e-10
    @test D.predict(fit, [10.5, 33.25]) ≈ 0.5 .- 0.1 .* [10.5, 33.25] .+ 0.002 .* [10.5, 33.25] .^ 2 atol = 1e-10
    # no extrapolation outside the data
    fit = D.loess(x[5:30], y[5:30])
    p = D.predict(fit, x)
    @test all(isnan, p[1:4]) && all(isnan, p[31:end]) && !any(isnan, p[5:30])
    # direct surface evaluates the local fit itself
    fd = D.loess(x, y; surface=:direct)
    @test D.predict(fd, [12.3])[1] ≈ 0.5 - 0.1 * 12.3 + 0.002 * 12.3^2 atol = 1e-10

    # errors: rates decrease with quality and rows sum to one
    trans = zeros(Int, 16, 41)
    rng = MersenneTwister(3)
    for q in 0:40, i in 0:3
        n = 10_000
        e = 10.0^(-q / 10)
        for j in 0:3
            trans[4i+j+1, q+1] = i == j ? round(Int, n * (1 - e)) : round(Int, n * e / 3 * (0.8 + 0.4rand(rng)))
        end
    end
    err = loessErrfun(trans)
    @test size(err) == (16, 41)
    @test all(0 .<= err .<= 1)
    for i in 0:3
        @test vec(sum(err[4i+1:4i+4, :]; dims=1)) ≈ ones(41)
    end
    @test err[2, 11] > err[2, 41]
    @test all(err[[2, 3, 4], :] .>= 1e-7)

    e2 = noqualErrfun(trans)
    @test all(e2[:, 1] .== e2[:, end])
    bq = makeBinnedQualErrfun([0, 11, 25, 40])
    e3 = bq(trans)
    @test size(e3) == (16, 41)
    @test inflateErr(err, 2)[2, 41] > err[2, 41]
end
