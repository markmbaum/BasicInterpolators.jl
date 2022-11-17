using Test
using BasicInterpolators, Random

@testset "base" begin
    
    @testset "findcell" begin
        V = 0:3 |> collect
        n = 4
        @test findcell(-1, V, n) == 1 #out of left boundary
        @test findcell(4, V, n) == 3 #out of right boundary
        @test findcell(1.5, V, n) == 2
    end

    @testset "boundaries" begin

        @test NoBoundaries()(nothing) === nothing

        B = WeakBoundaries()
        @test B(0, -1, 1) === nothing
        @test B(-1 - 1e-10, -1, 1) === nothing
        @test B( 1 + 1e-10, -1, 1) === nothing
        @test_throws ErrorException B(-2, -1, 1)
        @test_throws ErrorException B( 2, -1, 1)

        B = StrictBoundaries()
        @test B(0, -1, 1) === nothing
        @test_throws ErrorException B(-1 - 1e-10, -1, 1)
        @test_throws ErrorException B( 1 + 1e-10, -1, 1)

    end

end

@testset "chebyshev" begin
    
    @test ischebygrid(chebygrid(-2, 2, 9))

    p = ChebyshevInterpolator(sin, -3, 3, 32)
    d = chebyderiv(p)
    x = LinRange(-3, 3, 128)
    @test all(d.(x) .≈ cos.(x))

end


@testset "polynomial" begin
    
    @test linear(0.25, [0, 1], [0, 1]) ≈ 0.25
    @test quadratic(0.5, [-1, 0, 1], [1, 0, 1]) ≈ 0.25
    @test cubic(1.5, [0,1,2,3], [1,0,-1,2]) ≈ -0.75
    c = [4, -1, 2, -5, 1]
    x = 0:4 |> collect
    y = [4, 1, -14, -35, -32]
    @test neville(2.5, x, y) ≈ -25.0625
    @test all(vandermonde(x, y) .≈ c)
    @test cubichermite(0, -1, 1, 1, -1, 1, -1) ≈ 0.5

end

@testset "1D" begin
    
    x = LinRange(-3, 3, 101)
    f(x) = sin(x^2)
    for ℐ ∈ [LinearInterpolator,
             CubicInterpolator,
             CubicSplineInterpolator,
             ChebyshevInterpolator]
        n = ℐ <: OneDimensionalInterpolator ? 256 : 32
        p = ℐ(f, -3, 3, n)
        y = p.(x)
        @test maximum(@. abs(y - f(x))) < 1e-2
    end
    p = CubicSplineInterpolator(x, sin.(x), cos(-3), cos(3)) #clamped spline
    y = p.(x)
    @test maximum(@. abs(y - sin(x)) < 1e-2)

end

@testset "2D" begin
    
    n = 51
    x = LinRange(-3, 3, n)
    y = LinRange(-3, 3, n)
    X = x' .* ones(n)
    Y = ones(n)' .* y
    f(x,y) = sin(x)*cos(y) + x*y/4
    for ℐ ∈ [BilinearInterpolator,
             BicubicInterpolator,
             BichebyshevInterpolator]
        n = ℐ <: TwoDimensionalInterpolator ? 128 : 32
        P = ℐ(f, -3, 3, n, -3, 3, n)
        Z = P.(X, Y)
        @test maximum(@. abs(Z - f(X, Y))) < 1e-2
        @test maximum(@. abs(Z - f(X .|> Float32, Y))) < 1e-2
    end

end

@testset "scattered" begin
    
    f(x) = @inbounds sin(x[1])*cos(x[2]) + x[3]
    rng = MersenneTwister(609)

    for (n, 𝒯, tol) ∈ ([10_000, ShepardInterpolator, 0.8], [1_000, RBFInterpolator, 1e-2])
        X = 2*rand(rng, n, 3) .- 1
        y = map(f, eachslice(X, dims=1))
        X_test = 2*rand(rng, n ÷ 4, 3) .- 1
        y_test = map(f, eachslice(X_test, dims=1))
        P = 𝒯(X, y)
        y = map(P, eachslice(X_test, dims=1))
        @test maximum(@. abs(y - y_test)) < tol
    end

    #n = 10_000
    #X = 2*rand(n, 3) .- 1
    #y = map(f, eachslice(X, dims=1))
    #X_test = 2*rand(n, 3) .- 1
    #y_test = map(f, eachslice(X_test, dims=1))
    #P = ShepardInterpolator(X, y)
    #y = map(P, eachslice(X_test, dims=1))
    #@test maximum(@. abs(y - y_test)) < 0.8
#
    #n = 1_000
    #X = 2*rand(n, 3) .- 1
    #y = map(f, eachslice(X, dims=1))
    #X_test = 2*rand(n, 3) .- 1
    #y_test = map(f, eachslice(X_test, dims=1))
    #P = RBFInterpolator(X, y, 1.0)
    #y = map(P, eachslice(X_test, dims=1))
    #@test maximum(@. abs(y - y_test)) < 1e-2

end