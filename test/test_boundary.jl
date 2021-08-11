T1 = [LinearInterpolator,
      CubicInterpolator,
      CubicSplineInterpolator]
T2 = [BilinearInterpolator,
      BicubicInterpolator,
      BicubicSplineInterpolator]

for T in T1
    ϕ = T(x->1.0, -1, 1, 5)
    @test_throws ErrorException ϕ(-2)
    @test_throws ErrorException ϕ(2)
    ϕ = T(x->1.0, -1, 1, 5, NoBoundaries())
    @test isreal(ϕ(-2))
    @test isreal(ϕ(2))
    ϕ = T(x->1.0, -1, 1, 5, WeakBoundaries())
    @test isreal(ϕ(-1 - 1e-10))
    @test isreal(ϕ(1 + 1e-10))
end

for T in T2
    Φ = T((x,y)->1.0, -1, 1, 5, -1, 1, 5)
    @test_throws ErrorException Φ(-2, 0)
    @test_throws ErrorException Φ(0, 2)
    @test_throws ErrorException Φ(2, 0)
    @test_throws ErrorException Φ(0, -2)
    Φ = T((x,y)->1.0, -1, 1, 5, -1, 1, 5, NoBoundaries())
    @test isreal(Φ(-2, 0))
    @test isreal(Φ(0, 2))
    @test isreal(Φ(2, 0))
    @test isreal(Φ(0, -2))
    Φ = T((x,y)->1.0, -1, 1, 5, -1, 1, 5, WeakBoundaries())
    @test isreal(Φ(-1 - 1e-10, 0))
    @test isreal(Φ(0, 1 + 1e-10))
    @test isreal(Φ(1 + 1e-10, 0))
    @test isreal(Φ(0, -1 - 1e-10))
end
