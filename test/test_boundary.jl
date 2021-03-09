T1 = [LinearInterpolator,
      CubicInterpolator,
      CubicSplineInterpolator,
      ChebyshevInterpolator]
T2 = [BilinearInterpolator,
      BicubicInterpolator,
      BicubicSplineInterpolator,
      BichebyshevInterpolator]

for T in T1
    ϕ = T(x->1, -1, 1, 5)
    @test_throws AssertionError ϕ(-2)
    @test_throws AssertionError ϕ(2)
end

for T in T2
    Φ = T((x,y)->1, -1, 1, 5, -1, 1, 5)
    @test_throws AssertionError Φ(-2, 0)
    @test_throws AssertionError Φ(0, 2)
    @test_throws AssertionError Φ(2, 0)
    @test_throws AssertionError Φ(0, -2)
end
