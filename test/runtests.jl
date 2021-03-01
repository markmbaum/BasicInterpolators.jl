using Test, BasicInterpolators

#arrays of the interpolator types for running tests
T1 = [LinearInterpolator,
      CubicInterpolator,
      CubicSplineInterpolator,
      ChebyshevInterpolator]
T2 = [BilinearInterpolator,
      BicubicInterpolator,
      BicubicSplineInterpolator,
      BichebyshevInterpolator]

@testset "Accuracy" begin include("test_accuracy.jl") end
@testset "Boundary" begin include("test_boundary.jl") end
@testset "Nevilles" begin include("test_nevilles.jl") end
