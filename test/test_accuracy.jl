function interpolator_error(T::Type, f::Function, xa::Real, xb::Real, n::Int)
    ϕ = T(f, xa, xb, n)
    x = collect(LinRange(xa, xb, 1000))
    fy = f.(x)
    ϕy = ϕ.(x)
    maximum(abs.(ϕy .- fy))
end

function interpolator_error(T::Type, f::Function,
               xa::Real, xb::Real, nx::Int,
               ya::Real, yb::Real, ny::Int)
    Φ = T(f, xa, xb, nx, ya, yb, ny)
    x = collect(LinRange(xa, xb, 250))
    y = collect(LinRange(ya, yb, 250))
    X = x .* ones(length(y))'
    Y = y' .* ones(length(x))
    fZ = f.(X, Y)
    ΦZ = Φ.(X, Y)
    maximum(abs.(ΦZ .- fZ))
end

T1 = [LinearInterpolator,
      CubicInterpolator,
      CubicSplineInterpolator,
      ChebyshevInterpolator]
T2 = [BilinearInterpolator,
      BicubicInterpolator,
      BicubicSplineInterpolator,
      BichebyshevInterpolator]

for T ∈ T1
    @test interpolator_error(T, x->sin(x^2), -3, 3, 128) < 1e-2
end

for T ∈ T2
    @test interpolator_error(T, (x,y)->sin(x)+cos(y), -3, 3, 128, -3, 3, 128) < 1e-2
end
