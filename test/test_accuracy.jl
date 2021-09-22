function interpolator_error(T::Type, f::Function, xa, xb, n::Int)
    ϕ = T(f, xa, xb, n)
    x = collect(LinRange(xa, xb, 100))
    fy = f.(x)
    ϕy = ϕ.(x)
    maximum(abs.(ϕy .- fy))
end

function interpolator_error(T::Type, f::Function, xa, xb, nx::Int, ya, yb, ny::Int)
    Φ = T(f, xa, xb, nx, ya, yb, ny)
    x = collect(LinRange(xa, xb, 50))
    y = collect(LinRange(ya, yb, 50))
    X = x .* ones(length(y))'
    Y = y' .* ones(length(x))
    maximum(@. abs(Φ(X, Y) - f(X, Y)) )
end

T1 = [LinearInterpolator,
      CubicInterpolator,
      CubicSplineInterpolator,
      ChebyshevInterpolator]
T2 = [BilinearInterpolator,
      BicubicInterpolator,
      BichebyshevInterpolator]

for T ∈ T1
    @test interpolator_error(T, x->sin(x^2), -π, π, 128) < 1e-2
end

for T ∈ T2
    @test interpolator_error(T, (x,y)->sin(x)+cos(y), -π, π, 128, -π, π, 128) < 1e-2
end
