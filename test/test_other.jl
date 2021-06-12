#make sure the clamped spline is working
x = 2*π*sort(rand(100)) .- 3
y = sin.(x)
ϕ = CubicSplineInterpolator(x, y, cos(x[1]), cos(x[end]))
xx = collect(LinRange(minimum(x), maximum(x), 1000))
err = abs.(ϕ.(xx) - sin.(xx))
@test maximum(err) < 0.001

#test the cubic hermite
f(x) = 3x^3 - 2x^2 + 5x - 1
df(x) = 9x^2 - 4x + 5
xa = -2
xb = 3.1
xx = collect(LinRange(xa, xb, 100))
ff = f.(xx)
yy = cubichermite.(xx, xa, xb, f(xa), f(xb), df(xa), df(xb))
err = abs.(yy .- ff)
@test maximum(err) < 1e-6
