#make sure the clamped spline is working
x = 2*π*sort(rand(100)) .- 3
y = sin.(x)
ϕ = CubicSplineInterpolator(x, y, cos(x[1]), cos(x[end]))
xx = LinRange(minimum(x), maximum(x), 1000)
err = abs.(ϕ.(xx) - sin.(xx))
@test maximum(err) < 0.001

#test cubic hermite
f(x) = 3x^3 - 2x^2 + 5x - 1
df(x) = 9x^2 - 4x + 5
xa = -2
xb = 3.1
xx = LinRange(xa, xb, 100)
yy = cubichermite.(xx, xa, xb, f(xa), f(xb), df(xa), df(xb))
@test maximum(@. abs(yy - f.(xx)) ) < 1e-6

#test the type-flexible Bichebyshev functor
P = BichebyshevInterpolator((x,y)->sin(x)*cos(x), -π, π, 14, -π, π, 16)
for i in 1:10
    a, b = 2π*rand(2) .- π
    @test P(a,b) ≈ P(BigFloat(a),b)
end

#test cheby derivative
p = ChebyshevInterpolator(x->sin(x^2), -2π, 2π, 128)
d = chebyderiv(p)
for i = 1:10
    a = 4π*rand() - 2π
    @test d(a) - 2*a*cos(a^2) < 1e-6
end