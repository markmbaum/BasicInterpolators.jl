F(x) = sin(x[1])*cos(x[2])
npts = 1000
rng = MersenneTwister(1)
X = 2π*rand(rng, npts, 2) .- π
y = [F(X[i,:]) for i = 1:npts]

P = ShepardInterpolator(X, y)
Q = RBFInterpolator(X, y, 1.0)

x = LinRange(-π, π, 10)
y = copy(x)

errp = zeros(length(x), length(x))
errq = zeros(length(x), length(x))
for i ∈ eachindex(x)
    for j ∈ eachindex(y)
        f = F((x[i], y[j]))
        errp[i,j] = abs(P(x[i], y[j]) - f)
        errq[i,j] = abs(Q(x[i], y[j]) - f)
    end
end

@test maximum(errp) < 1
@test maximum(errq) < 0.1
