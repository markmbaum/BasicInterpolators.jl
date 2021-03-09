F(x) = sin(x[1]) + cos(x[2]) + sin(x[3])
npts = 10000
X = 6*rand(npts, 3) .- 3
y = [F(X[i,:]) for i = 1:npts]

P = ShepardInterpolator(X, y)
Q = RBFInterpolator(X, y, 1.0)

x = collect(LinRange(-3, 3, 50))
y = collect(LinRange(-3, 3, 50))
z = collect(LinRange(-3, 3, 50))

errp = zeros(50, 50, 50)
errq = zeros(50, 50, 50)
for i = 1:length(x)
    for j = 1:length(y)
        for k = 1:length(z)
            f = F([x[i], y[j], z[k]])
            errp[i,j,k] = abs(P(x[i], y[j], z[k]) - f)
            errq[i,j,k] = abs(Q(x[i], y[j], z[k]) - f)
        end
    end
end

@test maximum(errp) < 2
@test maximum(errq) < 0.1
