function vandermonde_error(n)
    a = collect(1:n)
    x = LinRange(-1, 1, n)
    y = zeros(n)
    for i = 1:n
        y[i] = sum(ntuple(j->a[j]*x[i]^(j-1), n))
    end
    c = vandermonde(x, y)
    return maximum(abs.(c - a))
end

for n = 2:5
    @test vandermonde_error(n) < 0.01
end
