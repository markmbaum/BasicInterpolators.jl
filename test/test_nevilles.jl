function neville_error(n::Int, f::Function)
    #random polynomial coefficents
    a = rand(n)
    #the polynomial
    p(x) = sum(a[i]*x^(i-1) for i=1:n)
    #samples
    x = rand(n)
    y = p.(x)
    #dense evaluations
    xx = collect(LinRange(minimum(x)-1, maximum(x)+1, 1000))
    pyy = p.(xx)
    fyy = [f(xx[i], x, y) for i=1:length(xx)]
    #error
    maximum(abs.(pyy - fyy))
end

@test neville_error(3, quadratic) < 1e-6
@test neville_error(4, cubic) < 1e-6
for i = 5:6
    @test neville_error(i, neville) < 1e-4
end
