function neville_error(n::Int, F::Function)
    #random polynomial coefficents
    a = rand(n)
    #the polynomial
    p(x) = sum(a[i]*x^(i-1) for i=1:n)
    #samples
    x = LinRange(-1, 1, n)
    y = p.(x)
    #closure
    f(ξ) = F(ξ, x, y)
    #maximum error in [-1,1]
    ξ = -1:0.001:1
    maximum( @. abs(p(ξ) - f(ξ)) )
end

#simple checks on polynomial interpolation accuracy
@test neville_error(3, quadratic) < 1e-3
@test neville_error(4, cubic) < 1e-3
@test neville_error(5, neville) < 1e-3
@test neville_error(6, neville) < 1e-3