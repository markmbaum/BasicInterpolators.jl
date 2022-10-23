export CubicSplineInterpolator, BicubicSplineInterpolator

#-------------------------------------------------------------------------------
# piecewise cubics with continuous derivatives (splines!)

struct CubicSplineInterpolator{T,B} <: OneDimensionalInterpolator
    r::InterpolatorRange{T}
    a::Vector{T}
    b::Vector{T}
    c::Vector{T}
    d::Vector{T}
    boundaries::B
end

#-------------------------------------------------------------------------------
# constructors

"""
    CubicSplineInterpolator(x, y, boundaries=StrictBoundaries())

Construct a `CubicSplineInterpolator` for the points defined by coordinates `x` and values `y`. This constructor creates a natural spline, where the second derivative is set to zero at the boundaries.
"""
function CubicSplineInterpolator(x, y, boundaries::AbstractBoundaries=StrictBoundaries())
    #construct the underlying range, triggering some checks
    T = promote_type(eltype(x), eltype(y))
    r = InterpolatorRange(x, y)
    n = r.n
    #compute coefficients
    #Burden, Richard L., and J. Douglas Faires. Numerical Analysis. 2011.
    a = collect(T, r.y)
    b = zeros(T, n - 1)
    d = zeros(T, n - 1)
    h = diff(r.x)
    α = zeros(T, n-1)
    for i = 2:n-1
        α[i] = 3*(a[i+1] - a[i])/h[i] - 3*(a[i] - a[i-1])/h[i-1]
    end
    c = zeros(T, n)
    l = ones(T, n)
    μ = zeros(T, n)
    z = zeros(T, n)
    l[1] = 1
    for i = 2:n-1
        l[i] = 2*(x[i+1] - x[i-1]) - h[i-1]*μ[i-1]
        μ[i] = h[i]/l[i]
        z[i] = (α[i] - h[i-1]*z[i-1])/l[i]
    end
    for j = n-1:-1:1
        c[j] = z[j] - μ[j]*c[j+1]
        b[j] = (a[j+1] - a[j])/h[j] - h[j]*(c[j+1] + 2*c[j])/3
        d[j] = (c[j+1] - c[j])/(3*h[j])
    end
    a = a[1:end-1]
    c = c[1:end-1]
    #static arrays
    #coef = Vector{NTuple{4,T}}(undef,n-1)
    #for i = 1:n-1
    #    coef[i] = (a[i], b[i], c[i], d[i])
    #end
    #construct the object
    CubicSplineInterpolator(r, a, b, c, d, boundaries)
end

"""
    CubicSplineInterpolator(x, y, dy₁, dyₙ, boundaries=StrictBoundaries())

Construct a `CubicSplineInterpolator` for the points defined by coordinates `x` and values `y`. This constructor creates a clamped spline, where the first derivatives at the boundaries are set by `dy₁` and `dyₙ`.
"""
function CubicSplineInterpolator(x,
                                 y,
                                 dy₁::Real,
                                 dyₙ::Real,
                                 boundaries::AbstractBoundaries=StrictBoundaries())
    #construct the underlying range, triggering some checks
    T = promote_type(eltype(x), eltype(y))
    r = InterpolatorRange(x, y)
    n = r.n
    #compute coefficients
    #Burden, Richard L., and J. Douglas Faires. Numerical Analysis. 2011.
    a = collect(T, r.y)
    b = zeros(T, n - 1)
    d = zeros(T, n - 1)
    h = diff(r.x)
    α = zeros(T, n)
    α[1] = 3*(a[2] - a[1])/h[1] - 3*dy₁
    for i = 2:n-1
        α[i] = 3*(a[i+1] - a[i])/h[i] - 3*(a[i] - a[i-1])/h[i-1]
    end
    α[n] = 3*dyₙ - 3*(a[n] - a[n - 1])/h[n-1]
    c = zeros(T, n)
    l = zeros(T, n)
    μ = zeros(T, n)
    z = zeros(T, n)
    l[1] = 2*h[1]
    μ[1] = 0.5
    z[1] = α[1]/l[1]
    for i = 2:n-1
        l[i] = 2*(x[i+1] - x[i-1]) - h[i-1]*μ[i-1]
        μ[i] = h[i]/l[i]
        z[i] = (α[i] - h[i-1]*z[i-1])/l[i]
    end
    l[n] = h[n-1]*(2 - μ[n-1])
    z[n] = (α[n] - h[n-1]*z[n-1])/l[n]
    c[n] = z[n]
    for j = n-1:-1:1
        c[j] = z[j] - μ[j]*c[j+1]
        b[j] = (a[j+1] - a[j])/h[j] - h[j]*(c[j+1] + 2*c[j])/3
        d[j] = (c[j+1] - c[j])/(3*h[j])
    end
    a = a[1:end-1]
    c = c[1:end-1]
    #static arrays
    #coef = Vector{NTuple{4,T}}(undef,n-1)
    #for i = 1:n-1
    #    coef[i] = (a[i], b[i], c[i], d[i])
    #end
    #construct the object
    CubicSplineInterpolator(r, a, b, c, d, boundaries)
end

"""
    CubicSplineInterpolator(f, xa, xb, n, boundaries=StrictBoundaries())

Construct a `CubicSplineInterpolator` for the function `f` using `n` evenly spaced function evaluations in the range [`xa`,`xb`]. A natural spline is created.
"""
function CubicSplineInterpolator(f,
                                 xa::Real,
                                 xb::Real,
                                 n::Int,
                                 boundaries::AbstractBoundaries=StrictBoundaries())
    linstruct(CubicSplineInterpolator, f, xa, xb, n, boundaries)
end

function (ϕ::CubicSplineInterpolator)(x)
    #enforce boundaries if desired
    ϕ.boundaries(x, ϕ.r.xa, ϕ.r.xb)
    #find the interpolation point
    i = findcell(x, ϕ)
    #offset from the nearest lower point
    @inbounds ξ = x - ϕ.r.x[i]
    #evaluate polynomial
    @inbounds ϕ.a[i] + ϕ.b[i]*ξ + ϕ.c[i]*ξ^2 + ϕ.d[i]*ξ^3
end

Base.getindex(ϕ::CubicSplineInterpolator, i) = ϕ.r.y[i]
function Base.copy(ϕ::CubicSplineInterpolator)
    CubicSplineInterpolator(ϕ.r, ϕ.coef, ϕ.boundaries)
end
