export quadratic, cubic, neville, vandermonde,
       LinearInterpolator, CubicInterpolator, CubicSplineInterpolator,
       ParametricCurveInterpolator,
       BilinearInterpolator, BicubicInterpolator, BicubicSplineInterpolator

"""
    quadratic(x, xp, yp)

Perform quadratic polynomial interpolation of the points defined by coordinates `xp` and values `yp`, at the coordinate `x`, using Neville's algorithm. `xp` and `yp` must both contain three points.
"""
function quadratic(x::Real,
                   xp::AbstractVector{<:Real},
                   yp::AbstractVector{<:Real})::Float64
    @assert length(xp) == length(yp) == 3 "need 3 points for quadratic interpolation"
    #first stage
    p12 = ((x - xp[2])*yp[1] + (xp[1] - x)*yp[2])/(xp[1] - xp[2])
    p23 = ((x - xp[3])*yp[2] + (xp[2] - x)*yp[3])/(xp[2] - xp[3])
    #final stage
    ((x - xp[3])*p12 + (xp[1] - x)*p23)/(xp[1] - xp[3])
end

"""
    cubic(x, xp, yp)

Perform cubic polynomial interpolation of the points defined by coordinates `xp` and values `yp`, at the coordinate `x`, using Neville's algorithm. `xp` and `yp` must both contain four points.
"""
function cubic(x::Real,
               xp::AbstractVector{<:Real},
               yp::AbstractVector{<:Real})::Float64
    @assert length(xp) == length(yp) == 4 "need 4 points for cubic interpolation"
    #first stage
    p12 = ((x - xp[2])*yp[1] + (xp[1] - x)*yp[2])/(xp[1] - xp[2])
    p23 = ((x - xp[3])*yp[2] + (xp[2] - x)*yp[3])/(xp[2] - xp[3])
    p34 = ((x - xp[4])*yp[3] + (xp[3] - x)*yp[4])/(xp[3] - xp[4])
    #second stage
    p123 = ((x - xp[3])*p12 + (xp[1] - x)*p23)/(xp[1] - xp[3])
    p234 = ((x - xp[4])*p23 + (xp[2] - x)*p34)/(xp[2] - xp[4])
    #final stage
    ((x - xp[4])*p123 + (xp[1] - x)*p234)/(xp[1] - xp[4])
end

"""
    neville(x, xp, yp)

Perform polynomial interpolation of the points defined by coordinates `xp` and values `yp`, at the coordinate `x`, using Neville's algorithm with as many points as are provided. `xp` and `yp` must have the same length. With only 3 or 4 points the [`quadratic`](@ref) and [`cubic`](@ref) functions will be faster.
"""
function neville(x::Real,
                 xp::AbstractVector{<:Real},
                 yp::AbstractVector{<:Real})::Float64
    @assert length(xp) == length(yp) "can't Neville with vectors of different lengths"
    n = length(xp)
    P = zeros(n, n)
    P[:,1] = yp
    for i = 2:n
        for j = i:n
            P[j,i] = ((x - xp[j-i+1])*P[j,i-1] - (x - xp[j])*P[j-1,i-1])/(xp[j] - xp[j-i+1])
        end
    end
    return P[n,n]
end

"""
    vandermonde(x, y)

Generate the coefficients of an arbitrary order polynomial passing through the ponts defined by coordinates `x` and value `y`. For n points, n coefficients ``[c_0, c_1, ..., c_{n-1}]`` are returned forming the polynomial ``c_0 + c_1x + ... + c_{n-1}x^{n-1}``

!!! warning

    Solving for the the coefficients of a high-order polynomial is a notoriously ill-conditioned problem. It is not recommended for orders greater than 5 or 6, although it depends on the application. If you must interpolate with a high-order polynomial, it's better to use the [`neville`](@ref) function instead of computing coefficients.
"""
function vandermonde(x::AbstractVector{<:Real},
                     y::AbstractVector{<:Real})::Vector{Float64}
    @assert length(x) == length(y) "length of x must equal length of y"
    n = length(x)
    s = zeros(n)
    c = zeros(n)
    s[n] = -x[1]
    for i = 2:n
        for j = n-i+1:n-1
            s[j] -= x[i]*s[j+1]
        end
        s[n] -= x[i]
    end
    for j = 1:n
        ϕ = n
        for k = n:-1:2
            ϕ = (k-1)*s[k] + x[j]*ϕ
        end
        f = y[j]/ϕ
        b = 1.0
        for k = n:-1:1
            c[k] += b*f
            b = s[k] + x[j]*b
        end
    end
    return c
end

#-------------------------------------------------------------------------------
# simple piecewise linear interpolator

struct LinearInterpolator
    r::InterpolatorRange
end

"""
    LinearInterpolator(x, y)

Construct a `LinearInterpolator` for the points defined by coordinates `x` and values `y`
"""
function LinearInterpolator(x::AbstractVector{<:Real}, y::AbstractVector{<:Real})
    LinearInterpolator(InterpolatorRange(x, y))
end

"""
    LinearInterpolator(f, xa, xb, n)

Construct a `LinearInterpolator` for the function `f` using `n` evenly spaced function evaluations in the range [`xa`,`xb`]
"""
function LinearInterpolator(f::Function, xa::Real, xb::Real, n::Int)
    linstruct(LinearInterpolator, f, xa, xb, n)
end

function (ϕ::LinearInterpolator)(x::Real, bounds::Bool=true)::Float64
    #enforce boundaries if desired
    enforcebounds(x, ϕ.r.xa, ϕ.r.xb, bounds)
    #find the interpolation point
    i = findcell(x, ϕ.r.x)
    #interpolate
    (x - ϕ.r.x[i])*(ϕ.r.y[i+1] - ϕ.r.y[i])/(ϕ.r.x[i+1] - ϕ.r.x[i]) + ϕ.r.y[i]
end

#-------------------------------------------------------------------------------
# piecewise cubic interpolator without continuous derivatives (not splines)

struct CubicInterpolator
    r::InterpolatorRange
end

"""
    CubicInterpolator(x, y)

Construct a `CubicInterpolator` for the points defined by coordinates `x` and values `y`
"""
function CubicInterpolator(x::AbstractVector{<:Real}, y::AbstractVector{<:Real})
    @assert length(x) > 3 "can't do cubic interpolation with <4 points"
    CubicInterpolator(InterpolatorRange(x, y))
end

"""
    CubicInterpolator(f, xa, xb, n)

Construct a `CubicInterpolator` for the function `f` using `n` evenly spaced function evaluations in the range [`xa`,`xb`]
"""
function CubicInterpolator(f::Function, xa::Real, xb::Real, n::Int)
    linstruct(CubicInterpolator, f, xa, xb, n)
end

function (ϕ::CubicInterpolator)(x::Real, bounds::Bool=true)::Float64
    #enforce boundaries if desired
    enforcebounds(x, ϕ.r.xa, ϕ.r.xb, bounds)
    #find the interpolation point
    i = findcell(x, ϕ.r.x)
    #determine which points to neville
    if i == 1
        I = 1:4
    elseif i == ϕ.r.n - 1
        I = ϕ.r.n-3:ϕ.r.n
    else
        I = i-1:i+2
    end
    #interpolate
    cubic(x, view(ϕ.r.x,I), view(ϕ.r.y,I))
end

#-------------------------------------------------------------------------------
# piecewise cubics with continuous derivatives (splines!)

struct CubicSplineInterpolator
    r::InterpolatorRange
    α::Array{Float64,2} #4xN array of cubic polynomial coefficients
end

"""
    CubicSplineInterpolator(x, y)

Construct a `CubicSplineInterpolator` for the points defined by coordinates `x` and values `y`. This constructor creates a natural spline, where the second derivative is set to zero at the boundaries.
"""
function CubicSplineInterpolator(x::AbstractVector{<:Real},
                                 y::AbstractVector{<:Real})
    #construct the underlying range, triggering some checks
    r = InterpolatorRange(x, y)
    n = r.n
    #compute coefficients
    #Burden, Richard L., and J. Douglas Faires. Numerical Analysis. 2011.
    a = collect(Float64, r.y)
    b = zeros(n - 1)
    d = zeros(n - 1)
    h = diff(r.x)
    α = zeros(n-1)
    for i = 2:n-1
        α[i] = 3*(a[i+1] - a[i])/h[i] - 3*(a[i] - a[i-1])/h[i-1]
    end
    c = zeros(n)
    l = ones(n)
    μ = zeros(n)
    z = zeros(n)
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
    #construct the object
    CubicSplineInterpolator(r, hcat(a, b, c, d)')
end

"""
    CubicSplineInterpolator(x, y, y₁′, yₙ′)

Construct a `CubicSplineInterpolator` for the points defined by coordinates `x` and values `y`. This constructor creates a clamped spline, where the first derivatives at the boundaries are set by `y₁′` and `yₙ′`.
"""
function CubicSplineInterpolator(x::AbstractVector{<:Real},
                                 y::AbstractVector{<:Real},
                                 y₁′::Real,
                                 yₙ′::Real)
    #construct the underlying range, triggering some checks
    r = InterpolatorRange(x, y)
    n = r.n
    #compute coefficients
    #Burden, Richard L., and J. Douglas Faires. Numerical Analysis. 2011.
    a = collect(Float64, r.y)
    b = zeros(n - 1)
    d = zeros(n - 1)
    h = diff(r.x)
    α = zeros(n)
    α[1] = 3*(a[2] - a[1])/h[1] - 3*y₁′
    for i = 2:n-1
        α[i] = 3*(a[i+1] - a[i])/h[i] - 3*(a[i] - a[i-1])/h[i-1]
    end
    α[n] = 3*yₙ′ - 3*(a[n] - a[n - 1])/h[n-1]
    c = zeros(n)
    l = zeros(n)
    μ = zeros(n)
    z = zeros(n)
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
    #coefficient table
    α = hcat(a[1:end-1], b, c[1:end-1], d)'
    #construct the object
    CubicSplineInterpolator(r, α)
end

"""
    CubicSplineInterpolator(f, xa, xb, n)

Construct a `CubicSplineInterpolator` for the function `f` using `n` evenly spaced function evaluations in the range [`xa`,`xb`]. A natural spline is created.
"""
function CubicSplineInterpolator(f::Function, xa::Real, xb::Real, n::Int)
    linstruct(CubicSplineInterpolator, f, xa, xb, n)
end

function (ϕ::CubicSplineInterpolator)(x::Real, bounds::Bool=true)::Float64
    #enforce boundaries if desired
    enforcebounds(x, ϕ.r.xa, ϕ.r.xb, bounds)
    #find the interpolation point
    i = findcell(x, ϕ.r.x)
    #offset from the nearest lower point
    ξ = x - ϕ.r.x[i]
    #evaluate polynomial
    ϕ.α[1,i] + ϕ.α[2,i]*ξ + ϕ.α[3,i]*ξ^2 + ϕ.α[4,i]*ξ^3
end

#-------------------------------------------------------------------------------
# parametric cubic splines for n dimensions

struct ParametricCurveInterpolator
    #interpolators for each dimension
    Φ::Vector{CubicSplineInterpolator}
    #number of dimensions
    N::Int64
end

"""
    ParametricCurveInterpolator(V...)

Construct an interpolator for a set of points in arbitrary dimensions that defines a one-dimensional curve, using natural cubic splines in each dimension.
"""
function ParametricCurveInterpolator(V::AbstractArray{<:Real}...)
    #standardize types
    V = [collect(Float64, v) for v ∈ V]
    #check lengths
    L = length(V[1])
    N = length(V)
    for i = 2:N
        @assert length(V[i]) == L "coordinate vectors must all be equal length"
    end
    #compute chordlengths along the curve
    s = zeros(L)
    for i = 2:L
        #segment distance
        d = 0
        for j = 1:N
            d += (V[j][i] - V[j][i-1])^2
        end
        d = sqrt(d)
        #incremental chord length
        s[i] = s[i-1] + d
    end
    #normalize s to have a maximum of 1
    s /= s[end]
    #construct interpolators
    Φ = [CubicSplineInterpolator(s, v) for v ∈ V]
    #construct global interpolator
    ParametricCurveInterpolator(Φ, N)
end

function (C::ParametricCurveInterpolator)(s::Real,
                                          bounds::Bool=true)::Vector{Float64}
    @assert 0 <= s <= 1 "interpolation coordinate must be in [0,1] for parametric curves"
    [ϕ(s, bounds) for ϕ ∈ C.Φ]
end

#-------------------------------------------------------------------------------
# bilinear interpolator

struct BilinearInterpolator
    G::InterpolatorGrid
end

"""
    BilinearInterpolator(x, y, Z)

Construct a `BilinearInterpolator` for the grid of points points defined by coordinates (`x`,`y`) and values `Z`.
"""
function BilinearInterpolator(x::AbstractVector{<:Real},
                              y::AbstractVector{<:Real},
                              Z::AbstractArray{<:Real,2})
    BilinearInterpolator(InterpolatorGrid(x, y, Z))
end

"""
    BilinearInterpolator(f, xa, xb, nx, ya, yb, ny)

Construct a `BilinearInterpolator` for the function `f` using a grid of `nx` points evenly spaced on the first axis in [`xa`,`xb`] and `ny` points evenly spaced on the second axis in [`ya`,`yb`].
"""
function BilinearInterpolator(f::Function,
                              xa::Real, xb::Real, nx::Int,
                              ya::Real, yb::Real, ny::Int)
    linstruct(BilinearInterpolator, f, xa, xb, nx, ya, yb, ny)
end

function (Φ::BilinearInterpolator)(x::Real,
                                   y::Real,
                                   bounds::Bool=true)::Float64
    #enforce boundaries if desired
    enforcebounds(x, Φ.G.xa, Φ.G.xb, y, Φ.G.ya, Φ.G.yb, bounds)
    #find the proper grid box to interpolate inside
    i = findcell(x, Φ.G.x)
    j = findcell(y, Φ.G.y)
    #name stuff for clarity
    xg, yg, Z = Φ.G.x, Φ.G.y, Φ.G.Z
    #interpolate along axis 1 first
    t = (x - xg[i])*(Z[i+1,j]   - Z[i,j]  )/(xg[i+1] - xg[i]) + Z[i,j]
    u = (x - xg[i])*(Z[i+1,j+1] - Z[i,j+1])/(xg[i+1] - xg[i]) + Z[i,j+1]
    #then finish by interpolating between the interpolated values along axis 2
    (y - yg[j])*(u - t)/(yg[j+1] - yg[j]) + t
end

#-------------------------------------------------------------------------------
# bicubic interpolator

struct BicubicInterpolator
    G::InterpolatorGrid
end

"""
    BicubicInterpolator(x, y, Z)

Construct a `BicubicInterpolator` for the grid of points points defined by coordinates (`x`,`y`) and values `Z`.
"""
function BicubicInterpolator(x::AbstractVector{<:Real},
                             y::AbstractVector{<:Real},
                             Z::AbstractArray{<:Real,2})
    #insist on at least 4 points in each dimension
    @assert (length(x) > 3) & (length(y) > 3) "bicubic interpolation requires at least 4 points in each dimension"
    BicubicInterpolator(InterpolatorGrid(x, y, Z))
end

"""
    BicubicInterpolator(f, xa, xb, nx, ya, yb, ny)

Construct a `BicubicInterpolator` for the function `f` using a grid of `nx` points evenly spaced on the first axis in [`xa`,`xb`] and `ny` points evenly spaced on the second axis in [`ya`,`yb`].
"""
function BicubicInterpolator(f::Function,
                             xa::Real, xb::Real, nx::Int,
                             ya::Real, yb::Real, ny::Int)
    linstruct(BicubicInterpolator, f, xa, xb, nx, ya, yb, ny)
end

function (Φ::BicubicInterpolator)(x::Real,
                                  y::Real,
                                  bounds::Bool=true)::Float64
    #enforce boundaries if desired
    enforcebounds(x, Φ.G.xa, Φ.G.xb, y, Φ.G.ya, Φ.G.yb, bounds)
    #find the proper grid box to interpolate inside
    i = findcell(x, Φ.G.x)
    j = findcell(y, Φ.G.y)
    #get indices of points along axis 1 to use for interpolation
    if i == 1
        I = 1:4
    elseif i == Φ.G.nx - 1
        I = Φ.G.nx-3:Φ.G.nx
    else
        I = i-1:i+2
    end
    #get indices along axis 2 where initial 4 interpolations occur
    if j == 1
        J = 1:4
    elseif j == Φ.G.ny - 1
        J = Φ.G.ny-3:Φ.G.ny
    else
        J = j-1:j+2
    end
    #perform initial 4 interpolations
    zx = [neville(x, view(Φ.G.x, I), view(Φ.G.Z, I, J[k])) for k=1:4]
    #final interpolation
    neville(y, view(Φ.G.y,J), zx)
end

#-------------------------------------------------------------------------------
# cubic splines on a regular grid

struct BicubicSplineInterpolator
    G::InterpolatorGrid
    α::Array{Array{Float64,2},2}
end

"""
    BicubicSplineInterpolator(x, y, Z)

Construct a `BicubicSplineInterpolator` for the grid of points points defined by coordinates (`x`,`y`) and values `Z`.
"""
function BicubicSplineInterpolator(x::AbstractVector{<:Real},
                                   y::AbstractVector{<:Real},
                                   Z::AbstractArray{<:Real,2})
    nx, ny = size(Z)
    #insist on at least 4 points in each dimension
    @assert nx >= 3 "bicubic interpolation requires at least 3 points along axis 1"
    @assert ny >= 3 "bicubic interpolation requires at least 3 points along axis 2"
    #insist on even spacing along both axes
    @assert all(diff(diff(x)) .< 1e-8*maximum(abs.(x))) "grid spacing along axis 1 must be uniform"
    @assert all(diff(diff(y)) .< 1e-8*maximum(abs.(y))) "grid spacing along axis 2 must be uniform"
    #space for derivatives
    dx = zeros(nx, ny)
    dy = zeros(nx, ny)
    dxy = zeros(nx, ny)
    #first derivatives for internal nodes
    for i = 2:nx-1
        for j = 1:ny
            dx[i,j] = (Z[i+1,j] - Z[i-1,j])/2
        end
    end
    for i = 1:nx
        for j = 2:ny-1
            dy[i,j] = (Z[i,j+1] - Z[i,j-1])/2
        end
    end
    #first derivatives for bounary nodes
    for i = 1:nx
        dy[i,1] = -3*Z[i,1]/2 + 2*Z[i,2] - Z[i,3]/2
        dy[i,ny] = Z[i,ny-2]/2 - 2*Z[i,ny-1] + 3*Z[i,ny]/2
    end
    for j = 1:ny
        dx[1,j] = -3*Z[1,j]/2 + 2*Z[2,j] - Z[3,j]/2
        dx[nx,j] = Z[nx-2,j]/2 - 2*Z[nx-1,j] + 3*Z[nx,j]/2
    end
    #mixed second order derivatives at internal nodes
    for i = 1:nx
        for j = 2:ny-1
            dxy[i,j] = (dx[i,j+1] - dx[i,j-1])/2
        end
    end
    #mixed second deriv along the sides
    for i = 1:nx
        dxy[i,1] = -3*dx[i,1]/2 + 2*dx[i,2] - dx[i,3]/2
        dxy[i,ny] = dx[i,ny-2]/2 - 2*dx[i,ny-1] + 3*dx[i,ny]/2
    end
    #matrix needed to compute coefficients and its transpose
    A = [1. 0. 0. 0.; 0. 0. 1. 0.; -3. 3. -2. -1.; 2. -2. 1. 1.]
    Aᵀ = transpose(A)
    #space for coefficients
    f = zeros(4, 4)
    α = Array{Array{Float64,2},2}(undef, nx-1, ny-1)
    for i = 1:nx-1
        for j = 1:ny-1
            #load the matrix of 16 values and derivatives
            f[1,1] = Z[i,j]
            f[1,2] = Z[i,j+1]
            f[1,3] = dy[i,j]
            f[1,4] = dy[i,j+1]
            f[2,1] = Z[i+1,j]
            f[2,2] = Z[i+1,j+1]
            f[2,3] = dy[i+1,j]
            f[2,4] = dy[i+1,j+1]
            f[3,1] = dx[i,j]
            f[3,2] = dx[i,j+1]
            f[3,3] = dxy[i,j]
            f[3,4] = dxy[i,j+1]
            f[4,1] = dx[i+1,j]
            f[4,2] = dx[i+1,j+1]
            f[4,3] = dxy[i+1,j]
            f[4,4] = dxy[i+1,j+1]
            #get the 16 double spline coefficients
            α[i,j] = A*f*Aᵀ
        end
    end
    BicubicSplineInterpolator(InterpolatorGrid(x, y, Z), α)
end

"""
    BicubicSplineInterpolator(f, xa, xb, nx, ya, yb, ny)

Construct a `BicubicSplineInterpolator` for the function `f` using a grid of `nx` points evenly spaced on the first axis in [`xa`,`xb`] and `ny` points evenly spaced on the second axis in [`ya`,`yb`].
"""
function BicubicSplineInterpolator(f::Function,
                                   xa::Real, xb::Real, nx::Int,
                                   ya::Real, yb::Real, ny::Int)
    linstruct(BicubicSplineInterpolator, f, xa, xb, nx, ya, yb, ny)
end

function (Φ::BicubicSplineInterpolator)(x::Real,
                                        y::Real,
                                        bounds::Bool=true)::Float64
    #enforce boundaries if desired
    enforcebounds(x, Φ.G.xa, Φ.G.xb, y, Φ.G.ya, Φ.G.yb, bounds)
    #find the proper grid box to interpolate inside
    i = findcell(x, Φ.G.x)
    j = findcell(y, Φ.G.y)
    #get the coefficients
    α = Φ.α[i,j]
    #offsets
    Δx = (x - Φ.G.x[i])/(Φ.G.x[i+1] - Φ.G.x[i])
    Δy = (y - Φ.G.y[j])/(Φ.G.y[j+1] - Φ.G.y[j])
    #powers of the offsets
    x1, x2, x3 = Δx, Δx^2, Δx^3
    y1, y2, y3 = Δy, Δy^2, Δy^3
    #final interpolation calculation (fastest written out like this)
    ( α[1,1]    + α[2,1]*x1    + α[3,1]*x2    + α[4,1]*x3
    + α[1,2]*y1 + α[2,2]*x1*y1 + α[3,2]*x2*y1 + α[4,2]*x3*y1
    + α[1,3]*y2 + α[2,3]*x1*y2 + α[3,3]*x2*y2 + α[4,3]*x3*y2
    + α[1,4]*y3 + α[2,4]*x1*y3 + α[3,4]*x2*y3 + α[4,4]*x3*y3)
end
