module BasicInterpolators

#for in-place matrix multiplication in chebyshev interpolators
using LinearAlgebra: mul!

#-------------------------------------------------------------------------------
#general things
export quadratic, cubic, neville

function findcell(q::Real, V::AbstractVector{Float64})::Int64
    n = length(V)
    #handle boundaries
    if q <= V[1]
        return 1
    end
    if q >= V[n]
        return n - 1 #last cell starts with second to last point
    end
    #bisection search for the containing cell
    ilo = 0
    ihi = n
    while ihi - ilo > 1
        imid = (ihi + ilo) ÷ 2
        if V[imid] > q
            ihi = imid
        else
            ilo = imid
        end
    end
    return ilo
end

"""
    quadratic(x::Real, xp::AbstractVector{<:Real}, yp::AbstractVector{<:Real})

Perform quadratic polynomial interpolation of `xp` and `yp` at the point `x`, using Neville's algorithm. `xp` and `yp` must both contain three points.
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
    cubic(x::Real, xp::AbstractVector{<:Real}, yp::AbstractVector{<:Real})

Perform cubic polynomial interpolation of `xp` and `yp` at the point `x`, using Neville's algorithm. `xp` and `yp` must both contain four points.
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
    neville(x::Real, xp::AbstractVector{<:Real}, yp::AbstractVector{<:Real})

Perform polynomial interpolation of `xp` and `yp` at the point `x`, using Neville's algorithm with as many points as are provided. `xp` and `yp` must have the same length. With only 3 or 4 points (quadratic or cubic polynomial interpolation) the [`quadratic`](@ref) and [`cubic`](@ref) functions will be faster.
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

#-------------------------------------------------------------------------------
# functions for handling boundaries

function enforcebounds(x, xa, xb, enforce::Bool)
    if enforce
        @assert xa <= x <= xb "Interpolation location $x out of range [$xa,$xb]"
    end
end

function enforcebounds(x, xa, xb, y, ya, yb, enforce::Bool)
    if enforce
        @assert xa <= x <= xb "Interpolation location $x out of range [$xa,$xb] on axis 1"
        @assert ya <= y <= yb "Interpolation location $y out of range [$ya,$yb] on axis 2"
    end
end

#-------------------------------------------------------------------------------
# all the basic chebyshev stuff
export chebygrid

ξ2x(ξ::Real, a::Real, b::Real)::Float64 = (ξ + 1)*((b - a)/2) + a

x2ξ(x::Real, a::Real, b::Real)::Float64 = 2*(x - a)/(b - a) - 1

x2θ(x::Real, a::Real, b::Real)::Float64 = acos(x2ξ(x, a, b))

"""
    chebygrid(n::Int)

Create an array of `n` chebyshev nodes in [-1,1]
"""
function chebygrid(n::Int)::Vector{Float64}
    cos.(π*(n-1:-1:0)/(n-1))
end

"""
    chebygrid(xa::Real, xb::Real, n::Int)

Create an array of `n` chebyshev nodes in [`xa`,`xb`]
"""
function chebygrid(xa::Real, xb::Real, n::Int)::Vector{Float64}
    ξ2x.(chebygrid(n), xa, xb)
end

"""
    chebygrid(xa::Real, xb::Real, nx::Int, ya::Real, yb::Real, ny::Int

Create a two-dimensional grid of chebyshev nodes using `nx` points along the first axis, in [`xa`,`xb`], and `ny` points along the second axis, in [`ya`,`yb`].
"""
function chebygrid(xa::Real, xb::Real, nx::Int, ya::Real, yb::Real, ny::Int)
    x = chebygrid(xa, xb, nx)
    y = chebygrid(ya, yb, ny)
    X = x .* ones(length(y))'
    Y = y' .* ones(length(x))
    return X, Y
end

function ischebygrid(x::AbstractVector{<:Real})::Bool
    all(x2ξ.(x, minimum(x), maximum(x)) .- chebygrid(length(x)) .< 1e-8)
end

cheby(ξ::Real, k::Int)::Float64 = cos(k*acos(ξ))

cheby(x::Real, k::Int, xa::Real, xb::Real)::Float64 = cheby(x2ξ(x, xa, xb), k)

function cheby(x::Real,
               a::AbstractVector{<:Real},
               xa::Real,
               xb::Real)::Float64
    ξ = x2ξ(x, xa, xb)
    y = a[1]
    for k = 2:length(a)
        y += a[k]*cheby(ξ, k-1)
    end
    return y
end

function chebymatrix(n::Int)::Array{Float64,2}
    @assert n > 1 "can't construct cheby matrix smaller than 2 x 2"
    A = zeros(n,n)
    ξ = chebygrid(n)
    for j = 1:n
        for k = 1:n
            A[j,k] = cheby(ξ[j], k-1)
        end
    end
    return A
end

function fce!(v::AbstractVector{<:Real}, θ::Real, n::Int)
    #only need to evaluate cosine directly at odd multiples of θ (even indices)
    for i = 2:2:n
        v[i] = cos((i-1)*θ)
        #then evaluate 2ⁿθ terms recursively
        j = 2*i - 1
        w = v[i]
        while j <= n
            v[j] = 2*w^2 - 1 #double-angle formula, cos(2θ) = 2cos²(θ) - 1
            w = v[j]
            j = 2*j - 1
        end
    end
end

#-------------------------------------------------------------------------------
# a base structure for 1d interpolations

struct InterpolatorRange
    #number of points
    n::Int64
    #grid/sample points
    x::Vector{Float64}
    #lowest sample value
    xa::Float64
    #highest sample value
    xb::Float64
    #values at sample points
    y::Vector{Float64}
end

function InterpolatorRange(x::AbstractVector{<:Real}, y::AbstractVector{<:Real})
    #check for basic problems
    rangecheck(x, y)
    #make sure types are uniform and make copies along the way
    x = convert.(Float64, x)
    y = convert.(Float64, y)
    #grid properties
    n = length(x)
    xa = minimum(x)
    xb = maximum(x)
    #construct
    InterpolatorRange(n, x, xa, xb, y)
end

function rangecheck(x::AbstractVector{<:Real}, y::AbstractVector{<:Real})
    @assert length(x) == length(y) "lengths mismatched"
    @assert all(diff(x) .> 0.0) "grid points must be in ascending order"
end

function linstruct(T::Type, f::Function, xa::Real, xb::Real, n::Int)
    x = collect(LinRange(xa, xb, n))
    y = f.(x)
    T(x, y)
end

#-------------------------------------------------------------------------------
# simple piecewise linear interpolator
export LinearInterpolator

struct LinearInterpolator
    r::InterpolatorRange
end

"""
    LinearInterpolator(x::AbstractVector{<:Real}, y::AbstractVector{<:Real})

Construct a `LinearInterpolator` for the points in `x` and `y`
"""
function LinearInterpolator(x::AbstractVector{<:Real}, y::AbstractVector{<:Real})
    LinearInterpolator(InterpolatorRange(x, y))
end

"""
    LinearInterpolator(f::Function, xa::Real, xb::Real, n::Int)

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
export CubicInterpolator

struct CubicInterpolator
    r::InterpolatorRange
end

"""
    CubicInterpolator(x::AbstractVector{<:Real}, y::AbstractVector{<:Real})

Construct a `CubicInterpolator` for the points in `x` and `y`
"""
function CubicInterpolator(x::AbstractVector{<:Real}, y::AbstractVector{<:Real})
    @assert length(x) > 3 "can't do cubic interpolation with <4 points"
    CubicInterpolator(InterpolatorRange(x, y))
end

"""
    CubicInterpolator(f::Function, xa::Real, xb::Real, n::Int)

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
export CubicSplineInterpolator

struct CubicSplineInterpolator
    r::InterpolatorRange
    α::Array{Float64,2} #4xN array of cubic polynomial coefficients
end

"""
    CubicSplineInterpolator(x::AbstractVector{<:Real}, y::AbstractVector{<:Real})

Construct a `CubicSplineInterpolator` for the points in `x` and `y`. This is a natural spline, where the second derivative is set to zero at the boundaries.
"""
function CubicSplineInterpolator(x::AbstractVector{<:Real},
                                 y::AbstractVector{<:Real})
    #construct the underlying range, triggering some checks
    r = InterpolatorRange(x, y)
    #compute coefficients
    #Burden, Richard L., and J. Douglas Faires. Numerical Analysis. 2011.
    a = y[:]
    b = zeros(r.n - 1)
    d = zeros(r.n - 1)
    h = diff(r.x)
    α = zeros(r.n-1)
    for i = 2:r.n-1
        α[i] = 3*(a[i+1] - a[i])/h[i] - 3*(a[i] - a[i-1])/h[i-1]
    end
    c = zeros(r.n)
    l = ones(r.n)
    μ = zeros(r.n)
    z = zeros(r.n)
    l[1] = 1
    for i = 2:r.n-1
        l[i] = 2*(x[i+1] - x[i-1]) - h[i-1]*μ[i-1]
        μ[i] = h[i]/l[i]
        z[i] = (α[i] - h[i-1]*z[i-1])/l[i]
    end
    for j = r.n-1:-1:1
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
    CubicSplineInterpolator(f::Function, xa::Real, xb::Real, n::Int)

Construct a `CubicSplineInterpolator` for the function `f` using `n` evenly spaced function evaluations in the range [`xa`,`xb`]
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
# one-dimensional chebyshev interpolation
export ChebyshevInterpolator

struct ChebyshevInterpolator
    #number of points
    n::Int64
    #lowest value in range
    xa::Float64
    #highest value in range
    xb::Float64
    #interpolation coefficents
    a::Vector{Float64}
    #vector for writing cosine expansion in place
    c::Vector{Float64}
end

"""
    ChebyshevInterpolator(x::AbstractVector{<:Real}, y::AbstractVector{<:Real})

Construct a `ChebyshevInterpolator` for the points in `x` and `y`. The `x` coordinates *must* be arranged on a chebyshev grid, which can be generated using the [`chebygrid`](@ref) function.
"""
function ChebyshevInterpolator(x::AbstractVector{<:Real},
                               y::AbstractVector{<:Real})
    #check for basic issues
    rangecheck(x, y)
    #demand that the input points have chebyshev spacing
    @assert ischebygrid(x) "points must be on a chebyshev grid"
    #number of points/coefficients
    n = length(x)
    #generate expansion coefficients
    a = chebymatrix(n)\y
    #construct
    ChebyshevInterpolator(n, minimum(x), maximum(x), a, ones(n))
end

"""
    ChebyshevInterpolator(f::Function, xa::Real, xb::Real, n::Int)

Construct a `ChebyshevInterpolator` for the function `f` using `n` function evaluations in the range [`xa`,`xb`]. The function evaluations will occur on the chebyshev nodes.
"""
function ChebyshevInterpolator(f::Function, xa::Real, xb::Real, n::Int)
    #set up the range coordinates
    x = chebygrid(xa, xb, n)
    #evaluate the function at those coordinates
    y = f.(x)
    #call the other constructor
    ChebyshevInterpolator(x, y)
end

function (ϕ::ChebyshevInterpolator)(x::Real)::Float64
    #always enforce boundaries
    enforcebounds(x, ϕ.xa, ϕ.xb, true)
    #get the target coordinate in theta space
    θ = x2θ(x, ϕ.xa, ϕ.xb)
    #evaluate the cosine expansion in-place
    fce!(ϕ.c, θ, ϕ.n)
    #then apply the coefficients
    ϕ.a'*ϕ.c
end

#-------------------------------------------------------------------------------
# a base structure for 2d grid interpolation and some general functions

struct InterpolatorGrid
    #number of points along axis 1
    nx::Int64
    #number of points along axis 2
    ny::Int64
    #grid points along axis 1
    x::Vector{Float64}
    #grid points along axis 2
    y::Vector{Float64}
    #lowest value on axis 1
    xa::Float64
    #highest value on axis 1
    xb::Float64
    #lowest value on axis 2
    ya::Float64
    #highest value on axis 2
    yb::Float64
    #values at grid points
    Z::Array{Float64,2}
end

function InterpolatorGrid(x::AbstractVector{<:Real},
                          y::AbstractVector{<:Real},
                          Z::AbstractArray{<:Real,2})
    #make sure types are uniform and make copies along the way
    x = convert.(Float64, x)
    y = convert.(Float64, y)
    Z = convert.(Float64, Z)
    #check for basic grid problems
    gridcheck(x, y, Z)
    #construct the object
    InterpolatorGrid(
        length(x),
        length(y),
        x,
        y,
        minimum(x),
        maximum(x),
        minimum(y),
        maximum(y),
        Z
    )
end

function gridcheck(x::AbstractVector{<:Real},
                   y::AbstractVector{<:Real},
                   Z::AbstractArray{<:Real,2})
    @assert (length(x) > 1) & (length(y) > 1) "grid must have more than one point in each direction"
    @assert size(Z) == (length(x), length(y)) "dimensions mismatched, size(Z) = $(size(Z)) but length(x) = $(length(x)) and length(y) = $(length(y))"
    @assert (all(diff(x) .> 0.0) & all(diff(y) .>= 0.0)) "grid coordinates must be monotonically increasing"
end

function linstruct(T::Type, f::Function,
                   xa::Real, xb::Real, nx::Int,
                   ya::Real, yb::Real, ny::Int)
    x = collect(LinRange(xa, xb, nx))
    y = collect(LinRange(ya, yb, ny))
    X = x .* ones(ny)'
    Y = y' .* ones(nx)
    Z = f.(X, Y)
    T(x, y, Z)
end

#-------------------------------------------------------------------------------
# bilinear interpolator
export BilinearInterpolator

struct BilinearInterpolator
    G::InterpolatorGrid
end

"""
    BilinearInterpolator(x::AbstractVector{<:Real},
                         y::AbstractVector{<:Real},
                         Z::AbstractArray{<:Real,2})

Construct a `BilinearInterpolator` for the grid of points with axis 1 coordinates in `x`, axis 2 coordinates in `y`, and grid values in `Z`.
"""
function BilinearInterpolator(x::AbstractVector{<:Real},
                              y::AbstractVector{<:Real},
                              Z::AbstractArray{<:Real,2})
    BilinearInterpolator(InterpolatorGrid(x, y, Z))
end

"""
    BilinearInterpolator(f::Function,
                         xa::Real, xb::Real, nx::Int,
                         ya::Real, yb::Real, ny::Int)

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
    zxa = (x - xg[i])*(Z[i+1,j]   - Z[i,j]  )/(xg[i+1] - xg[i]) + Z[i,j]
    zxb = (x - xg[i])*(Z[i+1,j+1] - Z[i,j+1])/(xg[i+1] - xg[i]) + Z[i,j+1]
    #then finish by interpolating between the interpolated values along axis 2
    (y - yg[j])*(zxb - zxa)/(yg[j+1] - yg[j]) + zxa
end

#-------------------------------------------------------------------------------
# bicubic interpolator
export BicubicInterpolator

struct BicubicInterpolator
    G::InterpolatorGrid
end

"""
    BicubicInterpolator(x::AbstractVector{<:Real},
                        y::AbstractVector{<:Real},
                        Z::AbstractArray{<:Real,2})

Construct a `BicubicInterpolator` for the grid of points with axis 1 coordinates in `x`, axis 2 coordinates in `y`, and grid values in `Z`.
"""
function BicubicInterpolator(x::AbstractVector{<:Real},
                             y::AbstractVector{<:Real},
                             Z::AbstractArray{<:Real,2})
    #insist on at least 4 points in each dimension
    @assert (length(x) > 3) & (length(y) > 3) "bicubic interpolation requires at least 4 points in each dimension"
    BicubicInterpolator(InterpolatorGrid(x, y, Z))
end

"""
    BicubicInterpolator(f::Function,
                        xa::Real, xb::Real, nx::Int,
                        ya::Real, yb::Real, ny::Int)

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
export BicubicSplineInterpolator

struct BicubicSplineInterpolator
    G::InterpolatorGrid
    α::Array{Array{Float64,2},2}
end

"""
    BicubicSplineInterpolator(x::AbstractVector{<:Real},
                              y::AbstractVector{<:Real},
                              Z::AbstractArray{<:Real,2})

Construct a `BicubicSplineInterpolator` for the grid of points with axis 1 coordinates in `x`, axis 2 coordinates in `y`, and grid values in `Z`.
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
    BicubicSplineInterpolator(f::Function,
                              xa::Real, xb::Real, nx::Int,
                              ya::Real, yb::Real, ny::Int)

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

#-------------------------------------------------------------------------------
# bichebyshev interpolation, a little tricky now!
export BichebyshevInterpolator

struct BichebyshevInterpolator
    #number of points along axis 1
    nx::Int64
    #number of points along axis 2
    ny::Int64
    #lowest value on axis 1
    xa::Float64
    #highest value on axis 1
    xb::Float64
    #lowest value on axis 2
    ya::Float64
    #highest value on axis 2
    yb::Float64
    #matrix and vectors for doing the interpolation
    M::Array{Float64,2} # size ny x nx
    a::Vector{Float64} # length ny for cosine expansion in θy
    b::Vector{Float64} # length nx for cosine expansion in θx
    c::Vector{Float64} # length ny for doing M*b in place
end

"""
    BichebyshevInterpolator(x::AbstractVector{<:Real},
                            y::AbstractVector{<:Real},
                            Z::AbstractArray{<:Real,2})

Construct a `BichebyshevInterpolator` for the grid of points with axis 1 coordinates in `x`, axis 2 coordinates in `y`, and grid values in `Z`. The given points must lie on a chebyshev grid in each direction. These can be generated with the [`chebygrid`](@ref) function or the interpolator can be constructed directly from a function using the method below.
"""
function BichebyshevInterpolator(x::AbstractVector{<:Real},
                                 y::AbstractVector{<:Real},
                                 Z::AbstractArray{<:Real,2})
    #check for basic grid problems
    gridcheck(x, y, Z)
    #grid properties
    nx, ny = length(x), length(y)
    xa, xb = minimum(x), maximum(x)
    ya, yb = minimum(y), maximum(y)
    #reject any non-cheby grid spacing
    @assert ischebygrid(x) "axis 1 coordinates must be on a chebyshev grid"
    @assert ischebygrid(y) "axis 2 coordinates must be on a chebyshev grid"
    #generate interpolation coefficients along axis 1 for each value of axis 2
    α = zeros(nx, ny)
    B = inv(chebymatrix(nx))
    for j = 1:ny
        α[:,j] = B*view(Z,:,j)
    end
    #generate the matrix needed to get coefficents in the other direction
    A = inv(chebymatrix(ny))
    #then combine α and A for efficiency
    M = A*α'
    #other vectors we need
    a = ones(ny)
    b = ones(nx)
    c = zeros(ny)
    #done
    BichebyshevInterpolator(nx, ny, xa, xb, ya, yb, M, a, b, c)
end

"""
    BichebyshevInterpolator(f::Function,
                            xa::Real, xb::Real, nx::Int,
                            ya::Real, yb::Real, ny::Int)

Construct a `BichebyshevInterpolator` for the function `f` using a grid of `nx` points on the first axis in [`xa`,`xb`] and `ny` points on the second axis in [`ya`,`yb`].
"""
function BichebyshevInterpolator(f::Function,
                                 xa::Real, xb::Real, nx::Int,
                                 ya::Real, yb::Real, ny::Int)
    #set up the grid
    X, Y = chebygrid(xa, xb, nx, ya, yb, ny)
    #evaluate the function at chebyshev grid points
    Z = f.(X, Y)
    #call the other constructor
    BichebyshevInterpolator(X[:,1], Y[1,:], Z)
end

function (Φ::BichebyshevInterpolator)(x::Real, y::Real)::Float64
    #always enforce boundaries
    enforcebounds(x, Φ.xa, Φ.xb, y, Φ.ya, Φ.yb, true)
    #get the target point in theta space
    θy = x2θ(y, Φ.ya, Φ.yb)
    θx = x2θ(x, Φ.xa, Φ.xb)
    #evaluate as few cosines directly as possible with no new allocations
    fce!(Φ.a, θy, Φ.ny)
    fce!(Φ.b, θx, Φ.nx)
    #perform M*b, which interpolates along the first axis, also in-place
    mul!(Φ.c, Φ.M, Φ.b)
    #then a'*c interpolates along the second axis
    Φ.a'*Φ.c
end

end
