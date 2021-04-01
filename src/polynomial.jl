export quadratic, cubic, neville, vandermonde, vandermonde!, cubichermite,
       LinearInterpolator, CubicInterpolator,
       BilinearInterpolator, BicubicInterpolator

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
    c = zeros(length(x))
    vandermonde!(c, x, y)
    return c
end

function vandermonde!(c::AbstractVector{<:Real},
                      x::AbstractVector{<:Real},
                      y::AbstractVector{<:Real})
    @assert length(c) == length(x) == length(y) "vandermonde! requires length(c) == length(x) == length(y)"
    n = length(x)
    s = zeros(n)
    c[:] .= 0.0
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
end

"""
    cubichermite(x, xa, xb, ya, yb, dya, dyb)

Interpolate a cubic polynomial between two points, given its values and first derivatives at the points.
"""
function cubichermite(x::Real,
                      xa::Real, xb::Real,
                      ya::Real, yb::Real,
                      dya::Real, dyb::Real)::Float64
    Δx = (xb - xa)
    ξ = (x - xa)/Δx
    u = Δx*dya
    v = Δx*dyb
    ξ^3*(2*ya + u - 2*yb + v) + ξ^2*(-3*ya - 2*u + 3*yb - v) + ξ*u + ya
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
# bilinear interpolator

struct BilinearInterpolator
    G::InterpolatorGrid
end

"""
    BilinearInterpolator(x, y, Z)

Construct a `BilinearInterpolator` for the grid of points points defined by coordinates `x`,`y` and values `Z`.
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

Construct a `BicubicInterpolator` for the grid of points points defined by coordinates `x`,`y` and values `Z`.
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
