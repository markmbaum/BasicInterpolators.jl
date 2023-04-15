export linear, quadratic, quaddiff, cubic,
       neville, vandermonde, vandermonde!, cubichermite,
       LinearInterpolator, CubicInterpolator,
       BilinearInterpolator, BicubicInterpolator


"""
    linear(x, xₚ, yₚ)

Perform simple linear interpolation of the points defined by coordinates `xₚ` and values `yₚ`, at the coordinate `x`. `xₚ` and `yₚ` must both contain two points.
"""
function linear(x, xₚ, yₚ)
    @assert length(xₚ) == length(yₚ) == 2 "Two (2) points/coordinates are required for linear interpolation. For more than two points, use a [`LinearInterpolator`](@ref)."
    @inbounds x₁, x₂ = xₚ[1], xₚ[2]
    @inbounds y₁, y₂ = yₚ[1], yₚ[2]
    return y₁ + (x - x₁)*(y₂ - y₁)/(x₂ - x₁)
end

"""
    quadratic(x, xₚ, yₚ)

Perform quadratic polynomial interpolation of the points defined by coordinates `xₚ` and values `yₚ`, at the coordinate `x`, using Neville's algorithm. `xₚ` and `yₚ` must both contain three points.
"""
function quadratic(x, xₚ, yₚ)
    @assert length(xₚ) == length(yₚ) == 3 "Three (3) points/coordinates are required for simple quadratic interpolation."
    #short names
    @inbounds x₁, x₂, x₃ = xₚ[1], xₚ[2], xₚ[3]
    @inbounds y₁, y₂, y₃ = yₚ[1], yₚ[2], yₚ[3]
    #first stage
    p₁₂ = ((x - x₂)*y₁ + (x₁ - x)*y₂)/(x₁ - x₂)
    p₂₃ = ((x - x₃)*y₂ + (x₂ - x)*y₃)/(x₂ - x₃)
    #final stage
    ((x - x₃)*p₁₂ + (x₁ - x)*p₂₃)/(x₁ - x₃)
end

"""
    cubic(x, xₚ, yₚ)

Perform cubic polynomial interpolation of the points defined by coordinates `xₚ` and values `yₚ`, at the coordinate `x`, using Neville's algorithm. `xₚ` and `yₚ` must both contain four points.
"""
function cubic(x, xₚ, yₚ)
    @assert length(xₚ) == length(yₚ) == 4 "Four (4) points/coordinates are required for simple cubic interpolation/"
    #short names
    @inbounds x₁, x₂, x₃, x₄ = xₚ[1], xₚ[2], xₚ[3], xₚ[4]
    @inbounds y₁, y₂, y₃, y₄ = yₚ[1], yₚ[2], yₚ[3], yₚ[4]
    #first stage
    p₁₂ = ((x - x₂)*y₁ + (x₁ - x)*y₂)/(x₁ - x₂)
    p₂₃ = ((x - x₃)*y₂ + (x₂ - x)*y₃)/(x₂ - x₃)
    p₃₄ = ((x - x₄)*y₃ + (x₃ - x)*y₄)/(x₃ - x₄)
    #second stage
    p₁₂₃ = ((x - x₃)*p₁₂ + (x₁ - x)*p₂₃)/(x₁ - x₃)
    p₂₃₄ = ((x - x₄)*p₂₃ + (x₂ - x)*p₃₄)/(x₂ - x₄)
    #final stage
    ((x - x₄)*p₁₂₃ + (x₁ - x)*p₂₃₄)/(x₁ - x₄)
end

"""
    neville(x, xₚ, yₚ)

Perform polynomial interpolation of the points defined by coordinates `xₚ` and values `yₚ`, at the coordinate `x`, using Neville's algorithm with as many points as are provided. `xₚ` and `yₚ` must have the same length. With only 3 or 4 points the [`quadratic`](@ref) and [`cubic`](@ref) functions will be considerably faster.
"""
function neville(x, xₚ, yₚ)
    @assert length(xₚ) == length(yₚ) "Can't perform Neville's algorithm with coordinate vectors of different lengths. length(xₚ) = $(length(xₚ)) but length(yₚ) == $(length(yₚ))"
    n = length(xₚ)
    P = zeros(n, n)
    P[:,1] = yₚ
    for i ∈ 2:n, j ∈ i:n
        @inbounds P[j,i] = ((x - xₚ[j-i+1])*P[j,i-1] - (x - xₚ[j])*P[j-1,i-1])/(xₚ[j] - xₚ[j-i+1])
    end
    return P[n,n]
end

"""
    vandermonde(x, y)

Generate the coefficients of an arbitrary order polynomial passing through the ponts defined by coordinates `x` and value `y`. For n points, n coefficients ``[c_0, c_1, ..., c_{n-1}]`` are returned forming the polynomial ``c_0 + c_1x + ... + c_{n-1}x^{n-1}``

!!! warning

    Solving for the the coefficients of a high-order polynomial is a notoriously ill-conditioned problem. It is not recommended for orders greater than 5 or 6, although it depends on the application. If you must interpolate with a high-order polynomial, it's better to use the [`neville`](@ref) function instead of computing coefficients.
"""
function vandermonde(x, y)
    @assert length(x) == length(y) "length of x must equal length of y"
    c = zeros(length(x))
    vandermonde!(c, x, y)
    return c
end

function vandermonde!(c, x, y)
    @assert length(c) == length(x) == length(y) "vandermonde! requires length(c) == length(x) == length(y)"
    n = length(x)
    s = similar(c)
    c[:] .= zero(eltype(c))
    s[n] = -x[1]
    @inbounds for i = 2:n
        for j = n-i+1:n-1
            s[j] -= x[i]*s[j+1]
        end
        s[n] -= x[i]
    end
    @inbounds for j = 1:n
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
    cubichermite(x, x₁, x₂, y₁, y₂, y₁′, y₂′)

Interpolate a cubic polynomial between two points, given its values and first derivatives at the points.
"""
function cubichermite(x, x₁, x₂, y₁, y₂, y₁′, y₂′)
    Δx = (x₂ - x₁)
    ξ = (x - x₁)/Δx
    u = Δx*y₁′
    v = Δx*y₂′
    ξ^3*(2*y₁ + u - 2*y₂ + v) + ξ^2*(-3*y₁ - 2*u + 3*y₂ - v) + ξ*u + y₁
end

#-------------------------------------------------------------------------------
# simple piecewise linear interpolator

struct LinearInterpolator{T,B} <: OneDimensionalInterpolator
    r::InterpolatorRange{T}
    boundaries::B
end

"""
    LinearInterpolator(x, y, boundaries=StrictBoundaries())

Construct a `LinearInterpolator` for the points defined by coordinates `x` and values `y`
"""
function LinearInterpolator(x, y, boundaries::AbstractBoundaries=StrictBoundaries())
    LinearInterpolator(InterpolatorRange(x, y), boundaries)
end

"""
    LinearInterpolator(f, xa, xb, n, boundaries=StrictBoundaries())

Construct a `LinearInterpolator` for the function `f` using `n` evenly spaced function evaluations in the range [`xa`,`xb`]
"""
function LinearInterpolator(f, xa, xb, n::Int, boundaries::AbstractBoundaries=StrictBoundaries())
    linstruct(LinearInterpolator, f, xa, xb, n, boundaries)
end

function (ϕ::LinearInterpolator)(x)
    #enforce boundaries
    ϕ.boundaries(x, ϕ.r.xa, ϕ.r.xb)
    #find the interpolation cell
    i = findcell(x, ϕ)
    #short names
    @inbounds x₁, x₂ = ϕ.r.x[i], ϕ.r.x[i+1]
    @inbounds y₁, y₂ = ϕ.r.y[i], ϕ.r.y[i+1]
    #interpolate
    (x - x₁)*(y₂ - y₁)/(x₂ - x₁) + y₁
end

Base.getindex(ϕ::LinearInterpolator, i) = ϕ.r.y[i]
Base.firstindex(::LinearInterpolator) = 1
Base.lastindex(ϕ::LinearInterpolator) = ϕ.r.n
function Base.setindex!(ϕ::LinearInterpolator, v, i)
    ϕ.r.y[i] = v
end
function Base.copy(ϕ::LinearInterpolator)
    LinearInterpolator(ϕ.r, ϕ.boundaries)
end

#-------------------------------------------------------------------------------
# piecewise cubic interpolator without continuous derivatives (not splines)

struct CubicInterpolator{T,B} <: OneDimensionalInterpolator
    r::InterpolatorRange{T}
    boundaries::B
end

"""
    CubicInterpolator(x, y, boundaries=StrictBoundaries())

Construct a `CubicInterpolator` for the points defined by coordinates `x` and values `y`
"""
function CubicInterpolator(x, y, boundaries::AbstractBoundaries=StrictBoundaries())
    @assert length(x) > 3 "can't do cubic interpolation with < 4 points"
    CubicInterpolator(InterpolatorRange(x, y), boundaries)
end

"""
    CubicInterpolator(f, xa, xb, n, boundaries=StrictBoundaries())

Construct a `CubicInterpolator` for the function `f` using `n` evenly spaced function evaluations in the range [`xa`,`xb`]
"""
function CubicInterpolator(f, xa, xb, n::Int, boundaries::AbstractBoundaries=StrictBoundaries())
    linstruct(CubicInterpolator, f, xa, xb, n, boundaries)
end

function (ϕ::CubicInterpolator)(x)
    #enforce boundaries if desired
    ϕ.boundaries(x, ϕ.r.xa, ϕ.r.xb)
    #find the interpolation point
    i = findcell(x, ϕ)
    #determine which points to neville
    if i == 1
        I = 1:4
    elseif i == ϕ.r.n - 1
        I = ϕ.r.n-3:ϕ.r.n
    else
        I = i-1:i+2
    end
    #interpolate
    @inbounds cubic(x, view(ϕ.r.x,I), view(ϕ.r.y,I))
end

Base.getindex(ϕ::CubicInterpolator, i) = ϕ.r.y[i]
Base.firstindex(::CubicInterpolator) = 1
Base.lastindex(ϕ::CubicInterpolator) = ϕ.r.n
function Base.setindex!(ϕ::CubicInterpolator, v, i)
    ϕ.r.y[i] = v
end
function Base.copy(ϕ::CubicInterpolator)
    CubicInterpolator(ϕ.r, ϕ.boundaries)
end

#-------------------------------------------------------------------------------
# bilinear interpolator

struct BilinearInterpolator{T,B} <: TwoDimensionalInterpolator
    G::InterpolatorGrid{T}
    boundaries::B
end

"""
    BilinearInterpolator(x, y, Z, boundaries=StrictBoundaries())

Construct a `BilinearInterpolator` for the grid of points points defined by coordinates `x`,`y` and values `Z`.
"""
function BilinearInterpolator(x, y, Z, boundaries::AbstractBoundaries=StrictBoundaries())
    BilinearInterpolator(InterpolatorGrid(x, y, Z), boundaries)
end

"""
    BilinearInterpolator(f, xa, xb, nx, ya, yb, ny, boundaries=StrictBoundaries())

Construct a `BilinearInterpolator` for the function `f` using a grid of `nx` points evenly spaced on the first axis in [`xa`,`xb`] and `ny` points evenly spaced on the second axis in [`ya`,`yb`].
"""
function BilinearInterpolator(f, xa, xb, nx::Int, ya, yb, ny::Int, boundaries::AbstractBoundaries=StrictBoundaries())
    linstruct(BilinearInterpolator, f, xa, xb, nx, ya, yb, ny, boundaries)
end

function (Φ::BilinearInterpolator)(x, y)
    #enforce boundaries if desired
    Φ.boundaries(x, Φ.G.xa, Φ.G.xb, y, Φ.G.ya, Φ.G.yb)
    #find the proper grid box to interpolate inside
    i, j = findcell(x, y, Φ)
    #clear names
    @inbounds x₁, x₂ = Φ.G.x[i], Φ.G.x[i+1]
    @inbounds y₁, y₂ = Φ.G.y[j], Φ.G.y[j+1]
    @inbounds Z₁₁, Z₂₁ = Φ.G.Z[i,j],     Φ.G.Z[i+1,j]
    @inbounds Z₂₂, Z₁₂ = Φ.G.Z[i+1,j+1], Φ.G.Z[i,j+1]
    #reused differences
    xx₁ = x - x₁
    x₂x₁ = x₂ - x₁
    #interpolate along axis 1 first
    t = xx₁*(Z₂₁ - Z₁₁)/x₂x₁ + Z₁₁
    u = xx₁*(Z₂₂ - Z₁₂)/x₂x₁ + Z₁₂
    #then finish by interpolating between the interpolated values along axis 2
    (y - y₁)*(u - t)/(y₂ - y₁) + t
end

Base.getindex(Φ::BilinearInterpolator, i, j) = Φ.G.Z[i,j]
function Base.setindex!(Φ::BilinearInterpolator, v, i, j)
    Φ.G.Z[i,j] = v
end
function Base.copy(Φ::BilinearInterpolator)
    BilinearInterpolator(Φ.G, Φ.boundaries)
end

#-------------------------------------------------------------------------------
# bicubic interpolator

struct BicubicInterpolator{T,B} <: TwoDimensionalInterpolator
    G::InterpolatorGrid{T}
    boundaries::B
end

"""
    BicubicInterpolator(x, y, Z, boundaries=StrictBoundaries())

Construct a `BicubicInterpolator` for the grid of points points defined by coordinates `x`,`y` and values `Z`.
"""
function BicubicInterpolator(x, y, Z, boundaries::AbstractBoundaries=StrictBoundaries())
    #insist on at least 4 points in each dimension
    @assert (length(x) > 3) & (length(y) > 3) "bicubic interpolation requires at least 4 points in each dimension"
    BicubicInterpolator(InterpolatorGrid(x, y, Z), boundaries)
end

"""
    BicubicInterpolator(f, xa, xb, nx, ya, yb, ny, boundaries=StrictBoundaries())

Construct a `BicubicInterpolator` for the function `f` using a grid of `nx` points evenly spaced on the first axis in [`xa`,`xb`] and `ny` points evenly spaced on the second axis in [`ya`,`yb`].
"""
function BicubicInterpolator(f, xa, xb, nx::Int, ya, yb, ny::Int, boundaries::AbstractBoundaries=StrictBoundaries())
    linstruct(BicubicInterpolator, f, xa, xb, nx, ya, yb, ny, boundaries)
end

function (Φ::BicubicInterpolator)(x, y)
    #enforce boundaries if desired
    Φ.boundaries(x, Φ.G.xa, Φ.G.xb, y, Φ.G.ya, Φ.G.yb)
    #find the proper grid box to interpolate inside
    i, j = findcell(x, y, Φ)
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
    #view the proper chunks of the arrays
    xᵣ = view(Φ.G.x, I)
    yᵣ = view(Φ.G.y, J)
    #perform initial 4 interpolations
    zₓ = (
        cubic(x, xᵣ, @view Φ.G.Z[I,J[1]]),
        cubic(x, xᵣ, @view Φ.G.Z[I,J[2]]),
        cubic(x, xᵣ, @view Φ.G.Z[I,J[3]]),
        cubic(x, xᵣ, @view Φ.G.Z[I,J[4]])
    )
    #final interpolation
    cubic(y, yᵣ, zₓ)
end

Base.getindex(Φ::BicubicInterpolator, i, j) = Φ.G.Z[i,j]
function Base.setindex!(Φ::BicubicInterpolator, v, i, j)
    Φ.G.Z[i,j] = v
end
function Base.copy(Φ::BicubicInterpolator)
    BicubicInterpolator(Φ.G, Φ.boundaries)
end