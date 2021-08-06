export chebygrid, ChebyshevInterpolator, BichebyshevInterpolator

ξ2x(ξ::Real, a::Real, b::Real) = (ξ + 1)*((b - a)/2) + a

x2ξ(x::Real, a::Real, b::Real) = 2*(x - a)/(b - a) - 1

x2θ(x::Real, a::Real, b::Real) = acos(x2ξ(x, a, b))

"""
    chebygrid(n)

Create an array of `n` chebyshev nodes in [-1,1]
"""
function chebygrid(n::Int)::Vector{Float64}
    cos.(π*(n-1:-1:0)/(n-1))
end

"""
    chebygrid(xa, xb, n)

Create an array of `n` chebyshev nodes in [`xa`,`xb`]
"""
function chebygrid(xa::Real, xb::Real, n::Int)::Vector{Float64}
    ξ2x.(chebygrid(n), xa, xb)
end

"""
    chebygrid(xa, xb, nx, ya, yb, ny)

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
    for j ∈ 1:n, k ∈ 1:n
        A[j,k] = cheby(ξ[j], k-1)
    end
    return A
end

function Tₖ!(T::Vector{Float64}, ξ::Float64, n::Int64)::Nothing
    Tₖ₋₂::Float64 = 1.0 
    Tₖ₋₁::Float64 = ξ
    @inbounds T[2] = Tₖ₋₁
    for k = 3:n
        #compute next value
        Tₖ::Float64 = 2*ξ*Tₖ₋₁ - Tₖ₋₂
        #set array value
        @inbounds T[k] = Tₖ
        #swaps
        Tₖ₋₂ = Tₖ₋₁
        Tₖ₋₁ = Tₖ
    end
    nothing
end

#-------------------------------------------------------------------------------
# cache for inverted cheby matrices needed for interpolator setup

#a cache for the inverted matrices needed for Bicheby setup
const invchebmat = Dict{Int64, Array{Float64,2}}()

function invertedchebymatrix(n::Int64)::Array{Float64,2}
    if !haskey(invchebmat, n)
        invchebmat[n] = inv(chebymatrix(n))
    end
    return invchebmat[n]
end

#-------------------------------------------------------------------------------
# one-dimensional interpolation

struct ChebyshevInterpolator{N} <: OneDimensionalInterpolator
    #number of points
    n::Int64
    #lowest value in range
    xa::Float64
    #highest value in range
    xb::Float64
    #interpolation coefficents
    a::NTuple{N,Float64}
    #vector for writing cosine expansion in place
    c::Vector{Float64}
    #must always have strict boundaries
    boundaries::StrictBoundaries
end

"""
    ChebyshevInterpolator(x, y)

Construct a `ChebyshevInterpolator` for the points defined by coordinates `x` and values `y`. The `x` coordinates *must* be arranged on a chebyshev grid, which can be generated using the [`chebygrid`](@ref) function.

!!! warning

    The Chebyshev interpolator is *not thread-safe*. It computes a cosine expansion in-place using an array stored with the object. A single `ChebyshevInterpolator` should never be called by multiple threads at once.
"""
function ChebyshevInterpolator(x::AbstractVector{<:Real},
                               y::AbstractVector{<:Real})
    #check for basic issues
    rangecheck(x, y)
    #demand that the input points have chebyshev spacing
    @assert ischebygrid(x) "points must be on a chebyshev grid"
    #number of points/coefficients
    n = length(x)
    #get the inverted cheby matrix, from cache or fresh
    A = invertedchebymatrix(n)
    #generate expansion coefficients
    a = Tuple(A*y)
    #construct
    ChebyshevInterpolator(n, minimum(x), maximum(x), a, ones(n), StrictBoundaries())
end

"""
    ChebyshevInterpolator(f, xa, xb, n)

Construct a `ChebyshevInterpolator` for the function `f` using `n` function evaluations in the range [`xa`,`xb`]. The function evaluations will occur on the chebyshev nodes.
"""
function ChebyshevInterpolator(f::F, xa::Real, xb::Real, n::Int) where {F}
    #set up the range coordinates
    x = chebygrid(xa, xb, n)
    #evaluate the function at those coordinates
    y = f.(x)
    #call the other constructor
    ChebyshevInterpolator(x, y)
end

function (ϕ::ChebyshevInterpolator)(x::Real)::Float64
    #always enforce boundaries
    ϕ.boundaries(x, ϕ.xa, ϕ.xb)
    #evaluate the cosine expansion in-place
    Tₖ!(ϕ.c, x2ξ(x, ϕ.xa, ϕ.xb), ϕ.n)
    #then apply the coefficients
    dot(ϕ.a, ϕ.c)
end

#-------------------------------------------------------------------------------
# bichebyshev interpolation, a little trickier now!

struct BichebyshevInterpolator <: TwoDimensionalInterpolator
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
    M::Matrix{Float64} # size ny x nx
    a::Vector{Float64} # length ny for cosine expansion in θy
    b::Vector{Float64} # length nx for cosine expansion in θx
    c::Vector{Float64} # length ny for doing M*b in place
    #must always use strict boundaries
    boundaries::StrictBoundaries
end

"""
    BichebyshevInterpolator(x, y, Z)

Construct a `BichebyshevInterpolator` for the grid of points defined by coordinates (`x`,`y`) and values `Z`. The given points must lie on a chebyshev grid in each direction. These can be generated with the [`chebygrid`](@ref) function or the interpolator can be constructed directly from a function using the method below.

!!! warning

    The Bichebyshev interpolator is *not thread-safe*. It computes a cosine expansion and does some linear algebra in-place using arrays stored with the object. A single `BichebyshevInterpolator` should never be called by multiple threads at once.
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
    #get inverted matrices from cache or generate them
    B = invertedchebymatrix(nx)
    A = invertedchebymatrix(ny)
    #generate interpolation coefficients along axis 1 for each value of axis 2
    α = zeros(nx, ny)
    for j = 1:ny
        mul!(view(α,:,j), B, view(Z,:,j))
    end
    #then combine α and A
    M = A*α'
    #other vectors we need for doing the actual interpolation
    a = ones(ny)
    b = ones(nx)
    c = zeros(ny)
    #done
    BichebyshevInterpolator(nx, ny, xa, xb, ya, yb, M, a, b, c, StrictBoundaries())
end

"""
    BichebyshevInterpolator(f, xa, xb, nx, ya, yb, ny)

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
    Φ.boundaries(x, Φ.xa, Φ.xb, y, Φ.ya, Φ.yb)
    #evaluate Chebyshev polys at the coordinates recursively and in-place
    Tₖ!(Φ.a, x2ξ(y, Φ.ya, Φ.yb), Φ.ny)
    Tₖ!(Φ.b, x2ξ(x, Φ.xa, Φ.xb), Φ.nx)
    #perform M*b, which interpolates along the first axis, also in-place
    mul!(Φ.c, Φ.M, Φ.b)
    #then a'*c interpolates along the second axis
    dot(Φ.a, Φ.c)
end
 