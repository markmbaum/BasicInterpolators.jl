#-------------------------------------------------------------------------------
# abstract classes for interpolators

export BasicInterpolator
export OneDimensionalInterpolator
export TwoDimensionalInterpolator

abstract type BasicInterpolator end
abstract type OneDimensionalInterpolator <: BasicInterpolator end
abstract type TwoDimensionalInterpolator <: BasicInterpolator end

function Base.show(io::IO, ϕ::P) where {P<:OneDimensionalInterpolator}
    print(io, "$P ∈ [$(ϕ.r.xa), $(ϕ.r.xb)]")
end

function Base.show(io::IO, ϕ::P) where {P<:TwoDimensionalInterpolator}
    print(io, "$P ∈ [$(ϕ.G.xa), $(ϕ.G.xb)], [$(ϕ.G.ya), $(ϕ.G.yb)]")
end

#-------------------------------------------------------------------------------
# cell finding stuff

export findcell

"""
    findcell(q, V, n)

Use bisection search to find the cell containing `q`, assuming `V` is a sorted vector of `n` coordinates. The returned integer is the index of the element in `V` immediately less than `q`. For example, if `findcell` returns 2, then `q ∈ [V[2],V[3])`. If `q` is less than every element in `V`, 1 is returned, indicating the first cell in `V`. If `q` is greater than every element in `V`, `length(V)-1` is returned, indicating the last cell in `V`.
"""
function findcell(q, V, n::Int64)::Int64
    #handle boundaries
    @inbounds (q <= V[1]) && return(1)
    @inbounds (q >= V[n]) && return(n - 1)
    #bisection search for the containing cell
    L = 0
    H = n
    while H - L > 1
        M = (H + L) ÷ 2
        @inbounds if V[M] > q
            H = M
        else
            L = M
        end
    end
    return L
end

findcell(x, ϕ::OneDimensionalInterpolator)::Int64 = findcell(x, ϕ.r.x, ϕ.r.n)

function findcell(x, y, Φ::TwoDimensionalInterpolator)::NTuple{2,Int64}
    i = findcell(x, Φ.G.x, Φ.G.nx)
    j = findcell(y, Φ.G.y, Φ.G.ny)
    return i, j
end

#-------------------------------------------------------------------------------
# types for handling boundaries

export AbstractBoundaries, NoBoundaries, WeakBoundaries, StrictBoundaries

upperbounderror(x, xb, axis::String) = error("Interpolation location $x outside upper interpolation limit $xb on the $axis axis")

lowerbounderror(x, xa, axis::String) = error("Interpolation location $x outside lower interpolation limit $xa on the $axis axis")

abstract type AbstractBoundaries end

#--------------------------------

"""
    NoBoundaries()

Performs no boundary checking. Allows interpolators to blindly extrapolate when possible.
"""
struct NoBoundaries <: AbstractBoundaries end

#totally empty function
function (B::NoBoundaries)(x...) end

#--------------------------------

"""
    WeakBoundaries()

Allows small overshoots at boundaries, but not large ones. Errors are only triggered when the interpolation coordinate is outside of a boundary and not close to it. At the lower boundary, for example, an error would be triggered when `(x < boundary) & !(x ≈ boundary)`.
"""
struct WeakBoundaries <: AbstractBoundaries end

function (B::WeakBoundaries)(x, xa, xb)::Nothing
    ((x < xa) && !(x ≈ xa)) && lowerbounderror(x, xa, "first")
    ((x > xb) && !(x ≈ xb)) && upperbounderror(x, xb, "first")
    nothing
end

function (B::WeakBoundaries)(x, xa, xb, y, ya, yb)::Nothing
    ((x < xa) && !(x ≈ xa)) && lowerbounderror(x, xa, "first")
    ((x > xb) && !(x ≈ xb)) && upperbounderror(x, xb, "first")
    ((y < ya) && !(y ≈ ya)) && lowerbounderror(y, ya, "second")
    ((y > yb) && !(y ≈ yb)) && upperbounderror(y, yb, "second")
    nothing
end

#--------------------------------

"""
    StrictBoundaries()

Triggers errors whenever interpolation coordinates are outside of the boundaries.
"""
struct StrictBoundaries <: AbstractBoundaries end

function (B::StrictBoundaries)(x, xa, xb)::Nothing
    (x < xa) && lowerbounderror(x, xa, "first")
    (x > xb) && upperbounderror(x, xb, "first")
    nothing
end

function (B::StrictBoundaries)(x, xa, xb, y, ya, yb)::Nothing
    (x < xa) && lowerbounderror(x, xa, "first")
    (x > xb) && upperbounderror(x, xb, "first")
    (y < ya) && lowerbounderror(y, ya, "second")
    (y > yb) && upperbounderror(y, yb, "second")
    nothing
end

#-------------------------------------------------------------------------------
# a base structure for 1d interpolations and some related functions

struct InterpolatorRange{T}
    #number of points
    n::Int64
    #grid/sample points
    x::Vector{T}
    #lowest sample value
    xa::T
    #highest sample value
    xb::T
    #values at sample points
    y::Vector{T}
end

function InterpolatorRange(x::AbstractVector, y::AbstractVector)
    #collect and copy
    T = promote_type(eltype(x), eltype(y))
    x = collect(T, x)
    y = collect(T, y)
    #check for basic problems
    rangecheck(x, y)
    #grid properties
    n = length(x)
    xa = minimum(x)
    xb = maximum(x)
    #construct
    InterpolatorRange{T}(n, x, xa, xb, y)
end

function rangecheck(x::AbstractVector, y::AbstractVector, minpts=2)
    nx, ny = length(x), length(y)
    @assert nx == ny "length of x ($nx) does not match length of y ($ny)"
    @assert nx >= minpts "must have at least $minpts points in the interpolation range"
    @assert issorted(x) "grid points must be in strictly ascending order"
end

function linstruct(I::Type, f::F, xa, xb, n::Int, boundaries::AbstractBoundaries) where {F}
    @assert n > 1 "cannot construct a length<=1 interpolation range"
    x = LinRange(xa, xb, n)
    y = f.(x)
    I(x, y, boundaries)
end

#-------------------------------------------------------------------------------
# a base structure for 2d grid interpolation and some related functions

struct InterpolatorGrid{T}
    #number of points along axis 1
    nx::Int64
    #number of points along axis 2
    ny::Int64
    #grid points along axis 1
    x::Vector{T}
    #grid points along axis 2
    y::Vector{T}
    #lowest value on axis 1
    xa::T
    #highest value on axis 1
    xb::T
    #lowest value on axis 2
    ya::T
    #highest value on axis 2
    yb::T
    #values at grid points
    Z::Array{T,2}
end

function InterpolatorGrid(x::AbstractVector, y::AbstractVector, Z::AbstractMatrix)
    #collect and copy
    T = promote_type(eltype(x), eltype(y), eltype(Z))
    x = collect(T, x)
    y = collect(T, y)
    Z = collect(T, Z)
    #check for basic grid problems
    gridcheck(x, y, Z)
    #boundaries
    xa, xb = minimum(x), maximum(x)
    ya, yb = minimum(y), maximum(y)
    #construct the object
    InterpolatorGrid{T}(length(x), length(y), x, y, xa, xb, ya, yb, Z)
end

function gridcheck(x::AbstractVector, y::AbstractVector, Z::AbstractMatrix, minpts=2)
    @assert (length(x) >= minpts) & (length(y) >= minpts) "must have at least $minpts points in each dimension of interpolation grid"
    @assert size(Z) == (length(x), length(y)) "dimensions mismatched, size(Z) = $(size(Z)) but length(x) = $(length(x)) and length(y) = $(length(y))"
    @assert (issorted(x) & issorted(y)) "grid coordinates must be monotonically increasing without duplicates"
end

function linstruct(I::Type, f::F,
                   xa, xb, nx::Int,
                   ya, yb, ny::Int,
                   boundaries::AbstractBoundaries) where {F}
    @assert (nx > 1) & (ny > 1) "cannot use a length<=1 grid axis"
    xa, xb, ya, yb = promote(xa, xb, ya, yb)
    T = typeof(xa)
    x = LinRange(xa, xb, nx)
    y = LinRange(ya, yb, ny)
    X = x .* ones(T, ny)'
    Y = y' .* ones(T, nx)
    Z = f.(X, Y)
    I(x, y, Z, boundaries)
end
