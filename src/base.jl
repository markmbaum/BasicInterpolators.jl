export findcell, NoBoundaries, WeakBoundaries, StrictBoundaries, values

#-------------------------------------------------------------------------------
# abstract classes for interpolators

export OneDimensionalInterpolator, TwoDimensionalInterpolator

abstract type OneDimensionalInterpolator end
abstract type TwoDimensionalInterpolator end

#-------------------------------------------------------------------------------
# cell finding stuff

"""
    findcell(q, V, n)

Use bisection search to find the cell containing `q`, assuming `V` is a sorted vector of `n` coordinates. The returned integer is the index of the element in `V` immediately less than `q`. For example, if `findcell` returns 2, then `q ∈ [V[2],V[3])`. If `q` is less than every element in `V`, 1 is returned, indicating the first cell in `V`. If `q` is greater than every element in `V`, `length(V)-1` is returned, indicating the last cell in `V`.
"""
function findcell(q::Real, V::AbstractVector, n::Int64)::Int64
    #handle boundaries
    @inbounds (q <= V[1]) && return(1)
    @inbounds (q >= V[n]) && return(n - 1)
    #bisection search for the containing cell
    ilo = 0
    ihi = n
    while ihi - ilo > 1
        imid = (ihi + ilo) ÷ 2
        @inbounds if V[imid] > q
            ihi = imid
        else
            ilo = imid
        end
    end
    return ilo
end

function findcell(x::Real, ϕ::OneDimensionalInterpolator)::Int64
    #check previous index used
    i::Int64 = ϕ.i[]
    @inbounds (x >= ϕ.r.x[i]) && (x <= ϕ.r.x[i+1]) && return(i)
    #didn't work, find a fresh cell
    i = findcell(x, ϕ.r.x, ϕ.r.n)
    #store the index
    ϕ.i[] = i
    return i
end

function findcell(x::Real, y::Real, Φ::TwoDimensionalInterpolator)::NTuple{2,Int64}
    #check previous indices used
    i::Int64 = Φ.i[]
    j::Int64 = Φ.j[]
    @inbounds if (x >= Φ.G.x[i]) && (x <= Φ.G.x[i+1])
        @inbounds if (y >= Φ.G.y[j]) & (y <= Φ.G.y[j+1])
            return i,j
        end
    end
    #didn't work, find a fresh cell
    i = findcell(x, Φ.G.x, Φ.G.nx)
    j = findcell(y, Φ.G.y, Φ.G.ny)
    #store the indices
    Φ.i[] = i
    Φ.j[] = j
    return i, j
end

#-------------------------------------------------------------------------------
# types for handling boundaries

function upperbounderror(x, xb, axis::String)
    error("Interpolation location $x outside upper interpolation limit $xb on the $axis axis")
end

function lowerbounderror(x, xa, axis::String)
    error("Interpolation location $x outside lower interpolation limit $xa on the $axis axis")
end

abstract type AbstractBoundaries end

#--------------------------------

"""
    NoBoundaries()

Allows interpolators to blindly extrapolate if possible
"""
struct NoBoundaries <: AbstractBoundaries end

function (B::NoBoundaries)(x...) end

#--------------------------------

"""
    WeakBoundaries()

Allows small overshoots at boundaries, but not large ones. Errors are only triggered when the interpolation coordinate is outside of a boundary and not close to it: `(x < boundary) & !(x ≈ boundary)`.
"""
struct WeakBoundaries <: AbstractBoundaries end

function (B::WeakBoundaries)(x, xa, xb)
    ((x < xa) && !(x ≈ xa)) && lowerbounderror(x, xa, "first")
    ((x > xb) && !(x ≈ xb)) && upperbounderror(x, xb, "first")
end

function (B::WeakBoundaries)(x, xa, xb, y, ya, yb)
    ((x < xa) && !(x ≈ xa)) && lowerbounderror(x, xa, "first")
    ((x > xb) && !(x ≈ xb)) && upperbounderror(x, xb, "first")
    ((y < ya) && !(y ≈ ya)) && lowerbounderror(y, ya, "second")
    ((y > yb) && !(y ≈ yb)) && upperbounderror(y, yb, "second")
end

#--------------------------------

"""
    StrictBoundaries()

Triggers errors whenever interpolation coordinates are outside of the boundaries.
"""
struct StrictBoundaries <: AbstractBoundaries end

function (B::StrictBoundaries)(x, xa, xb)
    (x < xa) && lowerbounderror(x, xa, "first")
    (x > xb) && upperbounderror(x, xb, "first")
end

function (B::StrictBoundaries)(x, xa, xb, y, ya, yb)
    (x < xa) && lowerbounderror(x, xa, "first")
    (x > xb) && upperbounderror(x, xb, "first")
    (y < ya) && lowerbounderror(y, ya, "second")
    (y > yb) && upperbounderror(y, yb, "second")
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
    #promote types
    x, y = promote(collect(x), collect(y))
    #check for basic problems
    rangecheck(x, y)
    #grid properties
    n = length(x)
    xa = minimum(x)
    xb = maximum(x)
    #construct
    InterpolatorRange(n, x, xa, xb, y)
end

function rangecheck(x::AbstractVector, y::AbstractVector)
    nx, ny = length(x), length(y)
    @assert nx == ny "length of x ($nx) does not match length of y ($ny)"
    @assert all(diff(x) .> 0.0) "grid points must be in strictly ascending order"
end

function linstruct(T::Type,
                   f::Function,
                   xa::Real,
                   xb::Real,
                   n::Int,
                   boundaries::AbstractBoundaries)
    x = LinRange(xa, xb, n)
    y = f.(x)
    T(x, y, boundaries)
end

values(ϕ::OneDimensionalInterpolator) = ϕ.r.y

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
    #promote types
    x, y = promote(collect(x), collect(y))
    Z = collect(typeof(x[1]), Z)
    #check for basic grid problems
    gridcheck(x, y, Z)
    #boundaries
    xa, xb = minimum(x), maximum(x)
    ya, yb = minimum(y), maximum(y)
    #construct the object
    InterpolatorGrid(length(x), length(y), x, y, xa, xb, ya, yb, Z)
end

function gridcheck(x::AbstractVector, y::AbstractVector, Z::AbstractMatrix)
    @assert (length(x) > 1) & (length(y) > 1) "grid must have more than one point in each direction"
    @assert size(Z) == (length(x), length(y)) "dimensions mismatched, size(Z) = $(size(Z)) but length(x) = $(length(x)) and length(y) = $(length(y))"
    @assert (all(diff(x) .> 0.0) & all(diff(y) .>= 0.0)) "grid coordinates must be monotonically increasing without duplicates"
end

function linstruct(T::Type, f::Function,
                   xa::Real, xb::Real, nx::Int,
                   ya::Real, yb::Real, ny::Int,
                   boundaries::AbstractBoundaries)
    x = LinRange(xa, xb, nx)
    y = LinRange(ya, yb, ny)
    X = x .* ones(ny)'
    Y = y' .* ones(nx)
    Z = f.(X, Y)
    T(x, y, Z, boundaries)
end

values(ϕ::TwoDimensionalInterpolator) = ϕ.G.Z