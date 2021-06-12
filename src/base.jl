export findcell, NoBoundaries, WeakBoundaries, StrictBoundaries

"""
    findcell(q, V)

Use bisection search to find the cell containing `q`, assuming `V` is a sorted vector of coordinates. The returned integer is the index of the element in `V` immediately less than `q`. For example, if findcell returns 2, then `q ∈ [V[2],V[3])`. If `q` is less than every element in `V`, 1 is returned, indicating the first cell in `V`. If `q` is greater than every element in `V`, `length(V)-1` is returned, indicating the last cell in `V`.
"""
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

#-------------------------------------------------------------------------------
# types for handling boundaries

function upperbounderror(x, xb, axis::Int)
    error("Interpolation location $x outside upper interpolation limit $xb on axis $axis")
end

function lowerbounderror(x, xa, axis::Int)
    error("Interpolation location $x outside lower interpolation limit $xa on axis $axis")
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
    if (x < xa) && !(x ≈ xa)
        lowerbounderror(x, xa, 1)
    end
    if (x > xb) && !(x ≈ xb)
        upperbounderror(x, xb, 1)
    end
end

function (B::WeakBoundaries)(x, xa, xb, y, ya, yb)
    if (x < xa) && !(x ≈ xa)
        lowerbounderror(x, xa, 1)
    end
    if (x > xb) && !(x ≈ xb)
        upperbounderror(x, xb, 1)
    end
    if (y < ya) && !(y ≈ ya)
        lowerbounderror(y, ya, 2)
    end
    if (y > yb) && !(y ≈ yb)
        upperbounderror(y, yb, 2)
    end
end

#--------------------------------

"""
    StrictBoundaries()

Triggers errors whenever interpolation coordinates are outside of the boundaries.
"""
struct StrictBoundaries <: AbstractBoundaries end

function (B::StrictBoundaries)(x, xa, xb)
    if (x < xa)
        lowerbounderror(x, xa, 1)
    end
    if (x > xb)
        upperbounderror(x, xb, 1)
    end
end

function (B::StrictBoundaries)(x, xa, xb, y, ya, yb)
    if (x < xa)
        lowerbounderror(x, xa, 1)
    end
    if (x > xb)
        upperbounderror(x, xb, 1)
    end
    if (y < ya)
        lowerbounderror(y, ya, 2)
    end
    if (y > yb)
        upperbounderror(y, yb, 2)
    end
end

#-------------------------------------------------------------------------------
# a base structure for 1d interpolations and some related functions

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
    x = collect(LinRange(xa, xb, n))
    y = f.(x)
    T(x, y, boundaries)
end

#-------------------------------------------------------------------------------
# a base structure for 2d grid interpolation and some related functions

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
    #boundaries
    xa, xb, ya, yb = minimum(x), maximum(x), minimum(y), maximum(y)
    #construct the object
    InterpolatorGrid(length(x), length(y), x, y, xa, xb, ya, yb, Z)
end

function gridcheck(x::AbstractVector{<:Real},
                   y::AbstractVector{<:Real},
                   Z::AbstractArray{<:Real,2})
    @assert (length(x) > 1) & (length(y) > 1) "grid must have more than one point in each direction"
    @assert size(Z) == (length(x), length(y)) "dimensions mismatched, size(Z) = $(size(Z)) but length(x) = $(length(x)) and length(y) = $(length(y))"
    @assert (all(diff(x) .> 0.0) & all(diff(y) .>= 0.0)) "grid coordinates must be monotonically increasing without duplicates"
end

function linstruct(T::Type, f::Function,
                   xa::Real, xb::Real, nx::Int,
                   ya::Real, yb::Real, ny::Int,
                   boundaries::AbstractBoundaries)
    x = collect(LinRange(xa, xb, nx))
    y = collect(LinRange(ya, yb, ny))
    X = x .* ones(ny)'
    Y = y' .* ones(nx)
    Z = f.(X, Y)
    T(x, y, Z, boundaries)
end
