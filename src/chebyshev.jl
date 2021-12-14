export chebygrid, ischebygrid, chebycoef, cheby, chebyderiv
export ChebyshevInterpolator, BichebyshevInterpolator

ξ2x(ξ, a, b) = (ξ + 1)*((b - a)/2) + a

x2ξ(x, a, b) = 2*(x - a)/(b - a) - 1

x2θ(x, a, b) = acos(x2ξ(x, a, b))

"""
    chebygrid(n)

Create an array of `n` chebyshev nodes in [-1,1]
"""
@memoize chebygrid(n::Int) = cos.(π*(n-1:-1:0)/(n-1))

"""
    chebygrid(xa, xb, n)

Create an array of `n` chebyshev nodes in [`xa`,`xb`]
"""
chebygrid(xa, xb, n::Int) = ξ2x.(chebygrid(n), xa, xb)

"""
    chebygrid(xa, xb, nx, ya, yb, ny)

Create a two-dimensional grid of chebyshev nodes using `nx` points along the first axis, in [`xa`,`xb`], and `ny` points along the second axis, in [`ya`,`yb`].
"""
function chebygrid(xa, xb, nx::Int, ya, yb, ny::Int)
    X = chebygrid(xa, xb, nx) .* ones(ny)'
    Y = chebygrid(ya, yb, ny)' .* ones(nx)
    return X, Y
end

function ischebygrid(x)::Bool
    n = length(x)
    c = chebygrid(n)
    xa, xb = minimum(x), maximum(x)
    @inbounds for i ∈ eachindex(x)
        ξ = x2ξ(x[i], xa, xb)
        r = abs((ξ - c[i])/c[i]) #relative difference
        r > 1e-3 && return false 
    end
    return true
end

@memoize function chebymatrix(n::Int)
    @assert n > 1 "can't construct cheby matrix smaller than 2 x 2"
    A = zeros(n,n)
    ξ = chebygrid(n)
    for j ∈ 1:n, k ∈ 1:n
        @inbounds A[j,k] = cos((k-1)*acos(ξ[j]))
    end
    return A
end

function chebyrecurrance!(T, ξ, L::Int)::Nothing
    Tₖ₋₂ = one(ξ)
    Tₖ₋₁ = ξ
    @inbounds T[2] = Tₖ₋₁
    for k = 3:L
        #compute next value
        Tₖ = 2ξ*Tₖ₋₁ - Tₖ₋₂
        #set array value
        @inbounds T[k] = Tₖ
        #swaps
        Tₖ₋₂ = Tₖ₋₁
        Tₖ₋₁ = Tₖ
    end
    nothing
end

function chebyrecurrance(ξ::U, L::Int) where {U}
    T = Vector{U}(undef, L)
    @inbounds T[1] = one(ξ)
    chebyrecurrance!(T, ξ, L)
    return T
end

"""
    chebycoef(y)

Compute the Chebyshev expansion coefficients for a set of points `y`, which are assumed to be located on the Chebyshev points for some interval.
"""
function chebycoef(y)
    n = length(y)
    @assert n > 1 "must have at least 2 points to form chebyshev coefficients"
    Aᴵ = invertedchebymatrix(n)
    return Aᴵ*y
end

"""
    cheby(coef, x, xa, xb)

Evaluates the Chebyshev expansion represented by the coefficients in `coef` and defined on the interval [`xa`,`xb`] at the point `x`.
"""
function cheby(coef, x, xa, xb)
    N = length(coef)
    @assert xa <= x <= xb "x must be inside [xa,xb]"
    #get coordinate in ξ space
    ξ = x2ξ(x, xa, xb)
    #first two elements of cheby recursion
    Tₖ₋₂ = one(ξ)
    Tₖ₋₁ = ξ
    #first two terms of dot product
    @inbounds y = Tₖ₋₂*coef[1] + Tₖ₋₁*coef[2]
    #cheby recursion and rest of terms in dot product, all at once
    for k = 3:N
        #next value in recursion
        Tₖ = 2*ξ*Tₖ₋₁ - Tₖ₋₂
        #next term in dot product
        @inbounds y += Tₖ*coef[k]
        #swaps
        Tₖ₋₂ = Tₖ₋₁
        Tₖ₋₁ = Tₖ
    end
    return y
end

"""
    chebyderiv(coef, xa, xb)

Generates the expansion coefficents for the derivative of a preexisting Chebyshev expansion defined on the interval [`xa`,`xb`].
"""
function chebyderiv(coef, xa, xb)
    n = length(coef)
    T = eltype(coef)
    @assert xa < xb "xa must be less than xb"
    #length 1 case is a trivial zero slope
    n == 1 && return zeros(T, 1)
    #start the array and check for length 2 case
    d = zeros(T, n-1)
    d[n-1] = 2*(n-2)*coef[n]
    n == 2 && return d
    #recurrance for all the rest
    d[n-2] = 2*(n-3)*coef[n-1]
    @inbounds for k ∈ n-2:-1:2
        d[k-1] = d[k+1] + 2*(k-1)*coef[k]
    end
    #interval width factor
    d .*= 2/(xb - xa)
    #extra DC scaling
    d[1] /= 2
    return d
end

#-------------------------------------------------------------------------------
# caching function for inverted cheby matrices, needed for interpolator setup

@memoize invertedchebymatrix(n::Int64)::Matrix{Float64} = inv(chebymatrix(n))

#-------------------------------------------------------------------------------
# one-dimensional interpolation

struct ChebyshevInterpolator{N,T}
    #lowest value in range
    xa::T
    #highest value in range
    xb::T
    #interpolation coefficents
    a::NTuple{N,T}
    #must always have strict boundaries
    boundaries::StrictBoundaries
end

function Base.show(io::IO, ϕ::ChebyshevInterpolator)
    print(io, "$(typeof(ϕ)) ∈ [$(ϕ.xa), $(ϕ.xb)]")
end

"""
    ChebyshevInterpolator(x, y)

Construct a `ChebyshevInterpolator` for the points defined by coordinates `x` and values `y`. The `x` coordinates *must* be arranged on a chebyshev grid, which can be generated using the [`chebygrid`](@ref) function.
"""
function ChebyshevInterpolator(x, y)
    #same types
    T = promote_type(eltype(x), eltype(y))
    x = collect(T, x)
    y = collect(T, y)
    #check for basic issues
    rangecheck(x, y, 3)
    #demand that the input points have chebyshev spacing
    @assert ischebygrid(x) "points must be on a chebyshev grid"
    #generate expansion coefficients
    a = Tuple(chebycoef(y))
    #construct
    ChebyshevInterpolator(minimum(x), maximum(x), a, StrictBoundaries())
end

"""
    ChebyshevInterpolator(f, xa, xb, n)

Construct a `ChebyshevInterpolator` for the function `f` using `n` function evaluations in the range [`xa`,`xb`]. The function evaluations will occur on the chebyshev nodes.
"""
function ChebyshevInterpolator(f::F, xa, xb, n::Int) where {F}
    #set up the range coordinates
    x = chebygrid(xa, xb, n)
    #evaluate the function at those coordinates
    y = f.(x)
    #call the other constructor
    ChebyshevInterpolator(x, y)
end

function (ϕ::ChebyshevInterpolator{N,U})(x) where {N,U}
    #always enforce boundaries
    ϕ.boundaries(x, ϕ.xa, ϕ.xb)
    #evaluate the expansion
    return cheby(ϕ.a, x, ϕ.xa, ϕ.xb)
end

"""
    chebyderiv(ϕ::ChebyshevInterpolator)

Construct a ChebyshevInterpolator representing the derivative of a preexisting interpolator.
"""
function chebyderiv(ϕ::ChebyshevInterpolator)
    #derivative expansion's coefficients
    coef = Tuple(chebyderiv(ϕ.a, ϕ.xa, ϕ.xb))
    #new interpolator
    ChebyshevInterpolator(ϕ.xa, ϕ.xb, coef, StrictBoundaries())
end

#-------------------------------------------------------------------------------
# bichebyshev interpolation, a little trickier now!

struct BichebyshevInterpolator{M,N,U}
    #lowest value on axis 1
    xa::U
    #highest value on axis 1
    xb::U
    #lowest value on axis 2
    ya::U
    #highest value on axis 2
    yb::U
    #matrix and vectors for doing the interpolation
    A::Matrix{U} # size (ny by nx) or (M by N)
    a::Vector{U} # length ny for cosine expansion in θy
    b::Vector{U} # length nx for cosine expansion in θx
    c::Vector{U} # length ny for doing M*b in place
    #must always use strict boundaries
    boundaries::StrictBoundaries
end

function Base.show(io::IO, ϕ::BichebyshevInterpolator)
    print(io, "$(typeof(ϕ)) ∈ [$(ϕ.xa), $(ϕ.xb)], [$(ϕ.ya), $(ϕ.yb)]")
end

"""
    BichebyshevInterpolator(x, y, Z)

Construct a `BichebyshevInterpolator` for the grid of points defined by coordinates (`x`,`y`) and values `Z`. The given points must lie on a chebyshev grid in each direction. These can be generated with the [`chebygrid`](@ref) function or the interpolator can be constructed directly from a function using the method below.

!!! warning

    The Bichebyshev interpolator is *not thread-safe*. It computes a cosine expansion and does some linear algebra in-place using arrays stored with the object. A single `BichebyshevInterpolator` should never be called by multiple threads at once.
"""
function BichebyshevInterpolator(x, y, Z)
    #same types
    T = promote_type(eltype(x), eltype(y), eltype(Z))
    x = collect(T, x)
    y = collect(T, y)
    Z = collect(T, Z)
    #check for basic grid problems
    gridcheck(x, y, Z, 3)
    #grid properties
    nx, ny = length(x), length(y)
    xa, xb = T(minimum(x)), T(maximum(x))
    ya, yb = T(minimum(y)), T(maximum(y))
    #reject any non-cheby grid spacing
    @assert ischebygrid(x) "axis 1 coordinates must be on a chebyshev grid"
    @assert ischebygrid(y) "axis 2 coordinates must be on a chebyshev grid"
    #get inverted matrices from cache or generate them
    B = invertedchebymatrix(nx)
    #generate interpolation coefficients along axis 1 for each value of axis 2
    α = zeros(T, nx, ny)
    for j = 1:ny
        mul!(view(α,:,j), B, view(Z,:,j))
    end
    #then combine α and A
    A = invertedchebymatrix(ny)*α'
    #other vectors we need for doing the actual interpolation
    a = ones(T, ny)
    b = ones(T, nx)
    c = zeros(T, ny)
    #done
    BichebyshevInterpolator{ny,nx,T}(xa, xb, ya, yb, A, a, b, c, StrictBoundaries())
end

"""
    BichebyshevInterpolator(f, xa, xb, nx, ya, yb, ny)

Construct a `BichebyshevInterpolator` for the function `f` using a grid of `nx` points on the first axis in [`xa`,`xb`] and `ny` points on the second axis in [`ya`,`yb`].
"""
function BichebyshevInterpolator(f::F, xa, xb, nx::Int, ya, yb, ny::Int) where {F}
    #set up the grid
    X, Y = chebygrid(xa, xb, nx, ya, yb, ny)
    #evaluate the function at chebyshev grid points
    Z = f.(X, Y)
    #call the other constructor
    BichebyshevInterpolator(X[:,1], Y[1,:], Z)
end

#=====
This is the fast implementation. It's executed when the types of the
input coordinates match the type of the stored coefficients in the
interpolator. When the types match, the Chebyshev expansions can be
evaluated in-place, using vectors pre-allocated in the interpolator.
See the a, b, and c fields of the struct. This method also guarantees
that the interpolator type and the coordinate types are <: AbstractFloat,
making the low-level linear algebra functions safe. Without such a
guarantee, there can be issues with, for example, the dual numbers in
FowardDiff routines.
=====#
function (Φ::BichebyshevInterpolator{M,N,U})(x::U, y::U) where {M,N,U<:AbstractFloat}
    #always enforce boundaries for Chebyshev
    Φ.boundaries(x, Φ.xa, Φ.xb, y, Φ.ya, Φ.yb)
    #evaluate Chebyshev polys at the coordinates recursively and in-place
    ξy = x2ξ(y, Φ.ya, Φ.yb)
    chebyrecurrance!(Φ.a, ξy, M)
    ξx = x2ξ(x, Φ.xa, Φ.xb)
    chebyrecurrance!(Φ.b, ξx, N)
    #perform M*b, which interpolates along the first axis, also in-place
    mul!(Φ.c, Φ.A, Φ.b)
    #then a'*c interpolates along the second axis
    return dot(Φ.a, Φ.c)
end

#=====
This is the slow, but type-flexible, implementation. It's
executed whenever the type of the interpolator's coefficients or
the coordinate types don't match OR are not AbstractFloats. The
price of this flexibility is allocations for the expansions and
loss of the low-level linear algebra routines 😭.
=====#
function (Φ::BichebyshevInterpolator{M,N,U})(x, y) where {M,N,U}
    #always enforce boundaries for Chebyshev
    Φ.boundaries(x, Φ.xa, Φ.xb, y, Φ.ya, Φ.yb)
    #coordinates in ξ space
    ξx, ξy = promote(x2ξ(x, Φ.xa, Φ.xb), x2ξ(y, Φ.ya, Φ.yb))
    #allocating expansion in the x direction
    b = chebyrecurrance(ξx, N)
    c = Φ.A * b
    #then perform the y-axis expansion and dot product simulataneously    
    Tₖ₋₂ = one(ξy)
    Tₖ₋₁ = ξy
    @inbounds z = c[1] + ξy*c[2]
    for i = 3:M
        #compute next value
        Tₖ = 2ξy*Tₖ₋₁ - Tₖ₋₂
        #update running dot product
        @inbounds z += c[i]*Tₖ
        #swaps
        Tₖ₋₂ = Tₖ₋₁
        Tₖ₋₁ = Tₖ
    end
    return z
end