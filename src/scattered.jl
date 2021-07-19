export ShepardInterpolator, RBFInterpolator
export multiquadratic, invmultiquadratic, gaussian

using LinearAlgebra: Symmetric

#-------------------------------------------------------------------------------

function distance(a::AbstractVector{Float64},
                  b::AbstractVector{Float64},
                  n::Int64)::Float64
    r = 0.0
    for i = 1:n
        @inbounds r += (a[i] - b[i])*(a[i] - b[i])
    end
    sqrt(r)
end

function checkcoords(X::AbstractArray{<:Real,2}, y::AbstractVector{<:Real})
    @assert size(X)[2] > 1 "generally it doesn't make sense to apply scattered interpolation methods in one dimension, consider using another interpolator"
    p = size(X)[1]
    n = length(y)
    @assert p == n "number of points in X ($p) does not match number of values in y ($n)"
end

function checkdims(x::AbstractVector{<:Real}, N::Int64)
    n = length(x)
    @assert n == N "Number of interpolation coordinates ($n) does not match number of interpolator dimensions ($N)"
end

struct ScatteredPoints
    #number of dimensions
    N::Int64
    #number of points
    p::Int64
    #table of points, p x N
    X::Array{Float64,2}
    #vector of values
    y::Vector{Float64}
end

function ScatteredPoints(X::AbstractArray{<:Real,2}, y::AbstractVector{<:Real})
    checkcoords(X, y)
    X = convert.(Float64, X)
    y = convert.(Float64, y)
    N = size(X)[2]
    p = size(X)[1]
    ScatteredPoints(N, p, X, y)
end

#-------------------------------------------------------------------------------

"""
    multiquadratic(r, ϵ)

``\\sqrt{1 + (r/ϵ)^2}``
"""
multiquadratic(r, ϵ) = sqrt(1 + (r/ϵ)^2)

"""
    invmultiquadratic(r, ϵ)

``\\frac{1}{\\sqrt{1 + (r/ϵ)^2}}``
"""
invmultiquadratic(r, ϵ) = 1/multiquadratic(r, ϵ)

"""
    gaussian(r, ϵ)

``e^{-r^2/ϵ^2}``
"""
gaussian(r, ϵ) = exp(-r^2/ϵ^2)

struct RBFInterpolator{T<:Function}
    #interpolation points
    S::ScatteredPoints
    #radial basis function
    rbf::T
    #distance scale factor
    ϵ::Float64
    #function weights
    w::Vector{Float64}
end

"""
    RBFInterpolator(X, y, ϵ, rbf=invmultiquadratic)

Construct a radial basis function (RBF) interpolator for an n-dimensional set of points with coordinates `X` and values `y`. `X` must be an p × N array, where p is the number of points and N is the number of dimensions. `y` must be a length p vector. The value of `ϵ` scales the radial basis function of choice, `f`, which is [`invmultiquadratic`](@ref) by default. Any function in the form ``ϕ(r,ϵ)`` can be passed to the `rbf` argument, where \$r\$ is the distance between points and ``ϵ`` is a scaling factor.
"""
function RBFInterpolator(X::AbstractArray{<:Real,2},
                         y::AbstractVector{<:Real},
                         ϵ::Real,
                         rbf::F=multiquadratic) where {F<:Function}
    S = ScatteredPoints(X, y)
    ϵ = Float64(ϵ)
    #matrix of basis function evaluations
    A = ones(S.p, S.p)
    for i = 1:S.p
        a = view(S.X,i,:)
        for j = i+1:S.p
            b = view(S.X,j,:)
            r = distance(a, b, S.N)
            A[i,j] = rbf(r, ϵ)
        end
    end
    #solve for weights
    w = Symmetric(A)\S.y
    #construct
    RBFInterpolator(S, rbf, ϵ, w)
end

function RBF(Φ::RBFInterpolator, x::Vector{Float64})::Float64
    checkdims(x, Φ.S.N)
    y = 0.0
    for i = 1:Φ.S.p
        #distance
        r = distance(x, view(Φ.S.X,i,:), Φ.S.N)
        #weighted evaluation
        y += Φ.w[i]*Φ.rbf(r, Φ.ϵ)
    end
    return y
end

function (Φ::RBFInterpolator)(x::AbstractVector{Float64})::Float64
    RBF(Φ, x)
end

function (Φ::RBFInterpolator)(x::AbstractVector{<:Real})::Float64
    RBF(Φ, collect(Float64, x))
end

function (Φ::RBFInterpolator)(x::Real...)::Float64
    RBF(Φ, collect(Float64, x))
end

#-------------------------------------------------------------------------------

struct ShepardInterpolator{T}
    #interpolation points
    S::ScatteredPoints
    #power-law exponent
    a::T
end

"""
    ShepardInterpolator(X, y, a=3)

Construct a `ShepardInterpolator` for an n-dimensional set of points with coordinates `X` and values `y`. `X` must be an p × N array, where p is the number of points and N is the number of dimensions. `y` must be a length p vector. The value of `a` defines the distance weighting function ``r^{-a}``.
"""
function ShepardInterpolator(X::AbstractArray{<:Real,2},
                             y::AbstractVector{<:Real},
                             a::Real=3)
    ShepardInterpolator(ScatteredPoints(X, y), a)
end

function shepard(Φ::ShepardInterpolator, x::Vector{Float64})::Float64
    checkdims(x, Φ.S.N)
    n = 0.0
    d = 0.0
    for i = 1:Φ.S.p
        #distance
        r = distance(x, view(Φ.S.X,i,:), Φ.S.N)
        #evaluate weighting function
        f = 1.0/r^Φ.a
        #weigted and unweighted contribution to numertor and denominator
        n += Φ.S.y[i]*f
        d += f
    end
    return n/d
end

function (Φ::ShepardInterpolator)(x::AbstractVector{Float64})::Float64
    shepard(Φ, x)
end

function (Φ::ShepardInterpolator)(x::AbstractVector{<:Real})::Float64
    shepard(Φ, collect(Float64, x))
end

function (Φ::ShepardInterpolator)(x::Real...)::Float64
    shepard(Φ, collect(Float64, x))
end
