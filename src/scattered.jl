export ShepardInterpolator, RBFInterpolator
export multiquadratic, invmultiquadratic, gaussian

using LinearAlgebra: Symmetric

#-------------------------------------------------------------------------------

function distance(a, b, n::Int)
    r = 0.0
    for i = 1:n
        @inbounds r += (a[i] - b[i])*(a[i] - b[i])
    end
    sqrt(r)
end

function checkcoords(X, y)
    @assert size(X)[2] > 1 "Generally it doesn't make sense to apply scattered interpolation methods in one dimension. Consider using a 1D interpolator."
    p = size(X)[1]
    n = length(y)
    @assert p == n "number of points in X ($p) does not match number of values in y ($n)"
end

function checkdims(x, N::Int)
    n = length(x)
    @assert n == N "Number of interpolation coordinates ($n) does not match number of interpolator dimensions ($N)"
end

struct ScatteredPoints{T}
    #number of dimensions
    N::Int64
    #number of points
    p::Int64
    #table of points, p x N
    X::Array{T,2}
    #vector of values
    y::Vector{T}
end

function ScatteredPoints(X, y)
    checkcoords(X, y)
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

struct RBFInterpolator{T,F}
    #interpolation points
    S::ScatteredPoints{T}
    #radial basis function
    rbf::F
    #distance scale factor
    ϵ::T
    #function weights
    w::Vector{T}
end

"""
    RBFInterpolator(X, y, ϵ, rbf=invmultiquadratic)

Construct a radial basis function (RBF) interpolator for an n-dimensional set of points with coordinates `X` and values `y`. `X` must be an p × N array, where p is the number of points and N is the number of dimensions. `y` must be a length p vector. The value of `ϵ` scales the radial basis function of choice, `f`, which is [`invmultiquadratic`](@ref) by default. Any function in the form ``ϕ(r,ϵ)`` can be passed to the `rbf` argument, where \$r\$ is the distance between points and ``ϵ`` is a scaling factor.
"""
function RBFInterpolator(X::AbstractMatrix,
                         y::AbstractVector,
                         ϵ,
                         rbf::F=multiquadratic) where {F}
    S = ScatteredPoints(X, y)
    #matrix of basis function evaluations
    A = ones(eltype(y), S.p, S.p)
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
    RBFInterpolator(S, rbf, convert(eltype(y), ϵ), w)
end

function RBF(Φ::RBFInterpolator, x)
    checkdims(x, Φ.S.N)
    y = 0.0
    for i = 1:Φ.S.p
        #distance
        r = distance(x, view(Φ.S.X,i,:), Φ.S.N)
        #weighted evaluation
        @inbounds y += Φ.w[i]*Φ.rbf(r, Φ.ϵ)
    end
    return y
end

function (Φ::RBFInterpolator)(x::Union{AbstractVector,Tuple})
    RBF(Φ, x)
end

function (Φ::RBFInterpolator)(x::Number...)
    RBF(Φ, collect(Float64, x))
end

#-------------------------------------------------------------------------------

struct ShepardInterpolator{T}
    #interpolation points
    S::ScatteredPoints{T}
    #power-law exponent
    a::Float64
end

"""
    ShepardInterpolator(X, y, a=3)

Construct a `ShepardInterpolator` for an n-dimensional set of points with coordinates `X` and values `y`. `X` must be an p × N array, where p is the number of points and N is the number of dimensions. `y` must be a length p vector. The value of `a` defines the distance weighting function ``r^{-a}``.
"""
function ShepardInterpolator(X::AbstractMatrix,
                             y::AbstractVector,
                             a::Real=3.0)
    ShepardInterpolator(ScatteredPoints(X, y), a)
end

function shepard(Φ::ShepardInterpolator, x)
    checkdims(x, Φ.S.N)
    n = 0.0
    d = 0.0
    for i = 1:Φ.S.p
        #distance
        r = distance(x, view(Φ.S.X,i,:), Φ.S.N)
        #evaluate weighting function
        f = 1.0/r^Φ.a
        #weigted and unweighted contribution to numertor and denominator
        @inbounds n += Φ.S.y[i]*f
        d += f
    end
    return n/d
end

(Φ::ShepardInterpolator)(x::Union{AbstractVector,Tuple}) = shepard(Φ, x)

(Φ::ShepardInterpolator)(x::Number...) = shepard(Φ, x)
