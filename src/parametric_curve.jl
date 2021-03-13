#-------------------------------------------------------------------------------
# parametric cubic splines for n dimensions

struct ParametricCurveInterpolator
    #interpolators for each dimension
    Φ::Vector{CubicSplineInterpolator}
    #number of dimensions
    N::Int64
end

"""
    ParametricCurveInterpolator(V...)

Construct an interpolator for a set of points in arbitrary dimensions that defines a one-dimensional curve, using natural cubic splines in each dimension.
"""
function ParametricCurveInterpolator(V::AbstractArray{<:Real}...; closed=false)
    #standardize types
    V = [collect(Float64, v) for v ∈ V]
    #check lengths
    L = length(V[1])
    N = length(V)
    for i = 2:N
        @assert length(V[i]) == L "coordinate vectors must all be equal length"
    end
    #compute chordlengths along the curve
    s = zeros(L)
    for i = 2:L
        #segment distance
        d = 0
        for j = 1:N
            d += (V[j][i] - V[j][i-1])^2
        end
        d = sqrt(d)
        #incremental chord length
        s[i] = s[i-1] + d
    end
    #normalize s to have a maximum of 1
    s /= s[end]
    #construct interpolators
    Φ = [CubicSplineInterpolator(s, v) for v ∈ V]
    #construct global interpolator
    ParametricCurveInterpolator(Φ, N)
end

function (C::ParametricCurveInterpolator)(s::Real,
                                          bounds::Bool=true)::Vector{Float64}
    @assert 0 <= s <= 1 "interpolation coordinate must be in [0,1] for parametric curves"
    [ϕ(s, bounds) for ϕ ∈ C.Φ]
end
