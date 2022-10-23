module BasicInterpolators

using LinearAlgebra: mul!, dot
using Memoize

include("base.jl")
include("polynomial.jl")
include("spline.jl")
include("chebyshev.jl")
include("scattered.jl")

end
