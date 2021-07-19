module BasicInterpolators

#for previous cell indices
using Base: RefValue
#to get value arrays of interpolators
import Base.values
#for in-place matrix multiplication in chebyshev interpolators
using LinearAlgebra: mul!, dot

include("base.jl")
include("polynomial.jl")
include("spline.jl")
include("chebyshev.jl")
include("scattered.jl")

end
