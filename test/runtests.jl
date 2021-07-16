pushfirst!(LOAD_PATH, "C:/Users/markm/Documents/code")

using Test, BasicInterpolators

@testset "Accuracy" begin include("test_accuracy.jl") end
@testset "Boundary" begin include("test_boundary.jl") end
@testset "Nevilles" begin include("test_nevilles.jl") end
@testset "Vandermonde" begin include("test_vandermonde.jl") end
@testset "Scattered" begin include("test_scattered.jl") end
@testset "Other" begin include("test_other.jl") end
