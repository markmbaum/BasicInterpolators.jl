# BasicInterpolators.jl

Welcome to the `BasicInterpolators` documentation. This package provides a collection of common interpolation recipes. Quick facts about the package:
* The `Interpolator` types are the main focus. To use them, you construct an interpolator object then call it like a function. There are also some functions for things like one-off polynomial interpolation.
* Interpolators can be constructed with any numeric type that supports the usual arithmetic operations necessary for interpolating. They will, however, insist that the element types of your grid(s) and values are uniform, often converting to consistent element types automatically.
* There are two-dimensional grid interpolators, but not higher.
* All the interpolators should be compatible with automatic differentiation, specifically [ForwardDiff.jl](https://github.com/JuliaDiff/ForwardDiff.jl).

--------------------------------------------------------------------------------

## Details

For a tutorial, see the [**tutorial page**](tutorial.md). For more information on the specific types and methods available, see these links:
* [One-dimensional interpolation](1d.md)
* [Two-dimensional interpolation](2d.md)
* [N-dimensional, scattered point interpolation](scattered.md)

For a little extra on the Chebyshev interpolators, look [**here**](chebyshev.md).

--------------------------------------------------------------------------------
## Other Packages

If you need more features/methods/dimensions, please look in other packages:
1. [Interpolations.jl](https://github.com/JuliaMath/Interpolations.jl)
2. [Dierckx.jl](https://github.com/kbarbary/Dierckx.jl)
3. [GridInterpolations.jl](https://github.com/sisl/GridInterpolations.jl)
4. [ApproXD](https://github.com/floswald/ApproXD.jl)
5. [FastChebInterp.jl](https://github.com/stevengj/FastChebInterp.jl)