# BasicInterpolators.jl
![downloads](https://shields.io/endpoint?url=https://pkgs.genieframework.com/api/v1/badge/BasicInterpolators)

*Interpolation methods with a simple interface*

*contributions welcome*

| Documentation | Status |
| :-----------: | :----: |
| [![docs](https://img.shields.io/badge/docs-latest-blue.svg)](https://markmbaum.github.io/BasicInterpolators.jl/dev)  | [![Build Status](https://github.com/markmbaum/BasicInterpolators.jl/workflows/CI/badge.svg)](https://github.com/markmbaum/BasicInterpolators.jl/actions) [![codecov](https://codecov.io/gh/markmbaum/BasicInterpolators.jl/branch/main/graph/badge.svg?token=yRg33tFcL3)](https://codecov.io/gh/markmbaum/BasicInterpolators.jl)  |


Use Julia's package manager to install
```
julia> ]add BasicInterpolators
```

-----

### Interpolation Methods

##### One Dimension

- [x] linear
- [x] piecewise cubic
- [x] cubic spline (natural or clamped)
- [x] Chebyshev
- [x] arbitrary order polynomials (Neville's method)
- [x] polynomial coefficients (efficient Vandermonde solver)
- [x] end-point cubic Hermite

##### Two Dimensions, Regular Grid

- [x] linear
- [x] piecewise cubic
- [x] Chebyshev

##### N-Dimensions, Scattered Points

- [x] radial basis functions (any choice of function)
- [x] Shepard

-----

### Basic Usage

See the [**tutorial**](https://markmbaum.github.io/BasicInterpolators.jl/dev/tutorial/) for more complete examples.

```julia
using BasicInterpolators, ForwardDiff

#some data to interpolate
x = [-1, 0.5, 2, 3]
y = [1, 3, -0.5, 0]

#a linear interpolation struct
p = LinearInterpolator(x, y)

#interpolate at one point
p(2.5)

#interpolate at lots of points
p.(LinRange(-1, 3, 100))

#compute the derivative dy/dx
ForwardDiff.derivative(p, 1.0)

#make an interpolator that doesn't check boundaries (allows extrapolation)
p = LinearInterpolator(x, y, NoBoundaries())
```

-----

### Other packages

Some notable packages with other/advanced methods:

1. [Interpolations.jl](https://github.com/JuliaMath/Interpolations.jl)
2. [Dierckx.jl](https://github.com/kbarbary/Dierckx.jl)
3. [GridInterpolations.jl](https://github.com/sisl/GridInterpolations.jl)
4. [ApproXD.jl](https://github.com/floswald/ApproXD.jl)
5. [FastChebInterp.jl](https://github.com/stevengj/FastChebInterp.jl)
6. [ApproxFun.jl](https://github.com/JuliaApproximation/ApproxFun.jl)

For a longer list, look [here](https://github.com/JuliaMath/Interpolations.jl#other-interpolation-packages).