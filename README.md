# BasicInterpolators.jl

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://wordsworthgroup.github.io/BasicInterpolators.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://wordsworthgroup.github.io/BasicInterpolators.jl/dev)
[![Build Status](https://github.com/wordsworthgroup/BasicInterpolators.jl/workflows/CI/badge.svg)](https://github.com/wordsworthgroup/BasicInterpolators.jl/actions)
[![codecov](https://codecov.io/gh/wordsworthgroup/BasicInterpolators.jl/branch/main/graph/badge.svg?token=yRg33tFcL3)](https://codecov.io/gh/wordsworthgroup/BasicInterpolators.jl)

-----

This Julia package provides straightforward access to basic interpolation algorithms in one and two dimensions. The methods are
* Linear
* Piecewise cubic (no smoothness guarantee)
* Cubic spline (guarantees smoothness)
* Chebyshev (extremely high accuracy for the right application)

There are also functions for arbitrary polynomial interpolation in one dimension.

See the [documentation](https://wordsworthgroup.github.io/BasicInterpolators.jl/stable) for examples, discussion, and details.

-----

More ambitions packages for more advanced applications include
1. [Interpolations.jl](https://github.com/JuliaMath/Interpolations.jl)
2. [Dierckx.jl](https://github.com/kbarbary/Dierckx.jl)
3. [GridInterpolations.jl](https://github.com/sisl/GridInterpolations.jl)
4. [ApproXD](https://github.com/floswald/ApproXD.jl)
