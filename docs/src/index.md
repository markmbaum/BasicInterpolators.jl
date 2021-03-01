# BasicInterpolators.jl

Welcome to the `BasicInterpolators` documentation. This package provides some basic and useful interpolation algorithms. Quick facts about the functionality:
* All types and functions assume the use of real numbers and will return double precision (`Float64`) numbers.
* There are four different interpolation methods, each implemented in one and two dimensions.
* The two-dimensional interpolation methods operate regular grids.
* Attempts to extrapolate can be set to throw an error or blindly use whatever curve/surface defines the nearest boundary cell.

--------------------------------------------------------------------------------
## Other Packages

If you need n-dimensional interpolation, irregular grids, complex numbers, or built-in treatment of different boundary behaviors, this package will not help you. But a lot of applications don't require those features and I hope this package will be useful and fast for basic applications. If you do need other features, however, please have a look at more ambitious packages,
1. [Interpolations.jl](https://github.com/JuliaMath/Interpolations.jl)
2. [Dierckx.jl](https://github.com/kbarbary/Dierckx.jl)
3. [GridInterpolations.jl](https://github.com/sisl/GridInterpolations.jl)
4. [ApproXD](https://github.com/floswald/ApproXD.jl)

--------------------------------------------------------------------------------

## Details

For an overview of how the package works, see the how-to section below. For more information on the specific types and methods available, follow these links:
* [Linear Interpolators](linear.md)
* [Cubic Interpolators](cubic.md)
* [Cubic Spline Interpolators](spline.md)
* [Chebyshev Interpolators](chebyshev.md)
* [General 1D Polynomial Interpolation](polynomial.md)

--------------------------------------------------------------------------------
## How to Start Interpolating Stuff

`BasicInterpolators` offers callable types, or [function-like objects](https://docs.julialang.org/en/v1/manual/methods/#Function-like-objects), for your interpolation needs. This simply means that you construct an `Interpolator` object, like a `LinearInterpolator` or a `BicubicInterpolator`, then call it like a function, passing it the interpolation coordinate.

#### Construct from points

 For example, to do linear interpolation on a range of three points, we construct a [`LinearInterpolator`](@ref) then call it on the interpolation coordinate:
```@repl
using BasicInterpolators: LinearInterpolator
x = [0, 1, 2];
y = [0, 0.5, 1];
p = LinearInterpolator(x, y);
p(0.5)
```
Pretty basic! Each of the types can be constructed by passing vectors/arrays with the points to interpolate. A nice result of this syntax is easy broadcasting. To interpolate a bunch of values in some array `x`, we simply call `p.(x)`.

#### Construct from a function

All the interpolators can also be constructed with a function to sample from. For example, to generate a [`ChebyshevInterpolator`](@ref) of the function $f(x) = \sin(2x) + x^2\exp(x)$, in the range $[-3,0]$, using 16 points:
```@repl
using BasicInterpolators: ChebyshevInterpolator
f(x) = sin(2x) + x^2*exp(x);
p = ChebyshevInterpolator(f, -3, 0, 16);
f(-1.5)
p(-1.5)
```
In the background, the constructor generates 16 points between -3 and 0, evaluates the function at those points, then sets up everything it needs to interpolate in that range.

#### Two-dimensional interpolation

The two-dimensional interpolators work the same way, but with another dimension (of course). You can supply the coordinates along each axis and a grid of points or you can supply a function along with where to evaluate it and how many points to use. For example, to make a `BilinearInterpolator` for a 10 by 10 grid of randomly spaced points:
```@repl
using BasicInterpolators: BilinearInterpolator
x = sort(rand(10));
y = sort(rand(10));
A = rand(10,10);
P = BilinearInterpolator(x, y, A);
P(0.5, 0.5)
```

#### Grid spacing

Unevenly spaced grids are acceptable for everything but the [chebyshev interpolators](chebyshev.md) and the [bivariate spline](spline.md), but you have to supply the coordinates directly instead of constructing with a function. When you do construct with a function, all the interpolators use evenly spaced points except for the chebyshev interpolators. The chebyshevs are a special case requiring a very specific grid spacing but they reward you with shockingly high accuracy as the number of points increases.

#### Boundaries

The [chebyshev interpolators](chebyshev.md) won't let you extrapolate, or pass a point outside the stated interpolation range, because they are undefined outside the boundaries. All the other ones will let you extrapolate if you turn off the boundary check, which is an optional argument whenever you call an interpolator. For example,
```@repl
using BasicInterpolators: CubicSplineInterpolator
x = [0, 1, 2];
y = rand(3);
p = CubicSplineInterpolator(x, y);
#trying to extrapolate will normally cause an error
p(-1)
#but you can override the boundary check by passing `false`
p(-1, false)
```
Extrapolation blindly uses the interpolating curve/surface of the nearest boundary cell in the grid, so it can easily return wild values, especially for the cubic interpolators. That's why it's not allowed by default. The same goes for two-dimensional interpolators. You can turn off the boundary check by passing `false` after the interpolation coordiantes, as in `P(x, y, false)`. In these cases, the `false` argument means "don't check the boundaries."

#### Methods

* [Linear Interpolation](linear.md)
* [Cubic Interpolation](cubic.md)
* [Cubic Spline Interpolation](spline.md)
* [Chebyshev Interpolation](chebyshev.md)
* [General Polynomial Interpolation in 1D](polynomial.md)
