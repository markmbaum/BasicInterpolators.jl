`BasicInterpolators` offers callable types, or [function-like objects](https://docs.julialang.org/en/v1/manual/methods/#Function-like-objects), for your interpolation needs. This simply means that you construct an `Interpolator` object, like a `LinearInterpolator` or a `BicubicInterpolator`, then call it like a function, passing it the coordinate where you'd like to interpolate.

## Construct from points

For example, to do linear interpolation on three points, we construct a [`LinearInterpolator`](@ref) then call it on the interpolation coordinate:
```@example
using BasicInterpolators: LinearInterpolator
x = [0.0, 1, 2];
y = [0, 0.5, 1];
p = LinearInterpolator(x, y);
p(0.5)
```
Pretty basic! A nice result of this syntax is easy broadcasting. To interpolate a bunch of values in some array `x`, you can simply use `p.(x)`.

## Construct from a function

Each of the types can be constructed by passing vectors/arrays with the points to interpolate. They can also be constructed with a function to sample from. This is probably most useful for the Chebyshev interpolators, which require a specific grid spacing. For example, to generate a [`ChebyshevInterpolator`](@ref) of the function $f(x) = \sin(2x) + x^2\exp(x)$, in the range $[-3,0]$, using 16 points:
```@example
using BasicInterpolators: ChebyshevInterpolator
f(x) = sin(2x) + x^2*exp(x);
p = ChebyshevInterpolator(f, -3, 0, 16);
(f(-1.5) - p(-1.5))/f(-1.5)
```
In the background, the constructor generates 16 points between -3 and 0, evaluates the function at those points, then sets up everything it needs to interpolate in that range.

## Two-dimensional interpolation

The two-dimensional interpolators work the same way, but with another dimension. You can supply the coordinates along each axis and a grid of points or you can supply a function along with where to evaluate it and how many points to use. For example, to make a `BilinearInterpolator` for a 10 by 10 grid of randomly spaced points with random values:
```@example
using BasicInterpolators: BilinearInterpolator
x = sort(rand(10));
y = sort(rand(10));
A = rand(10,10);
P = BilinearInterpolator(x, y, A);
P(0.5, 0.5)
```

## Grid spacing

Unevenly spaced grids are acceptable for all interpolators except the chebyshev interpolators, but you have to supply the coordinates directly instead of constructing with a function. When you do construct with a function, all the interpolators use evenly spaced points except for the chebyshev interpolators.

## Boundaries

The [chebyshev interpolators](chebyshev.md) will never let you pass a point outside the stated interpolation range (extrapolate), because they are undefined outside the boundaries.

For all the other interpolators, there are three options.

1. [`StrictBoundaries`](@ref) - error for any point outside the interpolation points, even *barely* outside them
2. [`WeakBoundaries`](@ref) - error only for points that are more than barely outside the boundaries, as evaluated by `≈`. For example, if `x` is the interpolation coordinate and `xa` is the lower boundary coordinate, an error is triggered only if `(x < xa) & !(x ≈ xa)`.
3. [`NoBoundaries`](@ref) - no boundary checks at all and no errors

To select one of these behaviors, pass one of the objects as the last argument of the constructor. The default behavior is `StrictBoundaries`.

```@repl
using BasicInterpolators
x = [0.0, 1, 2];
y = rand(3);
p = CubicSplineInterpolator(x, y);
#trying to extrapolate will normally cause an error
p(2 + 1e-10)
#unless the boundary behavior is changed
p = CubicSplineInterpolator(x, y, WeakBoundaries());
p(2 + 1e-10)
#but going way outside the boundary is only ok with NoBoundaries
p(3)
p = CubicSplineInterpolator(x, y, NoBoundaries());
p(3)
```
Extrapolation blindly uses the interpolating curve/surface of the nearest boundary cell in the grid, so it can easily return wild values, especially for the cubic interpolators. That's why it's not allowed by default. All of the above applies to the two-dimensional interpolators as well.

## Viewing/Changing Data

The interpolators store a copy of the coordinates and values. When possible, you can get set the values inside an interpolator with regular indexing (`getindex` and `setindex!`).

```@repl
using BasicInterpolators
p = LinearInterpolator([0, 1], [0, 1]);
p(0.5)
p[1]
p[1] = 2
p[1]
p(0.5)
```