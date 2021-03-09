# Parametric Curve Interpolation

It is sometimes useful to interpolate a set of points in $n$ dimensions that lie on a one-dimensional curve, a parametric curve. The curve is like the trajectory of a point moving through space. To accomplish this, one option is to interpolate each dimension of coordinates independently using the distance between the points as the independent variable. The [`ParametricCurveInterpolator`](@ref) does this using natural cubic splines. It's constructed from arrays representing the coordinates in each dimension, then interpolated with a normalized coordinate in [0,1]. For example, to interpolate the curve in 3D space represented by
```math
\begin{aligned}
x(t) &= t\sin(t) \\
y(t) &= \exp(t/2) \\
z(t) &= t^3
\end{aligned}
```
over $t ∈ [-2,2]$ with 25 points, we would do the following
```@repl
using BasicInterpolators: ParametricCurveInterpolator
t = collect(LinRange(-2, 2, 25));
x = t.*sin.(t);
y = exp.(t/2);
z = t.^3;
P = ParametricCurveInterpolator(x, y, z);
```

Then `P` can be called like a function on any coordinate in [0,1] to draw a smooth curve through the points. Of course, you wouldn't need to interpolate these points if you knew the underlying function, but it's an example.

Like the other interpolators, extrapolation will cause an error by default. To override the boundary check, pass `false` after the interpolation coordinate, as in `P(0.1, false)`. Each underlying spline will then extrapolate. As always, you're likely to get wild results when extrapolating far from the boundaries.

------

```@docs
ParametricCurveInterpolator
```
