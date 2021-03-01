# Cubic Spline Interpolation

The types below implement piecewise cubic *spline* interpolation in one and two dimensions. Both types guarantee smoothness at the interpolation nodes. The `CubicSplineInterpolator` is a natural spline, with the second derivative set to zero at the boundaries, and it can be used with unevenly spaced points. The `BicubicSplineInterpolator` requires points to be evenly spaced in both dimensions and patches together cubic polynomials using finite-difference approximations of the surface's derivatives at the interpolation nodes.

--------------------------------------------------------------------------------

```@docs
CubicSplineInterpolator
BicubicSplineInterpolator
```
