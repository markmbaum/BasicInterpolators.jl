# Cubic Interpolation

The types below implement piecewise cubic interpolation in one and two dimensions. These are not typical cubic splines, which are also implemented in this package. The `CubicInterpolator` and `BicubicInterpolator` do not guarantee smoothness at the interpolation nodes. However, they can be used with unevenly spaced grids and perform well when interpolating smooth curves/surfaces. For cubic splines, look [here](spline.md).

--------------------------------------------------------------------------------

```@docs
CubicInterpolator
BicubicInterpolator
```
