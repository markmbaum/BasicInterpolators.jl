`BasicInterpolators` implements two types for interpolating scattered points in any number of dimensions.
1. [radial basis function (RBF) interpolation](https://en.wikipedia.org/wiki/Radial_basis_function_interpolation) using *any* radial basis function in the form ``ϕ(r,ϵ)``, where ``r`` is the distance between points and ``ϵ`` is a scaling factor that can be used in any way. Several common functions are provided, but using your own function will not slow down the interpolator unless the function itself is slow. Constructing an RBF interpolator for ``p`` points requires the solution of a symmetric ``p × p`` matrix equation, which will require a noticeable amount of time and memory with a lot of points. Three common radial basis functions are available for convenience:
    * [`invmultiquadratic`](@ref) (default)
    * [`multiquadratic`](@ref)
    * [`gaussian`](@ref)
2. [Shepard interpolation](https://en.wikipedia.org/wiki/Inverse_distance_weighting), which is a special case of RBF interpolation where the radial basis function is ``ϕ(r) = r^{-a}``, where ``r`` is the distance between points and ``a`` is a user-specified exponent. Shepard interpolation is usually less accurate than RBF interpolation, but it requires negligible construction time and memory.

For $p$ points in $N$ dimensions, construct these two interpolators with a $p×N$ array of coordinates and a length $p$ vector of values to interpolate. For example, constructing an `RBFInterpolator` for 100 random points in three dimensions (with ``ϵ=1``) then interpolating at a point would look like
```@repl
using BasicInterpolators: RBFInterpolator
X = rand(100,3);
y = rand(100);
P = RBFInterpolator(X, y, 1);
P(0.5, 0.5, 0.5)
```
The [`invmultiquadratic`](@ref) RBF is used here by default. To use any other RBF, pass a function in the form ``ϕ(r,ϵ)`` to the constructor after the scale factor, like `P = RBFInterpolator(X, y, 1, myrbf)`.

The interpolators can be called with an array of coordinates, like `P([1,2,3])`, or with individual coordinates, like `P(1,2,3)`. Either way is fine.

There is no perfect rule for choosing the scale factor, ``ϵ``. For the functions provided here, it should generally be larger than the spacing between points but smaller than variations in the interpolated points. For example, if you were interpolating the function ``f(x,y) = \sin(x) + \sin(y)``, you might choose ``ϵ≈1`` because that's a little smaller than the width of the bumps in a sine curve. The best results may require some tinkering.

-----

```@docs
RBFInterpolator
invmultiquadratic
multiquadratic
gaussian
```
-----
```@docs
ShepardInterpolator
```
