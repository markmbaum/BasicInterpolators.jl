# General Polynomial Interpolation in One Dimension

For general polynomial interpolation in one dimension, `BasicInterpolators` provides a few implementations of [Neville's Algorithm](https://en.wikipedia.org/wiki/Neville%27s_algorithm).
* For quadratic interpolation of any three points, use the [`quadratic`](@ref) function.
* For cubic interpolation of any four points, use the [`cubic`](@ref) function.
* For nth order polynomial interpolation of arbitrary points, the [`neville`](@ref) function is available.
All three functions use the same basic algorithm, but the quadratic and cubic versions are written out explicitly for those cases and will be faster than the general neville function.

--------------------------------------------------------------------------------

```@docs
quadratic
cubic
neville
```
