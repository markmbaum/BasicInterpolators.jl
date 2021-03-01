# Chebyshev Interpolation

Chebyshev interpolation is an extremely accurate method for interpolating smooth functions.  However, the extremely high accuracy is achieved by careful placement of the interpolation nodes at the [Chebyshev nodes](https://en.wikipedia.org/wiki/Chebyshev_nodes). Chebyshev interpolation is a good way to approximate a smooth function that is very expensive to evaluate. By evaluating the function at the Chebyshev nodes, which are given by the [`chebygrid`](@ref) function, highly accurate and efficient interpolators can be generated to speed things up. This is the only part of `BasicInterpolators` that is not really so basic, but chebyshev interpolation is a powerful tool in the right context and it's worth providing easy access to it.

To create a `ChebyshevInterpolator` or `BichebyshevInterpolator`, it is probably easiest to pass a function to the constructor, which will evaluate the function at the Chebyshev nodes and generate interpolating coefficients. For example, for two-dimensional interpolation of $f(x,y) = \sin(x) + \cos(y)$ with $x,y ∈ [-3,3]$
```@repl
using BasicInterpolators: BichebyshevInterpolator
f(x,y) = sin(x) + cos(y)
P = BichebyshevInterpolator(f, -3, 3, 16, -3, 3, 16);
P(1, 2) - f(1, 2)
```
With 16 points in each dimension, the error is already about one part per trillion. This is a contrived example and the incredible accuracy of Chebyshev interpolation can be demonstrated more thoroughly, but the important fact is that Chebyshev approximations to smooth functions have exponential, or "infinite order" convergence. As soon as there are enough points for the Chebyshev approximation to capture the function's basic shape, the error will plummet to machine precision as more points are added. For more on this and related topics, I recommend the following book
* Boyd, John P. *Chebyshev and Fourier spectral methods*. Courier Corporation, 2001.

--------------------------------------------------------------------------------

!!! warning

    The Chebyshev interpolators are **not thread-safe**. They compute cosine expansions and, in the 2D case, do some linear algebra in-place using arrays stored with the objects. A single chebyshev interpolator should never be called by multiple threads at once.

As the warning above indicates, the Chebyshev interpolators do some things in the background that the other interpolators don't. I put some effort into making them fast (at least once constructed) by avoiding as many direct `cos` evaluations as possible and doing the unavoidable calculations in-place using fast linear algebra functions ([`mul!`](https://docs.julialang.org/en/v1/stdlib/LinearAlgebra/#LinearAlgebra.mul!)). I think the speed of the Chebyshev interpolators here should compare well with any other implementations, but they should be considerably easier to use.

--------------------------------------------------------------------------------

```@docs
chebygrid
ChebyshevInterpolator
BichebyshevInterpolator
```
