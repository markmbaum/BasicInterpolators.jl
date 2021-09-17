Chebyshev interpolation is an extremely accurate method for interpolating smooth functions. It's a good way to replace a very slow/costly function with a high accuracy approximation (if the function is smooth!). The extremely high accuracy is achieved by careful placement of the interpolation nodes at the [Chebyshev nodes](https://en.wikipedia.org/wiki/Chebyshev_nodes). This part of `BasicInterpolators` that is not really so basic, but chebyshev interpolation is a powerful tool in the right context and it's worth providing easy access to it.

To easily create a [`ChebyshevInterpolator`](@ref) or [`BichebyshevInterpolator`](@ref) without dealing with the specific grid spacing they require, pass a function to the constructor. The function will be evaluated at the appropriate locations to generate an interpolating approximation. For example, for two-dimensional interpolation of $f(x,y) = \sin(x) + \cos(y)$ with $x,y ∈ [-3,3]$
```@repl
using BasicInterpolators: BichebyshevInterpolator
f(x,y) = sin(x) + cos(y)
P = BichebyshevInterpolator(f, -3, 3, 16, -3, 3, 16);
P(1, 2) - f(1, 2)
```
With 16 points in each dimension, the error is already about one part per trillion. This is a contrived example and the incredible accuracy of Chebyshev interpolation can be demonstrated more thoroughly, but the important fact is that Chebyshev approximations to smooth functions have exponential, or "infinite order" convergence. As soon as there are enough points for the Chebyshev approximation to capture the function's basic shape, the error will plummet to machine precision as more points are added. For more on this and related topics, I recommend the following book
* Boyd, John P. *Chebyshev and Fourier spectral methods*. Courier Corporation, 2001.

The Chebyshev interpolators do some things in the background that the other interpolators don't. Construction requires some matrix math, but constructing many interpolators with the same dimensions will be fast because inverted matrices are cached. In the two-dimensional case, evaluation uses fast in-place linear algebra functions where appropriate ([`mul!`](https://docs.julialang.org/en/v1/stdlib/LinearAlgebra/#LinearAlgebra.mul!)). Because of the in-place operations, a single [`BichebyshevInterpolator`](@ref) interpolator is not thread-safe.
