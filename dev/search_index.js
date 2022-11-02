var documenterSearchIndex = {"docs":
[{"location":"scattered/","page":"Scattered, N Dimensions","title":"Scattered, N Dimensions","text":"BasicInterpolators implements two types for interpolating scattered points in any number of dimensions.","category":"page"},{"location":"scattered/","page":"Scattered, N Dimensions","title":"Scattered, N Dimensions","text":"radial basis function (RBF) interpolation using any radial basis function in the form ϕ(rϵ), where r is the distance between points and ϵ is a scaling factor that can be used in any way. Several common functions are provided, but using your own function will not slow down the interpolator unless the function itself is slow. Constructing an RBF interpolator for p points requires the solution of a symmetric p  p matrix equation, which will require a noticeable amount of time and memory with a lot of points. Three common radial basis functions are available for convenience:\ninvmultiquadratic (default)\nmultiquadratic\ngaussian\nShepard interpolation, which is a special case of RBF interpolation where the radial basis function is ϕ(r) = r^-a, where r is the distance between points and a is a user-specified exponent. Shepard interpolation is usually less accurate than RBF interpolation, but it requires negligible construction time and memory.","category":"page"},{"location":"scattered/","page":"Scattered, N Dimensions","title":"Scattered, N Dimensions","text":"For p points in N dimensions, construct these two interpolators with a pN array of coordinates and a length p vector of values to interpolate. For example, constructing an RBFInterpolator for 100 random points in three dimensions (with ϵ=1) then interpolating at a point would look like","category":"page"},{"location":"scattered/","page":"Scattered, N Dimensions","title":"Scattered, N Dimensions","text":"using BasicInterpolators: RBFInterpolator\nX = rand(100,3);\ny = rand(100);\nP = RBFInterpolator(X, y, 1);\nP(0.5, 0.5, 0.5)","category":"page"},{"location":"scattered/","page":"Scattered, N Dimensions","title":"Scattered, N Dimensions","text":"The invmultiquadratic RBF is used here by default. To use any other RBF, pass a function in the form ϕ(rϵ) to the constructor after the scale factor, like P = RBFInterpolator(X, y, 1, myrbf).","category":"page"},{"location":"scattered/","page":"Scattered, N Dimensions","title":"Scattered, N Dimensions","text":"The interpolators can be called with an array of coordinates, like P([1,2,3]), or with individual coordinates, like P(1,2,3). Either way is fine.","category":"page"},{"location":"scattered/","page":"Scattered, N Dimensions","title":"Scattered, N Dimensions","text":"There is no perfect rule for choosing the scale factor, ϵ. For the functions provided here, it should generally be larger than the spacing between points but smaller than variations in the interpolated points. For example, if you were interpolating the function f(xy) = sin(x) + sin(y), you might choose ϵ1 because that's a little smaller than the width of the bumps in a sine curve. The best results may require some tinkering.","category":"page"},{"location":"scattered/","page":"Scattered, N Dimensions","title":"Scattered, N Dimensions","text":"","category":"page"},{"location":"scattered/","page":"Scattered, N Dimensions","title":"Scattered, N Dimensions","text":"RBFInterpolator\ninvmultiquadratic\nmultiquadratic\ngaussian","category":"page"},{"location":"scattered/#BasicInterpolators.RBFInterpolator","page":"Scattered, N Dimensions","title":"BasicInterpolators.RBFInterpolator","text":"RBFInterpolator(X, y, ϵ, rbf=invmultiquadratic)\n\nConstruct a radial basis function (RBF) interpolator for an n-dimensional set of points with coordinates X and values y. X must be an p × N array, where p is the number of points and N is the number of dimensions. y must be a length p vector. The value of ϵ scales the radial basis function of choice, f, which is invmultiquadratic by default. Any function in the form ϕ(rϵ) can be passed to the rbf argument, where r is the distance between points and ϵ is a scaling factor.\n\n\n\n\n\n","category":"type"},{"location":"scattered/#BasicInterpolators.invmultiquadratic","page":"Scattered, N Dimensions","title":"BasicInterpolators.invmultiquadratic","text":"invmultiquadratic(r, ϵ)\n\nfrac1sqrt1 + (rϵ)^2\n\n\n\n\n\n","category":"function"},{"location":"scattered/#BasicInterpolators.multiquadratic","page":"Scattered, N Dimensions","title":"BasicInterpolators.multiquadratic","text":"multiquadratic(r, ϵ)\n\nsqrt1 + (rϵ)^2\n\n\n\n\n\n","category":"function"},{"location":"scattered/#BasicInterpolators.gaussian","page":"Scattered, N Dimensions","title":"BasicInterpolators.gaussian","text":"gaussian(r, ϵ)\n\ne^-r^2ϵ^2\n\n\n\n\n\n","category":"function"},{"location":"scattered/","page":"Scattered, N Dimensions","title":"Scattered, N Dimensions","text":"","category":"page"},{"location":"scattered/","page":"Scattered, N Dimensions","title":"Scattered, N Dimensions","text":"ShepardInterpolator","category":"page"},{"location":"scattered/#BasicInterpolators.ShepardInterpolator","page":"Scattered, N Dimensions","title":"BasicInterpolators.ShepardInterpolator","text":"ShepardInterpolator(X, y, a=3)\n\nConstruct a ShepardInterpolator for an n-dimensional set of points with coordinates X and values y. X must be an p × N array, where p is the number of points and N is the number of dimensions. y must be a length p vector. The value of a defines the distance weighting function r^-a.\n\n\n\n\n\n","category":"type"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"BasicInterpolators offers callable types, or function-like objects, for your interpolation needs. This simply means that you construct an Interpolator object, like a LinearInterpolator or a BicubicInterpolator, then call it like a function, passing it the coordinate where you'd like to interpolate.","category":"page"},{"location":"tutorial/#Construct-from-points","page":"Tutorial","title":"Construct from points","text":"","category":"section"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"For example, to do linear interpolation on three points, we construct a LinearInterpolator then call it on the interpolation coordinate:","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"using BasicInterpolators: LinearInterpolator\nx = [0.0, 1, 2];\ny = [0, 0.5, 1];\np = LinearInterpolator(x, y);\np(0.5)","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"Pretty basic! A nice result of this syntax is easy broadcasting. To interpolate a bunch of values in some array x, you can simply use p.(x).","category":"page"},{"location":"tutorial/#Construct-from-a-function","page":"Tutorial","title":"Construct from a function","text":"","category":"section"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"Each of the types can be constructed by passing vectors/arrays with the points to interpolate. They can also be constructed with a function to sample from. This is probably most useful for the Chebyshev interpolators, which require a specific grid spacing. For example, to generate a ChebyshevInterpolator of the function f(x) = sin(2x) + x^2exp(x), in the range -30, using 16 points:","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"using BasicInterpolators: ChebyshevInterpolator\nf(x) = sin(2x) + x^2*exp(x);\np = ChebyshevInterpolator(f, -3, 0, 16);\n(f(-1.5) - p(-1.5))/f(-1.5)","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"In the background, the constructor generates 16 points between -3 and 0, evaluates the function at those points, then sets up everything it needs to interpolate in that range.","category":"page"},{"location":"tutorial/#Two-dimensional-interpolation","page":"Tutorial","title":"Two-dimensional interpolation","text":"","category":"section"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"The two-dimensional interpolators work the same way, but with another dimension. You can supply the coordinates along each axis and a grid of points or you can supply a function along with where to evaluate it and how many points to use. For example, to make a BilinearInterpolator for a 10 by 10 grid of randomly spaced points with random values:","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"using BasicInterpolators: BilinearInterpolator\nx = sort(rand(10));\ny = sort(rand(10));\nA = rand(10,10);\nP = BilinearInterpolator(x, y, A);\nP(0.5, 0.5)","category":"page"},{"location":"tutorial/#Grid-spacing","page":"Tutorial","title":"Grid spacing","text":"","category":"section"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"Unevenly spaced grids are acceptable for all interpolators except the chebyshev interpolators and the bivariate spline, but you have to supply the coordinates directly instead of constructing with a function. When you do construct with a function, all the interpolators use evenly spaced points except for the chebyshev interpolators.","category":"page"},{"location":"tutorial/#Boundaries","page":"Tutorial","title":"Boundaries","text":"","category":"section"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"The chebyshev interpolators will never let you pass a point outside the stated interpolation range (extrapolate), because they are undefined outside the boundaries.","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"For all the other interpolators, there are three options.","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"StrictBoundaries - error for any point outside the interpolation points, even barely outside them\nWeakBoundaries - error only for points that are more than barely outside the boundaries, as evaluated by ≈. For example, if x is the interpolation coordinate and xa is the lower boundary coordinate, an error is triggered only if (x < xa) & !(x ≈ xa).\nNoBoundaries - no boundary checks at all and no errors","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"To select one of these behaviors, pass one of the objects as the last argument of the constructor. The default behavior is StrictBoundaries.","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"using BasicInterpolators\nx = [0.0, 1, 2];\ny = rand(3);\np = CubicSplineInterpolator(x, y);\n#trying to extrapolate will normally cause an error\np(2 + 1e-10)\n#unless the boundary behavior is changed\np = CubicSplineInterpolator(x, y, WeakBoundaries());\np(2 + 1e-10)\n#but going way outside the boundary is only ok with NoBoundaries\np(3)\np = CubicSplineInterpolator(x, y, NoBoundaries());\np(3)","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"Extrapolation blindly uses the interpolating curve/surface of the nearest boundary cell in the grid, so it can easily return wild values, especially for the cubic interpolators. That's why it's not allowed by default. All of the above applies to the two-dimensional interpolators as well.","category":"page"},{"location":"tutorial/#Viewing/Changing-Data","page":"Tutorial","title":"Viewing/Changing Data","text":"","category":"section"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"The interpolators store a copy of the coordinates and values. When possible, you can get set the values inside an interpolator with regular indexing (getindex and setindex!).","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"using BasicInterpolators\np = LinearInterpolator([0, 1], [0, 1]);\np(0.5)\np[1]\np[1] = 2\np[1]\np(0.5)","category":"page"},{"location":"1d/","page":"One Dimension","title":"One Dimension","text":"Here is a table of the types/constructors available for one-dimensional interpolation.","category":"page"},{"location":"1d/","page":"One Dimension","title":"One Dimension","text":"Interpolators Method\nLinearInterpolator piecewise linear\nCubicInterpolator piecewise cubic (no smoothness guarantee)\nCubicSplineInterpolator cubic spline, natural or clamped\nChebyshevInterpolator Chebyshev expansion","category":"page"},{"location":"1d/","page":"One Dimension","title":"One Dimension","text":"For the derivative of a ChebyshevInterpolator, use chebyderiv.","category":"page"},{"location":"1d/","page":"One Dimension","title":"One Dimension","text":"","category":"page"},{"location":"1d/","page":"One Dimension","title":"One Dimension","text":"Here is a table of the functions available for one-dimensional interpolation.","category":"page"},{"location":"1d/","page":"One Dimension","title":"One Dimension","text":"Functions Use\nquadratic quadratic interpolation of any 3 points\ncubic cubic interpolation of any 4 points\nneville (n-1)th order polynomial interpolation of any n points\nvandermonde coefficients of (n-1)th order polynomial passing through n points\ncubichermite cubic interpolation using two points with first derivatives","category":"page"},{"location":"1d/","page":"One Dimension","title":"One Dimension","text":"","category":"page"},{"location":"1d/","page":"One Dimension","title":"One Dimension","text":"LinearInterpolator\nCubicInterpolator\nCubicSplineInterpolator\nChebyshevInterpolator\nquadratic\ncubic\nneville\nvandermonde\ncubichermite","category":"page"},{"location":"1d/#BasicInterpolators.LinearInterpolator","page":"One Dimension","title":"BasicInterpolators.LinearInterpolator","text":"LinearInterpolator(x, y, boundaries=StrictBoundaries())\n\nConstruct a LinearInterpolator for the points defined by coordinates x and values y\n\n\n\n\n\nLinearInterpolator(f, xa, xb, n, boundaries=StrictBoundaries())\n\nConstruct a LinearInterpolator for the function f using n evenly spaced function evaluations in the range [xa,xb]\n\n\n\n\n\n","category":"type"},{"location":"1d/#BasicInterpolators.CubicInterpolator","page":"One Dimension","title":"BasicInterpolators.CubicInterpolator","text":"CubicInterpolator(x, y, boundaries=StrictBoundaries())\n\nConstruct a CubicInterpolator for the points defined by coordinates x and values y\n\n\n\n\n\nCubicInterpolator(f, xa, xb, n, boundaries=StrictBoundaries())\n\nConstruct a CubicInterpolator for the function f using n evenly spaced function evaluations in the range [xa,xb]\n\n\n\n\n\n","category":"type"},{"location":"1d/#BasicInterpolators.CubicSplineInterpolator","page":"One Dimension","title":"BasicInterpolators.CubicSplineInterpolator","text":"CubicSplineInterpolator(x, y, boundaries=StrictBoundaries())\n\nConstruct a CubicSplineInterpolator for the points defined by coordinates x and values y. This constructor creates a natural spline, where the second derivative is set to zero at the boundaries.\n\n\n\n\n\nCubicSplineInterpolator(x, y, dy₁, dyₙ, boundaries=StrictBoundaries())\n\nConstruct a CubicSplineInterpolator for the points defined by coordinates x and values y. This constructor creates a clamped spline, where the first derivatives at the boundaries are set by dy₁ and dyₙ.\n\n\n\n\n\nCubicSplineInterpolator(f, xa, xb, n, boundaries=StrictBoundaries())\n\nConstruct a CubicSplineInterpolator for the function f using n evenly spaced function evaluations in the range [xa,xb]. A natural spline is created.\n\n\n\n\n\n","category":"type"},{"location":"1d/#BasicInterpolators.ChebyshevInterpolator","page":"One Dimension","title":"BasicInterpolators.ChebyshevInterpolator","text":"ChebyshevInterpolator(x, y)\n\nConstruct a ChebyshevInterpolator for the points defined by coordinates x and values y. The x coordinates must be arranged on a chebyshev grid, which can be generated using the chebygrid function.\n\n\n\n\n\nChebyshevInterpolator(f, xa, xb, n)\n\nConstruct a ChebyshevInterpolator for the function f using n function evaluations in the range [xa,xb]. The function evaluations will occur on the chebyshev nodes.\n\n\n\n\n\n","category":"type"},{"location":"1d/#BasicInterpolators.quadratic","page":"One Dimension","title":"BasicInterpolators.quadratic","text":"quadratic(x, xₚ, yₚ)\n\nPerform quadratic polynomial interpolation of the points defined by coordinates xₚ and values yₚ, at the coordinate x, using Neville's algorithm. xₚ and yₚ must both contain three points.\n\n\n\n\n\n","category":"function"},{"location":"1d/#BasicInterpolators.cubic","page":"One Dimension","title":"BasicInterpolators.cubic","text":"cubic(x, xₚ, yₚ)\n\nPerform cubic polynomial interpolation of the points defined by coordinates xₚ and values yₚ, at the coordinate x, using Neville's algorithm. xₚ and yₚ must both contain four points.\n\n\n\n\n\n","category":"function"},{"location":"1d/#BasicInterpolators.neville","page":"One Dimension","title":"BasicInterpolators.neville","text":"neville(x, xₚ, yₚ)\n\nPerform polynomial interpolation of the points defined by coordinates xₚ and values yₚ, at the coordinate x, using Neville's algorithm with as many points as are provided. xₚ and yₚ must have the same length. With only 3 or 4 points the quadratic and cubic functions will be considerably faster.\n\n\n\n\n\n","category":"function"},{"location":"1d/#BasicInterpolators.vandermonde","page":"One Dimension","title":"BasicInterpolators.vandermonde","text":"vandermonde(x, y)\n\nGenerate the coefficients of an arbitrary order polynomial passing through the ponts defined by coordinates x and value y. For n points, n coefficients c_0 c_1  c_n-1 are returned forming the polynomial c_0 + c_1x +  + c_n-1x^n-1\n\nwarning: Warning\nSolving for the the coefficients of a high-order polynomial is a notoriously ill-conditioned problem. It is not recommended for orders greater than 5 or 6, although it depends on the application. If you must interpolate with a high-order polynomial, it's better to use the neville function instead of computing coefficients.\n\n\n\n\n\n","category":"function"},{"location":"1d/#BasicInterpolators.cubichermite","page":"One Dimension","title":"BasicInterpolators.cubichermite","text":"cubichermite(x, x₁, x₂, y₁, y₂, y₁′, y₂′)\n\nInterpolate a cubic polynomial between two points, given its values and first derivatives at the points.\n\n\n\n\n\n","category":"function"},{"location":"misc/","page":"Miscellaneous","title":"Miscellaneous","text":"Pieces of the package that might be useful 🤷","category":"page"},{"location":"misc/","page":"Miscellaneous","title":"Miscellaneous","text":"findcell\ncheby\nchebycoef\nchebygrid\nchebyderiv","category":"page"},{"location":"misc/#BasicInterpolators.findcell","page":"Miscellaneous","title":"BasicInterpolators.findcell","text":"findcell(q, V, n)\n\nUse bisection search to find the cell containing q, assuming V is a sorted vector of n coordinates. The returned integer is the index of the element in V immediately less than q. For example, if findcell returns 2, then q ∈ [V[2],V[3]). If q is less than every element in V, 1 is returned, indicating the first cell in V. If q is greater than every element in V, length(V)-1 is returned, indicating the last cell in V.\n\n\n\n\n\n","category":"function"},{"location":"misc/#BasicInterpolators.cheby","page":"Miscellaneous","title":"BasicInterpolators.cheby","text":"cheby(coef, x, xa, xb)\n\nEvaluates the Chebyshev expansion represented by the coefficients in coef and defined on the interval [xa,xb] at the point x.\n\n\n\n\n\n","category":"function"},{"location":"misc/#BasicInterpolators.chebycoef","page":"Miscellaneous","title":"BasicInterpolators.chebycoef","text":"chebycoef(y)\n\nCompute the Chebyshev expansion coefficients for a set of points y, which are assumed to be located on the Chebyshev points for some interval.\n\n\n\n\n\n","category":"function"},{"location":"misc/#BasicInterpolators.chebygrid","page":"Miscellaneous","title":"BasicInterpolators.chebygrid","text":"chebygrid(n)\n\nCreate an array of n chebyshev nodes in [-1,1]\n\n\n\n\n\nchebygrid(xa, xb, n)\n\nCreate an array of n chebyshev nodes in [xa,xb]\n\n\n\n\n\nchebygrid(xa, xb, nx, ya, yb, ny)\n\nCreate a two-dimensional grid of chebyshev nodes using nx points along the first axis, in [xa,xb], and ny points along the second axis, in [ya,yb].\n\n\n\n\n\n","category":"function"},{"location":"misc/#BasicInterpolators.chebyderiv","page":"Miscellaneous","title":"BasicInterpolators.chebyderiv","text":"chebyderiv(coef, xa, xb)\n\nGenerates the expansion coefficents for the derivative of a preexisting Chebyshev expansion defined on the interval [xa,xb].\n\n\n\n\n\nchebyderiv(ϕ::ChebyshevInterpolator)\n\nConstruct a ChebyshevInterpolator representing the derivative of a preexisting interpolator.\n\n\n\n\n\n","category":"function"},{"location":"chebyshev/","page":"-","title":"-","text":"Chebyshev interpolation is an extremely accurate method for interpolating smooth functions. It's a good way to replace a very slow/costly function with a high accuracy approximation (if the function is smooth!). The extremely high accuracy is achieved by careful placement of the interpolation nodes at the Chebyshev nodes. This part of BasicInterpolators that is not really so basic, but chebyshev interpolation is a powerful tool in the right context and it's worth providing easy access to it.","category":"page"},{"location":"chebyshev/","page":"-","title":"-","text":"To easily create a ChebyshevInterpolator or BichebyshevInterpolator without dealing with the specific grid spacing they require, pass a function to the constructor. The function will be evaluated at the appropriate locations to generate an interpolating approximation. For example, for two-dimensional interpolation of f(xy) = sin(x) + cos(y) with xy  -33","category":"page"},{"location":"chebyshev/","page":"-","title":"-","text":"using BasicInterpolators: BichebyshevInterpolator\nf(x,y) = sin(x) + cos(y)\nP = BichebyshevInterpolator(f, -3, 3, 16, -3, 3, 16);\nP(1, 2) - f(1, 2)","category":"page"},{"location":"chebyshev/","page":"-","title":"-","text":"With 16 points in each dimension, the error is already about one part per trillion. This is a contrived example and the incredible accuracy of Chebyshev interpolation can be demonstrated more thoroughly, but the important fact is that Chebyshev approximations to smooth functions have exponential, or \"infinite order\" convergence. As soon as there are enough points for the Chebyshev approximation to capture the function's basic shape, the error will plummet to machine precision as more points are added. For more on this and related topics, I recommend the following book","category":"page"},{"location":"chebyshev/","page":"-","title":"-","text":"Boyd, John P. Chebyshev and Fourier spectral methods. Courier Corporation, 2001.","category":"page"},{"location":"chebyshev/","page":"-","title":"-","text":"The Chebyshev interpolators do some things in the background that the other interpolators don't. Construction requires some matrix math, but constructing many interpolators with the same dimensions will be fast because inverted matrices are cached. In the two-dimensional case, evaluation uses fast in-place linear algebra functions when possible (mul!). Because of the in-place operations, a single BichebyshevInterpolator interpolator is generally not thread-safe.","category":"page"},{"location":"chebyshev/","page":"-","title":"-","text":"It's easy to evaluate the derivative of a Chebyshev expansion. See chebyderiv.","category":"page"},{"location":"2d/","page":"Two Dimensions","title":"Two Dimensions","text":"The following interpolator types all operate on regular grids.","category":"page"},{"location":"2d/","page":"Two Dimensions","title":"Two Dimensions","text":"The BilinearInterpolator and BicubicInterpolator will work with any regular, unevenly spaced grids.\nThe BicubicSplineInterpolator expects uniform spacing in each dimension, although the spacing on one axis doesn't have to match the other. This interpolator estimates derivatives using finite differences to construct the splines.\nThe BichebyshevInterpolator must have the usual Chebyshev spacing in each direction, but it can have a different number of points on each axis.","category":"page"},{"location":"2d/","page":"Two Dimensions","title":"Two Dimensions","text":"","category":"page"},{"location":"2d/","page":"Two Dimensions","title":"Two Dimensions","text":"Interpolators Method\nBilinearInterpolator piecewise linear\nBicubicInterpolator piecewise cubic (no smoothness guarantee)\nBichebyshevInterpolator Chebyshev expansions","category":"page"},{"location":"2d/","page":"Two Dimensions","title":"Two Dimensions","text":"","category":"page"},{"location":"2d/","page":"Two Dimensions","title":"Two Dimensions","text":"BilinearInterpolator\nBicubicInterpolator\nBichebyshevInterpolator","category":"page"},{"location":"2d/#BasicInterpolators.BilinearInterpolator","page":"Two Dimensions","title":"BasicInterpolators.BilinearInterpolator","text":"BilinearInterpolator(x, y, Z, boundaries=StrictBoundaries())\n\nConstruct a BilinearInterpolator for the grid of points points defined by coordinates x,y and values Z.\n\n\n\n\n\nBilinearInterpolator(f, xa, xb, nx, ya, yb, ny, boundaries=StrictBoundaries())\n\nConstruct a BilinearInterpolator for the function f using a grid of nx points evenly spaced on the first axis in [xa,xb] and ny points evenly spaced on the second axis in [ya,yb].\n\n\n\n\n\n","category":"type"},{"location":"2d/#BasicInterpolators.BicubicInterpolator","page":"Two Dimensions","title":"BasicInterpolators.BicubicInterpolator","text":"BicubicInterpolator(x, y, Z, boundaries=StrictBoundaries())\n\nConstruct a BicubicInterpolator for the grid of points points defined by coordinates x,y and values Z.\n\n\n\n\n\nBicubicInterpolator(f, xa, xb, nx, ya, yb, ny, boundaries=StrictBoundaries())\n\nConstruct a BicubicInterpolator for the function f using a grid of nx points evenly spaced on the first axis in [xa,xb] and ny points evenly spaced on the second axis in [ya,yb].\n\n\n\n\n\n","category":"type"},{"location":"2d/#BasicInterpolators.BichebyshevInterpolator","page":"Two Dimensions","title":"BasicInterpolators.BichebyshevInterpolator","text":"BichebyshevInterpolator(x, y, Z)\n\nConstruct a BichebyshevInterpolator for the grid of points defined by coordinates (x,y) and values Z. The given points must lie on a chebyshev grid in each direction. These can be generated with the chebygrid function or the interpolator can be constructed directly from a function using the method below.\n\nwarning: Warning\nThe Bichebyshev interpolator is not thread-safe. It computes a cosine expansion and does some linear algebra in-place using arrays stored with the object. A single BichebyshevInterpolator should never be called by multiple threads at once.\n\n\n\n\n\nBichebyshevInterpolator(f, xa, xb, nx, ya, yb, ny)\n\nConstruct a BichebyshevInterpolator for the function f using a grid of nx points on the first axis in [xa,xb] and ny points on the second axis in [ya,yb].\n\n\n\n\n\n","category":"type"},{"location":"#BasicInterpolators.jl","page":"Home","title":"BasicInterpolators.jl","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Welcome to the BasicInterpolators documentation. This package provides a collection of common interpolation recipes. Quick facts about the package:","category":"page"},{"location":"","page":"Home","title":"Home","text":"The Interpolator types are the main focus. To use them, you construct an interpolator object then call it like a function. There are also some functions for things like one-off polynomial interpolation.\nInterpolators can be constructed with any numeric type that supports the usual arithmetic operations necessary for interpolating. They will, however, insist that the element types of your grid(s) and values are uniform, often converting to consistent element types automatically.\nThere are two-dimensional grid interpolators, but not higher.\nAll the interpolators should be compatible with automatic differentiation, specifically ForwardDiff.jl.","category":"page"},{"location":"","page":"Home","title":"Home","text":"","category":"page"},{"location":"#Details","page":"Home","title":"Details","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"For a tutorial, see the tutorial page. For more information on the specific types and methods available, see these links:","category":"page"},{"location":"","page":"Home","title":"Home","text":"One-dimensional interpolation\nTwo-dimensional interpolation\nN-dimensional, scattered point interpolation","category":"page"},{"location":"","page":"Home","title":"Home","text":"For a little extra on the Chebyshev interpolators, look here.","category":"page"},{"location":"","page":"Home","title":"Home","text":"","category":"page"},{"location":"#Other-Packages","page":"Home","title":"Other Packages","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"If you need more features/methods/dimensions, please look in other packages:","category":"page"},{"location":"","page":"Home","title":"Home","text":"Interpolations.jl\nDierckx.jl\nGridInterpolations.jl\nApproXD\nFastChebInterp.jl","category":"page"},{"location":"boundaries/","page":"Boundary Behavior","title":"Boundary Behavior","text":"The following types can be passed to interpolator constructors to modify how points are handled at/outside the boundaries of the interpolation nodes. See the \"Boundaries\" section of the tutorial for examples.","category":"page"},{"location":"boundaries/","page":"Boundary Behavior","title":"Boundary Behavior","text":"StrictBoundaries\nWeakBoundaries\nNoBoundaries","category":"page"},{"location":"boundaries/#BasicInterpolators.StrictBoundaries","page":"Boundary Behavior","title":"BasicInterpolators.StrictBoundaries","text":"StrictBoundaries()\n\nTriggers errors whenever interpolation coordinates are outside of the boundaries.\n\n\n\n\n\n","category":"type"},{"location":"boundaries/#BasicInterpolators.WeakBoundaries","page":"Boundary Behavior","title":"BasicInterpolators.WeakBoundaries","text":"WeakBoundaries()\n\nAllows small overshoots at boundaries, but not large ones. Errors are only triggered when the interpolation coordinate is outside of a boundary and not close to it. At the lower boundary, for example, an error would be triggered when (x < boundary) & !(x ≈ boundary).\n\n\n\n\n\n","category":"type"},{"location":"boundaries/#BasicInterpolators.NoBoundaries","page":"Boundary Behavior","title":"BasicInterpolators.NoBoundaries","text":"NoBoundaries()\n\nPerforms no boundary checking. Allows interpolators to blindly extrapolate when possible.\n\n\n\n\n\n","category":"type"}]
}
