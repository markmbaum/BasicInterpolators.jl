export CubicSplineInterpolator,
       ParametricCurveInterpolator,
       BicubicSplineInterpolator

#-------------------------------------------------------------------------------
# piecewise cubics with continuous derivatives (splines!)

struct CubicSplineInterpolator{B}
    r::InterpolatorRange
    coef::Vector{NTuple{4,Float64}}
    boundaries::B
end

"""
    CubicSplineInterpolator(x, y, boundaries=StrictBoundaries())

Construct a `CubicSplineInterpolator` for the points defined by coordinates `x` and values `y`. This constructor creates a natural spline, where the second derivative is set to zero at the boundaries.
"""
function CubicSplineInterpolator(x::AbstractVector{<:Real},
                                 y::AbstractVector{<:Real},
                                 boundaries::AbstractBoundaries=StrictBoundaries())
    #construct the underlying range, triggering some checks
    r = InterpolatorRange(x, y)
    n = r.n
    #compute coefficients
    #Burden, Richard L., and J. Douglas Faires. Numerical Analysis. 2011.
    a = collect(Float64, r.y)
    b = zeros(n - 1)
    d = zeros(n - 1)
    h = diff(r.x)
    α = zeros(n-1)
    for i = 2:n-1
        α[i] = 3*(a[i+1] - a[i])/h[i] - 3*(a[i] - a[i-1])/h[i-1]
    end
    c = zeros(n)
    l = ones(n)
    μ = zeros(n)
    z = zeros(n)
    l[1] = 1
    for i = 2:n-1
        l[i] = 2*(x[i+1] - x[i-1]) - h[i-1]*μ[i-1]
        μ[i] = h[i]/l[i]
        z[i] = (α[i] - h[i-1]*z[i-1])/l[i]
    end
    for j = n-1:-1:1
        c[j] = z[j] - μ[j]*c[j+1]
        b[j] = (a[j+1] - a[j])/h[j] - h[j]*(c[j+1] + 2*c[j])/3
        d[j] = (c[j+1] - c[j])/(3*h[j])
    end
    a = a[1:end-1]
    c = c[1:end-1]
    #static arrays
    coef = Vector{NTuple{4,Float64}}(undef,n-1)
    for i = 1:n-1
        coef[i] = (a[i], b[i], c[i], d[i])
    end
    #construct the object
    CubicSplineInterpolator(r, coef, boundaries)
end

"""
    CubicSplineInterpolator(x, y, dy₁, dyₙ, boundaries=StrictBoundaries())

Construct a `CubicSplineInterpolator` for the points defined by coordinates `x` and values `y`. This constructor creates a clamped spline, where the first derivatives at the boundaries are set by `dy₁` and `dyₙ`.
"""
function CubicSplineInterpolator(x::AbstractVector{<:Real},
                                 y::AbstractVector{<:Real},
                                 dy₁::Real,
                                 dyₙ::Real,
                                 boundaries::AbstractBoundaries=StrictBoundaries())
    #construct the underlying range, triggering some checks
    r = InterpolatorRange(x, y)
    n = r.n
    #compute coefficients
    #Burden, Richard L., and J. Douglas Faires. Numerical Analysis. 2011.
    a = collect(Float64, r.y)
    b = zeros(n - 1)
    d = zeros(n - 1)
    h = diff(r.x)
    α = zeros(n)
    α[1] = 3*(a[2] - a[1])/h[1] - 3*dy₁
    for i = 2:n-1
        α[i] = 3*(a[i+1] - a[i])/h[i] - 3*(a[i] - a[i-1])/h[i-1]
    end
    α[n] = 3*dyₙ - 3*(a[n] - a[n - 1])/h[n-1]
    c = zeros(n)
    l = zeros(n)
    μ = zeros(n)
    z = zeros(n)
    l[1] = 2*h[1]
    μ[1] = 0.5
    z[1] = α[1]/l[1]
    for i = 2:n-1
        l[i] = 2*(x[i+1] - x[i-1]) - h[i-1]*μ[i-1]
        μ[i] = h[i]/l[i]
        z[i] = (α[i] - h[i-1]*z[i-1])/l[i]
    end
    l[n] = h[n-1]*(2 - μ[n-1])
    z[n] = (α[n] - h[n-1]*z[n-1])/l[n]
    c[n] = z[n]
    for j = n-1:-1:1
        c[j] = z[j] - μ[j]*c[j+1]
        b[j] = (a[j+1] - a[j])/h[j] - h[j]*(c[j+1] + 2*c[j])/3
        d[j] = (c[j+1] - c[j])/(3*h[j])
    end
    a = a[1:end-1]
    c = c[1:end-1]
    #static arrays
    coef = Vector{NTuple{4,Float64}}(undef,n-1)
    for i = 1:n-1
        coef[i] = (a[i], b[i], c[i], d[i])
    end
    #construct the object
    CubicSplineInterpolator(r, coef, boundaries)
end

"""
    CubicSplineInterpolator(f, xa, xb, n, boundaries=StrictBoundaries())

Construct a `CubicSplineInterpolator` for the function `f` using `n` evenly spaced function evaluations in the range [`xa`,`xb`]. A natural spline is created.
"""
function CubicSplineInterpolator(f::Function,
                                 xa::Real,
                                 xb::Real,
                                 n::Int,
                                 boundaries::AbstractBoundaries=StrictBoundaries())
    linstruct(CubicSplineInterpolator, f, xa, xb, n, boundaries)
end

function (ϕ::CubicSplineInterpolator)(x::Real)::Float64
    #enforce boundaries if desired
    ϕ.boundaries(x, ϕ.r.xa, ϕ.r.xb)
    #find the interpolation point
    i = findcell(x, ϕ.r.x)
    #offset from the nearest lower point
    ξ = x - ϕ.r.x[i]
    #evaluate polynomial
    ϕ.coef[i][1] + ϕ.coef[i][2]*ξ + ϕ.coef[i][3]*ξ^2 + ϕ.coef[i][4]*ξ^3
end

#-------------------------------------------------------------------------------
# cubic splines on a regular grid

struct BicubicSplineInterpolator{B}
    G::InterpolatorGrid
    coef::Array{NTuple{16,Float64},2}
    boundaries::B
end

"""
    BicubicSplineInterpolator(x, y, Z, boundaries=StrictBoundaries())

Construct a `BicubicSplineInterpolator` for the grid of points points defined by coordinates `x`,`y` and values `Z`.
"""
function BicubicSplineInterpolator(x::AbstractVector{<:Real},
                                   y::AbstractVector{<:Real},
                                   Z::AbstractArray{<:Real,2},
                                   boundaries::AbstractBoundaries=StrictBoundaries())
    nx, ny = size(Z)
    #insist on at least 4 points in each dimension
    @assert nx >= 3 "bicubic interpolation requires at least 3 points along axis 1"
    @assert ny >= 3 "bicubic interpolation requires at least 3 points along axis 2"
    #insist on even spacing along both axes
    @assert all(diff(diff(x)) .< 1e-8*maximum(abs.(x))) "grid spacing along axis 1 must be uniform"
    @assert all(diff(diff(y)) .< 1e-8*maximum(abs.(y))) "grid spacing along axis 2 must be uniform"
    #space for derivatives
    dx = zeros(nx, ny)
    dy = zeros(nx, ny)
    dxy = zeros(nx, ny)
    #first derivatives for internal nodes
    for i = 2:nx-1
        for j = 1:ny
            dx[i,j] = (Z[i+1,j] - Z[i-1,j])/2
        end
    end
    for i = 1:nx
        for j = 2:ny-1
            dy[i,j] = (Z[i,j+1] - Z[i,j-1])/2
        end
    end
    #first derivatives for bounary nodes
    for i = 1:nx
        dy[i,1] = -3*Z[i,1]/2 + 2*Z[i,2] - Z[i,3]/2
        dy[i,ny] = Z[i,ny-2]/2 - 2*Z[i,ny-1] + 3*Z[i,ny]/2
    end
    for j = 1:ny
        dx[1,j] = -3*Z[1,j]/2 + 2*Z[2,j] - Z[3,j]/2
        dx[nx,j] = Z[nx-2,j]/2 - 2*Z[nx-1,j] + 3*Z[nx,j]/2
    end
    #mixed second order derivatives at internal nodes
    for i = 1:nx
        for j = 2:ny-1
            dxy[i,j] = (dx[i,j+1] - dx[i,j-1])/2
        end
    end
    #mixed second deriv along the sides
    for i = 1:nx
        dxy[i,1] = -3*dx[i,1]/2 + 2*dx[i,2] - dx[i,3]/2
        dxy[i,ny] = dx[i,ny-2]/2 - 2*dx[i,ny-1] + 3*dx[i,ny]/2
    end
    #matrix needed to compute coefficients and its transpose
    A = [1. 0. 0. 0.; 0. 0. 1. 0.; -3. 3. -2. -1.; 2. -2. 1. 1.]
    Aᵀ = transpose(A)
    #space for coefficients
    f = zeros(4, 4)
    α = Array{Array{Float64,2},2}(undef, nx-1, ny-1)
    for i = 1:nx-1
        for j = 1:ny-1
            #load the matrix of 16 values and derivatives
            f[1,1] = Z[i,j]
            f[1,2] = Z[i,j+1]
            f[1,3] = dy[i,j]
            f[1,4] = dy[i,j+1]
            f[2,1] = Z[i+1,j]
            f[2,2] = Z[i+1,j+1]
            f[2,3] = dy[i+1,j]
            f[2,4] = dy[i+1,j+1]
            f[3,1] = dx[i,j]
            f[3,2] = dx[i,j+1]
            f[3,3] = dxy[i,j]
            f[3,4] = dxy[i,j+1]
            f[4,1] = dx[i+1,j]
            f[4,2] = dx[i+1,j+1]
            f[4,3] = dxy[i+1,j]
            f[4,4] = dxy[i+1,j+1]
            #get the 16 double spline coefficients
            α[i,j] = A*f*Aᵀ
        end
    end
    #static coefficients
    coef = Array{NTuple{16,Float64},2}(undef,nx-1,ny-1)
    for i = 1:nx-1
        for j = 1:ny-1
            coef[i,j] = Tuple(vec(α[i,j]))
        end
    end
    BicubicSplineInterpolator(InterpolatorGrid(x, y, Z), coef, boundaries)
end

"""
    BicubicSplineInterpolator(f, xa, xb, nx, ya, yb, ny, boundaries=StrictBoundaries())

Construct a `BicubicSplineInterpolator` for the function `f` using a grid of `nx` points evenly spaced on the first axis in [`xa`,`xb`] and `ny` points evenly spaced on the second axis in [`ya`,`yb`].
"""
function BicubicSplineInterpolator(f::Function,
                                   xa::Real, xb::Real, nx::Int,
                                   ya::Real, yb::Real, ny::Int,
                                   boundaries::AbstractBoundaries=StrictBoundaries())
    linstruct(BicubicSplineInterpolator, f, xa, xb, nx, ya, yb, ny, boundaries)
end

function (Φ::BicubicSplineInterpolator)(x::Real, y::Real)::Float64
    #enforce boundaries if desired
    Φ.boundaries(x, Φ.G.xa, Φ.G.xb, y, Φ.G.ya, Φ.G.yb)
    #find the proper grid box to interpolate inside
    i = findcell(x, Φ.G.x)
    j = findcell(y, Φ.G.y)
    #get the coefficients
    α = Φ.coef[i,j]
    #offsets
    Δx = (x - Φ.G.x[i])/(Φ.G.x[i+1] - Φ.G.x[i])
    Δy = (y - Φ.G.y[j])/(Φ.G.y[j+1] - Φ.G.y[j])
    #powers of the offsets
    Δx², Δx³ = Δx^2, Δx^3
    Δy², Δy³ = Δy^2, Δy^3
    #final interpolation calculation
    ( α[1]      + α[2]*Δx      + α[3]*Δx²      + α[4]*Δx³
    + α[5]*Δy   + α[6]*Δx*Δy   + α[7]*Δx²*Δy   + α[8]*Δx³*Δy
    + α[9]*Δy²  + α[10]*Δx*Δy² + α[11]*Δx²*Δy² + α[12]*Δx³*Δy²
    + α[13]*Δy³ + α[14]*Δx*Δy³ + α[15]*Δx²*Δy³ + α[16]*Δx³*Δy³ )
end
