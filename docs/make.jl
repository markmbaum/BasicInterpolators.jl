using BasicInterpolators
using Documenter

DocMeta.setdocmeta!(BasicInterpolators, :DocTestSetup, :(using BasicInterpolators); recursive=true)

makedocs(;
    modules=[BasicInterpolators],
    authors="Mark Baum <markmbaum@protonmail.com> and contributors",
    repo="https://github.com/wordsworthgroup/BasicInterpolators.jl/blob/{commit}{path}#{line}",
    sitename="BasicInterpolators.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://wordsworthgroup.github.io/BasicInterpolators.jl",
        assets=String[],
    ),
    pages=[
        "Home"=>"index.md",
        "Linear Interpolators"=>"linear.md",
        "Cubic Interpolators"=>"cubic.md",
        "Cubic Spline Interpolators"=>"spline.md",
        "Chebyshev Interpolators"=>"chebyshev.md",
        "General Polynomial Interpolation"=>"polynomial.md"
    ],
)

deploydocs(;
    repo="github.com/wordsworthgroup/BasicInterpolators.jl",
    devbranch="main",
    versions=["stable"=>"v^"]
)
