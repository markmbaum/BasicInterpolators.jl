using BasicInterpolators
using Documenter

DocMeta.setdocmeta!(BasicInterpolators, :DocTestSetup, :(using BasicInterpolators); recursive=true)

makedocs(;
    modules=[BasicInterpolators],
    authors="Mark Baum <markmbaum@protonmail.com>",
    repo="https://github.com/markmbaum/BasicInterpolators.jl/blob/{commit}{path}#{line}",
    sitename="BasicInterpolators.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://markmbaum.github.io/BasicInterpolators.jl",
        assets=String[],
    ),
    pages=[
        "Home"=>"index.md",
        "Tutorial"=>"tutorial.md",
        "One Dimension"=>"1d.md",
        "Two Dimensions"=>"2d.md",
        "Scattered, N Dimensions"=>"scattered.md",
        "Boundary Behavior"=>"boundaries.md",
        "Miscellaneous"=>"misc.md"
    ],
)

deploydocs(;
    repo="github.com/markmbaum/BasicInterpolators.jl",
    devbranch="main",
    versions=["stable"=>"v#.#"]
)
