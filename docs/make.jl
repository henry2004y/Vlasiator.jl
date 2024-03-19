using Documenter, DemoCards
using Vlasiator

branch = "master"
# generate demo files
demos, postprocess_cb, demo_assets = makedemos("examples"; branch)
# if there are generated css assets, pass it to Documenter.HTML
assets = String[]
isnothing(demo_assets) || (push!(assets, demo_assets))

makedocs(;
    modules=[Vlasiator],
    authors="Hongyang Zhou <hongyang.zhou@helsinki.fi> and contributors",
    sitename="Vlasiator.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://henry2004y.github.io/Vlasiator.jl",
        assets,
        sidebar_sitename=false
    ),
    pages=[
        "Home" => "index.md",
        "User Guide" => "manual.md",
        "Examples" => demos,
        "Calling from Python" => "python.md",
        "Gallery" => "gallery.md",
        "Benchmarks" => "benchmark.md",
        "API Reference" => "internal.md",
        "Contributing" => "contributing.md",
        "FAQ" => "log.md",
    ],
)

# postprocess after makedocs
postprocess_cb()

deploydocs(;
    repo="github.com/henry2004y/Vlasiator.jl",
)