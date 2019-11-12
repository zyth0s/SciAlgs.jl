using Documenter, SciAlgs

makedocs(;
    modules=[SciAlgs],
    format=:html,
    pages=[
        "Home" => "index.md",
    ],
    repo="https://github.com/zyth0s/SciAlgs.jl/blob/{commit}{path}#L{line}",
    sitename="SciAlgs.jl",
    authors="Daniel Menendez Crespo, MPI CPFS, Dresden",
    assets=[],
)

deploydocs(;
    repo="github.com/zyth0s/SciAlgs.jl",
    target="build",
    julia="1.0",
    deps=nothing,
    make=nothing,
)
