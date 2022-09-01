using Documenter, KomaMRI

makedocs(
    modules = [KomaMRI],
    sitename = "KomaMRI.jl: General MRI simulation framework",
    authors = "Carlos Castillo Passi",
    pages = [
        "About KomaMRI" => "index.md",
        "Getting Started" => "getting-started.md",
        "Simulation Examples" => "simulation-examples.md",
        "API Documentation" => "api.md"
    ]
)

deploydocs(
    repo = "github.com/cncastillo/KomaMRI.jl.git",
    # versions = ["stable" => "v^", "v#", "dev" => "dev"],
)

#julia --color=yes make.jl
