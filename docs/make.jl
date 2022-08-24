using Documenter, KomaMRI

makedocs(
    modules = [KomaMRI],
    sitename = "KomaMRI.jl: General MRI simulation framework",
    authors = "Carlos Castillo Passi",
    pages = [
        "Home" => "index.md",
        "Examples" => "examples.md"
    ]
)

deploydocs(
    repo = "github.com/cncastillo/KomaMRI.jl.git",
    # versions = ["stable" => "v^", "v#", "dev" => "dev"],
)

#julia --color=yes make.jl
