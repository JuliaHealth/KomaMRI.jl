using Documenter, KomaMRI

makedocs(
    modules = [KomaMRI],
    sitename = "KomaMRI.jl: General MRI simulation framework",
    authors = "Carlos Castillo Passi",
    pages = [
        "Getting Started" => "index.md",
        "CLI Examples" => "cli-examples.md",
        "Software" => "software.md"
    ]
)

deploydocs(
    repo = "github.com/cncastillo/KomaMRI.jl.git",
    # versions = ["stable" => "v^", "v#", "dev" => "dev"],
)

#julia --color=yes make.jl
