using Documenter, KomaMRI

makedocs(
    modules = [KomaMRI],
    sitename = "KomaMRI.jl: General MRI simulation framework",
    authors = "Boris Orostica Navarrete and Carlos Castillo Passi",
    pages = [
        "Home" => "index.md",
        "Getting Started" => "getting-started.md",
        "Sequence" => "sequence.md",
        "Simulation Examples" => "simulation-examples.md",
        "API Documentation" => "api.md"
    ],
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true",
        sidebar_sitename = false,
        assets = ["assets/extra-styles.css"]
    )
)

deploydocs(
    repo = "github.com/cncastillo/KomaMRI.jl.git",
    # versions = ["stable" => "v^", "v#", "dev" => "dev"],
)

#julia --color=yes make.jl
