using Documenter, Koma

makedocs(sitename="Koma.jl")

makedocs(
    modules = [Koma],
    sitename = "Koma.jl: General MRI simulator",
    authors = "Carlos Castillo",
    pages = [
        "Home" => "index.md",
    ],
)

deploydocs(
    repo = "github.com/cncastillo/Koma.jl.git",
)

#makedocs(format = LaTeX())
#julia --color=yes make.jl
