using Documenter, Koma

makedocs(sitename="Koma.jl")

deploydocs(
    repo = "github.com/cncastillo/Koma.jl.git",
)

#makedocs(format = LaTeX())
#julia --color=yes make.jl
