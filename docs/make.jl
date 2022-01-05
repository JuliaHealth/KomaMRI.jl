using Documenter, Koma

makedocs(sitename="Koma.jl")

deploydocs(
    repo = "github.com/cncastillo/MRIsim.jl.git",
)

#makedocs(format = LaTeX())
#julia --color=yes make.jl
