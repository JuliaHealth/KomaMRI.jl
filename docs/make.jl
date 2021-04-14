using Documenter, MRIsim

makedocs(sitename="MRIsim.jl")

deploydocs(
    repo = "github.com/cncastillo/MRIsim.jl.git",
)

#makedocs(format = LaTeX())
#julia --color=yes make.jl
