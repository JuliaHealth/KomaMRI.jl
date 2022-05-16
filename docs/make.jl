using Documenter, Koma

makedocs(
    modules = [Koma],
    sitename = "Koma.jl :General MRI simulation framework",
    authors = "Carlos Castillo Passi",
    pages = [
        "Home" => "index.md",
        "API" => "API.md"
    ]
)

# makedocs(
#     modules = [Koma],
#     sitename = "Koma.jl :General MRI simulation framework",
#     authors = "Carlos Castillo Passi",
#     pages = [
#         "Home" => "index.md",
#     ],
#     format = LaTeX()
# )

deploydocs(
    repo = "github.com/cncastillo/Koma.jl.git",
    versions = ["stable" => "v^", "v#", "dev" => "dev"],
)

#julia --color=yes make.jl
