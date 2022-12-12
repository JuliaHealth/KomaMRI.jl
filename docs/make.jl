using Documenter, Literate, KomaMRI

lit = joinpath(@__DIR__, "lit")
src = joinpath(@__DIR__, "src")

for (root, _, files) ∈ walkdir(lit), file ∈ files
    splitext(file)[2] == ".jl" || continue
    ipath = joinpath(root, file)
    opath = splitdir(replace(ipath, lit=>src))[1]
    Literate.markdown(ipath, opath)
    #Literate.notebook(ipath, opath; execute = false)
end

makedocs(
    modules = [KomaMRI],
    sitename = "KomaMRI.jl: General MRI simulation framework",
    authors = "Boris Orostica Navarrete and Carlos Castillo Passi",
    pages = [
        "Home" => "index.md",
        "Getting Started" => "getting-started.md",
        # "Sequence" => "sequence.md",
        # "Phantom" => "phantom.md",
        # "Scanner" => "scanner.md",
        # "Raw Signal" => "raw-signal.md",
        "Graphical User Interface" => "ui-details.md",
        "Examples" => "simulation-examples.md",
        "Simulation Method" => "mri-theory.md",
        "API Documentation" => "api.md",
        "Literate Examples" => "examples.md"
    ],
    format = Documenter.HTML(
        prettyurls = true, #get(ENV, "CI", nothing) == "true",
        sidebar_sitename = false,
        assets = ["assets/extra-styles.css"]
    )
)

deploydocs(
    repo = "github.com/cncastillo/KomaMRI.jl.git",
)
