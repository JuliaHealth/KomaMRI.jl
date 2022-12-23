using Documenter, Literate, KomaMRI

exa = joinpath(@__DIR__, "../examples")
src = joinpath(@__DIR__, "src")
gen = joinpath(@__DIR__, "src/generated")

for (root, _, files) ∈ walkdir(exa), file ∈ files
    #splitext(file)[2] == ".jl" || continue
    cmp(last(splitpath(root)), "examples") == 0 && cmp(first(split(file, "-")), "lit") == 0 || continue
    ipath = joinpath(root, file)
    opath = splitdir(replace(ipath, exa=>gen))[1]
    Literate.markdown(ipath, opath)
    #Literate.notebook(ipath, opath; execute = false)
end

# Documentation structure
ismd(f) = splitext(f)[2] == ".md"
pages() = [joinpath("generated", f) for f in readdir(gen) if ismd(f)]
#pages(folder) = [joinpath("generated", folder, f) for f in readdir(joinpath(gen, folder)) if ismd(f)]

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
        "Examples" => pages(),
        "Simulation Method" => "mri-theory.md",
        "API Documentation" => "api.md"
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
