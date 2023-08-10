using Documenter, Literate, KomaMRI

org, reps = :cncastillo, :KomaMRI

base = "$org/$reps.jl"
repo_root_url = "https://github.com/$base/blob/master"

exa = joinpath(@__DIR__, "../examples")
src = joinpath(@__DIR__, "src")
gen = joinpath(@__DIR__, "src/generated")
foldernames = ["basic"; "sequence-design"]
folderpaths = [joinpath(exa, foldername) for foldername in foldernames]

function preprocess_constants(file, folder)
    return content -> content = replace(content, "FILE_NAME" => file, "FOLDER_NAME" => folder)
end

for i ∈ eachindex(foldernames)
    for (root, _, files) ∈ walkdir(folderpaths[i]), file ∈ files
        file_minus_ext = split(file, ".")[1]
        ipath, opath = joinpath(root, file), joinpath(gen, foldernames[i])
        Literate.markdown(ipath, opath; repo_root_url, preprocess=preprocess_constants(file_minus_ext, foldernames[i]) )
        Literate.script(ipath, opath; repo_root_url)
        Literate.notebook(ipath, opath; execute=false)
    end
end

# Documentation structure
ismd(f) = splitext(f)[2] == ".md"
pages(folder) = [joinpath("generated", folder, f) for f in readdir(joinpath(gen, folder)) if ismd(f)]

makedocs(
    modules = [KomaMRI, KomaMRICore, KomaMRIPlots],
    sitename = "KomaMRI.jl: General MRI simulation framework",
    authors = "Boris Orostica Navarrete and Carlos Castillo Passi",
    pages = [
        "Home" => "index.md",
        #"Getting Started" => "getting-started.md",
        #"Graphical User Interface" => "ui-details.md",
        "Examples" => pages(foldernames[1]),
        "Sequence Design" => pages(foldernames[2]),
        #"Simulation Method" => "mri-theory.md",
        #"API Documentation" => "api.md"
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
