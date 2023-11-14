using Documenter, Literate, KomaMRI

org, reps = :cncastillo, :KomaMRI

base = "$org/$reps.jl"
repo_root_url = "https://github.com/$base/blob/master"

# User can add some directories from the literate folder
foldernames = ["examples"]

# Define some paths
com = joinpath(dirname(@__DIR__), "assets")
exa = joinpath(dirname(@__DIR__), "examples")
src = joinpath(@__DIR__, "src")
lit = joinpath(exa, "literate")
gen = joinpath(src, "generated")
assets = joinpath(src, "assets")
folderpaths = [joinpath(lit, foldername) for foldername in foldernames]

# Auxiliar function to replace texts into literate .jl files previous processing
function preprocess_constants(file, folder)
    return content -> content = replace(content, "FILE_NAME" => file, "FOLDER_NAME" => folder)
end

# Auxiliar function to check if a file has .md extension
function is_md_file(filename)
    return splitext(filename)[2] == ".md"
end

# Auxiliar function that returns an array with the subseccions of a literate seccion
function pages_for_literate_seccion(folder)
    return [joinpath("generated", folder, f) for f in readdir(joinpath(gen, folder)) if is_md_file(f)]
end

# Copy the logos into the assets folder if necessary
logo, logo_dark = joinpath(assets, "logo.svg"), joinpath(assets, "logo-dark.svg")
(!isfile(logo)) && cp(joinpath(com, "logo.svg"), logo)
(!isfile(logo_dark)) && cp(joinpath(com, "logo-dark.svg"), logo_dark)

# Create empty folders for in the assets directory for the literate generated sections
[mkpath(joinpath(assets, foldername)) for foldername in foldernames]
[mkpath(joinpath(gen, foldername)) for foldername in foldernames]

# Generate markdown, script and notebook for from the source literate file
for i ∈ eachindex(foldernames)
    subsection = joinpath(lit, foldernames[i])
    for (root, _, files) ∈ walkdir(subsection), file ∈ files
        file_minus_ext = split(file, ".")[1]
        ipath, opath = joinpath(root, file), joinpath(gen, foldernames[i])
        Literate.markdown(ipath, opath; repo_root_url, preprocess=preprocess_constants(file_minus_ext, foldernames[i]) )
        Literate.script(ipath, opath; repo_root_url)
        Literate.notebook(ipath, opath; execute=false)
    end
end

literate_seccions = [titlecase(replace(foldername, "-"=>" ")) => pages_for_literate_seccion(foldername) for foldername in foldernames]

# Documentation structure
makedocs(
    modules = [KomaMRI, KomaMRICore, KomaMRIPlots],
    sitename = "KomaMRI.jl: General MRI simulation framework",
    authors = "Boris Orostica Navarrete and Carlos Castillo Passi",
    pages = [
        "Home" => "index.md";
        "Getting Started" => "getting-started.md";
        "Simulation with User Interface" => "ui-details.md";
        "Simulation with Scripts" => "programming-workflow.md";
        "Create Your Own Phantom" => "create-your-own-phantom.md";
        "Create Your Own Sequence" => "create-your-own-sequence.md";
        "Sequence Definition" => "sequence.md";
        "Events Definition" => "events.md";
        literate_seccions;
        "Simulation Method" => "mri-theory.md";
        "API Documentation" => "api.md";
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
