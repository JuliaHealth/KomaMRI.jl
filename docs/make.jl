using Documenter, Literate, KomaMRI

org, reps = :cncastillo, :KomaMRI

base = "$org/$reps.jl"
repo_root_url = "https://github.com/$base/blob/master"

# Define some paths
com = joinpath(dirname(@__DIR__), "assets")
exa = joinpath(dirname(@__DIR__), "examples")
src = joinpath(@__DIR__, "src")
edu = joinpath(exa, "6.educational_pluto_notebook")
lit = joinpath(exa, "literate")
gen = joinpath(src, "generated")
assets = joinpath(src, "assets")
exaname = "examples"
eduname = "educational-examples"

# Create empty folders for in the assets directory for the literate generated sections
mkpath(joinpath(assets, exaname))
mkpath(joinpath(gen, exaname))
mkpath(joinpath(gen, eduname))

# Copy the logos into the assets folder if necessary
logo, logo_dark = joinpath(assets, "logo.svg"), joinpath(assets, "logo-dark.svg")
(!isfile(logo)) && cp(joinpath(com, "logo.svg"), logo)
(!isfile(logo_dark)) && cp(joinpath(com, "logo-dark.svg"), logo_dark)

# Auxiliar function to check if a file has .html extension
function is_html_file(filename)
    return splitext(filename)[2] == ".html"
end

# Copy html files for Pluto educational example
for (root, _, files) ∈ walkdir(edu), file ∈ files
    if is_html_file(file)
        ipath, opath = joinpath(root, file), joinpath(gen, eduname, file)
        cp(ipath, opath; force=true)
    end
end

# Auxiliar function to replace texts into literate .jl files previous processing
function preprocess_constants(file, folder)
    return content -> content = replace(content, "FILE_NAME" => file, "FOLDER_NAME" => folder)
end

# Auxiliar function to check if a file has .md extension
function is_md_file(filename)
    return splitext(filename)[2] == ".md"
end

# Generate markdown, script and notebook for from the source literate file
for (root, _, files) ∈ walkdir(joinpath(lit, exaname)), file ∈ files
    file_minus_ext = split(file, ".")[1]
    ipath, opath = joinpath(root, file), joinpath(gen, exaname)
    Literate.markdown(ipath, opath; repo_root_url, preprocess=preprocess_constants(file_minus_ext, exaname) )
    Literate.script(ipath, opath; repo_root_url)
    Literate.notebook(ipath, opath; execute=false)
end

# Create some subsections
ways_of_using_koma = ["ui-details.md", "programming-workflow.md", "notebooks.md"]
create_your_own_sequence = ["sequence.md", "events.md"]
literate_examples = [joinpath("generated", exaname, f) for f in readdir(joinpath(gen, exaname)) if is_md_file(f)]

# Documentation structure
makedocs(
    modules = [KomaMRI, KomaMRIBase, KomaMRICore, KomaMRIFiles, KomaMRIPlots],
    sitename = "KomaMRI.jl: General MRI simulation framework",
    authors = "Boris Orostica Navarrete and Carlos Castillo Passi",
    checkdocs = :exports,
    pages = [
        "Home" => "index.md";
        "Getting Started" => "getting-started.md";
        "Ways of Using KomaMRI" => ways_of_using_koma;
        #"Simulation with User Interface" => "ui-details.md";
        #"Simulation with Scripts" => "programming-workflow.md";
        #"Notebooks" => "notebooks.md";
        "Create Your Own Phantom" => "create-your-own-phantom.md";
        "Create Your Own Sequence" => create_your_own_sequence;
        #"Sequence Definition" => "sequence.md";
        #"Events Definition" => "events.md";
        "Examples" => literate_examples;
        "Educational Material 📚" => "educational-1d-simulation.md";
        "Simulation" => "mri-theory.md";
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
