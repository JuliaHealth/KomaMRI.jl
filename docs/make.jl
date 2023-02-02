using Documenter, Literate, KomaMRI

org, reps = :cncastillo, :KomaMRI

base = "$org/$reps.jl"
repo_root_url = "https://github.com/$base/blob/master"

exa = joinpath(@__DIR__, "../examples")
src = joinpath(@__DIR__, "src")
gen = joinpath(@__DIR__, "src/generated")

function update_file_name(file)
    return content -> content = replace(content, "FILE_NAME" => file)
end

for (root, _, files) ∈ walkdir(exa), file ∈ files
    #splitext(file)[2] == ".jl" || continue
    cmp(last(splitpath(root)), "examples") == 0 && cmp(first(split(file, "-")), "lit") == 0 || continue
    file_minus_ext = split(file, ".")[1] #lit-02-example.jl => lit-02-example
    ipath = joinpath(root, file)
    Literate.markdown(ipath, gen; repo_root_url, preprocess=update_file_name(file_minus_ext) )
    Literate.script(ipath, gen; repo_root_url)
    Literate.notebook(ipath, gen; execute=false)
end

# Documentation structure
ismd(f) = splitext(f)[2] == ".md"
pages() = [joinpath("generated", f) for f in readdir(gen) if ismd(f)]
#pages(folder) = [joinpath("generated", folder, f) for f in readdir(joinpath(gen, folder)) if ismd(f)]

makedocs(
    modules = [KomaMRI, KomaMRICore, KomaMRIPlots],
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
