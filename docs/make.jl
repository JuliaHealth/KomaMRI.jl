using Documenter, Literate, KomaMRI, PlutoSliderServer

org, reps = :JuliaHealth, :KomaMRI
base = "$org/$reps.jl"
repo_root_url = "https://github.com/$base/blob/master"

# Documentation folders
doc_assets         = joinpath(@__DIR__, "src/assets")
doc_tutorial       = joinpath(@__DIR__, "src/tutorial")
doc_tutorial_pluto = joinpath(@__DIR__, "src/tutorial-pluto")
doc_howto          = joinpath(@__DIR__, "src/how-to")
doc_explanation    = joinpath(@__DIR__, "src/explanation")
doc_reference      = joinpath(@__DIR__, "src/reference")
# External files
koma_assets = joinpath(dirname(@__DIR__), "assets")
koma_tutorials_lit = joinpath(dirname(@__DIR__), "examples/3.tutorials")
koma_tutorials_plu = joinpath(dirname(@__DIR__), "examples/4.reproducible_notebooks")

# Copying external files to the documentation folder
# docs/assets
cp(joinpath(koma_assets, "logo.svg"), joinpath(doc_assets, "logo.svg"); force=true)
cp(
    joinpath(koma_assets, "logo-dark.svg"),
    joinpath(doc_assets, "logo-dark.svg");
    force=true,
)
# docs/src
# Literate
lit_start_pattern = "lit-"
for (root, _, files) in walkdir(koma_tutorials_lit)
    for filename in filter(startswith(lit_start_pattern), files)
        if !endswith(filename, ".jl")
            continue
        end
        cp(joinpath(root, filename), joinpath(doc_tutorial, filename); force=true)
    end
end
# Pluto
pluto_start_pattern = "pluto-"
for (root, _, files) in walkdir(koma_tutorials_plu)
    for filename in filter(startswith(pluto_start_pattern), files)
        if !endswith(filename, ".jl")
            continue
        end
        filename_gen = filename[(length(pluto_start_pattern) + 1):end] # removes "pluto-"
        cp(joinpath(root, filename), joinpath(doc_tutorial_pluto, filename_gen); force=true)
    end
end

# Based on https://github.com/jump-dev/JuMP.jl/blob/master/docs/make.jl
# Add Literate links after the title
function _link_example(filename)
    function _link_example_for_filename(content)
        title_line = findfirst(r"\n# .+?\n", content)
        line = content[title_line]
        badges = """

        #md # [![](https://img.shields.io/badge/julia-script-9558B2?logo=julia)](./$filename.jl)
        #md # [![](https://img.shields.io/badge/jupyter-notebook-blue?logo=jupyter)](./$filename.ipynb)

        """
        return replace(content, line => line * badges)
    end
    return _link_example_for_filename
end

# TUTORIALS: Generate markdown, script and notebook from the source literate files
tutorial_list = []
for (root, _, files) in walkdir(doc_tutorial)
    for filename in filter(startswith(lit_start_pattern), files)
        filename_gen = splitext(filename)[1][(length(lit_start_pattern) + 1):end] # removes "lit-"
        tutorial_src = joinpath(doc_tutorial, filename)
        # tutorial_out = joinpath(doc_tutorial, filename_gen)
        tutorial_md  = "tutorial/$filename_gen.md"
        Literate.markdown(
            tutorial_src,
            doc_tutorial;
            repo_root_url,
            preprocess=_link_example(filename_gen),
            name=filename_gen,
        )
        Literate.script(tutorial_src, doc_tutorial; name=filename_gen, repo_root_url)
        Literate.notebook(tutorial_src, doc_tutorial; name=filename_gen, execute=false)
        push!(tutorial_list, tutorial_md)
    end
end
sort!(tutorial_list)

# REPRODUCIBLE TUTORIALS
reproducible_list = String[]
for (root, _, files) in walkdir(doc_tutorial_pluto)
    for filename in filter(endswith("jl"), files)
        filename_gen  = splitext(filename)[1]
        tutorial_src  = joinpath(doc_tutorial_pluto, filename)
        tutorial_md   = joinpath(doc_tutorial_pluto, "$filename_gen.md")
        # HTML to Markdown
        frontmatter = PlutoSliderServer.Pluto.frontmatter(tutorial_src)
        iframe = """
        # $(frontmatter["title"])

        ```@raw html
        <iframe type="text/html" src="../$filename_gen.html" style="height:100vh;width:100%;"></iframe>
        ```
        """
        open(tutorial_md, "w") do file
            write(file, iframe)
        end
        push!(reproducible_list, "tutorial-pluto/$filename_gen.md")
    end
end
sort!(reproducible_list)
# Pluto to HTML
PlutoSliderServer.export_directory(doc_tutorial_pluto)

# HOW-TO GUIDES
howto_list = String[]
for (root, _, files) in walkdir(doc_howto)
    for filename in filter(endswith(".md"), files)
        if occursin("getting-started", filename)
            continue
        end
        push!(howto_list, "how-to/$filename")
    end
end
sort!(howto_list)

# EXPLANATIONS
explanation_list = String[]
for (root, _, files) in walkdir(doc_explanation)
    for filename in filter(endswith(".md"), files)
        push!(explanation_list, "explanation/$filename")
    end
end
sort!(explanation_list)

# REFERENCE GUIDES
reference_list = String[]
for (root, _, files) in walkdir(doc_reference)
    for filename in filter(endswith(".md"), files)
        push!(reference_list, "reference/$filename")
    end
end
sort!(reference_list)

# Documentation structure
makedocs(;
    modules=[KomaMRI, KomaMRIBase, KomaMRICore, KomaMRIFiles, KomaMRIPlots],
    sitename="KomaMRI.jl",
    authors="Boris Orostica Navarrete and Carlos Castillo Passi",
    checkdocs=:exports,
    pages=[
        "🏠 Home" => "index.md",
        "🏃 Getting Started" => "how-to/1-getting-started.md",
        "🏋️ Tutorials" => tutorial_list,
        "🧑‍🔬 Reproducible Tutorials" => reproducible_list,
        "👨‍🍳 How to" => howto_list,
        "🤔 Explanations" => explanation_list,
        "👨‍💻 Reference Guides" => reference_list,
    ],
    format=Documenter.HTML(;
        prettyurls=true, #get(ENV, "CI", nothing) == "true",
        sidebar_sitename=false,
        collapselevel=1,
        assets=["assets/extra-styles.css"],
    ),
)

deploydocs(; repo="github.com/JuliaHealth/KomaMRI.jl.git", push_preview=true)
