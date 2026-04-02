using Documenter, DocumenterVitepress, Literate, KomaMRI, PlutoSliderServer

# Setup for Literate and Pluto
repo_base = "JuliaHealth/KomaMRI.jl"
repo_root_url = "https://github.com/$repo_base/blob/master"
lit_pattern = "lit-"
plu_pattern = "pluto-"
include("utils.jl")
include("EmbeddPlotlyJSSyncPlotLiterate.jl")

# Documentation folders KomaMRI.jl/docs/
doc_tutorial       = joinpath(dirname(@__DIR__), "docs/src/tutorial")
doc_tutorial_rep   = joinpath(dirname(@__DIR__), "docs/src/tutorial/pluto")
doc_howto          = joinpath(dirname(@__DIR__), "docs/src/how-to")
doc_explanation    = joinpath(dirname(@__DIR__), "docs/src/explanation")
doc_reference      = joinpath(dirname(@__DIR__), "docs/src/reference")
# For Tutorials: Literate and Pluto
koma_assets        = joinpath(dirname(@__DIR__), "assets")
doc_assets         = joinpath(dirname(@__DIR__), "docs/src/public/assets")
doc_assets_tmp     = joinpath(dirname(@__DIR__), "docs/src/assets")
koma_tutorials_lit = joinpath(dirname(@__DIR__), "examples/3.tutorials")
koma_tutorials_plu = joinpath(dirname(@__DIR__), "examples/4.reproducible_notebooks")

# Copying files from KomaMRI.jl/ to the documentation folder KomaMRI.jl/docs/
# Assets
mkpath(doc_assets_tmp)
mkpath(doc_assets)
cp(joinpath(koma_assets, "logo.svg"), joinpath(doc_assets, "logo.svg"); force=true)
cp(joinpath(koma_assets, "logo-dark.svg"), joinpath(doc_assets, "logo-dark.svg"); force=true)
# Tutorials: Literate and Pluto
move_examples_to_docs!(koma_tutorials_lit, doc_tutorial, lit_pattern)
move_examples_to_docs!(koma_tutorials_plu, doc_tutorial_rep, plu_pattern; remove_pattern=true)

## DOCUMENTATION GENERATION
# Get list of documentation md files from docs/src/section
howto_list       = list_md_not_lit(doc_howto, "how-to"; lit_pattern)
explanation_list = list_md_not_lit(doc_explanation, "explanation"; lit_pattern)
reference_list   = list_md_not_lit(doc_reference, "reference"; lit_pattern)
# Add literate examples strarting with "lit-" from docs/src/section
lit_howto_list       = literate_doc_folder(doc_howto, "how-to")
lit_explanation_list = literate_doc_folder(doc_explanation, "explanation")
 lit_reference_list   = literate_doc_folder(doc_reference, "reference")
# Tutorials (Literate only), and reproducible tutorials (Pluto only)
tutorial_list     = literate_doc_folder(doc_tutorial, "tutorial"; lit_pattern)
reproducible_list = pluto_directory_to_html(doc_tutorial_rep, "tutorial/pluto"; plu_pattern)

# Combine md files in docs/src/section with Literate/Pluto-generated md files
append!(howto_list, lit_howto_list)
append!(explanation_list, lit_explanation_list)
append!(reference_list, lit_reference_list)

# Generate documentation
makedocs(;
    modules=[KomaMRI, KomaMRIBase, KomaMRICore, KomaMRIFiles, KomaMRIPlots],
    sitename="KomaMRI.jl",
    authors="Carlos Castillo Passi and collaborators",
    checkdocs=:exports,
    pages=[
        "Getting Started" => "introduction/1-getting-started.md",
        "Tutorials" => [
            "" => sort(tutorial_list),
            "Reproducible Tutorials" => sort(reproducible_list),
        ],
        "How to" => sort(howto_list),
        "Explanations" => sort(explanation_list),
        "Reference" => sort(reference_list),
    ],
    format=DocumenterVitepress.MarkdownVitepress(;
        repo = "github.com/JuliaHealth/KomaMRI.jl",
        devbranch = "master", 
        devurl = "dev",
        assets=["public/assets/hide-documenter-example-output.css","public/assets/center-images.css"],
    ),
    clean=false,
)

DocumenterVitepress.deploydocs(;
    repo = "github.com/JuliaHealth/KomaMRI.jl",
    target = joinpath(@__DIR__, "build"),
    branch = "gh-pages",
    devbranch = "master",
    push_preview = true,
)
