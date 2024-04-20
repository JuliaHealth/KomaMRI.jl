using Documenter, Literate, KomaMRI, PlutoSliderServer

org, reps = :cncastillo, :KomaMRI
base = "$org/$reps.jl"
repo_root_url = "https://github.com/$base/blob/master"
include("utils.jl")

# Documentation folders
doc_assets         = joinpath(dirname(@__DIR__), "docs/src/assets")
doc_tutorial       = joinpath(dirname(@__DIR__), "docs/src/tutorial")
doc_tutorial_pluto = joinpath(dirname(@__DIR__), "docs/src/tutorial-pluto")
doc_howto          = joinpath(dirname(@__DIR__), "docs/src/how-to")
doc_explanation    = joinpath(dirname(@__DIR__), "docs/src/explanation")
doc_reference      = joinpath(dirname(@__DIR__), "docs/src/reference")
# External files
koma_assets = joinpath(dirname(@__DIR__), "assets")
koma_tutorials_lit = joinpath(dirname(@__DIR__), "examples/3.tutorials")
koma_tutorials_plu = joinpath(dirname(@__DIR__), "examples/4.reproducible_notebooks")
lit_start_pattern = "lit-"
pluto_start_pattern = "pluto-"

# Copying external files to the documentation folder
cp(joinpath(koma_assets, "logo.svg"), joinpath(doc_assets, "logo.svg"); force=true)
cp(joinpath(koma_assets, "logo-dark.svg"), joinpath(doc_assets, "logo-dark.svg"); force=true)
move_pattern_from_examples_to_docs!(koma_tutorials_lit, doc_tutorial, lit_start_pattern)
move_pattern_from_examples_to_docs!(koma_tutorials_plu, doc_tutorial_pluto, pluto_start_pattern)

## DOCUMENTATION GENERATION
# Get list of doc md files
howto_list = list_md_not_lit(doc_howto, "how-to"; exclude="getting-started", lit_start_pattern)
explanation_list = list_md_not_lit(doc_explanation, "explanation"; lit_start_pattern)
reference_list = list_md_not_lit(doc_reference, "reference"; lit_start_pattern)
# Add literate examples strarting with "lit-"
lit_howto_list = literate_doc_folder(doc_howto, "how-to")
lit_explanation_list = literate_doc_folder(doc_explanation, "explanation")
lit_reference_list = literate_doc_folder(doc_reference, "reference")
# Tutorials can only be literate, and reproducible tutorials are pluto
tutorial_list = literate_doc_folder(doc_tutorial, "tutorial"; lit_start_pattern)
reproducible_list =  pluto_directory_to_html_pattern(doc_tutorial_pluto, "tutorial-pluto"; pluto_start_pattern)
# Combine md files in docs with lit and pluto (copied from koma examples)
append!(howto_list, lit_howto_list)
append!(explanation_list, lit_explanation_list)
append!(reference_list, lit_reference_list)
# Sorting files to order documentation sections
sort!(tutorial_list)
sort!(reproducible_list)
sort!(howto_list)
sort!(explanation_list)
sort!(reference_list)

# Generate documentation
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
