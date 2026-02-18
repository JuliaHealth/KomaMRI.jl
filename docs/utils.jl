function move_examples_to_docs!(src_folder, dst_folder, start_pattern; remove_pattern=false)
    for (_, _, files) in walkdir(src_folder)
        for filename in filter(startswith(start_pattern), files)
            if !endswith(filename, ".jl")
                continue
            end
            # removes "pluto-"/"lit-" from the filename
            if remove_pattern
                filename_gen = filename[(length(start_pattern) + 1):end]
            else
                filename_gen = filename
            end
            cp(joinpath(src_folder, filename), joinpath(dst_folder, filename_gen); force=true)
        end
    end
end

# Based on https://github.com/jump-dev/JuMP.jl/blob/master/docs/make.jl
# Add Literate links after the title
function _link_example(filename)
    function _link_example_for_filename(content)
        title_line = findfirst(r"\n# .+?\n", content)
        line = content[title_line]
        koma_version = "dev"
        binder_link = "https://mybinder.org/v2/gh/$repo_base/master?urlpath=git-pull"
        binder_gitpull = "?repo=https://github.com/$repo_base&urlpath=lab/tree/KomaMRI.jl/$koma_version/tutorial/$filename.ipynb&branch=gh-pages"
        binder_gitpull = replace(binder_gitpull, "?"=>"%3F", "="=>"%3D", ":"=>"%253A", "/"=>"%252F", "&"=>"%26")
        badges = """

        #md # [![](https://img.shields.io/badge/julia-script-9558B2?logo=julia)](./$filename.jl) [![](https://img.shields.io/badge/jupyter-notebook-blue?logo=jupyter)](./$filename.ipynb) [![](https://mybinder.org/badge_logo.svg)]($(binder_link)$(binder_gitpull))

        """
        return replace(content, line => badges * line)
    end
    return _link_example_for_filename
end

function list_md_not_lit(input_folder, output_doc_section; exclude="-----------------------", lit_pattern="lit-")
    md_list = String[]
    for (_, _, files) in walkdir(input_folder)
        for filename in filter(endswith(".md"), files)
            if startswith(lit_pattern, filename)
                continue
            end
            if occursin(exclude, filename)
                continue
            end
            push!(md_list, "$output_doc_section/$filename")
        end
    end
    return md_list
end

function literate_doc_folder(input_folder, output_doc_section; lit_pattern="lit-")
    tutorial_list = []
    for (_, _, files) in walkdir(input_folder)
        for filename in filter(startswith(lit_pattern), files)
            filename_gen = splitext(filename)[1][(length(lit_pattern) + 1):end] # removes "lit-"
            tutorial_src = joinpath(input_folder, filename)
            tutorial_md  = "$output_doc_section/$filename_gen.md"
            Literate.markdown(
                tutorial_src,
                input_folder;
                repo_root_url,
                preprocess=_link_example(filename_gen),
                name=filename_gen,
                execute=true
            )
            Literate.script(tutorial_src, input_folder; name=filename_gen, repo_root_url)
            Literate.notebook(tutorial_src, input_folder; name=filename_gen, execute=false)
            push!(tutorial_list, tutorial_md)
        end
    end
    return tutorial_list
end

# TODO: copy files with "pluto-" to docs, and remove for generated html and md
function pluto_directory_to_html(doc_tutorial_pluto, doc_output_section; plu_pattern="pluto-")
    reproducible_list = String[]
    for (_, _, files) in walkdir(doc_tutorial_pluto)
        for filename in filter(endswith("jl"), files)
            # if !startswith(plu_pattern, filename)
            #     continue
            # end
            filename_gen  = splitext(filename)[1]
            tutorial_src  = joinpath(doc_tutorial_pluto, filename)
            tutorial_md   = joinpath(doc_tutorial_pluto, "$filename_gen.md")
            # HTML to Markdown
            frontmatter = PlutoSliderServer.Pluto.frontmatter(tutorial_src)
            koma_version = "dev"
            binder_link = "https://mybinder.org/v2/gh/$repo_base/master?urlpath=git-pull"
            binder_gitpull = "?repo=https://github.com/$repo_base&urlpath=pluto/open?path=KomaMRI.jl/$koma_version/tutorial-pluto/$filename&branch=gh-pages"
            binder_gitpull = replace(binder_gitpull, "?"=>"%3F", "="=>"%3D", ":"=>"%253A", "/"=>"%252F", "&"=>"%26")

            iframe = """
            # $(frontmatter["title"])

            [![](https://img.shields.io/badge/julia-script-9558B2?logo=julia)](./$filename) [![](https://mybinder.org/badge_logo.svg)]($(binder_link)$(binder_gitpull))
            
            ```@raw html
            <iframe type="text/html" src="../$filename_gen.html" style="height:100vh;width:100%;"></iframe>
            ```
            """
            open(tutorial_md, "w") do file
                write(file, iframe)
            end
            push!(reproducible_list, "$doc_output_section/$filename_gen.md")
        end
    end
    PlutoSliderServer.export_directory(doc_tutorial_pluto)
    return reproducible_list
end