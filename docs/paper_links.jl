function write_paper_redirects!(build_dir, papers)
    stable_url = "https://juliahealth.org/KomaMRI.jl/stable/"
    marker = "Generated from docs/paper-links.toml"

    for links in values(papers), (source, target) in links
        startswith(source, stable_url) || continue
        startswith(target, stable_url) ||
            error("Paper link target must remain under $stable_url: $target")

        source_path = strip(chopprefix(source, stable_url), '/')
        target_parts = split(chopprefix(target, stable_url), '#'; limit=2)
        target_path = strip(first(target_parts), '/')
        fragment = length(target_parts) == 2 ? "#$(last(target_parts))" : ""
        relative_target =
            joinpath(relpath(dirname(target_path), source_path), basename(target_path))
        relative_target = replace(relative_target, '\\' => '/') * fragment

        redirect_file = joinpath(build_dir, split(source_path, '/')..., "index.html")
        if isfile(redirect_file) && !occursin(marker, read(redirect_file, String))
            error("Paper link conflicts with generated documentation: $source")
        end
        mkpath(dirname(redirect_file))
        write(
            redirect_file,
            """<!doctype html>
            <!-- $marker -->
            <meta charset="utf-8">
            <meta http-equiv="refresh" content="0; url=$relative_target">
            <link rel="canonical" href="$target">
            <title>Page moved</title>
            <p>Page moved to <a href="$relative_target">$target</a>.</p>
            """,
        )
    end
    return nothing
end
