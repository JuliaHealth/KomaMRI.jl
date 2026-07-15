function select_export_folder(w::KomaWindow)
    isnothing(w.display[]) && return tempdir()
    folders = Base.run(w.display[].window, """
        electron.dialog.showOpenDialogSync(
            electron.BrowserWindow.getFocusedWindow(),
            {
                title: "Export .mat files",
                buttonLabel: "Export here",
                properties: ["openDirectory", "createDirectory"]
            }
        )
    """)
    return isnothing(folders) || isempty(folders) ? nothing : String(first(folders))
end

function save_ui!(
    w::KomaWindow, seq::Sequence, obj::Phantom, sys::Scanner, raw, img, rec_params; type="all"
)
    mat_folder = select_export_folder(w)
    isnothing(mat_folder) && return nothing
    message = if type == "all"
        export_2_mat(seq, obj, sys, raw, rec_params, img, mat_folder; type, matfilename="")
    else
        export_2_mat(seq, obj, sys, raw, rec_params, img, mat_folder; type)
    end
    toast!(w, 1, "Saved .mat files", message)
    suffix = Dict(
        "all" => "",
        "sequence" => "seq",
        "phantom" => "pha",
        "scanner" => "sca",
        "raw" => "raw",
        "image" => "ima",
    )[type]
    state = "matfolder$suffix"
    w.state[] = state
    evaljs(w, js"document.getElementById('content').dataset.content = $(state);")
    return nothing
end
