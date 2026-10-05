"""
    export_results(w, output; path)

Save the current `:raw_data_mrd`, `:raw_data_mat`, or `:image_mat` output to an
explicit file path, without a picker, download, or view change. The output key
determines the format, not the filename extension. MAT exports contain one file
without parameter sidecars.

The path is on the machine running Julia; its parent directory must exist.
Existing files are overwritten. Returns `(; path, url)` with the absolute path
and a download URL in browser mode (`nothing` otherwise). The URL remains valid
while the UI server runs and the file exists. For a remote server, retrieve it
through the SSH tunnel, substituting its local forwarded port if different.
"""
function export_results(w, output; path)
    w.operation[].state === :running && error("KomaUI is busy with $(w.operation[].action).")
    path = abspath(path)
    if output === :raw_data_mrd
        save(ISMRMRDFile(path), w.raw[])
    elseif output === :raw_data_mat
        export_2_mat_raw(w.raw[], dirname(path); matfilename=basename(path), include_params=false)
    elseif output === :image_mat
        export_2_mat_image(w.img[], w.rec_params[], dirname(path); matfilename=basename(path), include_params=false)
    else
        throw(ArgumentError("Unknown output $output; use :raw_data_mrd, :raw_data_mat, or :image_mat."))
    end
    url = nothing
    if !isnothing(w.server[])
        if isnothing(w.export_assets[])
            w.export_assets[] = Bonito.HTTPAssetServer(w.server[])
        end
        url = Bonito.url(w.export_assets[], Asset(path))
    end
    return (; path, url)
end

function setup_exports!(w)
    push!(w.on_render, session -> Bonito.on_document_load(session, js"""(() => {
        const handles = new Map();
        const buttons = document.querySelectorAll('#export .dropdown-item');
        const reset = () => {
            handles.clear();
            buttons.forEach(button => { button.disabled = false; });
        };
        const toast = (id, title, text) => {
            const message = document.createElement('span');
            message.textContent = text;
            document.querySelector('#liveToast' + id + ' rect')
                .setAttribute('fill', id === 0 ? '#dc3545' : '#2a7fb8');
            Toasty(id, title, message.innerHTML);
        };
        window.KomaUIExport = {
            reset,
            async select(action, filename) {
                if (!$(!isnothing(w.server[])) || typeof window.showSaveFilePicker !== 'function') {
                    KomaUI.notify(action);
                    return;
                }
                buttons.forEach(button => { button.disabled = true; });
                try {
                    handles.set(action, await window.showSaveFilePicker({suggestedName: filename}));
                    KomaUI.notify(action);
                } catch (error) {
                    reset();
                    if (error.name !== 'AbortError') toast(0, 'Export failed', error.message);
                }
            },
            async save(action, files) {
                const handle = handles.get(action);
                if (!handle) {
                    for (const file of files) {
                        const link = document.createElement('a');
                        link.href = file.url;
                        link.download = file.name;
                        document.body.appendChild(link);
                        link.click();
                        link.remove();
                    }
                    toast(1, 'Download requested', files.map(file => file.name).join(', '));
                    return;
                }
                try {
                    const response = await fetch(files[0].url);
                    if (!response.ok) throw new Error('Export download failed: HTTP ' + response.status);
                    const writable = await handle.createWritable();
                    await response.body.pipeTo(writable);
                    toast(1, 'Saved results', handle.name);
                } catch (error) {
                    toast(0, 'Export failed', error.message);
                } finally {
                    reset();
                }
            }
        };
    })()"""))
    return nothing
end

function select_export_folder(w::KomaWindow)
    isnothing(w.display[]) && return isnothing(w.server[]) ? tempdir() : mktempdir()
    window = w.window[]
    folders = Base.run(window.app, """
        electron.dialog.showOpenDialogSync(
            electron.BrowserWindow.fromId($(window.id)),
            {
                title: "Export results",
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
    output = get(Dict(
        "raw" => (:raw_data_mat, "data_raw.mat"),
        "image" => (:image_mat, "data_image.mat"),
        "mrd" => (:raw_data_mrd, "data_raw.mrd"),
    ), type, nothing)
    result = nothing
    message = try
        if !isnothing(output)
            key, filename = output
            result = export_results(w, key; path=joinpath(mat_folder, filename))
            "<ul><li><b>Name:</b> $filename</li><li><b>Path:</b> $mat_folder</li></ul>"
        elseif type == "all"
            export_2_mat(seq, obj, sys, raw, rec_params, img, mat_folder; type, matfilename="")
        else
            export_2_mat(seq, obj, sys, raw, rec_params, img, mat_folder; type, include_params=false)
        end
    catch
        evaljs(w, js"window.KomaUIExport?.reset()")
        rethrow()
    end
    suffix = Dict(
        "all" => "",
        "sequence" => "seq",
        "phantom" => "pha",
        "scanner" => "sca",
        "raw" => "raw",
        "image" => "ima",
        "mrd" => "mrd",
    )[type]
    state = "matfolder$suffix"
    if isnothing(w.server[])
        toast!(w, 1, "Saved results", message)
    else
        downloads = isnothing(result) ?
            [(url=Bonito.url(w.session[], Asset(file)), name=basename(file)) for file in readdir(mat_folder; join=true)] :
            [(url=result.url, name=basename(result.path))]
        action = type == "raw" ? "export_raw_data" : "export_$type"
        evaljs(w, js"KomaUIExport.save($(action), $(downloads))")
    end
    w.state[] = state
    evaljs(w, js"document.getElementById('content').dataset.content = $(state);")
    return nothing
end
