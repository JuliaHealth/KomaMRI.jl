function loaded_data_toast(w, name, view, next, noun)
    short_name = html_escape(display_filename(name, MAX_TOAST_FILENAME_CHARS))
    full_name = html_escape(name)
    body = """
    <ul class="list-unstyled mb-0">
        <li>
            <button type="button" class="btn btn-dark btn-circle btn-circle-sm m-1" title="View $noun" aria-label="View $noun" onclick="KomaUI.notify('$view')"><i class="bi bi-search"></i></button>
            Updating <b>$noun</b> plots ...
        </li>
        <li>
            <button type="button" class="btn btn-primary btn-circle btn-circle-sm m-1" title="$next" aria-label="$next" onclick="KomaUI.notify('$next')"><i class="bi bi-caret-right-fill"></i></button>
            Ready to <b>$next</b>?
        </li>
    </ul>
    """
    return toast!(w, 0, "Loaded <b title=\"$full_name\">$short_name</b>", body)
end

function callback_filepicker(filename::String, w::KomaWindow, seq::Sequence)
    ext = splitext(filename)[2]
    if ext == ".seqk"
        seq = JLD2.load(FileIO.File{FileIO.DataFormat{:JLD2}}(filename), "seq")
    elseif ext == ".seq"
        seq = read_seq(filename)
    end
    name = basename(filename)
    update_filename!(w, "seqname", name)
    loaded_data_toast(w, name, "pulses_seq", "simulate", "Sequence")
    return seq
end

function callback_filepicker(filename::String, w::KomaWindow, obj::Phantom)
    ext = splitext(filename)[2]
    if ext == ".phantom"
        obj = read_phantom(filename)
    elseif ext == ".h5"
        obj = read_phantom_jemris(filename)
    end
    name = basename(filename)
    update_filename!(w, "phaname", name)
    loaded_data_toast(w, name, "phantom", "simulate", "Phantom")
    return obj
end

function callback_filepicker(filename::String, w::KomaWindow, raw::RawAcquisitionData)
    raw = RawAcquisitionData(ISMRMRDFile(filename))
    if raw.params["systemVendor"] != "KomaMRI.jl"
        @warn "ISMRMRD files generated externally could cause problems during the reconstruction. We are currently improving compatibility."
    end
    name = basename(filename)
    update_filename!(w, "rawname", name)
    loaded_data_toast(w, name, "sig", "recon", "Raw data")
    return raw
end

function setup_filepickers!(w::KomaWindow; seq_file=Ref(""))
    setup_filepicker!(
        w,
        "#seqfilepicker",
        "#seqname",
        ".seq (Pulseq)",
        seq_ui;
        accept=".seq,.seqk",
        selected_file=seq_file,
    )
    setup_filepicker!(
        w,
        "#phafilepicker",
        "#phaname",
        ".phantom (Koma)/.h5 (JEMRIS)",
        obj_ui;
        accept=".phantom,.h5",
    )
    setup_filepicker!(
        w,
        "#sigfilepicker",
        "#rawname",
        ".h5/.mrd (ISMRMRD)",
        raw_ui;
        accept=".h5,.mrd",
    )
    return nothing
end

upload_bytes(data) = UInt8.(data)
upload_bytes(data::MsgPack.Extension) = data.data

function setup_filepicker!(
    w::KomaWindow,
    selector::String,
    current::String,
    label::String,
    output;
    accept,
    selected_file=nothing,
)
    upload = Observable{Any}(nothing)
    push!(w.listeners, on(upload) do file
        isnothing(file) && return nothing
        name = String(get(file, "name", get(file, :name, "upload")))
        data = get(file, "data", get(file, :data, UInt8[]))
        bytes = upload_bytes(data)
        directory = mktempdir()
        filename = joinpath(directory, basename(name))
        open(filename, "w") do io
            write(io, bytes)
        end
        isnothing(selected_file) || (selected_file[] = filename)
        output[] = callback_filepicker(filename, w, output[])
        return nothing
    end)
    push!(w.on_render, session -> Bonito.on_document_load(session, js"""
            const mount = () => {
                const host = document.querySelector($(selector));
                if (!host) {
                    requestAnimationFrame(mount);
                    return;
                }
                const input = document.createElement('input');
                input.type = 'file';
                input.accept = $(accept);
                input.style.display = 'none';
                const button = document.createElement('label');
                button.className = 'koma-file-input';
                button.title = $(label);
                const icon = document.createElement('i');
                icon.className = 'bi bi-folder2-open';
                const caption = document.createElement('span');
                caption.className = 'koma-file-name';
                const setCaption = name => {
                    caption.textContent = name.length > $(MAX_UI_FILENAME_CHARS)
                        ? name.slice(0, $(MAX_UI_FILENAME_CHARS - 3)) + '...'
                        : name;
                    caption.title = name;
                    button.title = name;
                };
                setCaption(document.querySelector($(current))?.innerText.trim() || $(label));
                button.append(icon, caption, input);
                host.replaceChildren(button);
                input.addEventListener('change', async event => {
                    const file = event.target.files[0];
                    if (!file) return;
                    setCaption(file.name);
                    const data = new Uint8Array(await file.arrayBuffer());
                    $(upload).notify({name: file.name, data: data});
                });
            };
            mount();
        """))
    return nothing
end
