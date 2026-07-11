struct BonitoWindow <: KomaWindow
    app::App
    display::Base.RefValue{Any}
    session::Base.RefValue{Union{Nothing,Session}}
    content::Observable{Any}
    state::Observable{String}
    events::Observable{String}
    home::Base.RefValue{Any}
    handlers::Dict{String,Function}
    on_render::Vector{Function}
    window_options::Dict{String,Any}
    dev_tools::Bool
end

function page_content(content, state)
    return DOM.div(
        content;
        id="content",
        class="koma-content",
        dataContent=state,
    )
end

function set_content!(w::BonitoWindow, content, state)
    w.state[] = state
    w.content[] = page_content(content, state)
    evaljs(w, js"""(() => {
        const section = {
            sequence: 'pulses', kspace: 'pulses', m0: 'pulses', m1: 'pulses', m2: 'pulses',
            phantom: 'phantom', scanneparams: 'scanner', sig: 'sig',
            absi: 'recon', angi: 'recon', absk: 'recon'
        }[$(state)];
        if ($(state) === 'loading') return;
        document.querySelectorAll('.koma-nav-link[aria-current="page"]')
            .forEach(link => link.removeAttribute('aria-current'));
        if (!section) return;
        document.querySelector('#' + section + ' > .koma-nav-link')
            ?.setAttribute('aria-current', 'page');
    })()""")
    return content
end

function evaljs(w::BonitoWindow, code)
    isnothing(w.session[]) && return nothing
    return Bonito.evaljs(w.session[], code)
end

function handle(callback::Function, w::BonitoWindow, event::String)
    w.handlers[event] = callback
    return nothing
end

function show!(w::BonitoWindow)
    display = Bonito.use_electron_display(; options=w.window_options, devtools=w.dev_tools)
    w.display[] = display
    Base.display(display, w.app)
    return w
end

function Base.close(w::BonitoWindow)
    isnothing(w.session[]) || close(w.session[])
    if !isnothing(w.display[])
        display = w.display[]
        popdisplay(display)
        close(display)
        w.display[] = nothing
    end
    return nothing
end

function asset_stylesheet(session, filename)
    css = replace(read(filename, String), r"url\(([^)]+)\)" => expression -> begin
        match = Base.match(r"url\(([^)]+)\)", expression)
        reference = strip(only(match.captures), ['\'', '"'])
        path = first(split(reference, '?'))
        if startswith(path, "data:") || occursin("://", path)
            expression
        else
            asset = Asset(normpath(joinpath(dirname(filename), path)))
            "url(\"$(Bonito.url(session, asset))\")"
        end
    end)
    return DOM.style(css)
end

function home_page(image)
    return DOM.div(
        DOM.div(
            DOM.div(; class="mb-auto"),
            DOM.div(
                DOM.h1("Welcome to Koma!"),
                DOM.p(
                    "KomaMRI.jl is a Julia package to simulate Magnetic Resonance Imaging (MRI) acquisitions. The main focus of this package is to simulate general scenarios that could arise in pulse sequence development.";
                    class="lead",
                    style="padding:32px clamp(24px,8%,128px);",
                ),
                DOM.img(; src=image, alt="Icon", style="padding:32px 128px;width:326px;");
                class="mt-auto text-white-75",
            ),
            DOM.footer(
                DOM.p(
                    "Koma was developed by ",
                    DOM.a(
                        DOM.u("Carlos Castillo-Passi");
                        target="_blank",
                        href="https://github.com/cncastillo",
                        style="color:#2a7fb8;",
                    ),
                    ", et al.";
                    style="padding-bottom:10px;",
                ),
                DOM.p(
                    "Cite us ",
                    DOM.a(
                        DOM.u("doi:10.1002/mrm.29635");
                        target="_blank",
                        href="https://onlinelibrary.wiley.com/doi/10.1002/mrm.29635",
                        style="color:#2a7fb8;",
                    );
                    style="padding-bottom:10px;",
                ),
                DOM.p("Thanks to all our contributors!"; class="text-white mb-1"),
                DOM.a(
                    DOM.img(
                        ;
                        src="https://contrib.rocks/image?repo=cncastillo/KomaMRI.jl",
                        style="height:60px;",
                    );
                    href="https://github.com/cncastillo/KomaMRI.jl/graphs/contributors",
                );
                class="mt-auto text-white-50",
            );
            class="cover-container d-flex w-100 h-100 p-3 mx-auto flex-column mask",
        );
        class="d-flex h-100 text-center text-white bg-image",
    )
end

table_entries(values) = values
table_entries(values::AbstractDict) = sort!(collect(values); by=pair -> string(first(pair)))

function dictionary_page(values, title)
    rows = map(enumerate(table_entries(values))) do (index, (name, value))
        DOM.tr(
            DOM.th(string(index); scope="row"),
            DOM.td(string(name)),
            DOM.td(string(value)),
        )
    end
    return DOM.div(
        DOM.h1(title; class="koma-table-title"),
        DOM.table(
            DOM.thead(
                DOM.tr(
                    DOM.th("#"; scope="col"),
                    DOM.th("Name"; scope="col"),
                    DOM.th("Value"; scope="col"),
                ),
            ),
            DOM.tbody(rows...),
            class="table table-dark table-striped",
        ),
        class="koma-table-view",
    )
end

function setup_bonito_window(; darkmode=true, frame=true, dev_tools=false, versions="")
    komamri_src_ui = @__DIR__
    komamri_root = dirname(dirname(komamri_src_ui))
    scripts = joinpath(komamri_src_ui, "scripts")
    css = joinpath(komamri_src_ui, "css")

    style_assets = Asset.([
        joinpath(css, "bootstrap.min.css"),
        joinpath(css, "custom.css"),
        joinpath(css, "sidebars.css"),
    ])
    font_stylesheets = [
        joinpath(css, "bootstrap-icons.css"),
        joinpath(css, "katex.min.css"),
        joinpath(css, "icons.css"),
    ]
    script_assets = Asset.([
        joinpath(scripts, "custom.js"),
        joinpath(scripts, "bootstrap.bundle.min.js"),
        joinpath(scripts, "katex.min.js"),
        joinpath(scripts, "auto-render.min.js"),
    ])

    display_ref = Ref{Any}(nothing)
    session_ref = Ref{Union{Nothing,Session}}(nothing)
    content = Observable{Any}(DOM.div())
    state = Observable("index")
    events = Observable("")
    home = Ref{Any}(DOM.div())
    handlers = Dict{String,Function}()
    render_callbacks = Function[]
    on(events) do event
        callback = get(handlers, event, nothing)
        isnothing(callback) || callback(event)
        return nothing
    end

    app = App(; title="KomaUI") do session
        session_ref[] = session
        font_styles = asset_stylesheet.(Ref(session), font_stylesheets)
        logo = Bonito.url(session, Asset(joinpath(komamri_root, "assets", "logo-dark.svg")))
        home_image = Bonito.url(session, Asset(joinpath(komamri_root, "assets", "home-image.svg")))

        sidebar = read(joinpath(komamri_src_ui, "html", "sidebar.html"), String)
        sidebar = replace(
            sidebar,
            "LOGO" => string(logo),
            "title=\"Hooray!\"" => "title=\"$(replace(versions, '\n' => "&#10;"))\"",
        )
        home[] = home_page(home_image)
        state[] = "index"
        content[] = page_content(home[], "index")

        event_bridge = js"document.title = 'KomaUI'; window.KomaUI = {notify: name => $(events).notify(name)};"
        Bonito.on_document_load(session, js"renderMathInElement(document.body);")
        foreach(callback -> callback(session), render_callbacks)
        background = darkmode ? "background-color:rgb(13,16,17);" : ""
        main = DOM.main(
            DOM.div_unesc(sidebar),
            content,
            event_bridge;
            id="main",
            style=background,
        )
        return DOM.div(style_assets..., font_styles..., script_assets..., main)
    end

    window_options = Dict{String,Any}(
        "title" => "KomaUI",
        "autoHideMenuBar" => true,
        "frame" => frame,
        "icon" => joinpath(komamri_root, "assets", "app-icon.png"),
        "width" => 1200,
        "height" => 800,
    )
    return BonitoWindow(
        app,
        display_ref,
        session_ref,
        content,
        state,
        events,
        home,
        handlers,
        render_callbacks,
        window_options,
        dev_tools,
    )
end

function display_loading!(w::BonitoWindow, msg::String)
    loading = DOM.div(
        DOM.div(; class="koma-loading-spinner", role="status"),
        DOM.div(msg; id="loaddes", class="koma-loading-message");
        class="koma-loading",
    )
    return set_content!(w, loading, "loading")
end

function toast!(w::BonitoWindow, id, title, body)
    title_bytes = collect(codeunits(title))
    body_bytes = collect(codeunits(body))
    return evaljs(w, js"""
        const showToast = () => {
            if (typeof Toasty !== 'function') return requestAnimationFrame(showToast);
            const decode = bytes => new TextDecoder().decode(new Uint8Array(bytes));
            Toasty($(string(id)), decode($(title_bytes)), decode($(body_bytes)));
        };
        showToast();
    """)
end

function update_filename!(w, id, name)
    return evaljs(w, js"""
        const current = document.getElementById($(id));
        const label = document.createElement('abbr');
        label.title = $(name);
        label.textContent = $(name);
        current.replaceChildren(label);
    """)
end

function loaded_data_toast(w, name, view, next, noun)
    short_name = replace(display_filename(name, MAX_TOAST_FILENAME_CHARS), '&' => "&amp;", '<' => "&lt;", '>' => "&gt;")
    full_name = replace(name, '&' => "&amp;", '<' => "&lt;", '>' => "&gt;", '"' => "&quot;")
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

function callback_filepicker(filename::String, w::BonitoWindow, seq::Sequence)
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

function callback_filepicker(filename::String, w::BonitoWindow, obj::Phantom)
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

function callback_filepicker(filename::String, w::BonitoWindow, raw::RawAcquisitionData)
    raw = RawAcquisitionData(ISMRMRDFile(filename))
    if raw.params["systemVendor"] != "KomaMRI.jl"
        @warn "ISMRMRD files generated externally could cause problems during the reconstruction. We are currently improving compatibility."
    end
    name = basename(filename)
    update_filename!(w, "rawname", name)
    loaded_data_toast(w, name, "sig", "recon", "Raw data")
    return raw
end

function show_sequence!(w, seq, view; darkmode=true)
    if view === :sequence
        display_loading!(w, "Plotting sequence ...")
        long_seq = length(seq) > 1_000
        time_end = long_seq ? dur(seq) * 1e3 : 30
        plot = plot_seq(seq; darkmode, range=[0 time_end], slider=!long_seq, gl=long_seq, show_adc=false)
        set_content!(w, plot_node(plot), "sequence")
    elseif view === :kspace
        display_loading!(w, "Plotting kspace ...")
        set_content!(w, plot_node(plot_kspace(seq; darkmode)), "kspace")
    elseif view === :moment0
        display_loading!(w, "Plotting moment 0 ...")
        set_content!(w, plot_node(plot_M0(seq; darkmode)), "m0")
    elseif view === :moment1
        display_loading!(w, "Plotting moment 1 ...")
        set_content!(w, plot_node(plot_M1(seq; darkmode)), "m1")
    elseif view === :moment2
        display_loading!(w, "Plotting moment 2 ...")
        set_content!(w, plot_node(plot_M2(seq; darkmode)), "m2")
    else
        throw(ArgumentError("Unsupported sequence view: $view"))
    end
    return nothing
end

function show_phantom!(w, obj, buttons; key=:ρ, darkmode=true)
    display_loading!(w, "Plotting phantom ...")
    plot = plot_phantom_map(obj, key; time_samples=5, darkmode)
    return set_content!(
        w,
        DOM.div(DOM.div(buttons...), plot_node(plot); class="koma-plot-stack"),
        "phantom",
    )
end

function show_scanner!(w, sys)
    display_loading!(w, "Displaying scanner parameters ...")
    values = [
        "B0" => sys.B0,
        "B1" => sys.B1,
        "Gmax" => sys.Gmax,
        "Smax" => sys.Smax,
        "ADC_dt" => sys.ADC_Δt,
        "DUR_dt" => sys.DUR_Δt,
        "GR_dt" => sys.GR_Δt,
        "RF_dt" => sys.RF_Δt,
        "RF_ring_down_time" => sys.RF_ring_down_time,
        "RF_dead_time" => sys.RF_dead_time,
        "ADC_dead_time" => sys.ADC_dead_time,
    ]
    return set_content!(w, dictionary_page(values, "Scanner parameters"), "scanneparams")
end

function show_parameters!(w, parameters, title, state)
    display_loading!(w, "Displaying $(lowercase(title)) ...")
    return set_content!(w, dictionary_page(parameters, title), state)
end

function show_signal!(w, raw; darkmode=true)
    display_loading!(w, "Plotting raw signal ...")
    return set_content!(w, plot_node(plot_signal(raw; darkmode)), "sig")
end

function image_plot(img, type, darkmode)
    data, zmin, zmax = if type === :absi
        values = abs.(img) * prod(size(img)[1:2])
        values, minimum(values), maximum(values)
    elseif type === :angi
        values = angle.(img)
        values, -π, π
    elseif type === :absk
        values = log.(abs.(fftc(img)) .+ 1)
        values, 0, 0.1 * maximum(values)
    else
        throw(ArgumentError("Unsupported image view: $type"))
    end
    slider = Slider(1:size(img, 3))
    plot = map(slider.value) do slice
        plot_node(plot_image(data[:, :, slice]; zmin, zmax, darkmode, title="Reconstruction ($slice/$(size(img, 3)))"))
    end
    return size(img, 3) == 1 ? plot : DOM.div(slider, plot; class="koma-plot-stack")
end

function show_image!(w, img, view; darkmode=true)
    messages = Dict(:absi => "Plotting image magnitude ...", :angi => "Plotting image phase ...", :absk => "Plotting image k ...")
    display_loading!(w, messages[view])
    return set_content!(w, image_plot(img, view, darkmode), string(view))
end

function select_export_folder(w::BonitoWindow)
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

function save_ui!(w::BonitoWindow, seq::Sequence, obj::Phantom, sys::Scanner, raw, img, rec_params; type="all")
    mat_folder = select_export_folder(w)
    isnothing(mat_folder) && return nothing
    message = if type == "all"
        export_2_mat(seq, obj, sys, raw, rec_params, img, mat_folder; type, matfilename="")
    else
        export_2_mat(seq, obj, sys, raw, rec_params, img, mat_folder; type)
    end
    toast!(w, 1, "Saved .mat files", message)
    suffix = Dict("all" => "", "sequence" => "seq", "phantom" => "pha", "scanner" => "sca", "raw" => "raw", "image" => "ima")[type]
    state = "matfolder$suffix"
    w.state[] = state
    evaljs(w, js"document.getElementById('content').dataset.content = $(state);")
    return nothing
end
