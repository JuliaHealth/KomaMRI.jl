struct KomaContent
    value::Observable{Any}
end

Base.getindex(content::KomaContent) = content.value[]
Base.setindex!(content::KomaContent, value) = content.value[] = value

function Bonito.jsrender(session::Session, content::KomaContent)
    pkgversion(Bonito) >= v"5" && return Bonito.jsrender(session, content.value)

    # Bonito 4 closes displaced subsessions before pending asset requests complete.
    root = DOM.div()
    current, html = Bonito.render_subsession(session, content[]; init=true)
    Bonito.mark_displayed!(current)
    previous = nothing
    Bonito.on(session, content.value) do value
        next = Bonito.update_session_dom!(
            session, Bonito.uuid(session, root), value; replace=false
        )
        if next !== current
            isnothing(previous) || close(previous)
            previous = current
            current = next
        end
        return nothing
    end
    push!(Bonito.children(root), html)
    return Bonito.jsrender(session, root)
end

Base.@kwdef struct UIOperation
    action::Symbol = :none
    state::Symbol = :idle
    progress::Union{Nothing,Float64} = nothing
    elapsed_seconds::Union{Nothing,Float64} = nothing
    eta_seconds::Union{Nothing,Float64} = nothing
    message::String = ""
end

struct KomaWindow
    app::App
    display::Base.RefValue{Any}
    window::Base.RefValue{Union{Nothing,Electron.Window}}
    session::Base.RefValue{Union{Nothing,Session}}
    content::KomaContent
    state::Observable{String}
    events::Observable{String}
    darkmode::Observable{Bool}
    home::Base.RefValue{Any}
    handlers::Dict{String,Function}
    files::Dict{Symbol,String}
    on_render::Vector{Function}
    listeners::Vector{ObserverFunction}
    window_options::Dict{String,Any}
    dev_tools::Bool
    sys::Observable{Scanner}
    seq::Observable{Sequence}
    obj::Observable{Phantom}
    physio::Observable{AbstractPhysioSignal}
    raw::Observable{RawAcquisitionData}
    img::Observable{Union{AbstractArray,ReconstructionResult}}
    sim_params::Observable{Dict{String,Any}}
    rec_params::Observable{Dict{Symbol,Any}}
    operation::Observable{UIOperation}
    operation_lock::ReentrantLock
    server::Base.RefValue{Union{Nothing,Bonito.Server}}
    export_assets::Base.RefValue{Union{Nothing,Bonito.ChildAssetServer}}
end

const WORKSPACE_ROUTES = IdDict{KomaWindow,String}()
window_route(w) = get(WORKSPACE_ROUTES, w, "/")

Base.show(io::IO, ::KomaWindow) = print(io, "KomaWindow()")

function page_content(content, state)
    return DOM.div(
        content;
        id="content",
        class="koma-content",
        dataContent=state,
    )
end

function set_content!(w::KomaWindow, content, state)
    w.state[] = state
    w.content[] = page_content(content, state)
    evaljs(w, js"""(() => {
        const section = {
            sequence: 'pulses', seqdefs: 'pulses', kspace: 'pulses', m0: 'pulses', m1: 'pulses', m2: 'pulses',
            slew_rate: 'pulses',
            phantom: 'phantom', coils: 'scanner', scanneparams: 'scanner', sig: 'sig',
            absi: 'recon', absk: 'recon'
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

function evaljs(w::KomaWindow, code)
    session = w.session[]
    (isnothing(session) || !isopen(session)) && return nothing
    return Bonito.evaljs(session, code)
end

function setup_host_status!(w)
    push!(w.on_render, session -> begin
        location = haskey(ENV, "SSH_CONNECTION") ? gethostname() : "Local"
        port = isnothing(w.server[]) ? "" : ":$(w.server[].port)"
        location *= port
        devices = Observable(get_devices())
        choice = Observable(0)
        on(session, choice) do index
            if index == 0
                devices[] = get_devices()
            else
                ui_action(w, :select_device) do
                    device = devices[][index]
                    click!(w, :select_device; gpu=device.backend != "CPU", id=device.id)
                end
            end
        end
        status = map(session, w.sim_params) do params
            "$(w.app.title)\nJulia host: $(gethostname())$port\nUsing $(simulation_device_label(params))"
        end
        options = map(session, w.sim_params, w.operation, devices) do params, operation, available
            backend = KomaMRICore.name(KomaMRICore.get_backend(get(params, "gpu", true); verbose=false))
            selected = get(params, "gpu_device", nothing)
            threads = get(params, "Nthreads", Threads.nthreads())
            rows = [(label=d.backend == "CPU" ? "CPU — $(d.name) · $threads thread$(threads == 1 ? "" : "s")" : "GPU $(d.id) — $(d.name)",
                short_label=d.backend == "CPU" ? "CPU" : "GPU$(d.id)",
                selected=d.backend == backend && (backend == "CPU" || (isnothing(selected) ? d.current : d.id == selected))) for d in available]
            (; rows, busy=operation.state === :running)
        end
        workspace_status = map(session, w.operation) do operation
            Dict{String,Any}(
                "action" => string(operation.action), "state" => string(operation.state),
                "progress" => operation.progress, "elapsed_seconds" => operation.elapsed_seconds,
                "eta_seconds" => operation.eta_seconds)
        end
        Bonito.on_document_load(session, js"""(() => {
            const host = document.getElementById('koma-host-status');
            const toggle = document.getElementById('koma-device-toggle');
            const menu = document.getElementById('koma-device-menu');
            const selectedDevice = document.getElementById('koma-selected-device');
            const report = operation => {
                if (window.parent === window) return;
                window.parent.postMessage({type: 'koma-workspace-status',
                    busy: $(options).value.busy,
                    connection: document.getElementById('main').dataset.connection,
                    ...operation}, window.location.origin);
            };
            $(workspace_status).on(report);
            new MutationObserver(() => report($(workspace_status).value))
                .observe(document.getElementById('main'),
                    {attributes: true, attributeFilter: ['data-connection']});
            new bootstrap.Dropdown(toggle, {offset: [0, 8], popperConfig: {placement: 'right-end'}});
            host.classList.toggle('is-ssh', $(haskey(ENV, "SSH_CONNECTION")));
            host.querySelector('.koma-nav-label').textContent = $(location);
            const updateDevices = value => {
                toggle.setAttribute('aria-disabled', String(value.busy));
                selectedDevice.textContent = value.rows.find(device => device.selected)?.short_label ?? '';
                document.title = [$(w.app.title), $(location), selectedDevice.textContent].filter(Boolean).join(' · ');
                report($(workspace_status).value);
                menu.replaceChildren(...value.rows.map((device, index) => {
                    const item = document.createElement('li');
                    const button = document.createElement('button');
                    button.type = 'button';
                    button.className = 'dropdown-item';
                    button.classList.toggle('active', device.selected);
                    button.setAttribute('aria-pressed', String(device.selected));
                    button.textContent = device.label;
                    button.disabled = value.busy;
                    button.onclick = () => {
                        bootstrap.Dropdown.getInstance(toggle).hide();
                        $(choice).notify(index + 1);
                    };
                    item.append(button);
                    return item;
                }));
            };
            $(options).on(updateDevices);
            updateDevices($(options).value);
            host.addEventListener('show.bs.dropdown', event => {
                if ($(options).value.busy || document.getElementById('main').dataset.connection !== 'connected') {
                    event.preventDefault();
                    return;
                }
                bootstrap.Tooltip.getInstance(host)?.hide();
                $(choice).notify(0);
            });
            const describe = value => {
                host.setAttribute('aria-label', value);
                toggle.setAttribute('aria-label', value + ' — Choose simulation device');
            };
            $(status).on(describe);
            describe($(status).value);
            host.addEventListener('show.bs.tooltip', event => {
                if (toggle.getAttribute('aria-expanded') === 'true') event.preventDefault();
            });
            new bootstrap.Tooltip(host, {
                title: () => $(status).value,
                trigger: 'hover focus', placement: 'right', container: host,
                customClass: 'koma-host-tooltip'
            });
        })()""")
    end)
    return nothing
end

function show!(w::KomaWindow)
    display = Bonito.use_electron_display(; options=w.window_options, devtools=w.dev_tools)
    w.display[] = display
    w.window[] = display.window.window
    run(w.window[].app, string(js"""
        electron.app.setName('KomaUI');
        if (process.platform === 'darwin') {
            const icon = electron.nativeImage.createFromPath($(w.window_options["icon"]));
            const {width, height} = icon.getSize();
            // A 10% transparent margin on each side matches macOS Dock icon sizing.
            const padding = Math.round(Math.max(width, height) / 8);
            const paddedWidth = width + 2 * padding;
            const paddedHeight = height + 2 * padding;
            const pixels = icon.toBitmap();
            const padded = Buffer.alloc(paddedWidth * paddedHeight * 4);
            for (let row = 0; row < height; row++) {
                pixels.copy(padded, ((row + padding) * paddedWidth + padding) * 4,
                    row * width * 4, (row + 1) * width * 4);
            }
            electron.app.dock.setIcon(electron.nativeImage.createFromBitmap(padded,
                {width: paddedWidth, height: paddedHeight}));
        }
    """))
    Base.display(display, w.app)
    return w
end

function Base.isopen(w::KomaWindow)
    isnothing(w.server[]) || return Bonito.HTTPServer.isrunning(w.server[])
    window = w.window[]
    return !isnothing(window) && isopen(window)
end

function Base.wait(w::KomaWindow)
    if haskey(WORKSPACE_ROUTES, w)
        while isopen(w)
            sleep(0.1)
        end
        return nothing
    end
    isnothing(w.server[]) || return wait(w.server[])
    window = w.window[]
    isnothing(window) && return nothing
    # Electron closes this channel when its window closes.
    foreach(_ -> nothing, Electron.msgchannel(window))
    return nothing
end

function Base.close(w::KomaWindow)
    window = w.window[]
    foreach(off, w.listeners)
    empty!(w.listeners)
    if !isnothing(w.session[])
        close(w.session[])
        w.session[] = nothing
    end
    if !isnothing(w.display[])
        display = w.display[]
        popdisplay(display)
        close(display)
        w.display[] = nothing
    end
    if !isnothing(window) && isopen(window.app)
        close(window.app)
    end
    w.window[] = nothing
    if !isnothing(w.export_assets[])
        close(w.export_assets[])
        w.export_assets[] = nothing
    end
    if !isnothing(w.server[])
        route = pop!(WORKSPACE_ROUTES, w, nothing)
        isnothing(route) ? close(w.server[]) : Bonito.delete_route!(w.server[], route)
        w.server[] = nothing
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
                class="mt-auto",
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
                    class="mb-1",
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
                DOM.p("Thanks to all our contributors!"; class="mb-1"),
                DOM.a(
                    DOM.img(
                        ;
                        src="https://contrib.rocks/image?repo=cncastillo/KomaMRI.jl",
                        style="height:58px;",
                    );
                    href="https://github.com/cncastillo/KomaMRI.jl/graphs/contributors",
                );
                class="mt-auto",
            );
            class="cover-container d-flex w-100 h-100 p-3 mx-auto flex-column mask",
        );
        class="koma-home d-flex h-100 text-center bg-image",
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
            class="table table-striped",
        ),
        class="koma-table-view",
    )
end

function setup_bonito_window(sys, seq, obj, sim_params, rec_params;
    darkmode=true, frame=true, dev_tools=false, versions="",
)
    title = "KomaUI"
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
    content = KomaContent(Observable{Any}(DOM.div()))
    state = Observable("index")
    events = Observable("")
    theme = Observable(darkmode)
    home = Ref{Any}(DOM.div())
    handlers = Dict{String,Function}()
    render_callbacks = Function[]
    listeners = ObserverFunction[]

    app = App(; title) do session
        isnothing(session_ref[]) || session_ref[] === session || close(session_ref[])
        session_ref[] = session
        session_theme = Bonito.bind_global(session, theme)
        # Plot modules must outlive views displaced while their imports are pending.
        foreach(asset -> Bonito.url(session, asset),
            (PLOTLY_ASSET, KomaMRIPlots.TIME_PLOT_JS, KomaMRIPlots.SPATIAL_PLOT_JS, KomaMRIPlots.SLICE_PLOT_JS))
        font_styles = asset_stylesheet.(Ref(session), font_stylesheets)
        logos = [Bonito.url(session, Asset(joinpath(komamri_root, "assets", name)))
            for name in ("logo.svg", "logo-dark.svg")]
        spin_icon = Bonito.url(session, Asset(joinpath(komamri_root, "assets", "app-icon.svg")))
        home_image = Bonito.url(session, Asset(joinpath(komamri_root, "assets", "home-image.svg")))

        sidebar = read(joinpath(komamri_src_ui, "html", "sidebar.html"), String)
        sidebar = replace(
            sidebar,
            "LOGO" => string(logos[theme[] ? 2 : 1]),
            "SPIN_ICON" => string(spin_icon),
            "CSS_VERSION" => string(hash(read(joinpath(css, "custom.css")))),
            "title=\"Hooray!\"" => "title=\"$(replace(versions, '\n' => "&#10;"))\"",
        )
        home[] = home_page(home_image)
        state[] == "index" && (content[] = page_content(home[], "index"))

        event_bridge = js"""(() => {
            window.KomaUI = {notify: name => $(events).notify(name)};
            const applyTheme = dark => {
                const main = document.getElementById('main');
                main.className = dark ? 'koma-dark' : 'koma-light';
                main.dataset.plotTheme = JSON.stringify($(KomaMRIPlots.theme_chooser.((false, true)))[dark ? 1 : 0]);
                document.querySelector('.koma-logo-full').src = $(logos)[dark ? 1 : 0];
                const toggle = document.getElementById('button_toggle_theme');
                toggle.checked = dark;
                toggle.parentElement.title = dark ? 'Switch to light mode' : 'Switch to dark mode';
            };
            $(session_theme).on(applyTheme);
            applyTheme($(session_theme).value);
        })()"""
        Bonito.on_document_load(session, js"renderMathInElement(document.body);")
        foreach(callback -> callback(session), render_callbacks)
        main = DOM.main(
            DOM.div_unesc(sidebar),
            DOM.div(
                DOM.header(; class="koma-window-titlebar"),
                DOM.div(content; class="koma-content-container");
                class="koma-window-panel"),
            event_bridge;
            id="main",
            dataPlotTheme=PlotlyBase.JSON.json(KomaMRIPlots.theme_chooser(theme[])),
            class=theme[] ? "koma-dark" : "koma-light",
        )
        return DOM.div(PLOT_THEME, style_assets..., font_styles..., script_assets..., main)
    end

    window_options = Dict{String,Any}(
        "title" => title,
        "autoHideMenuBar" => true,
        "frame" => frame,
        "icon" => joinpath(komamri_root, "assets", "app-icon.png"),
        "width" => 1200,
        "height" => 800,
    )
    w = KomaWindow(
        app,
        display_ref,
        Ref{Union{Nothing,Electron.Window}}(nothing),
        session_ref,
        content,
        state,
        events,
        theme,
        home,
        handlers,
        Dict{Symbol,String}(),
        render_callbacks,
        listeners,
        window_options,
        dev_tools,
        Observable{Scanner}(sys),
        Observable{Sequence}(seq),
        Observable{Phantom}(obj),
        Observable{AbstractPhysioSignal}(default_physio_signal(seq)),
        Observable{RawAcquisitionData}(setup_raw()),
        Observable{Union{AbstractArray,ReconstructionResult}}([0.0im 0.; 0. 0.]),
        Observable{Dict{String,Any}}(sim_params),
        Observable{Dict{Symbol,Any}}(rec_params),
        Observable(UIOperation()),
        ReentrantLock(),
        Ref{Union{Nothing,Bonito.Server}}(nothing),
        Ref{Union{Nothing,Bonito.ChildAssetServer}}(nothing),
    )
    push!(listeners, on(events) do event
        @async ui_action(w, event) do
            click!(w, event)
        end
        return nothing
    end)
    return w
end

function display_loading!(w::KomaWindow, msg::String; details="")
    loading = DOM.div(
        DOM.div(; class="koma-loading-spinner", role="status"),
        DOM.div(
            DOM.div(msg),
            DOM.div(details; id="loadstatus", class="koma-loading-status");
            id="loaddes",
            class="koma-loading-message",
        );
        class="koma-loading",
    )
    return set_content!(w, loading, "loading")
end

function restore_content!(w, content, state)
    w.state[] = state
    w.content[] = content
    return nothing
end

html_escape(text) = replace(
    string(text), '&' => "&amp;", '<' => "&lt;", '>' => "&gt;", '"' => "&quot;"
)

function toast!(w::KomaWindow, id, title, body; level=:info)
    color = level == :warning ? "#ffc107" : level == :error ? "#dc3545" :
        id == 0 ? "currentColor" : "#2a7fb8"
    title_bytes = collect(codeunits(title))
    body_bytes = collect(codeunits(body))
    return evaljs(w, js"""
        const showToast = () => {
            if (typeof Toasty !== 'function') return requestAnimationFrame(showToast);
            const decode = bytes => new TextDecoder().decode(new Uint8Array(bytes));
            document.querySelector('#liveToast' + $(string(id)) + ' rect').setAttribute('fill', $(color));
            Toasty($(string(id)), decode($(title_bytes)), decode($(body_bytes)));
        };
        showToast();
    """)
end

function failure_toast!(w, id, operation, error)
    message = html_escape(sprint(showerror, error))
    return toast!(w, id, "$operation failed", "<div class=\"text-break\" style=\"white-space: pre-wrap;\">$message</div>"; level=:error)
end

function update_filename!(w, id, name)
    return evaljs(w, js"""
        const current = document.getElementById($(id));
        const label = document.createElement('abbr');
        label.title = $(name);
        label.textContent = $(name);
        current.replaceChildren(label);
        const section = current.previousElementSibling;
        const reload = section.querySelector('.file-reload-button');
        if (reload) reload.disabled = false;
        const toggle = section.querySelector('.koma-nav-link');
        toggle.title = toggle.getAttribute('aria-label') + ': ' + $(name);
        const picker = section.querySelector('.koma-file-input');
        if (picker) {
            const caption = picker.querySelector('.koma-file-name');
            caption.textContent = $(display_filename(name, MAX_UI_FILENAME_CHARS));
            caption.title = picker.title = $(name);
        }
    """)
end
