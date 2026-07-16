struct KomaWindow
    app::App
    display::Base.RefValue{Any}
    window::Base.RefValue{Union{Nothing,Electron.Window}}
    session::Base.RefValue{Union{Nothing,Session}}
    content::Observable{Any}
    state::Observable{String}
    events::Observable{String}
    home::Base.RefValue{Any}
    handlers::Dict{String,Function}
    on_render::Vector{Function}
    listeners::Vector{ObserverFunction}
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

function set_content!(w::KomaWindow, content, state)
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

function evaljs(w::KomaWindow, code)
    session = w.session[]
    (isnothing(session) || !isopen(session)) && return nothing
    return Bonito.evaljs(session, code)
end

function handle(callback::Function, w::KomaWindow, event::String)
    w.handlers[event] = callback
    return nothing
end

function show!(w::KomaWindow)
    windows = Set(Iterators.flatten(Electron.windows.(Electron.applications())))
    display = Bonito.use_electron_display(; options=w.window_options, devtools=w.dev_tools)
    w.display[] = display
    w.window[] = only(
        setdiff(Iterators.flatten(Electron.windows.(Electron.applications())), windows)
    )
    Base.display(display, w.app)
    return w
end

function Base.isopen(w::KomaWindow)
    window = w.window[]
    return !isnothing(window) && isopen(window)
end

function Base.wait(w::KomaWindow)
    window = w.window[]
    isnothing(window) && return nothing
    for _ in Electron.msgchannel(window)
    end
    return nothing
end

function Base.close(w::KomaWindow)
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
    w.window[] = nothing
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
    listeners = ObserverFunction[]
    push!(listeners, on(events) do event
        callback = get(handlers, event, nothing)
        isnothing(callback) || callback(event)
        return nothing
    end)

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
    return KomaWindow(
        app,
        display_ref,
        Ref{Union{Nothing,Electron.Window}}(nothing),
        session_ref,
        content,
        state,
        events,
        home,
        handlers,
        render_callbacks,
        listeners,
        window_options,
        dev_tools,
    )
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

function toast!(w::KomaWindow, id, title, body)
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

function failure_toast!(w, id, operation, error)
    message = html_escape(sprint(showerror, error))
    return toast!(w, id, "$operation failed", "<pre class=\"mb-0\">$message</pre>")
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
