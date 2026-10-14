"""
    KomaWorkspace([w1, w2, ...]; host="127.0.0.1", port=9384)
    KomaWorkspace(; host="127.0.0.1", port=9384)

Show independent KomaUI windows in a browser grid on one server. Tile numbers
belong to this workspace and stay stable when a tile is expanded. Each member
keeps its own data, device, controls, and operation status.

Create members with `KomaUI(; show_window=false, return_window=true)`, or adopt
one existing browser window without restarting its server. Desktop windows and
members already attached to a workspace cannot be adopted. The no-argument
form starts with one window; use Add window or `add_window!` for more.

`workspace.windows` contains the Julia handles. Use the tile's Close button to
close its window and remove it from the workspace. `close(w)` closes only that window;
`close(workspace)` closes all members and the server. Forward the workspace's
single printed port when using SSH.
"""
struct KomaWorkspace
    app::App
    windows::Vector{KomaWindow}
    server::Bonito.Server
    entries::Observable{Vector{Tuple{Int,String}}}
    darkmode::Observable{Bool}
end

const WORKSPACE_SCRIPT = Bonito.ES6Module(joinpath(@__DIR__, "scripts", "workspace.js"))

function KomaWorkspace(windows::AbstractVector{KomaWindow}; host="127.0.0.1", port=9384)
    allunique(windows) || throw(ArgumentError("Each workspace tile needs a different KomaWindow."))
    all(w -> isnothing(w.window[]) && !haskey(WORKSPACE_ROUTES, w), windows) ||
        throw(ArgumentError("Use unshown or standalone browser windows for workspace members."))
    servers = unique([w.server[] for w in windows if !isnothing(w.server[])])
    length(servers) <= 1 || throw(ArgumentError("Adopt at most one existing browser server; create other members with show_window=false."))
    entries = Observable(Tuple{Int,String}[])
    add = Observable(0)
    close_window = Observable(0)
    darkmode = Observable(isempty(windows) || first(windows).darkmode[])
    workspace_ref = Ref{KomaWorkspace}()
    session_ref = Ref{Union{Nothing,Session}}(nothing)
    app = App(; title="Koma workspace") do session
        isnothing(session_ref[]) || session_ref[] === session || close(session_ref[])
        session_ref[] = session
        on(session, add) do value
            iszero(value) && return
            add[] = 0
            add_window!(workspace_ref[])
        end
        on(session, close_window) do id
            iszero(id) && return
            close_window[] = 0
            workspace = workspace_ref[]
            index = findfirst(entry -> first(entry) == id, workspace.entries[])
            isnothing(index) && return
            w = workspace.windows[index]
            trylock(w.operation_lock) || return
            try
                w.operation[].state === :running && return
                close(w)
                deleteat!(workspace.windows, index)
                workspace.entries[] = filter(entry -> first(entry) != id, workspace.entries[])
            finally
                unlock(w.operation_lock)
            end
        end
        on(session, darkmode) do dark
            foreach(w -> w.darkmode[] = dark, workspace_ref[].windows)
        end
        stylesheet = read(joinpath(@__DIR__, "css", "workspace.css"), String)
        return DOM.div(
            asset_stylesheet(session, joinpath(@__DIR__, "css", "bootstrap-icons.css")),
            DOM.style(stylesheet),
            DOM.header(
                DOM.div(
                    DOM.img(; src=Bonito.url(session, Asset(joinpath(pkgdir(KomaMRI), "assets", "app-icon.svg"))),
                        width="24", height="24", alt=""),
                    DOM.strong("Koma workspace"); class="workspace-brand"),
                DOM.div(
                    DOM.label("Columns ", DOM.select(
                        DOM.option("1"; value="1"), DOM.option("2"; value="2", selected=true),
                        DOM.option("4"; value="4"); id="workspace-columns")),
                    DOM.button("+ Add window"; id="workspace-add", type="button");
                    class="workspace-layout"),
                DOM.div(
                    DOM.details(
                        DOM.summary(DOM.i(; class="bi bi-github", ariaHidden="true");
                            ariaLabel="KomaMRI links", title="KomaMRI links"),
                        DOM.div(
                            DOM.a("Go to Koma's GitHub"; href="https://github.com/cncastillo/KomaMRI.jl/", target="_blank", rel="noopener"),
                            DOM.a("Documentation"; href="https://cncastillo.github.io/KomaMRI.jl/", target="_blank", rel="noopener"),
                            DOM.a("Post an issue/bug"; href="https://github.com/cncastillo/KomaMRI.jl/issues/new", target="_blank", rel="noopener");
                            class="workspace-links-menu"); class="workspace-links"),
                    DOM.label(
                        DOM.input(; type="checkbox", role="switch", ariaLabel="Dark mode", id="workspace-theme"),
                        DOM.i(; class="bi bi-moon-stars-fill", ariaHidden="true");
                        class="workspace-theme"),
                    DOM.span(; id="workspace-connection", role="status");
                    class="workspace-tools");
                class="workspace-toolbar"),
            DOM.div(; id="workspace-grid", class="workspace-grid"),
            js"""$(WORKSPACE_SCRIPT).then(module => module.mount($(entries), $(add), $(close_window), $(darkmode)))""";
            class="koma-workspace")
    end
    server = isempty(servers) ? Bonito.Server(host, port) : only(servers)
    workspace = KomaWorkspace(app, KomaWindow[], server, entries, darkmode)
    workspace_ref[] = workspace
    for w in windows
        register_workspace_window!(workspace, w)
    end
    Bonito.route!(server, "/" => app)
    @info "Open Koma workspace in your browser" url=Bonito.HTTPServer.online_url(server, "/")
    return workspace
end

KomaWorkspace(; kwargs...) = KomaWorkspace(
    [KomaUI(; show_window=false, return_window=true)]; kwargs...)

function register_workspace_window!(workspace, w)
    w.darkmode[] = workspace.darkmode[]
    id = isempty(workspace.entries[]) ? 1 : maximum(first, workspace.entries[]) + 1
    route = "/window/$id"
    w.server[] = workspace.server
    WORKSPACE_ROUTES[w] = route
    Bonito.route!(workspace.server, route => w.app)
    push!(workspace.windows, w)
    workspace.entries[] = [workspace.entries[]; [(id, route)]]
    return nothing
end

"""
    add_window!(workspace; kwargs...)

Create an independent workspace member and return its `KomaWindow`. Keywords
configure its inputs and options as in `KomaUI`; the workspace owns its display.
Adding a member does not reload existing tiles.
"""
function add_window!(workspace::KomaWorkspace; kwargs...)
    isopen(workspace) || error("Koma workspace is closed.")
    w = KomaUI(; show_window=false, return_window=true, kwargs...)
    register_workspace_window!(workspace, w)
    return w
end

Base.isopen(workspace::KomaWorkspace) = Bonito.HTTPServer.isrunning(workspace.server)
Base.wait(workspace::KomaWorkspace) = wait(workspace.server)

function Base.close(workspace::KomaWorkspace)
    foreach(close, workspace.windows)
    isnothing(workspace.app.session[]) || close(workspace.app.session[])
    close(workspace.app)
    close(workspace.server)
    return nothing
end
