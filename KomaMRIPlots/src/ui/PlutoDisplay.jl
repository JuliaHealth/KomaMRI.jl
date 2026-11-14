# Preserve Plotly's return type; only Koma-created figures opt into Pluto rendering.
const KOMA_PLOTS = WeakKeyDict{PlotlyBase.Plot,Union{Nothing,KspaceSpatialSource}}()
function pluto_plot(plot, source=nothing)
    KOMA_PLOTS[plot] = source
    return plot
end

plotly_loader(asset) = """
window.komaPlotly ??= new Promise((resolve, reject) => {
    if (window.Plotly) return resolve(window.Plotly);
    const url = URL.createObjectURL(new Blob([$asset], {type:'text/javascript'}));
    const script = document.createElement('script');
    script.src = url;
    script.onload = () => { URL.revokeObjectURL(url); resolve(window.Plotly); };
    script.onerror = reject;
    document.head.append(script);
});
"""

function Base.show(io::Union{IOContext,IOStream,IOBuffer}, mime::MIME"text/html", plot::PlotlyBase.Plot; kwargs...)
    if haskey(KOMA_PLOTS, plot) && !isnothing(get(io, :pluto_published_to_js, nothing))
        return show_pluto_plot(io, plot)
    end
    source = get(KOMA_PLOTS, plot, nothing)
    if !isnothing(source) && !isnothing(source.overview)
        return show_live_plot(io, KspacePlot(source, plot))
    end
    return invoke(show, Tuple{IO,MIME"text/html",PlotlyBase.Plot}, io, mime, plot; kwargs...)
end

function PlotlyBase.to_html(io::Union{IOContext,IOStream,IOBuffer}, plot::PlotlyBase.Plot;
    include_plotlyjs="cdn", post_script=missing, kwargs...)
    source = get(KOMA_PLOTS, plot, nothing)
    if !isnothing(source) && !isnothing(source.overview)
        return show_live_plot(io, KspacePlot(source, plot))
    end
    if !haskey(KOMA_PLOTS, plot)
        return invoke(PlotlyBase.to_html, Tuple{IO,PlotlyBase.Plot}, io, plot;
            include_plotlyjs, post_script, kwargs...)
    end
    theme = PlotlyBase.JSON.json(read(joinpath(@__DIR__, "PlotTheme.js"), String))
    slices = any(key -> haskey(get(plot.layout,:meta,Dict()),key), ("phantom_slices","image_slices")) ?
        PlotlyBase.JSON.json(read(joinpath(@__DIR__, "PhantomSlices.js"), String)) : "null"
    controls = """
    {
    const themeURL = URL.createObjectURL(new Blob([$theme], {type:'text/javascript'}));
    import(themeURL).then(theme => {
        const graph = document.getElementById('$(plot.divid)');
        for (const dimension of ['width', 'height'])
            if (graph.layout[dimension]) graph.style[dimension] = graph.layout[dimension] + 'px';
        theme.observeControls(graph, Plotly)();
        theme.observeResize(graph, Plotly);
        if ($slices) {
            const url = URL.createObjectURL(new Blob([$slices], {type:'text/javascript'}));
            import(url).then(viewer => viewer.mountSlices(graph, Plotly))
                .finally(() => URL.revokeObjectURL(url));
        }
    }).finally(() => URL.revokeObjectURL(themeURL));
    }
    """
    bundled_plotly = isequal(include_plotlyjs, "cdn")
    if bundled_plotly
        asset = PlotlyBase.JSON.json(read(TIME_PLOTLY.local_path, String))
        print(io, "<script>", plotly_loader(asset), "</script>")
    end
    html = sprint(buffer -> invoke(PlotlyBase.to_html, Tuple{IO,PlotlyBase.Plot}, buffer, plot;
        include_plotlyjs=bundled_plotly ? "require-loaded" : include_plotlyjs,
        post_script=missing, kwargs...))
    # PlotlyBase's post_script hook references an undefined plot_id.
    hook = ismissing(post_script) ? controls : controls * replace(post_script, "{plot_id}" => plot.divid)
    html = replace(html, "\n    )" => "\n    ).then(function() { $hook })"; count=1)
    # Use PlotlyBase's require-loaded callback with the shared bundled asset promise.
    if bundled_plotly
        html = replace(html, "require([\"plotly\"], function(Plotly) {" =>
            "window.komaPlotly.then(function(Plotly) {"; count=1)
    end
    print(io, html)
    return nothing
end

function Bonito.jsrender(session::Bonito.Session, plot::PlotlyBase.Plot)
    source = get(KOMA_PLOTS, plot, nothing)
    return !isnothing(source) && !isnothing(source.overview) ?
        Bonito.jsrender(session, KspacePlot(source, plot)) :
        Bonito.render_mime(session, MIME"text/html"(), plot)
end

live_plot_app(plot) = Bonito.App(session -> Bonito.jsrender(session, plot))
const NOTEBOOK_PLOT_BOOTSTRAPS = WeakKeyDict{Bonito.Session,String}()
struct NotebookPlotCleanupPolicy{P<:Bonito.CleanupPolicy} <: Bonito.CleanupPolicy
    fallback::P
end
NotebookPlotCleanupPolicy(policy::NotebookPlotCleanupPolicy) = policy
Bonito.should_cleanup(policy::NotebookPlotCleanupPolicy, session) =
    !haskey(NOTEBOOK_PLOT_BOOTSTRAPS, session) && Bonito.should_cleanup(policy.fallback, session)
Bonito.allow_soft_close(policy::NotebookPlotCleanupPolicy) = Bonito.allow_soft_close(policy.fallback)
Bonito.set_cleanup_time!(policy::NotebookPlotCleanupPolicy, hours) =
    Bonito.set_cleanup_time!(policy.fallback, hours)
plot_layout(plot::TimePlot) = plot_template(plot.source).layout
plot_layout(plot::SpatialPlot) = plot.template.layout
plot_layout(plot::KspacePlot) = plot.plot.layout

Base.showable(::MIME"text/html", ::Union{TimePlot,SpatialPlot,KspacePlot}) = true
Base.show(io::IO, ::MIME"text/html", plot::Union{TimePlot,SpatialPlot,KspacePlot}) = show_live_plot(io, plot)
Base.display(display::Bonito.HTTPServer.BrowserDisplay, plot::Union{TimePlot,SpatialPlot,KspacePlot}) =
    Base.display(display, live_plot_app(plot))

function show_live_plot(io, plot)
    !isnothing(get(io, :pluto_published_to_js, nothing)) && return show_pluto_plot(io, plot)
    if haskey(io, :QuartoNotebookRunner)
        server = Bonito.get_server()
        route = "/koma-plot-$(Bonito.uuid4())"
        layout = plot_layout(plot)
        background = get(layout, :paper_bgcolor, "white")
        scheme = background == "white" ? "light" : "dark"
        app = Bonito.App(session -> Bonito.DOM.div(
            Bonito.DOM.style(":root { color-scheme:$scheme; } body { margin:0; background:$background; }"),
            Bonito.jsrender(session, plot)))
        Bonito.HTTPServer.route!(server, route => app)
        height = get(layout, :height, 400)
        print(io, "<iframe src=\"", Bonito.HTTPServer.online_url(server, route),
            "\" style=\"width:100%;height:$(height)px;border:0\"></iframe>")
        return nothing
    end
    app = live_plot_app(plot)
    parent = Bonito.CURRENT_SESSION[]
    if isnothing(parent) || Bonito.isclosed(parent)
        parent = Bonito.CURRENT_SESSION[] = Bonito.Session()
    end
    # Saved notebook outputs reconnect for as long as their Julia kernel lives.
    Bonito.set_cleanup_policy!(NotebookPlotCleanupPolicy(Bonito.CLEANUP_POLICY[]))
    bootstrap = get!(NOTEBOOK_PLOT_BOOTSTRAPS, parent) do
        reconnect = isready(parent)
        dom = Bonito.session_dom(parent, Bonito.App(nothing; indicator=app.indicator))
        if reconnect
            connection = sprint(Bonito.print_js_code, Bonito.setup_connection(parent), Bonito.JSSourceContext(parent))
            dom = Bonito.DOM.div(dom, Bonito.inline_code(parent, parent.asset_server, "if (!window.WEBSOCKET) { $connection }"))
        end
        Bonito.mark_displayed!(parent)
        sprint(show, MIME"text/html"(), dom)
    end
    # Every saved output must initialize the shared connection after a page reload.
    print(io, bootstrap)
    return Bonito.show_html(io, app; parent)
end

const PLUTO_PLOT_ASSETS = Dict{Symbol,Vector{UInt8}}()
const PLUTO_PLOT_SESSION = Ref{Union{Nothing,Bonito.Session{Bonito.PlutoConnection}}}(nothing)
const PLUTO_PLOT_BOOTSTRAP = Ref("")
function pluto_plot_assets()
    if isempty(PLUTO_PLOT_ASSETS)
        Bonito.bundle!(TIME_PLOT_JS)
        Bonito.bundle!(SPATIAL_PLOT_JS)
        merge!(PLUTO_PLOT_ASSETS, Dict(
            :plotly => read(TIME_PLOTLY.local_path),
            :time => read(joinpath(@__DIR__, "TimePlots.bundled.js")),
            :spatial => read(joinpath(@__DIR__, "SpatialPlots.bundled.js")),
            :theme => read(joinpath(@__DIR__, "PlotTheme.js")),
            :slices => read(joinpath(@__DIR__, "PhantomSlices.js")),
        ))
    end
    return PLUTO_PLOT_ASSETS
end

pluto_plot_spec(plot::PlotlyBase.Plot) = pluto_plot_spec(plot, get(KOMA_PLOTS, plot, nothing))

pluto_plot_spec(plot, ::Nothing) = (; mode="static", initial=PlotlyBase.JSON.json(plot),
    template="", overview="", query=nothing, layout=plot.layout, navigator=nothing, kspace="")

pluto_plot_spec(plot, source::KspaceSpatialSource) =
    (; pluto_plot_spec(plot, nothing)..., navigator=pluto_plot_spec(source.overview),
        kspace=PlotlyBase.JSON.json(kspace_controls(KspacePlot(source, plot))))

pluto_plot_spec(plot::KspacePlot) = pluto_plot_spec(plot.plot, plot.source)

pluto_plot_spec(plot::TimePlot) = (; plot_payload(plot)..., mode="time",
    template="", navigator=nothing, kspace="")

function pluto_plot_spec(plot::SpatialPlot)
    return (; plot_payload(plot)..., mode="spatial", overview="",
        navigator=isnothing(spatial_controls(plot.source)) ? nothing : pluto_plot_spec(plot.source.overview), kspace="")
end

function pluto_plot_dimensions(value)
    width = get(value.layout, :width, nothing)
    height = get(value.layout, :height, isnothing(value.navigator) ? 400 : 600)
    return (; value.mode, value.initial, value.template, value.overview, value.kspace,
        width=isnothing(width) ? "100%" : "$(width)px", height="$(height)px",
        navigator=isnothing(value.navigator) ? nothing : pluto_plot_dimensions(value.navigator))
end

function show_pluto_plot(io, plot)
    publish = get(io, :pluto_published_to_js, nothing)
    value = pluto_plot_spec(plot)
    transport_id = string(Bonito.uuid4())
    if !isnothing(value.query) || !isnothing(value.navigator)
        app = Bonito.App(; indicator=nothing) do session
            request, response = isnothing(value.query) ? (nothing, nothing) : plot_observables(session, value)
            navigator_request, navigator_response = isnothing(value.navigator) ?
                (nothing, nothing) : plot_observables(session, value.navigator)
            node = Bonito.DOM.div(; id=transport_id, style="display:none")
            Bonito.onload(session, node, js"""node => {
                node.transport = {request:$(request), response:$(response),
                    navigator:{request:$(navigator_request), response:$(navigator_response)}};
                node.dispatchEvent(new Event('koma-transport-ready'));
            }""")
            return node
        end
        parent = PLUTO_PLOT_SESSION[]
        # Cached cell outputs must initialize the shared connection again after a page reload.
        if isnothing(parent) || Bonito.isclosed(parent)
            parent = PLUTO_PLOT_SESSION[] = Bonito.Session(Bonito.PlutoConnection())
            dom = Bonito.session_dom(parent, Bonito.App(nothing; indicator=nothing))
            PLUTO_PLOT_BOOTSTRAP[] = sprint(show, MIME"text/html"(), dom)
            Bonito.mark_displayed!(parent)
        end
        print(io, "<script id='koma-pluto-transport'>\nconst html = ")
        publish(io, PLUTO_PLOT_BOOTSTRAP[])
        print(io, ";\nconst id = '", parent.id, "';", raw"""
const connections = window.komaPlutoConnections ??= {};
connections[id] ??= (async () => {
    const container = document.createElement('div');
    container.hidden = true;
    container.innerHTML = html;
    const scripts = [...container.querySelectorAll('script')];
    document.body.appendChild(container);
    for (const original of scripts) {
        const script = document.createElement('script');
        for (const attribute of original.attributes) script.setAttribute(attribute.name, attribute.value);
        script.textContent = original.textContent;
        const loaded = new Promise((resolve, reject) => {script.onload = resolve; script.onerror = reject;});
        original.replaceWith(script);
        if (script.src || script.type === 'module') await loaded;
    }
})();
await connections[id];
</script>
""")
        Bonito.show_html(io, app; parent)
    end
    print(io, "<script id='koma-pluto-display'>\nconst assets = ")
    publish(io, pluto_plot_assets())
    print(io, ";\nconst spec = ")
    publish(io, pluto_plot_dimensions(value))
    print(io, ";\nconst transportNode = document.getElementById('", transport_id, "')")
    print(io, raw""";
const kind = spec.mode;
const root = this?.dataset.komaPlotKind === kind ? this : document.createElement('div');
this?.bridge?.cancel?.();
this?.navigator?.bridge?.cancel?.();
const revision = root.revision = (root.revision ?? 0) + 1;
let cancelled = false;
invalidation.then(() => { cancelled = true; });
const current = () => !cancelled && revision === root.revision;
root.dataset.komaPlotKind = kind;
if (!root.isConnected) currentScript.before(root);
const transport = transportNode && await Promise.race([
    transportNode.transport ?? new Promise(resolve => transportNode.addEventListener(
        'koma-transport-ready', () => resolve(transportNode.transport), {once:true})),
    invalidation.then(() => null)
]);
if (!current()) return root;
""", plotly_loader("assets.plotly"), raw"""
const Plotly = await window.komaPlotly;
await root.renderTask;
if (!current()) return root;
const cache = window.komaPlutoAssets ??= {};
function viewerModule(mode) {
    return cache[mode] ??= (async () => {
        const url = URL.createObjectURL(new Blob([assets[mode]], {type:'text/javascript'}));
        try { return await import(url); }
        finally { URL.revokeObjectURL(url); }
    })();
}
async function render(root, spec, transport) {
root.className = `koma-pluto-display${spec.mode === 'time' ? ' koma-time-plot' : ''}`;
Object.assign(root.style, {width:spec.width, height:spec.height, minHeight:'0', position:'relative'});
if (spec.navigator) {
    root.dataset.linkedTimePlots = 'true';
    Object.assign(root.style, {display:'flex', flexDirection:'column'});
    const navigator = root.navigator ??= root.appendChild(document.createElement('div'));
    await render(navigator, spec.navigator, transport.navigator);
    navigator.style.cssText += ';max-height:45%;flex:none';
} else if (root.navigator) {
    root.navigator.remove();
    root.navigator = null;
    delete root.dataset.linkedTimePlots;
    root.style.display = 'block';
}
if (spec.mode === 'static') {
    const graph = root.graph ??= document.createElement('div');
    if (!graph.isConnected) root.prepend(graph);
    graph.style.cssText = spec.navigator ? 'width:100%;flex:1;min-height:0' : 'width:100%;height:100%';
    const plot = JSON.parse(spec.initial);
    const theme = await viewerModule('theme');
    plot.layout.width = graph.clientWidth;
    plot.layout.height = graph.clientHeight;
    root.definition = await theme.renderPlot(graph, Plotly, plot, root.definition);
    root.controls ??= theme.observeControls(graph, Plotly);
    await root.controls();
    const sliceModel = plot.layout.meta?.phantom_slices ?? plot.layout.meta?.image_slices;
    if (sliceModel) {
        const viewer = await viewerModule('slices');
        if (!root.slices) root.slices = viewer.mountSlices(graph, Plotly);
        else await root.slices.rebind(sliceModel);
    } else if (root.slices) {
        root.slices.dispose(); root.slices = null;
    }
    if (spec.kspace) {
        const viewer = await viewerModule('spatial');
        if (!root.kspace) root.kspace = await viewer.mountKspace(root, graph, Plotly, JSON.parse(spec.kspace));
        else await root.kspace.rebind(JSON.parse(spec.kspace));
    } else if (root.kspace) {
        root.kspace.dispose();
        root.kspace = null;
    }
    root.resize ??= theme.observeResize(graph, Plotly);
} else {
    const viewer = await viewerModule(spec.mode);
    if (!current()) return;
    const bridge = root.bridge ??= {generation:0, request:null, callbacks:[]};
    const generation = ++bridge.generation;
    bridge.request = transport.request;
    const cancel = () => {
        if (generation !== bridge.generation) return;
        bridge.generation++;
        bridge.request = null;
    };
    bridge.cancel = cancel;
    transport.response.on(value => {
        if (generation !== bridge.generation) return false;
        bridge.callbacks.forEach(callback => callback(value));
    });
    invalidation.then(cancel);
    const response = {value:spec.initial, on:callback => bridge.callbacks.push(callback)};
    const request = {notify:q => bridge.request?.notify(q)};
    if (!root.viewer) root.viewer = await viewer.mount(root, Plotly, request, response,
        spec.mode === 'time' ? spec.overview : spec.template);
    else await root.viewer.rebind(spec.initial, spec.mode === 'time' ? spec.overview : spec.template);
}
}
root.renderTask = render(root, spec, transport);
await root.renderTask;
return root;
</script>
""")
    return nothing
end
