const PLOTLY_ASSET = Asset(joinpath(artifact"plotly-artifacts", "plotly.min.js"); name="Plotly")
const PLOT_THEME = Bonito.ES6Module(joinpath(pkgdir(KomaMRIPlots), "src", "ui", "PlotTheme.js"))

plot_node(plot::Union{KomaMRIPlots.TimePlot,KomaMRIPlots.SpatialPlot,KomaMRIPlots.KspacePlot}) =
    DOM.div(plot; style="width:100%;height:100%;padding:8px;box-sizing:border-box;")

plot_node(plot::PlotlyBase.Plot; fit_colorbar=false) =
    plot_node(plot, get(KomaMRIPlots.KOMA_PLOTS, plot, nothing); fit_colorbar)

plot_node(plot, ::Nothing; fit_colorbar=false) = plot_node(Observable(plot); fit_colorbar)
plot_node(plot, source::KomaMRIPlots.KspaceSpatialSource; fit_colorbar=false) =
    plot_node(KomaMRIPlots.KspacePlot(source, plot))

function plot_node(plot::Observable; fit_colorbar=false)
    div = DOM.div(; style="width:100%;height:100%;")
    spec = map(plot) do value
        sprint(show, MIME"application/vnd.plotly.v1+json"(), value)
    end
    source = js"""
        Promise.all([$(PLOTLY_ASSET), $(PLOT_THEME), $(KomaMRIPlots.SLICE_PLOT_JS)]).then(([Plotly, theme, viewer]) => {
            const graph = $(div);
            if (!graph?.isConnected) return;
            const themeObserver = theme.observeTheme(graph, Plotly);
            const spec = $(spec);
            let pending = null;
            let alignPending = false;
            let rendering = false;
            let slices = null;
            let definition = null;
            let controls = null;
            const dispose = () => {
                if (graph.isConnected) return;
                observer.disconnect();
                themeObserver.disconnect();
                slices?.dispose();
                pending = null;
                alignPending = false;
                // An in-flight draw must finish before its graph state is purged.
                if (!rendering) Plotly.purge(graph);
            };
            const observer = new MutationObserver(dispose);
            observer.observe(document, {childList: true, subtree: true});
            const render = async () => {
                rendering = true;
                try {
                    // Keep only the latest selection while a draw is in progress.
                    while ((pending !== null || alignPending) && graph.isConnected) {
                        if (pending !== null) {
                            const figure = JSON.parse(pending);
                            pending = null;
                            delete figure.layout.width;
                            delete figure.layout.height;
                            figure.layout.autosize = true;
                            definition = await theme.renderPlot(graph, Plotly, figure, definition);
                            controls ??= theme.observeControls(graph, Plotly);
                            await controls();
                            const model=figure.layout.meta?.image_slices;
                            if (model) {
                                alignPending=false;
                                if (slices) await slices.rebind(model);
                                else slices=viewer.mountSlices(graph,Plotly);
                            } else { slices?.dispose(); slices=null; }
                            if ($(fit_colorbar) && !model && graph.isConnected) {
                                graph.removeListener('plotly_afterplot', alignColorbar);
                                graph.on('plotly_afterplot', alignColorbar);
                                alignPending = true;
                            }
                        }
                        if (alignPending && graph.isConnected) {
                            alignPending = false;
                            const [bottom, top] = graph._fullLayout.yaxis.domain;
                            const len = top - bottom;
                            const y = (top + bottom) / 2;
                            for (const [index, trace] of graph._fullData.entries()) {
                                const colorbar = trace.colorbar;
                                if (!trace.visible || !colorbar) continue;
                                if (colorbar.len !== len || colorbar.y !== y) {
                                    await Plotly.restyle(graph, {'colorbar.len': len, 'colorbar.y': y}, [index]);
                                }
                            }
                        }
                    }
                } finally {
                    rendering = false;
                    dispose();
                }
            };
            const requestRender = () => {
                if (!rendering) render().catch(error => Bonito.send_error("Plot update failed", error));
            };
            const alignColorbar = () => {
                if (!graph.isConnected || graph.layout.meta?.image_slices) return;
                alignPending = true;
                requestRender();
            };
            const update = value => {
                if (!graph.isConnected) return false;
                pending = value;
                requestRender();
            };
            spec.on(update);
            update(spec.value);
        }).catch(error => Bonito.send_error("Plot initialization failed", error));
    """
    return DOM.div(PLOTLY_ASSET, div, source; style="width:100%;height:100%;padding:8px;box-sizing:border-box;")
end
