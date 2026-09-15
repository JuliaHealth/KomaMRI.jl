const PLOTLY_ASSET = Asset(joinpath(artifact"plotly-artifacts", "plotly.min.js"); name="Plotly")

plot_node(plot::PlotlyBase.Plot; fit_colorbar=false) =
    plot_node(Observable(plot); fit_colorbar)

function plot_node(plot::Observable; fit_colorbar=false)
    div = DOM.div(; style="width:100%;height:100%;")
    spec = map(plot) do value
        sprint(show, MIME"application/vnd.plotly.v1+json"(), value)
    end
    source = js"""
        $(PLOTLY_ASSET).then(Plotly => {
            const graph = $(div);
            if (!graph.isConnected) return;
            const spec = $(spec);
            let pending = null;
            let alignPending = false;
            let rendering = false;
            const dispose = () => {
                if (graph.isConnected) return;
                observer.disconnect();
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
                            await Plotly.react(
                                graph, figure.data, {...figure.layout, autosize: true},
                                {...figure.config, responsive: true}
                            );
                            if ($(fit_colorbar) && graph.isConnected) {
                                graph.removeListener('plotly_afterplot', alignColorbar);
                                graph.on('plotly_afterplot', alignColorbar);
                                alignPending = true;
                            }
                        }
                        if (alignPending && graph.isConnected) {
                            alignPending = false;
                            const colorbar = graph._fullData[0].colorbar;
                            if (!colorbar) continue;
                            const [bottom, top] = graph._fullLayout.yaxis.domain;
                            const len = top - bottom;
                            const y = (top + bottom) / 2;
                            if (colorbar.len !== len || colorbar.y !== y) {
                                await Plotly.restyle(graph, {'colorbar.len': len, 'colorbar.y': y}, [0]);
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
                if (!graph.isConnected) return;
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
        });
    """
    return DOM.div(PLOTLY_ASSET, div, source; style="width:100%;height:100%;padding:8px;box-sizing:border-box;")
end
