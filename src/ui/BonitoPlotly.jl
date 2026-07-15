const PLOTLY_ASSET = Asset(joinpath(artifact"plotly-artifacts", "plotly.min.js"); name="Plotly")

function plot_node(plot::PlotlyBase.Plot)
    div = DOM.div(; style="width:100%;height:100%;")
    spec = sprint(show, MIME"application/vnd.plotly.v1+json"(), plot)
    source = js"""
        const figure = JSON.parse($(spec));
        delete figure.layout.width;
        delete figure.layout.height;
        const horizontalLegend = figure.layout.legend?.orientation === 'h';
        const topMargin = figure.layout.margin?.t || 0;
        if (horizontalLegend) {
            figure.layout.margin = {...figure.layout.margin, autoexpand: false};
        }
        $(PLOTLY_ASSET).then(Plotly => Plotly.newPlot(
            $(div), figure.data, {...figure.layout, autosize: true},
            {...figure.config, responsive: true}
        ).then(() => {
            let resize = 0;
            const observer = new ResizeObserver(async () => {
                const current = ++resize;
                if (!$(div).isConnected) return observer.disconnect();
                await Plotly.Plots.resize($(div));
                await new Promise(requestAnimationFrame);
                if (current !== resize) return;
                const legend = $(div).querySelector('.legend');
                if (!legend || !horizontalLegend) return;
                const margin = {
                    ...$(div).layout.margin,
                    t: Math.max(topMargin, Math.ceil(legend.getBBox().height)),
                };
                await Plotly.relayout($(div), {margin});
            });
            observer.observe($(div).parentElement);
        }));
    """
    return DOM.div(PLOTLY_ASSET, div, source; style="width:100%;height:100%;padding:8px;box-sizing:border-box;")
end
