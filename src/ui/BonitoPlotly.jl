const PLOTLY_ASSET = Asset(PlotlyBase.cdn_url(); name="Plotly")

function plot_node(plot::PlotlyBase.Plot)
    div = DOM.div(; style="width:100%;height:100%;padding:8px;box-sizing:border-box;")
    spec = sprint(show, MIME"application/vnd.plotly.v1+json"(), plot)
    source = js"""
        const figure = JSON.parse($(spec));
        delete figure.layout.width;
        delete figure.layout.height;
        $(PLOTLY_ASSET).then(Plotly => Plotly.newPlot(
            $(div), figure.data, {...figure.layout, autosize: true},
            {...figure.config, responsive: true}
        ).then(() => {
            const observer = new ResizeObserver(() => {
                if (!$(div).isConnected) return observer.disconnect();
                Plotly.Plots.resize($(div));
            });
            observer.observe($(div).parentElement);
        }));
    """
    return DOM.div(PLOTLY_ASSET, div, source; style="width:100%;height:100%;")
end
