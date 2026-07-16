const PLOTLY_ASSET = Asset(joinpath(artifact"plotly-artifacts", "plotly.min.js"); name="Plotly")

function plot_node(plot::PlotlyBase.Plot)
    div = DOM.div(; style="width:100%;height:100%;")
    spec = sprint(show, MIME"application/vnd.plotly.v1+json"(), plot)
    source = js"""
        const figure = JSON.parse($(spec));
        delete figure.layout.width;
        delete figure.layout.height;
        $(PLOTLY_ASSET).then(Plotly => Plotly.newPlot(
            $(div), figure.data, {...figure.layout, autosize: true},
            {...figure.config, responsive: true}
        ));
    """
    return DOM.div(PLOTLY_ASSET, div, source; style="width:100%;height:100%;padding:8px;box-sizing:border-box;")
end
