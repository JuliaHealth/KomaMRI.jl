// An embedding UI may supply a theme; standalone plots keep their own colors.
function themeColors(root, layout) {
    const host = root.closest('[data-plot-theme]');
    if (!host) return {};
    const [background, text, plot, grid] = JSON.parse(host.dataset.plotTheme);
    const colors = {paper_bgcolor: background, plot_bgcolor: plot, 'font.color': text,
        'modebar.bgcolor': background, 'modebar.color': text, 'modebar.activecolor': plot};
    const axisColors = path => {
        colors[`${path}.gridcolor`] = grid;
        colors[`${path}.zerolinecolor`] = grid;
        colors[`${path}.backgroundcolor`] = plot;
    };
    for (const key of Object.keys(layout)) {
        if (/^[xy]axis\d*$/.test(key)) axisColors(key);
        if (/^scene\d*$/.test(key))
            for (const axis of ['xaxis', 'yaxis', 'zaxis']) axisColors(`${key}.${axis}`);
    }
    return colors;
}

export function themeLayout(root, layout) {
    for (const [path, color] of Object.entries(themeColors(root, layout))) {
        const keys = path.split('.'), key = keys.pop();
        const target = keys.reduce((value, name) => value[name] ??= {}, layout);
        target[key] = color;
    }
    return layout;
}

export function themeTraces(root, data) {
    const host = root.closest('[data-plot-theme]');
    if (host) for (const trace of data)
        if (trace.meta?.theme === 'contrast') trace.marker.color = JSON.parse(host.dataset.plotTheme)[4];
    return data;
}

export function observeTheme(graph, Plotly) {
    const host = graph.closest('[data-plot-theme]');
    const observer = new MutationObserver(() => {
        if (!graph.isConnected || !graph._fullLayout) return;
        const update = themeColors(graph, graph.layout);
        // A dragged camera can be newer than the camera saved in layout.
        for (const key of Object.keys(graph.layout))
            if (/^scene\d*$/.test(key))
                update[`${key}.camera`] = graph._fullLayout[key]._scene.getCamera();
        const indices = graph.data.flatMap((trace, i) => trace.meta?.theme === 'contrast' ? [i] : []);
        if (indices.length) {
            const colors = themeTraces(graph, indices.map(i => graph.data[i])).map(trace => trace.marker.color);
            Plotly.update(graph, {'marker.color': colors}, update, indices);
        } else Plotly.relayout(graph, update);
    });
    if (host) observer.observe(host, {attributes: true, attributeFilter: ['data-plot-theme']});
    return observer;
}
