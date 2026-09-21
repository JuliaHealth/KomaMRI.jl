// The same component is embedded by Bonito in notebooks, browsers and KomaUI.
export async function mount(root, Plotly, request, response, overviewJSON) {
    const overview = JSON.parse(overviewJSON);
    const fullRange = overview.full_range;
    const graph = document.createElement('div');
    const clip = document.createElement('div');
    const navigator = document.createElement('div');
    root.style.display = 'flex';
    root.style.flexDirection = 'column';
    graph.style.cssText = 'flex:1 1 0;min-height:260px;width:100%';
    const group = root.closest('[data-linked-time-plots]');
    if (group) {
        root.style.minHeight = '0';
        graph.style.minHeight = '0';
    }
    clip.style.cssText = 'height:80px;overflow:hidden;flex:none';
    navigator.style.cssText = 'height:200px;position:relative';
    clip.append(navigator);
    root.append(graph, clip);

    let range, visibility, revision = 0, rendering = false, pending = null;
    let ready = false, disposed = false, navigatorReady = false;
    const first = JSON.parse(response.value);
    range = first.layout.xaxis?.range?.length === 2 ? [...first.layout.xaxis.range] : [...fullRange];
    visibility = first.data.map(trace => trace.visible ?? true);
    const hasSlider = first.layout.xaxis?.rangeslider?.visible ?? false;
    if (!hasSlider) clip.style.display = 'none';

    function fetchDetail() {
        if (disposed) return;
        request.notify({id: ++revision, range: [...range], width: graph.clientWidth, visibility});
    }
    function setRange(next) {
        if (next.every((t, i) => t === range[i])) return;
        range = [...next];
        for (const plot of navigatorReady ? [graph, navigator] : [graph]) {
            if (!range.every((t, i) => t === plot.layout.xaxis.range[i])) {
                const update = {'xaxis.range': [...range], 'xaxis.autorange': false};
                if (plot === navigator) update['xaxis.rangeslider.range'] = [...fullRange];
                Plotly.relayout(plot, update);
            }
        }
        fetchDetail();
        group?.dispatchEvent(new CustomEvent('koma-time-range', {detail: {source: root, range}}));
    }
    function annotation(payload) {
        if (payload.shown >= payload.available) return [];
        const fraction = 100 * payload.shown / payload.available;
        const percent = fraction >= 99.995 ? '&lt;100' : fraction < 0.01 ? '&lt;0.01' : fraction.toFixed(2);
        return [{name: 'sampling-detail', text: `Downsampled · ${percent}%`,
            xref: 'paper', yref: 'paper', x: 1, y: 0, xanchor: 'right', yanchor: 'bottom',
            xshift: -2, yshift: 2, showarrow: false, borderpad: 0, font: {size: 7.2},
            hovertext: `${payload.shown.toLocaleString()} / ${payload.available.toLocaleString()} samples shown. Zoom in for full detail.`}];
    }
    async function draw(payload) {
        const layout = payload.layout;
        layout.xaxis ??= {};
        layout.xaxis.domain = [0, 1];
        layout.xaxis.range = [...range];
        layout.xaxis.autorange = false;
        layout.xaxis.tickformat = '~f';
        layout.xaxis.nticks = 4;
        layout.xaxis.rangeslider = {visible: false};
        delete layout.xaxis.rangeselector;
        layout.height = graph.clientHeight;
        layout.width = graph.clientWidth;
        layout.uirevision = 'time-plot';
        for (const [key, axis] of Object.entries(graph._fullLayout ?? {})) {
            if (!/^yaxis\d*$/.test(key)) continue;
            layout[key] = {...layout[key], range: [...axis.range], autorange: false};
        }
        layout.annotations = [...(layout.annotations ?? []), ...annotation(payload)];
        // Plotly controls remain native; preserve their selection when detail arrives.
        for (const key of ['updatemenus', 'sliders']) {
            layout[key]?.forEach((control, i) => {
                control.active = graph.layout?.[key]?.[i]?.active ?? control.active;
            });
        }
        await Plotly.react(graph, payload.data, layout, {...payload.config, responsive: true});
    }
    async function render() {
        if (rendering || disposed) return;
        rendering = true;
        try {
            while (pending && !disposed) {
                const value = pending;
                pending = null;
                if (value.id >= revision) await draw(value);
            }
        } finally {
            rendering = false;
            if (disposed) { Plotly.purge(graph); Plotly.purge(navigator); }
        }
    }
    function update(serialized) {
        if (disposed) return false;
        pending = JSON.parse(serialized);
        render().catch(error => Bonito.send_error('Time plot update failed', error));
    }
    await draw({...first, data: overview.data.map((trace, i) => ({...trace, visible: visibility[i]}))});
    const fullYRanges = Object.fromEntries(Object.entries(graph._fullLayout)
        .filter(([key]) => /^yaxis\d*$/.test(key))
        .map(([key, axis]) => [key, [...axis.range]]));
    await draw(first);
    response.on(update);

    if (hasSlider) {
        const layout = first.layout;
        await Plotly.newPlot(navigator, overview.data.map((trace, i) => ({...trace,
            visible: visibility[i], showlegend: false, hoverinfo: 'skip', hovertemplate: null})), {
            height: 200, margin: {l: graph._fullLayout.margin.l, r: graph._fullLayout.margin.r, t: 0, b: 0},
            template: layout.template, colorway: layout.colorway,
            paper_bgcolor: layout.paper_bgcolor, plot_bgcolor: layout.plot_bgcolor, font: layout.font,
            showlegend: false, hovermode: false,
            xaxis: {anchor: layout.xaxis.anchor ?? 'y', range: [...range], autorange: false,
                showticklabels: false, ticks: '', showgrid: false, zeroline: false,
                rangeslider: {visible: true, autorange: false, range: [...fullRange], thickness: 0.4}},
            yaxis: {...layout.yaxis, fixedrange: true, visible: false},
            ...(layout.yaxis2 ? {yaxis2: {...layout.yaxis2, visible: false}} : {}),
            ...(layout.yaxis3 ? {yaxis3: {...layout.yaxis3, visible: false}} : {}),
        }, {displayModeBar: false, responsive: true});
        const fitSlider = () => {
            if (disposed) return;
            const slider = navigator.querySelector('.rangeslider-container').getBoundingClientRect();
            navigator.style.top = `-${slider.top - navigator.getBoundingClientRect().top}px`;
            clip.style.height = `${slider.height}px`;
        };
        navigator.on('plotly_afterplot', fitSlider);
        fitSlider();
        navigator.on('plotly_relayout', event => {
            if (!Object.keys(event).some(key => key.startsWith('xaxis.range'))) return;
            setRange(navigator.layout.xaxis.range);
        });
        navigatorReady = true;
    }
    graph.on('plotly_relayout', async event => {
        const reset = {};
        for (const [axis, limits] of Object.entries(fullYRanges)) {
            if (!event[`${axis}.autorange`]) continue;
            reset[`${axis}.range`] = [...limits];
            reset[`${axis}.autorange`] = false;
        }
        if (Object.keys(reset).length) await Plotly.relayout(graph, reset);
        let next;
        if (event['xaxis.autorange']) next = fullRange;
        else if (event['xaxis.range']) next = event['xaxis.range'];
        else if (event['xaxis.range[0]'] !== undefined || event['xaxis.range[1]'] !== undefined)
            next = [event['xaxis.range[0]'] ?? range[0], event['xaxis.range[1]'] ?? range[1]];
        if (!next || next.every((t, i) => t === range[i])) return;
        setRange(next);
    });
    graph.on('plotly_restyle', () => {
        const next = graph.data.map(trace => trace.visible ?? true);
        if (next.every((v, i) => v === visibility[i])) return;
        visibility = next;
        if (navigatorReady) Plotly.restyle(navigator, {visible: visibility});
        fetchDetail();
    });
    const syncRange = event => {
        if (event.detail.source !== root) setRange(event.detail.range);
    };
    group?.addEventListener('koma-time-range', syncRange);
    let size = [root.clientWidth, root.clientHeight];
    const resize = new ResizeObserver(() => {
        const next = [root.clientWidth, root.clientHeight];
        if (!ready || next.every((value, i) => Math.abs(value - size[i]) < 2)) return;
        size = next;
        if (navigatorReady) Plotly.Plots.resize(navigator);
        fetchDetail();
    });
    const removal = new MutationObserver(() => {
        if (root.isConnected) return;
        disposed = true; pending = null; resize.disconnect(); removal.disconnect();
        group?.removeEventListener('koma-time-range', syncRange);
        if (!rendering) { Plotly.purge(graph); Plotly.purge(navigator); }
    });
    removal.observe(document, {subtree: true, childList: true});
    resize.observe(root);
    ready = true;
    fetchDetail();
}
