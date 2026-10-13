// deno-fmt-ignore-file
// deno-lint-ignore-file
// This code was bundled using `deno bundle` and it's not recommended to edit it manually

const same = (a, b)=>JSON.stringify(a) === JSON.stringify(b);
const choices = (control)=>control.steps ?? control.buttons ?? [];
const controlKey = (control)=>JSON.stringify([
        control.name,
        control.currentvalue?.prefix,
        choices(control)[0]?.label
    ]);
function reconcileControls(layout, current, previous) {
    for (const key of [
        'sliders',
        'updatemenus'
    ]){
        const controls = (layout[key] ?? []).map((control)=>[
                control,
                previous?.[key]?.find((item)=>controlKey(item) === controlKey(control)),
                current?.[key]?.find((item)=>controlKey(item) === controlKey(control))
            ]);
        const overrides = new Set(controls.filter(([control, old])=>old && !same(control.active ?? 0, old.active ?? 0)).map(([control])=>control.name ?? controlKey(control)));
        for (const [control, old, live] of controls){
            if (!old || !live || overrides.has(control.name ?? controlKey(control))) continue;
            const selected = choices(live)[live.active ?? 0]?.label;
            const index = live.active === -1 ? -1 : choices(control).findIndex((choice)=>same(choice.label, selected));
            if (live.active === -1 || index >= 0) control.active = index;
        }
    }
}
async function applyControlSelections(graph, Plotly) {
    const hidden = new Set(graph.data.filter((trace)=>trace.visible === 'legendonly').map((trace)=>trace.uid));
    const selections = [
        'sliders',
        'updatemenus'
    ].flatMap((key)=>(graph.layout[key] ?? []).flatMap((control)=>{
            const selection = choices(control)[control.active ?? 0];
            return selection && [
                'restyle',
                'relayout',
                'update'
            ].includes(selection.method) ? [
                selection
            ] : [];
        }));
    for (const selection of selections)await Plotly[selection.method](graph, ...selection.args);
    if (selections.length && hidden.size) await Plotly.restyle(graph, {
        visible: graph.data.map((trace)=>hidden.has(trace.uid) ? 'legendonly' : trace.visible ?? true)
    });
}
function identifyTraces(data) {
    const counts = new Map();
    for (const trace of data){
        const key = JSON.stringify([
            trace.name,
            trace.type,
            trace.mode,
            trace.xaxis,
            trace.yaxis,
            trace.legendgroup
        ]);
        const count = (counts.get(key) ?? 0) + 1;
        counts.set(key, count);
        trace.uid ??= `koma-${Array.from(key, (character)=>character.codePointAt(0).toString(16)).join('-')}-${count}`;
    }
    return data;
}
function reconcilePlot(graph, data, layout, previous) {
    identifyTraces(data);
    const definition = {
        layout: structuredClone(layout),
        traces: new Map(data.map((trace)=>[
                trace.uid,
                trace.visible ?? true
            ])),
        revision: (previous?.revision ?? 0) + 1
    };
    if (previous && graph._fullLayout) {
        const oldTraces = new Map(graph.data.map((trace)=>[
                trace.uid,
                trace
            ]));
        for (const trace of data){
            const live = oldTraces.get(trace.uid);
            if (live?.visible !== false && previous.traces.has(trace.uid) && same(trace.visible ?? true, previous.traces.get(trace.uid))) trace.visible = live?.visible ?? trace.visible;
        }
        for (const [key, axis] of Object.entries(layout)){
            if (!/^[xy]axis\d*$/.test(key)) continue;
            const old = previous.layout[key], current = graph._fullLayout[key];
            if (old && current && !axis.fixedrange && same(axis.range, old.range) && same(axis.autorange, old.autorange) && same(axis.title, old.title)) Object.assign(axis, {
                range: [
                    ...current.range
                ],
                autorange: false
            });
        }
        for (const [key, scene] of Object.entries(layout)){
            if (!/^scene\d*$/.test(key) || !previous.layout[key] || !graph._fullLayout[key]) continue;
            if (same(scene.camera, previous.layout[key].camera) && [
                'xaxis',
                'yaxis',
                'zaxis'
            ].every((axis)=>same(scene[axis]?.title, previous.layout[key][axis]?.title))) scene.camera = graph._fullLayout[key]._scene.getCamera();
        }
        reconcileControls(layout, graph.layout, previous.layout);
    }
    layout.uirevision = `${layout.uirevision ?? 'koma'}-${definition.revision}`;
    return definition;
}
function themeColors(root, layout) {
    const host = root.closest('[data-plot-theme]');
    if (!host) return {};
    const [background, text, plot, grid] = JSON.parse(host.dataset.plotTheme);
    const colors = {
        paper_bgcolor: background,
        plot_bgcolor: plot,
        'font.color': text,
        'modebar.bgcolor': background,
        'modebar.color': text,
        'modebar.activecolor': plot
    };
    const axisColors = (path)=>{
        colors[`${path}.gridcolor`] = grid;
        colors[`${path}.zerolinecolor`] = grid;
        colors[`${path}.backgroundcolor`] = plot;
    };
    for (const key of Object.keys(layout)){
        if (/^[xy]axis\d*$/.test(key)) axisColors(key);
        if (/^scene\d*$/.test(key)) for (const axis of [
            'xaxis',
            'yaxis',
            'zaxis'
        ])axisColors(`${key}.${axis}`);
    }
    return colors;
}
function themeLayout(root, layout) {
    for (const [path, color] of Object.entries(themeColors(root, layout))){
        const keys = path.split('.'), key = keys.pop();
        const target = keys.reduce((value, name)=>value[name] ??= {}, layout);
        target[key] = color;
    }
    return layout;
}
function themeTraces(root, data) {
    const host = root.closest('[data-plot-theme]');
    if (host) {
        for (const trace of data)if (trace.meta?.theme === 'contrast') trace.marker.color = JSON.parse(host.dataset.plotTheme)[4];
    }
    return data;
}
function observeTheme(graph, Plotly) {
    const host = graph.closest('[data-plot-theme]');
    const observer = new MutationObserver(()=>{
        if (!graph.isConnected || !graph._fullLayout) return;
        const update = themeColors(graph, graph.layout);
        for (const key of Object.keys(graph.layout))if (/^scene\d*$/.test(key)) update[`${key}.camera`] = graph._fullLayout[key]._scene.getCamera();
        const indices = graph.data.flatMap((trace, i)=>trace.meta?.theme === 'contrast' ? [
                i
            ] : []);
        if (indices.length) {
            const colors = themeTraces(graph, indices.map((i)=>graph.data[i])).map((trace)=>trace.marker.color);
            Plotly.update(graph, {
                'marker.color': colors
            }, update, indices);
        } else Plotly.relayout(graph, update);
    });
    if (host) observer.observe(host, {
        attributes: true,
        attributeFilter: [
            'data-plot-theme'
        ]
    });
    return observer;
}
async function mount(root, Plotly, request, response, overviewJSON) {
    let overview = JSON.parse(overviewJSON);
    let fullRange = overview.full_range;
    const graph = document.createElement('div');
    const clip = document.createElement('div');
    const navigator = document.createElement('div');
    root.style.display = 'flex';
    root.style.flexDirection = 'column';
    graph.style.cssText = 'flex:1 1 0;min-height:0;width:100%';
    const group = root.closest('[data-linked-time-plots]');
    if (group) {
        root.style.minHeight = '0';
    }
    clip.style.cssText = 'height:80px;overflow:hidden;flex:none';
    navigator.style.cssText = 'height:200px;position:relative';
    clip.append(navigator);
    root.append(graph, clip);
    const theme = observeTheme(graph, Plotly);
    const navigatorTheme = observeTheme(navigator, Plotly);
    let range, visibility, revision = 0, rendering = false, pending = null;
    let ready = false, disposed = false, navigatorReady = false, renderTask;
    const first = JSON.parse(response.value);
    let definition = reconcilePlot(graph, first.data, first.layout);
    range = first.layout.xaxis?.range?.length === 2 ? [
        ...first.layout.xaxis.range
    ] : [
        ...fullRange
    ];
    visibility = first.data.map((trace)=>trace.visible ?? true);
    let hasSlider = first.layout.xaxis?.rangeslider?.visible ?? false, fullYRanges;
    function fetchDetail() {
        if (disposed) return;
        request.notify({
            id: ++revision,
            range: [
                ...range
            ],
            width: graph.clientWidth,
            visibility
        });
    }
    function setRange(next) {
        if (next.every((t, i)=>t === range[i])) return;
        range = [
            ...next
        ];
        for (const plot of hasSlider && navigatorReady ? [
            graph,
            navigator
        ] : [
            graph
        ]){
            if (!range.every((t, i)=>t === plot.layout.xaxis.range[i])) {
                const update = {
                    'xaxis.range': [
                        ...range
                    ],
                    'xaxis.autorange': false
                };
                if (plot === navigator) update['xaxis.rangeslider.range'] = [
                    ...fullRange
                ];
                Plotly.relayout(plot, update);
            }
        }
        fetchDetail();
        group?.dispatchEvent(new CustomEvent('koma-time-range', {
            detail: {
                source: root,
                range
            }
        }));
    }
    function annotation(payload) {
        if (payload.shown >= payload.available) return [];
        const fraction = 100 * payload.shown / payload.available;
        const percent = fraction >= 99.995 ? '&lt;100' : fraction < 0.01 ? '&lt;0.01' : fraction.toFixed(2);
        return [
            {
                name: 'sampling-detail',
                text: `Downsampled · ${percent}%`,
                xref: 'paper',
                yref: 'paper',
                x: 1,
                y: 0,
                xanchor: 'right',
                yanchor: 'bottom',
                xshift: -2,
                yshift: 2,
                showarrow: false,
                borderpad: 0,
                font: {
                    size: 7.2
                },
                hovertext: `${payload.shown.toLocaleString()} / ${payload.available.toLocaleString()} samples shown. Zoom in for full detail.`
            }
        ];
    }
    async function draw(payload, preserveY = true) {
        const layout = themeLayout(root, structuredClone(payload.layout));
        layout.xaxis ??= {};
        layout.xaxis.domain = [
            0,
            1
        ];
        layout.xaxis.range = [
            ...range
        ];
        layout.xaxis.autorange = false;
        layout.xaxis.tickformat = '~f';
        layout.xaxis.nticks = 4;
        layout.xaxis.rangeslider = {
            visible: false
        };
        delete layout.xaxis.rangeselector;
        layout.height = graph.clientHeight;
        layout.width = graph.clientWidth;
        if (graph.layout && (graph.layout.width !== layout.width || graph.layout.height !== layout.height)) {
            await Plotly.relayout(graph, {
                width: layout.width,
                height: layout.height
            });
        }
        layout.uirevision = `time-plot-${definition.revision}`;
        for (const [key, axis] of Object.entries(preserveY ? graph._fullLayout ?? {} : {})){
            if (!/^yaxis\d*$/.test(key)) continue;
            layout[key] = {
                ...layout[key],
                range: [
                    ...axis.range
                ],
                autorange: false
            };
        }
        layout.annotations = [
            ...layout.annotations ?? [],
            ...annotation(payload)
        ];
        if (preserveY) reconcileControls(layout, graph.layout, definition.layout);
        await Plotly.react(graph, themeTraces(root, identifyTraces(payload.data)), layout, {
            ...payload.config,
            responsive: true
        });
    }
    async function render() {
        if (rendering || disposed) return;
        rendering = true;
        try {
            while(pending && !disposed){
                const value = pending;
                pending = null;
                if (value.id >= revision) await draw(value);
            }
        } finally{
            rendering = false;
            if (disposed) {
                Plotly.purge(graph);
                Plotly.purge(navigator);
            }
        }
    }
    function update(serialized) {
        if (disposed) return false;
        pending = JSON.parse(serialized);
        if (!rendering) renderTask = render().catch((error)=>console.error('Time plot update failed', error));
    }
    await draw({
        ...first,
        data: overview.data.map((trace, i)=>({
                ...trace,
                visible: visibility[i]
            }))
    });
    const yRanges = ()=>Object.fromEntries(Object.entries(graph._fullLayout).filter(([key])=>/^yaxis\d*$/.test(key)).map(([key, axis])=>[
                key,
                [
                    ...axis.range
                ]
            ]));
    fullYRanges = yRanges();
    await draw(first);
    response.on(update);
    async function updateNavigator(layout) {
        clip.style.display = hasSlider ? '' : 'none';
        if (!hasSlider) return;
        await Plotly.react(navigator, themeTraces(root, identifyTraces(overview.data).map((trace, i)=>({
                ...trace,
                visible: visibility[i],
                showlegend: false,
                hoverinfo: 'skip',
                hovertemplate: null
            }))), {
            height: 200,
            margin: {
                l: graph._fullLayout.margin.l,
                r: graph._fullLayout.margin.r,
                t: 0,
                b: 0
            },
            template: layout.template,
            colorway: layout.colorway,
            paper_bgcolor: layout.paper_bgcolor,
            plot_bgcolor: layout.plot_bgcolor,
            font: layout.font,
            showlegend: false,
            hovermode: false,
            xaxis: {
                anchor: layout.xaxis.anchor ?? 'y',
                range: [
                    ...range
                ],
                autorange: false,
                showticklabels: false,
                ticks: '',
                showgrid: false,
                zeroline: false,
                rangeslider: {
                    visible: true,
                    autorange: false,
                    range: [
                        ...fullRange
                    ],
                    thickness: 0.4
                }
            },
            yaxis: {
                ...layout.yaxis,
                fixedrange: true,
                visible: false
            },
            ...layout.yaxis2 ? {
                yaxis2: {
                    ...layout.yaxis2,
                    visible: false
                }
            } : {},
            ...layout.yaxis3 ? {
                yaxis3: {
                    ...layout.yaxis3,
                    visible: false
                }
            } : {}
        }, {
            displayModeBar: false,
            responsive: true
        });
        if (navigatorReady) return;
        const fitSlider = ()=>{
            if (disposed) return;
            const slider = navigator.querySelector('.rangeslider-container').getBoundingClientRect();
            navigator.style.top = `-${slider.top - navigator.getBoundingClientRect().top}px`;
            clip.style.height = `${slider.height}px`;
        };
        navigator.on('plotly_afterplot', fitSlider);
        fitSlider();
        navigator.on('plotly_relayout', (event)=>{
            if (rendering || !Object.keys(event).some((key)=>key.startsWith('xaxis.range'))) return;
            setRange(navigator.layout.xaxis.range);
        });
        navigatorReady = true;
    }
    await updateNavigator(first.layout);
    graph.on('plotly_relayout', async (event)=>{
        if (rendering) return;
        const reset = {};
        for (const [axis, limits] of Object.entries(fullYRanges)){
            if (!event[`${axis}.autorange`]) continue;
            reset[`${axis}.range`] = [
                ...limits
            ];
            reset[`${axis}.autorange`] = false;
        }
        if (Object.keys(reset).length) await Plotly.relayout(graph, reset);
        let next;
        if (event['xaxis.autorange']) next = fullRange;
        else if (event['xaxis.range']) next = event['xaxis.range'];
        else if (event['xaxis.range[0]'] !== undefined || event['xaxis.range[1]'] !== undefined) next = [
            event['xaxis.range[0]'] ?? range[0],
            event['xaxis.range[1]'] ?? range[1]
        ];
        if (!next || next.every((t, i)=>t === range[i])) return;
        setRange(next);
    });
    graph.on('plotly_restyle', ()=>{
        if (rendering) return;
        const next = graph.data.map((trace)=>trace.visible ?? true);
        if (next.length === visibility.length && next.every((v, i)=>v === visibility[i])) return;
        visibility = next;
        if (hasSlider && navigatorReady) Plotly.restyle(navigator, {
            visible: visibility
        });
        fetchDetail();
    });
    const syncRange = (event)=>{
        if (event.detail.source !== root) setRange(event.detail.range);
    };
    group?.addEventListener('koma-time-range', syncRange);
    let size = [
        root.clientWidth,
        root.clientHeight
    ];
    const resize = new ResizeObserver(()=>{
        const next = [
            root.clientWidth,
            root.clientHeight
        ];
        if (!ready || next.every((value, i)=>Math.abs(value - size[i]) < 2)) return;
        size = next;
        if (hasSlider && navigatorReady) Plotly.Plots.resize(navigator);
        fetchDetail();
    });
    const removal = new MutationObserver(()=>{
        if (root.isConnected) return;
        disposed = true;
        pending = null;
        resize.disconnect();
        removal.disconnect();
        theme.disconnect();
        navigatorTheme.disconnect();
        group?.removeEventListener('koma-time-range', syncRange);
        if (!rendering) {
            Plotly.purge(graph);
            Plotly.purge(navigator);
        }
    });
    removal.observe(document, {
        subtree: true,
        childList: true
    });
    resize.observe(root);
    ready = true;
    fetchDetail();
    return {
        graph,
        async rebind (initialJSON, overviewJSON) {
            revision++;
            pending = null;
            await renderTask;
            const next = JSON.parse(initialJSON);
            const originalLayout = structuredClone(next.layout);
            const oldYRanges = yRanges();
            definition = reconcilePlot(graph, next.data, next.layout, definition);
            overview = JSON.parse(overviewJSON);
            fullRange = overview.full_range;
            range = next.layout.xaxis?.range?.length === 2 ? [
                ...next.layout.xaxis.range
            ] : [
                ...fullRange
            ];
            if (range[1] < fullRange[0] || range[0] > fullRange[1]) range = [
                ...fullRange
            ];
            visibility = next.data.map((trace)=>trace.visible ?? true);
            hasSlider = originalLayout.xaxis?.rangeslider?.visible ?? false;
            ready = false;
            rendering = true;
            try {
                await draw({
                    ...next,
                    layout: originalLayout,
                    data: overview.data.map((trace, i)=>({
                            ...trace,
                            visible: visibility[i]
                        }))
                }, false);
                const nextYRanges = yRanges();
                for (const [key, limits] of Object.entries(nextYRanges)){
                    if (JSON.stringify(oldYRanges[key]) === JSON.stringify(fullYRanges[key])) Object.assign(next.layout[key] ??= {}, {
                        range: limits,
                        autorange: false
                    });
                }
                fullYRanges = nextYRanges;
                await draw(next, false);
                await applyControlSelections(graph, Plotly);
                visibility = graph.data.map((trace)=>trace.visible ?? true);
                await updateNavigator(next.layout);
            } finally{
                ready = true;
                rendering = false;
            }
            fetchDetail();
            group?.dispatchEvent(new CustomEvent('koma-time-range', {
                detail: {
                    source: root,
                    range
                }
            }));
        }
    };
}
export { mount as mount };
