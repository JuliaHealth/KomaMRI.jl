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
function styleTimeTicks(graph) {
    graph.querySelectorAll('.slider-group').forEach((group)=>{
        const slider = group.__data__;
        if (slider.name !== 'time') return;
        const labels = [
            ...group.querySelectorAll('.slider-label-group')
        ];
        if (!labels.length) return;
        const ticks = [
            ...group.querySelectorAll('.slider-tick-rect')
        ];
        if (labels.length === 1) labels.push(labels[0].parentNode.appendChild(labels[0].cloneNode(true)));
        for (const tick of ticks)if (Number(tick.getAttribute('height')) === slider.ticklen) tick.setAttribute('height', slider.minorticklen);
        labels.forEach((label, i)=>{
            const index = Math.round(i * (ticks.length - 1) / (labels.length - 1));
            const tick = ticks[index], value = slider._visibleSteps[index].label;
            const x = tick.transform.baseVal.consolidate().matrix.e + slider.tickwidth / 2;
            const y = label.transform.baseVal.consolidate().matrix.f;
            label.setAttribute('transform', `translate(${x},${y})`);
            const text = label.querySelector('.slider-label');
            text.textContent = value;
            text.setAttribute('data-unformatted', value);
            text.setAttribute('text-anchor', i === 0 ? 'start' : i === labels.length - 1 ? 'end' : 'middle');
            tick.setAttribute('height', slider.ticklen);
        });
    });
}
async function styleControls(graph, Plotly, selected = (menu, index)=>index === menu.active) {
    graph.querySelector('.menulayer')?.classList.toggle('koma-slice-controls', Boolean(graph.layout.meta?.phantom_slices));
    graph.querySelectorAll('.slider-group').forEach((group)=>group.classList.toggle('koma-slice-position', Boolean(group.__data__.name?.startsWith('slice-'))));
    if (!graph.querySelector('.koma-control-style')) {
        const style = document.createElement('style');
        style.className = 'koma-control-style';
        style.textContent = `
.slider-grip-rect:focus { outline:none; }
.slider-grip-rect:focus-visible { stroke:#2a7fb8 !important; stroke-width:2px !important; }
.koma-slice-position .slider-labels { display:none; }
.koma-component-button .updatemenu-item-text,
.koma-slice-controls :is(.updatemenu-header, .updatemenu-dropdown-button) .updatemenu-item-text { fill:white !important; }
.koma-component-button .updatemenu-item-rect,
.koma-slice-controls :is(.updatemenu-header, .updatemenu-dropdown-button) .updatemenu-item-rect { fill:#2a7fb8 !important; fill-opacity:1 !important; }
.koma-component-button[aria-pressed="true"] .updatemenu-item-rect,
.koma-component-button:hover .updatemenu-item-rect,
.koma-slice-controls :is(.updatemenu-header, .updatemenu-dropdown-button):hover .updatemenu-item-rect { fill:#1e5c85 !important; stroke:#1e5c85 !important; }`;
        graph.append(style);
        graph.on('plotly_afterplot', ()=>styleTimeTicks(graph));
    }
    const menus = [
        ...graph.querySelectorAll('.updatemenu-header-group')
    ];
    menus.forEach((group)=>{
        const menu = group.__data__;
        if (![
            'phantom-property',
            'coil-component',
            'image-component',
            'motion-playback'
        ].includes(menu.name)) return;
        group.querySelectorAll('.updatemenu-button').forEach((button, index)=>{
            button.classList.add('koma-component-button');
            button.setAttribute('role', 'button');
            button.setAttribute('aria-pressed', String(selected(menu, index)));
        });
    });
    const groups = menus.filter((group)=>group.__data__.name === 'phantom-property');
    if (!groups.length) return styleTimeTicks(graph);
    let bounds = groups.map((group)=>group.getBoundingClientRect());
    const plotBounds = graph.getBoundingClientRect();
    const modebar = graph.querySelector('.modebar')?.getBoundingClientRect();
    const thickness = menus.find((group)=>group.__data__.name === 'slab-thickness');
    if (thickness) {
        const header = thickness.querySelector('.updatemenu-header');
        const text = header.querySelector('.updatemenu-item-text');
        const arrow = thickness.querySelector('.updatemenu-header-arrow');
        const inset = 6, width = text.getBBox().width + arrow.getBBox().width + 3 * inset;
        header.querySelector('rect').style.width = `${width}px`;
        text.style.transform = `translateX(${inset - Number(text.getAttribute('x'))}px)`;
        arrow.style.transform = `translateX(${width - inset - Number(arrow.getAttribute('x'))}px)`;
    }
    const thicknessBounds = thickness?.getBoundingClientRect();
    const reserve = thickness ? thicknessBounds.width + 12 : 0;
    const controlWidth = bounds.reduce((width, bound, i)=>width + bound.width + (i ? groups[i - 1].__data__.pad.r + groups[i].__data__.pad.l : 0), reserve);
    const right = modebar?.width && (plotBounds.left + plotBounds.right + controlWidth) / 2 > modebar.left ? modebar.left - 12 : plotBounds.right;
    const viewport = {
        left: plotBounds.left,
        right,
        width: right - plotBounds.left
    };
    const legend = graph.querySelector('.legend');
    async function shiftControls(controls, delta) {
        if (!delta) return;
        const update = {
            'margin.t': graph.layout.margin.t + delta
        };
        for (const group of controls)update[`updatemenus[${group.__data__._index}].pad.b`] = group.__data__.pad.b + delta;
        await Plotly.relayout(graph, update);
    }
    if (thickness) {
        await shiftControls([
            ...groups,
            thickness
        ], 24 - thickness.__data__.pad.b);
    }
    let stacked = false;
    if (groups.length === 2) {
        const [intrinsicIndex, velocityIndex] = groups.map((group)=>group.__data__._index);
        const [intrinsic, velocity] = [
            intrinsicIndex,
            velocityIndex
        ].map((index)=>graph.layout.updatemenus[index]);
        const gap = intrinsic.pad.r + velocity.pad.l;
        stacked = bounds[0].width + bounds[1].width + gap + reserve > viewport.width;
        const padding = velocity.pad.b + (stacked ? Math.max(...bounds.map((bound)=>bound.height)) + gap / 2 : 0);
        if (intrinsic.pad.b !== padding) {
            await Plotly.relayout(graph, {
                'margin.t': graph.layout.margin.t + padding - intrinsic.pad.b,
                [`updatemenus[${intrinsicIndex}].pad.b`]: padding,
                [`updatemenus[${intrinsicIndex}].xanchor`]: stacked ? 'center' : 'right',
                [`updatemenus[${velocityIndex}].xanchor`]: stacked ? 'center' : 'left'
            });
        }
    }
    groups.forEach((group)=>group.style.transform = '');
    bounds = groups.map((group)=>group.getBoundingClientRect());
    const rowWidth = stacked ? bounds.at(-1).width : bounds.at(-1).right - bounds[0].left;
    const alongside = rowWidth + reserve <= viewport.width;
    if (thickness) {
        const index = thickness.__data__._index;
        const dropdown = graph.layout.updatemenus[index];
        const separation = Math.min(...groups.map((group)=>group.__data__.pad.b)) - dropdown.pad.b;
        const delta = (alongside ? 0 : thicknessBounds.height + 12) - separation;
        await shiftControls(groups, delta);
        const bound = thickness.getBoundingClientRect();
        const offset = (viewport.left + viewport.right - bound.left - bound.right) / 2 + (alongside ? (rowWidth + 12) / 2 : 0);
        if (Math.round(offset)) await Plotly.relayout(graph, {
            [`updatemenus[${index}].x`]: dropdown.x + offset / graph._fullLayout._size.w
        });
    }
    groups.forEach((group)=>group.style.transform = '');
    bounds = groups.map((group)=>group.getBoundingClientRect());
    const offset = (viewport.left + viewport.right - bounds[0].left - bounds.at(-1).right - (alongside ? reserve : 0)) / 2;
    groups.forEach((group, i)=>group.style.transform = `translateX(${stacked ? (viewport.left + viewport.right - bounds[i].left - bounds[i].right - (alongside && i === groups.length - 1 ? reserve : 0)) / 2 : offset}px)`);
    bounds = groups.map((group)=>group.getBoundingClientRect());
    graph.querySelectorAll('.annotation').forEach((annotation)=>{
        if (annotation.textContent !== '—') return;
        annotation.style.transform = '';
        const bound = annotation.getBoundingClientRect();
        annotation.style.transform = `translate(${(bounds[0].right + bounds.at(-1).left - bound.left - bound.right) / 2}px,` + `${(bounds[0].top + bounds[0].bottom - bound.top - bound.bottom) / 2}px)`;
        annotation.style.visibility = stacked ? 'hidden' : '';
    });
    if (legend) {
        legend.style.translate = '';
        legend.style.translate = `0px ${Math.max(...bounds.map((bound)=>bound.bottom), thickness?.getBoundingClientRect().bottom ?? 0) - legend.getBoundingClientRect().top + groups[0].__data__.pad.r}px`;
    }
    styleTimeTicks(graph);
}
async function resetSpatialAxes(graph, Plotly, event) {
    const range = graph.layout.meta?.spatial_range;
    if (!range) return;
    const reset = {};
    for (const axis of Object.keys(graph.layout).filter((key)=>/^[xy]axis\d*$/.test(key)))if (event[`${axis}.autorange`]) {
        const reversed = graph.layout.meta?.phantom_slices?.planes.some((plane)=>plane.axes[plane.dimensions.indexOf(1)] === axis);
        reset[`${axis}.range`] = reversed ? [
            ...range
        ].reverse() : [
            ...range
        ];
        reset[`${axis}.autorange`] = false;
    }
    if (Object.keys(reset).length) await Plotly.relayout(graph, reset);
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
const FILTER_SETTLE_MS = 150;
function mountKspaceRange(root, graph, config, onChange) {
    root.style.display = 'flex';
    root.style.flexDirection = 'column';
    graph.style.cssText = 'width:100%;flex:1;min-height:0';
    const navigator = root.querySelector('.koma-time-plot');
    let pending, timer, dragging = false;
    const control = {
        query: {},
        dispose,
        rebind (next) {
            clearTimeout(timer);
            pending = null;
            dragging = false;
            config = next;
            const interval = navigator.querySelector('.js-plotly-plot')?.layout?.xaxis?.range;
            if (interval) control.query = {
                time_range: interval.map((t)=>Math.max(0, Math.min(config.time_max, t)))
            };
        }
    };
    control.rebind(config);
    function apply() {
        if (dragging || !pending) return;
        control.query = {
            time_range: pending
        };
        pending = null;
        onChange();
    }
    function schedule() {
        clearTimeout(timer);
        if (!dragging) timer = setTimeout(apply, FILTER_SETTLE_MS);
    }
    function change(event) {
        if (event.detail.source !== navigator) return;
        pending = event.detail.range.map((t)=>Math.max(0, Math.min(config.time_max, t)));
        schedule();
    }
    function begin(event) {
        if (!event.target.closest('.rangeslider-container,.draglayer')) return;
        dragging = true;
        clearTimeout(timer);
    }
    function finish() {
        if (!dragging) return;
        dragging = false;
        if (pending) schedule();
    }
    function dispose() {
        clearTimeout(timer);
        root.removeEventListener('koma-time-range', change);
        navigator.removeEventListener('pointerdown', begin);
        document.removeEventListener('pointerup', finish);
        document.removeEventListener('pointercancel', finish);
    }
    root.addEventListener('koma-time-range', change);
    navigator.addEventListener('pointerdown', begin);
    document.addEventListener('pointerup', finish);
    document.addEventListener('pointercancel', finish);
    return control;
}
function sliceAnnotations(model, counts, time) {
    return model.planes.flatMap((plane, panel)=>{
        const { shown , available  } = counts[panel];
        if (shown >= available) return [];
        return {
            name: 'slice-sampling',
            text: `Downsampled · ${(100 * shown / available).toFixed(2)}%`,
            xref: `${plane.refs[0]} domain`,
            yref: `${plane.refs[1]} domain`,
            x: 1,
            y: 0,
            xanchor: 'right',
            yanchor: 'bottom',
            xshift: -4,
            yshift: 4,
            showarrow: false,
            borderpad: 0,
            font: {
                size: 10
            },
            hovertext: `${shown.toLocaleString()} / ${available.toLocaleString()} spins in this slab's viewport ` + `at t = ${time} ms. Zoom for more detail.`
        };
    });
}
function mountSlices(graph, Plotly, onChange) {
    let model = graph.layout.meta?.phantom_slices ?? graph.layout.meta?.image_slices;
    if (!model) return null;
    let image = Boolean(model.shape);
    let center = [
        ...model.center
    ], thickness = model.thickness;
    let frame = 0, property = graph.layout.meta.phantom_property;
    let disposed = false, rendering = false, pending = false, updateFrame, playTimer, clickTimer, playing = false, playStarted;
    const positionSliders = ()=>graph.layout.sliders.filter((slider)=>slider.name?.startsWith('slice-'));
    const clamp = (value, axis)=>Math.max(model.bounds[axis][0], Math.min(model.bounds[axis][1], image ? Math.round(value) : value));
    const nearest = (steps, value)=>steps.reduce((best, step, i)=>Math.abs(Number(step.value) - value) < Math.abs(Number(steps[best].value) - value) ? i : best, 0);
    function applyLayout(layout) {
        for (const plane of model.planes){
            const position = layout.annotations?.find((item)=>item.name === `slice-position-${plane.normal}`);
            if (position) position.text = image ? `${'xyz'[plane.normal - 1]} = ${plane.normal === 3 ? model.zlabels[center[2]] : center[plane.normal - 1]}` : `${'xyz'[plane.normal - 1]} = ${Math.round(1000 * center[plane.normal - 1]) / 100} mm`;
        }
        layout.shapes = model.planes.flatMap((plane)=>{
            const [x, y] = plane.refs, [horizontal, vertical] = plane.dimensions.map((axis)=>center[axis - 1]);
            const tickLength = 12, tick = {
                type: 'line',
                line: {
                    color: '#2a7fb8',
                    width: 1
                },
                opacity: 0.8
            };
            return (image ? [
                0
            ] : [
                -thickness / 2,
                0,
                thickness / 2
            ]).flatMap((offset)=>[
                    0,
                    1
                ].flatMap((edge)=>{
                    const start = edge ? -tickLength : 0, end = edge ? 0 : tickLength;
                    const marker = {
                        ...tick,
                        line: {
                            ...tick.line,
                            dash: offset ? '3px,3px' : 'solid'
                        }
                    };
                    return [
                        {
                            ...marker,
                            xref: x,
                            x0: horizontal + offset,
                            x1: horizontal + offset,
                            yref: `${y} domain`,
                            ysizemode: 'pixel',
                            yanchor: edge,
                            y0: start,
                            y1: end
                        },
                        {
                            ...marker,
                            yref: y,
                            y0: vertical + offset,
                            y1: vertical + offset,
                            xref: `${x} domain`,
                            xsizemode: 'pixel',
                            xanchor: edge,
                            x0: start,
                            x1: end
                        }
                    ];
                }));
        });
        for (const slider of layout.sliders ?? [])if (slider.name?.startsWith('slice-')) slider.active = nearest(slider.steps, center[Number(slider.name.slice(6)) - 1]);
        const menu = layout.updatemenus?.find((item)=>item.name === 'slab-thickness');
        if (menu) menu.active = menu.buttons.findIndex((button)=>button.args[0] === thickness);
    }
    function query() {
        return {
            center: [
                ...center
            ],
            thickness,
            ranges: model.planes.map((plane)=>plane.axes.map((axis)=>[
                        ...graph._fullLayout[axis].range
                    ])),
            sizes: model.planes.map((plane)=>plane.axes.map((axis)=>graph._fullLayout[axis]._length))
        };
    }
    function localData() {
        if (image) return graph.data.map((trace, index)=>{
            const plane = model.planes[index % model.planes.length], [a, b] = plane.dimensions.map((d)=>d - 1);
            const values = model.values[Math.floor(index / model.planes.length)], z = [], customdata = [];
            for(let row = 0; row < model.shape[b]; row++){
                const line = [], points = [];
                for(let column = 0; column < model.shape[a]; column++){
                    const point = [
                        ...center
                    ];
                    point[a] = column;
                    point[b] = row;
                    line.push(values[point[0] + model.shape[0] * (point[1] + model.shape[1] * point[2])]);
                    points.push([
                        point[0],
                        point[1],
                        model.zlabels[point[2]]
                    ]);
                }
                z.push(line);
                customdata.push(points);
            }
            return {
                z,
                customdata
            };
        });
        const current = model.frames[frame], coordinates = current.coordinates;
        const entry = new Set(current.entry), exit = new Set(current.exit);
        const velocity = [
            'speed',
            'vx',
            'vy',
            'vz'
        ].includes(property);
        const otherFrame = frame === 0 ? Math.min(1, model.frames.length - 1) : frame - 1;
        const other = model.frames[otherFrame].coordinates, dt = model.times[frame] - model.times[otherFrame];
        const reset = frame === 0 ? new Set(model.frames[otherFrame].entry) : entry;
        const value = (i)=>{
            if (!velocity) return model.values[property][i];
            if (!dt || reset.has(i + 1)) return 0;
            if (property === 'speed') return Math.hypot(...coordinates.map((v, axis)=>v[i] - other[axis][i])) / Math.abs(dt);
            const axis = [
                'vx',
                'vy',
                'vz'
            ].indexOf(property);
            return (coordinates[axis][i] - other[axis][i]) / dt;
        };
        const data = model.planes.flatMap((plane, panel)=>{
            const groups = [
                [],
                [],
                []
            ], normal = plane.normal - 1;
            for(let i = 0; i < coordinates[normal].length; i++){
                if (Math.abs(coordinates[normal][i] - center[normal]) > thickness / 2) continue;
                const group = entry.has(i + 1) ? 2 : exit.has(i + 1) ? 1 : 0;
                groups[group].push(i);
            }
            return groups.map((ids, group)=>{
                const trace = structuredClone(model.properties[property]);
                const values = ids.map(value);
                trace.uid = `phantom-slice-${panel + 1}-${group + 1}`;
                [trace.xaxis, trace.yaxis] = plane.refs;
                [trace.x, trace.y] = plane.dimensions.map((axis)=>ids.map((i)=>coordinates[axis - 1][i]));
                trace.customdata = ids.map((i)=>i + 1);
                trace.ids = trace.customdata.map(String);
                trace.text = values.map((v)=>Math.round(v * 1e4) / 1e4);
                trace.marker.color = values;
                trace.marker.showscale = false;
                trace.showlegend = velocity && panel === 2 && group > 0;
                const colorbar = trace.marker.colorbar;
                trace.hovertemplate = `Spin %{customdata}<br>${'xyz'[plane.dimensions[0] - 1]}: %{x:.2f} cm<br>` + `${'xyz'[plane.dimensions[1] - 1]}: %{y:.2f} cm<br>` + `<b>${colorbar.title}</b>: %{text}${colorbar.ticksuffix}<extra></extra>`;
                if (group > 0) {
                    trace.name = group === 1 ? 'Exit' : 'Re-entry';
                    if (velocity) {
                        trace.marker.color = group === 1 ? '#ff40df' : '#00e5ff';
                        trace.marker.size *= 2;
                    }
                    if (!ids.length) trace.x = trace.y = [
                        null
                    ];
                }
                return trace;
            });
        });
        const marker = structuredClone(model.properties[property].marker);
        marker.colorbar.title = {
            text: marker.colorbar.title
        };
        marker.color = [
            marker.cmin,
            marker.cmax
        ];
        marker.size = marker.opacity = 0;
        data.push({
            type: 'scatter',
            uid: 'phantom-slice-colorbar',
            x: [
                0,
                0
            ],
            y: [
                0,
                0
            ],
            mode: 'markers',
            marker,
            showlegend: false,
            hoverinfo: 'skip'
        });
        return data;
    }
    async function update() {
        if (disposed) return;
        if (!image && !onChange && model.frames.length > 1 && !graph.layout.updatemenus.some((menu)=>menu.name === 'motion-playback')) {
            const timeSlider = graph.layout.sliders.find((item)=>item.name === 'time');
            graph.layout.updatemenus.push({
                name: 'motion-playback',
                type: 'buttons',
                showactive: false,
                x: 0,
                xanchor: 'left',
                y: 0,
                yanchor: 'top',
                pad: {
                    t: timeSlider.pad.t + 4
                },
                buttons: [
                    {
                        name: 'motion-playback',
                        label: '▶',
                        method: 'skip',
                        args: []
                    }
                ]
            });
        }
        const controls = ()=>{
            const layout = {
                sliders: graph.layout.sliders.map((control)=>({
                        ...control
                    })),
                updatemenus: (graph.layout.updatemenus ?? []).map((control)=>({
                        ...control
                    })),
                annotations: graph.layout.annotations.map((annotation)=>({
                        ...annotation
                    }))
            };
            applyLayout(layout);
            const changes = {
                shapes: layout.shapes,
                annotations: layout.annotations
            };
            if (image) {
                const size = graph._fullLayout._size, [x, y, z] = model.shape.map((n, i)=>n * model.spacing[i]);
                const minimumPanelWidth = Math.min(128, (size.w - 24) / 2), minimumPanelHeight = Math.min(100, (size.h - 60) / 2);
                const scale = Math.min((size.w - 24) / (x + z), (size.h - 60) / (y + z), (size.w - 24 - minimumPanelWidth) / Math.max(x, z), (size.h - 60 - minimumPanelHeight) / Math.max(y, z));
                const widths = [
                    x,
                    z
                ].map((value)=>Math.max(minimumPanelWidth, value * scale));
                const heights = [
                    y,
                    z
                ].map((value)=>Math.max(minimumPanelHeight, value * scale));
                const left = (size.w - widths[0] - widths[1] - 24) / 2, bottom = (size.h - heights[0] - heights[1] - 60) / 2;
                const right = left + widths[0] + 24, top = bottom + heights[1] + 60;
                const panels = [
                    [
                        left,
                        top,
                        widths[0],
                        heights[0]
                    ],
                    [
                        right,
                        top,
                        widths[1],
                        heights[0]
                    ],
                    [
                        left,
                        bottom,
                        widths[0],
                        heights[1]
                    ]
                ];
                model.planes.forEach((plane, i)=>{
                    const [a, b, w, h] = panels[i], [imageWidth, imageHeight] = plane.dimensions.map((axis)=>model.shape[axis - 1] * model.spacing[axis - 1] * scale);
                    const imageLeft = a + (w - imageWidth) / 2, imageBottom = b + (h - imageHeight) / 2;
                    changes[`${plane.axes[0]}.domain`] = [
                        imageLeft / size.w,
                        (imageLeft + imageWidth) / size.w
                    ];
                    changes[`${plane.axes[1]}.domain`] = [
                        imageBottom / size.h,
                        (imageBottom + imageHeight) / size.h
                    ];
                    changes[`sliders[${i}].x`] = a / size.w;
                    changes[`sliders[${i}].len`] = w / size.w;
                    changes[`sliders[${i}].y`] = imageBottom / size.h;
                    for (const [name, edge] of [
                        [
                            'slice-title',
                            0
                        ],
                        [
                            'slice-position',
                            1
                        ]
                    ]){
                        const annotation = layout.annotations.find((item)=>item.name === `${name}-${plane.normal}`);
                        if (annotation) Object.assign(annotation, {
                            xref: 'paper',
                            yref: 'paper',
                            yanchor: 'bottom',
                            yshift: 4,
                            x: (a + edge * w) / size.w,
                            y: (imageBottom + imageHeight) / size.h
                        });
                    }
                });
                changes.imageColorbar = {
                    x: (right + widths[1] / 2) / size.w,
                    y: (bottom + heights[1] / 2) / size.h,
                    len: heights[1] / size.h
                };
            }
            for (const key of [
                'sliders',
                'updatemenus'
            ])layout[key]?.forEach((control, i)=>changes[`${key}[${i}].active`] = control.active);
            return changes;
        };
        if (onChange) {
            await Plotly.relayout(graph, controls());
            onChange();
            return;
        }
        if (rendering) {
            pending = true;
            return;
        }
        rendering = true;
        try {
            do {
                pending = false;
                const data = localData(), attributes = {};
                for (const key of image ? [
                    'z',
                    'customdata'
                ] : [
                    'x',
                    'y',
                    'customdata',
                    'ids',
                    'text',
                    'name',
                    'hovertemplate',
                    'showlegend'
                ])attributes[key] = data.map((trace)=>trace[key]);
                for (const key of image ? [] : [
                    'color',
                    'colorscale',
                    'colorbar',
                    'cmin',
                    'cmax',
                    'size',
                    'showscale'
                ])attributes[`marker.${key}`] = data.map((trace)=>trace.marker[key]);
                const changes = controls(), colorbar = changes.imageColorbar;
                delete changes.imageColorbar;
                if (colorbar) for (const [key, value] of Object.entries(colorbar))attributes[`colorbar.${key}`] = data.map(()=>value);
                await Plotly.update(graph, attributes, changes, data.map((_, i)=>i));
            }while (pending && !disposed)
        } finally{
            rendering = false;
        }
    }
    function changed() {
        cancelAnimationFrame(updateFrame);
        updateFrame = requestAnimationFrame(update);
    }
    function wheel(event) {
        if (event.ctrlKey || event.metaKey || !event.deltaY) return;
        const rect = graph.getBoundingClientRect(), x = event.clientX - rect.left, y = event.clientY - rect.top;
        const panel = model.planes.findIndex((plane)=>{
            const [a, b] = plane.axes.map((axis)=>graph._fullLayout[axis]);
            return a._offset <= x && x <= a._offset + a._length && b._offset <= y && y <= b._offset + b._length;
        });
        if (panel < 0) return;
        event.preventDefault();
        event.stopImmediatePropagation();
        const axis = model.planes[panel].normal - 1, slider = positionSliders()[panel];
        const index = Math.max(0, Math.min(slider.steps.length - 1, nearest(slider.steps, center[axis]) + Math.sign(event.deltaY)));
        center[axis] = Number(slider.steps[index].value);
        changed();
    }
    function click(event) {
        if (image && event.event?.detail > 1) return;
        const point = event.points?.[0];
        const plane = model.planes.find((item)=>item.refs[0] === point?.xaxis?._id);
        if (!plane) return;
        const select = ()=>{
            plane.dimensions.forEach((axis, row)=>center[axis - 1] = clamp(row === 0 ? point.x : point.y, axis - 1));
            changed();
        };
        if (!image) return select();
        clearTimeout(clickTimer);
        clickTimer = setTimeout(select, graph._context.doubleClickDelay);
    }
    const doubleClick = ()=>clearTimeout(clickTimer);
    function slider(event) {
        if (event.slider.name?.startsWith('slice-')) {
            center[Number(event.slider.name.slice(6)) - 1] = Number(event.step.value);
            changed();
        } else if (!onChange && event.slider.name === 'time') {
            frame = Number(event.step.value) - 1;
            changed();
        }
    }
    function labelSliders() {
        graph.querySelectorAll('.slider-group').forEach((group)=>{
            const control = group.__data__;
            if (!control.name?.startsWith('slice-')) return;
            const handle = group.querySelector('.slider-grip-rect');
            handle.setAttribute('tabindex', '0');
            handle.setAttribute('role', 'slider');
            handle.setAttribute('aria-label', `${'xyz'[Number(control.name.slice(6)) - 1]} ${image ? 'image' : 'slice'} position`);
            handle.setAttribute('aria-orientation', 'horizontal');
            handle.setAttribute('aria-valuemin', control.steps[0].label);
            handle.setAttribute('aria-valuemax', control.steps.at(-1).label);
            handle.setAttribute('aria-valuenow', control.steps[control.active].label);
            handle.setAttribute('aria-valuetext', `${control.steps[control.active].label}${image ? '' : ' mm'}`);
        });
    }
    function focusSlider(event) {
        const group = event.target.closest('.slider-group');
        if (group?.__data__?.name?.startsWith('slice-')) group.querySelector('.slider-grip-rect').focus({
            preventScroll: true
        });
    }
    function key(event) {
        const control = event.target.closest('.slider-group')?.__data__;
        if (!control?.name?.startsWith('slice-') || event.altKey || event.ctrlKey || event.metaKey) return;
        if (event.key !== 'ArrowLeft' && event.key !== 'ArrowRight') return;
        event.preventDefault();
        event.stopPropagation();
        const axis = Number(control.name.slice(6)) - 1;
        const index = Math.max(0, Math.min(control.steps.length - 1, nearest(control.steps, center[axis]) + (event.key === 'ArrowRight' ? 1 : -1)));
        center[axis] = Number(control.steps[index].value);
        changed();
    }
    function button(event) {
        if (event.button.name === 'slab-thickness') {
            thickness = event.button.args[0];
            changed();
        }
        if (!onChange && event.button.name === 'phantom-property') {
            property = event.button.args[0];
            for (const menu of graph.layout.updatemenus)if (menu.name === 'phantom-property') menu.active = menu.buttons.findIndex((item)=>item.args[0] === property);
            changed();
        }
        if (!onChange && event.button.name === 'motion-playback') {
            playing = !playing;
            playStarted = performance.now() - 5000 * (model.times[frame] - model.times[0]) / (model.times.at(-1) - model.times[0]);
            clearTimeout(playTimer);
            const index = graph.layout.updatemenus.findIndex((menu)=>menu.name === 'motion-playback');
            Plotly.relayout(graph, {
                [`updatemenus[${index}].buttons[0].label`]: playing ? '❚❚' : '▶'
            });
            if (playing) play();
        }
    }
    function play() {
        if (!playing || disposed) return;
        const time = model.times[0] + (performance.now() - playStarted) % 5000 / 5000 * (model.times.at(-1) - model.times[0]);
        frame = Math.max(0, model.times.findLastIndex((t)=>t <= time));
        const timeSlider = graph.layout.sliders.find((item)=>item.name === 'time');
        timeSlider.active = frame;
        changed();
        playTimer = setTimeout(play, 30);
    }
    graph.addEventListener('wheel', wheel, {
        capture: true,
        passive: false
    });
    graph.addEventListener('pointerdown', focusSlider, true);
    graph.addEventListener('keydown', key);
    graph.on('plotly_click', click);
    graph.on('plotly_sliderchange', slider);
    graph.on('plotly_buttonclicked', button);
    graph.on('plotly_doubleclick', doubleClick);
    graph.on('plotly_afterplot', labelSliders);
    function resized(event) {
        if (image && [
            'width',
            'height',
            'autosize'
        ].some((key)=>key in event)) changed();
    }
    graph.on('plotly_relayout', resized);
    labelSliders();
    const removal = new MutationObserver(()=>{
        if (!graph.isConnected) dispose();
    });
    removal.observe(document, {
        childList: true,
        subtree: true
    });
    function dispose() {
        disposed = true;
        clearTimeout(playTimer);
        clearTimeout(clickTimer);
        cancelAnimationFrame(updateFrame);
        removal.disconnect();
        graph.removeEventListener('wheel', wheel, true);
        graph.removeEventListener('pointerdown', focusSlider, true);
        graph.removeEventListener('keydown', key);
        graph.removeListener?.('plotly_click', click);
        graph.removeListener?.('plotly_sliderchange', slider);
        graph.removeListener?.('plotly_buttonclicked', button);
        graph.removeListener?.('plotly_afterplot', labelSliders);
        graph.removeListener?.('plotly_relayout', resized);
        graph.removeListener?.('plotly_doubleclick', doubleClick);
    }
    if (!onChange) changed();
    return {
        query,
        applyLayout,
        dispose,
        async rebind (next) {
            playing = false;
            clearTimeout(playTimer);
            clearTimeout(clickTimer);
            if (JSON.stringify(next.center) !== JSON.stringify(model.center)) center = [
                ...next.center
            ];
            else center = center.map((value, axis)=>Math.max(next.bounds[axis][0], Math.min(next.bounds[axis][1], value)));
            if (next.thickness !== model.thickness) thickness = next.thickness;
            model = next;
            image = Boolean(model.shape);
            if (!onChange) {
                property = (graph.layout.updatemenus ?? []).filter((menu)=>menu.name === 'phantom-property').map((menu)=>menu.buttons[menu.active ?? 0]?.args[0]).find((value)=>value !== undefined) ?? graph.layout.meta.phantom_property;
                frame = graph.layout.sliders.find((slider)=>slider.name === 'time')?.active ?? 0;
                await update();
            }
        }
    };
}
function multiply(a, b) {
    return Array.from({
        length: 16
    }, (_, i)=>{
        const row = i % 4, col = Math.floor(i / 4);
        return [
            0,
            1,
            2,
            3
        ].reduce((sum, k)=>sum + a[4 * k + row] * b[4 * col + k], 0);
    });
}
function projection(graph, dimensions) {
    const layout = graph._fullLayout;
    if (dimensions.length === 3) {
        const scene = layout.scene._scene;
        const camera = scene.glplot.cameraParams;
        scene.getCamera();
        const matrix = multiply(camera.projection, multiply(scene.camera.matrix, camera.model));
        for(let col = 0; col < 3; col++)for(let row = 0; row < 4; row++)matrix[4 * col + row] *= scene.dataScale[col];
        return matrix;
    }
    const matrix = Array(16).fill(0);
    matrix[15] = 1;
    [
        'xaxis',
        'yaxis'
    ].forEach((name, row)=>{
        const [lo, hi] = layout[name].range;
        matrix[4 * (dimensions[row] - 1) + row] = 2 / (hi - lo);
        matrix[12 + row] = -(hi + lo) / (hi - lo);
    });
    return matrix;
}
const PLAYBACK_DURATION_MS = 5000;
function playbackFrame(steps, elapsed) {
    const start = Number(steps[0].label), end = Number(steps.at(-1).label);
    const time = start + elapsed % 5000 / 5000 * (end - start);
    let lo = 0, hi = steps.length - 1;
    while(lo < hi){
        const mid = Math.ceil((lo + hi) / 2);
        if (Number(steps[mid].label) <= time) lo = mid;
        else hi = mid - 1;
    }
    return lo + 1;
}
async function mountKspace(root, graph, Plotly, config) {
    let data = structuredClone(graph.data);
    if (!root.style.height) root.style.height = `${graph.clientHeight}px`;
    if (!root.style.width && graph.layout.width) root.style.width = `${graph.layout.width}px`;
    const range = mountKspaceRange(root, graph, config, filter);
    const theme = observeTheme(graph, Plotly);
    let selection = [];
    function bound(times, value, upper) {
        let lo = 0, hi = times.length;
        while(lo < hi){
            const mid = lo + hi >>> 1;
            if (times[mid] < value || upper && times[mid] === value) lo = mid + 1;
            else hi = mid;
        }
        return lo;
    }
    async function filter() {
        const time = range.query.time_range;
        const selections = [
            config.times.trajectory,
            config.times.adc
        ].map((times)=>[
                bound(times, time[0], false),
                bound(times, time[1], true)
            ]);
        const indices = selections.flat();
        if (indices.every((value, i)=>value === selection[i])) return;
        const update = {};
        for (const key of [
            'x',
            'y',
            'z',
            'text'
        ])if (data.some((trace)=>key in trace)) update[key] = data.slice(0, 2).map((trace, i)=>trace[key]?.slice(...selections[i]) ?? []);
        update['marker.color'] = [
            null,
            data[1].marker.color.slice(...selections[1])
        ];
        const view = {};
        if (graph.layout.scene) {
            const scene = graph._fullLayout.scene._scene;
            scene.glplot.selection.object = null;
            Plotly.Fx.unhover(graph);
            view['scene.camera'] = scene.getCamera();
        } else {
            for (const axis of [
                'xaxis',
                'yaxis'
            ]){
                view[`${axis}.range`] = [
                    ...graph._fullLayout[axis].range
                ];
                view[`${axis}.autorange`] = false;
                view[`${axis}.constrain`] = 'domain';
            }
        }
        await Plotly.update(graph, update, view, [
            0,
            1
        ]);
        selection = indices;
    }
    const fit = ()=>{
        if (graph.layout.width !== graph.clientWidth || graph.layout.height !== graph.clientHeight) Plotly.relayout(graph, {
            width: graph.clientWidth,
            height: graph.clientHeight
        });
    };
    fit();
    const resize = new ResizeObserver(fit);
    resize.observe(root);
    const removed = new MutationObserver(()=>{
        if (root.isConnected) return;
        removed.disconnect();
        resize.disconnect();
        theme.disconnect();
        range.dispose();
        Plotly.purge(graph);
    });
    removed.observe(document, {
        childList: true,
        subtree: true
    });
    return {
        dispose () {
            removed.disconnect();
            resize.disconnect();
            theme.disconnect();
            range.dispose();
        },
        async rebind (nextConfig) {
            config = nextConfig;
            range.rebind(config);
            data = structuredClone(graph.data);
            selection = [];
            if (range.query.time_range) await filter();
        }
    };
}
async function mount(root, Plotly, request, response, templateJSON) {
    if (!root.isConnected) return;
    let template = JSON.parse(templateJSON);
    const graph = document.createElement('div');
    graph.style.cssText = 'width:100%;height:100%';
    root.prepend(graph);
    let frame = 1, coil = 1, component = 1, revision = 0;
    let disposed = false, rendering = false, inFlight = false, wanted = false, timer;
    let preview = false, dragging = false, interacting = false;
    let pending = null, renderTask;
    let playing = false, playTimer, cameraFrame, cameraTimer, playStarted, drawnFrame = 1;
    const initial = JSON.parse(response.value);
    let definition = reconcilePlot(graph, initial.data, template.layout);
    let property = initial.property;
    let range = template.kspace ? mountKspaceRange(root, graph, template.kspace, fetchDetail) : null;
    let slices;
    const plotSize = ()=>range ? [
            graph.clientWidth,
            graph.clientHeight
        ] : [
            root.clientWidth,
            root.clientHeight
        ];
    let timeSlider, playbackMenu;
    function playbackControls() {
        playbackMenu = undefined;
        timeSlider = template.layout.sliders?.find((slider)=>slider.name === 'time' && slider.visible !== false && slider.steps?.length > 1);
        if (!timeSlider) return;
        timeSlider.pad = {
            ...timeSlider.pad,
            l: 55
        };
        template.layout.updatemenus ??= [];
        playbackMenu = template.layout.updatemenus.length;
        template.layout.updatemenus.push({
            name: 'motion-playback',
            type: 'buttons',
            showactive: false,
            x: 0,
            xanchor: 'left',
            y: 0,
            yanchor: 'top',
            pad: {
                t: 34
            },
            buttons: [
                {
                    name: 'motion-playback',
                    label: '▶',
                    method: 'skip',
                    args: []
                }
            ]
        });
        template.layout.updatemenus[playbackMenu].pad.t = (timeSlider.pad.t ?? 30) + 4;
    }
    playbackControls();
    function schedulePlayback() {
        clearTimeout(playTimer);
        if (!playing || disposed || rendering || inFlight || wanted || pending || interacting) return;
        playTimer = setTimeout(()=>{
            frame = playbackFrame(timeSlider.steps, performance.now() - playStarted);
            if (frame !== drawnFrame) fetchPreview();
            else schedulePlayback();
        }, 16);
    }
    async function setPlaying(value, refine = true) {
        playing = value;
        clearTimeout(playTimer);
        clearTimeout(timer);
        revision++;
        wanted = false;
        pending = null;
        if (playing) {
            const start = Number(timeSlider.steps[0].label), end = Number(timeSlider.steps.at(-1).label);
            playStarted = performance.now() - PLAYBACK_DURATION_MS * (Number(timeSlider.steps[drawnFrame - 1].label) - start) / (end - start);
        }
        const label = playing ? '❚❚' : '▶';
        template.layout.updatemenus[playbackMenu].buttons[0].label = label;
        await Plotly.relayout(graph, {
            [`updatemenus[${playbackMenu}].buttons[0].label`]: label
        });
        await styleControls1();
        if (!playing && refine) {
            frame = drawnFrame;
            fetchDetail();
        }
        schedulePlayback();
    }
    function send() {
        if (disposed || inFlight || rendering || interacting || !wanted) return;
        wanted = false;
        inFlight = true;
        const size = graph._fullLayout._size;
        request.notify({
            id: revision,
            frame,
            coil,
            component,
            preview,
            property,
            ...range?.query,
            ...slices ? {
                slices: slices.query()
            } : {},
            visible: graph.data.map((trace)=>trace.visible !== false && trace.visible !== 'legendonly'),
            width: size.w,
            height: size.h,
            projection: projection(graph, template.dimensions)
        });
    }
    function fetchDetail() {
        if (disposed) return;
        revision++;
        preview = false;
        wanted = true;
        clearTimeout(playTimer);
        clearTimeout(timer);
        send();
    }
    function fetchPreview(refineDelay = 120) {
        if (disposed) return;
        if (!preview) revision++;
        preview = true;
        wanted = true;
        clearTimeout(playTimer);
        clearTimeout(timer);
        send();
        if (!playing && !dragging) timer = setTimeout(()=>fetchDetail(), refineDelay);
    }
    async function styleControls1() {
        await styleControls(graph, Plotly, (menu, index)=>menu.name === 'motion-playback' ? playing : menu.name === 'phantom-property' ? menu.buttons[index].args[0] === property : index === component - 1);
        graph.querySelectorAll('.updatemenu-header-group').forEach((group)=>{
            const playback = group.__data__.name === 'motion-playback';
            if (playback) group.querySelectorAll('.updatemenu-button').forEach((button)=>button.setAttribute('aria-label', playing ? 'Pause motion' : 'Play motion'));
        });
        graph.querySelectorAll('.slider-group').forEach((group)=>{
            if (group.__data__.name !== 'time') return;
            const slider = graph.layout.sliders[group.__data__._index];
            const handle = group.querySelector('.slider-grip-rect');
            handle.setAttribute('tabindex', '0');
            handle.setAttribute('role', 'slider');
            handle.setAttribute('aria-label', 'Motion time');
            handle.setAttribute('aria-orientation', 'horizontal');
            handle.setAttribute('aria-valuemin', '1');
            handle.setAttribute('aria-valuemax', String(slider.steps.length));
            handle.setAttribute('aria-valuenow', String(slider.active + 1));
            handle.setAttribute('aria-valuetext', `${slider.steps[slider.active].label} ms`);
        });
    }
    async function draw(payload, rebind = false) {
        const layout = themeLayout(root, structuredClone(template.layout));
        [layout.width, layout.height] = plotSize();
        layout.uirevision = `spatial-plot-${definition.revision}`;
        (layout.sliders ?? []).forEach((slider, i)=>{
            if (slider.name !== 'time' && slider.name !== 'coil') return;
            slider.active = slider.name === 'coil' ? coil - 1 : payload.preview && !playing ? graph.layout?.sliders?.find((item)=>item.name === 'time')?.active ?? frame - 1 : (payload.frame ?? frame) - 1;
        });
        for (const menu of layout.updatemenus ?? []){
            if (menu.name === 'phantom-property') menu.active = menu.buttons.findIndex((button)=>button.args[0] === property);
            if (menu.name === 'coil-component') menu.active = component - 1;
        }
        slices?.applyLayout(layout);
        if (payload.coloraxis) layout.coloraxis = payload.coloraxis;
        const reduced = payload.coarse || payload.shown < payload.available;
        const counts = payload.slice_counts ? sliceAnnotations(layout.meta.phantom_slices, payload.slice_counts, layout.sliders.find((slider)=>slider.name === 'time').steps[(payload.frame ?? frame) - 1].label) : [];
        layout.annotations = [
            ...layout.annotations ?? [],
            ...counts,
            ...!payload.slice_counts && reduced ? [
                {
                    text: payload.coarse ? 'Coarser grid' : `Downsampled · ${(100 * payload.shown / payload.available).toFixed(2)}%`,
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
                        size: 10
                    },
                    hovertext: payload.detail ?? (payload.coarse ? `Zoom to refine the grid to ${1000 * payload.spacing} mm spacing.` : `${payload.shown.toLocaleString()} / ${payload.available.toLocaleString()} ` + 'visible spins shown. ' + (payload.preview ? 'Motion preview; pause to refine.' : 'Zoom for more detail.'))
                }
            ] : []
        ];
        identifyTraces(payload.data);
        if (graph._fullLayout && !rebind) {
            const old = new Map(graph.data.map((trace)=>[
                    trace.uid,
                    trace
                ]));
            payload.data = payload.data.map((trace)=>old.get(trace.uid)?.visible === 'legendonly' ? old.get(trace.uid) : trace);
        }
        if (rebind) {
            await Plotly.react(graph, payload.data, layout, {
                ...template.config,
                responsive: true
            });
        } else if (!graph._fullLayout) {
            await Plotly.newPlot(graph, payload.data, layout, {
                ...template.config,
                responsive: true
            });
        } else {
            const data = {}, changes = {};
            const equal = (a, b)=>a === b || Array.isArray(a) && Array.isArray(b) && a.length === b.length && a.every((value, i)=>value === b[i]);
            for (const key of [
                'x',
                'y',
                'z',
                'ids',
                'text',
                'customdata',
                'name',
                'hovertemplate',
                'showlegend'
            ])if (payload.data.some((trace, i)=>!equal(trace[key], graph.data[i][key]))) data[key] = payload.data.map((trace)=>trace[key]);
            for (const key of [
                'color',
                'colorscale',
                'colorbar',
                'cmin',
                'cmax',
                'size'
            ])if (payload.data.some((trace, i)=>key === 'color' ? !equal(trace.marker?.[key], graph.data[i].marker?.[key]) : JSON.stringify(trace.marker?.[key]) !== JSON.stringify(graph.data[i].marker?.[key]))) data[`marker.${key}`] = payload.data.map((trace)=>trace.marker?.[key]);
            for (const key of [
                'annotations',
                'coloraxis',
                'shapes'
            ])if (JSON.stringify(layout[key]) !== JSON.stringify(graph.layout[key])) changes[key] = layout[key];
            for (const key of [
                'sliders',
                'updatemenus'
            ])layout[key]?.forEach((control, i)=>{
                if (control.active !== graph.layout[key]?.[i]?.active) changes[`${key}[${i}].active`] = control.active;
            });
            if (Object.keys(data).length || Object.keys(changes).length) {
                if (template.dimensions.length === 3) {
                    const scene = graph._fullLayout.scene._scene;
                    scene.glplot.selection.object = null;
                    Plotly.Fx.unhover(graph);
                    changes['scene.camera'] = scene.getCamera();
                }
                await Plotly.update(graph, data, changes, payload.data.map((_, i)=>i));
            }
        }
        drawnFrame = payload.frame ?? frame;
        await styleControls1();
    }
    async function render() {
        if (rendering || disposed || interacting) return;
        rendering = true;
        try {
            while(pending && !disposed && !interacting){
                const payload = pending;
                pending = null;
                if (payload.id === revision) {
                    await draw(payload);
                }
            }
        } finally{
            rendering = false;
            if (disposed) Plotly.purge(graph);
            else {
                send();
                schedulePlayback();
            }
        }
    }
    const theme = observeTheme(graph, Plotly);
    await draw(initial);
    if (!root.isConnected) {
        Plotly.purge(graph);
        return;
    }
    slices = mountSlices(graph, Plotly, fetchPreview);
    response.on((serialized)=>{
        if (disposed) return false;
        inFlight = false;
        const payload = JSON.parse(serialized);
        if (payload.id === revision) pending = payload;
        if (!rendering) renderTask = render().catch((error)=>console.error('Spatial plot update failed', error));
    });
    function refreshView() {
        cancelAnimationFrame(cameraFrame);
        cameraFrame = requestAnimationFrame(()=>{
            if (playing) {
                revision++;
                frame = playbackFrame(timeSlider.steps, performance.now() - playStarted);
                fetchPreview();
            } else fetchDetail();
        });
    }
    graph.on('plotly_relayout', async (event)=>{
        if ((!rendering || interacting) && Object.keys(event).some((key)=>/^(scene\.|[xy]axis\d*\.)/.test(key))) {
            await resetSpatialAxes(graph, Plotly, event);
            refreshView();
        }
    });
    graph.addEventListener('wheel', (event)=>{
        if (!event.target.closest('.gl-container, .nsewdrag')) return;
        beginInteraction();
        cameraTimer = setTimeout(finishInteraction, 120);
    }, true);
    graph.on('plotly_restyle', ([changes])=>{
        if (!rendering && Object.hasOwn(changes, 'visible')) fetchDetail();
    });
    graph.on('plotly_sliderchange', (event)=>{
        if (event.slider.name === 'time') {
            if (playing) setPlaying(false, false);
            frame = Number(event.step.value);
            fetchPreview();
        } else if (event.slider.name === 'coil') {
            coil = Number(event.step.value);
            fetchDetail();
        }
    });
    graph.on('plotly_sliderstart', (event)=>{
        if (event.slider.name !== 'time') return;
        dragging = true;
        if (playing) setPlaying(false, false);
        clearTimeout(timer);
    });
    graph.on('plotly_sliderend', (event)=>{
        if (event.slider.name !== 'time') return;
        dragging = false;
        frame = Number(event.step.value);
        fetchDetail();
    });
    graph.addEventListener('keydown', (event)=>{
        const slider = event.target.closest('.slider-group');
        if (slider?.__data__?.name !== 'time' || event.altKey || event.ctrlKey || event.metaKey) return;
        if (event.key !== 'ArrowLeft' && event.key !== 'ArrowRight') return;
        event.preventDefault();
        event.stopPropagation();
        if (playing) {
            frame = drawnFrame;
            setPlaying(false, false);
        }
        frame = Math.max(1, Math.min(timeSlider.steps.length, frame + (event.key === 'ArrowRight' ? 1 : -1)));
        Plotly.relayout(graph, {
            [`sliders[${slider.__data__._index}].active`]: frame - 1
        }).then(styleControls1);
        fetchPreview(300);
    });
    graph.on('plotly_buttonclicked', (event)=>{
        if (event.button.name === 'motion-playback') return setPlaying(!playing);
        if (event.button.name === 'phantom-property') {
            property = event.button.args[0];
            revision++;
        } else if (event.button.name === 'coil-component') component = event.button.args[0];
        else return;
        styleControls1();
        if (playing) fetchPreview();
        else fetchDetail();
    });
    function beginInteraction() {
        interacting = true;
        revision++;
        pending = null;
        wanted = false;
        clearTimeout(playTimer);
        clearTimeout(cameraTimer);
    }
    graph.addEventListener('pointerdown', (event)=>{
        const slider = event.target.closest('.slider-group');
        if (slider?.__data__?.name === 'time') slider.querySelector('.slider-grip-rect').focus({
            preventScroll: true
        });
        if (!event.target.closest('.gl-container, .nsewdrag')) return;
        beginInteraction();
    }, true);
    const finishInteraction = ()=>{
        if (!interacting) return;
        interacting = false;
        clearTimeout(cameraTimer);
        pending = null;
        revision++;
        refreshView();
    };
    document.addEventListener('pointerup', finishInteraction, true);
    document.addEventListener('pointercancel', finishInteraction, true);
    let size = plotSize();
    const resize = new ResizeObserver(()=>{
        const next = plotSize();
        if (next.every((value, i)=>Math.abs(value - size[i]) < 2)) return;
        size = next;
        Plotly.relayout(graph, {
            width: next[0],
            height: next[1]
        }).then(async ()=>{
            await styleControls1();
            playing ? fetchPreview() : fetchDetail();
        });
    });
    const removal = new MutationObserver(()=>{
        if (root.isConnected) return;
        disposed = true;
        pending = null;
        clearTimeout(timer);
        clearTimeout(playTimer);
        clearTimeout(cameraTimer);
        cancelAnimationFrame(cameraFrame);
        resize.disconnect();
        removal.disconnect();
        theme.disconnect();
        range?.dispose();
        slices?.dispose();
        document.removeEventListener('pointerup', finishInteraction, true);
        document.removeEventListener('pointercancel', finishInteraction, true);
        if (!rendering) Plotly.purge(graph);
    });
    removal.observe(document, {
        subtree: true,
        childList: true
    });
    resize.observe(root);
    if (range) resize.observe(graph);
    fetchDetail();
    return {
        graph,
        async rebind (initialJSON, templateJSON) {
            revision++;
            pending = null;
            wanted = false;
            clearTimeout(timer);
            clearTimeout(playTimer);
            clearTimeout(cameraTimer);
            cancelAnimationFrame(cameraFrame);
            playing = dragging = interacting = false;
            await renderTask;
            const next = JSON.parse(initialJSON);
            template = JSON.parse(templateJSON);
            for (const menu of graph.layout.updatemenus ?? []){
                if (menu.name === 'phantom-property') menu.active = menu.buttons.findIndex((button)=>button.args[0] === property);
                if (menu.name === 'coil-component') menu.active = component - 1;
            }
            for (const slider of graph.layout.sliders ?? [])if (slider.name === 'time' || slider.name === 'coil') slider.active = (slider.name === 'time' ? frame : coil) - 1;
            definition = reconcilePlot(graph, next.data, template.layout, definition);
            property = template.layout.updatemenus?.filter((menu)=>menu.name === 'phantom-property').map((menu)=>menu.buttons[menu.active ?? 0]?.args[0]).find((value)=>value !== undefined) ?? next.property;
            next.property = property;
            const menus = template.layout.updatemenus ?? [];
            component = (menus.find((menu)=>menu.name === 'coil-component')?.active ?? 0) + 1;
            coil = (template.layout.sliders?.find((slider)=>slider.name === 'coil')?.active ?? 0) + 1;
            frame = (template.layout.sliders?.find((slider)=>slider.name === 'time')?.active ?? 0) + 1;
            playbackControls();
            if (slices && template.layout.meta?.phantom_slices) await slices.rebind(template.layout.meta.phantom_slices);
            else {
                slices?.dispose();
                slices = null;
            }
            if (range && template.kspace) range.rebind(template.kspace);
            else {
                range?.dispose();
                graph.style.cssText = 'width:100%;height:100%';
                range = template.kspace ? mountKspaceRange(root, graph, template.kspace, fetchDetail) : null;
            }
            inFlight = false;
            rendering = true;
            try {
                await draw(next, true);
            } finally{
                rendering = false;
            }
            slices ??= mountSlices(graph, Plotly, fetchPreview);
            fetchDetail();
        }
    };
}
export { projection as projection };
export { playbackFrame as playbackFrame };
export { mountKspace as mountKspace };
export { mount as mount };
