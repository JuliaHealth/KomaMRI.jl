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
async function renderPlot(graph, Plotly, plot, previous) {
    const layout = themeLayout(graph, plot.layout);
    const definition = reconcilePlot(graph, plot.data, layout, previous);
    graph._komaRendering = true;
    try {
        await graph._komaWindowDraw;
        await Plotly.react(graph, plot.data, layout, {
            ...plot.config,
            responsive: true
        });
        observeImageWindow(graph, Plotly);
        if (graph._transitionData?._frames?.length) await Plotly.deleteFrames(graph);
        if (plot.frames?.length) await Plotly.addFrames(graph, plot.frames);
        await applyControlSelections(graph, Plotly);
    } finally{
        graph._komaRendering = false;
    }
    await graph._komaImageWindow?.();
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
function observeControls(graph, Plotly) {
    observeImageWindow(graph, Plotly);
    const update = ()=>styleControls(graph, Plotly);
    for (const event of [
        'plotly_afterplot',
        'plotly_buttonclicked'
    ])graph.on(event, update);
    graph.on('plotly_relayout', (event)=>resetSpatialAxes(graph, Plotly, event));
    return update;
}
function observeImageWindow(graph, Plotly) {
    if (graph._komaImageWindow || !graph.data.some((trace)=>trace.meta?.image_window)) return;
    graph._komaImageWindow = update;
    let window = null, drag = null, drawing = false;
    const images = ()=>graph.data.flatMap((trace, index)=>trace.meta?.image_window ? [
                index
            ] : []);
    async function update() {
        if (!graph.isConnected || graph._komaRendering || drawing) return;
        const indices = images().filter((index)=>{
            const bounds = window ?? graph.data[index].meta.image_window;
            return graph.data[index].zmin !== bounds[0] || graph.data[index].zmax !== bounds[1];
        });
        if (!indices.length) return;
        drawing = true;
        const bounds = indices.map((index)=>window ?? graph.data[index].meta.image_window);
        try {
            graph._komaWindowDraw = Plotly.restyle(graph, {
                zmin: bounds.map((value)=>value[0]),
                zmax: bounds.map((value)=>value[1]),
                zauto: false
            }, indices);
            await graph._komaWindowDraw;
        } finally{
            drawing = false;
        }
        if (graph.isConnected) update();
    }
    function down(event) {
        if (event.button !== 1 || !event.target.closest('.draglayer')) return;
        const rect = graph.getBoundingClientRect();
        for (const index of images()){
            const trace = graph._fullData[index];
            if (trace.visible !== true) continue;
            const x = graph._fullLayout[trace.xaxis.replace('x', 'xaxis')];
            const y = graph._fullLayout[trace.yaxis.replace('y', 'yaxis')];
            const px = event.clientX - rect.left, py = event.clientY - rect.top;
            if (px < x._offset || px > x._offset + x._length || py < y._offset || py > y._offset + y._length) continue;
            event.preventDefault();
            event.stopImmediatePropagation();
            graph.setPointerCapture(event.pointerId);
            drag = {
                x: event.clientX,
                y: event.clientY,
                width: Math.max(trace.zmax - trace.zmin, Number.EPSILON),
                level: (trace.zmax + trace.zmin) / 2,
                horizontal: x._length,
                vertical: y._length
            };
            break;
        }
    }
    function move(event) {
        if (!drag) return;
        event.preventDefault();
        event.stopImmediatePropagation();
        const width = drag.width * Math.exp((event.clientX - drag.x) / drag.horizontal);
        const level = drag.level + drag.width * (drag.y - event.clientY) / drag.vertical;
        window = [
            level - width / 2,
            level + width / 2
        ];
        update();
    }
    function up(event) {
        if (!drag) return;
        drag = null;
        if (graph.hasPointerCapture(event.pointerId)) graph.releasePointerCapture(event.pointerId);
    }
    function reset() {
        window = null;
        update();
    }
    graph.addEventListener('pointerdown', down, true);
    graph.addEventListener('pointermove', move, true);
    for (const event of [
        'pointerup',
        'pointercancel',
        'lostpointercapture'
    ])graph.addEventListener(event, up);
    graph.on('plotly_doubleclick', reset);
}
function observeResize(graph, Plotly) {
    const resize = new ResizeObserver(()=>{
        if (!graph.isConnected || !graph._fullLayout) return;
        const width = graph.clientWidth, height = graph.clientHeight;
        if (!width || !height || width === graph._fullLayout.width && height === graph._fullLayout.height) return;
        Plotly.relayout(graph, {
            width,
            height
        });
    });
    const removal = new MutationObserver(()=>{
        if (graph.isConnected) return;
        resize.disconnect();
        removal.disconnect();
        Plotly.purge(graph);
    });
    resize.observe(graph);
    removal.observe(document, {
        childList: true,
        subtree: true
    });
    return resize;
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
export { reconcileControls as reconcileControls };
export { applyControlSelections as applyControlSelections };
export { identifyTraces as identifyTraces };
export { reconcilePlot as reconcilePlot };
export { renderPlot as renderPlot };
export { styleControls as styleControls };
export { resetSpatialAxes as resetSpatialAxes };
export { observeControls as observeControls };
export { observeResize as observeResize };
export { themeLayout as themeLayout };
export { themeTraces as themeTraces };
export { observeTheme as observeTheme };
