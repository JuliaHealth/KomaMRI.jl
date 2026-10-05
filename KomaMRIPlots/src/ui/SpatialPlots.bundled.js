// deno-fmt-ignore-file
// deno-lint-ignore-file
// This code was bundled using `deno bundle` and it's not recommended to edit it manually

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
function mountKspaceRange(root, graph, config, initial, onChange) {
    root.style.display = 'flex';
    root.style.flexDirection = 'column';
    graph.style.cssText = 'width:100%;flex:1;min-height:0';
    const footer = document.createElement('div');
    footer.className = 'koma-kspace-range';
    footer.innerHTML = `
        <style>
        .koma-kspace-range {flex:none;padding:8px 24px 12px;font:13px sans-serif;color:var(--range-text)}
        .koma-kspace-range .range-fields {display:flex;align-items:center;gap:12px;flex-wrap:wrap}
        .koma-kspace-range label {display:flex;align-items:center;gap:6px}
        .koma-kspace-range select,.koma-kspace-range input[type=number] {
            font:inherit;color:inherit;background:var(--range-bg);border:1px solid var(--range-border);
            border-radius:4px;padding:4px 6px;box-sizing:border-box}
        .koma-kspace-range input[type=number] {width:100px}
        .koma-kspace-range output {margin-left:auto;font-size:12px}
        .koma-kspace-range .range-track {position:relative;height:28px;margin:10px 8px 0}
        .koma-kspace-range .range-gradient {position:absolute;inset:10px 0;height:8px;border-radius:4px}
        .koma-kspace-range .range-mask {position:absolute;top:10px;height:8px;background:var(--range-bg);opacity:.75}
        .koma-kspace-range input[type=range] {position:absolute;inset:0;width:100%;height:28px;
            margin:0;background:none;pointer-events:none;appearance:none;-webkit-appearance:none}
        .koma-kspace-range input[type=range]::-webkit-slider-runnable-track {height:8px;background:transparent}
        .koma-kspace-range input[type=range]::-moz-range-track {height:8px;background:transparent}
        .koma-kspace-range input[type=range]::-webkit-slider-thumb {appearance:none;-webkit-appearance:none;
            width:16px;height:20px;margin-top:-6px;border:2px solid var(--range-text);border-radius:4px;
            background:var(--range-bg);pointer-events:auto;cursor:ew-resize}
        .koma-kspace-range input[type=range]::-moz-range-thumb {width:12px;height:16px;
            border:2px solid var(--range-text);border-radius:4px;background:var(--range-bg);pointer-events:auto;cursor:ew-resize}
        .koma-kspace-range .range-limits {display:flex;justify-content:space-between;font-size:11px}
        </style>
        <div class="range-fields">
            <label>Filter <select aria-label="K-space filter units"><option value="adc">ADC order</option><option value="time">Time (ms)</option></select></label>
            <label>From <input type="number" aria-label="K-space range start"></label>
            <label>To <input type="number" aria-label="K-space range end"></label>
            <output aria-live="polite"></output>
        </div>
        <div class="range-track">
            <div class="range-gradient"></div><div class="range-mask range-before"></div><div class="range-mask range-after"></div>
            <input type="range" aria-label="ADC order lower bound"><input type="range" aria-label="ADC order upper bound">
        </div>
        <div class="range-limits"><span></span><span></span></div>`;
    root.append(footer);
    const unit = footer.querySelector('select');
    const numbers = [
        ...footer.querySelectorAll('input[type=number]')
    ];
    const sliders = [
        ...footer.querySelectorAll('input[type=range]')
    ];
    const limits = [
        ...footer.querySelectorAll('.range-limits span')
    ];
    const before = footer.querySelector('.range-before'), after = footer.querySelector('.range-after');
    const ranges = structuredClone(initial.filter_ranges);
    const control = {
        query: {},
        update,
        theme
    };
    let colors;
    function display() {
        const max = unit.value === 'adc' ? config.adc_count - 1 : config.time_max;
        const step = unit.value === 'adc' ? 1 : 'any';
        const bounds = ranges[unit.value];
        numbers.forEach((input, i)=>{
            input.min = 0;
            input.max = max;
            input.step = step;
            if (input !== document.activeElement) input.value = unit.value === 'adc' ? bounds[i] : Number(bounds[i].toFixed(3));
        });
        sliders.forEach((input, i)=>{
            input.min = 0;
            input.max = config.adc_count - 1;
            input.step = 1;
            input.value = ranges.adc[i];
        });
        const adcMax = config.adc_count - 1;
        const percentage = (value)=>`${100 * Math.max(0, Math.min(adcMax, value)) / (adcMax || 1)}%`;
        before.style.cssText = `left:0;width:${percentage(ranges.adc[0])}`;
        after.style.cssText = `right:0;width:${percentage(adcMax - ranges.adc[1])}`;
        limits[0].textContent = 'ADC order: 0';
        limits[1].textContent = adcMax.toLocaleString();
    }
    function change(input, i, selectedUnit) {
        const max = selectedUnit === 'adc' ? config.adc_count - 1 : config.time_max;
        const bounds = ranges[selectedUnit].map((value)=>Math.max(0, Math.min(max, value)));
        let value = Math.max(0, Math.min(max, Number(input.value)));
        if (selectedUnit === 'adc') value = Math.round(value);
        bounds[i] = i === 0 ? Math.min(value, bounds[1]) : Math.max(value, bounds[0]);
        ranges[selectedUnit] = bounds;
        control.query = {
            [`${selectedUnit}_range`]: [
                ...bounds
            ]
        };
        display();
        onChange();
    }
    function update(payload) {
        Object.assign(ranges, payload.filter_ranges);
        footer.querySelector('output').textContent = `${payload.selected_adc.toLocaleString()} / ${config.adc_count.toLocaleString()} samples`;
        display();
    }
    function theme() {
        footer.style.setProperty('--range-text', graph.layout.font.color);
        footer.style.setProperty('--range-bg', graph.layout.paper_bgcolor);
        footer.style.setProperty('--range-border', graph.layout.xaxis?.gridcolor ?? graph.layout.scene.xaxis.gridcolor);
        colors ??= graph._fullData[1].marker?.colorscale;
        if (colors) footer.querySelector('.range-gradient').style.background = `linear-gradient(to right, ${colors.map(([position, color])=>`${color} ${100 * position}%`).join(',')})`;
    }
    numbers.forEach((input, i)=>input.addEventListener('change', ()=>change(input, i, unit.value)));
    sliders.forEach((input, i)=>input.addEventListener('input', ()=>change(input, i, 'adc')));
    unit.addEventListener('change', display);
    update(initial);
    return control;
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
async function mount(root, Plotly, request, response, templateJSON) {
    if (!root.isConnected) return;
    const template = JSON.parse(templateJSON);
    const graph = document.createElement('div');
    graph.style.cssText = 'width:100%;height:100%';
    root.append(graph);
    let frame = 1, coil = 1, component = 1, revision = 0;
    let disposed = false, rendering = false, inFlight = false, wanted = false, timer;
    let preview = false, dragging = false, interacting = false;
    let pending = null;
    let playing = false, playTimer, cameraFrame, cameraTimer, playStarted, drawnFrame = 1;
    const initial = JSON.parse(response.value);
    let property = initial.property;
    const range = template.kspace ? mountKspaceRange(root, graph, template.kspace, initial, fetchDetail) : null;
    const plotSize = ()=>range ? [
            graph.clientWidth,
            graph.clientHeight
        ] : [
            root.clientWidth,
            root.clientHeight
        ];
    const timeSlider = template.layout.sliders?.find((slider)=>slider.name === 'time' && slider.visible !== false && slider.steps?.length > 1);
    let playbackMenu;
    if (timeSlider) {
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
    }
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
        styleControls();
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
    function styleControls() {
        const propertyGroups = [];
        graph.querySelectorAll('.updatemenu-header-group').forEach((group)=>{
            const playback = group.__data__.name === 'motion-playback';
            const phantom = group.__data__.name === 'phantom-property';
            if (phantom) propertyGroups.push(group);
            if (!playback && !phantom && group.__data__.name !== 'coil-component') return;
            group.querySelectorAll('.updatemenu-button').forEach((button, index)=>{
                button.classList.add('koma-component-button');
                button.setAttribute('role', 'button');
                button.setAttribute('aria-pressed', String(playback ? playing : phantom ? group.__data__.buttons[index].args[0] === property : index === component - 1));
                if (playback) button.setAttribute('aria-label', playing ? 'Pause motion' : 'Play motion');
            });
        });
        if (propertyGroups.length) {
            propertyGroups.forEach((group)=>group.style.transform = '');
            const bounds = propertyGroups.map((group)=>group.getBoundingClientRect());
            const viewport = graph.getBoundingClientRect();
            const offset = (viewport.left + viewport.right - bounds[0].left - bounds.at(-1).right) / 2;
            propertyGroups.forEach((group)=>group.style.transform = `translateX(${offset}px)`);
            graph.querySelectorAll('.annotation').forEach((annotation)=>{
                if (annotation.textContent === '—') annotation.style.transform = `translateX(${offset}px)`;
            });
        }
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
    async function draw(payload) {
        if (range) {
            range.update(payload);
            payload.data[1].marker.showscale = false;
        }
        const layout = themeLayout(root, structuredClone(template.layout));
        [layout.width, layout.height] = plotSize();
        layout.uirevision = 'spatial-plot';
        (layout.sliders ?? []).forEach((slider, i)=>{
            slider.active = slider.name !== 'time' ? coil - 1 : payload.preview && !playing ? graph.layout?.sliders[i].active ?? frame - 1 : (payload.frame ?? frame) - 1;
        });
        for (const menu of layout.updatemenus ?? [])menu.active = menu.name === 'phantom-property' ? menu.buttons.findIndex((button)=>button.args[0] === property) : component - 1;
        if (payload.coloraxis) layout.coloraxis = payload.coloraxis;
        const reduced = payload.coarse || payload.shown < payload.available;
        layout.annotations = [
            ...layout.annotations ?? [],
            ...reduced ? [
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
        if (!graph._fullLayout) {
            await Plotly.newPlot(graph, payload.data, layout, {
                ...template.config,
                responsive: true
            });
        } else {
            payload.data = payload.data.map((trace, i)=>graph.data[i].visible === 'legendonly' ? graph.data[i] : trace);
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
                'coloraxis'
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
        styleControls();
        range?.theme();
    }
    async function render() {
        if (rendering || disposed || interacting) return;
        rendering = true;
        try {
            while(pending && !disposed && !interacting){
                const payload = pending;
                pending = null;
                if (payload.id === revision) await draw(payload);
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
    response.on((serialized)=>{
        if (disposed) return false;
        inFlight = false;
        const payload = JSON.parse(serialized);
        if (payload.id === revision) pending = payload;
        render().catch((error)=>Bonito.send_error('Spatial plot update failed', error));
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
    graph.on('plotly_relayout', (event)=>{
        if ((!rendering || interacting) && Object.keys(event).some((key)=>/^(scene\.|[xy]axis\.)/.test(key))) {
            refreshView();
        }
    });
    if (range) graph.on('plotly_afterplot', ()=>range.theme());
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
        } else {
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
        }).then(styleControls);
        fetchPreview(300);
    });
    graph.on('plotly_buttonclicked', (event)=>{
        if (event.button.name === 'motion-playback') return setPlaying(!playing);
        if (event.button.name === 'phantom-property') {
            property = event.button.args[0];
            revision++;
        } else component = event.button.args[0];
        styleControls();
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
        }).then(()=>{
            styleControls();
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
}
export { projection as projection };
export { playbackFrame as playbackFrame };
export { mount as mount };
