import {themeLayout, observeTheme, styleControls as stylePlotControls, identifyTraces, reconcilePlot, resetSpatialAxes} from './PlotTheme.js';
import {mountKspaceRange} from './KspaceRange.js';
import {mountSlices, sliceAnnotations} from './PhantomSlices.js';

function multiply(a, b) {
    return Array.from({length: 16}, (_, i) => {
        const row = i % 4, col = Math.floor(i / 4);
        return [0, 1, 2, 3].reduce((sum, k) => sum + a[4 * k + row] * b[4 * col + k], 0);
    });
}

// Plotly's pinned gl3d renderer projects scaled coordinates with these matrices.
// Keep this private-API boundary here, rather than approximating its camera in Julia.
export function projection(graph, dimensions) {
    const layout = graph._fullLayout;
    if (dimensions.length === 3) {
        const scene = layout.scene._scene;
        const camera = scene.glplot.cameraParams;
        scene.getCamera();
        const matrix = multiply(camera.projection, multiply(scene.camera.matrix, camera.model));
        for (let col = 0; col < 3; col++)
            for (let row = 0; row < 4; row++) matrix[4 * col + row] *= scene.dataScale[col];
        return matrix;
    }
    const matrix = Array(16).fill(0);
    matrix[15] = 1;
    ['xaxis', 'yaxis'].forEach((name, row) => {
        const [lo, hi] = layout[name].range;
        matrix[4 * (dimensions[row] - 1) + row] = 2 / (hi - lo);
        matrix[12 + row] = -(hi + lo) / (hi - lo);
    });
    return matrix;
}

const PLAYBACK_DURATION_MS = 5000;
const KEYBOARD_REFINE_DELAY_MS = 300;

export function playbackFrame(steps, elapsed) {
    const start = Number(steps[0].label), end = Number(steps.at(-1).label);
    const time = start + (elapsed % PLAYBACK_DURATION_MS) / PLAYBACK_DURATION_MS * (end - start);
    let lo = 0, hi = steps.length - 1;
    while (lo < hi) {
        const mid = Math.ceil((lo + hi) / 2);
        if (Number(steps[mid].label) <= time) lo = mid;
        else hi = mid - 1;
    }
    return lo + 1;
}

export async function mountKspace(root, graph, Plotly, config) {
    let data = structuredClone(graph.data);
    if (!root.style.height) root.style.height = `${graph.clientHeight}px`;
    if (!root.style.width && graph.layout.width) root.style.width = `${graph.layout.width}px`;
    const range = mountKspaceRange(root, graph, config, filter);
    const theme = observeTheme(graph, Plotly);
    let selection = [];
    function bound(times, value, upper) {
        let lo = 0, hi = times.length;
        while (lo < hi) {
            const mid = (lo + hi) >>> 1;
            if (times[mid] < value || (upper && times[mid] === value)) lo = mid + 1;
            else hi = mid;
        }
        return lo;
    }
    async function filter() {
        const time = range.query.time_range;
        const selections = [config.times.trajectory, config.times.adc].map(times =>
            [bound(times,time[0],false), bound(times,time[1],true)]);
        const indices = selections.flat();
        if (indices.every((value, i) => value === selection[i])) return;
        const update = {};
        for (const key of ['x','y','z','text'])
            if (data.some(trace => key in trace)) update[key] = data.slice(0,2).map((trace,i) =>
                trace[key]?.slice(...selections[i]) ?? []);
        update['marker.color'] = [null, data[1].marker.color.slice(...selections[1])];
        const view = {};
        if (graph.layout.scene) {
            const scene = graph._fullLayout.scene._scene;
            scene.glplot.selection.object = null;
            Plotly.Fx.unhover(graph);
            view['scene.camera'] = scene.getCamera();
        } else {
            for (const axis of ['xaxis','yaxis']) {
                view[`${axis}.range`] = [...graph._fullLayout[axis].range];
                view[`${axis}.autorange`] = false;
                view[`${axis}.constrain`] = 'domain';
            }
        }
        await Plotly.update(graph, update, view, [0,1]);
        selection = indices;
    }
    const fit = () => {
        if (graph.layout.width !== graph.clientWidth || graph.layout.height !== graph.clientHeight)
            Plotly.relayout(graph, {width:graph.clientWidth, height:graph.clientHeight});
    };
    fit();
    const resize = new ResizeObserver(fit);
    resize.observe(root);
    const removed = new MutationObserver(() => {
        if (root.isConnected) return;
        removed.disconnect(); resize.disconnect(); theme.disconnect(); range.dispose(); Plotly.purge(graph);
    });
    removed.observe(document, {childList:true, subtree:true});
    return {dispose() {
        removed.disconnect(); resize.disconnect(); theme.disconnect(); range.dispose();
    }, async rebind(nextConfig) {
        config = nextConfig;
        range.rebind(config);
        data = structuredClone(graph.data);
        selection = [];
        if (range.query.time_range) await filter();
    }};
}

export async function mount(root, Plotly, request, response, templateJSON) {
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
    const plotSize = () => range ? [graph.clientWidth, graph.clientHeight] : [root.clientWidth, root.clientHeight];
    let timeSlider, playbackMenu;
    function playbackControls() {
        playbackMenu = undefined;
        timeSlider = template.layout.sliders?.find(slider => slider.name === 'time' &&
            slider.visible !== false && slider.steps?.length > 1);
        if (!timeSlider) return;
        timeSlider.pad = {...timeSlider.pad, l: 55};
        template.layout.updatemenus ??= [];
        playbackMenu = template.layout.updatemenus.length;
        template.layout.updatemenus.push({name: 'motion-playback', type: 'buttons',
            showactive: false, x: 0, xanchor: 'left', y: 0, yanchor: 'top', pad: {t: 34},
            buttons: [{name: 'motion-playback', label: '▶', method: 'skip', args: []}]});
        template.layout.updatemenus[playbackMenu].pad.t = (timeSlider.pad.t ?? 30) + 4;
    }
    playbackControls();

    function schedulePlayback() {
        clearTimeout(playTimer);
        if (!playing || disposed || rendering || inFlight || wanted || pending || interacting) return;
        playTimer = setTimeout(() => {
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
            playStarted = performance.now() - PLAYBACK_DURATION_MS *
                (Number(timeSlider.steps[drawnFrame - 1].label) - start) / (end - start);
        }
        const label = playing ? '❚❚' : '▶';
        template.layout.updatemenus[playbackMenu].buttons[0].label = label;
        await Plotly.relayout(graph, {[`updatemenus[${playbackMenu}].buttons[0].label`]: label});
        await styleControls();
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
        request.notify({id: revision, frame, coil, component, preview, property, ...range?.query,
            ...(slices ? {slices:slices.query()} : {}),
            visible: graph.data.map(trace => trace.visible !== false && trace.visible !== 'legendonly'),
            width: size.w, height: size.h, projection: projection(graph, template.dimensions)});
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
        // Accept intermediate motion frames while coalescing the next requested time.
        if (!preview) revision++;
        preview = true;
        wanted = true;
        clearTimeout(playTimer);
        clearTimeout(timer);
        send();
        // Only refinement is debounced; previews never wait for the pointer to stop.
        if (!playing && !dragging) timer = setTimeout(() => fetchDetail(), refineDelay);
    }
    async function styleControls() {
        await stylePlotControls(graph, Plotly, (menu, index) => menu.name === 'motion-playback' ? playing :
            menu.name === 'phantom-property' ? menu.buttons[index].args[0] === property : index === component - 1);
        graph.querySelectorAll('.updatemenu-header-group').forEach(group => {
            const playback = group.__data__.name === 'motion-playback';
            if (playback) group.querySelectorAll('.updatemenu-button').forEach(button =>
                button.setAttribute('aria-label', playing ? 'Pause motion' : 'Play motion'));
        });
        graph.querySelectorAll('.slider-group').forEach(group => {
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
        (layout.sliders ?? []).forEach((slider, i) => {
            // Do not pull a dragging thumb back to an intermediate preview frame.
            if (slider.name !== 'time' && slider.name !== 'coil') return;
            slider.active = slider.name === 'coil' ? coil - 1 :
                payload.preview && !playing ? (graph.layout?.sliders?.find(item => item.name === 'time')?.active ?? frame - 1) :
                (payload.frame ?? frame) - 1;
        });
        for (const menu of layout.updatemenus ?? []) {
            if (menu.name === 'phantom-property') menu.active = menu.buttons.findIndex(button => button.args[0] === property);
            if (menu.name === 'coil-component') menu.active = component - 1;
        }
        slices?.applyLayout(layout);
        if (payload.coloraxis) layout.coloraxis = payload.coloraxis;
        const reduced = payload.coarse || payload.shown < payload.available;
        const counts = payload.slice_counts ? sliceAnnotations(layout.meta.phantom_slices, payload.slice_counts,
            layout.sliders.find(slider => slider.name === 'time').steps[(payload.frame ?? frame)-1].label) : [];
        layout.annotations = [...(layout.annotations ?? []), ...counts, ...(!payload.slice_counts && reduced ? [{
            text: payload.coarse ? 'Coarser grid' : `Downsampled · ${(100 * payload.shown / payload.available).toFixed(2)}%`,
            xref: 'paper', yref: 'paper', x: 1, y: 0, xanchor: 'right', yanchor: 'bottom',
            xshift: -2, yshift: 2, showarrow: false, borderpad: 0, font: {size: 10},
            hovertext: payload.detail ?? (payload.coarse ? `Zoom to refine the grid to ${1000 * payload.spacing} mm spacing.` :
                `${payload.shown.toLocaleString()} / ${payload.available.toLocaleString()} ` +
                'visible spins shown. ' +
                (payload.preview ? 'Motion preview; pause to refine.' : 'Zoom for more detail.')),
        }] : [])];
        // Emptying a hidden trace removes its legend entry in Plotly.
        // Retain its bounded browser data until the user shows it again.
        identifyTraces(payload.data);
        if (graph._fullLayout && !rebind) {
            const old = new Map(graph.data.map(trace => [trace.uid, trace]));
            payload.data = payload.data.map(trace => old.get(trace.uid)?.visible === 'legendonly' ? old.get(trace.uid) : trace);
        }
        if (rebind) {
            await Plotly.react(graph, payload.data, layout, {...template.config, responsive: true});
        } else if (!graph._fullLayout) {
            await Plotly.newPlot(graph, payload.data, layout, {...template.config, responsive: true});
        } else {
            const data = {}, changes = {};
            const equal = (a, b) => a === b || (Array.isArray(a) && Array.isArray(b) &&
                a.length === b.length && a.every((value, i) => value === b[i]));
            for (const key of ['x', 'y', 'z', 'ids', 'text', 'customdata', 'name', 'hovertemplate', 'showlegend'])
                if (payload.data.some((trace, i) => !equal(trace[key], graph.data[i][key])))
                    data[key] = payload.data.map(trace => trace[key]);
            for (const key of ['color', 'colorscale', 'colorbar', 'cmin', 'cmax', 'size'])
                if (payload.data.some((trace, i) => key === 'color' ?
                    !equal(trace.marker?.[key], graph.data[i].marker?.[key]) :
                    JSON.stringify(trace.marker?.[key]) !== JSON.stringify(graph.data[i].marker?.[key])))
                    data[`marker.${key}`] = payload.data.map(trace => trace.marker?.[key]);
            for (const key of ['annotations', 'coloraxis', 'shapes'])
                if (JSON.stringify(layout[key]) !== JSON.stringify(graph.layout[key])) changes[key] = layout[key];
            for (const key of ['sliders', 'updatemenus']) layout[key]?.forEach((control, i) => {
                if (control.active !== graph.layout[key]?.[i]?.active)
                    changes[`${key}[${i}].active`] = control.active;
            });
            if (Object.keys(data).length || Object.keys(changes).length) {
                if (template.dimensions.length === 3) {
                    const scene = graph._fullLayout.scene._scene;
                    // A pick indexes the old WebGL buffers, not the replacement points.
                    scene.glplot.selection.object = null;
                    Plotly.Fx.unhover(graph);
                    // Preserve the live camera, which can be newer than gd.layout.
                    changes['scene.camera'] = scene.getCamera();
                }
                await Plotly.update(graph, data, changes, payload.data.map((_, i) => i));
            }
        }
        drawnFrame = payload.frame ?? frame;
        await styleControls();
    }
    async function render() {
        if (rendering || disposed || interacting) return;
        rendering = true;
        try {
            while (pending && !disposed && !interacting) {
                const payload = pending;
                pending = null;
                if (payload.id === revision) {
                    await draw(payload);
                }
            }
        } finally {
            rendering = false;
            if (disposed) Plotly.purge(graph);
            else { send(); schedulePlayback(); }
        }
    }
    const theme = observeTheme(graph, Plotly);
    await draw(initial);
    if (!root.isConnected) { Plotly.purge(graph); return; }
    slices = mountSlices(graph, Plotly, fetchPreview);
    response.on(serialized => {
        if (disposed) return false;
        inFlight = false;
        const payload = JSON.parse(serialized);
        if (payload.id === revision) pending = payload;
        if (!rendering) renderTask = render().catch(error => console.error('Spatial plot update failed', error));
    });
    function refreshView() {
        cancelAnimationFrame(cameraFrame);
        // gl3d publishes camera matrices on its next render, after the input event.
        cameraFrame = requestAnimationFrame(() => {
            if (playing) {
                revision++;
                frame = playbackFrame(timeSlider.steps, performance.now() - playStarted);
                fetchPreview();
            } else fetchDetail();
        });
    }
    graph.on('plotly_relayout', async event => {
        if ((!rendering || interacting) && Object.keys(event).some(key => /^(scene\.|[xy]axis\d*\.)/.test(key))) {
            await resetSpatialAxes(graph, Plotly, event);
            refreshView();
        }
    });
    graph.addEventListener('wheel', event => {
        if (!event.target.closest('.gl-container, .nsewdrag')) return;
        beginInteraction();
        cameraTimer = setTimeout(finishInteraction, 120);
    }, true);
    graph.on('plotly_restyle', ([changes]) => {
        if (!rendering && Object.hasOwn(changes, 'visible')) fetchDetail();
    });
    graph.on('plotly_sliderchange', event => {
        if (event.slider.name === 'time') {
            if (playing) setPlaying(false, false);
            frame = Number(event.step.value);
            fetchPreview();
        }
        else if (event.slider.name === 'coil') { coil = Number(event.step.value); fetchDetail(); }
    });
    graph.on('plotly_sliderstart', event => {
        if (event.slider.name !== 'time') return;
        dragging = true;
        if (playing) setPlaying(false, false);
        clearTimeout(timer);
    });
    graph.on('plotly_sliderend', event => {
        if (event.slider.name !== 'time') return;
        dragging = false;
        frame = Number(event.step.value);
        fetchDetail();
    });
    graph.addEventListener('keydown', event => {
        const slider = event.target.closest('.slider-group');
        if (slider?.__data__?.name !== 'time' || event.altKey || event.ctrlKey || event.metaKey) return;
        if (event.key !== 'ArrowLeft' && event.key !== 'ArrowRight') return;
        event.preventDefault();
        event.stopPropagation();
        if (playing) {
            frame = drawnFrame;
            setPlaying(false, false);
        }
        frame = Math.max(1, Math.min(timeSlider.steps.length,
            frame + (event.key === 'ArrowRight' ? 1 : -1)));
        Plotly.relayout(graph, {[`sliders[${slider.__data__._index}].active`]: frame - 1})
            .then(styleControls);
        fetchPreview(KEYBOARD_REFINE_DELAY_MS);
    });
    graph.on('plotly_buttonclicked', event => {
        if (event.button.name === 'motion-playback') return setPlaying(!playing);
        if (event.button.name === 'phantom-property') {
            property = event.button.args[0];
            revision++;
        } else if (event.button.name === 'coil-component') component = event.button.args[0];
        else return;
        styleControls();
        if (playing) fetchPreview();
        else fetchDetail();
    });
    // Let native camera gestures finish without replacing their WebGL buffers mid-drag.
    function beginInteraction() {
        interacting = true;
        revision++;
        pending = null;
        wanted = false;
        clearTimeout(playTimer);
        clearTimeout(cameraTimer);
    }
    graph.addEventListener('pointerdown', event => {
        const slider = event.target.closest('.slider-group');
        if (slider?.__data__?.name === 'time')
            slider.querySelector('.slider-grip-rect').focus({preventScroll: true});
        if (!event.target.closest('.gl-container, .nsewdrag')) return;
        beginInteraction();
    }, true);
    const finishInteraction = () => {
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
    const resize = new ResizeObserver(() => {
        const next = plotSize();
        if (next.every((value, i) => Math.abs(value - size[i]) < 2)) return;
        size = next;
        Plotly.relayout(graph, {width: next[0], height: next[1]}).then(async () => {
            await styleControls();
            playing ? fetchPreview() : fetchDetail();
        });
    });
    const removal = new MutationObserver(() => {
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
    removal.observe(document, {subtree: true, childList: true});
    resize.observe(root);
    if (range) resize.observe(graph);
    fetchDetail();
    return {graph, async rebind(initialJSON, templateJSON) {
        revision++;
        pending = null;
        wanted = false;
        clearTimeout(timer); clearTimeout(playTimer); clearTimeout(cameraTimer);
        cancelAnimationFrame(cameraFrame);
        playing = dragging = interacting = false;
        await renderTask;
        const next = JSON.parse(initialJSON);
        template = JSON.parse(templateJSON);
        for (const menu of graph.layout.updatemenus ?? []) {
            if (menu.name === 'phantom-property')
                menu.active = menu.buttons.findIndex(button => button.args[0] === property);
            if (menu.name === 'coil-component') menu.active = component - 1;
        }
        for (const slider of graph.layout.sliders ?? [])
            if (slider.name === 'time' || slider.name === 'coil') slider.active = (slider.name === 'time' ? frame : coil) - 1;
        definition = reconcilePlot(graph, next.data, template.layout, definition);
        property = template.layout.updatemenus?.filter(menu => menu.name === 'phantom-property')
            .map(menu => menu.buttons[menu.active ?? 0]?.args[0]).find(value => value !== undefined) ?? next.property;
        next.property = property;
        const menus = template.layout.updatemenus ?? [];
        component = (menus.find(menu => menu.name === 'coil-component')?.active ?? 0) + 1;
        coil = (template.layout.sliders?.find(slider => slider.name === 'coil')?.active ?? 0) + 1;
        frame = (template.layout.sliders?.find(slider => slider.name === 'time')?.active ?? 0) + 1;
        playbackControls();
        if (slices && template.layout.meta?.phantom_slices) await slices.rebind(template.layout.meta.phantom_slices);
        else { slices?.dispose(); slices = null; }
        if (range && template.kspace) range.rebind(template.kspace);
        else {
            range?.dispose();
            graph.style.cssText = 'width:100%;height:100%';
            range = template.kspace ? mountKspaceRange(root, graph, template.kspace, fetchDetail) : null;
        }
        inFlight = false;
        rendering = true;
        try { await draw(next, true); }
        finally { rendering = false; }
        slices ??= mountSlices(graph, Plotly, fetchPreview);
        fetchDetail();
    }};
}
