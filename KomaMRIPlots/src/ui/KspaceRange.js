export function mountKspaceRange(root, graph, config, initial, onChange) {
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
    const numbers = [...footer.querySelectorAll('input[type=number]')];
    const sliders = [...footer.querySelectorAll('input[type=range]')];
    const limits = [...footer.querySelectorAll('.range-limits span')];
    const before = footer.querySelector('.range-before'), after = footer.querySelector('.range-after');
    const ranges = structuredClone(initial.filter_ranges);
    const control = {query: {}, update, theme};
    let colors;

    function display() {
        const max = unit.value === 'adc' ? config.adc_count - 1 : config.time_max;
        const step = unit.value === 'adc' ? 1 : 'any';
        const bounds = ranges[unit.value];
        numbers.forEach((input, i) => {
            input.min = 0; input.max = max; input.step = step;
            if (input !== document.activeElement)
                input.value = unit.value === 'adc' ? bounds[i] : Number(bounds[i].toFixed(3));
        });
        sliders.forEach((input, i) => {
            input.min = 0; input.max = config.adc_count - 1; input.step = 1; input.value = ranges.adc[i];
        });
        const adcMax = config.adc_count - 1;
        const percentage = value => `${100 * Math.max(0, Math.min(adcMax, value)) / (adcMax || 1)}%`;
        before.style.cssText = `left:0;width:${percentage(ranges.adc[0])}`;
        after.style.cssText = `right:0;width:${percentage(adcMax - ranges.adc[1])}`;
        limits[0].textContent = 'ADC order: 0';
        limits[1].textContent = adcMax.toLocaleString();
    }
    function change(input, i, selectedUnit) {
        const max = selectedUnit === 'adc' ? config.adc_count - 1 : config.time_max;
        const bounds = ranges[selectedUnit].map(value => Math.max(0, Math.min(max, value)));
        let value = Math.max(0, Math.min(max, Number(input.value)));
        if (selectedUnit === 'adc') value = Math.round(value);
        bounds[i] = i === 0 ? Math.min(value, bounds[1]) : Math.max(value, bounds[0]);
        ranges[selectedUnit] = bounds;
        control.query = {[`${selectedUnit}_range`]: [...bounds]};
        display();
        onChange();
    }
    function update(payload) {
        Object.assign(ranges, payload.filter_ranges);
        footer.querySelector('output').textContent =
            `${payload.selected_adc.toLocaleString()} / ${config.adc_count.toLocaleString()} samples`;
        display();
    }
    function theme() {
        footer.style.setProperty('--range-text', graph.layout.font.color);
        footer.style.setProperty('--range-bg', graph.layout.paper_bgcolor);
        footer.style.setProperty('--range-border', graph.layout.xaxis?.gridcolor ?? graph.layout.scene.xaxis.gridcolor);
        colors ??= graph._fullData[1].marker?.colorscale;
        if (colors) footer.querySelector('.range-gradient').style.background =
            `linear-gradient(to right, ${colors.map(([position, color]) => `${color} ${100 * position}%`).join(',')})`;
    }
    numbers.forEach((input, i) => input.addEventListener('change', () => change(input, i, unit.value)));
    sliders.forEach((input, i) => input.addEventListener('input', () => change(input, i, 'adc')));
    unit.addEventListener('change', display);
    update(initial);
    return control;
}
