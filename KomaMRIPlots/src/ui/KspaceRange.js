const FILTER_SETTLE_MS = 150;

export function mountKspaceRange(root, graph, config, onChange) {
    root.style.display = 'flex';
    root.style.flexDirection = 'column';
    graph.style.cssText = 'width:100%;flex:1;min-height:0';
    const navigator = root.querySelector('.koma-time-plot');
    let pending, timer, dragging = false;
    const control = {query: {}, dispose, rebind(next) {
        clearTimeout(timer);
        pending = null;
        dragging = false;
        config = next;
        const interval = navigator.querySelector('.js-plotly-plot')?.layout?.xaxis?.range;
        if (interval) control.query = {time_range: interval.map(t => Math.max(0, Math.min(config.time_max, t)))};
    }};
    control.rebind(config);

    function apply() {
        if (dragging || !pending) return;
        control.query = {time_range: pending};
        pending = null;
        onChange();
    }
    function schedule() {
        clearTimeout(timer);
        if (!dragging) timer = setTimeout(apply, FILTER_SETTLE_MS);
    }
    function change(event) {
        if (event.detail.source !== navigator) return;
        pending = event.detail.range.map(t => Math.max(0, Math.min(config.time_max, t)));
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
