// Run from the repository root: node KomaMRIPlots/test/time-viewer.mjs
import assert from 'node:assert/strict';
import {readFile} from 'node:fs/promises';

class Element {
    style = {}; children = []; handlers = {};
    clientWidth = 640; clientHeight = 400; isConnected = true;
    append(...children) { this.children.push(...children); }
    on(event, callback) { this.handlers[event] = callback; }
    closest() { return null; }
    querySelector() { return this; }
    getBoundingClientRect() { return {top: 0, height: 80}; }
}
globalThis.document = {createElement: () => new Element()};
globalThis.ResizeObserver = globalThis.MutationObserver = class {
    observe() {}
    disconnect() {}
};
const Plotly = {
    async react(graph, data, layout) {
        Object.assign(graph, {data, layout, _fullLayout: {
            margin: {l: 30, r: 10}, yaxis: {range: [-1, 1]}}});
    },
    async newPlot(...args) { return this.react(...args); },
    async relayout(graph, event) {
        const axis = graph.layout.xaxis;
        if (event['xaxis.range']) axis.range = [...event['xaxis.range']];
        if (event['xaxis.rangeslider.range'])
            axis.rangeslider.range = [...event['xaxis.rangeslider.range']];
        // Plotly expands the overview to contain a selection outside its previous bounds.
        if (axis.rangeslider.visible) {
            axis.rangeslider.range = [Math.min(axis.rangeslider.range[0], axis.range[0]),
                Math.max(axis.rangeslider.range[1], axis.range[1])];
        }
        await graph.handlers.plotly_relayout?.(event);
    },
    purge() {}, Plots: {resize() {}}
};
const source = await readFile(new URL('../src/ui/TimePlots.js', import.meta.url), 'utf8');
const {mount} = await import(`data:text/javascript;base64,${Buffer.from(source).toString('base64')}`);
const fullRange = [0, 100];
const payload = {id: 0, full_range: fullRange, shown: 2, available: 2,
    data: [{x: fullRange, y: [0, 1]}], config: {}, layout: {
        xaxis: {range: fullRange, rangeslider: {visible: true}}, yaxis: {}}};
const root = new Element(), requests = [];
await mount(root, Plotly, {notify: value => requests.push(value)},
    {value: JSON.stringify(payload), on() {}}, JSON.stringify(payload));
const [graph, clip] = root.children, navigator = clip.children[0];

// Reset after zooming beyond the sequence restores both the selected interval and overview.
for (const outside of [[-20, 120], [-50, 100]]) {
    await Plotly.relayout(graph, {'xaxis.range': outside});
    assert.deepEqual(navigator.layout.xaxis.rangeslider.range, outside);
    await Plotly.relayout(graph, {'xaxis.autorange': true});
    assert.deepEqual(graph.layout.xaxis.range, fullRange);
    assert.deepEqual(navigator.layout.xaxis.range, fullRange);
    assert.deepEqual(navigator.layout.xaxis.rangeslider.range, fullRange);
    assert.deepEqual(requests.at(-1).range, fullRange);
}

// Selecting a smaller interval still retains the full-sequence overview and updates the plot.
await Plotly.relayout(navigator, {'xaxis.range': [30, 40]});
assert.deepEqual(graph.layout.xaxis.range, [30, 40]);
assert.deepEqual(navigator.layout.xaxis.rangeslider.range, fullRange);
console.log('PASS: time slider selection and overview stay synchronized after reset');
