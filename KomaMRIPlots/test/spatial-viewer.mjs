// Run from the repository root: node KomaMRIPlots/test/spatial-viewer.mjs
import assert from 'node:assert/strict';
import {readFile} from 'node:fs/promises';

class Element {
    style = {}; children = []; handlers = {};
    clientWidth = 640; clientHeight = 400; isConnected = true;
    append(child) { this.children.push(child); }
    on(event, callback) { this.handlers[event] = callback; }
    addEventListener(event, callback) { this.handlers[event] = callback; }
    removeEventListener(event) { delete this.handlers[event]; }
    querySelectorAll() { return []; }
}
const observers = [];
globalThis.document = Object.assign(new Element(), {createElement: () => new Element()});
let now = 0;
globalThis.performance = {now: () => now};
globalThis.requestAnimationFrame = callback => setTimeout(callback, 0);
globalThis.cancelAnimationFrame = clearTimeout;
globalThis.ResizeObserver = globalThis.MutationObserver = class {
    constructor(callback) { observers.push(callback); }
    observe() {}
    disconnect() {}
};
const source = await readFile(new URL('../src/ui/SpatialPlots.js', import.meta.url), 'utf8');
const {mount, projection, playbackFrame} = await import(`data:text/javascript;base64,${Buffer.from(source).toString('base64')}`);
const root = new Element(), requests = [];
let respond, purged = false, updates = 0;
const initial = {id: 0, data: [{x: [0], y: [0], marker: {color: [1]}}], shown: 1, available: 1};
const template = {dimensions: [1, 2], config: {}, layout: {
    xaxis: {range: [-1, 1]}, yaxis: {range: [-1, 1]},
    sliders: [{name: 'time', active: 0, steps: Array.from({length: 5}, (_, i) =>
        ({value: String(i + 1), label: 2 * i}))}], updatemenus: []}};
const Plotly = {
    async newPlot(graph, data, layout) {
        Object.assign(graph, {data, layout, _fullLayout: {...layout, _size: {w: 600, h: 300}}});
    },
    async update(graph, data, changes) {
        updates++;
        for (const [key, value] of Object.entries(data)) {
            if (key === 'marker.color') graph.data[0].marker.color = value[0];
            else graph.data[0][key] = value[0];
        }
        for (const [key, value] of Object.entries(changes)) {
            const control = key.match(/^(sliders|updatemenus)\[(\d+)\]\.active$/);
            if (control) graph.layout[control[1]][control[2]].active = value;
            else graph.layout[key] = value;
        }
    },
    async relayout(graph, changes) {
        for (const [key, value] of Object.entries(changes)) {
            const menu = key.match(/^updatemenus\[(\d+)\]\.buttons\[0\]\.label$/);
            if (menu) graph.layout.updatemenus[menu[1]].buttons[0].label = value;
            else graph.layout[key] = value;
        }
    },
    purge() { purged = true; }, Plots: {async resize() {}}
};
await mount(root, Plotly, {notify: value => requests.push(value)},
    {value: JSON.stringify(initial), on: callback => { respond = callback; }}, JSON.stringify(template));
const graph = root.children[0];
const settle = () => new Promise(resolve => setTimeout(resolve, 160));
const tick = () => new Promise(resolve => setTimeout(resolve, 0));
const reply = (request, extra = {}) => respond(JSON.stringify({...initial,
    id: request.id, frame: request.frame, preview: request.preview,
    available: request.preview ? 2 : 1,
    data: [{...initial.data[0], x: [request.frame]}], ...extra}));
const scrub = frame => {
    graph.layout.sliders[0].active = frame - 1;
    graph.handlers.plotly_sliderchange({slider: {name: 'time'}, step: {value: String(frame)}});
};
reply(requests.at(-1));
await tick();

// Continuous scrubbing displays intermediate previews instead of waiting for an idle pointer.
for (const frame of [2, 3, 4, 5]) {
    scrub(frame);
    assert.equal(requests.at(-1).frame, frame);
    assert.equal(requests.at(-1).preview, true);
    reply(requests.at(-1));
    await tick();
    assert.deepEqual(graph.data[0].x, [frame]);
    assert.match(graph.layout.annotations[0].text, /50.00%/);
}

// A slow response still advances the picture during a drag, without queuing every pointer event.
scrub(2);
const intermediate = requests.at(-1), beforeCoalescing = requests.length;
scrub(3);
scrub(4);
assert.equal(requests.length, beforeCoalescing);
reply(intermediate);
await tick();
assert.deepEqual(graph.data[0].x, [2]);
assert.equal(graph.layout.sliders[0].active, 3);
assert.equal(requests.at(-1).frame, 4);
reply(requests.at(-1));
await tick();

// A delayed preview cannot replace the final selection; stopping requests full detail.
scrub(2);
const delayed = requests.at(-1), beforeRelease = requests.length;
scrub(5);
await settle();
assert.equal(requests.length, beforeRelease);
reply(delayed);
await tick();
assert.deepEqual(graph.data[0].x, [4]);
assert.equal(requests.at(-1).frame, 5);
assert.equal(requests.at(-1).preview, false);
reply(requests.at(-1));
await tick();
assert.deepEqual(graph.data[0].x, [5]);
assert.equal(graph.layout.sliders[0].active, 4);
assert.deepEqual(graph.layout.annotations, []);

// Holding the thumb still must not trigger an expensive full cloud before pointer release.
graph.handlers.plotly_sliderstart({slider:{name:'time'}});
scrub(3);
reply(requests.at(-1));
await tick();
const holding = requests.length;
await settle();
assert.equal(requests.length, holding);
assert.equal(requests.at(-1).preview, true);
graph.handlers.plotly_sliderend({slider:{name:'time'},step:{value:'3'}});
assert.equal(requests.at(-1).preview, false);
reply(requests.at(-1));
await tick();

// Detail replacement preserves a selected viewport and uses the displayed physical axes.
graph._fullLayout.xaxis.range = [2, 4];
graph.handlers.plotly_relayout({'xaxis.range': [2, 4]});
await settle();
assert.equal(requests.at(-1).projection[12], -3);
respond(JSON.stringify({...initial, id: requests.at(-1).id}));
await settle();
assert.deepEqual(graph.layout.xaxis.range, [2, 4]);
assert.deepEqual(projection(graph, [3, 2]).slice(0, 4), [0, 0, 0, 0]);
assert.equal(projection(graph, [3, 2])[8], 1);

// Even a fully loaded cloud is reselected on zoom so offscreen spins cannot consume its budget.
respond(JSON.stringify({...initial, id: requests.at(-1).id, complete: true}));
await settle();
const requestCount = requests.length;
const updateCount = updates;
respond(JSON.stringify({...initial, id: requests.at(-1).id, complete: true}));
await settle();
assert.equal(updates, updateCount);
graph.handlers.plotly_relayout({'xaxis.range': [2.5, 3.5]});
await settle();
assert.equal(requests.length, requestCount + 1);
reply(requests.at(-1));
await tick();
scrub(2);
assert.equal(requests.length, requestCount + 2);

// Playback follows physical time on a five-second clock, skipping frames instead of slowing the loop.
reply(requests.at(-1));
await settle();
reply(requests.at(-1));
await tick();
const playback = {button: {name: 'motion-playback'}};
await graph.handlers.plotly_buttonclicked(playback);
for (const expected of [3, 4, 1, 2]) {
    now += 1250;
    await settle();
    assert.equal(requests.at(-1).frame, expected);
    assert.equal(requests.at(-1).preview, true);
    const waiting = requests.length;
    await settle();
    assert.equal(requests.length, waiting);
    reply(requests.at(-1));
    await tick();
}
assert.equal(playbackFrame([{label: 0}, {label: 1}, {label: 9}, {label: 10}], 4500), 3);
assert.equal(playbackFrame(template.layout.sliders[0].steps, 5000), 1);
await graph.handlers.plotly_buttonclicked(playback);
now += 1250;
assert.equal(requests.at(-1).preview, false);
reply(requests.at(-1));
await tick();
const paused = requests.length;
await settle();
assert.equal(requests.length, paused);
assert.equal(graph.layout.updatemenus[0].buttons[0].label, '▶');

// Manual scrubbing pauses playback; an in-flight frame cannot overwrite the selected time.
await graph.handlers.plotly_buttonclicked(playback);
await settle();
const playingRequest = requests.at(-1);
scrub(4);
reply(playingRequest);
await settle();
const stoppedPreview = requests.at(-1);
reply(stoppedPreview);
await tick();
assert.equal(requests.at(-1).frame, 4);
assert.equal(requests.at(-1).preview, false);
reply(requests.at(-1));
await tick();
assert.equal(graph.layout.updatemenus[0].buttons[0].label, '▶');
assert.equal(graph.layout.sliders[0].active, 3);

// Camera gestures retain native control during playback; release requests the current time and viewport.
await graph.handlers.plotly_buttonclicked(playback);
graph.handlers.pointerdown({target:{closest:()=>true}});
const beforeGesture = requests.length;
now += 2800;
graph._fullLayout.xaxis.range = [3, 4];
graph.handlers.plotly_relayout({'xaxis.range':[3, 4]});
await settle();
assert.equal(requests.length, beforeGesture);
document.handlers.pointerup();
await tick();
assert.equal(requests.at(-1).projection[12], -7);
assert.equal(requests.at(-1).preview, true);
reply(requests.at(-1));
await tick();
await graph.handlers.plotly_buttonclicked(playback);
reply(requests.at(-1));
await tick();

// Collapsing the sidebar gives its freed width to the existing plot without resetting the time.
const selectedFrame = graph.layout.sliders[0].active;
root.clientWidth += 128;
observers[0]();
await tick();
assert.equal(graph.layout.width, root.clientWidth);
assert.equal(graph.layout.height, root.clientHeight);
assert.equal(graph.layout.sliders[0].active, selectedFrame);
reply(requests.at(-1));
await tick();

// The gl3d projection composes data scaling, model translation, view and projection in order.
const identity = [1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1];
const model = [...identity]; model[12] = 7;
const perspective = [...identity]; perspective[0] = 3;
graph._fullLayout.scene = {_scene: {dataScale: [2, 4, 8], getCamera() {}, camera:{matrix:identity},
    glplot: {cameraParams: {model, view: identity, projection: perspective}}}};
const matrix = projection(graph, [1, 2, 3]);
assert.deepEqual([matrix[0], matrix[5], matrix[10], matrix[12]], [6, 4, 8, 21]);

// A native 3D gesture can advance the GL camera before Plotly saves gd.layout; data must not rewind it.
const sceneRoot = new Element(), sceneRequests = [];
let sceneRespond;
const liveCamera = {eye:{x:2, y:-1, z:0.5}};
await mount(sceneRoot, {...Plotly, async newPlot(g, data, layout) {
    await Plotly.newPlot(g, data, layout);
    g._fullLayout.scene = {_scene:{getCamera:()=>liveCamera, camera:{matrix:identity}, dataScale:[1,1,1],
        glplot:{cameraParams:{model:identity,view:identity,projection:identity}}}};
}}, {notify:q=>sceneRequests.push(q)},
{value:JSON.stringify(initial),on:callback=>{sceneRespond=callback;}},
JSON.stringify({...template,dimensions:[1,2,3],layout:{...template.layout,scene:{camera:{eye:{x:1,y:1,z:1}}}}}));
sceneRespond(JSON.stringify({...initial,id:sceneRequests.at(-1).id,data:[{...initial.data[0],x:[2]}]}));
await tick();
assert.deepEqual(sceneRoot.children[0].layout['scene.camera'], liveCamera);
sceneRoot.isConnected = false;
observers.at(-1)();

// Removing a plot cancels queued work and ignores late backend responses.
root.isConnected = false;
observers[1]();
assert.equal(respond(JSON.stringify(initial)), false);
assert.equal(purged, true);
const detached = new Element();
detached.isConnected = false;
await mount(detached, Plotly, {}, {}, '');
assert.equal(detached.children.length, 0);
console.log('PASS: spatial requests, playback, pause, looping, viewport and disposal');
