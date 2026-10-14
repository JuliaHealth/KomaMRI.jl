export function sliceAnnotations(model, counts, time) {
    return model.planes.flatMap((plane, panel) => {
        const {shown, available} = counts[panel];
        if (shown >= available) return [];
        return {name:'slice-sampling', text:`Downsampled · ${(100*shown/available).toFixed(2)}%`,
            xref:`${plane.refs[0]} domain`, yref:`${plane.refs[1]} domain`, x:1, y:0,
            xanchor:'right', yanchor:'bottom', xshift:-4, yshift:4,
            showarrow:false, borderpad:0, font:{size:10},
            hovertext:`${shown.toLocaleString()} / ${available.toLocaleString()} spins in this slab's viewport ` +
                `at t = ${time} ms. Zoom for more detail.`};
    });
}

export function mountSlices(graph, Plotly, onChange) {
    let model = graph.layout.meta?.phantom_slices ?? graph.layout.meta?.image_slices;
    if (!model) return null;
    let image = Boolean(model.shape);
    let center = [...model.center], thickness = model.thickness;
    let frame = 0, property = graph.layout.meta.phantom_property;
    let disposed = false, rendering = false, pending = false, updateFrame, playTimer, clickTimer, playing = false, playStarted;
    const positionSliders = () => graph.layout.sliders.filter(slider => slider.name?.startsWith('slice-'));
    const clamp = (value, axis) => Math.max(model.bounds[axis][0], Math.min(model.bounds[axis][1], image ? Math.round(value) : value));
    const nearest = (steps, value) => steps.reduce((best, step, i) =>
        Math.abs(Number(step.value) - value) < Math.abs(Number(steps[best].value) - value) ? i : best, 0);
    function applyLayout(layout) {
        for (const plane of model.planes) {
            const position = layout.annotations?.find(item => item.name === `slice-position-${plane.normal}`);
            if (position) position.text = image ? `${'xyz'[plane.normal-1]} = ${plane.normal===3 ? model.zlabels[center[2]] : center[plane.normal-1]}` :
                `${'xyz'[plane.normal-1]} = ${Math.round(1000*center[plane.normal-1])/100} mm`;
        }
        layout.shapes = model.planes.flatMap(plane => {
            const [x, y] = plane.refs, [horizontal, vertical] = plane.dimensions.map(axis => center[axis - 1]);
            const tickLength = 12, tick = {type:'line', line:{color:'#2a7fb8', width:1}, opacity:0.8};
            return (image ? [0] : [-thickness/2,0,thickness/2]).flatMap(offset => [0,1].flatMap(edge => {
                const start = edge ? -tickLength : 0, end = edge ? 0 : tickLength;
                const marker = {...tick, line:{...tick.line, dash:offset ? '3px,3px' : 'solid'}};
                return [{...marker, xref:x, x0:horizontal+offset, x1:horizontal+offset,
                    yref:`${y} domain`, ysizemode:'pixel', yanchor:edge, y0:start, y1:end},
                    {...marker, yref:y, y0:vertical+offset, y1:vertical+offset,
                    xref:`${x} domain`, xsizemode:'pixel', xanchor:edge, x0:start, x1:end}];
            }));
        });
        for (const slider of layout.sliders ?? []) if (slider.name?.startsWith('slice-'))
            slider.active = nearest(slider.steps, center[Number(slider.name.slice(6)) - 1]);
        const menu = layout.updatemenus?.find(item => item.name === 'slab-thickness');
        if (menu) menu.active = menu.buttons.findIndex(button => button.args[0] === thickness);
    }
    function query() {
        return {center:[...center], thickness, ranges:model.planes.map(plane =>
            plane.axes.map(axis => [...graph._fullLayout[axis].range])), sizes:model.planes.map(plane =>
            plane.axes.map(axis => graph._fullLayout[axis]._length))};
    }
    function localData() {
        if (image) return graph.data.map((trace,index) => {
            const plane=model.planes[index%model.planes.length], [a,b]=plane.dimensions.map(d=>d-1);
            const values=model.values[Math.floor(index/model.planes.length)], z=[], customdata=[];
            for (let row=0;row<model.shape[b];row++) {
                const line=[],points=[];
                for (let column=0;column<model.shape[a];column++) {
                    const point=[...center]; point[a]=column; point[b]=row;
                    line.push(values[point[0]+model.shape[0]*(point[1]+model.shape[1]*point[2])]);
                    points.push([point[0],point[1],model.zlabels[point[2]]]);
                }
                z.push(line); customdata.push(points);
            }
            return {z,customdata};
        });
        const current = model.frames[frame], coordinates = current.coordinates;
        const entry = new Set(current.entry), exit = new Set(current.exit);
        const velocity = ['speed','vx','vy','vz'].includes(property);
        const otherFrame = frame === 0 ? Math.min(1, model.frames.length - 1) : frame - 1;
        const other = model.frames[otherFrame].coordinates, dt = model.times[frame] - model.times[otherFrame];
        const reset = frame === 0 ? new Set(model.frames[otherFrame].entry) : entry;
        const value = i => {
            if (!velocity) return model.values[property][i];
            if (!dt || reset.has(i + 1)) return 0;
            if (property === 'speed') return Math.hypot(...coordinates.map((v, axis) => v[i] - other[axis][i])) / Math.abs(dt);
            const axis = ['vx','vy','vz'].indexOf(property);
            return (coordinates[axis][i] - other[axis][i]) / dt;
        };
        const data = model.planes.flatMap((plane, panel) => {
            const groups = [[],[],[]], normal = plane.normal - 1;
            for (let i = 0; i < coordinates[normal].length; i++) {
                if (Math.abs(coordinates[normal][i] - center[normal]) > thickness / 2) continue;
                const group = entry.has(i+1) ? 2 : exit.has(i+1) ? 1 : 0;
                groups[group].push(i);
            }
            return groups.map((ids, group) => {
                const trace = structuredClone(model.properties[property]);
                const values = ids.map(value);
                trace.uid = `phantom-slice-${panel+1}-${group+1}`;
                [trace.xaxis, trace.yaxis] = plane.refs;
                [trace.x, trace.y] = plane.dimensions.map(axis => ids.map(i => coordinates[axis-1][i]));
                trace.customdata = ids.map(i => i+1);
                trace.ids = trace.customdata.map(String);
                trace.text = values.map(v => Math.round(v * 1e4) / 1e4);
                trace.marker.color = values;
                trace.marker.showscale = false;
                trace.showlegend = velocity && panel === 2 && group > 0;
                const colorbar = trace.marker.colorbar;
                trace.hovertemplate = `Spin %{customdata}<br>${'xyz'[plane.dimensions[0]-1]}: %{x:.2f} cm<br>` +
                    `${'xyz'[plane.dimensions[1]-1]}: %{y:.2f} cm<br>` +
                    `<b>${colorbar.title}</b>: %{text}${colorbar.ticksuffix}<extra></extra>`;
                if (group > 0) {
                    trace.name = group === 1 ? 'Exit' : 'Re-entry';
                    if (velocity) {
                        trace.marker.color = group === 1 ? '#ff40df' : '#00e5ff';
                        trace.marker.size *= 2;
                    }
                    if (!ids.length) trace.x = trace.y = [null];
                }
                return trace;
            });
        });
        const marker = structuredClone(model.properties[property].marker);
        marker.colorbar.title = {text:marker.colorbar.title};
        marker.color = [marker.cmin, marker.cmax];
        marker.size = marker.opacity = 0;
        data.push({type:'scatter', uid:'phantom-slice-colorbar', x:[0,0], y:[0,0],
            mode:'markers', marker, showlegend:false, hoverinfo:'skip'});
        return data;
    }
    async function update() {
        if (disposed) return;
        if (!image && !onChange && model.frames.length > 1 && !graph.layout.updatemenus.some(menu => menu.name === 'motion-playback')) {
            const timeSlider = graph.layout.sliders.find(item => item.name === 'time');
            graph.layout.updatemenus.push({name:'motion-playback', type:'buttons', showactive:false, x:0, xanchor:'left',
                y:0, yanchor:'top', pad:{t:timeSlider.pad.t+4}, buttons:[{name:'motion-playback', label:'▶', method:'skip', args:[]}]});
        }
        const controls = () => {
            const layout = {sliders:graph.layout.sliders.map(control => ({...control})),
                updatemenus:(graph.layout.updatemenus ?? []).map(control => ({...control})),
                annotations:graph.layout.annotations.map(annotation => ({...annotation}))};
            applyLayout(layout);
            const changes = {shapes:layout.shapes, annotations:layout.annotations};
            if (image) {
                const size=graph._fullLayout._size, [x,y,z]=model.shape.map((n,i)=>n*model.spacing[i]);
                const gap=24,rowGap=60;
                const minimumPanelWidth=Math.min(128,(size.w-gap)/2), minimumPanelHeight=Math.min(100,(size.h-rowGap)/2);
                const scale=Math.min((size.w-gap)/(x+z),(size.h-rowGap)/(y+z),
                    (size.w-gap-minimumPanelWidth)/Math.max(x,z),
                    (size.h-rowGap-minimumPanelHeight)/Math.max(y,z));
                const widths=[x,z].map(value=>Math.max(minimumPanelWidth,value*scale));
                const heights=[y,z].map(value=>Math.max(minimumPanelHeight,value*scale));
                const left=(size.w-widths[0]-widths[1]-gap)/2, bottom=(size.h-heights[0]-heights[1]-rowGap)/2;
                const right=left+widths[0]+gap, top=bottom+heights[1]+rowGap;
                const panels=[[left,top,widths[0],heights[0]],[right,top,widths[1],heights[0]],[left,bottom,widths[0],heights[1]]];
                model.planes.forEach((plane,i)=>{
                    const [a,b,w,h]=panels[i], [imageWidth,imageHeight]=plane.dimensions.map(axis=>model.shape[axis-1]*model.spacing[axis-1]*scale);
                    const imageLeft=a+(w-imageWidth)/2, imageBottom=b+(h-imageHeight)/2;
                    changes[`${plane.axes[0]}.domain`]=[imageLeft/size.w,(imageLeft+imageWidth)/size.w];
                    changes[`${plane.axes[1]}.domain`]=[imageBottom/size.h,(imageBottom+imageHeight)/size.h];
                    changes[`sliders[${i}].x`]=a/size.w;
                    changes[`sliders[${i}].len`]=w/size.w;
                    changes[`sliders[${i}].y`]=imageBottom/size.h;
                    for (const [name,edge] of [['slice-title',0],['slice-position',1]]) {
                        const annotation=layout.annotations.find(item=>item.name===`${name}-${plane.normal}`);
                        if (annotation) Object.assign(annotation,{xref:'paper',yref:'paper',yanchor:'bottom',yshift:4,x:(a+edge*w)/size.w,y:(imageBottom+imageHeight)/size.h});
                    }
                });
                changes.imageColorbar={x:(right+widths[1]/2)/size.w,y:(bottom+heights[1]/2)/size.h,len:heights[1]/size.h};
            }
            for (const key of ['sliders','updatemenus']) layout[key]?.forEach((control, i) =>
                changes[`${key}[${i}].active`] = control.active);
            return changes;
        };
        if (onChange) {
            await Plotly.relayout(graph, controls());
            onChange();
            return;
        }
        if (rendering) { pending = true; return; }
        rendering = true;
        try {
            do {
                pending = false;
                const data = localData(), attributes = {};
                for (const key of image ? ['z','customdata'] : ['x','y','customdata','ids','text','name','hovertemplate','showlegend'])
                    attributes[key] = data.map(trace => trace[key]);
                for (const key of image ? [] : ['color','colorscale','colorbar','cmin','cmax','size','showscale'])
                    attributes[`marker.${key}`] = data.map(trace => trace.marker[key]);
                const changes=controls(), colorbar=changes.imageColorbar;
                delete changes.imageColorbar;
                if (colorbar) for (const [key,value] of Object.entries(colorbar))
                    attributes[`colorbar.${key}`]=data.map(()=>value);
                await Plotly.update(graph, attributes, changes, data.map((_, i) => i));
            } while (pending && !disposed);
        } finally { rendering = false; }
    }
    function changed() {
        cancelAnimationFrame(updateFrame);
        updateFrame = requestAnimationFrame(update);
    }
    function wheel(event) {
        if (event.ctrlKey || event.metaKey || !event.deltaY) return;
        const rect = graph.getBoundingClientRect(), x = event.clientX - rect.left, y = event.clientY - rect.top;
        const panel = model.planes.findIndex(plane => {
            const [a,b] = plane.axes.map(axis => graph._fullLayout[axis]);
            return a._offset <= x && x <= a._offset+a._length && b._offset <= y && y <= b._offset+b._length;
        });
        if (panel < 0) return;
        event.preventDefault(); event.stopImmediatePropagation();
        const axis = model.planes[panel].normal - 1, slider = positionSliders()[panel];
        const index = Math.max(0, Math.min(slider.steps.length-1, nearest(slider.steps, center[axis]) + Math.sign(event.deltaY)));
        center[axis] = Number(slider.steps[index].value);
        changed();
    }
    function click(event) {
        if (image && event.event?.detail > 1) return;
        const point = event.points?.[0];
        const plane = model.planes.find(item => item.refs[0] === point?.xaxis?._id);
        if (!plane) return;
        const select = () => {
            plane.dimensions.forEach((axis, row) => center[axis-1] = clamp(row === 0 ? point.x : point.y, axis-1));
            changed();
        };
        if (!image) return select();
        clearTimeout(clickTimer);
        clickTimer = setTimeout(select, graph._context.doubleClickDelay);
    }
    const doubleClick = () => clearTimeout(clickTimer);
    function slider(event) {
        if (event.slider.name?.startsWith('slice-')) {
            center[Number(event.slider.name.slice(6))-1] = Number(event.step.value);
            changed();
        } else if (!onChange && event.slider.name === 'time') { frame = Number(event.step.value)-1; changed(); }
    }
    function labelSliders() {
        graph.querySelectorAll('.slider-group').forEach(group => {
            const control = group.__data__;
            if (!control.name?.startsWith('slice-')) return;
            const handle = group.querySelector('.slider-grip-rect');
            handle.setAttribute('tabindex', '0');
            handle.setAttribute('role', 'slider');
            handle.setAttribute('aria-label', `${'xyz'[Number(control.name.slice(6))-1]} ${image ? 'image' : 'slice'} position`);
            handle.setAttribute('aria-orientation', 'horizontal');
            handle.setAttribute('aria-valuemin', control.steps[0].label);
            handle.setAttribute('aria-valuemax', control.steps.at(-1).label);
            handle.setAttribute('aria-valuenow', control.steps[control.active].label);
            handle.setAttribute('aria-valuetext', `${control.steps[control.active].label}${image ? '' : ' mm'}`);
        });
    }
    function focusSlider(event) {
        const group = event.target.closest('.slider-group');
        if (group?.__data__?.name?.startsWith('slice-'))
            group.querySelector('.slider-grip-rect').focus({preventScroll:true});
    }
    function key(event) {
        const control = event.target.closest('.slider-group')?.__data__;
        if (!control?.name?.startsWith('slice-') || event.altKey || event.ctrlKey || event.metaKey) return;
        if (event.key !== 'ArrowLeft' && event.key !== 'ArrowRight') return;
        event.preventDefault(); event.stopPropagation();
        const axis = Number(control.name.slice(6))-1;
        const index = Math.max(0, Math.min(control.steps.length-1,
            nearest(control.steps, center[axis]) + (event.key === 'ArrowRight' ? 1 : -1)));
        center[axis] = Number(control.steps[index].value);
        changed();
    }
    function button(event) {
        if (event.button.name === 'slab-thickness') { thickness = event.button.args[0]; changed(); }
        if (!onChange && event.button.name === 'phantom-property') {
            property = event.button.args[0];
            for (const menu of graph.layout.updatemenus) if (menu.name === 'phantom-property')
                menu.active = menu.buttons.findIndex(item => item.args[0] === property);
            changed();
        }
        if (!onChange && event.button.name === 'motion-playback') {
            playing = !playing;
            playStarted = performance.now() - 5000 * (model.times[frame]-model.times[0]) / (model.times.at(-1)-model.times[0]);
            clearTimeout(playTimer);
            const index = graph.layout.updatemenus.findIndex(menu => menu.name === 'motion-playback');
            Plotly.relayout(graph, {[`updatemenus[${index}].buttons[0].label`]:playing ? '❚❚' : '▶'});
            if (playing) play();
        }
    }
    function play() {
        if (!playing || disposed) return;
        const time = model.times[0] + ((performance.now()-playStarted)%5000)/5000 * (model.times.at(-1)-model.times[0]);
        frame = Math.max(0, model.times.findLastIndex(t => t <= time));
        const timeSlider = graph.layout.sliders.find(item => item.name === 'time');
        timeSlider.active = frame;
        changed();
        playTimer = setTimeout(play, 30);
    }
    graph.addEventListener('wheel', wheel, {capture:true, passive:false});
    graph.addEventListener('pointerdown', focusSlider, true);
    graph.addEventListener('keydown', key);
    graph.on('plotly_click', click); graph.on('plotly_sliderchange', slider); graph.on('plotly_buttonclicked', button);
    graph.on('plotly_doubleclick', doubleClick);
    graph.on('plotly_afterplot', labelSliders);
    function resized(event) { if (image && ['width','height','autosize'].some(key=>key in event)) changed(); }
    graph.on('plotly_relayout',resized);
    labelSliders();
    const removal = new MutationObserver(() => { if (!graph.isConnected) dispose(); });
    removal.observe(document, {childList:true, subtree:true});
    function dispose() {
        disposed = true; clearTimeout(playTimer); clearTimeout(clickTimer); cancelAnimationFrame(updateFrame); removal.disconnect();
        graph.removeEventListener('wheel', wheel, true);
        graph.removeEventListener('pointerdown', focusSlider, true);
        graph.removeEventListener('keydown', key);
        graph.removeListener?.('plotly_click', click); graph.removeListener?.('plotly_sliderchange', slider);
        graph.removeListener?.('plotly_buttonclicked', button);
        graph.removeListener?.('plotly_afterplot', labelSliders);
        graph.removeListener?.('plotly_relayout',resized);
        graph.removeListener?.('plotly_doubleclick', doubleClick);
    }
    if (!onChange) changed();
    return {query, applyLayout, dispose, async rebind(next) {
        playing = false;
        clearTimeout(playTimer);
        clearTimeout(clickTimer);
        if (JSON.stringify(next.center) !== JSON.stringify(model.center)) center = [...next.center];
        else center = center.map((value, axis) => Math.max(next.bounds[axis][0], Math.min(next.bounds[axis][1], value)));
        if (next.thickness !== model.thickness) thickness = next.thickness;
        model = next;
        image = Boolean(model.shape);
        if (!onChange) {
            property = (graph.layout.updatemenus ?? []).filter(menu => menu.name === 'phantom-property')
                .map(menu => menu.buttons[menu.active ?? 0]?.args[0]).find(value => value !== undefined) ??
                graph.layout.meta.phantom_property;
            frame = graph.layout.sliders.find(slider => slider.name === 'time')?.active ?? 0;
            await update();
        }
    }};
}
