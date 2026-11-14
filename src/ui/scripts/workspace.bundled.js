// deno-fmt-ignore-file
// deno-lint-ignore-file
// This code was bundled using `deno bundle` and it's not recommended to edit it manually

function mount(entries, add, closeWindow, darkmode) {
    const grid = document.getElementById('workspace-grid');
    const columns = document.getElementById('workspace-columns');
    const addButton = document.getElementById('workspace-add');
    const theme = document.getElementById('workspace-theme');
    const workspace = document.querySelector('.koma-workspace');
    const links = document.querySelector('.workspace-links');
    const tiles = new Map();
    let focused = null;
    let connected = false;
    const expandIcon = '<svg viewBox="0 0 24 24" aria-hidden="true"><path d="M8 3H3v5m13-5h5v5M3 16v5h5m13-5v5h-5"/></svg>';
    const restoreIcon = '<svg viewBox="0 0 24 24" aria-hidden="true"><path d="M8 8V3h13v13h-5M3 8h13v13H3z"/></svg>';
    columns.onchange = ()=>grid.style.setProperty('--columns', columns.value);
    addButton.onclick = ()=>add.notify(1);
    theme.onchange = ()=>darkmode.notify(theme.checked);
    const applyTheme = (dark)=>{
        workspace.dataset.theme = dark ? 'dark' : 'light';
        theme.checked = dark;
        theme.parentElement.title = dark ? 'Switch to light mode' : 'Switch to dark mode';
        theme.nextElementSibling.className = dark ? 'bi bi-moon-stars-fill' : 'bi bi-sun-fill';
    };
    darkmode.on(applyTheme);
    applyTheme(darkmode.value);
    const focus = (id)=>{
        focused = focused === id ? null : id;
        grid.classList.toggle('is-focused', focused !== null);
        for (const [number, tile] of tiles){
            tile.element.classList.toggle('is-focused', focused === number);
            tile.expand.innerHTML = focused === number ? restoreIcon : expandIcon;
            tile.expand.title = focused === number ? 'Restore grid' : 'Expand window';
            tile.expand.setAttribute('aria-label', tile.expand.title);
            tile.expand.setAttribute('aria-expanded', String(focused === number));
        }
    };
    const update = (rows)=>{
        const ids = new Set(rows.map(([id])=>id));
        for (const [id, tile] of tiles){
            if (ids.has(id)) continue;
            if (focused === id) focus(id);
            tile.element.remove();
            tiles.delete(id);
        }
        for (const [id, url] of rows){
            if (tiles.has(id)) continue;
            const element = document.createElement('section');
            element.className = 'workspace-tile';
            element.dataset.window = id;
            const home = document.createElement('button');
            home.type = 'button';
            home.className = 'koma-titlebar-home';
            home.title = 'Home';
            home.setAttribute('aria-label', 'Home');
            const number = document.createElement('strong');
            number.textContent = '#' + id;
            home.append(number);
            const status = document.createElement('span');
            status.className = 'koma-titlebar-status';
            status.setAttribute('role', 'status');
            const expand = document.createElement('button');
            expand.type = 'button';
            expand.className = 'koma-titlebar-expand';
            expand.innerHTML = expandIcon;
            expand.title = 'Expand window';
            expand.setAttribute('aria-label', expand.title);
            expand.setAttribute('aria-expanded', 'false');
            expand.onclick = ()=>focus(id);
            const closeButton = document.createElement('button');
            closeButton.type = 'button';
            closeButton.className = 'koma-titlebar-close';
            closeButton.innerHTML = '<svg viewBox="0 0 24 24" aria-hidden="true"><path d="M6 6l12 12M18 6L6 18"/></svg>';
            closeButton.title = 'Close window';
            closeButton.setAttribute('aria-label', 'Close window #' + id);
            closeButton.disabled = true;
            closeButton.onclick = ()=>closeWindow.notify(id);
            const frame = document.createElement('iframe');
            frame.dataset.komaWorkspace = 'true';
            frame.src = url;
            frame.title = 'KomaUI workspace window #' + id;
            home.onclick = ()=>frame.contentWindow.KomaUI?.notify('home');
            const identity = document.createElement('div');
            identity.className = 'koma-titlebar-title';
            identity.append(home);
            const controls = document.createElement('div');
            controls.className = 'koma-titlebar-actions';
            controls.append(status, expand, closeButton);
            frame.addEventListener('load', ()=>{
                frame.contentDocument.querySelector('.koma-window-titlebar').replaceChildren(identity, controls);
                frame.contentWindow.dispatchEvent(new frame.contentWindow.Event('resize'));
            });
            element.append(frame);
            tiles.set(id, {
                element,
                home,
                status,
                expand,
                closeButton,
                frame,
                busy: true
            });
            grid.append(element);
        }
    };
    document.addEventListener('click', (event)=>{
        if (!links.contains(event.target)) links.open = false;
    });
    window.addEventListener('message', (event)=>{
        if (event.origin !== window.location.origin || event.data?.type !== 'koma-workspace-status') return;
        const tile = [
            ...tiles.values()
        ].find((tile)=>tile.frame.contentWindow === event.source);
        if (!tile) return;
        const data = event.data;
        tile.busy = data.busy;
        tile.closeButton.disabled = !connected || data.busy;
        const disconnected = data.connection && data.connection !== 'connected';
        tile.home.disabled = data.busy || disconnected;
        tile.status.dataset.state = disconnected ? 'disconnected' : data.state;
        tile.status.textContent = disconnected ? 'Disconnected' : data.state === 'failed' ? 'Failed' : '';
        tile.status.title = tile.status.textContent;
    });
    entries.on(update);
    update(entries.value);
    Bonito.register_connection_indicator({
        onStatusChange (connection) {
            connected = connection === 'connected';
            for (const tile of tiles.values())tile.closeButton.disabled = !connected || tile.busy;
            addButton.disabled = connection !== 'connected';
            theme.disabled = connection !== 'connected';
            document.getElementById('workspace-connection').textContent = connection === 'connected' ? '' : 'Disconnected';
        }
    });
}
export { mount as mount };
