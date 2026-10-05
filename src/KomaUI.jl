include("ui/BonitoPlotly.jl")
include("ui/WindowUI.jl")
include("ui/FilePickerUI.jl")
include("ui/ViewUI.jl")
include("ui/ExportUI.jl")
include("ui/SimulationUI.jl")
include("ui/UIActions.jl")
include("ui/BonitoUI.jl")
include("ui/WorkspaceUI.jl")

"""
    KomaUI(; kwargs...)

Open Koma's desktop UI, or serve the same UI with `browser=true`.

Browser mode starts a server at `host="127.0.0.1", port=9384` and logs its URL;
it does not launch Electron or a browser on the server. Use `ui_status(w).url`
to retrieve the actual address. For a remote host, connect through an SSH tunnel.
The Julia process must remain running. The server is for one user's UI, with no
built-in authentication; keep the default loopback host when tunneling.

Dark mode is the default. Switch themes at the bottom of the sidebar or start
with `darkmode=false` for the light theme.

Use `KomaWorkspace` to show independent windows together in a browser grid.
Numbers belong to workspace tiles, not standalone windows.

Set `return_window=true` to return the `KomaWindow`, and `show_window=false` to
build the UI without opening a window or starting a server. `close(w)` stops the
window/server and releases its listeners; `wait(w)` waits for its lifetime, not
for an operation. `isopen(w)` reports whether the window/server is running.

Each window owns its observables: `w.seq`, `w.obj`, `w.sys`, `w.physio`, `w.raw`,
`w.img`, `w.sim_params`, and `w.rec_params`. Assign values with `[]` (for example,
`w.seq[] = seq`) to update that window's data and view.

Use `load_file!`, `simulate!`, and `reconstruct!` for synchronous operations with
returned results and propagated errors. `ui_status(w)` reports operation progress
and failures, and `ui_actions(w)` lists navigation/button actions for `click!`.

See the [browser-based UI guide](https://juliahealth.org/KomaMRI.jl/stable/how-to/1-1-use-koma-ui#Browser-based-KomaUI-easier-for-AI-agents)
and [Julia control examples](https://juliahealth.org/KomaMRI.jl/stable/how-to/1-1-use-koma-ui#Controlling-the-UI-from-Julia-easier-for-AI-agents).
"""
function KomaUI(; kwargs...)
    return launch_ui(; kwargs...)
end
