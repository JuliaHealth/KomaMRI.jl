include("ui/BonitoPlotly.jl")
include("ui/WindowUI.jl")
include("ui/FilePickerUI.jl")
include("ui/ViewUI.jl")
include("ui/ExportUI.jl")
include("ui/SimulationUI.jl")
include("ui/UIActions.jl")
include("ui/BonitoUI.jl")

"""
    KomaUI(; kwargs...)

Open Koma's desktop UI.

Set `return_window=true` to return the `KomaWindow`, and `show_window=false` to
build the UI without opening its window.

Each window owns its observables: `w.seq`, `w.obj`, `w.sys`, `w.physio`, `w.raw`,
`w.img`, `w.sim_params`, and `w.rec_params`. Assign values with `[]` (for example,
`w.seq[] = seq`) to update that window's data and view.
"""
function KomaUI(; kwargs...)
    return launch_ui(; kwargs...)
end
