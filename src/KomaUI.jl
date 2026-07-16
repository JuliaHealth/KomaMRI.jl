include("ui/BonitoPlotly.jl")
include("ui/WindowUI.jl")
include("ui/FilePickerUI.jl")
include("ui/ViewUI.jl")
include("ui/ExportUI.jl")
include("ui/SimulationUI.jl")
include("ui/BonitoUI.jl")

"""
    KomaUI(; kwargs...)

Open Koma's desktop UI.

Set `return_window=true` to return the `KomaWindow`, and `show_window=false` to
build the UI without opening its window.
"""
function KomaUI(; kwargs...)
    return launch_ui(; kwargs...)
end
