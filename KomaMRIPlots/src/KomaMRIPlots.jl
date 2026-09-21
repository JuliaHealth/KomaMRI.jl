module KomaMRIPlots

using KomaMRIBase
using MAT, Interpolations, PlotlyBase
import PlotlyKaleido
using QMRIColors
using Artifacts
import Bonito
using Bonito: @js_str

include("ui/DisplayFunctions.jl")
include("ui/LegacyDisplayFunctions.jl")
include("ui/TimePlotSampling.jl")
include("ui/SequencePlotSource.jl")
include("ui/TimePlots.jl")
include("ui/CoilSensitivities.jl")
include("ui/SpatialPlots.jl")

"""Save a Plotly figure, starting Kaleido when needed."""
function savefig(args...; kwargs...)
    PlotlyKaleido.start()
    return PlotlyKaleido.savefig(args...; kwargs...)
end

export plot_seq,
    plot_M0,
    plot_M1,
    plot_M2,
    plot_eddy_currents,
    plot_seqd,
    plot_slew_rate,
    plot_kspace,
    plot_phantom_map,
    plot_phantom,
    plot_coil_sens,
    get_coil_sens_fov,
    plot_signal,
    plot_image,
    plot_dict,
    savefig

end
