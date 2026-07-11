module KomaMRIPlots

using KomaMRIBase
using MAT, Interpolations, PlotlyBase
import PlotlyKaleido
using QMRIColors

include("ui/PlotBackends.jl")
include("ui/DisplayFunctions.jl")

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
    plot_signal,
    plot_image,
    plot_dict,
    savefig

end
