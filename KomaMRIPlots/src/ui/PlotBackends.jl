# Available backends for plots
const PLOT_BACKENDS = ("PlotlyBase", "PlutoPlotly")
const PLOT_BACKEND = Ref{String}("PlotlyBase")

# Choose plot backend
function plot_backend!(backend::String)
    if backend == PLOT_BACKEND[]
        @info """
        KomaMRIPlots backend is already set to: $backend.
        No need to do anything else.
        """
        return
    end
    backend in PLOT_BACKENDS || throw(ArgumentError("""
    Unsupported KomaMRIPlots backend: $backend.
    Supported backends are: $PLOT_BACKENDS.
    """))
    PLOT_BACKEND[] = backend
    @info """
    New KomaMRIPlots backend set: $backend.
    """
end

# Backends
function plot_koma(args...; kwargs...)
    if PLOT_BACKEND[] == "PlotlyBase"
        plot_koma(KomaPlotlyBaseBackend(), args...; kwargs...)
    elseif PLOT_BACKEND[] == "PlutoPlotly"
        plot_koma(KomaPlutoPlotlyBackend(), args...; kwargs...)
    else
        error("""
        Unsupported KomaMRIPlots backend: $PLOT_BACKEND.
        Supported backends are: $PLOT_BACKENDS.
        """)
    end
end

# PlotlyBase
struct KomaPlotlyBaseBackend end
function plot_koma(::KomaPlotlyBaseBackend, args...; kwargs...)
    return PlotlyBase.Plot(args...; kwargs...)
end

# PlutoPlotly
struct KomaPlutoPlotlyBackend end
const PLUTOPLOTLY_LOADED = Ref{Bool}(false)
function plot_koma(::KomaPlutoPlotlyBackend, args...; kwargs...)
    if PLUTOPLOTLY_LOADED[]
        return _plutoplotly_plot(args...; kwargs...)
    else
        @info """
        The PlutoPlotly functionality is being called but
        `PlutoPlotly.jl` must be loaded to access it.
        Add `using PlutoPlotly` or `import PlutoPlotly` to your code.
        """ maxlog=1
        return nothing
    end
end
function _plutoplotly_plot end
