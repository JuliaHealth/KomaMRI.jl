module KomaPlotsPlutoPlotlyExt

    using KomaMRIPlots
    using PlutoPlotly

    function __init__()
        KomaMRIPlots.PLUTOPLOTLY_LOADED[] = true
        KomaMRIPlots.plot_backend!("PlutoPlotly")
    end

    # Define plot
    function KomaMRIPlots._plutoplotly_plot(args...; kwargs...)
        return PlutoPlotly.plot(args...; kwargs...)
    end

    # savefig
    function KomaMRIPlots.savefig(p::PlutoPlotly.PlutoPlot, fn::AbstractString; kwargs...)
        return savefig(p.Plot, fn; kwargs...)
    end

end
