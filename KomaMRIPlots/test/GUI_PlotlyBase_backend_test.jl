@testitem "PlotlyBase backend" tags=[:plots] begin
    using KomaMRIBase, MRIFiles, PlotlyBase
    import KomaMRIPlots: plot_backend!

    plot_backend!("PlotlyBase")

    function plot_json(plot)
        @test plot isa Plot
        @test !isempty(plot.data)
        json = sprint(show, MIME"application/vnd.plotly.v1+json"(), plot)
        @test startswith(json, '{')
        @test !occursin(r"\b(?:NaN|Inf)\b", json)
        return json
    end

    @testset "Phantom" begin
        phantom = brain_phantom2D()
        density = plot_phantom_map(phantom, :ρ; width=800, height=600)
        plot_json(density)
        @test density.layout[:width] == 800
        @test density.layout[:height] == 600

        phantom.motion = MotionList(
            translate(0.1, 0.1, 0.1, TimeRange(0.0, 1.0), SpinRange(1:1000)),
            rotate(0.0, 0.0, 90.0, TimeRange(; t_start=0.05, t_end=0.5), SpinRange(1:1000)),
        )
        plot_json(plot_phantom_map(phantom, :T1; max_spins=1_000))
        plot_json(plot_phantom_map(phantom, :ρ; view_2d=true, max_spins=1_000))
    end

    @testset "Sequence" begin
        sys = Scanner()
        excitation = PulseDesigner.RF_hard(sys.B1, π / 2 / (2π * γ * sys.B1), sys; G=[0, 0, 0])
        sequence = excitation + PulseDesigner.EPI(23e-2, 101, sys)

        waveform = plot_seq(sequence; slider=true, show_seq_blocks=true)
        plot_json(waveform)
        trace_names = [get(trace, :name, nothing) for trace in waveform.data]
        @test all(in(trace_names), ("Gx", "Gy", "Gz", "ADC"))

        labeled = Sequence()
        @addblock labeled += (ADC(1, 1e-6), LabelSet(0, "LIN"), LabelSet(1, "REV"))
        @addblock labeled += (ADC(1, 1e-6), LabelSet(3, "LIN"), LabelSet(0, "REV"))
        labels = plot_seq(labeled)
        plot_json(labels)
        trace_names = [get(trace, :name, nothing) for trace in labels.data]
        @test all(in(trace_names), ("LIN", "REV"))
        @test labels.layout[:updatemenus][1][:yanchor] == "bottom"

        for plot in (plot_M0(sequence), plot_M1(sequence), plot_M2(sequence))
            plot_json(plot)
            @test isempty(plot.layout[:updatemenus])
        end

        plot_json(plot_kspace(sequence))
        plot_json(plot_eddy_currents(sequence, 80e-3))
        plot_json(plot_slew_rate(sequence))
        plot_json(plot_seqd(sequence))
    end

    @testset "Raw signal and table" begin
        raw = RawAcquisitionData(ISMRMRDFile(joinpath(@__DIR__, "test_files", "Koma_signal.mrd")))
        signal = plot_signal(raw; width=800, height=600)
        plot_json(signal)
        @test signal.layout[:width] == 800
        @test signal.layout[:height] == 600

        filename = tempname() * ".png"
        try
            @test savefig(signal, filename) == filename
            @test filesize(filename) > 0
        finally
            rm(filename; force=true)
        end

        table = plot_dict(Dict("B0" => 1.5, "Gmax" => 0.06))
        @test occursin("<table", table)
        @test occursin("B0", table)
        @test occursin("Gmax", table)
    end
end
