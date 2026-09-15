@testitem "PlotlyBase backend" tags=[:plots] begin
    using KomaMRIBase, MRIFiles, PlotlyBase
    test_plot(plot) = @test plot isa Plot

    @testset "Phantom" begin
        phantom = brain_phantom2D()
        for key in (:ρ, :T1, :T2, :x, :Δw)
            test_plot(plot_phantom_map(phantom, key))
        end
        test_plot(plot_phantom_map(phantom, :ρ; width=800, height=600, view_2d=true))

        phantom.motion = MotionList(
            translate(0.1, 0.1, 0.1, TimeRange(0.0, 1.0), SpinRange(1:1000)),
            rotate(
                0.0,
                0.0,
                90.0,
                TimeRange(; t_start=0.05, t_end=0.5),
                SpinRange(1:1000),
            ),
        )
        for key in (:ρ, :T1, :T2, :x, :Δw)
            test_plot(plot_phantom_map(phantom, key; max_spins=1_000))
        end
    end

    @testset "Sequence" begin
        sys = Scanner()
        excitation = PulseDesigner.RF_hard(
            sys.limits.B1, π / 2 / (2π * γ * sys.limits.B1), sys; G=[0, 0, 0]
        )
        sequence = excitation + PulseDesigner.EPI(23e-2, 101, sys)

        test_plot(plot_seq(sequence))
        test_plot(
            plot_seq(sequence; width=800, height=600, slider=true, show_seq_blocks=true)
        )

        labeled = Sequence()
        @addblock labeled += (ADC(1, 1e-6), LabelSet(0, "LIN"), LabelSet(1, "REV"))
        @addblock labeled += (ADC(1, 1e-6), LabelSet(3, "LIN"), LabelSet(0, "REV"))
        labels = plot_seq(labeled)
        test_plot(labels)
        trace_names = [get(trace, :name, nothing) for trace in labels.data]
        @test "LIN" in trace_names
        @test "REV" in trace_names

        triggered = Sequence()
        preceding_delay = Delay(1e-3)
        trigger = PulseDesigner.make_trigger(:physio2; delay=2e-3, duration=1e-3)
        @addblock triggered += preceding_delay
        @addblock triggered += trigger
        trigger_time = (dur(preceding_delay) + trigger.delay) * 1e3

        trigger_plot = plot_seq(triggered)
        @test only(trigger_plot.layout[:shapes])[:x0] ≈ trigger_time
        @test any(trace -> get(trace, :name, nothing) == "Trigger", trigger_plot.data)

        r_peak = 10e-3
        physio_plot = plot_seq(triggered; physio=CardiacSignal(; r_peaks=[r_peak]))
        @test any(trace -> get(trace, :name, nothing) == "ECG", physio_plot.data)
        @test only(physio_plot.layout[:shapes])[:x0] ≈ r_peak * 1e3

        for plot in (plot_M0(sequence), plot_M1(sequence), plot_M2(sequence))
            test_plot(plot)
        end

        test_plot(plot_kspace(sequence; width=800, height=600))
        test_plot(plot_eddy_currents(sequence, 80e-3))
        test_plot(plot_slew_rate(sequence))
        test_plot(plot_seqd(sequence))
    end

    @testset "Raw signal and table" begin
        raw = RawAcquisitionData(
            ISMRMRDFile(joinpath(@__DIR__, "test_files", "Koma_signal.mrd"))
        )
        test_plot(plot_signal(raw; width=800, height=600))

        @test plot_dict(Dict("B0" => 1.5, "Gmax" => 0.06)) isa AbstractString
    end
end
