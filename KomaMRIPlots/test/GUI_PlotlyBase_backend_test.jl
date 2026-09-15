@testitem "PlotlyBase backend" tags=[:plots] begin
    using KomaMRIBase, MRIFiles, PlotlyBase
    test_plot(plot) = @test plot isa Plot

    @testset "Image components" begin
        # Real images retain signed values; only complex images expose component controls.
        real_image = [-2 3; 4 -5]
        real_plot = plot_image(real_image)
        @test only(real_plot.data)[:z] == real_image
        @test isempty(get(real_plot.layout, :updatemenus, []))

        # Known phasors select magnitude or radians without applying magnitude limits to phase.
        image = ComplexF64[2 3im; -4 -5im]
        plot = plot_image(image; zmin=0, zmax=6)
        buttons = only(plot.layout[:updatemenus])[:buttons]
        magnitude, phase = [only(plot.data[button[:args][1][:visible]]) for button in buttons]
        @test magnitude[:z] == [2 3; 4 5]
        @test phase[:z] ≈ [0 π/2; π -π/2]
        @test (magnitude[:zmin], magnitude[:zmax]) == (0, 6)
        @test (phase[:zmin], phase[:zmax]) == (-π, π)
        @test first(phase[:colorscale])[2] == last(phase[:colorscale])[2]

        # A navigator uses the same component controls, with headroom only for magnitude.
        profile = plot_image(reshape(ComplexF64[2, -4], :, 1); zmin=0, zmax=6)
        buttons = only(profile.layout[:updatemenus])[:buttons]
        @test [trace[:y] for trace in profile.data] == [[2, 4], [0, π]]
        @test buttons[1][:args][2][:yaxis][:range] == [0, 1.1 * 6]
        @test buttons[2][:args][2][:yaxis][:range] == [-π, π]
    end

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

    @testset "Receive sensitivities" begin
        # Uniform reception stays unity and has no meaningless single-coil selector.
        uniform = plot_coil_sens(Scanner())
        @test all(==(1), only(uniform.data)[:marker][:color])
        @test isempty(uniform.layout[:sliders])

        # Known sampled phasors retain gain and phase in each physical plane.
        for plane in ((1, 2), (1, 3), (2, 3))
            coordinates = ntuple(d -> d in plane ? [-0.01, 0.0, 0.01] : [0.0], 3)
            values = ones(ComplexF64, length.(coordinates)..., 2)
            values[:, :, :, 1] .= 3im
            values[:, :, :, 2] .= -4
            plot = plot_coil_sens(Scanner(receiver=ArbitraryCoilSens(coordinates..., values)))
            selector = only(plot.layout[:sliders])
            @test [only(plot.data[step[:args][1][:visible]])[:marker][:color]
                for step in selector[:steps]] == [fill(3, 9), fill(4, 9)]
            buttons = only(plot.layout[:updatemenus])[:buttons]
            @test buttons[2][:args][1]["marker.color"] == [fill(π/2, 9), fill(π, 9)]
            @test plot.layout[:coloraxis][:cmax] == 4
            @test (plot.layout[:xaxis][:title], plot.layout[:yaxis][:title]) ==
                map(d -> ("x", "y", "z")[d], plane)
        end

        # A short birdcage changes the sampled extent, not the one-centimetre grid spacing.
        short = plot_coil_sens(Scanner(receiver=BirdcageCoilSens(L=0.01)))
        @test sort!(unique(first(short.data)[:z])) == [-1, 0, 1]
        @test all(isfinite, first(short.data)[:marker][:color])
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
        signal_plot = plot_signal(raw; width=800, height=600)
        test_plot(signal_plot)

        # Single-coil data has no coil selector.
        @test isempty(signal_plot.layout[:sliders])
        profile = first(raw.profiles)
        head = deepcopy(profile.head)
        head.available_channels = head.active_channels = head.number_of_samples = 2
        data = ComplexF32[3+4im -5-12im; 0 0]
        coils = Profile(head, profile.traj[:, 1:2], data)
        multi_coil_plot = plot_signal(RawAcquisitionData(raw.params, [coils]))

        # Each coil selects its known phasor components; the hidden coil also sets the range.
        selector = only(multi_coil_plot.layout[:sliders])
        @test [
            [collect(skipmissing(trace[:y])) for trace in multi_coil_plot.data[step[:args][1][:visible]]]
            for step in selector[:steps]
        ] == [[[5, 0], [3, 0], [4, 0]], [[13, 0], [-5, 0], [-12, 0]]]
        signal_limits = multi_coil_plot.layout[:yaxis][:range]
        @test first(signal_limits) <= -12 && last(signal_limits) >= 13

        @test plot_dict(Dict("B0" => 1.5, "Gmax" => 0.06)) isa AbstractString
    end
end
