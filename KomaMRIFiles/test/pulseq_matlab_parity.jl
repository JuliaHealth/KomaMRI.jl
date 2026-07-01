using KomaMRIBase
import KomaMRIBase: PulseDesigner as PD
using KomaMRIFiles
using Test
using Unitful

const MATLAB = get(ENV, "MATLAB", "matlab")
const PULSEQ_MATLAB = get(ENV, "PULSEQ_MATLAB", "")
const PULSEQ_MATLAB_SEQ = get(ENV, "PULSEQ_MATLAB_SEQ", "")
const MATLAB_PULSEQ_EVENT_SIGNIFICANT_DIGITS = 6
const MATLAB_PULSEQ_SHAPE_SIGNIFICANT_DIGITS = 9

matlab_string(s) = replace(s, "'" => "''")

function pulseq_block(event::Union{RF,ADC}, sys)
    seq = Sequence(sys)
    addblock!(seq, event)
    seq.DUR[end] = PD.ceil_to_raster(dur(seq[end], sys), sys.DUR_Δt)
    return seq
end

function run_pulseq_matlab_parity()
    sys = Scanner(
        B1=20u"μT", Gmax=40u"mT/m", Smax=170u"T/m/s", ADC_Δt=100u"ns",
        DUR_Δt=10u"μs", GR_Δt=10u"μs", RF_Δt=1u"μs",
        RF_ring_down_time=30u"μs", RF_dead_time=80u"μs",
        ADC_dead_time=12u"μs",
    )
    # Tiny-dwell RF edge cases use a low flip angle to stay inside the test scanner B1 limit.
    edge_rf_flip = 0.1u"deg"
    cases = [
        "trap_area_short_triangle" =>
            PD.build_trapezoid(:x; area=2e-8u"T*s/m", sys),
        "trap_area_short_maxgrad" =>
            PD.build_trapezoid(:x; area=30e-6u"T*s/m", sys),
        "trap_area_duration" =>
            PD.build_trapezoid(:x; area=8e-6u"T*s/m", duration=1u"ms", sys),
        "trap_area_duration_rise" =>
            PD.build_trapezoid(:x; area=8e-6u"T*s/m", duration=1u"ms", rise_time=120u"μs", sys),
        "trap_area_duration_rise_fall" =>
            PD.build_trapezoid(:x; area=8e-6u"T*s/m", duration=1u"ms", rise_time=100u"μs", fall_time=130u"μs", sys),
        "trap_flat_area_flat_time" =>
            PD.build_trapezoid(:x; flat_area=10e-6u"T*s/m", flat_time=0.8u"ms", sys),
        "trap_flat_area_flat_time_rise" =>
            PD.build_trapezoid(:x; flat_area=10e-6u"T*s/m", flat_time=0.8u"ms", rise_time=90u"μs", sys),
        "trap_flat_area_flat_time_rise_fall" =>
            PD.build_trapezoid(:x; flat_area=10e-6u"T*s/m", flat_time=0.8u"ms", rise_time=90u"μs", fall_time=110u"μs", sys),
        "trap_amplitude_flat_time" =>
            PD.build_trapezoid(:x; amplitude=20u"mT/m", flat_time=0.8u"ms", sys),
        "trap_amplitude_flat_time_rise" =>
            PD.build_trapezoid(:x; amplitude=20u"mT/m", flat_time=0.8u"ms", rise_time=120u"μs", sys),
        "trap_amplitude_flat_time_rise_fall" =>
            PD.build_trapezoid(:x; amplitude=20u"mT/m", flat_time=0.8u"ms", rise_time=120u"μs", fall_time=130u"μs", sys),
        "trap_amplitude_duration" =>
            PD.build_trapezoid(:x; amplitude=20u"mT/m", duration=1.2u"ms", sys),
        "trap_amplitude_duration_rise" =>
            PD.build_trapezoid(:x; amplitude=20u"mT/m", duration=1.2u"ms", rise_time=120u"μs", sys),
        "trap_amplitude_duration_rise_fall" =>
            PD.build_trapezoid(:x; amplitude=20u"mT/m", duration=1.2u"ms", rise_time=120u"μs", fall_time=130u"μs", sys),
        "trap_negative_area" =>
            PD.build_trapezoid(:x; area=-8e-6u"T*s/m", sys),
        "trap_delay" =>
            PD.build_trapezoid(:x; area=8e-6u"T*s/m", delay=40u"μs", sys),
        "trap_overrides" =>
            PD.build_trapezoid(:x; flat_area=6e-6u"T*s/m", flat_time=0.8u"ms", max_grad=30u"mT/m", max_slew=120u"T/m/s", sys),
        "block_duration" =>
            PD.build_block_pulse(90u"deg"; duration=1u"ms", freq_offset=123u"Hz", phase_offset=0.3u"rad", sys, use=Excitation()),
        "block_bandwidth" =>
            PD.build_block_pulse(90u"deg"; bandwidth=250u"Hz", sys, use=Excitation()),
        "block_time_bw_product" =>
            PD.build_block_pulse(90u"deg"; bandwidth=500u"Hz", time_bw_product=3.0, sys, use=Excitation()),
        "block_deadtime" =>
            PD.build_block_pulse(90u"deg"; duration=1u"ms", delay=20u"μs", sys, use=Excitation()),
        "sinc_default" =>
            PD.build_sinc_pulse(90u"deg"; duration=2u"ms", sys, use=Excitation()),
        "sinc_custom" =>
            PD.build_sinc_pulse(90u"deg"; duration=2.4u"ms", freq_offset=123u"Hz", phase_offset=0.4u"rad", time_bw_product=6.0, apodization=0.46, center_pos=0.43, dwell=2u"μs", sys, use=Excitation()),
        "sinc_slice" =>
            PD.build_sinc_pulse(90u"deg"; duration=2u"ms", slice_thickness=5u"mm", sys, use=Excitation()),
        "sinc_slice_deadtime_overrides" =>
            PD.build_sinc_pulse(90u"deg"; duration=2u"ms", slice_thickness=5u"mm", delay=20u"μs", max_grad=30u"mT/m", max_slew=120u"T/m/s", sys, use=Excitation()),
        "arbitrary_rf" =>
            PD.build_arbitrary_rf([1, 2, 1], 45u"deg"; dwell=200u"μs", sys, use=Excitation()),
        "arbitrary_rf_custom" =>
            PD.build_arbitrary_rf([1, 1im, 1, -1im], 45u"deg"; dwell=100u"μs", center=150u"μs", freq_offset=77u"Hz", phase_offset=0.2u"rad", delay=20u"μs", sys, use=Excitation()),
        "rf_time_shape_default_raster" =>
            begin
                rf_samples = [1, 2, 1]
                rf0, _, _ = PD.make_arbitrary_rf(rf_samples, edge_rf_flip; dwell=sys.RF_Δt, sys, use=Excitation())
                rf_delay = sys.RF_dead_time + sys.RF_Δt / 2
                pulseq_block(
                    RF(
                        rf0.A, fill(sys.RF_Δt, length(rf_samples) - 1), rf0.Δf,
                        rf_delay; center=sys.RF_Δt, ϕ=rf0.ϕ, use=Excitation(),
                    ),
                    sys,
                )
            end,
        "rf_uniform_half_raster" =>
            PD.build_arbitrary_rf([1, 2, 1], edge_rf_flip; dwell=sys.RF_Δt / 2, sys, use=Excitation()),
        "rf_time_shape_start_offset" =>
            begin
                rf_samples = [1, 2, 1]
                rf0, _, _ = PD.make_arbitrary_rf(rf_samples, edge_rf_flip; dwell=sys.RF_Δt, sys, use=Excitation())
                pulseq_delay = sys.RF_dead_time
                time_shape_start = 2sys.RF_Δt
                time_shape_intervals = [3sys.RF_Δt, 4sys.RF_Δt]
                time_shape_center = time_shape_start + first(time_shape_intervals)
                pulseq_block(
                    RF(
                        rf0.A, time_shape_intervals, rf0.Δf,
                        pulseq_delay + time_shape_start;
                        center=time_shape_center - time_shape_start,
                        ϕ=rf0.ϕ, use=Excitation(),
                    ),
                    sys,
                )
            end,
        "rf_time_shape_irregular" =>
            begin
                rf_samples = [1, 2, 1]
                rf0, _, _ = PD.make_arbitrary_rf(rf_samples, edge_rf_flip; dwell=sys.RF_Δt, sys, use=Excitation())
                pulseq_delay = sys.RF_dead_time
                time_shape_start = sys.RF_Δt / 2
                time_shape_intervals = [3sys.RF_Δt / 2, 5sys.RF_Δt]
                time_shape_center = 2sys.RF_Δt
                pulseq_block(
                    RF(
                        rf0.A, time_shape_intervals, rf0.Δf,
                        pulseq_delay + time_shape_start;
                        center=time_shape_center - time_shape_start,
                        ϕ=rf0.ϕ, use=Excitation(),
                    ),
                    sys,
                )
            end,
        "arbitrary_rf_slice_bandwidth" =>
            PD.build_arbitrary_rf([1, 1, 1, 1], 90u"deg"; dwell=100u"μs", bandwidth=2u"kHz", slice_thickness=5u"mm", sys, use=Excitation()),
        "arbitrary_rf_slice_time_bw_product" =>
            PD.build_arbitrary_rf([1, 1, 1, 1], 90u"deg"; dwell=100u"μs", time_bw_product=2.0, slice_thickness=5u"mm", sys, use=Excitation()),
        "gauss_default" =>
            PD.build_gauss_pulse(90u"deg"; duration=2u"ms", sys, use=Excitation()),
        "gauss_bandwidth" =>
            PD.build_gauss_pulse(90u"deg"; duration=2u"ms", bandwidth=2u"kHz", sys, use=Excitation()),
        "gauss_custom" =>
            PD.build_gauss_pulse(90u"deg"; duration=2.4u"ms", bandwidth=1.8u"kHz", apodization=0.2, center_pos=0.43, dwell=2u"μs", freq_offset=123u"Hz", phase_offset=0.4u"rad", delay=20u"μs", sys, use=Excitation()),
        "gauss_slice" =>
            PD.build_gauss_pulse(90u"deg"; duration=2u"ms", bandwidth=2u"kHz", slice_thickness=5u"mm", sys, use=Excitation()),
        "gauss_slice_overrides" =>
            PD.build_gauss_pulse(90u"deg"; duration=2u"ms", bandwidth=2u"kHz", slice_thickness=5u"mm", max_grad=30u"mT/m", max_slew=120u"T/m/s", sys, use=Excitation()),
        "adc_positional_dwell" =>
            PD.build_adc(16, 2u"μs"; sys),
        "adc_dwell" =>
            PD.build_adc(16; dwell=2u"μs", freq_offset=77u"Hz", phase_offset=0.2u"rad", sys),
        "adc_duration" =>
            PD.build_adc(16; duration=32u"μs", sys),
        "adc_deadtime" =>
            PD.build_adc(16; dwell=2u"μs", delay=4u"μs", sys),
        "adc_one_sample" =>
            PD.build_adc(1; duration=3u"μs", sys),
        "delay_zero" =>
            PD.build_delay(0u"s"; sys),
        "delay_positive" =>
            PD.build_delay(1.5u"ms"; sys),
        "label_set_lin" =>
            PD.build_label(:SET, :LIN, 3; sys),
        "label_inc_slc" =>
            PD.build_label(:INC, :SLC, -2; sys),
        "label_set_ref" =>
            PD.build_label(:SET, :REF, true; sys),
        "rotation_phi" =>
            PD.build_rotation(60u"deg"; sys),
        "rotation_phi_theta_negative" =>
            PD.build_rotation(-45u"deg", 30u"deg"; sys),
        "rotation_axis" =>
            PD.build_rotation([1, 2, 3], 36u"deg"; sys),
        "rotation_negative_axis_angle" =>
            PD.build_rotation([-1, 0, 0], -30u"deg"; sys),
        "rotation_quaternion" =>
            PD.build_rotation([0.5, 0.5, -0.5, 0.5]; sys),
        "rotation_matrix" =>
            PD.build_rotation([0 -1 0; 1 0 0; 0 0 1]; sys),
        "trigger_physio1" =>
            PD.build_trigger(:physio1; sys),
        "trigger_physio2_delay" =>
            PD.build_trigger(:physio2; delay=20u"μs", duration=100u"μs", sys),
        "digital_osc0" =>
            PD.build_digital_output_pulse(:osc0; sys),
        "digital_osc1_delay" =>
            PD.build_digital_output_pulse(:osc1; delay=20u"μs", duration=100u"μs", sys),
        "digital_ext1" =>
            PD.build_digital_output_pulse(:ext1; duration=200u"μs", sys),
        "arbitrary_default" =>
            PD.build_arbitrary_grad(:x, [0, 1, 0] .* u"mT/m"; sys),
        "arbitrary_first_last_delay" =>
            PD.build_arbitrary_grad(:x, [0, 0.5, 1, 0.5, 0] .* u"mT/m"; first=-0.25u"mT/m", last=-0.25u"mT/m", delay=20u"μs", sys),
        "arbitrary_oversampling" =>
            PD.build_arbitrary_grad(:x, [0, 0.2, 0.4, 0.2, 0] .* u"mT/m"; oversampling=true, first=-0.2u"mT/m", last=-0.2u"mT/m", sys),
        "extended_raster" =>
            PD.build_extended_trapezoid(:x, [0, 0.5, 1] .* u"ms", [0, 1, 0] .* u"mT/m"; sys),
        "extended_delay_zero_start" =>
            PD.build_extended_trapezoid(:x, [20, 500, 1000] .* u"μs", [0, 1, 0] .* u"mT/m"; sys),
        "extended_skip_check" =>
            PD.build_extended_trapezoid(:x, [20, 500, 1000] .* u"μs", [1, 2, 0] .* u"mT/m"; skip_check=true, sys),
        "extended_convert2arbitrary" =>
            PD.build_extended_trapezoid(:x, [0, 0.355, 1] .* u"ms", [0, 1, 0] .* u"mT/m"; convert2arbitrary=true, sys),
        "extended_area_zero_edges" =>
            PD.build_extended_trapezoid_area(:x, 0u"mT/m", 0u"mT/m", 10e-6u"T*s/m"; sys),
        "extended_area_nonzero_edges" =>
            begin
                seq = PD.build_extended_trapezoid(:x, [0, 100] .* u"μs", [0, 10] .* u"mT/m"; sys)
                append!(seq, PD.build_extended_trapezoid_area(:x, 10u"mT/m", 5u"mT/m", 2e-6u"T*s/m"; sys))
                seq
            end,
        "extended_area_negative" =>
            PD.build_extended_trapezoid_area(:x, 0u"mT/m", 0u"mT/m", -5e-6u"T*s/m"; sys),
    ]
    koma_seq = Sequence(sys)
    ranges = Pair{String,UnitRange{Int}}[]
    block = 1
    for (name, case) in cases
        append!(koma_seq, case)
        range = block:(block + length(case) - 1)
        push!(ranges, name => range)
        block = last(range) + 1
    end
    mktempdir() do dir
        koma_path = joinpath(dir, "koma_pulseq_parity.seq")
        matlab_path = isempty(PULSEQ_MATLAB_SEQ) ?
            joinpath(dir, "pulseq_matlab_parity.seq") : PULSEQ_MATLAB_SEQ
        if isempty(PULSEQ_MATLAB_SEQ)
            matlab = Sys.which(MATLAB)
            matlab === nothing && isfile(MATLAB) && (matlab = MATLAB)
            matlab === nothing && error("MATLAB not found. Set MATLAB or put matlab on PATH.")
            isempty(PULSEQ_MATLAB) && error("Set PULSEQ_MATLAB to the Pulseq MATLAB source directory.")
            isdir(PULSEQ_MATLAB) || error("Pulseq MATLAB path not found at $PULSEQ_MATLAB")
            script = joinpath(@__DIR__, "write_pulseq_matlab_parity_sequences.m")
            withenv(
                "PULSEQ_MATLAB_SEQ" => matlab_path,
                "KOMA_GAMMA" => string(γ),
                "PULSEQ_MATLAB" => PULSEQ_MATLAB,
            ) do
                run(`$matlab -batch $("run('$(matlab_string(script))')")`)
            end
        else
            isfile(matlab_path) || error("Pulseq MATLAB reference sequence not found at $matlab_path")
        end
        write_seq(
            koma_seq, koma_path; verbose=false,
            significant_digits=MATLAB_PULSEQ_EVENT_SIGNIFICANT_DIGITS,
            shape_significant_digits=MATLAB_PULSEQ_SHAPE_SIGNIFICANT_DIGITS,
        )
        koma_seq = read_seq(koma_path; verbose=false)
        matlab_seq = read_seq(matlab_path; verbose=false)
        @test length(koma_seq) == length(matlab_seq)
        for (name, range) in ranges
            @testset "$name" begin
                @test matlab_seq[range] ≈ koma_seq[range]
            end
        end
    end
end

@testset "MATLAB Pulseq parity" begin
    run_pulseq_matlab_parity()
end
