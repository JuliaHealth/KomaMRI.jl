@with_kw mutable struct HardwareLimits
    #Main
    B0::Float64 = 1.5
    B1::Float64 = 10e-6
    Gmax::Float64 = 60e-3
    Smax::Float64 = 500.0
    #Sampling
    ADC_Δt::Float64 = 2e-6
    DUR_Δt::Float64 = 1e-5
    GR_Δt::Float64 = 1e-5
    RF_Δt::Float64 = 1e-6
    #Secondary
    RF_ring_down_time::Float64 = 20e-6
    RF_dead_time::Float64 = 100e-6
    ADC_dead_time::Float64 = 10e-6

    function HardwareLimits(
        B0, B1, Gmax, Smax, ADC_Δt, DUR_Δt, GR_Δt, RF_Δt,
        RF_ring_down_time, RF_dead_time, ADC_dead_time,
    )
        return new(
            scanner_field_strength(B0), scanner_rf_amplitude(B1),
            scanner_gradient_amplitude(Gmax), scanner_slew_rate(Smax),
            scanner_time(ADC_Δt), scanner_time(DUR_Δt), scanner_time(GR_Δt),
            scanner_time(RF_Δt), scanner_time(RF_ring_down_time),
            scanner_time(RF_dead_time), scanner_time(ADC_dead_time),
        )
    end
end

scanner_field_strength(x) = Float64(x)
scanner_rf_amplitude(x) = Float64(x)
scanner_gradient_amplitude(x) = Float64(x)
scanner_slew_rate(x) = Float64(x)
scanner_time(x) = Float64(x)