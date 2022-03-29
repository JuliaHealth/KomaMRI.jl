@with_kw mutable struct Scanner
    #Main
    B0::Real=1.5      # Main magnetic field [T]
    B1::Real=10e-6     # Max RF amplitude [T]
    Gmax::Real=60e-3  # Max Gradient [T/m]
    Smax::Real=500    # Max Slew-rate [mT/m/ms] or [T/m/s]
    #Sampling
    ADC_Δt::Real=2e-6 # ADC raster time
    seq_Δt::Real=1e-5 # Seq-block raster time 
    GR_Δt::Real=1e-5  # GR raster time
    RF_Δt::Real=1e-6  # RF raster time 
    #Secondary
    RF_ring_down_T::Real=20e-6  # RF ring down time [s]
    RF_dead_time_T::Real=100e-6 # RF dead tim [s]
    ADC_dead_time_T::Real=10e-6 # ADC dead time [s]
end

#Functions that check that Sequence satisfy hardware requirements:
#check_sys_req(seq::Sequence,sys::Scanner)