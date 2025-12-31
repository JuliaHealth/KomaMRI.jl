mutable struct Trigger <: Extension 
    type::Int # Type of trigger (system dependent). 0: undefined / unused
    channel::Int # channel of trigger (system dependent). 0: undefined / unused
    d1::Float64 # Delay prior to the trigger event [s]
    d2::Float64 # Duration of trigger event [s]
end

get_scale(::Type{Trigger}) = [1 1 1e-6 1e-6]