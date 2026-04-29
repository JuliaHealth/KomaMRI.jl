struct Trigger <: Extension
    type::Int # Type of trigger (system dependent). 0: undefined / unused
    channel::Int # channel of trigger (system dependent). 0: undefined / unused
    delay::Float64 # Delay prior to the trigger event [s]
    duration::Float64 # Duration of trigger event [s]
end

dur(t::Trigger) = t.delay + t.duration

get_scale(::Type{Trigger}) = [1 1 1e-6 1e-6]
get_pulseq_format(::Type{Trigger}) = "%i %i %i %i"
get_EXT_type_from_symbol(::Val{:TRIGGERS}) = Trigger
get_symbol_from_EXT_type(::Type{Trigger})  = "TRIGGERS"
extension_type_header(::Type{Trigger}) = "# Extension specification for triggers:\n# id type channel d1 [us] d2 [us]\n"
