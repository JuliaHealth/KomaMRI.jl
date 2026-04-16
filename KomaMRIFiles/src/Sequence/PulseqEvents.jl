abstract type GradEvent end
struct TrapGradEvent <: GradEvent
    amplitude::Float64
    rise::Float64
    flat::Float64
    fall::Float64
    delay::Float64
end

mutable struct ArbGradEvent <: GradEvent
    amplitude::Float64
    first::Float64
    last::Float64
    amp_shape_id::Int
    time_shape_id::Int
    delay::Float64
end

struct RFEvent
    amplitude::Float64
    mag_id::Int
    phase_id::Int
    time_shape_id::Int
    center::Union{Nothing, Float64}
    delay::Float64
    freq_ppm::Float64
    phase_ppm::Float64
    freq::Float64
    phase::Float64
    use::Char
end

struct ADCEvent
    num::Int
    dwell::Float64
    delay::Float64
    freq_ppm::Float64
    phase_ppm::Float64
    freq::Float64
    phase::Float64
    phase_id::Int
end

struct ExtensionInstanceEvent
    type::Int
    ref::Int
    next_id::Int
end

const ShapeLibrary = Dict{Int, Tuple{Int, Vector{Float64}}}

struct PulseqEventLibraries
    grad_library::Dict{Int, GradEvent}
    rf_library::Dict{Int, RFEvent}
    adc_library::Dict{Int, ADCEvent}
    tmp_delay_library::Dict{Int, Float64}
    shape_library::ShapeLibrary
    extension_instance_library::Dict{Int, ExtensionInstanceEvent}
    extension_type_library::Dict{Int, Type{<:Extension}}
    extension_spec_library::Dict{Int, Dict{Int, Extension}}
    definitions::Dict{String, Any}
    block_duration_raster::Float64
    gradient_raster_time::Float64
    radiofrequency_raster_time::Float64
    adc_raster_time::Float64
end

function PulseqEventLibraries(
    grad_library,
    rf_library,
    adc_library,
    tmp_delay_library,
    shape_library,
    extension_instance_library,
    extension_type_library,
    extension_spec_library,
    definitions,
)
    defs = Dict{String, Any}(definitions)
    return PulseqEventLibraries(
        grad_library,
        rf_library,
        adc_library,
        tmp_delay_library,
        shape_library,
        extension_instance_library,
        extension_type_library,
        extension_spec_library,
        defs,
        defs["BlockDurationRaster"],
        defs["GradientRasterTime"],
        defs["RadiofrequencyRasterTime"],
        defs["AdcRasterTime"],
    )
end

struct PulseqBlock
    id::Int
    duration::Int
    rf::Int
    gx::Int
    gy::Int 
    gz::Int
    adc::Int
    ext::Int 
end

struct PulseqRaster
    BlockDurationRaster::Float64
    GradientRasterTime::Float64
    RadiofrequencyRasterTime::Float64
    AdcRasterTime::Float64
    PulseqRaster(seq, sys) = new(
        get_raster_time("BlockDurationRaster", seq, sys.seq_Δt),
        get_raster_time("GradientRasterTime", seq, sys.GR_Δt),
        get_raster_time("RadiofrequencyRasterTime", seq, sys.RF_Δt),
        get_raster_time("AdcRasterTime", seq, sys.ADC_Δt)
    )
end

function get_raster_time(key::String, seq::Sequence, scanner_default)
    haskey(seq.DEF, key) || return scanner_default
    seq.DEF[key] == scanner_default || @warn "Sequence and Scanner definition for $key do not match (($(seq.DEF[key]) != $(scanner_default))). Using the Scanner definition ($key = $scanner_default)."
    return scanner_default
end

function merge_definitions_with_raster!(definitions::AbstractDict{String,Any}, raster::PulseqRaster)
    definitions["BlockDurationRaster"] = raster.BlockDurationRaster
    definitions["GradientRasterTime"] = raster.GradientRasterTime
    definitions["RadiofrequencyRasterTime"] = raster.RadiofrequencyRasterTime
    definitions["AdcRasterTime"] = raster.AdcRasterTime
    return definitions
end

TrapGradEvent(data) = TrapGradEvent(data...)
ArbGradEvent(data::NTuple{3, T}) where {T<:Real} = ArbGradEvent(data[1], 0.0, 0.0, Int(data[2]), 0, data[3])
ArbGradEvent(data::NTuple{4, T}) where {T<:Real} = ArbGradEvent(data[1], 0.0, 0.0, Int(data[2]), Int(data[3]), data[4])
ArbGradEvent(data::NTuple{6, T}) where {T<:Real} = ArbGradEvent(data[1], data[2], data[3], Int(data[4]), Int(data[5]), data[6])

RFEvent(data::NTuple{6, T}) where {T<:Real} = RFEvent(data[1], Int(data[2]), Int(data[3]), 0, nothing, data[4], 0.0, 0.0, data[5], data[6], 'u')
RFEvent(data::NTuple{7, T}) where {T<:Real} = RFEvent(data[1], Int(data[2]), Int(data[3]), Int(data[4]), nothing, data[5], 0.0, 0.0, data[6], data[7], 'u')
RFEvent(data::Tuple{Vararg{Any, 11}}) = RFEvent(data[1], Int(data[2]), Int(data[3]), Int(data[4]), data[5], data[6], data[7], data[8], data[9], data[10], data[11])

ADCEvent(data::NTuple{5, T}) where {T<:Real} = ADCEvent(Int(data[1]), data[2], data[3], 0.0, 0.0, data[4], data[5], 0)
ADCEvent(data::NTuple{8, T}) where {T<:Real} = ADCEvent(Int(data[1]), data[2], data[3], data[4], data[5], data[6], data[7], Int(data[8]))

ExtensionInstanceEvent(data) = ExtensionInstanceEvent(Int(data[1]), Int(data[2]), Int(data[3]))
