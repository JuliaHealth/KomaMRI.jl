abstract type PulseqGradEvent end
struct PulseqTrapGradEvent <: PulseqGradEvent
    amplitude::Float64
    rise::Float64
    flat::Float64
    fall::Float64
    delay::Float64
end

struct PulseqArbGradEvent <: PulseqGradEvent
    amplitude::Float64
    first::Float64
    last::Float64
    amp_shape_id::Int
    time_shape_id::Int
    delay::Float64
end

struct PulseqRFEvent
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

struct PulseqADCEvent
    num::Int
    dwell::Float64
    delay::Float64
    freq_ppm::Float64
    phase_ppm::Float64
    freq::Float64
    phase::Float64
    phase_id::Int
end

struct PulseqExtensionInstanceEvent
    type::Int
    ref::Int
    next_id::Int
end

const ShapeLibrary = Dict{Int, Tuple{Int, Vector{Float64}}}

# Exact contents of the [DEFINITIONS] section plus the active raster values.
struct PulseqDefinitions
    entries::Vector{Pair{String,String}}
    block_duration_raster::Float64
    gradient_raster_time::Float64
    radiofrequency_raster_time::Float64
    adc_raster_time::Float64
    required_extensions::Vector{String}
end

PulseqDefinitions() = PulseqDefinitions(
    Pair{String,String}[],
    DEFAULT_RASTER.BlockDurationRaster,
    DEFAULT_RASTER.GradientRasterTime,
    DEFAULT_RASTER.RadiofrequencyRasterTime,
    DEFAULT_RASTER.AdcRasterTime,
    String[],
)

# Pulseq file-format libraries keyed by the integer ids used in [BLOCKS].
struct PulseqFileEventLibraries
    grad_library::Dict{Int,PulseqGradEvent}
    rf_library::Dict{Int,PulseqRFEvent}
    adc_library::Dict{Int,PulseqADCEvent}
    tmp_delay_library::Dict{Int,Float64}
    shape_library::ShapeLibrary
    extension_instance_library::Dict{Int,PulseqExtensionInstanceEvent}
    extension_type_library::Dict{Int,Type{<:Extension}}
    extension_spec_library::Dict{Int,Dict{Int,Extension}}
    definitions::PulseqDefinitions
end

# One serialized [BLOCKS] row.
struct PulseqBlockEventIDs
    block_id::Int
    duration_ticks::Int
    rf_id::Int
    gx_id::Int
    gy_id::Int
    gz_id::Int
    adc_id::Int
    ext_id::Int
end
struct PulseqRaster
    BlockDurationRaster::Float64
    GradientRasterTime::Float64
    RadiofrequencyRasterTime::Float64
    AdcRasterTime::Float64
    PulseqRaster(seq) = new(
        get_raster_time("BlockDurationRaster", seq),
        get_raster_time("GradientRasterTime", seq),
        get_raster_time("RadiofrequencyRasterTime", seq),
        get_raster_time("AdcRasterTime", seq),
    )
    PulseqRaster(sys::Scanner) = new(
        sys.DUR_Δt,
        sys.GR_Δt,
        sys.RF_Δt,
        sys.ADC_Δt,
    )
    PulseqRaster(seq, sys) = new(
        get_raster_time("BlockDurationRaster", seq, sys.DUR_Δt),
        get_raster_time("GradientRasterTime", seq, sys.GR_Δt),
        get_raster_time("RadiofrequencyRasterTime", seq, sys.RF_Δt),
        get_raster_time("AdcRasterTime", seq, sys.ADC_Δt),
    )
end

function get_raster_time(key::String, seq::KomaMRIBase.Sequence)
    haskey(seq.DEF, key) || error("Sequence has no Pulseq raster definition for `$key`.")
    return seq.DEF[key]
end

function get_raster_time(key::String, seq::KomaMRIBase.Sequence, scanner_value)
    haskey(seq.DEF, key) || return scanner_value
    value = seq.DEF[key]
    value == scanner_value || @warn "Sequence and Scanner definition for $key do not match (($(value) != $(scanner_value))). Using the Scanner definition ($key = $scanner_value)."
    return scanner_value
end

PulseqTrapGradEvent(data) = PulseqTrapGradEvent(data...)
PulseqArbGradEvent(data::NTuple{3, T}) where {T<:Real} = PulseqArbGradEvent(data[1], 0.0, 0.0, Int(data[2]), 0, data[3])
PulseqArbGradEvent(data::NTuple{4, T}) where {T<:Real} = PulseqArbGradEvent(data[1], 0.0, 0.0, Int(data[2]), Int(data[3]), data[4])
PulseqArbGradEvent(data::NTuple{6, T}) where {T<:Real} = PulseqArbGradEvent(data[1], data[2], data[3], Int(data[4]), Int(data[5]), data[6])

PulseqRFEvent(data::NTuple{6, T}) where {T<:Real} = PulseqRFEvent(data[1], Int(data[2]), Int(data[3]), 0, nothing, data[4], 0.0, 0.0, data[5], data[6], 'u')
PulseqRFEvent(data::NTuple{7, T}) where {T<:Real} = PulseqRFEvent(data[1], Int(data[2]), Int(data[3]), Int(data[4]), nothing, data[5], 0.0, 0.0, data[6], data[7], 'u')
PulseqRFEvent(data::Tuple{Vararg{Any, 11}}) = PulseqRFEvent(data[1], Int(data[2]), Int(data[3]), Int(data[4]), data[5], data[6], data[7], data[8], data[9], data[10], data[11])

PulseqADCEvent(data::NTuple{5, T}) where {T<:Real} = PulseqADCEvent(Int(data[1]), data[2], data[3], 0.0, 0.0, data[4], data[5], 0)
PulseqADCEvent(data::NTuple{8, T}) where {T<:Real} = PulseqADCEvent(Int(data[1]), data[2], data[3], data[4], data[5], data[6], data[7], Int(data[8]))

PulseqExtensionInstanceEvent(data) = PulseqExtensionInstanceEvent(Int(data[1]), Int(data[2]), Int(data[3]))
