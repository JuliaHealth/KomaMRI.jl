module KomaMRIBase

#IMPORT PACKAGES
import Base.*, Base.+, Base.-, Base./, Base.vcat, Base.size, Base.abs, Base.getproperty
#General
using Reexport
#Datatypes
using Parameters
#Simulation
using Interpolations
#Reconstruction
using MRIBase
@reexport using MRIBase:
    Profile, RawAcquisitionData, AcquisitionData, AcquisitionHeader, EncodingCounters, Limit
using MAT   # For loading example phantoms

const global γ = 42.5774688e6 # Hz/T gyromagnetic constant for H1, JEMRIS uses 42.5756 MHz/T

_deepcopy_fields(x) = ntuple(i -> deepcopy(getfield(x, i)), fieldcount(typeof(x)))

fields_equal(x, y) =
    typeof(x) === typeof(y) &&
    all(i -> getfield(x, i) == getfield(y, i), 1:fieldcount(typeof(x)))

field_isapprox(x, y; kwargs...) = x == y
field_isapprox(x::Number, y::Number; kwargs...) = isapprox(x, y; kwargs...)
field_isapprox(x::AbstractArray{<:Number}, y::AbstractArray{<:Number}; kwargs...) = isapprox(x, y; kwargs...)
function field_isapprox(x::AbstractArray, y::AbstractArray; kwargs...)
    axes(x) == axes(y) || return false
    return all(eachindex(x, y)) do i
        field_isapprox(x[i], y[i]; kwargs...)
    end
end

fields_isapprox(x, y; kwargs...) =
    typeof(x) === typeof(y) &&
    all(i -> field_isapprox(getfield(x, i), getfield(y, i); kwargs...), 1:fieldcount(typeof(x)))

_has_negative_timings(T::Real) = T < 0
_has_negative_timings(T::AbstractVector{<:Real}) = any(<(0), T)

_shape_samples(x::Number) = [zero(x), x, x, zero(x)]
_shape_samples(x::AbstractVector{<:Number}) = [zero(eltype(x)); x; zero(eltype(x))]

_shape_times(::Number, T::Real) = [zero(T), T]
_shape_times(::Number, T::AbstractVector{<:Number}) = [zero(eltype(T)), sum(T)]
_shape_times(x::AbstractVector, T::Real) = range(zero(T), T; length=length(x))
function _shape_times(x::AbstractVector, T::AbstractVector{<:Number})
    n = length(x)
    n == 0 && return similar(T, 0)
    if length(T) == n - 1
        return cumsum([zero(eltype(T)); T])
    elseif length(T) == n
        return cumsum([zero(eltype(T)); T[1:(end - 1)]])
    end
    throw(DimensionMismatch("Expected time vector of length $(n - 1) or $n for $n samples, got $(length(T))."))
end

# Hardware
include("datatypes/Scanner.jl")
# Sequence
include("datatypes/sequence/Grad.jl")
include("datatypes/sequence/RF.jl")
include("datatypes/sequence/ADC.jl")
include("datatypes/sequence/EXT.jl")
include("timing/KeyValuesCalculation.jl")
include("datatypes/Sequence.jl")
include("datatypes/sequence/RotationExtensions.jl")
include("datatypes/sequence/TimingChecks.jl")
include("datatypes/sequence/HardwareChecks.jl")
include("datatypes/sequence/AddBlockMacro.jl")
include("datatypes/sequence/DelayDuration.jl")
# Motion
include("motion/MotionList.jl")
include("motion/NoMotion.jl")
# Phantom
include("datatypes/Phantom.jl")
# Simulator
include("datatypes/simulation/DiscreteSequence.jl")
include("timing/TimeStepCalculation.jl")
include("timing/TrapezoidalIntegration.jl")

# Main
export γ    # gyro-magnetic ratio [Hz/T]
export Scanner, Sequence, Phantom
export addblock!, @addblock, @addblocks
export Grad, RF, ADC, Delay, Duration, QuaternionRot
export dur, get_block_start_times, get_samples
export RFuse, Excitation, Refocusing, Inversion, Saturation, Preparation, Other, Undefined
export DiscreteSequence
export discretize, get_adc_phase_compensation, get_adc_sampling_times
export is_Gx_on, is_Gy_on, is_Gz_on, is_RF_on, is_ADC_on
export times, ampls, freqs
# These are also used for simulation
export kfoldperm, trapz, cumtrapz
# Phantom
export brain_phantom2D, brain_phantom3D, pelvis_phantom2D, heart_phantom
# Motion
export MotionList, NoMotion, Motion
export translate, rotate, heartbeat, path, flowpath
export Translate, TranslateX, TranslateY, TranslateZ
export Rotate, RotateX, RotateY, RotateZ, CenterOfMass
export HeartBeat, Path, FlowPath
export TimeRange, Periodic, TimeCurve
export SpinRange, AllSpins
export get_spin_coords
# Secondary
export get_kspace, rotx, roty, rotz
# Additionals
export get_flip_angles, is_RF_on, is_GR_on, is_ADC_on
# Sequence related
export get_Mk, get_kspace, get_M0, get_M1, get_M2, get_labels
export check_timing, check_hw_limits, apply_rotations

# PulseDesigner submodule
include("sequences/PulseDesigner.jl")
export PulseDesigner

end # module KomaMRIBase
