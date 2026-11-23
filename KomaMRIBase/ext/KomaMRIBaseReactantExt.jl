module KomaMRIBaseReactantExt

using KomaMRIBase
using Reactant

import KomaMRIBase: Phantom, RF, _has_waveform, VECTOR_PHANTOM_FIELDS,
    phantom_from_vectors, phantom_with_concrete_vectors
import KomaMRIBase: DiscreteSequence, sampled_sequence_storage, same_boundary_sample,
    append_waveform_samples!
import KomaMRIBase: linear_interpolate_samples, _foreach_interpolation_sample

# Keep native Phantom assignment semantics and select concrete storage at upload.
function Reactant.to_rarray(obj::Phantom; kwargs...)
    fields = Tuple(VECTOR_PHANTOM_FIELDS)
    vectors = map(field -> getfield(obj, field), fields)
    concrete = phantom_with_concrete_vectors(obj.name, vectors, obj.motion)
    return invoke(Reactant.to_rarray, Tuple{Any}, concrete; kwargs...)
end

function phantom_from_vectors(name, vectors::Tuple{Reactant.AnyTracedRVector,Vararg}, motion)
    return phantom_with_concrete_vectors(name, vectors, motion)
end

# A traced RF waveform has a statically known shape, but its values cannot be
# inspected on the host to determine whether the event is active.
_has_waveform(rf::RF{<:Reactant.AnyTracedRVector}) = !isempty(rf.A)

function sampled_sequence_storage(out, B1::Vector, ::Reactant.AnyTracedRVector)
    return DiscreteSequence(out.Gx, out.Gy, out.Gz, Reactant.promote_to(Reactant.TracedRArray, B1),
        out.Δf, out.ψ, out.ADC, out.excitation_bool, out.t, out.Δt)
end

same_boundary_sample(out, values, B1::Reactant.AnyTracedRVector) = false

append_waveform_samples!(dst::Reactant.TracedRArray, src, rows) =
    append!(dst, Reactant.promote_to(dst, src[rows]))

_gather_samples(A, indices) = A[indices]
# Avoid malformed reverse IR for complex gathers followed by strided slices.
_gather_samples(A::Reactant.AnyTracedRVector{<:Complex}, indices) =
    complex.(real.(A)[indices], imag.(A)[indices])

function linear_interpolate_samples(
    samples::NamedTuple{(:t, :A), <:Tuple{Any, Reactant.AnyTracedRVector}}, t;
    default=zero(eltype(samples.A)), interpolate=true,
)
    if isempty(samples.A) || isempty(samples.t)
        return fill!(similar(samples.A, typeof(default), length(t)), default)
    end
    lo, hi = ones(Int, length(t)), ones(Int, length(t))
    weights = zeros(float(promote_type(eltype(t), eltype(samples.t))), length(t))
    valid = fill(false, length(t))
    _foreach_interpolation_sample(samples, t, interpolate) do indices, l, h, w
        lo[indices] .= l
        hi[indices] .= h
        weights[indices] .= w
        valid[indices] .= true
    end
    lower, upper = _gather_samples(samples.A, lo), _gather_samples(samples.A, hi)
    values = lower .+ (upper .- lower) .* weights
    values = ifelse.(lo .== hi, lower, values)
    out = similar(samples.A, typeof(default), length(t))
    out .= ifelse.(valid, values, default)
    return out
end

end
