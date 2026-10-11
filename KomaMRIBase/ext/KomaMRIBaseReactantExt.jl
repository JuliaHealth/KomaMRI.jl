module KomaMRIBaseReactantExt

using KomaMRIBase
using Reactant

import KomaMRIBase: RF, is_on
import KomaMRIBase: linear_interpolate_samples, _foreach_interpolation_sample

# A traced RF waveform has a statically known shape, but its values cannot be
# inspected on the host to determine whether the event is active.
is_on(rf::RF{<:Reactant.AnyTracedRVector}) = !isempty(rf.A)

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
