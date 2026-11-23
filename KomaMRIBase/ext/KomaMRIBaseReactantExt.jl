module KomaMRIBaseReactantExt

using KomaMRIBase
using Reactant

import KomaMRIBase: Phantom, RF, _has_waveform, VECTOR_PHANTOM_FIELDS,
    phantom_from_vectors, phantom_with_concrete_vectors
import KomaMRIBase: DiscreteSequence, sampled_sequence_storage, same_boundary_sample,
    append_waveform_samples!

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

end
