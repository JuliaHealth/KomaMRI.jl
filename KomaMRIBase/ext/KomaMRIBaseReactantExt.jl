module KomaMRIBaseReactantExt

using KomaMRIBase
using Reactant

import KomaMRIBase: Phantom, VECTOR_PHANTOM_FIELDS,
    phantom_from_vectors, phantom_with_concrete_vectors

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

end
