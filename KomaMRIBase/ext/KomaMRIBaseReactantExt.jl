module KomaMRIBaseReactantExt

using KomaMRIBase
using Reactant

import KomaMRIBase: RF, _scalar_getindex, _scalar_setindex!, is_on

const ReactantRVector = Union{
    Reactant.AnyTracedRVector,
    Reactant.AbstractConcreteArray{<:Number,1},
}

# A traced RF waveform has a statically known shape, but its values cannot be
# inspected on the host to determine whether the event is active.
is_on(rf::RF{<:ReactantRVector}) = !isempty(rf.A)

function _scalar_getindex(x::ReactantRVector, i)
    return Reactant.allowscalar() do
        x[i]
    end
end

function _scalar_setindex!(x::ReactantRVector, value, i)
    Reactant.allowscalar() do
        x[i] = value
    end
    return x
end

end
