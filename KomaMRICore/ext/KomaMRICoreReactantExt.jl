module KomaMRICoreReactantExt

using Reactant

import KomaMRICore: _scalar_getindex, _scalar_setindex!

const TracedArray = AbstractArray{<:Reactant.TracedRNumber}

function _scalar_getindex(x::TracedArray, i)
    return Reactant.allowscalar() do
        x[i]
    end
end

function _scalar_setindex!(x::TracedArray, value, i)
    Reactant.allowscalar() do
        x[i] = value
    end
    return x
end

end
