module KomaMRIBaseReactantExt

using KomaMRIBase
using Reactant

import KomaMRIBase: RF, _has_waveform

# A traced RF waveform has a statically known shape, but its values cannot be
# inspected on the host to determine whether the event is active.
_has_waveform(rf::RF{<:Reactant.AnyTracedRVector}) = !isempty(rf.A)

end
