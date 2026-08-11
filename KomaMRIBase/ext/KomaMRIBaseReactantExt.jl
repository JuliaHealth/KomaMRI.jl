module KomaMRIBaseReactantExt

using KomaMRIBase
using Reactant

import KomaMRIBase: RF, is_on

# A traced RF waveform has a statically known shape, but its values cannot be
# inspected on the host to determine whether the event is active.
is_on(rf::RF{<:Reactant.AnyTracedRVector}) = !isempty(rf.A)

end
