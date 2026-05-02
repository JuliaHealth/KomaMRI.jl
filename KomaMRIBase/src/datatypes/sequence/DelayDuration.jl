"""
    delay = Delay(T)

The Delay struct adds a delay block, or sets a minimum block duration when used
inside `addblock!`.

# Arguments
- `T`: (`::Real`, `[s]`) time delay value

# Returns
- `delay`: (`::Delay`) delay struct

# Examples
```julia-repl
julia> delay = Delay(0.5)

julia> s = Sequence([Grad(1, 1, 0.1)])

julia> seq = delay + s; plot_seq(seq)
```
"""
struct Delay
    T::Float64
    function Delay(T)
        T < 0 ? error("Delays must be positive.") : new(Float64(T))
    end
end

"""
    duration = Duration(T)

Strict block-duration marker. When used inside `addblock!`, the block duration
is set to `T` and an error is thrown if any event is longer.
"""
struct Duration
    T::Float64
    function Duration(T)
        T < 0 ? error("Durations must be positive.") : new(Float64(T))
    end
end

"""
    str = show(io::IO, s::Delay)

Displays the delay time in m[s] of the delay struct `s` in the julia REPL.

# Arguments
- `s`: (`::Delay`) delay struct

# Returns
- `str`: (`::String`) output string message
"""
Base.show(io::IO, s::Delay) = begin
    print(io, "Delay($(s.T*1e3)ms)")
end
Base.show(io::IO, s::Duration) = begin
    print(io, "Duration($(s.T*1e3)ms)")
end

dur(d::Delay) = d.T
dur(d::Duration) = d.T

function empty_seq(T)
    seq = Sequence([Grad(0., 0.);;])
    seq.DUR[1] = T
    return seq
end

const _DelayBlockEvent = Union{_BlockEvent,Delay,Duration}
const _DelayBlockEventTuple = Tuple{_DelayBlockEvent,Vararg{_DelayBlockEvent}}

_check_block_event(::Delay) = nothing
_check_block_event(::Duration) = nothing
_block_duration(block_duration, ::Duration) = block_duration
_block_duration_constraint(d::Duration) = dur(d)

+(s::Sequence, d::Delay) = s + empty_seq(dur(d))
+(d::Delay, s::Sequence) = empty_seq(dur(d)) + s
+(s::Sequence, d::Duration) = s + empty_seq(dur(d))
+(d::Duration, s::Sequence) = empty_seq(dur(d)) + s
addblock!(seq::Sequence, d::Delay) = append!(seq, empty_seq(dur(d)))
addblock!(seq::Sequence, d::Duration) = append!(seq, empty_seq(dur(d)))
function addblock!(seq::Sequence, events::_DelayBlockEventTuple; x=nothing, y=nothing, z=nothing)
    return append!(seq, _block_sequence(events; x, y, z))
end
Base.append!(seq::Sequence, d::Delay) = addblock!(seq, d)
Base.append!(seq::Sequence, d::Duration) = addblock!(seq, d)
+(s::Sequence, events::_DelayBlockEventTuple) = s + _block_sequence(events)
+(events::_DelayBlockEventTuple, s::Sequence) = _block_sequence(events) + s
_addblock_term!(seq::Sequence, events::_DelayBlockEventTuple) = addblock!(seq, events)
_addblock_term!(seq::Sequence, events::_DelayBlockEvent...) = addblock!(seq, events...)
