"""
    delay = Delay(T)

The Delay struct is meant to add a delay to a sequence by using a sum operator.

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
    T::Real
    function Delay(T)
		T < 0 ? error("Delays must be positive.") : new(T)
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

"""
    seq = +(s::Sequence, d::Delay)
    seq = +(d::Delay, s::Sequence)

Add a delay to sequence struct. It ultimately affects to the duration of the gradients of a
sequence.

# Arguments
- `s`: (`::Sequence`) sequence struct
- `d`: (`::Delay`) delay struct

# Returns
- `seq`: (`::Sequence`) delayed sequence
"""
+(s::Sequence, d::Delay) = s + empty_seq(d.T)
+(d::Delay, s::Sequence) = empty_seq(d.T) + s
function empty_seq(T)
    seq = Sequence([Grad(0., 0.);;])
    seq.DUR[1] = T
    return seq
end
