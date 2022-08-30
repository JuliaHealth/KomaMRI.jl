"""
    delay = Delay(T)

The Delay struct. The input delay time `T` must be non-negative.

!!! note
    This struct is meant to add delays to a sequence struct that ultimately affects to the
    duration of the gradients of a sequence.

# Arguments
- `T`: (`::Real`, `[s]`) the time delay value

# Returns
- `delay`: (`::Delay`) the Delay struct
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
- `s`: (`::Delay`) the delay struct

# Returns
- `str`: (`::String`) the output string message

# Examples
``` julia-repl
julia> x = Delay(1)
Delay(1000.0ms)

julia> x
Delay(1000.0ms)

julia> show(x)
Delay(1000.0ms)

julia> display(x)
Delay(1000.0ms)
```
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
- `s`: (`::Sequence`) the sequence struct
- `d`: (`::Delay`) the delay struct

# Returns
- `seq`: (`::Sequence`) the delayed sequence
"""
+(s::Sequence, d::Delay) = s + Sequence([Grad(0.,d.T)])
+(d::Delay, s::Sequence) = Sequence([Grad(0.,d.T)]) + s
