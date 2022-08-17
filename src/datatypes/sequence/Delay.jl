"""
    Delay(T)

The Delay object. The input delay time `T` must be non-negative.

!!! note
    This object is meant to add delays to a sequence object that ultimately affects to the
    gradients of a sequence.

# Arguments
- `T::Real`: the time delay value in [s]
"""
struct Delay
    T::Real #Time in [s]
    function Delay(T)
		T < 0 ? error("Delays must be positive.") : new(T)
    end
end

#Interactions with other objects are context aware object
#Sequence contatenation
+(s::Sequence, d::Delay) = s + Sequence([Grad(0.,d.T)])
+(d::Delay, s::Sequence) = Sequence([Grad(0.,d.T)]) + s

#aux
"""
    str = show(io::IO, s::Delay)

Displays delay time in m[s] in the julia REPL.

# Arguments
- `s::Delay`: the delay object

# Returns
- `str` (::String) the output string message

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
