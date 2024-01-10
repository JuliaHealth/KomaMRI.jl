"""
    delay = Delay(T)

The Delay struct is designed to introduce a delay into a sequence through the utilization of
the sum operator.

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

Displays the delay time in milliseconds of a Delay struct in the Julia REPL.

# Arguments
- `s`: (`::Delay`) delay struct

# Returns
- `str`: (`::String`) output string message
"""
Base.show(io::IO, d::Delay) = begin
	print(io, "Delay($(d.T*1e3)ms)")
end

"""
    seq = +(s::Sequence, d::Delay)
    seq = +(d::Delay, s::Sequence)

Introduces a delay to a sequence struct, ultimately influencing the duration of the
gradients within the sequence.

# Arguments
- `s`: (`::Sequence`) sequence struct
- `d`: (`::Delay`) delay struct

# Returns
- `seq`: (`::Sequence`) delayed sequence
"""
+(s::Sequence, d::Delay) = s + Sequence([Grad(0.,d.T)])
+(d::Delay, s::Sequence) = Sequence([Grad(0.,d.T)]) + s
