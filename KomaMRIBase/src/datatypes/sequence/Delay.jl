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

# Display on the REPL
Base.show(io::IO, d::Delay) = begin
	print(io, "Delay($(d.T*1e3)ms)")
end

"""
    seq = +(s::Sequence, delay::Delay)
    seq = +(delay::Delay, s::Sequence)

Introduces a delay to a sequence struct, ultimately influencing the duration of the
gradients within the sequence.

# Arguments
- `s`: (`::Sequence`) sequence struct
- `delay`: (`::Delay`) delay struct

# Returns
- `seq`: (`::Sequence`) delayed sequence
"""
+(s::Sequence, d::Delay) = s + Sequence([Grad(0.,d.T)])
+(d::Delay, s::Sequence) = Sequence([Grad(0.,d.T)]) + s
