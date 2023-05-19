"""
    delay = Delay(T)

The Delay struct.

# Arguments
- `T`: (`::Real`, `[s]`) time delay value

# Returns
- `delay`: (`::Delay`) delay struct
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
+(s::Sequence, d::Delay) = s + Sequence([Grad(0.,d.T)])
+(d::Delay, s::Sequence) = Sequence([Grad(0.,d.T)]) + s
