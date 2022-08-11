"""
    Delay

The Delay object.
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
Base.show(io::IO,s::Delay) = begin
	print(io, "Delay($(s.T*1e3)ms)")
end
