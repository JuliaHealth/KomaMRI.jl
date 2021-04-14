using MRIsim
spin = Phantom()

mutable struct Spin
    x
    y
    z
    T1
    T2
end
s = Spin(0,0,0,1,0.08)

Base.show(io::IO,s::Spin) = begin
    print(io,"Spin[ T1 = $(s.T1*1e3) ms | T2 = $(s.T2*1e3) ms ]")
end