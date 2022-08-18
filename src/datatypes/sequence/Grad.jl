"""
    Rx = rotx(θ::Real)

Rotates vector counter-clockwise with respect to the x-axis.

# Arguments
- `θ::Real`: rotation angle

# Returns
- `Rx::Matrix{Int64}`: rotation matrix
"""
rotx(θ::Real) = [1		0		0
				 0	cos(θ)	-sin(θ);
				 0	sin(θ)	cos(θ)]

"""
    Ry = roty(θ::Real)

Rotates vector counter-clockwise with respect to the y-axis.

# Arguments
- `θ::Real`: rotation angle

# Returns
- `Ry::Matrix{Int64}`: rotation matrix
"""
roty(θ::Real) = [cos(θ) 0	sin(θ);
				 0		1		0;
				-sin(θ)	0 	cos(θ)]

"""
    Rz = rotz(θ::Real)

Rotates vector counter-clockwise with respect to the z-axis.

# Arguments
- `θ::Real`: rotation angle

# Returns
- `Rz::Matrix{Int64}`: rotation matrix
"""
rotz(θ::Real) = [cos(θ) -sin(θ)	0;
				 sin(θ) cos(θ)	0;
				 0 		0 		1]

"""
    Grad(A, T)
    Grad(A, T, rise)
    Grad(A, T, rise, delay)
    Grad(A, T, rise, fall, delay)

The Gradient object.

# Arguments
- `A`: (::Float64 or ::Vector{Float64}) amplitude in [T], if this is a scalar it makes a
    trapezoid
- `T`: (::Float64 or ::Vector{Float64}) duration of the flat-top [s]
- `rise::Real`: duration of the rise [s]
- `fall::Real`: duration of the fall [s]
- `delay::Real`: duration of the delay [s]
"""
mutable struct Grad
	A           #Amplitud [T], if this is a scalar it makes a trapezoid
	T           #Duration of flat-top [s]
	rise::Real  #Duration of rise [s]
	fall::Real  #Duration of fall [s]
	delay::Real #Duration of fall [s]
    function Grad(A, T, rise, fall, delay)
		all(T .< 0) || rise < 0 || fall < 0 || delay < 0 ? error("Gradient timings must be positive.") : new(A, T, rise, fall, delay)
    end
	function Grad(A, T, rise, delay)
		all(T .< 0) < 0 || rise < 0 || delay < 0 ? error("Gradient timings must be positive.") : new(A, T, rise, rise, delay)
    end
	function Grad(A, T, rise)
		all(T .< 0) < 0 || rise < 0 ? error("Gradient timings must be positive.") : new(A, T, rise, rise, 0)
    end
	function Grad(A, T)
		all(T .< 0) < 0 ? error("Gradient timings must be positive.") : new(A, T, 0, 0, 0)
    end
end

"""
    Grad(f::Function, T::Real, N::Int64)

Generates an arbitrary gradient waveform defined by function `f` in the interval t ∈
[0,`T`]. It uses `N` square gradients uniformly spaced in the interval.

# Arguments
- `f::Function`: the gradient waveform
- `T::Real`: the gradient duration
- `N::Int64`: the number of samples of the gradient
"""
Grad(f::Function, T::Real, N::Int64=300) = begin
	t = range(0, T; length=N)
	G = f.(t)
	Grad(G, T)
end

"""
    str = show(io::IO, x::Grad)

Displays information about the Grad object `x` in the julia REPL.

# Arguments
- `x::Grad`: the Grad object

# Returns
- `str` (::String) the output string message
"""
Base.show(io::IO, x::Grad) = begin
	r(x) = round.(x,digits=4)
	compact = get(io, :compact, false)
	if !compact
		wave = length(x.A) == 1 ? r(x.A*1e3) : "∿"
		if x.rise == x.fall == 0.
			print(io, (x.delay>0 ? "←$(r(x.delay*1e3)) ms→ " : "")*"Grad($(wave) mT, $(r(x.T*1e3)) ms)")
		else
			print(io, (x.delay>0 ? "←$(r(x.delay*1e3)) ms→ " : "")*"Grad($(wave) mT, $(r(x.T*1e3)) ms, ↑$(r(x.rise*1e3)) ms, ↓$(r(x.fall*1e3)) ms)")
		end
	else
		wave = length(x.A) == 1 ? "⊓" : "∿"
		print(io, (sum(abs.(x.A)) > 0 ? wave : "⇿")*"($(r((x.delay+x.rise+x.fall+x.T)*1e3)) ms)")
	end
end

"""
    y = getproperty(x::Vector{Grad}, f::Symbol)
    y = getproperty(x::Matrix{Grad}, f::Symbol)

Overchages Base.getproperty(). It is meant to access properties of the Grad vector `x`
directly without the need to iterate elementwise.

# Arguments
- `x::Vector{Grad}`: the vector of Grad objects
- `x::Matrix{Grad}`: the matrix of Grad objects
- `f::Symbol`: custom options are the `:x`, `:y`, `:z`, `:T`, `:delay`, `:rise`, `:delay`
    and `:dur` symbols

# Returns
- `y`: (::Vector{Any} or ::Matrix{Any}) the vector with the property defined by the symbol
    `f` for all elements of the Grad vector or matrix `x`
"""
getproperty(x::Vector{Grad}, f::Symbol) = getproperty.(x,f)
getproperty(x::Matrix{Grad}, f::Symbol) = begin
	if f == :x
		x[1,:]
	elseif f == :y && size(x,1) >= 2
		x[2,:]
	elseif f == :z && size(x,1) >= 3
		x[3,:]
	elseif f == :T || f == :rise || f == :fall || f == :delay
		getproperty.(x,f)
	elseif f == :dur
		T, ζ1, ζ2, delay = x.T, x.rise, x.fall, x.delay
		ΔT = [sum(t) for t=T] .+ ζ1 .+ ζ2 .+ delay
		maximum(ΔT,dims=1)[:]
	else
		getproperty.(x,f)
	end
end

# Gradient operations
*(x::Grad,α::Real) = Grad(α*x.A,x.T,x.rise,x.fall,x.delay)
*(α::Real,x::Grad) = Grad(α*x.A,x.T,x.rise,x.fall,x.delay)
*(A::Matrix{Float64},GR::Matrix{Grad}) = begin
	N, M = size(GR)
	[sum(A[i,1:N] .* GR[:,j]) for i=1:N, j=1:M]
end
zero(::Grad) = Grad(0,0)
size(g::Grad, i::Int64) = 1 #To fix [g;g;;] concatenation of Julia 1.7.3
/(x::Grad,α::Real) = Grad(x.A/α,x.T,x.rise,x.fall,x.delay)
+(x::Grad,y::Grad) = Grad(x.A+y.A,x.T,x.rise,x.fall,x.delay)
+(x::Array{Grad,1}, y::Array{Grad,1}) = [x[i]+y[i] for i=1:length(x)]
-(x::Grad) = -1*x
-(x::Grad,y::Grad) = Grad(x.A-y.A,x.T,x.rise,x.fall,x.delay)

# Gradient functions
vcat(x::Array{Grad,1},y::Array{Grad,1}) = [i==1 ? x[j] : y[j] for i=1:2,j=1:length(x)]
vcat(x::Array{Grad,1},y::Array{Grad,1},z::Array{Grad,1}) = [i==1 ? x[j] : i==2 ? y[j] : z[j] for i=1:3,j=1:length(x)]

"""
   y = dur(x::Grad)
   y = dur(x::Vector{Grad})

Duration time in [s] of Grad object or Grad array.

# Arguments
- `x`: (::Grad or ::Vector{Grad}) the RF object or RF array

# Returns
- `y`: (::Float64) the duration of the RF object or RF array in [s]
"""
dur(x::Grad) = x.delay + x.rise + sum(x.T) + x.fall
dur(x::Vector{Grad}) = maximum(dur.(x), dims=1)[:]
