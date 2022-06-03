"Rotates vector counter-clockwise with respecto to the x-axis."
rotx(θ::Real) = [1		0		0
				 0	cos(θ)	-sin(θ);
				 0	sin(θ)	cos(θ)]
"Rotates vector counter-clockwise with respecto to the y-axis."
roty(θ::Real) = [cos(θ) 0	sin(θ);
				 0		1		0;
				-sin(θ)	0 	cos(θ)]
"Rotates vector counter-clockwise with respecto to the z-axis."
rotz(θ::Real) = [cos(θ) -sin(θ)	0;
				 sin(θ) cos(θ)	0;
				 0 		0 		1]
#####################
## Gradient OBJECT ##
#####################
# @everywhere begin
"""
    Grad(A,T)

Gradient object.

# Arguments
- `A::Real` := Gradient amplitude [T].
- `T::Real` := Gradient duration [s].

# Examples
```@example
Grad(1,2)
```
"""
mutable struct Grad
	A #Amplitud [T], if this is a scalar it makes a trapezoid
	T     #Duration of flat-top [s]
	rise::Real  #Duration of rise [s]
	fall::Real  #Duration of fall [s]
	delay::Real #Duration of fall [s]
    function Grad(A,T,rise,fall,delay)
		all(T .< 0) || rise < 0 || fall < 0 || delay < 0 ? error("Gradient timings must be positive.") : new(A, T, rise, fall, delay)
    end
	function Grad(A,T,rise,delay)
		all(T .< 0) < 0 || rise < 0 || delay < 0 ? error("Gradient timings must be positive.") : new(A, T, rise, rise, delay)
    end
	function Grad(A,T,rise)
		all(T .< 0) < 0 || rise < 0 ? error("Gradient timings must be positive.") : new(A, T, rise, rise, 0)
    end
	function Grad(A,T)
		all(T .< 0) < 0 ? error("Gradient timings must be positive.") : new(A, T, 0, 0, 0)
    end
end
#Gradient operations
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
#Gradient functions
vcat(x::Array{Grad,1},y::Array{Grad,1}) = [i==1 ? x[j] : y[j] for i=1:2,j=1:length(x)]
vcat(x::Array{Grad,1},y::Array{Grad,1},z::Array{Grad,1}) = [i==1 ? x[j] : i==2 ? y[j] : z[j] for i=1:3,j=1:length(x)]
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

#TIMINGS
"Duration `T` [s] of Grad datatype."
dur(x::Grad) = x.delay + x.rise + sum(x.T) + x.fall
"Duration `T` [s] of Vector{Grad}."
dur(x::Vector{Grad}) = maximum(dur.(x), dims=1)[:]

"""
	Grad_fun(f,T,N)

Generates **arbitrary gradient waveform** defined by function `f` in the interval t∈[0,`T`].
It uses `N` square gradients uniformly spaced in the interval.

# Arguments
 - `f::Function` := Gradient waveform.
 - `T::Real`     := Gradient duration.
 - `N::Integer=1`:= Number of sample points.

# Examples
```julia-repl
julia> Grad_fun(x-> sin(π*x),1,4)
1×4 Array{Grad,2}:
 Grad(0.0, 0.333333)  Grad(0.866025, 0.333333)  Grad(0.866025, 0.333333)  Grad(1.22465e-16, 0.333333)
```
"""
Grad(f::Function,T::Real,N::Int64=300) = begin
	t = range(0,T,N)
	G = f.(t')
	Grad(G,T)
end
#aux
Base.show(io::IO,x::Grad) = begin
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