"""
    Rx = rotx(θ::Real)

Rotates a three-dimensional vector or matrix with three rows counter-clockwise with respect
to the x-axis.

# Arguments
- `θ`: (`::Real`, `[rad]`) rotation angle

# Returns
- `Rx`: (`::Matrix{Real}`) rotation matrix
"""
rotx(θ::Real) = [1		0		0
				 0	cos(θ)	-sin(θ);
				 0	sin(θ)	cos(θ)]

"""
    Ry = roty(θ::Real)

Rotates a three-dimensional vector or matrix with three rows counter-clockwise with respect
to the y-axis.

# Arguments
- `θ`: (`::Real`, `[rad]`) rotation angle

# Returns
- `Ry`: (`::Matrix{Real}`) rotation matrix
"""
roty(θ::Real) = [cos(θ) 0	sin(θ);
				 0		1		0;
				-sin(θ)	0 	cos(θ)]

"""
    Rz = rotz(θ::Real)

Rotates a three-dimensional vector or matrix with three rows counter-clockwise with respect
to the z-axis.

# Arguments
- `θ`: (`::Real`, `[rad]`) rotation angle

# Returns
- `Rz`: (`::Matrix{Real}`) rotation matrix
"""
rotz(θ::Real) = [cos(θ) -sin(θ)	0;
				 sin(θ) cos(θ)	0;
				 0 		0 		1]

"""
    gr = Grad(A, T)
    gr = Grad(A, T, rise)
    gr = Grad(A, T, rise, delay)
    gr = Grad(A, T, rise, fall, delay)

The Grad struct represents a gradient event of a sequence.

# Arguments
- `A`: (`::Real` or `::Vector`, `[T/m]`) amplitude of the gradient
- `T`: (`::Real` or `::Vector`, `[s]`) duration of the flat-top
- `rise`: (`::Real`, `[s]`) duration of the rise
- `fall`: (`::Real`, `[s]`) duration of the fall
- `delay`: (`::Real`, `[s]`) duration of the delay

# Returns
- `gr`: (`::Grad`) gradient struct

# Examples
```julia-repl
julia> gr = Grad(1, 1, 0.1, 0.1, 0.2)

julia> seq = Sequence([gr]); plot_seq(seq)
```
"""
mutable struct Grad
	A
	T
	rise::Real
	fall::Real
	delay::Real
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
    gr = Grad(f::Function, T::Real, N::Integer; delay::Real)

Generates an arbitrary gradient waveform defined by the function `f` in the interval t ∈
[0,`T`]. The time separation between two consecutive samples is given by T/(N-1).

# Arguments
- `f`: (`::Function`) function that describes the gradient waveform
- `T`: (`::Real`, `[s]`) duration of the gradient waveform
- `N`: (`::Integer`, `=300`) number of samples of the gradient waveform

# Keywords
- `delay`: (`::Real`, `=0`, `[s]`) delay time of the waveform

# Returns
- `gr`: (`::Grad`) gradient struct

# Examples
```julia-repl
julia> gx = Grad(t -> sin(π*t / 0.8), 0.8)

julia> seq = Sequence([gx]); plot_seq(seq)
```
"""
Grad(f::Function, T::Real, N::Integer=300; delay::Real=0) = begin
	t = range(0, T; length=N)
	G = f.(t)
	return Grad(G, T, 0, 0, delay)
end


# Display on the REPL
Base.show(io::IO, x::Grad) = begin
	r(x) = round.(x,digits=4)
	compact = get(io, :compact, false)
	if !compact
		wave = length(x.A) == 1 ? r(x.A*1e3) : "∿"
		if x.rise == x.fall == 0.
			print(io, (x.delay>0 ? "←$(r(x.delay*1e3)) ms→ " : "")*"Grad($(wave) mT, $(r(sum(x.T)*1e3)) ms)")
		else
			print(io, (x.delay>0 ? "←$(r(x.delay*1e3)) ms→ " : "")*"Grad($(wave) mT, $(r(sum(x.T)*1e3)) ms, ↑$(r(x.rise*1e3)) ms, ↓$(r(x.fall*1e3)) ms)")
		end
	else
		wave = length(x.A) == 1 ? "⊓" : "∿"
		print(io, (sum(abs.(x.A)) > 0 ? wave : "⇿")*"($(r((x.delay+x.rise+x.fall+sum(x.T))*1e3)) ms)")
	end
end

"""
    y = getproperty(x::Vector{Grad}, f::Symbol)
    y = getproperty(x::Matrix{Grad}, f::Symbol)

Overloads Base.getproperty(). It is meant to access properties of the Grad vector `x`
directly without the need to iterate elementwise.

# Arguments
- `x`: (`::Vector{Grad}` or `::Matrix{Grad}`) vector or matrix of Grad structs
- `f`: (`::Symbol`, opts: [`:x`, `:y`, `:z`, `:T`, `:delay`, `:rise`, `:delay`, `:dur`,
    `:A`, `f`]) input symbol that represents a property of the vector or matrix of Grad
    structs

# Returns
- `y`: (`::Vector{Any}` or `::Matrix{Any}`) vector or matrix with the property defined
    by the symbol `f` for all elements of the Grad vector or matrix `x`
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

# Gradient comparison
Base.isapprox(gr1::Grad, gr2::Grad) = begin
    return all(length(getfield(gr1, k)) ≈ length(getfield(gr2, k)) for k ∈ fieldnames(Grad)) &&
        all(getfield(gr1, k) ≈ getfield(gr2, k) for k ∈ fieldnames(Grad))
end

# Gradient operations
*(x::Grad,α::Real) = Grad(α*x.A,x.T,x.rise,x.fall,x.delay)
*(α::Real,x::Grad) = Grad(α*x.A,x.T,x.rise,x.fall,x.delay)
*(A::Matrix{Float64},GR::Matrix{Grad}) = begin
	N, M = size(GR)
	[sum(A[i,1:N] .* GR[:,j]) for i=1:N, j=1:M]
end
Base.zero(::Grad) = Grad(0,0)
size(g::Grad, i::Int64) = 1 #To fix [g;g;;] concatenation of Julia 1.7.3
/(x::Grad,α::Real) = Grad(x.A/α,x.T,x.rise,x.fall,x.delay)
+(x::Grad,y::Grad) = Grad(x.A.+y.A,x.T,x.rise,x.fall,x.delay)
+(x::Array{Grad,1}, y::Array{Grad,1}) = [x[i]+y[i] for i=1:length(x)]
-(x::Grad) = -1*x
-(x::Grad,y::Grad) = Grad(x.A.-y.A,x.T,x.rise,x.fall,x.delay)

# Gradient functions
vcat(x::Array{Grad,1},y::Array{Grad,1}) = [i==1 ? x[j] : y[j] for i=1:2,j=1:length(x)]
vcat(x::Array{Grad,1},y::Array{Grad,1},z::Array{Grad,1}) = [i==1 ? x[j] : i==2 ? y[j] : z[j] for i=1:3,j=1:length(x)]

"""
    time = dur(gr::Grad)
    time = dur(gr::Vector{Grad})

Duration time in seconds of a Grad struct or Grad array. When the input is a gradient vector,
the duration is the maximum duration among all the elements of the gradient vector.

# Arguments
- `gr`: (`::Grad` or `::Vector{Grad}`) Grad struct or Grad array

# Returns
- `time`: (`::Real`, `[s]`) duration of the Grad struct or Grad array
"""
dur(x::Grad) = x.delay + x.rise + sum(x.T) + x.fall
dur(x::Vector{Grad}) = maximum(dur.(x), dims=1)[:]
