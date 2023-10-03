"""
    Rx = rotx(θ::Real)

Rotates vector counter-clockwise with respect to the x-axis.

# Arguments
- `θ`: (`::Real`, `[rad]`) rotation angle

# Returns
- `Rx`: (`::Matrix{Int64}`) rotation matrix
"""
rotx(θ::Real) = [1		0		0
				 0	cos(θ)	-sin(θ);
				 0	sin(θ)	cos(θ)]

"""
    Ry = roty(θ::Real)

Rotates vector counter-clockwise with respect to the y-axis.

# Arguments
- `θ`: (`::Real`, `[rad]`) rotation angle

# Returns
- `Ry`: (`::Matrix{Int64}`) rotation matrix
"""
roty(θ::Real) = [cos(θ) 0	sin(θ);
				 0		1		0;
				-sin(θ)	0 	cos(θ)]

"""
    Rz = rotz(θ::Real)

Rotates vector counter-clockwise with respect to the z-axis.

# Arguments
- `θ`: (`::Real`, `[rad]`) rotation angle

# Returns
- `Rz`: (`::Matrix{Int64}`) rotation matrix
"""
rotz(θ::Real) = [cos(θ) -sin(θ)	0;
				 sin(θ) cos(θ)	0;
				 0 		0 		1]

"""
    grad = Grad(A, T)
    grad = Grad(A, T, rise)
    grad = Grad(A, T, rise, delay)
    grad = Grad(A, T, rise, fall, delay)

The Gradient struct.

# Arguments
- `A`: (`::Float64`, `[T/m]`) amplitude of the gradient
- `T`: (`::Float64`, `[s]`) duration of the flat-top
- `rise`: (`::Real`, `[s]`) duration of the rise
- `fall`: (`::Real`, `[s]`) duration of the fall
- `delay`: (`::Real`, `[s]`) duration of the delay

# Returns
- `grad`: (`::Grad`) gradient struct
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
    grad = Grad(f::Function, T::Real, N::Int64; delay::Real)

Generates an arbitrary gradient waveform defined by function `f` in the interval t ∈
[0,`T`]. It uses `N` square gradients uniformly spaced in the interval.

# Arguments
- `f`: (`::Function`) function that describes the gradient waveform
- `T`: (`::Real`, `[s]`) duration of the gradient waveform
- `N`: (`::Int64`) number of samples of the gradient waveform

# Keywords
- `delay`: (`::Real`, `=0`, `[s]`) starting delay for the waveform

# Returns
- `grad`: (`::Grad`) gradient struct

# Examples
```julia-repl
julia> f1 = t -> sin(π*t / 0.8)

julia> seq = Sequence([Grad(f1, 0.8)])

julia> plot_seq(seq)
```
"""
Grad(f::Function, T::Real, N::Int64=300; delay::Real=0) = begin
	t = range(0,T;length=N)
	G = f.(t)
	Grad(G,T,0,0,delay)
end


"""
    str = show(io::IO, x::Grad)

Displays information about the Grad struct `x` in the julia REPL.

# Arguments
- `x`: (`::Grad`) Grad struct

# Returns
- `str` (`::String`) output string message
"""
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
    y = dur(x::Grad)
    y = dur(x::Vector{Grad})

Duration time in [s] of Grad struct or Grad array. When the input is a gradient vector, then
the duration is the maximum duration of all the elements of the gradient vector.

# Arguments
- `x`: (`::Grad` or `::Vector{Grad}`) RF struct or RF array

# Returns
- `y`: (`::Float64`, `[s]`) duration of the RF struct or RF array
"""
dur(x::Grad) = x.delay + x.rise + sum(x.T) + x.fall
dur(x::Vector{Grad}) = maximum(dur.(x), dims=1)[:]



############################################################################################
############################################################################################
############################################################################################
"""
For detecting if the gradient event is on
"""
function ison(gr::Grad)
    return (sum(abs.(gr.A)) != 0.)
end

"""
For getting the time and amplitude samples of the gradient event
"""
function samples(gr::Grad)
    Δt, t, a = Float64[], Float64[], Float64[]
    if ison(gr)
        NT, NA = length(gr.T), length(gr.A)
        if (NA == 1 && NT == 1)     # Trapezoid
            Δt = [gr.delay; gr.rise; gr.T; gr.fall]
            t = cumsum(Δt)
            a = [0.; gr.A; gr.A; 0.]
        elseif (NA > 1 && NT == 1)  # Uniformly-Sampled waveform
            Δt = [gr.delay; gr.rise; (ones(NA-1).*(gr.T/(NA-1))); gr.fall]
            t = cumsum(Δt)
            a = [0.; gr.A; 0.]
        elseif (NA > 1 && NT > 1)   # Time-Shaped waveform
            Δt = [gr.delay; gr.rise; (ones(NT) .* gr.T); gr.fall]
            t = cumsum(Δt)
            a = [0.; gr.A; 0.]
        end
    end
    return (Δt = Δt, t = t, a = a)
end

"""
For getting the time and amplitude samples of the gradient event at times given
by the "ts" vector which must be increasing and with unique time points
"ts" must have at least 2 samples.
"""
function samples(gr::Grad, ts::Vector{Float64})
    Δt, t, a, ionfirst, ionlast = (ts[2:end]-ts[1:end-1]), ts, zeros(length(ts)), 0, 0
    if ison(gr)
        gre = samples(gr)                       # Get the samples of the event
        grs = interpolate(gre.t, gre.a, ts)     # Get interpolated and extrapolated values
        t, a = grs.t, grs.a                     # Assign returned values
        Δt = t[2:end] - t[1:end-1]              # Assign returned values
        ionfirst, ionlast = grs.ion[1], grs.ion[2]  # Assign returned values for on indexes
    end
    return (Δt = Δt, t = t, a = a, ion = (ionfirst, ionlast))
end

"""
For getting critical times that must be simulated and considered for the block-samples
For gradients the times of the trapezoid are critical
"""
function criticaltimes(gr::Grad)
    tc = Float64[]
    if ison(gr)
        t = samples(gr).t
        tc = [t[1]; t[2]; t[end-1]; t[end]]
    end
    return tc
end

"""
For adding properties to the struct (gradient event) calculated dynamically
"""
function Base.getproperty(gr::Grad, sym::Symbol)
    (sym == :ison) && return ison(gr)
    (sym == :tc  ) && return criticaltimes(gr)
    (sym == :Δt  ) && return samples(gr).Δt
    (sym == :t   ) && return samples(gr).t
    (sym == :a   ) && return samples(gr).a
    return getfield(gr, sym)
end

"""
For displaying additional properties of the struct (gradient event) when a user
press TAB twice in the REPL. This additional properties must be added previously in the
Base.getproperty() function
"""
function Base.propertynames(::Grad)
    return (:A, :T, :rise, :fall, :delay,
            :ison, :tc, :Δt, :t, :a)
end
