"""
    rf = RF(A, T)
    rf = RF(A, T, Δf)
    rf = RF(A, T, Δf, delay)

The RF struct.

# Arguments
- `A`: (`::Complex{Int64}`, `[T]`) the amplitud-phase B1x + i B1y
- `T`: (`::Int64`, [`s`]) the duration of the RF
- `Δf`: (`::Float64`, [`Hz`]) the frequency offset of the RF
- `delay`: (`::Float64`, [`s`]) the delay time of the RF

# Returns
- `rf`: (`::RF`) the RF struct
"""
mutable struct RF
	A
	T
	Δf
	delay::Float64
	function RF(A, T, Δf, delay)
        any(T .< 0) || delay < 0 ? error("RF timings must be non-negative.") : new(A, T, Δf, delay)
    end
	function RF(A, T, Δf)
        any(T .< 0) ? error("RF timings must be non-negative.") : new(A, T, Δf, 0.)
    end
	function RF(A, T)
        any(T .< 0) ? error("RF timings must be non-negative.") : new(A, T, 0., 0.)
    end
end

"""
    str = show(io::IO, x::RF)

Displays information about the RF struct `x` in the julia REPL.

# Arguments
- `x`: (`::RF`) RF struct

# Returns
- `str`: (`::String`) output string message
"""
Base.show(io::IO, x::RF) = begin
	r(x) = round.(x,digits=4)
	compact = get(io, :compact, false)
	if !compact
		wave = length(x.A) == 1 ? r(x.A*1e6) : "∿"
		print(io, (x.delay>0 ? "←$(r(x.delay*1e3)) ms→ " : "")*"RF($(wave) uT, $(r(sum(x.T)*1e3)) ms, $(r(x.Δf)) Hz)")
	else
		wave = length(x.A) == 1 ? "⊓" : "∿"
		print(io, (sum(abs.(x.A)) > 0 ? wave : "⇿")*"($(r((x.delay+sum(x.T))*1e3)) ms)")
	end
end

"""
    y = getproperty(x::Vector{RF}, f::Symbol)
    y = getproperty(x::Matrix{RF}, f::Symbol)

Overloads Base.getproperty(). It is meant to access properties of the RF vector `x`
directly without the need to iterate elementwise.

# Arguments
- `x`: (`::Vector{RF}` or `::Matrix{RF}`) vector or matrix of RF structs
- `f`: (`::Symbol`, opts: [`:A`, `:Bx`, `:By`, `:T`, `:Δf`, `:delay` and `:dur`]) input
    symbol that represents a property of the vector or matrix of RF structs

# Returns
- `y`: (`::Vector{Any}` or `::Matrix{Any}`) vector with the property defined by the
    symbol `f` for all elements of the RF vector or matrix `x`
"""
getproperty(x::Vector{RF}, f::Symbol) = getproperty.(x,f)
getproperty(x::Matrix{RF}, f::Symbol) = begin
	if     f == :Bx
		real.(getproperty.(x,:A))
	elseif f == :By
		imag.(getproperty.(x,:A))
	elseif f == :Δf
		getproperty.(x,:Δf)
	elseif f == :T || f == :delay
		getproperty.(x[1,:],f)
	elseif f == :dur
		T, delay = x.T, x.delay
		ΔT = [sum(t) for t=T] .+ delay
		ΔT
	else
		getproperty.(x,f)
	end
end

# RF comparison
Base.isapprox(rf1::RF, rf2::RF) = begin
    return all(length(getfield(rf1, k)) == length(getfield(rf2, k)) for k ∈ fieldnames(RF))
        all(≈(getfield(rf1, k), getfield(rf2, k), atol=1e-9) for k ∈ fieldnames(RF))
end

# Properties
size(r::RF, i::Int64) = 1 #To fix [r;r;;] concatenation of Julia 1.7.3
*(α::Complex{T}, x::RF) where {T<:Real} = RF(α*x.A,x.T,x.Δf,x.delay)

"""
    y = dur(x::RF)
    y = dur(x::Array{RF,1})
    y = dur(x::Array{RF,2})

Duration time in [s] of RF struct or RF array.

# Arguments
- `x`: (`::RF` or `::Array{RF,1}` or `::Array{RF,2}`) RF struct or RF array

# Returns
- `y`: (`::Float64`, [`s`]) duration of the RF struct or RF array
"""
dur(x::RF) = sum(x.T)
dur(x::Array{RF,1}) = sum(sum(x[i].T) for i=1:size(x,1))
dur(x::Array{RF,2}) = maximum(sum([sum(x[i,j].T) for i=1:size(x,1),j=1:size(x,2)],dims=2))

"""
    rf = RF_fun(f::Function, T::Real, N::Int64)

Generate an RF sequence with amplitudes sampled from a function waveform.

!!! note
    This function is not being used in this KomaMRI version.

# Arguments
- `f`: (`::Function`, [`T`]) function for the RF amplitud waveform
- `T`: (`::Real`, [`s`]) duration of the RF pulse
- `N`: (`::Int64`) number of samples of the RF pulse

# Returns
- `rf`:(`::RF`) RF struct with amplitud defined by the function `f`
"""
RF(f::Function, T::Real, N::Int64=301; delay::Real=0, Δf=0) = begin
	t = range(0,T;length=N)
	A = f.(t)
	RF(A,T,Δf,delay)
end

# DEPRECATED?
"""
    α = get_flip_angle(x::RF)

Calculates the flip angle α [deg] of an RF struct. α = γ ∫ B1(τ) dτ

# Arguments
- `x`: (`::RF`) RF struct

# Returns
- `α`: (`::Int64`, `[deg]`) flip angle RF struct `x`
"""
get_flip_angle(x::RF) = begin
	A, NA, T, NT = x.A, length(x.A), x.T, length(x.T)
	dT = T / NA * NT
	α = round(360 * γ * abs(sum(A .* dT)), digits=3) #Pulseq
	return α
end

# DEPRECATED?
"""
    t = get_RF_center(x::RF)

Calculates the time where is the center of the RF pulse `x`. This calculation includes the
RF delay.

# Arguments
- `x`: (`::RF`) RF struct

# Returns
- `t`: (`::Int64`, `[s]`) time where is the center of the RF pulse `x`
"""
get_RF_center(x::RF) = begin
	A, NA, T, NT, delay = x.A, length(x.A), x.T, length(x.T), x.delay
	dT = T / NA * NT .* ones(NA)
	t = cumsum([0; dT])[1:end-1]
	i_center = argmax(abs.(A))
	t_center = t[i_center] + dT[i_center]/2
	return t_center + delay
end



############################################################################################
############################################################################################
############################################################################################
"""
For detecting if the rf event is on
"""
function ison(rf::RF)
    return ((sum(abs.(rf.A)) != 0.) && (sum(rf.T) != 0.))
end

"""
For getting the time, amplitude and carrier-frequency-difference samples of the rf event
"""
function samples(rf::RF)
    Δt, t, a, Δfc = Float64[], Float64[], ComplexF64[], Float64[]
    if ison(rf)
        NT, NA, NΔf = length(rf.T), length(rf.A), length(rf.Δf)
        if (NA == 1 && NT == 1)         # Block Pulse
            Δt = [rf.delay; 0.; rf.T; 0.]
            t  = cumsum(Δt)
            a  = [0.; rf.A; rf.A; 0.]
            Δfc = [0.; rf.Δf; rf.Δf; 0.]
        elseif (NA > 1 && NT == 1)
            Δt = [rf.delay; 0.; (ones(NA-1).*(rf.T/(NA-1))); 0.]
            t = cumsum(Δt)
            a = [0.; rf.A; 0.]
            Δfc = ((NΔf == 1) ? (rf.Δf .* ones(NA+2)) : ([rf.Δf[1]; rf.Δf; rf.Δf[end]]))
        elseif (NA > 1 && NT > 1)
            Δt = [rf.delay; 0.; (ones(NT) .* rf.T); 0.]
            t = cumsum(Δt)
            a = [0.; rf.A; 0.]
            Δfc = ((NΔf == 1) ? (rf.Δf .* ones(NA+2)) : ([rf.Δf[1]; rf.Δf; rf.Δf[end]]))
        end
    end
    return (Δt = Δt, t = t, a = a, Δfc = Δfc)
end

"""
For getting the time and amplitude samples of the rf event at times given
by the "ts" vector which must be increasing and with unique time points
"ts" must have at least 2 samples.
"""
function samples(rf::RF, ts::Vector{Float64})
    Δt, t, a, Δfc, ionfirst, ionlast, ix = (ts[2:end]-ts[1:end-1]), ts, zeros(length(ts)), zeros(length(ts)), 0, 0, 0
    if ison(rf)
        rfe = samples(rf)                       # Get the samples of the event
        rfs = interpolate(rfe.t, rfe.a, ts)     # Get interpolated and extrapolated values
        rfΔfs = interpolate(rfe.t, rfe.Δfc, ts) # Get interpolated and extrapolated values (the user must put the same samples for A and Δf)
        t, a, Δfc = rfs.t, rfs.a, rfΔfs.a       # Assign returned values
        Δt = t[2:end] - t[1:end-1]              # Assign returned values
        ionfirst, ionlast = rfs.ion[1], rfs.ion[2]   # Assign returned values for on indexes
        ix = argmin(abs.(center(rf).t .- t))    # Index of the RF center (must be present in "ts" in order to get the exact point)
    end
    return (Δt = Δt, t = t, a = a, Δfc = Δfc, ion = (ionfirst, ionlast), ix = ix)
end

"""
For getting the flipangle and type (excitation or refocusing)
"""
function flipangle(rf::RF)
    α, type = NaN, false
    if ison(rf)
        rfs = samples(rf)
        α = 180. * γ * abs(sum((rfs.a[2:end] + rfs.a[1:end-1]) .* rfs.Δt[2:end]))
        type = (α <= 90.01)
    end
    return (α = α, type = type)
end

"""
For getting the center
"""
function center(rf::RF)
    tx, ax = NaN, NaN
    if ison(rf)
        if length(rf.A) == 1
            tx, ax = rf.delay + rf.T/2, rf.A
        else
            rfs = samples(rf)
            ix = argmax(abs.(rfs.a))
            tx, ax = rfs.t[ix], rfs.a[ix]
        end
    end
    return (t = tx, a = ax)
end

"""
For getting critical times that must be simulated and considered for the block-samples
For rfs the extremes and the center are considered critical times
"""
function criticaltimes(rf::RF)
    tc = Float64[]
    if ison(rf)
        t = samples(rf).t
        tx = center(rf).t
        tc = [t[1]; tx; t[end]]
    end
    return tc
end

"""
For adding properties to the struct (rf event) calculated dynamically
"""
function Base.getproperty(rf::RF, sym::Symbol)
    (sym == :ison  ) && return ison(rf)
    (sym == :tc    ) && return criticaltimes(rf)
    (sym == :Δt    ) && return samples(rf).Δt
    (sym == :t     ) && return samples(rf).t
    (sym == :a     ) && return samples(rf).a
    (sym == :Δfc   ) && return samples(rf).Δfc
    (sym == :α     ) && return flipangle(rf).α
    (sym == :type  ) && return flipangle(rf).type
    (sym == :center) && return center(rf)
    return getfield(rf, sym)
end

"""
For displaying additional properties of the struct (rf event) when a user
press TAB twice in the REPL. This additional properties must be added previously in the
Base.getproperty() function
"""
function Base.propertynames(::RF)
    return (:A, :T, :rise, :fall, :delay,
            :ison, :tc, :Δt, :t, :a, :Δfc,
            :α, :type, :center)
end
