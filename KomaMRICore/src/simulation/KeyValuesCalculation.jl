"""
    A = get_theo_A(g::Grad; off_val=0)
    A = get_theo_A(r::RF; off_val=0, max_rf_samples=Inf)
    A = get_theo_A(d::ADC; off_val=0)

Get the theoretical amplitudes of a rectangle waveform for Grad, RF or ADC structs. This are
5 points: delay, start, rise, stop and fall.

!!! note
    In some cases the array result can have duplicated points, so it is necessary to remove
    them whenever necessary.

# Arguments
- `g`: (`::Grad`) Gradient struct
- `r`: (`::RF`) RF struct
- `d`: (`::ADC`) ADC truct

# Keywords
- `off_val`: (`::Float64`, `=0`) offset value for amplitude. In general, it is used for
    not showing some points in plots by giving an `Inf` value
- `max_rf_samples`: (`::Float64`, `=Inf`) number of maximum samples for the RF struct.
    In general, this parameter is not necessary to set

# Returns
- `A`: (`::Vector{Float64}`) vector with the amplitude key points of the rectangle
    waveform
"""
get_theo_A(g::Grad; off_val=0) = begin
	A = g.A
	if length(A) == 1
		aux = [off_val; 0; A; A; 0]
	else
		aux = [off_val; 0; A; 0]
	end
	if sum(abs.(g.A)) == 0
		aux *= off_val
	end
	aux
end
#function get_theo_A(gr::Grad)
#    ((sum(abs.(gr.A)) == 0) || (sum(gr.T) == 0)) && return []
#    (length(gr.A) == 1) && return [0; gr.A; gr.A; 0]
#    return [0; gr.A; 0]
#end

get_theo_A(r::RF; off_val=0, max_rf_samples=Inf) = begin
	A = [r.A r.A]'; A = A[:]
	if sum(abs.(r.A)) == 0
		aux = off_val * [1; 1; ones(length(A)); 1]
	else
		#Calculating frequency waveform
		NA = length(r.A)
		NT = length(r.T)
		if NA > 1 && NT == 1
			dt = repeat([r.T], outer=length(r.A))
		elseif NA > 1 && NT > 1
			dt = r.T[:]
		elseif NA == 1 && NT == 1
			dt = r.T
		end
		freq = exp.(1im*2π*cumsum(r.Δf.*dt)) # exp(i ∫ᵗw(t)dt)
		added_phase = transpose([freq freq])[:]
		A = A .* added_phase
		A[abs.(A).<=1e-10] .= 0 #Remove small values
		#Output
		aux = [off_val; 0; A; 0]
	end
	#Subsample
	if length(aux) > max_rf_samples
		n = floor(Int, length(aux) / max_rf_samples)
		aux = aux[[1; 2:n:end-1; end]]
	end
	aux
end

get_theo_A(d::ADC; off_val=0) = begin
	N = [1 1]'; N = N[:]
	if d.N == 0
		aux = off_val * [1; 1; ones(length(N)); 1]
	else
		aux = [off_val; 0; N; 0]
	end
	aux
end


"""
    t = get_theo_t(g::Grad)
    t = get_theo_t(r::RF; max_rf_samples=Inf)
    t = get_theo_t(d::ADC)

Get the theoretical times of a rectangle waveform for Grad, RF or ADC structs. This are
5 points: delay, start, rise, stop and fall.

!!! note
    In some cases the array result can have duplicated points, so it is necessary to remove
    them whenever necessary.

# Arguments
- `g`: (`::Grad`) Gradient struct
- `r`: (`::RF`) RF struct
- `d`: (`::ADC`) ADC truct

# Keywords
- `max_rf_samples`: (`::Float64`, `=Inf`) number of maximum samples for the RF struct.
    In general, this parameter is not necessary to set

# Returns
- `t`: (`::Vector{Float64}`) vector with the time key points of the rectangle waveform
"""
get_theo_t(g::Grad) = begin
	NT, T, NA = length(g.T), g.T, length(g.A)
	if sum(T) != 0
		if NA == 1 && NT == 1
			trf = [0; T]
		elseif NA>1 && NT == 1
			trf = cumsum([0; T/(NA-1).*ones(NA-1)])
		elseif NA>1 && NT>1
			trf = cumsum([0; T.*ones(NT)])
		end
		t = g.delay .+ g.rise .+ trf
	else
		t = [g.delay+g.rise; g.delay+g.rise]
	end
	aux = [0; g.delay; t; g.delay+g.rise+sum(g.T)+g.fall]
	aux
end
#function get_theo_t(gr::Grad)
#    NT, NA = length(gr.T), length(gr.A)
#    ((sum(abs.(gr.A)) == 0) || (sum(gr.T) == 0)) && return []
#    (NA == 1 && NT == 1) && return cumsum([gr.delay; gr.rise; gr.T; gr.fall])
#    (NA > 1 && NT == 1) && return cumsum([gr.delay; gr.rise; (ones(NA-1).*(gr.T/(NA-1))); gr.fall])
#	return cumsum([gr.delay; gr.rise; (ones(NT) .* T); gr.fall])
#end

get_theo_t(r::RF; max_rf_samples=Inf) = begin
	if sum(r.T) != 0
		NA, T, NT = length(r.A), r.T, length(r.T)
		dT = T / NA * NT
		trf = cumsum([0; dT.*ones(NA)])
		t = r.delay .+ trf
	else
		t = [r.delay; r.delay]
	end
	t = [t t]'; t = t[:]
	t = [0; t]
	#Subsample
	if length(t) > max_rf_samples
		n = floor(Int, length(t) / max_rf_samples)
		t = t[[1; 2:n:end-1; end]]
	end
	t
end

get_theo_t(d::ADC) = begin
	t = [d.delay; d.delay+d.T]
	t = [t t]'; t = t[1:end-1]
	[0; t; d.delay+d.T]
end


"""
    t, g = get_theo_Gi(seq, idx)

Get the theoretical gradient for a sequence in a defined axis.

# Arguments
- `seq`: (`::Sequence`) Sequence struct
- `idx`: (`::Int64`, opts=[1, 2, 3]) axis x, y or z for the gradient

# Returns
- `t`: (`::Vector{Float64}`) time key points
- `g`: (`::Vector{Float64}`) amplitude key points
"""
get_theo_Gi(seq, idx) = begin
	ΔT, N = durs(seq), length(seq)
	T0 = cumsum([0; ΔT], dims=1)
	t = vcat([get_theo_t(seq.GR[idx,i]) .+ T0[i] for i=1:N]...)
	G = vcat([get_theo_A(seq.GR[idx,i]) for i=1:N]...) #; off_val=0 <---potential solution
	#Removing duplicated points
	#TODO: do this properly. As it is now it generates a bug for slew rates that are too high
	# mask = (G .== 0) #<---potential solution
	# t = t[mask]
	# G = G[mask]
	Interpolations.deduplicate_knots!(t; move_knots=true)
	return (t, G)
end

##################
### DEPRECATED ###
##################
#"""
#    t, r = get_theo_rf(seq, idx)
#
#Get the theoretical RF timings and amplitude of a sequence.
#
#!!! note
#    Experimental, not being used yet.
#
## Arguments
#- `seq`: (`::Sequence`) Sequence struct
#- `idx`: (`::Int64`, opts=[1, 2]) it selects 1:B1 or 2:Δf of the RF
#
## Returns
#- `t`: (`::Vector{Float64}`) time key points
#- `r`: (`::Vector{Float64}`) amplitude key points
#"""
#get_theo_RF(seq, idx) = begin
#	ΔT, N = durs(seq), length(seq)
#	T0 = cumsum([0; ΔT], dims=1)
#	t = vcat([get_theo_t(seq.RF[i]) .+ T0[i] for i=1:N]...)
#	R = vcat([get_theo_A(seq.RF[i]; off_val=0) for i=1:N]...)
#	#Removing duplicated points
#	#TODO: do this properly. As it is now it generates a bug for slew rates that are too high
#	mask = (G .== 0)
#	# Interpolations.deduplicate_knots!(t; move_knots=true)
#	return (t[mask], R[mask])
#end


############################################################################################
### For Grad ###############################################################################
############################################################################################
"""
Get values from object definition
Returns a vector of unique times and a vector of vectors with up to two amplitudes at a unique times
Note that length(tu) = length(a)
It also return if the gradient is on
"""
function event_values(gr::Grad)

    # Define empty vectors to be filled and some params
    amps, times = Float64[], Float64[]
    NT, NA = length(gr.T), length(gr.A)
    ison = (sum(abs.(gr.A)) != 0.)

    # Get the amplitudes from object parameters definitions
    if !ison
        return (ison = ison, tu = Float64[], a = [Float64[]])    # empty values definition
    elseif (NA == 1 && NT == 1)
        amps, times = [0.; gr.A; gr.A; 0.], cumsum([gr.delay; gr.rise; gr.T; gr.fall])
    elseif (NA > 1 && NT == 1)
        amps, times = [0.; gr.A; 0.], cumsum([gr.delay; gr.rise; (ones(NA-1).*(gr.T/(NA-1))); gr.fall])
    elseif (NA > 1 && NT > 1)
        amps, times = [0.; gr.A; 0.], cumsum([gr.delay; gr.rise; (ones(NT) .* gr.T); gr.fall])
    end

    # Here, at the same time can be multiple amplitude values,
    # so we keep up to 2 amplitudes at the same time
    tu =  sort(unique(times))
    Ntu = length(tu)
    a = Vector{Vector{Float64}}(undef, Ntu)
    for i in eachindex(tu)
        mask = findall(times .== tu[i])
        a[i] = (length(mask) == 1) ? amps[mask] : [amps[mask][1]; amps[mask][end]]
    end

    return (ison = ison, tu = tu, a = a)
end

"""
Get amplitude values from object at certain time
It returns a vector of amplitudes with up to two samples
"""
function event_values(gr::Grad, t::Float64)
    # Get the values of the object definition
    egr = event_values(gr)
    ison, tu, amps = egr.ison, egr.tu, egr.a
    # Return zero when zero-signal or outside the critical-times
    (!ison || (t < tu[1]) || (tu[end] < t)) && return (a = [0.],)
    # Return the interpolated value when is in between two points of the signal
    if !(t in tu)
        i = findlast(tu .< t)
        m = (amps[i+1][1] - amps[i][end]) / (tu[i+1] - tu[i])
        return (a = [m * (t - tu[i]) + amps[i][end]],)
    end
    # Return the exact values when is an exact point of the signal
    return (a = amps[t .== tu][1],)
end

"""
Get multiple amplitude values at certain time points for the event
Returns a vector of length(t), each element of the vector is a vector with up to two amplitudes
"""
function event_values(gr::Grad, t::Vector{Float64})
    Nt = length(t)
    a = Vector{Vector{Float64}}(undef, Nt)
    for i in eachindex(t)
        grt = event_values(gr, t[i])
        a[i] = grt.a
    end
    return (a = a,)
end

"""
Returns a vector of times and a vector of gradient amplitudes
The times are not necessary unique, at the same time can be up to two amplitudes
"""
function event_samples(gr::Grad)
    egr = event_values(gr)
    ison = egr.ison
    t, a = Float64[], Float64[]
    for i in eachindex(egr.tu)
        append!(t, fill(egr.tu[i], length(egr.a[i])))
        append!(a, egr.a[i])
    end
    return (ison = ison, t = t, a = a)
end


############################################################################################
### For RF #################################################################################
############################################################################################
"""
Get values from object definition
Returns a vector of unique times and a vector of vectors with up to two amplitudes at a unique times
Note that length(tu) = length(a)
It also return if the gradient is on
"""
function event_values(rf::RF)

    # Define empty vectors to be filled and some params
    amps, Δfs, times = Float64[], Float64[], Float64[]
    NT, NA = length(rf.T), length(rf.A)
    ison = (sum(abs.(rf.A)) != 0) && (sum(rf.T) != 0)

    # Get the amplitudes from object parameters definitions
    if !ison
        return (ison = ison, tu = Float64[], a = [Float64[]], Δf = [Float64[]])    # empty values definition
    elseif (NA == 1 && NT == 1)
        amps, Δfs, times = [0.; rf.A; rf.A; 0.], [0.; 0.; rf.Δf; 0.], cumsum([rf.delay; 0.; rf.T; 0.])
    elseif (NA > 1 && NT == 1)
        amps, Δfs, times = [0.; rf.A; 0.], (rf.Δf .* ones(NA+2)), cumsum([rf.delay; 0.; (ones(NA-1).*(rf.T/(NA-1))); 0.])
    elseif (NA > 1 && NT > 1)
        amps, Δfs, times = [0.; rf.A; 0.], [rf.Δf[1]; rf.Δf; rf.Δf[end]], cumsum([rf.delay; 0.; (ones(NT) .* rf.T); 0.])
    end

    # Here, at the same time can be multiple amplitude values,
    # so we keep up to 2 amplitudes at the same time
    tu =  sort(unique(times))
    Ntu = length(tu)
    a = Vector{Vector{Float64}}(undef, Ntu)
    Δf = Vector{Vector{Float64}}(undef, Ntu)
    for i in eachindex(tu)
        mask = findall(times .== tu[i])
        Nmask = length(mask)
        a[i]  = (Nmask == 1) ? amps[mask] : [amps[mask][1]; amps[mask][end]]
        Δf[i] = (Nmask == 1) ? Δfs[mask] : [Δfs[mask][1]; Δfs[mask][end]]
    end

    # Return the object values
    return (ison = ison, tu = tu, a = a, Δf = Δf)
end

"""
Get amplitude and Δf values from object at certain time
It returns a vector of amplitudes with up to two samples for a and Δf
"""
function event_values(rf::RF, t::Float64)
    # Get the values of the object definition
    erf = event_values(rf)
    ison, tu, amps, Δfs = erf.ison, erf.tu, erf.a, erf.Δf
    # Return zero when zero-signal or outside the critical-times
    (!ison || (t < tu[1]) || (tu[end] < t)) && return (a = [0.], Δf = [0.])
    # Return the interpolated value when is in between two points of the signal
    if !(t in tu)
        i = findlast(tu .< t)
        ma = (amps[i+1][1] - amps[i][end]) / (tu[i+1] - tu[i])
        mΔf = (Δfs[i+1][1] - Δfs[i][end]) / (tu[i+1] - tu[i])
        return (a = [ma * (t - tu[i]) + amps[i][end]], Δf = [mΔf * (t - tu[i]) + Δfs[i][end]])
    end
    # Return the exact values when is an exact point of the signal
    mask = (t .== tu)
    return (a = amps[mask][1], Δf = Δfs[mask][1])
end

"""
Get multiple amplitude and Δf values at certain time points for the event
Returns a vector of length(t), each element of the vector is a vector with up to two amplitudes
"""
function event_values(rf::RF, t::Vector{Float64})
    Nt = length(t)
    a, Δf = Vector{Vector{Float64}}(undef, Nt), Vector{Vector{Float64}}(undef, Nt)
    for i in eachindex(t)
        rft = event_values(rf, t[i])
        a[i], Δf[i] = rft.a, rft.Δf
    end
    return (a = a, Δf = Δf)
end

"""
Returns a vector of times and a vector of rf amplitudes
The times are not necessary unique, at the same time can be up to two amplitudes
"""
function event_samples(rf::RF)
    erf = event_values(rf)
    ison = erf.ison
    t, a, Δf = Float64[], Float64[], Float64[]
    for i in eachindex(erf.tu)
        append!(t, fill(erf.tu[i], length(erf.a[i])))
        append!(a, erf.a[i])
        append!(Δf, erf.Δf[i])
    end
    return (ison = ison, t = t, a = a, Δf = Δf)
end

############################################################################################
### For ADC ################################################################################
############################################################################################
"""
Return adc values from object definition
This should have unique time points
It returns a vector
"""
function event_values(adc::ADC)
    # Return for zero adc
    ison = (adc.N > 0)
    (!ison) && return (ison = ison, t = Float64[])
    # Return for adc with samples
    Nsamp, T = adc.N, adc.T
    t = adc.delay .+ ((Nsamp == 1) ? ([T/2]) : (range(0, T; length=Nsamp)))
    return (ison = ison, t = t)
end

"""
This should have unique time points
It returns a vector
"""
function event_samples(adc::ADC)
    return event_values(adc)
end
