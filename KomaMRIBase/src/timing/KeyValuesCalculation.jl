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
function get_theo_A(gr::Grad)
    if !(gr.A isa Vector) && !(gr.T isa Vector)     # trapezoidal
        return [0; gr.A; gr.A; 0]
	end
	return [0; gr.A; 0]     # uniformly-sampled and time-shaped
end
function get_theo_A(rf::RF, key::Symbol)
	rfA = getproperty(rf, key)
    if !(rfA isa Vector) && !(rf.T isa Vector)      # pulse
        return [0; rfA; rfA; 0]
	end
	return [0; rfA; 0]     # uniformly-sampled and time-shaped
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
    t = get_theo_t(gr::Grad)
    t = get_theo_t(rf::RF)
    t = get_theo_t(adc::ADC)

Get the theoretical times for Grad, RF or ADC structs.

# Arguments
- `gr`: (`::Grad`) Gradient struct
- `rf`: (`::RF`) RF struct
- `adc`: (`::ADC`) ADC truct

# Returns
- `t`: (`::Vector{Number}`) vector with the time theoretical points
"""
function get_theo_t(gr::Grad)
    if !(gr.A isa Vector) && !(gr.T isa Vector)
        return cumsum([gr.delay; gr.rise; gr.T; gr.fall])   # trapezoidal
    elseif gr.A isa Vector && gr.T isa Vector
        NT = length(gr.T)
        return cumsum([gr.delay; gr.rise; gr.T.*ones(NT); gr.fall])    # time-shaped
    end
    NA = length(gr.A)
    return cumsum([gr.delay; gr.rise; gr.T/(NA-1).*ones(NA-1); gr.fall])    # uniformly-sampled
end
function get_theo_t(rf::RF, key::Symbol)
	NT, T, NA = length(g.T), g.T, length(getproperty(g, key))
	if sum(T) != 0
		if !(getproperty(g, key) isa Vector) && !(g.T isa Vector)
			trf = [0; T]
        elseif getproperty(g, key) isa Vector && g.T isa Vector
			trf = cumsum([0; T.*ones(NT)])
		elseif getproperty(g, key) isa Vector
			trf = cumsum([0; T/(NA-1).*ones(NA-1)])
		end
		t = g.delay .+ trf
	else
		t = [g.delay; g.delay]
	end
    t[end] -= MIN_RISE_TIME #previous float, avoids incorrect block category
	aux = [0; g.delay; t; g.delay+sum(g.T)]
	aux
end
function get_theo_t(adc::ADC)
	t = adc.delay .+ [0; adc.T]
	t = t[1:end-1]
	return [0; t; adc.delay + adc.T]
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
	N = length(seq)
	T0 = get_block_start_times(seq)
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


get_theo_RF(seq, key::Symbol) = begin
	N = length(seq)
	T0 = get_block_start_times(seq)
	t = vcat([get_theo_t(seq.RF[i], key) .+ T0[i] for i=1:N]...)
	A = vcat([get_theo_A(seq.RF[i], key) for i=1:N]...)
	Interpolations.deduplicate_knots!(t; move_knots=true)
	return (t, A)
end
