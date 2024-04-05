"""
    A = get_theo_A(g::Grad)
    A = get_theo_A(r::RF)
    A = get_theo_A(d::ADC)

Get the theoretical amplitudes for Grad, RF or ADC structs.

# Arguments
- `gr`: (`::Grad`) Gradient struct
- `rf`: (`::RF`) RF struct
- `adc`: (`::ADC`) ADC truct

# Returns
- `A`: (`::Vector{Number}`) vector with the amplitude theoretical points
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
function get_theo_A(adc::ADC)
    return ones(adc.N)
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
    rfA = getproperty(rf, key)
    if !(rfA isa Vector) && !(rf.T isa Vector)
        return cumsum([rf.delay; 0; rf.T; 0])         # pulse
    elseif rfA isa Vector && rf.T isa Vector
        NT = length(rf.T)
        return cumsum([rf.delay; 0; rf.T.*ones(NT); 0])    # time-shaped
    end
    NA = length(rf.A)
    return cumsum([rf.delay; 0; rf.T/(NA-1).*ones(NA-1); 0])    # uniformly-sampled
end
function get_theo_t(adc::ADC)
	return adc.delay .+ ((adc.N == 1) ? [adc.T/2] : [range(0, adc.T; length=adc.N)...])
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
