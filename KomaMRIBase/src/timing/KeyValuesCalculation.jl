"""
    A = ampls(g::Grad)
    A = ampls(r::RF)
    A = ampls(d::ADC)

Get amplitude samples of MRI sequence event.

# Arguments
- `gr`: (`::Grad`) Gradient struct
- `rf`: (`::RF`) RF struct
- `adc`: (`::ADC`) ADC struct

# Returns
- `A`: (`::Vector{Number}`) vector with amplitude samples
"""
function ampls(gr::Grad)
    if !is_on(gr)
        return Float64[]
    end
    if !(gr.A isa Vector) && !(gr.T isa Vector)     # trapezoidal
        return [gr.first; gr.A; gr.A; gr.last]
	end
	return [gr.first; gr.A; gr.last]     # uniformly-sampled and time-shaped
end
function ampls(rf::RF; freq_in_phase=false)
    if !is_on(rf)
        return ComplexF64[]
    end
    if !(rf.A isa Vector)
        A = ComplexF64[0.0; rf.A; rf.A; 0.0]
	else
        A = ComplexF64[0.0; rf.A; 0.0]     # uniformly-sampled and time-shaped
    end
    if freq_in_phase
        t0 = times(rf, :Δf)
        t  = times(rf)
        Interpolations.deduplicate_knots!(t0; move_knots=true)
        Interpolations.deduplicate_knots!(t; move_knots=true)
        f = linear_interpolation(t0, freqs(rf)).(t)
        A = A .* exp.(im*2π*[0; cumtrapz(diff(t), f)])
    end
	return A
end
"""
    f = freqs(r::RF)

Get frequency samples of MRI sequence event.

# Arguments
- `rf`: (`::RF`) RF struct

# Returns
- `f`: (`::Vector{Number}`) vector with frequency samples
"""
function freqs(rf::RF)
    if !is_on(rf)
        return Float64[]
    end
    if !(rf.Δf isa Vector)
        df = [0.0; rf.Δf; rf.Δf; 0.0]
	else
        df = [0.0; rf.Δf; 0.0]     # uniformly-sampled and time-shaped
    end
	return df
end
function ampls(adc::ADC)
    if !is_on(adc)
        return Bool[]
    end
    return ones(Bool, adc.N)
end

"""
    c = cents(rf::RF)

Get center time of RF pulse. It includes the RF delay.

# Arguments
- `rf`: (`::RF`) RF struct

# Returns
- `c`: (`::Vector{Number}`) vector with center time of the RF pulse `rf`
"""
function cents(rf::RF)
    if !is_on(rf)
        return Float64[]
    end
    return [rf.center + rf.delay]
end

"""
    t = times(gr::Grad)
    t = times(rf::RF)
    t = times(adc::ADC)

Get time samples of MRI sequence event.

# Arguments
- `gr`: (`::Grad`) Gradient struct
- `rf`: (`::RF`) RF struct
- `adc`: (`::ADC`) ADC struct

# Returns
- `t`: (`::Vector{Number}`) vector with time samples
"""
function times(gr::Grad)
    if !is_on(gr)
        return Float64[]
    end
    if !(gr.A isa Vector) && !(gr.T isa Vector)
        t = cumsum([gr.delay; gr.rise; gr.T; gr.fall])   # trapezoidal
    elseif gr.A isa Vector && gr.T isa Vector
        t =  cumsum([gr.delay; gr.rise; gr.T; gr.fall])    # time-shaped
    else
        NA = length(gr.A)
        t = cumsum([gr.delay; gr.rise; gr.T/(NA-1).*ones(NA-1); gr.fall]) # uniformly-sampled
    end
    t[end-1] -= MIN_RISE_TIME #Fixes incorrect block interpretation
    return t
end
function times(rf::RF, key::Symbol)
    if !is_on(rf)
        return Float64[]
    end
    rfA = getproperty(rf, key)
    if !(rfA isa Vector) && !(rf.T isa Vector)
        t =  cumsum([rf.delay; 0.0; rf.T; 0.0])         # pulse
    elseif rfA isa Vector && rf.T isa Vector
        t =  cumsum([rf.delay; 0.0; rf.T; 0.0])    # time-shaped
    elseif !(rfA isa Vector)
        t =  cumsum([rf.delay; 0.0; sum(rf.T); 0.0]) # df constant
    else
        NA = length(rf.A)
        t = cumsum([rf.delay; 0.0; rf.T/(NA-1).*ones(NA-1); 0.0])    # uniformly-sampled
    end
    t[end-1] -= MIN_RISE_TIME #Fixes incorrect block interpretation
    return t
end
times(rf::RF) = times(rf, :A)
function times(adc::ADC)
    if !is_on(adc)
        return range(0.0,0.0,0) #empty
    end
	if adc.N != 1
        t = range(0.0, adc.T; length=adc.N) .+ adc.delay
    else
        t = range(adc.T/2, adc.T/2, 1) .+ adc.delay
    end
    return t
end

is_on(x::Grad) = any(x->abs.(x) > .0, x.A)
is_on(x::RF)   = any(x->abs.(x) > .0, x.A)
is_on(x::ADC)  = x.N > 0
