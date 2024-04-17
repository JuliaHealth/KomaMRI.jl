"""
    A = ampl(g::Grad)
    A = ampl(r::RF)
    A = ampl(d::ADC)

Get the theoretical amplitudes for Grad, RF or ADC structs.

# Arguments
- `gr`: (`::Grad`) Gradient struct
- `rf`: (`::RF`) RF struct
- `adc`: (`::ADC`) ADC truct

# Returns
- `A`: (`::Vector{Number}`) vector with the amplitude theoretical points
"""
function ampl(gr::Grad)
    if !is_on(gr)
        return Float64[]
    end
    if !(gr.A isa Vector) && !(gr.T isa Vector)     # trapezoidal
        return [gr.first; gr.A; gr.A; gr.last]
	end
	return [gr.first; gr.A; gr.last]     # uniformly-sampled and time-shaped
end
function ampl(rf::RF; freq_in_phase=false)
    if !is_on(rf)
        return ComplexF64[]
    end
    if !(rf.A isa Vector)
        A = ComplexF64[0.0; rf.A; rf.A; 0.0]
	else
        A = ComplexF64[0.0; rf.A; 0.0]     # uniformly-sampled and time-shaped
    end
    if freq_in_phase
        t0 = time(rf, :Δf)
        t  = time(rf)
        Interpolations.deduplicate_knots!(t0; move_knots=true)
        Interpolations.deduplicate_knots!(t; move_knots=true)
        f = linear_interpolation(t0, freq(rf); extrapolation_bc=0.0).(t)
        A = A .* exp.(im*2π*[0; cumtrapz(diff(t), f)])
    end
	return A
end
function freq(rf::RF)
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
function ampl(adc::ADC)
    if !is_on(adc)
        return Bool[]
    end
    return ones(Bool, adc.N)
end


"""
    t = time(gr::Grad)
    t = time(rf::RF)
    t = time(adc::ADC)

Get the theoretical times for Grad, RF or ADC structs.

# Arguments
- `gr`: (`::Grad`) Gradient struct
- `rf`: (`::RF`) RF struct
- `adc`: (`::ADC`) ADC truct

# Returns
- `t`: (`::Vector{Number}`) vector with the time theoretical points
"""
function time(gr::Grad)
    if !is_on(gr)
        return Float64[]
    end
    if !(gr.A isa Vector) && !(gr.T isa Vector)
        t = cumsum([gr.delay; gr.rise; gr.T; gr.fall])   # trapezoidal
    elseif gr.A isa Vector && gr.T isa Vector
        NT = length(gr.T)
        t =  cumsum([gr.delay; gr.rise; gr.T.*ones(NT); gr.fall])    # time-shaped
    else
        NA = length(gr.A)
        t = cumsum([gr.delay; gr.rise; gr.T/(NA-1).*ones(NA-1); gr.fall]) # uniformly-sampled
    end
    t[end-1] -= MIN_RISE_TIME #Fixes incorrect block interpretation
    return t
end
function time(rf::RF, key::Symbol)
    if !is_on(rf)
        return Float64[]
    end
    rfA = getproperty(rf, key)
    if !(rfA isa Vector) && !(rf.T isa Vector)
        t =  cumsum([rf.delay; 0.0; rf.T; 0.0])         # pulse
    elseif rfA isa Vector && rf.T isa Vector
        NT = length(rf.T)
        t =  cumsum([rf.delay; 0.0; rf.T.*ones(NT); 0.0])    # time-shaped
    elseif !(rfA isa Vector)
        t =  cumsum([rf.delay; 0.0; sum(rf.T); 0.0]) # df constant
    else
        NA = length(rf.A)
        t = cumsum([rf.delay; 0.0; rf.T/(NA-1).*ones(NA-1); 0.0])    # uniformly-sampled
    end
    t[end-1] -= MIN_RISE_TIME #Fixes incorrect block interpretation
    return t
end
time(rf::RF) = time(rf, :A)
function time(adc::ADC)
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
