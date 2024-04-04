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
    if !(g.A isa Vector) && !(g.T isa Vector)
        aux = [off_val; 0; A; A; 0]
	else
		aux = [off_val; g.first; A; g.last]
	end
    # If no signal, then don't show samples
	if sum(abs.(g.A)) == 0
		aux *= off_val
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
		if !(g.A isa Vector) && !(g.T isa Vector)
			trf = [0; T]
        elseif g.A isa Vector && g.T isa Vector
			trf = cumsum([0; T.*ones(NT)])
		elseif g.A isa Vector
			trf = cumsum([0; T/(NA-1).*ones(NA-1)])
		end
		t = g.delay .+ g.rise .+ trf
	else
		t = [g.delay+g.rise; g.delay+g.rise]
	end
	aux = [0; g.delay; t; g.delay+g.rise+sum(g.T)+g.fall]
	aux
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


get_theo_A(g::RF, key::Symbol; off_val=0, max_rf_samples=1) = begin
	A = getproperty(g, key)
    if !(g.A isa Vector) && !(g.T isa Vector)
        aux = [off_val; 0; A; A; 0]
	else
		aux = [off_val; 0; A; 0]
	end
    # If no signal, then don't show samples
	if sum(abs.(g.A)) == 0
		aux *= off_val
	end
	aux
end

get_theo_t(g::RF, key::Symbol,) = begin
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
	aux = [0; g.delay; t; g.delay+sum(g.T)]
	aux
end

get_theo_RF(seq, key::Symbol) = begin
	N = length(seq)
	T0 = get_block_start_times(seq)
	t = vcat([get_theo_t(seq.RF[i], key) .+ T0[i] for i=1:N]...)
	A = vcat([get_theo_A(seq.RF[i], key) for i=1:N]...) #; off_val=0 <---potential solution
	#Removing duplicated points
	#TODO: do this properly. As it is now it generates a bug for slew rates that are too high
	# mask = (G .== 0) #<---potential solution
	# t = t[mask]
	# G = G[mask]
	Interpolations.deduplicate_knots!(t; move_knots=true)
	return (t, A)
end
