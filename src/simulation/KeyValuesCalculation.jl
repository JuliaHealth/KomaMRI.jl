"""
    get_theo_A(g::Grad; off_val=0)
    get_theo_A(r::RF; off_val=0, max_rf_samples=Inf)
    get_theo_A(d::ADC; off_val=0)

Get the theoretical amplitude of a signal.

# Arguments
- `g::Grad`: the gradient
- `r::RF`: the RF signal
- `d::ADC`: the ADC

# Keywords
- `off_val`: offset value
- `max_rf_samples`: number of maximum samples for RF signal

# Returns
- `aux`: vector with information about the amplitude of the input signal
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

get_theo_A(r::RF; off_val=0, max_rf_samples=Inf) = begin
	A = [r.A r.A]'; A = A[:]
	if sum(abs.(r.A)) == 0
		aux = off_val * [1; 1; ones(length(A)); 1]
	else
		#Calculating frequency waveform
		if length(r.A) != 1
			tf = range(0,sum(r.T),length(r.A))
		else
			tf = [0]
		end
		freq = exp.(1im*2π*r.Δf.*tf)
		freq = transpose([freq freq])[:]
		A = A .* freq
		A[abs.(A).<=1e-14] .= 0 #Remove small values
		#Output
		aux = [off_val; 0; A; 0]
	end
	#Subsample
	if length(aux) > max_rf_samples
		n = floor(Int, length(aux) / max_rf_samples)
		aux[1:n:end]
	else
		aux
	end
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
    get_theo_t(g::Grad)
    get_theo_t(r::RF; max_rf_samples=Inf)
    get_theo_t(d::ADC)

Get the theoretical time of a signal.

# Arguments
- `g::Grad`: the gradient
- `r::RF`: the RF signal
- `d::ADC`: the ADC

# Keywords
- `max_rf_samples`: number of maximum samples for RF signal

# Returns
- `aux`: vector with information about the time of the input signal
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
		t[1:n:end]
	else
		t
	end
end

get_theo_t(d::ADC) = begin
	t = [d.delay; d.delay+d.T]
	t = [t t]'; t = t[1:end-1]
	[0; t; d.delay+d.T]
end


"""
    get_theo_Gi(seq, idx)

Get the theoretical gradient for a sequence in a defined axis.

# Arguments
- `seq`: the sequence
- `idx`: the axis

# Returns
- `(t, G)`: tuple:
    - `t`: the time
    - `G`: the theoretical gradient
"""
get_theo_Gi(seq, idx) = begin
	ΔT, N = durs(seq), length(seq)
	T0 = cumsum([0; ΔT], dims=1)
	t = vcat([get_theo_t(seq.GR[idx,i]) .+ T0[i] for i=1:N]...)
	G = vcat([get_theo_A(seq.GR[idx,i]) for i=1:N]...)
	#Removing duplicated points
	#TODO: do this properly. As it is now it generates a bug for slew rates that are too high
	Interpolations.deduplicate_knots!(t; move_knots=true)
	return (t, G)
end
