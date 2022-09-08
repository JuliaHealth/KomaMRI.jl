"""Divides a list of indices 1:N in k groups"""
function kfoldperm(N,k; type="random", breaks=[])
	k = min(N,k)
	n,r = divrem(N, k) #N >= k, N < k
	b = collect(1:n:N+1)
	Nb = length(b)
	for i in 1:Nb
		b[i] += i > r ? r : i-1
	end
	b = sort(unique(append!(b, breaks)))
	Nbreaks = length(b) - Nb
	if type=="random"
		p = randperm(N)
	elseif type=="ordered"
		p = 1:N
	end
	[p[r] for r in [b[i]:b[i+1]-1 for i=1:k+Nbreaks]] #TODO: use RF starts and ends differently to remove PATCH in run_sim_time_iter
end

"""Calcualtes key time points that can not be missed during the simulation"""
function points_from_key_times(times ; dt)
	t = Float64[]
	for i = 1:length(times)-1
		if dt < times[i+1] - times[i]
			taux = collect(range(times[i],times[i+1];step=dt))
		else
			taux = [times[i], times[i+1]]
		end
		append!(t, taux)
	end
	t
end
"""
Uniform time-step calculation
"""
function get_uniform_times(seq,Δt;Δt_rf=1e-4)
	t, Δt = get_variable_times(seq; dt=Δt, dt_rf=Δt_rf)
end

"""
Variable time-step calculation
"""
function get_variable_times(seq; dt=1, dt_rf=1e-4)
	t = Float64[]
	ΔT = durs(seq) #Duration of sequence block
	T0 = cumsum([0; ΔT[:]]) #Start time of each block
	for i = 1:length(seq)
		s = seq[i] #Current sequence block
		t0 = T0[i]
		if is_ADC_on(s)
			ts = get_sample_times(s) .+ t0
			taux = points_from_key_times(ts; dt) # ADC sampling
			append!(t, taux)
		end
		if is_RF_on(s)
			y = s.RF[1]
			delay, T = y.delay, y.T
			t1 = t0 + delay
			t2 = t1 + sum(T)
			rf0 = t0 + get_RF_center(y) #get_RF_center includes delays
			taux = points_from_key_times([t1,rf0,t2]; dt=dt_rf) # Arbitrary RF
			append!(t, taux)
		end
		if is_GR_on(s)
			active_gradients = []
			if is_Gx_on(s) append!(active_gradients, s.GR.x) end
			if is_Gy_on(s) append!(active_gradients, s.GR.y) end
			if is_Gz_on(s) append!(active_gradients, s.GR.z) end
			for y = active_gradients
				ts = get_theo_t(y) .+ t0
				taux = points_from_key_times(ts; dt)
				append!(t, taux)
			end
		end
	end
	t = sort(unique(x -> round(x*1e8), t)) #Removing repeated points
	Δt = t[2:end] .- t[1:end-1]
	t = t[1:end-1]
	t, Δt
end

"""Calculate RF key points (start-end) to split simulation in RF and non-RF parts. """
function get_breaks_in_RF_key_points(seq::Sequence, t)
	T0 = cumsum([0; durs(seq)[:]])
	# Identifying RF key points
	key_points = Float64[]
	key_idxs = Int[]
	for (i, s) = enumerate(seq)
		if is_RF_on(s)
			t0 = T0[i] + s.RF.delay[1]	#start of RF waverform
			tf = T0[i] + s.RF.dur[1]	#end of RF waveform
			append!(key_points, [t0; tf])
			idx0 = argmin(abs.(t.-t0))
			idxf = argmin(abs.(t.-tf))
			append!(key_idxs,   [idx0; idxf])
		end
	end
	key_idxs
end