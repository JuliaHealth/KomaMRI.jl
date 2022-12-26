using KomaMRI, Suppressor, MAT

obj = @suppress read_phantom_jemris("../../../2.phantoms/brain.h5")
seq = @suppress read_seq("../sequences/Spiral/spiral_100x100_FOV230_SPZ_INTER1.seq")
sys = Scanner()
#Fix slight differences in the sequence (reading Pulseq 1.2.1 files is depricated)
adc_dwell_time_s = seq[2].ADC[1].T / seq[2].ADC[1].N
rf_dwell_time_s = seq.DEF["RadiofrequencyRasterTime"]
gr_dwell_time_s = seq.DEF["GradientRasterTime"]

seq[1].RF[1].delay = 0
seq[1].RF[1].T    += rf_dwell_time_s

seq[2].GR[1].rise  = 0
seq[2].GR[2].rise  = 0
seq[2].GR[1].delay = rf_dwell_time_s/2
seq[2].GR[2].delay = rf_dwell_time_s/2

seq[2].ADC[1].delay = rf_dwell_time_s/2 + adc_dwell_time_s/2
seq[2].ADC[1].T    += adc_dwell_time_s
#Time
simParams = Dict{String, Any}("return_type"=>"mat")
warmup = @suppress simulate(obj,seq,sys; simParams) #warmup, to precompile
signal = simulate(obj,seq,sys; simParams)
#Export
matwrite("./signal_koma.mat", Dict("signal" => signal ./ prod(size(obj))); compress = true)
