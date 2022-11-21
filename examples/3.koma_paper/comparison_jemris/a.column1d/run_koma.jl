using KomaMRI, Suppressor

obj = @suppress read_phantom_jemris("../../../2.phantoms/column1d.h5")
seq = @suppress read_seq("../sequences/EPI/epi_100x100_TE100_FOV230.seq")
sys = Scanner()
#Time
warmup = @suppress simulate(obj,seq,sys) #warmup, to precompile
signal = simulate(obj,seq,sys)
