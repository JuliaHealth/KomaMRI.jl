using KomaMRI

obj = read_phantom_jemris("koma/column1d_phantom.h5")
seq = read_seq("koma/epi_100x100_TE100_FOV230.seq")
sys = Scanner()
#Time
warmup = simulate(obj,seq[1:8],sys) #warmup, to precompile
signal = simulate(obj,seq,sys)
