using Koma
obj = read_phantom_jemris("column1d_phantom.h5")
seq = read_seq("epi_100x100_TE100_FOV230.seq")
sys = Scanner()
#Time
signal = simulate(obj,seq[1:2],sys) #warmup
signal = simulate(obj,seq,sys)
#Output
nothing
