using KomaMRI, Suppressor

obj = @suppress read_phantom_jemris("../../../2.phantoms/brain.h5")
seq = @suppress read_seq("../sequences/Spiral/spiral_100x100_FOV230_SPZ_INTER1.seq")
sys = Scanner()
#Time
warmup = @suppress simulate(obj,seq,sys) #warmup, to precompile
signal = simulate(obj,seq,sys)
