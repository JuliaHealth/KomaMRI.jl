using KomaMRICore
obj = brain_phantom2D();
sys = Scanner();
seq = read_seq(joinpath(dirname(pathof(KomaMRICore)), "../../examples/3.koma_paper/comparison_accuracy/sequences/EPI/epi_100x100_TE100_FOV230.seq"));
simParams = KomaMRICore.default_sim_params();
raw = simulate(obj, seq, sys; simParams);

using KomaMRICore, CUDA
obj = brain_phantom2D();
sys = Scanner();
seq = read_seq(joinpath(dirname(pathof(KomaMRICore)), "../../examples/3.koma_paper/comparison_accuracy/sequences/EPI/epi_100x100_TE100_FOV230.seq"));
simParams = KomaMRICore.default_sim_params();
raw = simulate(obj, seq, sys; simParams);