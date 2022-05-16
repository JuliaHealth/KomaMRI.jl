#!/bin/bash
mpirun -np 4 pjemris run_jemris.xml > sim_jemris_koma_out
#julia run_koma.jl > sim_jemris_koma_out
