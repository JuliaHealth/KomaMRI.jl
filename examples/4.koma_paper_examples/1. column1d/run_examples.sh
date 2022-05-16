#!/bin/bash
mpirun -np 4 pjemris ./jemris/run.xml > out_jemris
julia run_koma.jl > out_koma
