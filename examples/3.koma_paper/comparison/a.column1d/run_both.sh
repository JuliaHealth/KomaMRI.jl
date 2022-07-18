#!/bin/bash
echo "### JEMRIS ###"
mpirun pjemris ./run_jemris.xml
rm signals.h5 signals_ismrmrd.h5
echo "### KOMA ###"
julia ./run_koma.jl
echo "### END ###"

