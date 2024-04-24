#!/bin/bash
if [ "$3" == "c" ] || [ -z "$3" ]
then

    if [ "$1" == "jemris" ] || [ -z "$1" ]
    then
        echo "### JEMRIS ###"
        mpirun pjemris ./run_jemris.xml
        rm signals_ismrmrd.h5
    fi

    if [ "$1" == "koma" ] || [ -z "$1" ]
    then
        echo "### KOMA ###"
        julia $2 ./run_koma.jl
    fi

    echo "### END ###"

fi

