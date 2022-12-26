#!/bin/bash
echo "### KOMA Multi-Shot Spiral CPU ###" 	| tee .out
julia $1 ./MRiLab_speed.jl cpu 0 8		| tee -a .out #CPU
echo "### KOMA Multi-Shot Spiral GPU0 ###" 	| tee -a .out
julia $1 ./MRiLab_speed.jl gpu 0 1		| tee -a .out #GPU0
echo "### KOMA Multi-Shot Spiral GPU1 ###" 	| tee -a .out
julia $1 ./MRiLab_speed.jl gpu 1 1		| tee -a .out #GPU1
