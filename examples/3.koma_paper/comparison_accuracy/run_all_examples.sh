#!/bin/bash
echo "### (A) Column 1D ###"                   | tee .out
cd a.column1d
./run_both.sh $1 $2 $3                         | tee -a ../.out
echo "### (B) Sphere CS ###"                   | tee -a ../.out
cd ../b.spheresCS/
./run_both.sh $1 $2 $3                         | tee -a ../.out
echo "### (C) Brain with Susceptibility ###"   | tee -a ../.out
cd ../c.brain_suscep/
./run_both.sh $1 $2 $3                         | tee -a ../.out
echo "### (D) Brain with Motion ###"           | tee -a ../.out
cd ../d.brain_motion/
./run_both.sh $1 $2 $3                         | tee -a ../.out
echo "### (E) Brain Spiral ###"                | tee -a ../.out
cd ../e.brain_spiral/
./run_both.sh $1 $2 $3                         | tee -a ../.out
cd ../
sed -i '/Preparation\|Simulating/d' .out #Removing unecessary lines in output file
