#!/bin/sh   
source /home/mayuan/intel/oneapi/setvars.sh
ulimit -s unlimited
mpirun -np 2 /home/mayuan/soft/qe-7.3.1/bin/matdyn.x -npool 1 <matdyn.in> matdyn.out 
