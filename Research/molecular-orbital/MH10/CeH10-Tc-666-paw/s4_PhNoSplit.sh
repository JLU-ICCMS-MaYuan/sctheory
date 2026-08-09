#!/bin/sh   
source /opt/intel/oneapi/setvars.sh --force  
ulimit -s unlimited
 /work/home/mayuan/software/qe-7.4.1/bin/ph.x -npool 1 <ph_no_split.in> ph_no_split.out 
