#!/bin/sh   
source /opt/intel/oneapi/setvars.sh --force  
ulimit -s unlimited
/work/home/mayuan/software/qe-7.4.1/bin/plotband.x <phonobanddata.in> phonobanddata.out 
