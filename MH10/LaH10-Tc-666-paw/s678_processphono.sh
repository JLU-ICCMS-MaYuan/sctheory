#!/bin/sh                           
#SBATCH  --job-name=mayuan
#SBATCH  --output=log.out                       
#SBATCH  --error=log.err                       
#SBATCH  --partition=public
#SBATCH  --nodes=1
#SBATCH  --ntasks=48
#SBATCH  --ntasks-per-node=48                          
#SBATCH  --cpus-per-task=1                         
#SBATCH  --exclude=node98

source /work/home/mayuan/intel/oneapi/setvars.sh --force      
ulimit -s unlimited
export I_MPI_ADJUST_REDUCE=3
export MPIR_CVAR_COLL_ALIAS_CHECK=0
#echo "run q2r"                   
#mpirun -np 48 /work/home/mayuan/software/qe-7.1/bin/q2r.x    -npool 4 <q2r.in> q2r.out 
#grep nqs q2r.out > nqs           
#echo "run phonoband"             
mpirun -np 48 /work/home/mayuan/software/qe-7.1/bin/matdyn.x -npool 4 <matdyn.in> matdyn.out 
#mkdir -p phonoband               
#cp La1H10.freq     phonoband         
#cp La1H10.freq.gp  phonoband         
#cp La1H10.dyn0     phonoband         
#cp gam.lines   phonoband         
#cp elph.gamma* phonoband         
#echo "run phonoband-data-process"
#/work/home/mayuan/software/qe-7.1/bin/plotband.x <phonobanddata.in> phonobanddata.out 
#echo "run phonodos"              
#mpirun -np 48 /work/home/mayuan/software/qe-7.1/bin/matdyn.x  <phonodos.in> phonodos.out  
#mkdir -p phonodos                
#cp La1H10"_phono.dos"  phonodos      
