#!/bin/sh                           
#SBATCH  --job-name=mayuan
#SBATCH  --output=log.out                       
#SBATCH  --error=log.err                       
#SBATCH  --partition=public
#SBATCH  --nodes=1
#SBATCH  --ntasks=48
#SBATCH  --ntasks-per-node=48                          
#SBATCH  --cpus-per-task=1                         
#SBATCH  --exclude=node49,node98

source /work/home/mayuan/intel/oneapi/setvars.sh --force      
ulimit -s unlimited
export I_MPI_ADJUST_REDUCE=3
export MPIR_CVAR_COLL_ALIAS_CHECK=0
mpirun -np 48 /work/home/mayuan/software/qe-7.4.1/bin/pw.x -npool 4 <relax.in> relax.out 
wait
Nat=$(grep 'number of k points' -B 2 relax.out |head -n 1|awk {'print($1)'})
StruLine=$(expr $Nat + 5)
grep 'CELL_' -A $StruLine relax.out |tail -n `expr $StruLine + 1` > new_structure.out
sed  -i '/^$/d' new_structure.out




cell_parameters=$(grep -n 'CELL_PARAMETERS' scffit.in | awk -F : '{print $1}')
k_points=$(grep -n 'K_POINTS' scffit.in | awk -F : '{print $1}')
insert_position=$(expr $cell_parameters - 1)
stop_delete_position=$(expr $k_points - 1)
sed -i "${cell_parameters}, ${stop_delete_position}d" scffit.in
sed -i "${insert_position}r new_structure.out" scffit.in
mpirun -np 48 /work/home/mayuan/software/qe-7.4.1/bin/pw.x -npool 4 <scffit.in> scffit.out 




cell_parameters=$(grep -n 'CELL_PARAMETERS' scf.in | awk -F : '{print $1}')
k_points=$(grep -n 'K_POINTS' scf.in | awk -F : '{print $1}')
insert_position=$(expr $cell_parameters - 1)
stop_delete_position=$(expr $k_points - 1)
sed -i "${cell_parameters}, ${stop_delete_position}d" scf.in
sed -i "${insert_position}r new_structure.out" scf.in
mpirun -np 48 /work/home/mayuan/software/qe-7.4.1/bin/pw.x -npool 4 <scf.in> scf.out 
