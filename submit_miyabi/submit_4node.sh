#!/bin/bash
#PBS -q debug-g
#PBS -l select=4:ompthreads=2
#PBS -l walltime=00:30:00
#PBS -W group_list=xg24i061
#PBS -j oe

module list

cd ${PBS_O_WORKDIR}

#echo "--------- this script ------------------"
#cat $0
#echo "--------- this script, done ------------"

export OMP_NUM_THREADS=2
echo OMP_NUM_THREAD: $OMP_NUM_THREADS
export OMP_PROC_BIND=spread

BIN=../benchmark_hmc_dwf/hmc_dwf_eo.elf

ls -l $BIN

mpirun $BIN
