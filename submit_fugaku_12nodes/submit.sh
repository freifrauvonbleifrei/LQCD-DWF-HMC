#!/bin/bash
#PJM -L rscgrp=small
#PJM -L node=12:torus
#PJM -L elapse=0:10:00
#PJM -g ra000001
#PJM --mpi proc=48
#PJM -j
#PJM -S
#PJM --rsc-list "freq=2200,eco_state=2"

date
echo "-------- printing this script ----------------"
cat $0
echo "-------- printing this script, done. ---------"
echo
echo "-------- printing main.yaml ------------------"
cat main.yaml
echo "-------- printing main.yaml, done ------------"
echo
echo "---------printing HMC_DWF_Nf2p1_eo.yaml ------"
cat HMC_DWF_Nf2p1_eo.yaml
echo "---------printing HMC_DWF_Nf2p1_eo.yaml ------"

module list

BIN=../benchmark_hmc_dwf/hmc_dwf_eo.elf
llio_transfer $BIN

export OMP_NUM_THREADS=12
export FLIB_BARRIER=HARD

mpiexec  $BIN

