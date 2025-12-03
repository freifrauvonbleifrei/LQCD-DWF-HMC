#!/bin/bash
#PJM -L node=1
#PJM -L elapse=0:10:00
#PJM -L rscgrp=small
#PJM -g ra000001
#PJM --mpi proc=4
#PJM -j
#PJM -S

export OMP_NUM_THREADS=12
export FLIB_BARRIER=HARD

BIN=../fs_hmc_dwf/hmc_dwf_eo.elf

fapp -C -d ./rep1  -Inompi -Hevent=pa1  mpiexec $BIN alt_qxs
fapp -C -d ./rep2  -Inompi -Hevent=pa2  mpiexec $BIN alt_qxs
fapp -C -d ./rep3  -Inompi -Hevent=pa3  mpiexec $BIN alt_qxs
fapp -C -d ./rep4  -Inompi -Hevent=pa4  mpiexec $BIN alt_qxs
fapp -C -d ./rep5  -Inompi -Hevent=pa5  mpiexec $BIN alt_qxs
fapp -C -d ./rep6  -Inompi -Hevent=pa6  mpiexec $BIN alt_qxs
fapp -C -d ./rep7  -Inompi -Hevent=pa7  mpiexec $BIN alt_qxs
fapp -C -d ./rep8  -Inompi -Hevent=pa8  mpiexec $BIN alt_qxs
fapp -C -d ./rep9  -Inompi -Hevent=pa9  mpiexec $BIN alt_qxs
fapp -C -d ./rep10 -Inompi -Hevent=pa10 mpiexec $BIN alt_qxs
fapp -C -d ./rep11 -Inompi -Hevent=pa11 mpiexec $BIN alt_qxs
fapp -C -d ./rep12 -Inompi -Hevent=pa12 mpiexec $BIN alt_qxs
fapp -C -d ./rep13 -Inompi -Hevent=pa13 mpiexec $BIN alt_qxs
fapp -C -d ./rep14 -Inompi -Hevent=pa14 mpiexec $BIN alt_qxs
fapp -C -d ./rep15 -Inompi -Hevent=pa15 mpiexec $BIN alt_qxs
fapp -C -d ./rep16 -Inompi -Hevent=pa16 mpiexec $BIN alt_qxs
fapp -C -d ./rep17 -Inompi -Hevent=pa17 mpiexec $BIN alt_qxs

