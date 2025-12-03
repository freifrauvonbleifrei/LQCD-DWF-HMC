# if compilied with use_mpi=no
#../benchmark_hmc_dwf/hmc_dwf_eo.elf >& run.log
mpirun -np 1 ../benchmark_hmc_dwf/hmc_dwf_eo.elf >& run.log
