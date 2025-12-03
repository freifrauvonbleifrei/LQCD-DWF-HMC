# Bridge++: LQCD-DWF-HMC

2025.12.3  Issaku Kanamori (kanamori-i@riken.jp)

## Description

Benchmark of HMC for Lattice QCD with Domainwall type fermion.
It is based on Bridge++ https://bridge.kek.jp/Lattice-code/ 
but some local extensions have been applied.
The version of the base Bridge++ is 2.1.0.

The corelib version does not work temporary [2025.12.3].

This benchmark program runs 2+1 flavor HMC.
  * 1 flavor part: Rational HMC with even-odd decomposed fermion operator, for which multi shift CG solver (double prec.) is applied.
  * 2 flavor part: HMC with even-odd decomposed fermion operator, for which mixed precision solver (double + single prec.) with Richardson iterations.
  * All fermion operators use stout smeared gauge field.

It supports three targets:
  * Fugaku version, uses ACLE (alt_qxs)
  * Fugaku version with ACLE emulation (alt_qxs with use_qxs_arch=general)
  * general version (corelib, a different data layout and does not use mixed prec.)
  * GPU version with OpenAcc (alt_accel)

License: GPLv3.

## Requirements

N/A

## Build and Run

edit the `Makefile` depending on the target. There are several samples so tha one can just copy it to `Makefile`:
* `Makefile_fugaku`:  uses ACLE and FJMPI,  the clang mode of Fujitsu compiler
* `Makefile_fugaku_general`: w/o ACLE or FJMPI but uses a data layout for Fugaku, the clang mode of Fujitsu compiler
* `Makefile_x86_acle_general`: the same as Makefile_fugaku_general but the target is x86 with gcc.
* `Makefile_corelib`: w/o data layout for Fugaku.
* `Makefile_openacc`: GPU version with OpenACC.

### important switches in the Makefile

* target:          specifies the compiler. See Makefile_target.inc for the details of each target
* use_mpi:         to use Fujitsu extension of the persistent communication, set fjmpi
* use_alternative: set yes to use architecture specific tuning
* use_alt_qxs:     set yes to use A64FX (Fugaku) specific tuning
* use_qxs_arch:    set general to emulate ACLE implementation
* use_alt_accel:   set yes to use GPU vesion

### How to make

```
# make library
make
# make the benchmark executable (Fugaku or GPU version)
cd benchmark_hmc_dwf
make
cd -

# make the benchmark executable (w/o using data layout change)
cd benchmark_hmc_dwf_corelib
make
```

## Execution samples

see `submit_*` and `sample_x86*` directories.  In short, just
```
mpiexec hmc_dwf_eo.elf

```


To change the target problem size, edit main.yaml
```
  lattice_size             : [32,8,8,12]
  grid_size                : [1,1,2,2]
  number_of_thread         : 12
```
from the above, total 4-dim lattice size, MPI process division, number of OMP threads for each process.
If `use_alt_qxs=yes`, `lattce_size` divided by `grid_size` must be dividable by [8,4,1,1].
The Flop count in the log may not be working.


### check of the result

`H(diff)` at the almost end of the log should be a small number, typically in order of 0.1--0.001 (it can be negative).  It becomes larger as the problem size becomes larger.  If the same parameters (including the problem size) is used as the sample provided, the same value of `H(diff)` should be reproduced.
Note that it is a difference of `H_total` (~10^6 in lattice_size = [8,8,8,8] case), so that the last digit(s) of `H(diff)` may differ from the log there.

If `H(diff)` is too large, e.g. 1 or larger, most likely something is wrong.
Another non-trivial check is change a parameter in `HMC_DWF_Nf2p1_eo.yaml`
```
  number_of_steps   : [2,1,3]
```
to
```
  number_of_steps   : [4,1,3]
```
that is, increase the first number (let us call it `N`).  The execution time is roughly proportional to `N`, and `H(diff)` should decrees as `1/N^2`.  If it does not decrees as `1/N^2` (typically stays in constant), the results are incorrect.
