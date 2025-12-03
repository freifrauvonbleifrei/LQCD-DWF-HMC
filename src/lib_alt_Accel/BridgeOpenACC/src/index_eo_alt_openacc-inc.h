/*!
      @file    index_eo_alt_openacc-inc.h
      @brief
      @author  Hideo Matsufuru (matufuru)
               $LastChangedBy: matufuru $
      @date    $LastChangedDate:: 2025-08-20 16:24:03 #$
      @version $LastChangedRevision: 2653 $
*/

#include "lib_alt_Accel/inline/define_params.h"
#include "lib_alt_Accel/inline/define_index.h"

//====================================================================
void split(real_t* RESTRICT ve, real_t* RESTRICT vo,
           real_t* RESTRICT w,
           int ieo_origin, int nin, int* Nsize)
{
  int Nx   = Nsize[0];
  int Ny   = Nsize[1];
  int Nz   = Nsize[2];
  int Nt   = Nsize[3];

  int Nvol  = Nx * Ny * Nz * Nt;
  int Nvol2 = Nvol/2;

  int Nvol_pad  = CEIL_NWP(Nvol);
  int Nvol2_pad = CEIL_NWP(Nvol2);
  
  int nv  = nin * Nvol_pad;
  int nv2 = nin * Nvol2_pad;

#pragma acc data present(ve[0:nv2], vo[0:nv2], w[0:nv]), \
                 copyin(ieo_origin, nin, Nx, Ny, Nz, Nt, Nvol, Nvol2_pad)
 {

#pragma acc parallel num_workers(NUM_WORKERS) vector_length(VECTOR_LENGTH)
  {
#pragma acc loop gang worker vector
    for(int ist2 = Nvol2_pad - NWP; ist2 < Nvol2_pad; ++ist2){
      for(int in = 0; in < nin; ++in){
        ve[IDX2(nin, in, ist2)] = 0.0;
        vo[IDX2(nin, in, ist2)] = 0.0;
      }
    }
  }

#pragma acc parallel num_workers(NUM_WORKERS) vector_length(VECTOR_LENGTH)
  {
#pragma acc loop gang worker vector
    for(int ist = 0; ist < Nvol; ++ist){
      int ix = ist % Nx;
      int iy = (ist/Nx) % Ny;
      int iz = (ist/(Nx * Ny)) % Nz;
      int it = ist/(Nx * Ny * Nz);
      int ieo = (ix + iy + iz + it + ieo_origin) % 2;
      int ist2 = ist/2;
      if(ieo == 0){
        for(int in = 0; in < nin; ++in){
          ve[IDX2(nin, in, ist2)] = w[IDX2(nin, in, ist)];
        }
      }else{
        for(int in = 0; in < nin; ++in){
          vo[IDX2(nin, in, ist2)] = w[IDX2(nin, in, ist)];
        }
      }
    }
  }

 } // acc data

}

//====================================================================
void merge(real_t* RESTRICT v,
           real_t* RESTRICT we, real_t* RESTRICT wo,
           int ieo_origin, int nin, int* Nsize)
{
  int Nx   = Nsize[0];
  int Ny   = Nsize[1];
  int Nz   = Nsize[2];
  int Nt   = Nsize[3];

  int Nvol  = Nx * Ny * Nz * Nt;
  int Nvol2 = Nvol/2;

  int Nvol_pad  = CEIL_NWP(Nvol);
  int Nvol2_pad = CEIL_NWP(Nvol2);

  int nv  = nin * Nvol_pad;
  int nv2 = nin * Nvol2_pad;

#pragma acc data present(we[0:nv2], wo[0:nv2], v[0:nv]), \
                 copyin(ieo_origin, nin, Nx, Ny, Nz, Nt, Nvol)
#pragma acc parallel num_workers(NUM_WORKERS) vector_length(VECTOR_LENGTH)
 {

#pragma acc loop gang worker vector
  for(int ist = 0; ist < Nvol_pad; ++ist){
    int ix = ist % Nx;
    int iy = (ist/Nx) % Ny;
    int iz = (ist/(Nx * Ny)) % Nz;
    int it = ist/(Nx * Ny * Nz);
    int ieo = (ix + iy + iz + it + ieo_origin) % 2;
    int ist2 = ist/2;
    if(ist < Nvol){
      if(ieo == 0){
        for(int in = 0; in < nin; ++in){
          v[IDX2(nin, in, ist)] = we[IDX2(nin, in, ist2)];
        }
      }else{
        for(int in = 0; in < nin; ++in){
          v[IDX2(nin, in, ist)] = wo[IDX2(nin, in, ist2)];
        }
      }
    }else{
      for(int in = 0; in < nin; ++in){
        v[IDX2(nin, in, ist)] = 0.0;
      }
    }
  }
 }

}

//============================================================END=====
