/*!
      @file    mult_Doainwall_5din_eo_4d_openacc-inc.h
      @brief
      @author  Hideo Matsufuru (matufuru)
               $LastChangedBy: matufuru $
      @date    $LastChangedDate: 2013-01-22 13:51:53 #$
      @version $LastChangedRevision: 2653 $
*/

#ifndef MULT_DOMAINWALL_5DIN_EO_4D_ACC_INCLUDED
#define MULT_DOMAINWALL_5DIN_EO_4D_ACC_INCLUDED

//====================================================================
void mult_domainwall_5din_eo_hopb_dirac_4d(
       real_t *RESTRICT vp, real_t *RESTRICT up, real_t *RESTRICT wp,
       int Ns, int *bc, int *Nsize, int *do_comm, int ieo, int jeo,
       int jgm5)
{
  int Nx  = Nsize[0];
  int Ny  = Nsize[1];
  int Nz  = Nsize[2];
  int Nt  = Nsize[3];
  int Nst = Nx * Ny * Nz * Nt;
  int Nst_pad = CEIL_NWP(Nst);

  int Nin5 = NVCD * Ns;

  int size   = Nin5 * Nst_pad;
  int size_u = NDF  * Nst_pad * 2 * NDIM;

#pragma acc data present(vp[0:size], up[0:size_u], wp[0:size]) \
                 copyin(Nst, Nst_pad, Nx, Ny, Nz, Nt, Nin5, Ns, \
                        bc[0:NDIM], do_comm[0:NDIM], ieo, jeo, jgm5)
 {
#pragma acc parallel num_workers(NUM_WORKERS) vector_length(VECTOR_LENGTH)
 {

#pragma acc loop gang worker vector
  for (int site = 0; site < Nst_pad; ++site) {
   if(site < Nst){
    int Nxy  = Nx  * Ny;
    int Nxyz = Nxy * Nz;

    int ix   = site % Nx;
    int iyzt = site / Nx;
    int ixy  = site % Nxy;
    int iy   = iyzt % Ny;
    int izt  = site / Nxy;
    int iz   = izt % Nz;
    int it   = izt / Nz;
    int ixyz = site % Nxyz;

    int keo  = (jeo + iy + iz + it) % 2;

    int idir;

    for(int is = 0; is < Ns; ++is){

      real_t vL[NVCD];

      for(int ivcd = 0; ivcd < NVCD; ++ivcd){
        vL[ivcd] = 0.0;
      }

      idir = 0;

      // mult_xp
      {
        int ix2 = (ix + keo) % Nx;
        int nei = ix2 + Nx * iyzt;
        real_t bc2 = 1.0;
        if(ix == Nx-1 && keo == 1) bc2 = bc[0];

        real_t ut[NDF];
        load_u(ut, up, site + Nst_pad * (ieo + 2 * idir));

        real_t wt[NVCD];
        if(jgm5 == 0){
          for(int ivcd = 0; ivcd < NVCD; ++ivcd){
            wt[ivcd] = bc2 * wp[IDX2(Nin5, (ivcd + NVCD * is), nei)];
          }
        }else{
          for(int ivc = 0; ivc < NVC; ++ivc){
            wt[ivc+ID1] = bc2 * wp[IDX2(Nin5, (ivc+ID3 + NVCD * is), nei)];
            wt[ivc+ID2] = bc2 * wp[IDX2(Nin5, (ivc+ID4 + NVCD * is), nei)];
            wt[ivc+ID3] = bc2 * wp[IDX2(Nin5, (ivc+ID1 + NVCD * is), nei)];
            wt[ivc+ID4] = bc2 * wp[IDX2(Nin5, (ivc+ID2 + NVCD * is), nei)];
          }
        }
        mult_wilson_xpb(vL, ut, wt);
      }
      // mult_xm
      {
        int ix2 = (ix - 1 + keo + Nx) % Nx;
        int nei  = ix2 + Nx * iyzt;
        real_t bc2 = 1.0;
        if(ix == 0 && keo == 0) bc2 = bc[0];

        real_t ut[NDF];
        load_u(ut, up, nei + Nst_pad * (1-ieo + 2 * idir));

        real_t wt[NVCD];
        if(jgm5 == 0){
          for(int ivcd = 0; ivcd < NVCD; ++ivcd){
            wt[ivcd] = bc2 * wp[IDX2(Nin5, (ivcd + NVCD * is), nei)];
          }
        }else{
          for(int ivc = 0; ivc < NVC; ++ivc){
            wt[ivc+ID1] = bc2 * wp[IDX2(Nin5, (ivc+ID3 + NVCD * is), nei)];
            wt[ivc+ID2] = bc2 * wp[IDX2(Nin5, (ivc+ID4 + NVCD * is), nei)];
            wt[ivc+ID3] = bc2 * wp[IDX2(Nin5, (ivc+ID1 + NVCD * is), nei)];
            wt[ivc+ID4] = bc2 * wp[IDX2(Nin5, (ivc+ID2 + NVCD * is), nei)];
          }
        }
        mult_wilson_xmb(vL, ut, wt);
      }

      idir = 1;

      if ((iy < Ny-1) || (do_comm[idir] == 0)) {
        int iy2 = (iy + 1) % Ny;
        int nei = ix + Nx * (iy2 + Ny * izt);
        real_t bc2 = 1.0;
        if(iy == Ny-1) bc2 = bc[1];

        real_t ut[NDF];
        load_u(ut, up, site + Nst_pad * (ieo + 2 * idir));

        real_t wt[NVCD];
        if(jgm5 == 0){
          for(int ivcd = 0; ivcd < NVCD; ++ivcd){
            wt[ivcd] = bc2 * wp[IDX2(Nin5, (ivcd + NVCD * is), nei)];
          }
        }else{
          for(int ivc = 0; ivc < NVC; ++ivc){
            wt[ivc+ID1] = bc2 * wp[IDX2(Nin5, (ivc+ID3 + NVCD * is), nei)];
            wt[ivc+ID2] = bc2 * wp[IDX2(Nin5, (ivc+ID4 + NVCD * is), nei)];
            wt[ivc+ID3] = bc2 * wp[IDX2(Nin5, (ivc+ID1 + NVCD * is), nei)];
            wt[ivc+ID4] = bc2 * wp[IDX2(Nin5, (ivc+ID2 + NVCD * is), nei)];
          }
        }
        mult_wilson_ypb(vL, ut, wt);
      }

      if ((iy > 0) || (do_comm[idir] == 0)) {
        int iy2 = (iy - 1 + Ny) % Ny;
        int nei = ix + Nx * (iy2 + Ny * izt);
        real_t bc2 = 1.0;
        if(iy == 0) bc2 = bc[1];

        real_t ut[NDF];
        load_u(ut, up, nei + Nst_pad * (1-ieo + 2 * idir));

        real_t wt[NVCD];
        if(jgm5 == 0){
          for(int ivcd = 0; ivcd < NVCD; ++ivcd){
            wt[ivcd] = bc2 * wp[IDX2(Nin5, (ivcd + NVCD * is), nei)];
          }
        }else{
          for(int ivc = 0; ivc < NVC; ++ivc){
            wt[ivc+ID1] = bc2 * wp[IDX2(Nin5, (ivc+ID3 + NVCD * is), nei)];
            wt[ivc+ID2] = bc2 * wp[IDX2(Nin5, (ivc+ID4 + NVCD * is), nei)];
            wt[ivc+ID3] = bc2 * wp[IDX2(Nin5, (ivc+ID1 + NVCD * is), nei)];
            wt[ivc+ID4] = bc2 * wp[IDX2(Nin5, (ivc+ID2 + NVCD * is), nei)];
          }
        }
        mult_wilson_ymb(vL, ut, wt);
      }

      idir = 2;

      if ((iz < Nz-1) || (do_comm[idir] == 0)) {
        int iz2 = (iz + 1) % Nz;
        int nei = ixy + Nxy * (iz2 + Nz * it);
        real_t bc2 = 1.0;
        if(iz == Nz-1) bc2 = bc[2];

        real_t ut[NDF];
        load_u(ut, up, site + Nst_pad * (ieo + 2 * idir));

        real_t wt[NVCD];
        if(jgm5 == 0){
          for(int ivcd = 0; ivcd < NVCD; ++ivcd){
            wt[ivcd] = bc2 * wp[IDX2(Nin5, (ivcd + NVCD * is), nei)];
          }
        }else{
          for(int ivc = 0; ivc < NVC; ++ivc){
            wt[ivc+ID1] = bc2 * wp[IDX2(Nin5, (ivc+ID3 + NVCD * is), nei)];
            wt[ivc+ID2] = bc2 * wp[IDX2(Nin5, (ivc+ID4 + NVCD * is), nei)];
            wt[ivc+ID3] = bc2 * wp[IDX2(Nin5, (ivc+ID1 + NVCD * is), nei)];
            wt[ivc+ID4] = bc2 * wp[IDX2(Nin5, (ivc+ID2 + NVCD * is), nei)];
          }
        }
        mult_wilson_zpb(vL, ut, wt);
      }

      if ((iz > 0) || (do_comm[idir] == 0)) {
        int iz2 = (iz - 1 + Nz) % Nz;
        int nei = ixy + Nxy * (iz2 + Nz * it);
        real_t bc2 = 1.0;
        if(iz == 0) bc2 = bc[2];

        real_t ut[NDF];
        load_u(ut, up, nei + Nst_pad * (1-ieo + 2 * idir));

        real_t wt[NVCD];
        if(jgm5 == 0){
          for(int ivcd = 0; ivcd < NVCD; ++ivcd){
            wt[ivcd] = bc2 * wp[IDX2(Nin5, (ivcd + NVCD * is), nei)];
          }
        }else{
          for(int ivc = 0; ivc < NVC; ++ivc){
            wt[ivc+ID1] = bc2 * wp[IDX2(Nin5, (ivc+ID3 + NVCD * is), nei)];
            wt[ivc+ID2] = bc2 * wp[IDX2(Nin5, (ivc+ID4 + NVCD * is), nei)];
            wt[ivc+ID3] = bc2 * wp[IDX2(Nin5, (ivc+ID1 + NVCD * is), nei)];
            wt[ivc+ID4] = bc2 * wp[IDX2(Nin5, (ivc+ID2 + NVCD * is), nei)];
          }
        }
        mult_wilson_zmb(vL, ut, wt);
      }

      idir = 3;

      if ((it < Nt-1) || (do_comm[idir] == 0)) {
        int it2 = (it + 1) % Nt;
        int nei = ixyz + Nxyz * it2;
        real_t bc2 = 1.0;
        if(it == Nt-1) bc2 = bc[3];

        real_t ut[NDF];
        load_u(ut, up, site + Nst_pad * (ieo + 2 * idir));

        real_t wt[NVCD];
        if(jgm5 == 0){
          for(int ivcd = 0; ivcd < NVCD; ++ivcd){
            wt[ivcd] = bc2 * wp[IDX2(Nin5, (ivcd + NVCD * is), nei)];
          }
        }else{
          for(int ivc = 0; ivc < NVC; ++ivc){
            wt[ivc+ID1] = bc2 * wp[IDX2(Nin5, (ivc+ID3 + NVCD * is), nei)];
            wt[ivc+ID2] = bc2 * wp[IDX2(Nin5, (ivc+ID4 + NVCD * is), nei)];
            wt[ivc+ID3] = bc2 * wp[IDX2(Nin5, (ivc+ID1 + NVCD * is), nei)];
            wt[ivc+ID4] = bc2 * wp[IDX2(Nin5, (ivc+ID2 + NVCD * is), nei)];
          }
        }
        mult_wilson_tpb_dirac(vL, ut, wt);
      }

      if ((it > 0) || (do_comm[idir] == 0)) {
        int it2 = (it - 1 + Nt) % Nt;
        int nei = ixyz + Nxyz * it2;
        real_t bc2 = 1.0;
        if(it == 0) bc2 = bc[3];

        real_t ut[NDF];
        load_u(ut, up, nei + Nst_pad * (1-ieo + 2 * idir));

        real_t wt[NVCD];
        if(jgm5 == 0){
          for(int ivcd = 0; ivcd < NVCD; ++ivcd){
            wt[ivcd] = bc2 * wp[IDX2(Nin5, (ivcd + NVCD * is), nei)];
          }
        }else{
          for(int ivc = 0; ivc < NVC; ++ivc){
            wt[ivc+ID1] = bc2 * wp[IDX2(Nin5, (ivc+ID3 + NVCD * is), nei)];
            wt[ivc+ID2] = bc2 * wp[IDX2(Nin5, (ivc+ID4 + NVCD * is), nei)];
            wt[ivc+ID3] = bc2 * wp[IDX2(Nin5, (ivc+ID1 + NVCD * is), nei)];
            wt[ivc+ID4] = bc2 * wp[IDX2(Nin5, (ivc+ID2 + NVCD * is), nei)];
          }
        }
        mult_wilson_tmb_dirac(vL, ut, wt);
      }

      for(int ivcd = 0; ivcd < NVCD; ++ivcd){
        vp[IDX2(Nin5, (ivcd + NVCD * is), site)] = vL[ivcd];
      }

    } // is loop

   }
  } // site loop

 }  // acc parallel
 }  // acc data

}

#endif
//============================================================END=====
