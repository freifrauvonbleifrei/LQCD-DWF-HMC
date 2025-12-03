/*!
      @file    mult_Doainwall_5din_LUinv_openacc-inc.h
      @brief
      @author  Hideo Matsufuru (matufuru)
               $LastChangedBy: matufuru $
      @date    $LastChangedDate: 2013-01-22 13:51:53 #$
      @version $LastChangedRevision: 2653 $
*/

#ifndef MULT_DOMAINWALL_5DIN_LUINV_ACC_INCLUDED
#define MULT_DOMAINWALL_5DIN_LUINV_ACC_INCLUDED

//====================================================================
void mult_domainwall_5din_LUinv_dirac(
                        real_t *RESTRICT vp, real_t *RESTRICT wp,
                        int Ns, int *Nsize,
                        real_t *e, real_t *f,
                        real_t *dpinv, real_t *dm, real_t alpha)
{
  int Nx  = Nsize[0];
  int Ny  = Nsize[1];
  int Nz  = Nsize[2];
  int Nt  = Nsize[3];
  int Nst = Nx * Ny * Nz * Nt;
  int Nst_pad = CEIL_NWP(Nst);

  int Nin5 = NVCD * Ns;
  int size = Nin5 * Nst_pad;

#pragma acc data present(vp[0:size], wp[0:size]) \
                 copyin(Nst, Nst_pad, Nin5, Ns, e[0:Ns-1], f[0:Ns-1], \
                        dpinv[0:Ns], dm[0:Ns], alpha)
 {

#pragma acc parallel num_workers(NUM_WORKERS) vector_length(VECTOR_LENGTH)
 {

#pragma acc loop gang worker vector
  for (int idx = 0; idx < NVC * Nst_pad; ++idx) {
    int idx2_wp = idx / NWP;
    int idx_in  = idx % NWP;
    int ivc = idx2_wp % NVC;
    int idx_out = idx2_wp / NVC;
    int site = idx_in + NWP * idx_out;
    if(site < Nst){

    real_t vt[ND], yt[ND], xt[ND];

    int is = 0;
    for(int id = 0; id < ND; ++id){
      int ivcd = ivc + NVC * id;
      vt[id] = wp[IDX2(Nin5, (ivcd + NVCD * is), site)];
    }
    for(int id = 0; id < ND; ++id){
      int ivcd = ivc + NVC * id;
      vp[IDX2(Nin5, (ivcd + NVCD * is), site)] = vt[id];
    }

    for(int id = 0; id < ND; ++id){
      yt[id] = e[0] * vt[id];
    }

    for (int is = 1; is < Ns-1; ++is) {

      for(int id = 0; id < ND; ++id){
        xt[id] = vt[id];
      }

      for(int id = 0; id < ND; ++id){
        int ivcd = ivc + NVC * id;
        vt[id] = wp[IDX2(Nin5, (ivcd + NVCD * is), site)];
      }

      real_t a = real_t(0.5) * dm[is] * dpinv[is-1];

      vt[0] += a * (xt[0] + xt[2]);
      vt[1] += a * (xt[1] + xt[3]);
      vt[2] += a * (xt[2] + xt[0]);
      vt[3] += a * (xt[3] + xt[1]);

      for(int id = 0; id < ND; ++id){
        int ivcd = ivc + NVC * id;
	vp[IDX2(Nin5, (ivcd + NVCD * is), site)] = vt[id];
      }

      for(int id = 0; id < ND; ++id){
        yt[id] += e[is] * vt[id];
      }

    }

    is = Ns-1;

    for(int id = 0; id < ND; ++id){
      xt[id] = vt[id];
    }

    for(int id = 0; id < ND; ++id){
        int ivcd = ivc + NVC * id;
      vt[id] = wp[IDX2(Nin5, (ivcd + NVCD * is), site)];
    }

    real_t a = real_t(0.5) * dm[is] * dpinv[is-1];

    vt[0] += a * (xt[0] + xt[2]);
    vt[1] += a * (xt[1] + xt[3]);
    vt[2] += a * (xt[2] + xt[0]);
    vt[3] += a * (xt[3] + xt[1]);

    vt[0] += -0.5 * (yt[0] - yt[2]);
    vt[1] += -0.5 * (yt[1] - yt[3]);
    vt[2] += -0.5 * (yt[2] - yt[0]);
    vt[3] += -0.5 * (yt[3] - yt[1]);

    for(int id = 0; id < ND; ++id){
        int ivcd = ivc + NVC * id;
      vp[IDX2(Nin5, (ivcd + NVCD * is), site)] = vt[id];
    }
    // L_inv completed

    is = Ns-1;

    a = dpinv[Ns-1];
    real_t f1 = 0.5 *  (1.0 + alpha);
    real_t f2 = 0.5 * (-1.0 + alpha);

    real_t vt1, vt2, vt3, vt4;
    vt1 = vp[IDX2(Nin5, (ID1 + ivc + NVCD * is), site)];
    vt2 = vp[IDX2(Nin5, (ID2 + ivc + NVCD * is), site)];
    vt3 = vp[IDX2(Nin5, (ID3 + ivc + NVCD * is), site)];
    vt4 = vp[IDX2(Nin5, (ID4 + ivc + NVCD * is), site)];
    vt[0] = a * (f1 * vt1 + f2 * vt3);
    vt[1] = a * (f1 * vt2 + f2 * vt4);
    vt[2] = a * (f1 * vt3 + f2 * vt1);
    vt[3] = a * (f1 * vt4 + f2 * vt2);

    for(int id = 0; id < ND; ++id){
        int ivcd = ivc + NVC * id;
      vp[IDX2(Nin5, (ivcd + NVCD * is), site)] = vt[id];
    }

    yt[0] = 0.5 * (vt[0] + vt[2]);
    yt[1] = 0.5 * (vt[1] + vt[3]);
    yt[2] = 0.5 * (vt[2] + vt[0]);
    yt[3] = 0.5 * (vt[3] + vt[1]);

    for (int is = Ns-2; is >= 0; --is) {

      for(int id = 0; id < ND; ++id){
        xt[id] = vt[id];
      }

      for(int id = 0; id < ND; ++id){
        int ivcd = ivc + NVC * id;
        vt[id] = vp[IDX2(Nin5, (ivcd + NVCD * is), site)];
      }

      real_t a = real_t(0.5) * dm[is];

      vt[0] += a * (xt[0] - xt[2]);
      vt[1] += a * (xt[1] - xt[3]);
      vt[2] += a * (xt[2] - xt[0]);
      vt[3] += a * (xt[3] - xt[1]);

      for(int id = 0; id < ND; ++id){
        vt[id] += - f[is] * yt[id];
      }

      real_t aa = dpinv[is];

      for(int id = 0; id < ND; ++id){
        vt[id] *= aa;
      }

      if(is == 0){
        real_t f1 = 0.5 * (1.0 + alpha);
        real_t f2 = 0.5 * (1.0 - alpha);
        vt1 = f1 * vt[0] + f2 * vt[2];
        vt2 = f1 * vt[1] + f2 * vt[3];
        vt3 = f1 * vt[2] + f2 * vt[0];
        vt4 = f1 * vt[3] + f2 * vt[1];
        vp[IDX2(Nin5, (ID1 + ivc + NVCD * is), site)] = vt1;
        vp[IDX2(Nin5, (ID2 + ivc + NVCD * is), site)] = vt2;
        vp[IDX2(Nin5, (ID3 + ivc + NVCD * is), site)] = vt3;
        vp[IDX2(Nin5, (ID4 + ivc + NVCD * is), site)] = vt4;
      }else{
        for(int id = 0; id < ND; ++id){
          int ivcd = ivc + NVC * id;
          vp[IDX2(Nin5, (ivcd + NVCD * is), site)] = vt[id];
        }
      }

    }

    }
  } // idx loop

 }
 }

}

//====================================================================
void mult_domainwall_5din_LUdaginv_dirac(
                        real_t *RESTRICT vp, real_t *RESTRICT wp,
                        int Ns, int *Nsize,
                        real_t *e, real_t *f,
                        real_t *dpinv, real_t *dm, real_t alpha)
{
  int Nx  = Nsize[0];
  int Ny  = Nsize[1];
  int Nz  = Nsize[2];
  int Nt  = Nsize[3];
  int Nst = Nx * Ny * Nz * Nt;
  int Nst_pad = CEIL_NWP(Nst);

  int Nin5 = NVCD * Ns;
  int size = Nin5 * Nst_pad;

#pragma acc data present(vp[0:size], wp[0:size]) \
                 copyin(Nst, Nst_pad, Nin5, Ns, e[0:Ns-1], f[0:Ns-1], \
                        dpinv[0:Ns], dm[0:Ns], alpha)
 {

#pragma acc parallel num_workers(NUM_WORKERS) vector_length(VECTOR_LENGTH)
 {

#pragma acc loop gang worker vector
  for (int idx = 0; idx < NVC * Nst_pad; ++idx) {
    int idx2_wp = idx / NWP;
    int idx_in  = idx % NWP;
    int ivc = idx2_wp % NVC;
    int idx_out = idx2_wp / NVC;
    int site = idx_in + NWP * idx_out;
    if(site < Nst){

    real_t vt[ND], yt[ND], xt[ND];

    int is = 0;

    real_t a = dpinv[0];
    real_t f1 = 0.5 * (1.0 + alpha);
    real_t f2 = 0.5 * (1.0 - alpha);
    {
    real_t vt1, vt2, vt3, vt4;
    vt1 = wp[IDX2(Nin5, (ID1 + ivc + NVCD * is), site)];
    vt2 = wp[IDX2(Nin5, (ID2 + ivc + NVCD * is), site)];
    vt3 = wp[IDX2(Nin5, (ID3 + ivc + NVCD * is), site)];
    vt4 = wp[IDX2(Nin5, (ID4 + ivc + NVCD * is), site)];
    vt[0] = a * (f1 * vt1 + f2 * vt3);
    vt[1] = a * (f1 * vt2 + f2 * vt4);
    vt[2] = a * (f1 * vt3 + f2 * vt1);
    vt[3] = a * (f1 * vt4 + f2 * vt2);
    }

    for(int id = 0; id < ND; ++id){
      int ivcd = ivc + NVC * id;
      vp[IDX2(Nin5, (ivcd + NVCD * is), site)] = vt[id];
    }

    for(int id = 0; id < ND; ++id){
      yt[id] = f[0] * vt[id];
    }

    for (int is = 1; is < Ns-1; ++is) {

      for(int id = 0; id < ND; ++id){
        xt[id] = vt[id];
      }

      for(int id = 0; id < ND; ++id){
        int ivcd = ivc + NVC * id;
        vt[id] = wp[IDX2(Nin5, (ivcd + NVCD * is), site)];
      }

      real_t a = real_t(0.5) * dm[is - 1];

      vt[0] += a * (xt[0] - xt[2]);
      vt[1] += a * (xt[1] - xt[3]);
      vt[2] += a * (xt[2] - xt[0]);
      vt[3] += a * (xt[3] - xt[1]);

      real_t aa = dpinv[is];

      for(int id = 0; id < ND; ++id){
        vt[id] *= aa;
      }

      for(int id = 0; id < ND; ++id){
        int ivcd = ivc + NVC * id;
        vp[IDX2(Nin5, (ivcd + NVCD * is), site)] = vt[id];
      }

      for(int id = 0; id < ND; ++id){
        yt[id] += f[is] * vt[id];
      }

    }

    is = Ns-1;

    for(int id = 0; id < ND; ++id){
      xt[id] = vt[id];
    }

    for(int id = 0; id < ND; ++id){
      int ivcd = ivc + NVC * id;
      vt[id] = wp[IDX2(Nin5, (ivcd + NVCD * is), site)];
    }

    a = real_t(0.5) * dm[is - 1];

    vt[0] += a * (xt[0] - xt[2]);
    vt[1] += a * (xt[1] - xt[3]);
    vt[2] += a * (xt[2] - xt[0]);
    vt[3] += a * (xt[3] - xt[1]);

    vt[0] += -0.5 * (yt[0] + yt[2]);
    vt[1] += -0.5 * (yt[1] + yt[3]);
    vt[2] += -0.5 * (yt[2] + yt[0]);
    vt[3] += -0.5 * (yt[3] + yt[1]);

    real_t aa = dpinv[is];

    for(int id = 0; id < ND; ++id){
      vt[id] *= aa;
    }

    real_t ff1 = 0.5 * ( 1.0 + alpha);
    real_t ff2 = 0.5 * (-1.0 + alpha);
    {
    real_t vt1, vt2, vt3, vt4;
    vt1 = ff1 * vt[0] + ff2 * vt[2];
    vt2 = ff1 * vt[1] + ff2 * vt[3];
    vt3 = ff1 * vt[2] + ff2 * vt[0];
    vt4 = ff1 * vt[3] + ff2 * vt[1];
    vp[IDX2(Nin5, (ID1 + ivc + NVCD*is), site)] = vt1;
    vp[IDX2(Nin5, (ID2 + ivc + NVCD*is), site)] = vt2;
    vp[IDX2(Nin5, (ID3 + ivc + NVCD*is), site)] = vt3;
    vp[IDX2(Nin5, (ID4 + ivc + NVCD*is), site)] = vt4;
    }
    // Udag_inv completed

    is = Ns-1;

    for(int id = 0; id < ND; ++id){
      int ivcd = ivc + NVC * id;
      vt[id] = vp[IDX2(Nin5, (ivcd + NVCD * is), site)];
    }

    for(int id = 0; id < ND; ++id){
      int ivcd = ivc + NVC * id;
      vp[IDX2(Nin5, (ivcd + NVCD * is), site)] = vt[id];
    }

    yt[0] = 0.5 * (vt[0] - vt[2]);
    yt[1] = 0.5 * (vt[1] - vt[3]);
    yt[2] = 0.5 * (vt[2] - vt[0]);
    yt[3] = 0.5 * (vt[3] - vt[1]);

    for (int is = Ns-2; is >= 0; --is) {

      for(int id = 0; id < ND; ++id){
        xt[id] = vt[id];
      }

      for(int id = 0; id < ND; ++id){
        int ivcd = ivc + NVC * id;
        vt[id] = vp[IDX2(Nin5, (ivcd + NVCD * is), site)];
      }

      real_t a = real_t(0.5) * dm[is + 1] * dpinv[is];

      vt[0] += a * (xt[0] + xt[2]);
      vt[1] += a * (xt[1] + xt[3]);
      vt[2] += a * (xt[2] + xt[0]);
      vt[3] += a * (xt[3] + xt[1]);

      for(int id = 0; id < ND; ++id){
        vt[id] += -e[is] * yt[id];
      }
      
      for(int id = 0; id < ND; ++id){
        int ivcd = ivc + NVC * id;
        vp[IDX2(Nin5, (ivcd + NVCD * is), site)] = vt[id];
      }
    }

    }
  } // idx loop end

 }
 }

}

#endif
//============================================================END=====
