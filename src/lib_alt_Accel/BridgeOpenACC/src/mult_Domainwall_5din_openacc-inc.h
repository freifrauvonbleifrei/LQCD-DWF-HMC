/*!
      @file    mult_Doainwall_5din_openacc-inc.h
      @brief
      @author  Hideo Matsufuru (matufuru)
               $LastChangedBy: matufuru $
      @date    $LastChangedDate: 2013-01-22 13:51:53 #$
      @version $LastChangedRevision: 2653 $
*/

#ifndef MULT_DOMAINWALL_5DIN_ACC_INCLUDED
#define MULT_DOMAINWALL_5DIN_ACC_INCLUDED

//====================================================================
void mult_domainwall_5din_5dir_dirac(
      real_t *RESTRICT vp, real_t *RESTRICT yp, real_t *RESTRICT wp,
      real_t mq, real_t M0, int Ns, real_t *b, real_t *c, real_t alpha,
      int *Nsize)
{
  int Nx  = Nsize[0];
  int Ny  = Nsize[1];
  int Nz  = Nsize[2];
  int Nt  = Nsize[3];
  int Nst = Nx * Ny * Nz * Nt;
  int Nst_pad = CEIL_NWP(Nst);

  int Nin5 = NVCD * Ns;
  int size = Nin5 * Nst;

#pragma acc data present(vp[0:size], yp[0:size], wp[0:size]) \
                 copyin(Nst, Nst_pad, Ns, Nin5, mq, M0, \
                        b[0:Ns], c[0:Ns], alpha)
 {
#pragma acc parallel num_workers(NUM_WORKERS) vector_length(VECTOR_LENGTH)
 {

#pragma acc loop gang worker vector
  for (int idx = 0; idx < NVC * Nst_pad; ++idx) {
    int idx2_wp = idx / NWP;
    int idx_in  = idx % NWP;
    int idx_out = idx2_wp / NVC;
    int site = idx_in + NWP * idx_out;
    int ivc  = idx2_wp % NVC;

    if(site < Nst){

    for (int is = 0; is < Ns; ++is) {

      real_t wt1, wt2, wt3, wt4;
      real_t vt1, vt2, vt3, vt4;
      real_t xt1, xt2, xt3, xt4;

      int is_up = (is+1) % Ns;
      real_t Fup   = 0.5 * alpha;
      if (is == Ns-1) Fup = -0.5 * mq;

      wt1 = wp[IDX2(Nin5, ID1 + ivc + NVCD * is_up, site)];
      wt2 = wp[IDX2(Nin5, ID2 + ivc + NVCD * is_up, site)];
      wt3 = wp[IDX2(Nin5, ID3 + ivc + NVCD * is_up, site)];
      wt4 = wp[IDX2(Nin5, ID4 + ivc + NVCD * is_up, site)];

      vt1 = Fup * (wt1 - wt3);
      vt2 = Fup * (wt2 - wt4);
      vt3 = Fup * (wt3 - wt1);
      vt4 = Fup * (wt4 - wt2);
      
      int is_dn = (is-1 + Ns) % Ns;
      real_t Fdn   = 0.5 * alpha;
      if (is == 0) Fdn = -0.5 * mq;

      wt1 = wp[IDX2(Nin5, ID1 + ivc + NVCD * is_dn, site)];
      wt2 = wp[IDX2(Nin5, ID2 + ivc + NVCD * is_dn, site)];
      wt3 = wp[IDX2(Nin5, ID3 + ivc + NVCD * is_dn, site)];
      wt4 = wp[IDX2(Nin5, ID4 + ivc + NVCD * is_dn, site)];

      vt1 += Fdn * (wt1 + wt3);
      vt2 += Fdn * (wt2 + wt4);
      vt3 += Fdn * (wt3 + wt1);
      vt4 += Fdn * (wt4 + wt2);
      
      wt1 = wp[IDX2(Nin5, ID1 + ivc + NVCD * is, site)];
      wt2 = wp[IDX2(Nin5, ID2 + ivc + NVCD * is, site)];
      wt3 = wp[IDX2(Nin5, ID3 + ivc + NVCD * is, site)];
      wt4 = wp[IDX2(Nin5, ID4 + ivc + NVCD * is, site)];

      real_t B1 = b[is] * (4.0 - M0) + 1.0;
      real_t C1 = c[is] * (4.0 - M0) - 1.0;
      real_t B2 = -0.5 * b[is];
      real_t C2 = -0.5 * c[is];

      if(alpha != 1.0){
        if(is == 0){
          real_t F1 = 0.5 * ( 1.0 + alpha);
          real_t F2 = 0.5 * (-1.0 + alpha);
          xt1 = F1 * wt1 + F2 * wt3;
          xt2 = F1 * wt2 + F2 * wt4;
          xt3 = F1 * wt3 + F2 * wt1;
          xt4 = F1 * wt4 + F2 * wt2;
          wt1 = xt1;
          wt2 = xt2;
          wt3 = xt3;
          wt4 = xt4;
        }else if(is == Ns-1){
          real_t F1 = 0.5 * (1.0 + alpha);
          real_t F2 = 0.5 * (1.0 - alpha);
          xt1 = F1 * wt1 + F2 * wt3;
          xt2 = F1 * wt2 + F2 * wt4;
          xt3 = F1 * wt3 + F2 * wt1;
          xt4 = F1 * wt4 + F2 * wt2;
          wt1 = xt1;
          wt2 = xt2;
          wt3 = xt3;
          wt4 = xt4;
        }else{
          B1 *= alpha;
          B2 *= alpha;
        }
      }

      int ivcs = ivc + NVCD * is;
      vp[IDX2(Nin5, ID1 + ivcs, site)] = B1 * wt1 + C1 * vt1;
      vp[IDX2(Nin5, ID2 + ivcs, site)] = B1 * wt2 + C1 * vt2;
      vp[IDX2(Nin5, ID3 + ivcs, site)] = B1 * wt3 + C1 * vt3;
      vp[IDX2(Nin5, ID4 + ivcs, site)] = B1 * wt4 + C1 * vt4;

      yp[IDX2(Nin5, ID1 + ivcs, site)] = B2 * wt1 + C2 * vt1;
      yp[IDX2(Nin5, ID2 + ivcs, site)] = B2 * wt2 + C2 * vt2;
      yp[IDX2(Nin5, ID3 + ivcs, site)] = B2 * wt3 + C2 * vt3;
      yp[IDX2(Nin5, ID4 + ivcs, site)] = B2 * wt4 + C2 * vt4;

    }
    }
  }

 }
 }

}

//====================================================================
void mult_domainwall_5din_5dirdag_dirac(
       real_t *RESTRICT vp, real_t *RESTRICT yp, real_t *RESTRICT wp,
       real_t mq, real_t M0, int Ns, real_t *b, real_t *c, real_t alpha,
       int *Nsize)
{
  int Nx  = Nsize[0];
  int Ny  = Nsize[1];
  int Nz  = Nsize[2];
  int Nt  = Nsize[3];
  int Nst = Nx * Ny * Nz * Nt;
  int Nst_pad = CEIL_NWP(Nst);

  int Nin5 = NVCD * Ns;
  int size = Nin5 * Nst_pad;

#pragma acc data present(vp[0:size], yp[0:size], wp[0:size]) \
                 copyin(Nst, Nst_pad, Ns, Nin5, mq, M0, \
                        b[0:Ns], c[0:Ns], alpha)
 {
#pragma acc parallel num_workers(NUM_WORKERS) vector_length(VECTOR_LENGTH)
 {

#pragma acc loop gang worker vector
  for (int idx = 0; idx < NVC * Nst_pad; ++idx) {
    int idx2_wp = idx / NWP;
    int idx_in  = idx % NWP;
    int idx_out = idx2_wp / NVC;
    int site = idx_in + NWP * idx_out;
    int ivc  = idx2_wp % NVC;

    if(site < Nst){

    for (int is = 0; is < Ns; ++is) {

      real_t yt1, yt2, yt3, yt4;
      real_t wt1, wt2, wt3, wt4;
      real_t vt1, vt2, vt3, vt4;
      real_t xt1, xt2, xt3, xt4;

      yt1 = yp[IDX2(Nin5, ID1 + ivc + NVCD * is, site)];
      yt2 = yp[IDX2(Nin5, ID2 + ivc + NVCD * is, site)];
      yt3 = yp[IDX2(Nin5, ID3 + ivc + NVCD * is, site)];
      yt4 = yp[IDX2(Nin5, ID4 + ivc + NVCD * is, site)];

      wt1 = wp[IDX2(Nin5, ID1 + ivc + NVCD * is, site)];
      wt2 = wp[IDX2(Nin5, ID2 + ivc + NVCD * is, site)];
      wt3 = wp[IDX2(Nin5, ID3 + ivc + NVCD * is, site)];
      wt4 = wp[IDX2(Nin5, ID4 + ivc + NVCD * is, site)];

      if(is == 0){
        real_t B1 = b[is] * (4.0 - M0) + 1.0;
        real_t a1 = -0.5 * b[is];
        xt1 = B1 * wt1 + a1 * yt3;
        xt2 = B1 * wt2 + a1 * yt4;
        xt3 = B1 * wt3 + a1 * yt1;
        xt4 = B1 * wt4 + a1 * yt2;

        real_t F1  = 0.5 * ( 1.0 + alpha);
        real_t F2  = 0.5 * (-1.0 + alpha);
        vt1 = F1 * xt1 + F2 * xt3;
        vt2 = F1 * xt2 + F2 * xt4;
        vt3 = F1 * xt3 + F2 * xt1;
        vt4 = F1 * xt4 + F2 * xt2;
      }else if(is == Ns-1){
        real_t B1 = b[is] * (4.0 - M0) + 1.0;
        real_t a1 = -0.5 * b[is];
        xt1 = B1 * wt1 + a1 * yt3;
        xt2 = B1 * wt2 + a1 * yt4;
        xt3 = B1 * wt3 + a1 * yt1;
        xt4 = B1 * wt4 + a1 * yt2;

        real_t F1  = 0.5 * (1.0 + alpha);
        real_t F2  = 0.5 * (1.0 - alpha);
        vt1 = F1 * xt1 + F2 * xt3;
        vt2 = F1 * xt2 + F2 * xt4;
        vt3 = F1 * xt3 + F2 * xt1;
        vt4 = F1 * xt4 + F2 * xt2;
      }else{
        real_t B1 = (b[is] * (4.0 - M0) + 1.0) * alpha;
        real_t a1 = -0.5 * b[is] * alpha;
        vt1 = B1 * wt1 + a1 * yt3;
        vt2 = B1 * wt2 + a1 * yt4;
        vt3 = B1 * wt3 + a1 * yt1;
        vt4 = B1 * wt4 + a1 * yt2;
      }

      int is_up = (is+1) % Ns;
      real_t C1 = c[is_up] * (4.0 - M0) - 1.0;
      real_t aup = -0.5 * c[is_up];

      yt1 = yp[IDX2(Nin5, ID1 + ivc + NVCD * is_up, site)];
      yt2 = yp[IDX2(Nin5, ID2 + ivc + NVCD * is_up, site)];
      yt3 = yp[IDX2(Nin5, ID3 + ivc + NVCD * is_up, site)];
      yt4 = yp[IDX2(Nin5, ID4 + ivc + NVCD * is_up, site)];

      wt1 = wp[IDX2(Nin5, ID1 + ivc + NVCD * is_up, site)];
      wt2 = wp[IDX2(Nin5, ID2 + ivc + NVCD * is_up, site)];
      wt3 = wp[IDX2(Nin5, ID3 + ivc + NVCD * is_up, site)];
      wt4 = wp[IDX2(Nin5, ID4 + ivc + NVCD * is_up, site)];

      xt1 = C1 * wt1 + aup * yt3;
      xt2 = C1 * wt2 + aup * yt4;
      xt3 = C1 * wt3 + aup * yt1;
      xt4 = C1 * wt4 + aup * yt2;

      real_t Fup = 0.5 * alpha;
      if (is == Ns-1) Fup = -0.5 * mq;

      vt1 += Fup * (xt1 + xt3);
      vt2 += Fup * (xt2 + xt4);
      vt3 += Fup * (xt3 + xt1);
      vt4 += Fup * (xt4 + xt2);

      int is_dn = (is-1 + Ns) % Ns;
      real_t C2 = c[is_dn] * (4.0 - M0) - 1.0;
      real_t adn = -0.5 * c[is_dn];

      yt1 = yp[IDX2(Nin5, ID1 + ivc + NVCD * is_dn, site)];
      yt2 = yp[IDX2(Nin5, ID2 + ivc + NVCD * is_dn, site)];
      yt3 = yp[IDX2(Nin5, ID3 + ivc + NVCD * is_dn, site)];
      yt4 = yp[IDX2(Nin5, ID4 + ivc + NVCD * is_dn, site)];

      wt1 = wp[IDX2(Nin5, ID1 + ivc + NVCD * is_dn, site)];
      wt2 = wp[IDX2(Nin5, ID2 + ivc + NVCD * is_dn, site)];
      wt3 = wp[IDX2(Nin5, ID3 + ivc + NVCD * is_dn, site)];
      wt4 = wp[IDX2(Nin5, ID4 + ivc + NVCD * is_dn, site)];

      xt1 = C2 * wt1 + adn * yt3;
      xt2 = C2 * wt2 + adn * yt4;
      xt3 = C2 * wt3 + adn * yt1;
      xt4 = C2 * wt4 + adn * yt2;

      real_t Fdn   = 0.5 * alpha;
      if (is == 0) Fdn = -0.5 * mq;

      vt1 += Fdn * (xt1 - xt3);
      vt2 += Fdn * (xt2 - xt4);
      vt3 += Fdn * (xt3 - xt1);
      vt4 += Fdn * (xt4 - xt2);

      vp[IDX2(Nin5, ID1 + ivc + NVCD * is, site)] = vt1;
      vp[IDX2(Nin5, ID2 + ivc + NVCD * is, site)] = vt2;
      vp[IDX2(Nin5, ID3 + ivc + NVCD * is, site)] = vt3;
      vp[IDX2(Nin5, ID4 + ivc + NVCD * is, site)] = vt4;

    }
    }
  }
 }
 }
}

//====================================================================
void mult_domainwall_5din_mult_gm5_dirac(
         real_t *RESTRICT vp, real_t *RESTRICT wp, int Ns, int *Nsize)
{
  int Nx  = Nsize[0];
  int Ny  = Nsize[1];
  int Nz  = Nsize[2];
  int Nt  = Nsize[3];
  int Nst = Nx * Ny * Nz * Nt;
  int Nst_pad = CEIL_NWP(Nst);

  int Nin5 = NVCD * Ns;
  int size = Nin5 * Nst_pad;

#pragma acc data present(vp[0:size], wp[0:size]) copyin(Nst, Nst_pad, Ns, Nin5)
 {
#pragma acc parallel num_workers(NUM_WORKERS) vector_length(VECTOR_LENGTH)
 {

#pragma acc loop gang worker vector
  for (int idx = 0; idx < Ns * Nst_pad; ++idx) {
    int idx2_wp = idx / NWP;
    int idx_in  = idx % NWP;
    int idx_out = idx2_wp / Ns;
    int site    = idx_in + NWP * idx_out;
    int is      = idx2_wp % Ns;
    if(site < Nst){

      for (int ivc = 0; ivc < NVC; ++ivc) {
        real_t vt1 = wp[IDX2(Nin5, ID1 + ivc + NVCD * is, site)];
        real_t vt2 = wp[IDX2(Nin5, ID2 + ivc + NVCD * is, site)];
        real_t vt3 = wp[IDX2(Nin5, ID3 + ivc + NVCD * is, site)];
        real_t vt4 = wp[IDX2(Nin5, ID4 + ivc + NVCD * is, site)];

	vp[IDX2(Nin5, ID1 + ivc + NVCD * is, site)] = vt3;
	vp[IDX2(Nin5, ID2 + ivc + NVCD * is, site)] = vt4;
	vp[IDX2(Nin5, ID3 + ivc + NVCD * is, site)] = vt1;
	vp[IDX2(Nin5, ID4 + ivc + NVCD * is, site)] = vt2;
      }

    }
  }

 }
 }

}

//====================================================================
void mult_domainwall_5din_mult_R(
         real_t *RESTRICT vp, real_t *RESTRICT wp, int Ns, int *Nsize)
{
  int Nx  = Nsize[0];
  int Ny  = Nsize[1];
  int Nz  = Nsize[2];
  int Nt  = Nsize[3];
  int Nst = Nx * Ny * Nz * Nt;
  int Nst_pad = CEIL_NWP(Nst);

  int Nin5 = NVCD * Ns;
  int size = Nin5 * Nst_pad;

#pragma acc data present(vp[0:size], wp[0:size]) copyin(Nst, Nst_pad, Ns, Nin5)
 {
#pragma acc parallel num_workers(NUM_WORKERS) vector_length(VECTOR_LENGTH)
 {

#pragma acc loop gang worker vector
  for (int site = 0; site < Nst_pad; ++site) {
   if(site < Nst){

    for (int is = 0; is < Ns; ++is) {

      real_t vt[NVCD];
      int isR = Ns-1 - is;

      for (int ivcd = 0; ivcd < NVCD; ++ivcd) {
        vt[ivcd] = wp[IDX2(Nin5, ivcd + NVCD * isR, site)];
      }

      for (int ivcd = 0; ivcd < NVCD; ++ivcd) {
        vp[IDX2(Nin5, ivcd + NVCD * is, site)] = vt[ivcd];
      }
    }
   }
  }

 }
 }

}

//====================================================================
void mult_domainwall_5din_mult_gm5R_dirac(
         real_t *RESTRICT vp, real_t *RESTRICT wp, int Ns, int *Nsize)
{
  int Nx  = Nsize[0];
  int Ny  = Nsize[1];
  int Nz  = Nsize[2];
  int Nt  = Nsize[3];
  int Nst = Nx * Ny * Nz * Nt;
  int Nst_pad = CEIL_NWP(Nst);

  int Nin5 = NVCD * Ns;
  int size = Nin5 * Nst_pad;

#pragma acc data present(vp[0:size], wp[0:size]) copyin(Nst, Nst_pad, Ns, Nin5)
 {
#pragma acc parallel num_workers(NUM_WORKERS) vector_length(VECTOR_LENGTH)
 {

#pragma acc loop gang worker vector
  for (int site = 0; site < Nst_pad; ++site) {
   if(site < Nst){

    for (int is = 0; is < Ns; ++is) {

      int isR = Ns-1 - is;
      real_t vt[NVCD], wt[NVCD];

      for (int ivcd = 0; ivcd < NVCD; ++ivcd) {
        wt[ivcd] = wp[IDX2(Nin5, ivcd + NVCD*isR, site)];
      }

      for (int ivc = 0; ivc < NVC; ++ivc) {
        vt[ID1 + ivc] = wt[ID3 + ivc];
        vt[ID2 + ivc] = wt[ID4 + ivc];
        vt[ID3 + ivc] = wt[ID1 + ivc];
        vt[ID4 + ivc] = wt[ID2 + ivc];
      }

      for (int ivcd = 0; ivcd < NVCD; ++ivcd) {
        vp[IDX2(Nin5, ivcd + NVCD*is, site)] = vt[ivcd];
      }

    }
   }
  }  // site loop
 }
 }

}

//====================================================================
void mult_domainwall_5din_hopb_dirac_5d(
     real_t *RESTRICT vp, real_t *RESTRICT up, real_t *RESTRICT wp,
     int Ns, int *bc, int *Nsize, int *do_comm, int flag)
{
  int Nx  = Nsize[0];
  int Ny  = Nsize[1];
  int Nz  = Nsize[2];
  int Nt  = Nsize[3];
  int Nst = Nx * Ny * Nz * Nt;
  int Nst_pad = CEIL_NWP(Nst);

  int Nin5 = NVCD * Ns;

  int size   = Nin5 * Nst_pad;
  int size_u = NDF * Nst_pad * NDIM;

#pragma acc data present(vp[0:size], up[0:size_u], wp[0:size]) \
                  copyin(Nst, Nst_pad, Nx, Ny, Nz, Nt, Nin5, Ns, \
                         bc[0:NDIM], do_comm[0:NDIM])
 {
#pragma acc parallel num_workers(NUM_WORKERS) vector_length(VECTOR_LENGTH)
 {

#pragma acc loop gang worker vector
  for (int idx = 0; idx < Ns * Nst_pad; ++idx) {
    int idx2_wp = idx / NWP;
    int idx_in  = idx % NWP;
    int is = idx2_wp % Ns;
    int idx_out = idx2_wp / Ns;

    int site = idx_in + NWP*idx_out;
    if(site < Nst) {

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

      int idir;

      real_t vL[NVCD];

      if(flag == 0){
        for(int ivcd = 0; ivcd < NVCD; ++ivcd){
          vL[ivcd] = 0.0;
        }
      }else{
        for(int ivcd = 0; ivcd < NVCD; ++ivcd){
          vL[ivcd] = vp[IDX2(Nin5, ivcd + NVCD * is, site)];
        }
      }

      idir = 0;

      if ((ix < Nx-1) || (do_comm[idir] == 0)) {
        int ix2 = (ix + 1) % Nx;
        int nei = ix2 + Nx * iyzt;
        real_t bc2 = 1.0;
        if(ix == Nx-1) bc2 = bc[0];

        real_t ut[NDF];
        load_u(ut, up, site + Nst_pad * idir);

        real_t wt[NVCD];
        for(int ivcd = 0; ivcd < NVCD; ++ivcd){
          wt[ivcd] = bc2 * wp[IDX2(Nin5, ivcd + NVCD * is, nei)];
	}
        mult_wilson_xpb(vL, ut, wt);
      }

      if ((ix > 0) || (do_comm[idir] == 0)) {
        int ix2 = (ix - 1 + Nx) % Nx;
        int nei = ix2 + Nx * iyzt;
        real_t bc2 = 1.0;
        if(ix == 0) bc2 = bc[0];

        real_t ut[NDF];
        load_u(ut, up, nei + Nst_pad * idir);

        real_t wt[NVCD];
        for(int ivcd = 0; ivcd < NVCD; ++ivcd){
          wt[ivcd] = bc2 * wp[IDX2(Nin5, ivcd + NVCD * is, nei)];
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
        load_u(ut, up, site + Nst_pad * idir);

        real_t wt[NVCD];
        for(int ivcd = 0; ivcd < NVCD; ++ivcd){
          wt[ivcd] = bc2 * wp[IDX2(Nin5, ivcd + NVCD * is, nei)];
	}
        mult_wilson_ypb(vL, ut, wt);
      }

      if ((iy > 0) || (do_comm[idir] == 0)) {
        int iy2 = (iy - 1 + Ny) % Ny;
        int nei = ix + Nx * (iy2 + Ny * izt);
        real_t bc2 = 1.0;
        if(iy == 0) bc2 = bc[1];

        real_t ut[NDF];
        load_u(ut, up, nei + Nst_pad * idir);

        real_t wt[NVCD];
        for(int ivcd = 0; ivcd < NVCD; ++ivcd){
          wt[ivcd] = bc2 * wp[IDX2(Nin5, ivcd + NVCD * is, nei)];
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
        load_u(ut, up, site + Nst_pad * idir);

        real_t wt[NVCD];
        for(int ivcd = 0; ivcd < NVCD; ++ivcd){
          wt[ivcd] = bc2 * wp[IDX2(Nin5, ivcd + NVCD * is, nei)];
	}
        mult_wilson_zpb(vL, ut, wt);
      }

      if ((iz > 0) || (do_comm[idir] == 0)) {
        int iz2 = (iz - 1 + Nz) % Nz;
        int nei = ixy + Nxy * (iz2 + Nz * it);
        real_t bc2 = 1.0;
        if(iz == 0) bc2 = bc[2];

        real_t ut[NDF];
        load_u(ut, up, nei + Nst_pad * idir);

        real_t wt[NVCD];
        for(int ivcd = 0; ivcd < NVCD; ++ivcd){
          wt[ivcd] = bc2 * wp[IDX2(Nin5, ivcd + NVCD * is, nei)];
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
        load_u(ut, up, site + Nst_pad * idir);

        real_t wt[NVCD];
        for(int ivcd = 0; ivcd < NVCD; ++ivcd){
          wt[ivcd] = bc2 * wp[IDX2(Nin5, ivcd + NVCD * is, nei)];
	}
        mult_wilson_tpb_dirac(vL, ut, wt);
      }

      if ((it > 0) || (do_comm[idir] == 0)) {
        int it2 = (it - 1 + Nt) % Nt;
        int nei = ixyz + Nxyz * it2;
        real_t bc2 = 1.0;
        if(it == 0) bc2 = bc[3];

        real_t ut[NDF];
        load_u(ut, up, nei + Nst_pad * idir);

        real_t wt[NVCD];
        for(int ivcd = 0; ivcd < NVCD; ++ivcd){
          wt[ivcd] = bc2 * wp[IDX2(Nin5, ivcd + NVCD * is, nei)];
	}
        mult_wilson_tmb_dirac(vL, ut, wt);
      }

      for(int ivcd = 0; ivcd < NVCD; ++ivcd){
        vp[IDX2(Nin5, ivcd + NVCD * is, site)] = vL[ivcd];
      }

   }
  } // idx loop

 }  // acc parallel
 }  // acc data

}


//====================================================================
void mult_domainwall_5din_hop1_dirac(
                  real_t *RESTRICT buf_xp, real_t *RESTRICT buf_xm,
                  real_t *RESTRICT buf_yp, real_t *RESTRICT buf_ym,
                  real_t *RESTRICT buf_zp, real_t *RESTRICT buf_zm,
                  real_t *RESTRICT buf_tp, real_t *RESTRICT buf_tm,
                  real_t *RESTRICT up, real_t *RESTRICT wp,
                  int Ns, int *bc, int *Nsize, int *do_comm)
{
  int Nx  = Nsize[0];
  int Ny  = Nsize[1];
  int Nz  = Nsize[2];
  int Nt  = Nsize[3];
  int Nst = Nx * Ny * Nz * Nt;
  int Nst_pad = CEIL_NWP(Nx * Ny * Nz * Nt);

  int Nin5   = NVCD * Ns;
  int Nin5bd = NVCD2 * Ns;

  int size   = Nin5 * Nst_pad;
  int size_u = NDF * Nst_pad * NDIM;

  int size_bx = Nin5bd * CEIL_NWP(Ny * Nz * Nt);
  int size_by = Nin5bd * CEIL_NWP(Nx * Nz * Nt);
  int size_bz = Nin5bd * CEIL_NWP(Nx * Ny * Nt);
  int size_bt = Nin5bd * CEIL_NWP(Nx * Ny * Nz);

#pragma acc data present(up[0:size_u], wp[0:size], \
                         buf_xp[0:size_bx], buf_xm[0:size_bx], \
                         buf_yp[0:size_by], buf_ym[0:size_by], \
                         buf_zp[0:size_bz], buf_zm[0:size_bz], \
                         buf_tp[0:size_bt], buf_tm[0:size_bt] ) \
                  copyin(Nst, Nst_pad, Nx, Ny, Nz, Nt, Nin5, Nin5bd, \
                         Ns, bc[0:NDIM], do_comm[0:4])
 {
  if (do_comm[0] > 0) {
    int idir = 0;
    int Nyzt = Ny * Nz * Nt;

#pragma acc parallel async \
            num_workers(NUM_WORKERS) vector_length(VECTOR_LENGTH)
   {
     //int Nyzt_pad = CEIL_NWP(Nyzt);

#pragma acc loop gang worker vector
     //for (int iyzt = 0; iyzt < Nyzt_pad; ++iyzt) {
    for (int iyzt = 0; iyzt < Nyzt; ++iyzt) {
     if(iyzt < Nyzt){
      int ix = 0;
      int site = ix + Nx * iyzt;
      real_t bc2 = bc[0];
      for (int is = 0; is < Ns; ++is) {
        real_t wt[NVCD], vt[NVCD2];
        for(int ivcd = 0; ivcd < NVCD; ++ivcd){
          wt[ivcd] = wp[IDX2(Nin5, ivcd + NVCD * is, site)];
	}
        mult_wilson_xp1(vt, wt);
        for(int ivcd = 0; ivcd < NVCD2; ++ivcd){
          buf_xp[IDX2(Nin5bd, ivcd + NVCD2 * is, iyzt)] = bc2 * vt[ivcd];
	}
      }

      ix = Nx-1;
      site = ix + Nx * iyzt;
      real_t ut[NDF];
      load_u(ut, up, site + Nst_pad * idir);
      for (int is = 0; is < Ns; ++is) {
        real_t wt[NVCD], vt[NVCD2];
        for(int ivcd = 0; ivcd < NVCD; ++ivcd){
          wt[ivcd] = wp[IDX2(Nin5, ivcd + NVCD * is, site)];
	}
        mult_wilson_xm1(vt, ut, wt);
        for(int ivcd = 0; ivcd < NVCD2; ++ivcd){
          buf_xm[IDX2(Nin5bd, ivcd + NVCD2 * is, iyzt)] = vt[ivcd];
	}
      }
     }
    }
   }
  } // do_comm[0]

  if (do_comm[1] > 0) {
    int idir = 1;
    int Nxzt = Nx * Nz * Nt;

#pragma acc parallel async \
            num_workers(NUM_WORKERS) vector_length(VECTOR_LENGTH)
   {
    int Nxzt_pad = CEIL_NWP(Nxzt);

#pragma acc loop gang worker vector
    for (int ixzt = 0; ixzt < Nxzt_pad; ++ixzt) {
     if(ixzt < Nxzt){
      int iy = 0;
      int ix  = ixzt % Nx;
      int izt = ixzt / Nx;
      int site = ix + Nx * (iy + Ny * izt);
      real_t bc2 = bc[1];
      for (int is = 0; is < Ns; ++is) {
        real_t wt[NVCD], vt[NVCD2];
        for(int ivcd = 0; ivcd < NVCD; ++ivcd){
          wt[ivcd] = wp[IDX2(Nin5, ivcd + NVCD * is, site)];
	}
        mult_wilson_yp1(vt, wt);
        for(int ivcd = 0; ivcd < NVCD2; ++ivcd){
          buf_yp[IDX2(Nin5bd, ivcd + NVCD2 * is, ixzt)] = bc2 * vt[ivcd];
	}
      }

      iy = Ny-1;
      site = ix + Nx * (iy + Ny * izt);
      real_t ut[NDF];
      load_u(ut, up, site + Nst_pad * idir);
      for (int is = 0; is < Ns; ++is) {
        real_t wt[NVCD], vt[NVCD2];
        for(int ivcd = 0; ivcd < NVCD; ++ivcd){
          wt[ivcd] = wp[IDX2(Nin5, ivcd + NVCD * is, site)];
	}
        mult_wilson_ym1(vt, ut, wt);
        for(int ivcd = 0; ivcd < NVCD2; ++ivcd){
          buf_ym[IDX2(Nin5bd, ivcd + NVCD2 * is, ixzt)] = vt[ivcd];
	}
      }
     }
    }
   }
  } // do_comm[1]

  if (do_comm[2] > 0) {
    int idir = 2;
    int Nxy  = Nx * Ny;
    int Nxyt = Nx * Ny * Nt;

#pragma acc parallel async \
            num_workers(NUM_WORKERS) vector_length(VECTOR_LENGTH)
   {
    int Nxyt_pad = CEIL_NWP(Nxyt);

#pragma acc loop gang worker vector
    for (int ixyt = 0; ixyt < Nxyt_pad; ++ixyt) {
     if(ixyt < Nxyt){
      int iz = 0;
      int ixy = ixyt % Nxy;
      int it  = ixyt / Nxy;
      int site = ixy + Nxy * (iz + Nz * it);
      real_t bc2 = bc[2];
      for (int is = 0; is < Ns; ++is) {
        real_t wt[NVCD], vt[NVCD2];
        for(int ivcd = 0; ivcd < NVCD; ++ivcd){
          wt[ivcd] = wp[IDX2(Nin5, ivcd + NVCD * is, site)];
	}
        mult_wilson_zp1(vt, wt);
        for(int ivcd = 0; ivcd < NVCD2; ++ivcd){
          buf_zp[IDX2(Nin5bd, ivcd + NVCD2 * is, ixyt)] = bc2 * vt[ivcd];
	}
      }

      iz = Nz-1;
      site = ixy + Nxy * (iz + Nz * it);
      real_t ut[NDF];
      load_u(ut, up, site + Nst_pad * idir);
      for (int is = 0; is < Ns; ++is) {
        real_t wt[NVCD], vt[NVCD2];
        for(int ivcd = 0; ivcd < NVCD; ++ivcd){
          wt[ivcd] = wp[IDX2(Nin5, ivcd + NVCD * is, site)];
	}
        mult_wilson_zm1(vt, ut, wt);
        for(int ivcd = 0; ivcd < NVCD2; ++ivcd){
          buf_zm[IDX2(Nin5bd, ivcd + NVCD2 * is, ixyt)] = vt[ivcd];
	}
      }
     }
    }
   }
  } // do_comm[2]

  if (do_comm[3] > 0) {
    int idir = 3;
    int Nxyz = Nx * Ny * Nz;

#pragma acc parallel async \
            num_workers(NUM_WORKERS) vector_length(VECTOR_LENGTH)
   {
    int Nxyz_pad = CEIL_NWP(Nxyz);

#pragma acc loop gang worker vector
    for (int ixyz = 0; ixyz < Nxyz_pad; ++ixyz) {
     if(ixyz < Nxyz){
      int it = 0;
      int site = ixyz + Nxyz * it;
      real_t bc2 = bc[3];
      for (int is = 0; is < Ns; ++is) {
        real_t wt[NVCD], vt[NVCD2];
        for(int ivcd = 0; ivcd < NVCD; ++ivcd){
          wt[ivcd] = wp[IDX2(Nin5, ivcd + NVCD * is, site)];
	}
        mult_wilson_tp1_dirac(vt, wt);
        for(int ivcd = 0; ivcd < NVCD2; ++ivcd){
          buf_tp[IDX2(Nin5bd, ivcd + NVCD2 * is, ixyz)] = bc2 * vt[ivcd];
	}
      }

      it = Nt-1;
      site = ixyz + Nxyz * it;
      real_t ut[NDF];
      load_u(ut, up, site + Nst_pad * idir);
      for (int is = 0; is < Ns; ++is) {
        real_t wt[NVCD], vt[NVCD2];
        for(int ivcd = 0; ivcd < NVCD; ++ivcd){
          wt[ivcd] = wp[IDX2(Nin5, ivcd + NVCD * is, site)];
	}
        mult_wilson_tm1_dirac(vt, ut, wt);
        for(int ivcd = 0; ivcd < NVCD2; ++ivcd){
          buf_tm[IDX2(Nin5bd, ivcd + NVCD2 * is, ixyz)] = vt[ivcd];
	}
      }
     }
    }
   }
  } // do_comm[3]

#pragma acc wait

 } // acc data

  if(do_comm[0] > 0){
    #pragma acc update async host (buf_xp[0:size_bx])
    #pragma acc update async host (buf_xm[0:size_bx])
  }
  if(do_comm[1] > 0){
    #pragma acc update async host (buf_yp[0:size_by])
    #pragma acc update async host (buf_ym[0:size_by])
  }
  if(do_comm[2] > 0){
    #pragma acc update async host (buf_zp[0:size_bz])
    #pragma acc update async host (buf_zm[0:size_bz])
  }
  if(do_comm[3] > 0){
    #pragma acc update async host (buf_tp[0:size_bt])
    #pragma acc update async host (buf_tm[0:size_bt])
  }

 #pragma acc wait

}

//====================================================================
void mult_domainwall_5din_hop2_dirac(
      real_t *RESTRICT vp, real_t *RESTRICT up, real_t *RESTRICT wp,
      real_t *RESTRICT buf_xp, real_t *RESTRICT buf_xm,
      real_t *RESTRICT buf_yp, real_t *RESTRICT buf_ym,
      real_t *RESTRICT buf_zp, real_t *RESTRICT buf_zm,
      real_t *RESTRICT buf_tp, real_t *RESTRICT buf_tm,
      int Ns, int *bc, int *Nsize, int *do_comm)
{
  int Nx  = Nsize[0];
  int Ny  = Nsize[1];
  int Nz  = Nsize[2];
  int Nt  = Nsize[3];
  int Nst = Nx * Ny * Nz * Nt;
  int Nst_pad = CEIL_NWP(Nst);

  int Nin5   = NVCD * Ns;
  int Nin5bd = (NVCD/2) * Ns;

  int size    = Nin5 * Nst_pad;
  int size_u  = NDF * Nst_pad * 4;
  int size_bx = Nin5bd * CEIL_NWP(Ny * Nz * Nt);
  int size_by = Nin5bd * CEIL_NWP(Nx * Nz * Nt);
  int size_bz = Nin5bd * CEIL_NWP(Nx * Ny * Nt);
  int size_bt = Nin5bd * CEIL_NWP(Nx * Ny * Nz);

  if(do_comm[0] > 0){
    #pragma acc update async device (buf_xp[0:size_bx])
    #pragma acc update async device (buf_xm[0:size_bx])
  }
  if(do_comm[1] > 0){
    #pragma acc update async device (buf_yp[0:size_by])
    #pragma acc update async device (buf_ym[0:size_by])
  }
  if(do_comm[2] > 0){
    #pragma acc update async device (buf_zp[0:size_bz])
    #pragma acc update async device (buf_zm[0:size_bz])
  }
  if(do_comm[3] > 0){
    #pragma acc update async device (buf_tp[0:size_bt])
    #pragma acc update async device (buf_tm[0:size_bt])
  }

#pragma acc wait

#pragma acc data present(buf_xp[0:size_bx], buf_xm[0:size_bx], \
                         buf_yp[0:size_by], buf_ym[0:size_by], \
                         buf_zp[0:size_bz], buf_zm[0:size_bz], \
                         buf_tp[0:size_bt], buf_tm[0:size_bt], \
                         vp[0:size], up[0:size_u]) \
                 copyin(Nst, Nst_pad, Nx, Ny, Nz, Nt, Nin5, Nin5bd, Ns, \
                        bc[0:NDIM], do_comm[0:4])
 {

#pragma acc parallel num_workers(NUM_WORKERS) vector_length(VECTOR_LENGTH)
 {
  int Nxy  = Nx * Ny;
  int Nxyz = Nx * Ny * Nz;

#pragma acc loop gang worker vector
  for (int site = 0; site < Nst_pad; ++site) {
   if(site < Nst){

    int ix   = site % Nx;
    int iyzt = site / Nx;
    int ixy  = site % Nxy;
    int iy   = iyzt % Ny;
    int izt  = site / Nxy;
    int iz   = izt % Nz;
    int it   = izt / Nz;
    int ixyz = site % Nxyz;

    int idir;

    for(int is = 0; is < Ns; ++is){

      real_t vL[NVCD];

      for(int ivcd = 0; ivcd < NVCD; ++ivcd){
        vL[ivcd] = 0.0;
      }

      int opr_any = 0;

      idir = 0;
      if (do_comm[idir] > 0) {

        if (ix == Nx-1) {
          real_t ut[NDF];
          load_u(ut, up, site + Nst_pad * idir);
          real_t wt[NVCD2];
          for(int ivcd = 0; ivcd < NVCD2; ++ivcd){
            wt[ivcd] = buf_xp[IDX2(Nin5bd, ivcd + NVCD2 * is, iyzt)];
          }
          mult_wilson_xp2(vL, ut, wt);
          ++opr_any;
        }

        if (ix == 0) {
          real_t bc2 = bc[0];
          real_t wt[NVCD2];
          for(int ivcd = 0; ivcd < NVCD2; ++ivcd){
            wt[ivcd] = bc2 * buf_xm[IDX2(Nin5bd, ivcd + NVCD2 * is, iyzt)];
          }
          mult_wilson_xm2(vL, wt);
          ++opr_any;
        }
      }

      idir = 1;
      if (do_comm[idir] > 0) {
        int ixzt = ix + Nx * izt;

        if (iy == Ny-1) {
          real_t ut[NDF];
          load_u(ut, up, site + Nst_pad * idir);
          real_t wt[NVCD2];
          for(int ivcd = 0; ivcd < NVCD2; ++ivcd){
            wt[ivcd] = buf_yp[IDX2(Nin5bd, ivcd + NVCD2 * is, ixzt)];
          }
          mult_wilson_yp2(vL, ut, wt);
          ++opr_any;
        }

        if (iy == 0) {
          real_t bc2 = bc[1];
          real_t wt[NVCD2];
          for(int ivcd = 0; ivcd < NVCD2; ++ivcd){
            wt[ivcd] = bc2 * buf_ym[IDX2(Nin5bd, ivcd + NVCD2 * is, ixzt)];
          }
          mult_wilson_ym2(vL, wt);
          ++opr_any;
        }
      }

      idir = 2;
      if (do_comm[idir] > 0) {
        int ixyt = ixy + Nxy * it;
        if (iz == Nz-1) {
          real_t ut[NDF];
          load_u(ut, up, site + Nst_pad * idir);
          real_t wt[NVCD2];
          for(int ivcd = 0; ivcd < NVCD2; ++ivcd){
            wt[ivcd] = buf_zp[IDX2(Nin5bd, ivcd + NVCD2 * is, ixyt)];
          }
          mult_wilson_zp2(vL, ut, wt);
          ++opr_any;
        }

        if (iz == 0) {
          real_t bc2 = bc[2];
          real_t wt[NVCD2];
          for(int ivcd = 0; ivcd < NVCD2; ++ivcd){
            wt[ivcd] = bc2 * buf_zm[IDX2(Nin5bd, ivcd + NVCD2 * is, ixyt)];
          }
          mult_wilson_zm2(vL, wt);
          ++opr_any;
        }
      }

      idir = 3;
      if (do_comm[idir] > 0) {
        if (it == Nt-1) {
          real_t ut[NDF];
          load_u(ut, up, site + Nst_pad * idir);
          real_t wt[NVCD2];
          for(int ivcd = 0; ivcd < NVCD2; ++ivcd){
            wt[ivcd] = buf_tp[IDX2(Nin5bd, ivcd + NVCD2 * is, ixyz)];
          }
          mult_wilson_tp2_dirac(vL, ut, wt);
          ++opr_any;
        }

        if (it == 0) {
          real_t bc2 = bc[3];
          real_t wt[NVCD2];
          for(int ivcd = 0; ivcd < NVCD2; ++ivcd){
            wt[ivcd] = bc2 * buf_tm[IDX2(Nin5bd, ivcd + NVCD2 * is, ixyz)];
          }
          mult_wilson_tm2_dirac(vL, wt);
          ++opr_any;
        }
      }

      if (opr_any > 0) {
        for(int ivcd = 0; ivcd < NVCD; ++ivcd){
          vp[IDX2(Nin5, ivcd + NVCD * is, site)] += vL[ivcd];
	}
      }

    }   // is loop
   }
  }   // site loop

 }  // acc parallel
 }  // acc data

}


#endif
//============================================================END=====
