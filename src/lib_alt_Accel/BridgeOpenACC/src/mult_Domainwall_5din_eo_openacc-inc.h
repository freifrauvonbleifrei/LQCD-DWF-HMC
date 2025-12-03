/*!
      @file    mult_Doainwall_5din_eo_openacc-inc.h
      @brief
      @author  Hideo Matsufuru (matufuru)
               $LastChangedBy: matufuru $
      @date    $LastChangedDate: 2013-01-22 13:51:53 #$
      @version $LastChangedRevision: 2653 $
*/

#ifndef MULT_DOMAINWALL_5DIN_EO_ACC_INCLUDED
#define MULT_DOMAINWALL_5DIN_EO_ACC_INCLUDED

//====================================================================
void mult_domainwall_5din_ee_5dir_dirac(
     real_t *RESTRICT vp, real_t *RESTRICT wp,
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

#pragma acc data present(vp[0:size], wp[0:size]) \
             copyin(Nst, Nst_pad, Ns, Nin5, mq, M0, b[0:Ns], c[0:Ns], alpha)
 {
#pragma acc parallel num_workers(NUM_WORKERS) vector_length(VECTOR_LENGTH)
 {

#pragma acc loop gang worker vector
  for (int site = 0; site < Nst_pad; ++site) {
   if(site < Nst){

    for (int is = 0; is < Ns; ++is) {

      real_t FF1 = b[is] * (4.0 - M0) + 1.0;
      real_t FF2 = c[is] * (4.0 - M0) - 1.0;

      real_t vt[NVCD], wt[NVCD];

      int is_up = (is+1) % Ns;
      real_t Fup   = 0.5 * FF2 * alpha;
      if (is == Ns-1) Fup = -0.5 * mq * FF2;

      for (int ivcd = 0; ivcd < NVCD; ++ivcd) {
        wt[ivcd]  = wp[IDX2(Nin5, (ivcd + NVCD * is_up), site)];
      }
      for (int ivc = 0; ivc < NVC; ++ivc) {
        vt[ID1 + ivc]  = Fup * (wt[ID1 + ivc] - wt[ID3 + ivc]);
        vt[ID2 + ivc]  = Fup * (wt[ID2 + ivc] - wt[ID4 + ivc]);
        vt[ID3 + ivc]  = Fup * (wt[ID3 + ivc] - wt[ID1 + ivc]);
        vt[ID4 + ivc]  = Fup * (wt[ID4 + ivc] - wt[ID2 + ivc]);
      }

      int is_dn = (is-1 + Ns) % Ns;
      real_t Fdn   = 0.5 * FF2 * alpha;
      if (is == 0) Fdn = -0.5 * mq * FF2;

      for (int ivcd = 0; ivcd < NVCD; ++ivcd) {
        wt[ivcd] = wp[IDX2(Nin5, (ivcd + NVCD * is_dn), site)];
      }
      for (int ivc = 0; ivc < NVC; ++ivc) {
        vt[ID1 + ivc] += Fdn * (wt[ID1 + ivc] + wt[ID3 + ivc]);
        vt[ID2 + ivc] += Fdn * (wt[ID2 + ivc] + wt[ID4 + ivc]);
        vt[ID3 + ivc] += Fdn * (wt[ID3 + ivc] + wt[ID1 + ivc]);
        vt[ID4 + ivc] += Fdn * (wt[ID4 + ivc] + wt[ID2 + ivc]);
      }

      if(is == 0){
	real_t fac1 = FF1 * 0.5 * ( 1.0 + alpha);
        real_t fac2 = FF1 * 0.5 * (-1.0 + alpha);
        for(int ivc = 0; ivc < NVC; ++ivc){
          real_t wt1, wt2, wt3, wt4;
          wt1 = wp[IDX2(Nin5, (ID1 + ivc + NVCD * is), site)];
          wt2 = wp[IDX2(Nin5, (ID2 + ivc + NVCD * is), site)];
          wt3 = wp[IDX2(Nin5, (ID3 + ivc + NVCD * is), site)];
          wt4 = wp[IDX2(Nin5, (ID4 + ivc + NVCD * is), site)];
          wt[ID1 + ivc] = fac1 * wt1 + fac2 * wt3;
          wt[ID2 + ivc] = fac1 * wt2 + fac2 * wt4;
          wt[ID3 + ivc] = fac1 * wt3 + fac2 * wt1;
          wt[ID4 + ivc] = fac1 * wt4 + fac2 * wt2;
	}
      }else if(is == Ns-1){
	real_t fac1 = FF1 * 0.5 * (1.0 + alpha);
        real_t fac2 = FF1 * 0.5 * (1.0 - alpha);
        for(int ivc = 0; ivc < NVC; ++ivc){
          real_t wt1, wt2, wt3, wt4;
          wt1 = wp[IDX2(Nin5, (ID1 + ivc + NVCD * is), site)];
          wt2 = wp[IDX2(Nin5, (ID2 + ivc + NVCD * is), site)];
          wt3 = wp[IDX2(Nin5, (ID3 + ivc + NVCD * is), site)];
          wt4 = wp[IDX2(Nin5, (ID4 + ivc + NVCD * is), site)];
          wt[ID1 + ivc] = fac1 * wt1 + fac2 * wt3;
          wt[ID2 + ivc] = fac1 * wt2 + fac2 * wt4;
          wt[ID3 + ivc] = fac1 * wt3 + fac2 * wt1;
          wt[ID4 + ivc] = fac1 * wt4 + fac2 * wt2;
	}
      }else{
        real_t f1 = FF1 * alpha;
	for (int ivcd = 0; ivcd < NVCD; ++ivcd) {
          wt[ivcd] = f1 * wp[IDX2(Nin5, (ivcd + NVCD * is), site)];
        }
      }

      for (int ivcd = 0; ivcd < NVCD; ++ivcd) {
        vp[IDX2(Nin5, (ivcd + NVCD * is), site)] = wt[ivcd] + vt[ivcd];
      }
    }
   }
  }  // site loop

 }
 }

}

//====================================================================
void mult_domainwall_5din_eo_5dir_dirac(
     real_t *RESTRICT yp, real_t *RESTRICT wp,
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

#pragma acc data present(yp[0:size], wp[0:size]) \
                 copyin(Nst, Nst_pad, Ns, Nin5, mq, M0, b[0:Ns], c[0:Ns], alpha)
 {
#pragma acc parallel num_workers(NUM_WORKERS) vector_length(VECTOR_LENGTH)
 {

#pragma acc loop gang worker vector
  for (int idx = 0; idx < Nst_pad * NVC; ++idx) {
    int idx2_wp = idx / NWP;
    int idx_in  = idx % NWP;
    int ivc = idx2_wp % NVC;
    int idx_out = idx2_wp / NVC;
    int site = idx_in + NWP*idx_out;
   if(site < Nst){

    for (int is = 0; is < Ns; ++is) {

      real_t vt1, vt2, vt3, vt4;
      real_t wt1, wt2, wt3, wt4;

      int is_up = (is+1) % Ns;
      real_t Fup   = 0.5 * c[is] * alpha;
      if (is == Ns-1) Fup = -0.5 * mq * c[is];
      wt1  = wp[IDX2(Nin5, (ID1 + ivc + NVCD * is_up), site)];
      wt2  = wp[IDX2(Nin5, (ID2 + ivc + NVCD * is_up), site)];
      wt3  = wp[IDX2(Nin5, (ID3 + ivc + NVCD * is_up), site)];
      wt4  = wp[IDX2(Nin5, (ID4 + ivc + NVCD * is_up), site)];

      vt1  = Fup * (wt1 - wt3);
      vt2  = Fup * (wt2 - wt4);
      vt3  = Fup * (wt3 - wt1);
      vt4  = Fup * (wt4 - wt2);

      int is_dn = (is-1 + Ns) % Ns;
      real_t Fdn   = 0.5 * c[is] * alpha;
      if (is == 0) Fdn = -0.5 * mq * c[is];
      wt1 = wp[IDX2(Nin5, (ID1 + ivc + NVCD * is_dn), site)];
      wt2 = wp[IDX2(Nin5, (ID2 + ivc + NVCD * is_dn), site)];
      wt3 = wp[IDX2(Nin5, (ID3 + ivc + NVCD * is_dn), site)];
      wt4 = wp[IDX2(Nin5, (ID4 + ivc + NVCD * is_dn), site)];

      vt1 += Fdn * (wt1+ wt3);
      vt2 += Fdn * (wt2+ wt4);
      vt3 += Fdn * (wt3+ wt1);
      vt4 += Fdn * (wt4+ wt2);

      if(is == 0){
        real_t b1 = b[is] * 0.5 * ( 1.0 + alpha);
        real_t b2 = b[is] * 0.5 * (-1.0 + alpha);
        real_t yt1, yt2, yt3, yt4;
        yt1 = wp[IDX2(Nin5, (ID1 + ivc + NVCD * is), site)];
        yt2 = wp[IDX2(Nin5, (ID2 + ivc + NVCD * is), site)];
        yt3 = wp[IDX2(Nin5, (ID3 + ivc + NVCD * is), site)];
        yt4 = wp[IDX2(Nin5, (ID4 + ivc + NVCD * is), site)];
        wt1 = b1 * yt1 + b2 * yt3;
        wt2 = b1 * yt2 + b2 * yt4;
        wt3 = b1 * yt3 + b2 * yt1;
        wt4 = b1 * yt4 + b2 * yt2;
      }else if(is == Ns-1){
        real_t b1 = b[is] * 0.5 * (1.0 + alpha);
        real_t b2 = b[is] * 0.5 * (1.0 - alpha);
        real_t yt1, yt2, yt3, yt4;
        yt1 = wp[IDX2(Nin5, (ID1 + ivc + NVCD * is), site)];
        yt2 = wp[IDX2(Nin5, (ID2 + ivc + NVCD * is), site)];
        yt3 = wp[IDX2(Nin5, (ID3 + ivc + NVCD * is), site)];
        yt4 = wp[IDX2(Nin5, (ID4 + ivc + NVCD * is), site)];
        wt1 = b1 * yt1 + b2 * yt3;
        wt2 = b1 * yt2 + b2 * yt4;
        wt3 = b1 * yt3 + b2 * yt1;
        wt4 = b1 * yt4 + b2 * yt2;
      }else{
        real_t bb = b[is] * alpha;
        wt1 = bb * wp[IDX2(Nin5, (ID1 + ivc + NVCD * is), site)];
        wt2 = bb * wp[IDX2(Nin5, (ID2 + ivc + NVCD * is), site)];
        wt3 = bb * wp[IDX2(Nin5, (ID3 + ivc + NVCD * is), site)];
        wt4 = bb * wp[IDX2(Nin5, (ID4 + ivc + NVCD * is), site)];
      }

      yp[IDX2(Nin5, (ID1 + ivc + NVCD * is), site)] = -0.5 * (wt1 + vt1);
      yp[IDX2(Nin5, (ID2 + ivc + NVCD * is), site)] = -0.5 * (wt2 + vt2);
      yp[IDX2(Nin5, (ID3 + ivc + NVCD * is), site)] = -0.5 * (wt3 + vt3);
      yp[IDX2(Nin5, (ID4 + ivc + NVCD * is), site)] = -0.5 * (wt4 + vt4);

    } //is

   } // site < Nst
  } // idx
 }

 }
}

//====================================================================
void mult_domainwall_5din_ee_5dirdag_dirac(
     real_t *RESTRICT vp, real_t *RESTRICT wp,
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

#pragma acc data present(vp[0:size], wp[0:size]) \
                 copyin(Nst, Nst_pad, Ns, Nin5, mq, M0, b[0:Ns], c[0:Ns], alpha)
 {
#pragma acc parallel num_workers(NUM_WORKERS) vector_length(VECTOR_LENGTH)
 {

#pragma acc loop gang worker vector
  for (int site = 0; site < Nst_pad; ++site) {
   if(site < Nst){

    for (int is = 0; is < Ns; ++is) {

      real_t vt[NVCD], xt[NVCD];

      real_t B1 = b[is] * (4.0 - M0) + 1.0;
      real_t a1 = -0.5 * b[is];

      if(is == 0){
	real_t fac1 = B1 * 0.5 * ( 1.0 + alpha);
        real_t fac2 = B1 * 0.5 * (-1.0 + alpha);
        for(int ivc = 0; ivc < NVC; ++ivc){
          real_t wt1, wt2, wt3, wt4;
          wt1 = wp[IDX2(Nin5, (ID1 + ivc + NVCD * is), site)];
          wt2 = wp[IDX2(Nin5, (ID2 + ivc + NVCD * is), site)];
          wt3 = wp[IDX2(Nin5, (ID3 + ivc + NVCD * is), site)];
          wt4 = wp[IDX2(Nin5, (ID4 + ivc + NVCD * is), site)];
          vt[ID1 + ivc] = fac1 * wt1 + fac2 * wt3;
          vt[ID2 + ivc] = fac1 * wt2 + fac2 * wt4;
          vt[ID3 + ivc] = fac1 * wt3 + fac2 * wt1;
          vt[ID4 + ivc] = fac1 * wt4 + fac2 * wt2;
	}
      }else if(is == Ns-1){
	real_t fac1 = B1 * 0.5 * (1.0 + alpha);
        real_t fac2 = B1 * 0.5 * (1.0 - alpha);
        for(int ivc = 0; ivc < NVC; ++ivc){
          real_t wt1, wt2, wt3, wt4;
          wt1 = wp[IDX2(Nin5, (ID1 + ivc + NVCD * is), site)];
          wt2 = wp[IDX2(Nin5, (ID2 + ivc + NVCD * is), site)];
          wt3 = wp[IDX2(Nin5, (ID3 + ivc + NVCD * is), site)];
          wt4 = wp[IDX2(Nin5, (ID4 + ivc + NVCD * is), site)];
          vt[ID1 + ivc] = fac1 * wt1 + fac2 * wt3;
          vt[ID2 + ivc] = fac1 * wt2 + fac2 * wt4;
          vt[ID3 + ivc] = fac1 * wt3 + fac2 * wt1;
          vt[ID4 + ivc] = fac1 * wt4 + fac2 * wt2;
	}
      }else{
        real_t bb = B1 * alpha;
        for (int ivcd = 0; ivcd < NVCD; ++ivcd) {
          vt[ivcd] = bb * wp[IDX2(Nin5, (ivcd + NVCD * is), site)];
        }
      }

      int is_up = (is+1) % Ns;
      for (int ivcd = 0; ivcd < NVCD; ++ivcd) {
        xt[ivcd] = wp[IDX2(Nin5, (ivcd + NVCD * is_up), site)];
      }

      real_t Fup = 0.5 * (c[is_up] * (4.0 - M0) - 1.0);
      if (is == Ns-1){
	Fup *= -mq;
      }else{
        Fup *= alpha;
      }
      for (int ivc = 0; ivc < NVC; ++ivc) {
        vt[ID1 + ivc] += Fup * (xt[ID1 + ivc] + xt[ID3 + ivc]);
        vt[ID2 + ivc] += Fup * (xt[ID2 + ivc] + xt[ID4 + ivc]);
        vt[ID3 + ivc] += Fup * (xt[ID3 + ivc] + xt[ID1 + ivc]);
        vt[ID4 + ivc] += Fup * (xt[ID4 + ivc] + xt[ID2 + ivc]);
      }

      int is_dn = (is-1 + Ns) % Ns;
      for (int ivcd = 0; ivcd < NVCD; ++ivcd) {
        xt[ivcd] = wp[IDX2(Nin5, (ivcd + NVCD * is_dn), site)];
      }

      real_t Fdn   = 0.5 * (c[is_dn] * (4.0 - M0) - 1.0);
      if (is == 0){
        Fdn *= -mq;
      }else{
        Fdn *= alpha;
      }
      for (int ivc = 0; ivc < NVC; ++ivc) {
        vt[ID1 + ivc] += Fdn * (xt[ID1 + ivc] - xt[ID3 + ivc]);
        vt[ID2 + ivc] += Fdn * (xt[ID2 + ivc] - xt[ID4 + ivc]);
        vt[ID3 + ivc] += Fdn * (xt[ID3 + ivc] - xt[ID1 + ivc]);
        vt[ID4 + ivc] += Fdn * (xt[ID4 + ivc] - xt[ID2 + ivc]);
      }

      for (int ivcd = 0; ivcd < NVCD; ++ivcd) {
        vp[IDX2(Nin5, (ivcd + NVCD * is), site)] = vt[ivcd];
      }
    }
   }
  }

 }
 }

}

//====================================================================
void mult_domainwall_5din_eo_5dirdag_dirac(
     real_t *RESTRICT vp, real_t *RESTRICT yp,
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

#pragma acc data present(vp[0:size], yp[0:size]) \
                 copyin(Nst, Nst_pad, Ns, Nin5, mq, M0, b[0:Ns], c[0:Ns], alpha)
 {
#pragma acc parallel num_workers(NUM_WORKERS) vector_length(VECTOR_LENGTH)
 {

#pragma acc loop gang worker vector
  for (int idx = 0; idx < NVC * Nst_pad; ++idx) {
    int idx2_wp = idx / NWP;
    int idx_in  = idx % NWP;
    int ivc = idx2_wp % NVC;
    int idx_out = idx2_wp / NVC;
    int site = idx_in + NWP*idx_out;
    if(site<Nst){

      for (int is = 0; is < Ns; ++is) {

      real_t vt1, vt2, vt3, vt4;
      real_t yt1, yt2, yt3, yt4;

      yt1  = yp[IDX2(Nin5, (ID1 + ivc + NVCD * is), site)];
      yt2  = yp[IDX2(Nin5, (ID2 + ivc + NVCD * is), site)];
      yt3  = yp[IDX2(Nin5, (ID3 + ivc + NVCD * is), site)];
      yt4  = yp[IDX2(Nin5, (ID4 + ivc + NVCD * is), site)];

      if(is == 0){
        real_t b1 = -0.5 * b[is] * 0.5 * ( 1.0 + alpha);
        real_t b2 = -0.5 * b[is] * 0.5 * (-1.0 + alpha);
        vt1  = b1 * yt3 + b2 * yt1;
        vt2  = b1 * yt4 + b2 * yt2;
        vt3  = b1 * yt1 + b2 * yt3;
        vt4  = b1 * yt2 + b2 * yt4;
      }else if(is == Ns-1){
        real_t b1 = -0.5 * b[is] * 0.5 * (1.0 + alpha);
        real_t b2 = -0.5 * b[is] * 0.5 * (1.0 - alpha);
        vt1  = b1 * yt3 + b2 * yt1;
        vt2  = b1 * yt4 + b2 * yt2;
        vt3  = b1 * yt1 + b2 * yt3;
        vt4  = b1 * yt2 + b2 * yt4;
      }else{
        real_t bb = -0.5 * b[is] * alpha;
        vt1  = bb * yt3;
        vt2  = bb * yt4;
        vt3  = bb * yt1;
        vt4  = bb * yt2;
      }

      int is_up = (is+1) % Ns;
      yt1  = yp[IDX2(Nin5, (ID1 + ivc + NVCD * is_up), site)];
      yt2  = yp[IDX2(Nin5, (ID2 + ivc + NVCD * is_up), site)];
      yt3  = yp[IDX2(Nin5, (ID3 + ivc + NVCD * is_up), site)];
      yt4  = yp[IDX2(Nin5, (ID4 + ivc + NVCD * is_up), site)];

      real_t Fup = -0.5 * c[is_up] * 0.5 * alpha;
      if (is == Ns-1) Fup = -0.5 * c[is_up] * (-0.5) * mq;

      vt1  += Fup * (yt3 + yt1);
      vt2  += Fup * (yt4 + yt2);
      vt3  += Fup * (yt1 + yt3);
      vt4  += Fup * (yt2 + yt4);

      int is_dn = (is-1 + Ns) % Ns;
      yt1  = yp[IDX2(Nin5, (ID1 + ivc + NVCD * is_dn), site)];
      yt2  = yp[IDX2(Nin5, (ID2 + ivc + NVCD * is_dn), site)];
      yt3  = yp[IDX2(Nin5, (ID3 + ivc + NVCD * is_dn), site)];
      yt4  = yp[IDX2(Nin5, (ID4 + ivc + NVCD * is_dn), site)];

      real_t Fdn = -0.5 * c[is_dn] * (0.5) * alpha;
      if (is == 0) Fdn = -0.5 * c[is_dn] * (-0.5) * mq;

      vt1  += Fdn * (yt3 - yt1);
      vt2  += Fdn * (yt4 - yt2);
      vt3  += Fdn * (yt1 - yt3);
      vt4  += Fdn * (yt2 - yt4);

      vp[IDX2(Nin5, (ID1 + ivc + NVCD * is), site)] = vt1;
      vp[IDX2(Nin5, (ID2 + ivc + NVCD * is), site)] = vt2;
      vp[IDX2(Nin5, (ID3 + ivc + NVCD * is), site)] = vt3;
      vp[IDX2(Nin5, (ID4 + ivc + NVCD * is), site)] = vt4;
      } // is

    } // site < Nst
  } // idx
 }

}
}

//====================================================================
void mult_domainwall_5din_eo_hopb_dirac_5d(
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

  int size = Nin5 * Nst_pad;
  int size_u = NDF * Nst_pad * 2 * NDIM;

#pragma acc data present(vp[0:size], up[0:size_u], wp[0:size])	\
                 copyin(Nst, Nst_pad, Nx, Ny, Nz, Nt, Nin5, Ns, \
                        bc[0:NDIM], do_comm[0:NDIM], ieo, jeo, jgm5)
 {

#pragma acc parallel num_workers(NUM_WORKERS) vector_length(VECTOR_LENGTH)
 {
  real_t *RESTRICT u_up = up;
  real_t *RESTRICT u_dn = up;

#pragma acc loop gang worker vector
  for (int idx = 0; idx < Ns * Nst_pad; ++idx) {
    int idx2_wp = idx / NWP;
    int idx_in  = idx % NWP;
    int is = idx2_wp % Ns;
    int idx_out = idx2_wp / Ns;

    int site = idx_in + NWP*idx_out;
    if(site < Nst) {

      int Nxy  = Nx  * Ny;
      int ix   = site % Nx;
      int iyzt = site / Nx;
      int iy   = iyzt % Ny;
      int izt  = site / Nxy;
      int iz   = izt % Nz;
      int it   = izt / Nz;
      int Nxyz = Nx * Ny * Nz;

      int keo  = (jeo + iy + iz + it) % 2;

      real_t u_0, u_1, u_2, u_3, u_4, u_5;
      real_t u_6, u_7, u_8, u_9, u10, u11;
      real_t u12, u13, u14, u15, u16, u17;
      real_t vt1_0, vt1_1, vt1_2, vt1_3, vt1_4, vt1_5;
      real_t vt2_0, vt2_1, vt2_2, vt2_3, vt2_4, vt2_5;
      real_t wt1r, wt1i, wt2r, wt2i;

      real_t v2_01, v2_11, v2_21, v2_31, v2_41, v2_51;
      real_t v2_02, v2_12, v2_22, v2_32, v2_42, v2_52;
      real_t v2_03, v2_13, v2_23, v2_33, v2_43, v2_53;
      real_t v2_04, v2_14, v2_24, v2_34, v2_44, v2_54;

      real_t *RESTRICT v1 = wp;
      real_t *RESTRICT v2 = vp;
 
#include "inc/mult_Domainwall_eo_xyz_openacc-inc.h"

#include "inc/mult_Domainwall_eo_t_dirac_openacc-inc.h"

      v2[IDX2_SP_5D_R(0,0,is,Ns,site)] = v2_01;
      v2[IDX2_SP_5D_I(0,0,is,Ns,site)] = v2_11;
      v2[IDX2_SP_5D_R(1,0,is,Ns,site)] = v2_21;
      v2[IDX2_SP_5D_I(1,0,is,Ns,site)] = v2_31;
      v2[IDX2_SP_5D_R(2,0,is,Ns,site)] = v2_41;
      v2[IDX2_SP_5D_I(2,0,is,Ns,site)] = v2_51;

      v2[IDX2_SP_5D_R(0,1,is,Ns,site)] = v2_02;
      v2[IDX2_SP_5D_I(0,1,is,Ns,site)] = v2_12;
      v2[IDX2_SP_5D_R(1,1,is,Ns,site)] = v2_22;
      v2[IDX2_SP_5D_I(1,1,is,Ns,site)] = v2_32;
      v2[IDX2_SP_5D_R(2,1,is,Ns,site)] = v2_42;
      v2[IDX2_SP_5D_I(2,1,is,Ns,site)] = v2_52;

      v2[IDX2_SP_5D_R(0,2,is,Ns,site)] = v2_03;
      v2[IDX2_SP_5D_I(0,2,is,Ns,site)] = v2_13;
      v2[IDX2_SP_5D_R(1,2,is,Ns,site)] = v2_23;
      v2[IDX2_SP_5D_I(1,2,is,Ns,site)] = v2_33;
      v2[IDX2_SP_5D_R(2,2,is,Ns,site)] = v2_43;
      v2[IDX2_SP_5D_I(2,2,is,Ns,site)] = v2_53;

      v2[IDX2_SP_5D_R(0,3,is,Ns,site)] = v2_04;
      v2[IDX2_SP_5D_I(0,3,is,Ns,site)] = v2_14;
      v2[IDX2_SP_5D_R(1,3,is,Ns,site)] = v2_24;
      v2[IDX2_SP_5D_I(1,3,is,Ns,site)] = v2_34;
      v2[IDX2_SP_5D_R(2,3,is,Ns,site)] = v2_44;
      v2[IDX2_SP_5D_I(2,3,is,Ns,site)] = v2_54;
     }

   } // site loop
 }  // acc parallel
 }  // acc data

}

//====================================================================
void mult_domainwall_5din_eo_hop1_dirac(
       real_t *RESTRICT buf_xp, real_t *RESTRICT buf_xm,
       real_t *RESTRICT buf_yp, real_t *RESTRICT buf_ym,
       real_t *RESTRICT buf_zp, real_t *RESTRICT buf_zm,
       real_t *RESTRICT buf_tp, real_t *RESTRICT buf_tm,
       real_t *RESTRICT up, real_t *RESTRICT wp,
       int Ns, int *bc, int *Nsize, int *do_comm, int ieo, int jeo,
       int jgm5)
{
  int Nx  = Nsize[0];
  int Ny  = Nsize[1];
  int Nz  = Nsize[2];
  int Nt  = Nsize[3];
  int Nst = Nx * Ny * Nz * Nt;
  int Nst_pad = CEIL_NWP(Nst);

  int Nin5   = NVCD * Ns;
  int Nin5bd = NVC * ND2 * Ns;

  int size   = Nin5 * Nst_pad;
  int size_u = NDF * Nst_pad * 2 * NDIM;

  int size_bx = Nin5bd * CEIL_NWP((Ny * Nz * Nt + 1)/2);
  int size_by = Nin5bd * CEIL_NWP(Nst/Ny);
  int size_bz = Nin5bd * CEIL_NWP(Nst/Nz);
  int size_bt = Nin5bd * CEIL_NWP(Nst/Nt);

#pragma acc data present(up[0:size_u], wp[0:size], \
                         buf_xp[0:size_bx], buf_xm[0:size_bx], \
                         buf_yp[0:size_by], buf_ym[0:size_by], \
                         buf_zp[0:size_bz], buf_zm[0:size_bz], \
                         buf_tp[0:size_bt], buf_tm[0:size_bt] ) \
                 copyin(Nst, Nst_pad, Nx, Ny, Nz, Nt, Nin5, Ns, \
                        bc[0:NDIM], do_comm[0:4], ieo, jeo, jgm5)
 {

  if (do_comm[0] > 0) {
    int idir = 0;
    int Nyzt = Ny * Nz * Nt;
    int Nyzt_pad = CEIL_NWP(Nyzt);

#pragma acc parallel async \
            num_workers(NUM_WORKERS) vector_length(VECTOR_LENGTH)
   {

#pragma acc loop gang worker vector
    for (int iyzt = 0; iyzt < Nyzt_pad; ++iyzt) {
     if(iyzt < Nyzt){

      int iy = iyzt % Ny;
      int iz = (iyzt/Ny) % Nz;
      int it = iyzt/(Ny * Nz);
      int keo = (jeo + iy + iz + it) % 2;

      if(keo == 1){
        int ix = 0;
        int iyzt2 = iyzt/2;
        int site = ix + Nx * iyzt;
        real_t bc2 = bc[0];
        for (int is = 0; is < Ns; ++is) {
          real_t wt[NVCD], vt[NVC * ND2];
	  if(jgm5 == 0){
            for(int ivcd = 0; ivcd < NVCD; ++ivcd){
              wt[ivcd] = wp[IDX2(Nin5, (ivcd + NVCD * is), site)];
            }
	  }else{
            for(int ivc = 0; ivc < NVC; ++ivc){
              wt[ivc+ID1] = wp[IDX2(Nin5, (ivc+ID3 + NVCD * is), site)];
              wt[ivc+ID2] = wp[IDX2(Nin5, (ivc+ID4 + NVCD * is), site)];
              wt[ivc+ID3] = wp[IDX2(Nin5, (ivc+ID1 + NVCD * is), site)];
              wt[ivc+ID4] = wp[IDX2(Nin5, (ivc+ID2 + NVCD * is), site)];
            }
	  }
          mult_wilson_xp1(vt, wt);
          for(int ivcd = 0; ivcd < NVC * ND2; ++ivcd){
            buf_xp[IDX2(Nin5bd, (ivcd + NVCD2 * is), iyzt2)] = bc2 * vt[ivcd];
          }
	}
      }

      if(keo == 0){
        int iyzt2 = iyzt/2;
        int ix = Nx-1;
        int site = ix + Nx * iyzt;
        real_t ut[NDF];
        load_u(ut, up, site + Nst_pad * (1-ieo + 2*idir));
        for (int is = 0; is < Ns; ++is) {
          real_t wt[NVCD], vt[NVC * ND2];
	  if(jgm5 == 0){
            for(int ivcd = 0; ivcd < NVCD; ++ivcd){
              wt[ivcd] = wp[IDX2(Nin5, (ivcd + NVCD * is), site)];
            }
	  }else{
            for(int ivc = 0; ivc < NVC; ++ivc){
              wt[ivc+ID1] = wp[IDX2(Nin5, (ivc+ID3 + NVCD * is), site)];
              wt[ivc+ID2] = wp[IDX2(Nin5, (ivc+ID4 + NVCD * is), site)];
              wt[ivc+ID3] = wp[IDX2(Nin5, (ivc+ID1 + NVCD * is), site)];
              wt[ivc+ID4] = wp[IDX2(Nin5, (ivc+ID2 + NVCD * is), site)];
            }
	  }
          mult_wilson_xm1(vt, ut, wt);
          for(int ivcd = 0; ivcd < NVC * ND2; ++ivcd){
            buf_xm[IDX2(Nin5bd, (ivcd + NVCD2 * is), iyzt2)] = vt[ivcd];
          }
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
        real_t wt[NVCD], vt[NVC * ND2];
        if(jgm5 == 0){
          for(int ivcd = 0; ivcd < NVCD; ++ivcd){
            wt[ivcd] = wp[IDX2(Nin5, (ivcd + NVCD * is), site)];
          }
        }else{
          for(int ivc = 0; ivc < NVC; ++ivc){
            wt[ivc+ID1] = wp[IDX2(Nin5, (ivc+ID3 + NVCD * is), site)];
            wt[ivc+ID2] = wp[IDX2(Nin5, (ivc+ID4 + NVCD * is), site)];
            wt[ivc+ID3] = wp[IDX2(Nin5, (ivc+ID1 + NVCD * is), site)];
            wt[ivc+ID4] = wp[IDX2(Nin5, (ivc+ID2 + NVCD * is), site)];
          }
        }
        mult_wilson_yp1(vt, wt);
        for(int ivcd = 0; ivcd < NVC * ND2; ++ivcd){
          buf_yp[IDX2(Nin5bd, (ivcd + NVCD2 * is), ixzt)] = bc2 * vt[ivcd];
	}
      }

      iy = Ny-1;
      site = ix + Nx * (iy + Ny * izt);
      real_t ut[NDF];
      load_u(ut, up, site + Nst_pad * (1-ieo + 2*idir));
      for (int is = 0; is < Ns; ++is) {
        real_t wt[NVCD], vt[NVC * ND2];
        if(jgm5 == 0){
          for(int ivcd = 0; ivcd < NVCD; ++ivcd){
            wt[ivcd] = wp[IDX2(Nin5, (ivcd + NVCD * is), site)];
          }
        }else{
          for(int ivc = 0; ivc < NVC; ++ivc){
            wt[ivc+ID1] = wp[IDX2(Nin5, (ivc+ID3 + NVCD * is), site)];
            wt[ivc+ID2] = wp[IDX2(Nin5, (ivc+ID4 + NVCD * is), site)];
            wt[ivc+ID3] = wp[IDX2(Nin5, (ivc+ID1 + NVCD * is), site)];
            wt[ivc+ID4] = wp[IDX2(Nin5, (ivc+ID2 + NVCD * is), site)];
          }
        }
        mult_wilson_ym1(vt, ut, wt);
        for(int ivcd = 0; ivcd < NVC * ND2; ++ivcd){
          buf_ym[IDX2(Nin5bd, (ivcd + NVCD2 * is), ixzt)] = vt[ivcd];
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
        real_t wt[NVCD], vt[NVC * ND2];
        if(jgm5 == 0){
          for(int ivcd = 0; ivcd < NVCD; ++ivcd){
            wt[ivcd] = wp[IDX2(Nin5, (ivcd + NVCD * is), site)];
          }
        }else{
          for(int ivc = 0; ivc < NVC; ++ivc){
            wt[ivc+ID1] = wp[IDX2(Nin5, (ivc+ID3 + NVCD * is), site)];
            wt[ivc+ID2] = wp[IDX2(Nin5, (ivc+ID4 + NVCD * is), site)];
            wt[ivc+ID3] = wp[IDX2(Nin5, (ivc+ID1 + NVCD * is), site)];
            wt[ivc+ID4] = wp[IDX2(Nin5, (ivc+ID2 + NVCD * is), site)];
          }
        }
        mult_wilson_zp1(vt, wt);
        for(int ivcd = 0; ivcd < NVC * ND2; ++ivcd){
          buf_zp[IDX2(Nin5bd, (ivcd + NVCD2 * is), ixyt)] = bc2 * vt[ivcd];
	}
      }

      iz = Nz-1;
      site = ixy + Nxy * (iz + Nz * it);
      real_t ut[NDF];
      load_u(ut, up, site + Nst_pad * (1-ieo + 2*idir));
      for (int is = 0; is < Ns; ++is) {
        real_t wt[NVCD], vt[NVC * ND2];
        if(jgm5 == 0){
          for(int ivcd = 0; ivcd < NVCD; ++ivcd){
            wt[ivcd] = wp[IDX2(Nin5, (ivcd + NVCD * is), site)];
          }
        }else{
          for(int ivc = 0; ivc < NVC; ++ivc){
            wt[ivc+ID1] = wp[IDX2(Nin5, (ivc+ID3 + NVCD * is), site)];
            wt[ivc+ID2] = wp[IDX2(Nin5, (ivc+ID4 + NVCD * is), site)];
            wt[ivc+ID3] = wp[IDX2(Nin5, (ivc+ID1 + NVCD * is), site)];
            wt[ivc+ID4] = wp[IDX2(Nin5, (ivc+ID2 + NVCD * is), site)];
          }
        }
        mult_wilson_zm1(vt, ut, wt);
        for(int ivcd = 0; ivcd < NVC * ND2; ++ivcd){
          buf_zm[IDX2(Nin5bd, (ivcd + NVCD2 * is), ixyt)] = vt[ivcd];
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
        real_t wt[NVCD], vt[NVC * ND2];
        if(jgm5 == 0){
          for(int ivcd = 0; ivcd < NVCD; ++ivcd){
            wt[ivcd] = wp[IDX2(Nin5, (ivcd + NVCD * is), site)];
          }
        }else{
          for(int ivc = 0; ivc < NVC; ++ivc){
            wt[ivc+ID1] = wp[IDX2(Nin5, (ivc+ID3 + NVCD * is), site)];
            wt[ivc+ID2] = wp[IDX2(Nin5, (ivc+ID4 + NVCD * is), site)];
            wt[ivc+ID3] = wp[IDX2(Nin5, (ivc+ID1 + NVCD * is), site)];
            wt[ivc+ID4] = wp[IDX2(Nin5, (ivc+ID2 + NVCD * is), site)];
          }
        }
        mult_wilson_tp1_dirac(vt, wt);
        for(int ivcd = 0; ivcd < NVC * ND2; ++ivcd){
          buf_tp[IDX2(Nin5bd, (ivcd + NVCD2 * is), ixyz)] = bc2 * vt[ivcd];
	}
      }

      it = Nt-1;
      site = ixyz + Nxyz * it;
      real_t ut[NDF];
      load_u(ut, up, site + Nst_pad * (1-ieo + 2*idir));
      for (int is = 0; is < Ns; ++is) {
        real_t wt[NVCD], vt[NVC * ND2];
        if(jgm5 == 0){
          for(int ivcd = 0; ivcd < NVCD; ++ivcd){
            wt[ivcd] = wp[IDX2(Nin5, (ivcd + NVCD * is), site)];
          }
        }else{
          for(int ivc = 0; ivc < NVC; ++ivc){
            wt[ivc+ID1] = wp[IDX2(Nin5, (ivc+ID3 + NVCD * is), site)];
            wt[ivc+ID2] = wp[IDX2(Nin5, (ivc+ID4 + NVCD * is), site)];
            wt[ivc+ID3] = wp[IDX2(Nin5, (ivc+ID1 + NVCD * is), site)];
            wt[ivc+ID4] = wp[IDX2(Nin5, (ivc+ID2 + NVCD * is), site)];
          }
        }
        mult_wilson_tm1_dirac(vt, ut, wt);
        for(int ivcd = 0; ivcd < NVC * ND2; ++ivcd){
          buf_tm[IDX2(Nin5bd, (ivcd + NVCD2 * is), ixyz)] = vt[ivcd];
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
void mult_domainwall_5din_eo_hop2_dirac(
      real_t *RESTRICT vp, real_t *RESTRICT up, real_t *RESTRICT wp,
      real_t *RESTRICT buf_xp, real_t *RESTRICT buf_xm,
      real_t *RESTRICT buf_yp, real_t *RESTRICT buf_ym,
      real_t *RESTRICT buf_zp, real_t *RESTRICT buf_zm,
      real_t *RESTRICT buf_tp, real_t *RESTRICT buf_tm,
      int Ns, int *bc, int *Nsize, int *do_comm, int ieo, int jeo)
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
  int size_u  = NDF  * Nst_pad * 2 * NDIM;
  int size_bx = Nin5bd * CEIL_NWP((Ny * Nz * Nt + 1)/2);
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
                 copyin(Nst, Nx, Ny, Nz, Nt, Nin5, Ns, \
                        bc[0:NDIM], do_comm[0:4], ieo, jeo)
 {

#pragma acc parallel num_workers(NUM_WORKERS) vector_length(VECTOR_LENGTH)
 {
  int Nxy  = Nx * Ny;
  int Nxyz = Nx * Ny * Nz;
  int Nst_pad = CEIL_NWP(Nst);

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
    int keo  = (jeo + iy + iz + it) % 2;

    int idir;

    for(int is = 0; is < Ns; ++is){

      real_t vL[NVCD];

      for(int ivcd = 0; ivcd < NVCD; ++ivcd){
        vL[ivcd] = 0.0;
      }

      int opr_any = 0;

      idir = 0;
      if (do_comm[idir] > 0) {

        if(ix == Nx-1 && keo == 1){
          real_t ut[NDF];
          load_u(ut, up, site + Nst_pad * (ieo + 2*idir));
          real_t wt[NVC * ND2];
          int iyzt2 = iyzt/2;
          for(int ivcd = 0; ivcd < NVCD2; ++ivcd){
            wt[ivcd] = buf_xp[IDX2(Nin5bd, ivcd + NVCD2 * is, iyzt2)];
          }
          mult_wilson_xp2(vL, ut, wt);
          ++opr_any;
        }

        if(ix == 0 && keo == 0){
          real_t bc2 = bc[0];
          int iyzt2  = iyzt/2;
          real_t wt[NVC * ND2];
          for(int ivcd = 0; ivcd < NVCD2; ++ivcd){
            wt[ivcd] = bc2 * buf_xm[IDX2(Nin5bd, ivcd + NVCD2 * is, iyzt2)];
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
          load_u(ut, up, site + Nst_pad * (ieo + 2*idir));
          real_t wt[NVC * ND2];
          for(int ivcd = 0; ivcd < NVCD2; ++ivcd){
            wt[ivcd] = buf_yp[IDX2(Nin5bd, ivcd + NVCD2 * is, ixzt)];
          }
          mult_wilson_yp2(vL, ut, wt);
          ++opr_any;
        }

        if (iy == 0) {
          real_t bc2 = bc[1];
          real_t wt[NVC * ND2];
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
          load_u(ut, up, site + Nst_pad * (ieo + 2*idir));
          real_t wt[NVC * ND2];
          for(int ivcd = 0; ivcd < NVCD2; ++ivcd){
            wt[ivcd] = buf_zp[IDX2(Nin5bd, ivcd + NVCD2 * is, ixyt)];
          }
          mult_wilson_zp2(vL, ut, wt);
          ++opr_any;
        }

        if (iz == 0) {
          real_t bc2 = bc[2];
          real_t wt[NVC * ND2];
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
          load_u(ut, up, site + Nst_pad * (ieo + 2*idir));
          real_t wt[NVC * ND2];
          for(int ivcd = 0; ivcd < NVCD2; ++ivcd){
            wt[ivcd] = buf_tp[IDX2(Nin5bd, ivcd + NVCD2 * is, ixyz)];
          }
          mult_wilson_tp2_dirac(vL, ut, wt);
          ++opr_any;
        }

        if (it == 0) {
          real_t bc2 = bc[3];
          real_t wt[NVC * ND2];
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
