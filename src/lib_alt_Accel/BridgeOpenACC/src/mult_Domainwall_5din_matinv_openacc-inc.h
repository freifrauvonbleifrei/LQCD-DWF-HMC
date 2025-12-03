/*!
      @file    mult_Doainwall_5din_matinv_openacc-inc.h
      @brief
      @author  Hideo Matsufuru (matufuru)
               $LastChangedBy: matufuru $
      @date    $LastChangedDate: 2013-01-22 13:51:53 #$
      @version $LastChangedRevision: 2653 $
*/

#ifndef MULT_DOMAINWALL_5DIN_MATINV_ACC_INCLUDED
#define MULT_DOMAINWALL_5DIN_MATINV_ACC_INCLUDED

//====================================================================
void mult_domainwall_5din_ee_inv_dirac_5d(
                        real_t *RESTRICT vp, real_t *RESTRICT wp,
			int jd, int Ns,
                        real_t *RESTRICT mat_inv, int *Nsize)
{
  int Nx  = Nsize[0];
  int Ny  = Nsize[1];
  int Nz  = Nsize[2];
  int Nt  = Nsize[3];
  int Nst = Nx * Ny * Nz * Nt;
  int Nst_pad = CEIL_NWP(Nst);

  int Nin5 = NVCD * Ns;
  int size = Nin5 * Nst_pad;
  int mat_size  = ND2 * Ns;
  int mat_size2 = mat_size * mat_size;

#pragma acc data present(vp[0:size], wp[0:size]) \
                 copyin(Ns, Nin5, Nst, Nst_pad, mat_inv[0:mat_size2])
 {
#pragma acc parallel num_workers(NUM_WORKERS) vector_length(VECTOR_LENGTH)
 {

#pragma acc loop gang worker vector
  for (int idx = 0; idx < Ns * NC * Nst_pad; ++idx) {
    int idx2_wp = idx / NWP;
    int idx_in  = idx % NWP;
    int ic      = idx2_wp % NC;
    int is1     = (idx2_wp/NC) % Ns;
    int idx_out = idx2_wp/(NC * Ns);
    int site = idx_in + NWP * idx_out;
    if(site < Nst){

      for (int id = 0; id < ND2; ++id) {

        real_t vt1r = 0.0;
        real_t vt1i = 0.0;
        real_t vt2r = 0.0;
        real_t vt2i = 0.0;
        int idx1 = id + ND2 * is1;
        for (int is2 = 0; is2 < Ns; ++is2) {
          int idx2, ivcd1, ivcd2;
          real_t mat1, mat2;
          if(jd == 1){
            idx2 = 0 + ND2 * is2;
            mat1 = mat_inv[idx1 + (ND2 * Ns) * idx2];
            idx2 = 1 + ND2 * is2;
            mat2 = mat_inv[idx1 + (ND2 * Ns) * idx2];
          }else{
            idx2 = 0 + ND2 * is2;
            mat1 = mat_inv[idx2 + (ND2 * Ns) * idx1];
            idx2 = 1 + ND2 * is2;
            mat2 = mat_inv[idx2 + (ND2 * Ns) * idx1];
          }
          ivcd1 = 0 + 2 * (ic + NC * 0);
          ivcd2 = 0 + 2 * (ic + NC * 2);
          vt1r +=   mat1 * wp[IDX2(Nin5, (ivcd1 + NVCD * is2), site)]
                  + mat2 * wp[IDX2(Nin5, (ivcd2 + NVCD * is2), site)];

          ivcd1 = 1 + 2 * (ic + NC * 0);
          ivcd2 = 1 + 2 * (ic + NC * 2);
          vt1i +=   mat1 * wp[IDX2(Nin5, (ivcd1 + NVCD * is2), site)]
                  + mat2 * wp[IDX2(Nin5, (ivcd2 + NVCD * is2), site)];

          ivcd1 = 0 + 2 * (ic + NC * 1);
          ivcd2 = 0 + 2 * (ic + NC * 3);
          vt2r +=   mat1 * wp[IDX2(Nin5, (ivcd1 + NVCD * is2), site)]
                  + mat2 * wp[IDX2(Nin5, (ivcd2 + NVCD * is2), site)];

          ivcd1 = 1 + 2 * (ic + NC * 1);
          ivcd2 = 1 + 2 * (ic + NC * 3);
          vt2i +=   mat1 * wp[IDX2(Nin5, (ivcd1 + NVCD * is2), site)]
                  + mat2 * wp[IDX2(Nin5, (ivcd2 + NVCD * is2), site)];
        }
        int id1 = 2 * id;
        int id2 = 2 * id + 1;
        vp[IDX2(Nin5, 0 + 2*(ic+NC*id1) + NVCD * is1, site)] = vt1r;
        vp[IDX2(Nin5, 1 + 2*(ic+NC*id1) + NVCD * is1, site)] = vt1i;
        vp[IDX2(Nin5, 0 + 2*(ic+NC*id2) + NVCD * is1, site)] = vt2r;
        vp[IDX2(Nin5, 1 + 2*(ic+NC*id2) + NVCD * is1, site)] = vt2i;
      }
    }
  }

 }
 }

}

//====================================================================
void mult_domainwall_5din_ee_inv_dirac_4d(
                        real_t *RESTRICT vp, real_t *RESTRICT wp,
			int jd, int Ns,
                        real_t *RESTRICT mat_inv, int *Nsize)
{
  int Nx  = Nsize[0];
  int Ny  = Nsize[1];
  int Nz  = Nsize[2];
  int Nt  = Nsize[3];
  int Nst = Nx * Ny * Nz * Nt;
  int Nst_pad = CEIL_NWP(Nst);

  int Nin5 = NVCD * Ns;
  int size = Nin5 * Nst_pad;
  int mat_size  = ND2 * Ns;
  int mat_size2 = mat_size * mat_size;

#pragma acc data present(vp[0:size], wp[0:size]) \
                 copyin(Ns, Nin5, Nst, Nst_pad, mat_inv[0:mat_size2])
 {
#pragma acc parallel num_workers(NUM_WORKERS) vector_length(VECTOR_LENGTH)
 {

#pragma acc loop gang worker vector
  for (int idx = 0; idx < NC * Nst_pad; ++idx) {
    int idx2_wp = idx / NWP;
    int idx_in  = idx % NWP;
    int ic      = idx2_wp % NC;
    int idx_out = idx2_wp / NC;
    int site    = idx_in + NWP * idx_out;
   if(site < Nst){

    for (int is1 = 0; is1 < Ns; ++is1) {
      {
        for (int id = 0; id < ND2; ++id) {

          real_t vt1r = 0.0;
          real_t vt1i = 0.0;
          real_t vt2r = 0.0;
          real_t vt2i = 0.0;
          int idx1 = id + ND2 * is1;
          for (int is2 = 0; is2 < Ns; ++is2) {
            int idx2, ivcd1, ivcd2;
            real_t mat1, mat2;
	    if(jd == 1){
              idx2 = 0 + ND2 * is2;
              mat1 = mat_inv[idx1 + (ND2 * Ns) * idx2];
              idx2 = 1 + ND2 * is2;
              mat2 = mat_inv[idx1 + (ND2 * Ns) * idx2];
	    }else{
              idx2 = 0 + ND2 * is2;
              mat1 = mat_inv[idx2 + (ND2 * Ns) * idx1];
              idx2 = 1 + ND2 * is2;
              mat2 = mat_inv[idx2 + (ND2 * Ns) * idx1];
	    }
            ivcd1 = 0 + 2 * (ic + NC * 0);
            ivcd2 = 0 + 2 * (ic + NC * 2);
            vt1r +=   mat1 * wp[IDX2(Nin5, (ivcd1 + NVCD * is2), site)]
                    + mat2 * wp[IDX2(Nin5, (ivcd2 + NVCD * is2), site)];

            ivcd1 = 1 + 2 * (ic + NC * 0);
            ivcd2 = 1 + 2 * (ic + NC * 2);
            vt1i +=   mat1 * wp[IDX2(Nin5, (ivcd1 + NVCD * is2), site)]
                    + mat2 * wp[IDX2(Nin5, (ivcd2 + NVCD * is2), site)];

            ivcd1 = 0 + 2 * (ic + NC * 1);
            ivcd2 = 0 + 2 * (ic + NC * 3);
            vt2r +=   mat1 * wp[IDX2(Nin5, (ivcd1 + NVCD * is2), site)]
                    + mat2 * wp[IDX2(Nin5, (ivcd2 + NVCD * is2), site)];

            ivcd1 = 1 + 2 * (ic + NC * 1);
            ivcd2 = 1 + 2 * (ic + NC * 3);
            vt2i +=   mat1 * wp[IDX2(Nin5, (ivcd1 + NVCD * is2), site)]
                    + mat2 * wp[IDX2(Nin5, (ivcd2 + NVCD * is2), site)];
	  }
          int id1 = 2 * id;
          int id2 = 2 * id + 1;
          vp[IDX2(Nin5, 0 + 2*(ic+NC*id1) + NVCD * is1, site)] = vt1r;
          vp[IDX2(Nin5, 1 + 2*(ic+NC*id1) + NVCD * is1, site)] = vt1i;
          vp[IDX2(Nin5, 0 + 2*(ic+NC*id2) + NVCD * is1, site)] = vt2r;
          vp[IDX2(Nin5, 1 + 2*(ic+NC*id2) + NVCD * is1, site)] = vt2i;
	}
      }
    }

   }
  }

 }
 }

 }

#endif
//============================================================END=====
