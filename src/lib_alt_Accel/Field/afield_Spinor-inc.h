/*!
        @file    afield_Spinor-inc.h
        @brief
        @author  Issaku Kanamori (kanamori)
                 $LastChangedBy: kanamori $
        @date    $LastChangedDate: 2025-12-03 19:35:35 +0900 (2025年12月03日 (水)) $
        @version $LastChangedRevision: 2668 $
*/

#ifndef ACCEL_AFIELD_SPINOR_INC_INCLUDED
#define ACCEL_AFIELD_SPINOR_INC_INCLUDED

#include <cstdlib> 

#include "lib/ResourceManager/threadManager.h"

#include "lib_alt_Accel/Field/afield.h"
#include "lib_alt_Accel/inline/afield_th-inc.h"


namespace {
  inline int idxGr(int ic1, int ic2){ return     2*(ic1 + NC * ic2); }
  inline int idxGi(int ic1, int ic2){ return 1 + 2*(ic1 + NC * ic2); }

  inline int idxSr(int ic, int id){ return     2*(id + ND * ic); }
  inline int idxSi(int ic, int id){ return 1 + 2*(id + ND * ic); }
}

namespace Accel_Spinor{

//====================================================================
template <typename REALTYPE>
void mult_Gnv(AField<REALTYPE,ACCEL>& v,       const int exv,
              const AField<REALTYPE,ACCEL>& u, const int exu,
              const AField<REALTYPE,ACCEL>& w, const int exw)
{
#pragma omp barrier
  vout.crucial("mult_Gnv is not yet ready\n");
  exit(EXIT_FILURE);

#pragma omp barrier
}

//====================================================================
template <typename REALTYPE>
void mult_Gdv(AField<REALTYPE,ACCEL>& v,       const int exv,
              const AField<REALTYPE,ACCEL>& u, const int exu,
              const AField<REALTYPE,ACCEL>& w, const int exw)
{
#pragma omp barrier
  vout.crucial("mult_Gdv is not yet ready\n");
  exit(EXIT_FILURE);

#pragma omp barrier
}

//====================================================================
template <typename REALTYPE>
void tensorProd(AField<REALTYPE,ACCEL>& u, const int mu,
                const AField<REALTYPE,ACCEL>& v,
                const AField<REALTYPE,ACCEL>& w)
{
#pragma omp barrier
  vout.crucial("tensorProd is not yet ready\n");
  exit(EXIT_FILURE);
#pragma omp barrier
}

//====================================================================
template <typename REALTYPE>
void tensorProd_5din(AField<REALTYPE,ACCEL>& u, const int mu,
                     const AField<REALTYPE,ACCEL>& v,
                     const AField<REALTYPE,ACCEL>& w)
{
#pragma omp barrier

  // u[mu]_{ab} = sum_spinor v(a, spin)^* w(b, spin)
  typedef AField<REALTYPE,ACCEL> AFIELD;
  typedef REALTYPE real_t;

  int ith = ThreadManager::get_thread_id();
  int nth = ThreadManager::get_num_threads();
  int ith_kernel = 0;

  int Nst  = w.nvol();
  int Nst_pad = CEIL_NWP(Nst);
  int NinF = w.nin();
  int Ns   = NinF/NVCD;
  assert(NinF == Ns * NVCD);

  int size = NinF * Nst_pad;

  real_t* up = u.ptr(NDF * Nst_pad * mu);
  real_t* vp = const_cast<AFIELD*>(&v)->ptr(0);
  real_t* wp = const_cast<AFIELD*>(&w)->ptr(0);

  if (ith == ith_kernel){

#pragma acc data present(vp[0:size], wp[0:size], up[0:NDF*Nst_pad])        \
  copyin(Nst, Ns)
    {
#pragma acc parallel num_workers(NUM_WORKERS) vector_length(VECTOR_LENGTH)
      {
#pragma acc loop gang worker vector
        for (int idx = 0; idx < Nst; ++idx) {
          int site = idx;

          for(int ic1 = 0; ic1 < NC; ++ic1) {
            for(int ic2 = 0; ic2 < NC; ++ic2){
              real_t ur = 0.0;
              real_t ui = 0.0;
              for(int is = 0; is < Ns; ++is) {
                real_t v1r =  vp[IDX2_SP_5D_R(ic1, 0, is, Ns, site) ];
                real_t v1i =  vp[IDX2_SP_5D_I(ic1, 0, is, Ns, site) ];
                real_t v2r =  vp[IDX2_SP_5D_R(ic1, 1, is, Ns, site) ];
                real_t v2i =  vp[IDX2_SP_5D_I(ic1, 1, is, Ns, site) ];
                real_t v3r =  vp[IDX2_SP_5D_R(ic1, 2, is, Ns, site) ];
                real_t v3i =  vp[IDX2_SP_5D_I(ic1, 2, is, Ns, site) ];
                real_t v4r =  vp[IDX2_SP_5D_R(ic1, 3, is, Ns, site) ];
                real_t v4i =  vp[IDX2_SP_5D_I(ic1, 3, is, Ns, site) ];

                real_t w1r =  wp[IDX2_SP_5D_R(ic2, 0, is, Ns, site) ];
                real_t w1i =  wp[IDX2_SP_5D_I(ic2, 0, is, Ns, site) ];
                real_t w2r =  wp[IDX2_SP_5D_R(ic2, 1, is, Ns, site) ];
                real_t w2i =  wp[IDX2_SP_5D_I(ic2, 1, is, Ns, site) ];
                real_t w3r =  wp[IDX2_SP_5D_R(ic2, 2, is, Ns, site) ];
                real_t w3i =  wp[IDX2_SP_5D_I(ic2, 2, is, Ns, site) ];
                real_t w4r =  wp[IDX2_SP_5D_R(ic2, 3, is, Ns, site) ];
                real_t w4i =  wp[IDX2_SP_5D_I(ic2, 3, is, Ns, site) ];

                ur += (v1r*w1r + v1i*w1i);
                ui += (v1r*w1i - v1i*w1r);
                ur += (v2r*w2r + v2i*w2i);
                ui += (v2r*w2i - v2i*w2r);
                ur += (v3r*w3r + v3i*w3i);
                ui += (v3r*w3i - v3i*w3r);
                ur += (v4r*w4r + v4i*w4i);
                ui += (v4r*w4i - v4i*w4r);
              } // is
              up[IDX2_G_R(ic1,ic2,site)] = ur;
              up[IDX2_G_I(ic1,ic2,site)] = ui;
            } // ic2
          } // ic1
        } // idx
      } // acc parallel
    }
  } // ith == ith_kernel

#pragma omp barrier

 }



} // namespace ACCEL_Gauge

//============================================================END=====
#endif
