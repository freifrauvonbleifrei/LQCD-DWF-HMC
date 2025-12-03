/*!
        @file    aforce_F_Wilson_Nf2.cpp
        @brief
        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: kanamori $
        @date    $LastChangedDate:: 2025-12-03 19:35:35 #$
        @version $LastChangedRevision: 2668 $
*/

#include "lib_alt/Force/Fermion/aforce_F_Wilson_Nf2.h"

#include<cassert>

// include files in core library
#include "lib/ResourceManager/threadManager.h"
#include "lib/IO/bridgeIO.h"
using Bridge::vout;

#include "lib_alt_QXS/inline/define_vlen.h"
#include "lib_alt_QXS/inline/define_params.h"

#define  VLEN     VLEND
#define  VLENX    VLENXD
#define  VLENY    VLENYD

typedef double real_t;

#include "lib_alt_QXS/inline/vsimd_double-inc.h"
#include "lib_alt_QXS/inline/vsimd_common_double-inc.h"
#include "lib_alt_QXS/inline/vsimd_Wilson_SU3_double-inc.h"

// include files in alt-code dorectories
#include "lib_alt_QXS/Field/afield.h"
#include "lib_alt_QXS/Field/aindex_lex.h"
#include "lib_alt_QXS/Field/afield-inc.h"
#include "lib_alt_QXS/Field/afield_Spinor-inc.h"
#include "lib_alt_QXS/Field/shiftAField_lex.h"
namespace Alt_Spinor = QXS_Spinor;
#include "lib_alt_QXS/Fopr/afopr_Wilson.h"

#include "lib_alt_QXS/BridgeQXS/bridgeQXS_Wilson.h"

#include "lib_alt/Force/Fermion/aforce_F_Wilson_Nf2-tmpl.h"

typedef AField<double,QXS> AField_d;
//template<> class ShiftAField_lex<AField_d>;

template<>
void AForce_F_Wilson_Nf2<AField_d>::init_impl()
{
  int Nc = CommonParameters::Nc();
  int Nd = CommonParameters::Nd();
  int NinF  = Nc * Nd * 2;
  int NinF2 = NinF/2;
  int NinG = Nc * Nc * 2;
  int Nvol = CommonParameters::Nvol();
  int Ndim = CommonParameters::Ndim();
  m_b1.reset(NinF2, Nvol, 1);
  m_b2.reset(NinF2, Nvol, 1);
}

template<>
void AForce_F_Wilson_Nf2<AField_d>::set_parameters_impl(const Parameters& )
{
  int ith = ThreadManager::get_thread_id();
  if (ith == 0) {
    int nin = 2*NC*ND2; // 2 spinor
    m_shift.reset(new ShiftAField_lex<AField_d>(nin, m_boundary));
  }
}


template<>
void AForce_F_Wilson_Nf2<AField_d>::force_udiv1_impl(AField_d& force,
                                                     const AField_d& zeta,
                                                     const AField_d& eta)
{
#pragma omp barrier

  const int Nst  = CommonParameters::Nvol();

  int Nstv = Nst / VLEN;
  int Nsize[4]={ CommonParameters::Nx()/VLENX,
                 CommonParameters::Ny()/VLENY,
                 CommonParameters::Nz(),
                 CommonParameters::Nt() };

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, Nstv);

  svbool_t pg = set_predicate();

  // m_eta3 = gm5 x zeta
  {
    real_t *in  = const_cast<real_t*>(zeta.ptr(0));
    real_t *out = m_eta3.ptr(0);
    BridgeQXS::mult_wilson_gm5_dirac(out, in, Nsize);
  }
  scal(m_eta3, -m_kappa);

  {// xp
    constexpr int mu=0;
    for(int site=is; site<ns; ++site){
      real_t* __restrict in  = const_cast<real_t*>(eta.ptr(0)) + site * VLEN * NVCD;
      real_t* __restrict out = m_b1.ptr(0) + site * VLEN * 2*NC*ND2;
      for(int ic=0; ic<NC; ++ic){
        svreal_t sp1r, sp1i, sp2r, sp2i;
        set_sp2_xp(pg, sp1r, sp1i, sp2r, sp2i, in, ic);
        save_vec(pg, out,          sp1r);
        save_vec(pg, out +   VLEN, sp1i);
        save_vec(pg, out + 2*VLEN, sp2r);
        save_vec(pg, out + 3*VLEN, sp2i);
        out += 4*VLEN;
      }
    }
    m_shift->backward(m_b2, m_b1, mu);
    for(int site=is; site<ns; ++site){
      real_t* __restrict in  = m_b2.ptr(0)  + site * VLEN * 2*NC*ND2;
      real_t* __restrict out = m_eta2.ptr(0)  + site * VLEN*NVCD;
      for(int ic=0; ic<NC; ++ic){
        int ic2 = ND * 2 * ic;
        svreal_t sp1r, sp1i, sp2r, sp2i;
        svreal_t msp1i, msp2i;
        load_vec(pg, sp1r, in  );
        load_vec(pg, sp1i, in +   VLEN);
        load_vec(pg, sp2r, in + 2*VLEN);
        load_vec(pg, sp2i, in + 3*VLEN);
        in  += 4*VLEN;

        flip_sign(pg, msp1i, sp1i);
        flip_sign(pg, msp2i, sp2i);
        save_vec(pg, out + VLEN*(ic2   + ID1),  sp1r);  // sp1
        save_vec(pg, out + VLEN*(ic2+1 + ID1),  sp1i);
        save_vec(pg, out + VLEN*(ic2   + ID2),  sp2r);  // sp2
        save_vec(pg, out + VLEN*(ic2+1 + ID2),  sp2i);
        save_vec(pg, out + VLEN*(ic2   + ID3), msp2i);  // -i x sp2
        save_vec(pg, out + VLEN*(ic2+1 + ID3),  sp2r);
        save_vec(pg, out + VLEN*(ic2   + ID4), msp1i);  // -i x sp1
        save_vec(pg, out + VLEN*(ic2+1 + ID4),  sp1r);
      }
    }
    Alt_Spinor::tensorProd(force, mu, m_eta3, m_eta2);
  } // xp

  {// yp
    constexpr int mu=1;
    for(int site=is; site<ns; ++site){
      real_t* __restrict in  = const_cast<real_t*>(eta.ptr(0)) + site * VLEN * NVCD;
      real_t* __restrict out = m_b1.ptr(0) + site * VLEN * 2*NC*ND2;
      for(int ic=0; ic<NC; ++ic){
        svreal_t sp1r, sp1i, sp2r, sp2i;
        set_sp2_yp(pg, sp1r, sp1i, sp2r, sp2i, in, ic);
        save_vec(pg, out,          sp1r);
        save_vec(pg, out +   VLEN, sp1i);
        save_vec(pg, out + 2*VLEN, sp2r);
        save_vec(pg, out + 3*VLEN, sp2i);
        out += 4*VLEN;
      }
    }
    m_shift->backward(m_b2, m_b1, mu);

    for(int site=is; site<ns; ++site){
      real_t* __restrict in  = m_b2.ptr(0)  + site * VLEN * 2*NC*ND2;
      real_t* __restrict out = m_eta2.ptr(0)  + site * VLEN*NVCD;
      for(int ic=0; ic<NC; ++ic){
        int ic2 = ND * 2 * ic;
        svreal_t sp1r, sp1i, sp2r, sp2i;
        svreal_t msp1r, msp1i;
        load_vec(pg, sp1r, in         );
        load_vec(pg, sp1i, in +   VLEN);
        load_vec(pg, sp2r, in + 2*VLEN);
        load_vec(pg, sp2i, in + 3*VLEN);
        in += 4*VLEN;

        flip_sign(pg, msp1r, sp1r);
        flip_sign(pg, msp1i, sp1i);
        save_vec(pg, out + VLEN*(ic2   + ID1),  sp1r);  // sp1
        save_vec(pg, out + VLEN*(ic2+1 + ID1),  sp1i);
        save_vec(pg, out + VLEN*(ic2   + ID2),  sp2r);  // sp2
        save_vec(pg, out + VLEN*(ic2+1 + ID2),  sp2i);
        save_vec(pg, out + VLEN*(ic2   + ID3),  sp2r);  // sp2
        save_vec(pg, out + VLEN*(ic2+1 + ID3),  sp2i);
        save_vec(pg, out + VLEN*(ic2   + ID4), msp1r);  // -sp1
        save_vec(pg, out + VLEN*(ic2+1 + ID4), msp1i);
      }
    }
    Alt_Spinor::tensorProd(force, mu, m_eta3, m_eta2);
  } // yp

  {// zp
    constexpr int mu=2;
    for(int site=is; site<ns; ++site){
      real_t* __restrict in  = const_cast<real_t*>(eta.ptr(0)) + site * VLEN * NVCD;
      real_t* __restrict out = m_b1.ptr(0) + site * VLEN * 2*NC*ND2;
      for(int ic=0; ic<NC; ++ic){
        svreal_t sp1r, sp1i, sp2r, sp2i;
        set_sp2_zp(pg, sp1r, sp1i, sp2r, sp2i, in, ic);
        save_vec(pg, out,          sp1r);
        save_vec(pg, out +   VLEN, sp1i);
        save_vec(pg, out + 2*VLEN, sp2r);
        save_vec(pg, out + 3*VLEN, sp2i);
        out += 4*VLEN;
      }
    }
    m_shift->backward(m_b2, m_b1, mu);

    for(int site=is; site<ns; ++site){
      real_t* __restrict in  = m_b2.ptr(0)  + site * VLEN * 2*NC*ND2;
      real_t* __restrict out = m_eta2.ptr(0)  + site * VLEN*NVCD;
      for(int ic=0; ic<NC; ++ic){
        int ic2 = ND * 2 * ic;
        svreal_t sp1r, sp1i, sp2r, sp2i;
        svreal_t msp1i, msp2r;
        load_vec(pg, sp1r, in         );
        load_vec(pg, sp1i, in +   VLEN);
        load_vec(pg, sp2r, in + 2*VLEN);
        load_vec(pg, sp2i, in + 3*VLEN);
        in += 4*VLEN;

        flip_sign(pg, msp1i, sp1i);
        flip_sign(pg, msp2r, sp2r);
        save_vec(pg, out + VLEN*(ic2   + ID1),  sp1r);  // sp1
        save_vec(pg, out + VLEN*(ic2+1 + ID1),  sp1i);
        save_vec(pg, out + VLEN*(ic2   + ID2),  sp2r);  // sp2
        save_vec(pg, out + VLEN*(ic2+1 + ID2),  sp2i);
        save_vec(pg, out + VLEN*(ic2   + ID3), msp1i);  // i x sp1
        save_vec(pg, out + VLEN*(ic2+1 + ID3),  sp1r);
        save_vec(pg, out + VLEN*(ic2   + ID4),  sp2i);  // -i x sp2
        save_vec(pg, out + VLEN*(ic2+1 + ID4), msp2r);
      }
    }
    Alt_Spinor::tensorProd(force, mu, m_eta3, m_eta2);
  } // zp

  {// tp (Dirac rep.)
    constexpr int mu=3;
    for(int site=is; site<ns; ++site){
      real_t* __restrict in  = const_cast<real_t*>(eta.ptr(0)) + site * VLEN * NVCD;
      real_t* __restrict out = m_b1.ptr(0) + site * VLEN * 2*NC*ND2;
      for(int ic=0; ic<NC; ++ic){
        svreal_t sp1r, sp1i, sp2r, sp2i;
        set_sp2_tp_dirac(pg, sp1r, sp1i, sp2r, sp2i, in, ic);
        save_vec(pg, out,          sp1r);
        save_vec(pg, out +   VLEN, sp1i);
        save_vec(pg, out + 2*VLEN, sp2r);
        save_vec(pg, out + 3*VLEN, sp2i);
        out += 4*VLEN;
      }
    }
    m_shift->backward(m_b2, m_b1, mu);

    for(int site=is; site<ns; ++site){
      real_t* __restrict in  = m_b2.ptr(0)  + site * VLEN * 2*NC*ND2;
      real_t* __restrict out = m_eta2.ptr(0)  + site * VLEN*NVCD;
      for(int ic=0; ic<NC; ++ic){
        int ic2 = ND * 2 * ic;
        svreal_t sp1r, sp1i, sp2r, sp2i;
        svreal_t vzero;

        load_vec(pg, sp1r, in         );
        load_vec(pg, sp1i, in +   VLEN);
        load_vec(pg, sp2r, in + 2*VLEN);
        load_vec(pg, sp2i, in + 3*VLEN);
        in += 4*VLEN;

        clear_vec(pg, vzero);

        save_vec(pg, out + VLEN*(ic2   + ID1), vzero);  // 0
        save_vec(pg, out + VLEN*(ic2+1 + ID1), vzero);
        save_vec(pg, out + VLEN*(ic2   + ID2), vzero);  // 0
        save_vec(pg, out + VLEN*(ic2+1 + ID2), vzero);
        save_vec(pg, out + VLEN*(ic2   + ID3),  sp1r);  // sp1
        save_vec(pg, out + VLEN*(ic2+1 + ID3),  sp1i);
        save_vec(pg, out + VLEN*(ic2   + ID4),  sp2r);  // sp2
        save_vec(pg, out + VLEN*(ic2+1 + ID4),  sp2i);
      }
    }
    Alt_Spinor::tensorProd(force, mu, m_eta3, m_eta2);
  } // tp
  //  scal(force, -m_kappa); --> multiplied to m_eta3 instead
#pragma omp barrier

}

// explicit instanciation for AField<double,QXS>.
template<>
const std::string AForce_F_Wilson_Nf2<AField_d >::class_name
                         = "AForce_F_Wilson_Nf2<AField<double,QXS> >";


template class AForce_F_Wilson_Nf2<AField_d>;

//============================================================END=====
