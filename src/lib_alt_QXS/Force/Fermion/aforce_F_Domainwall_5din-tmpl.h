/*!
        @file    aforce_F_Domainwall_5din-tmpl.h

        @brief

        @author  Issaku Kanamori (kanamori)
                 $LastChangedBy: kanamori $

        @date    $LastChangedDate:: 2025-12-03 19:35:35 #$

        @version $LastChangedRevision: 2668 $
*/

#include "lib_alt_QXS/Force/Fermion/aforce_F_Domainwall_5din.h"
#include "lib/Force/Fermion/force_F_Domainwall.h"
template<typename AFIELD>
const std::string AForce_F_Domainwall_5din<AFIELD>::class_name
                               = "AForce_F_Domainwall_5din<AFIELD>";
//====================================================================
template<typename AFIELD>
void AForce_F_Domainwall_5din<AFIELD>::init(const Parameters& params)
{
  ThreadManager::assert_single_thread(class_name);

  std::string vlevel;
  if (!params.fetch_string("verbose_level", vlevel)) {
    m_vl = vout.set_verbose_level(vlevel);
  } else {
    m_vl = CommonParameters::Vlevel();
  }

  vout.general(m_vl, "%s: construction\n", class_name.c_str());
  vout.increase_indent();

  int Nc   = CommonParameters::Nc();
  int Nd   = CommonParameters::Nd();
  int Nvol = CommonParameters::Nvol();
  int Ndim = CommonParameters::Ndim();
  int NinG = Nc * Nc * 2;

  m_Nvcd = 2 * Nc * Nd;
  m_boundary.resize(Ndim);

  set_parameters(params);

  int err = 0;

  Parameters params_kernel = params;
  //double kappa = 1.0 / (8.0 - 2.0 * m_M0);
  //params_kernel.set_double("hopping_parameter", kappa);


  m_v1.reset(m_NinF, Nvol, 1);
  m_v2.reset(m_NinF, Nvol, 1);
  m_h1.reset(m_NinF/2, Nvol, 1); // 2-spin
  m_h2.reset(m_NinF/2, Nvol, 1);
  m_z1.reset(m_NinF, Nvol, 1);
  m_z2.reset(m_NinF, Nvol, 1);

  m_force1.reset(NinG, Nvol, Ndim);
  m_force2.reset(NinG, Nvol, Ndim);

  vout.decrease_indent();

  vout.general(m_vl, "%s: construction finished.\n",
               class_name.c_str());

}

//====================================================================
template<typename AFIELD>
void AForce_F_Domainwall_5din<AFIELD>::tidyup()
{
  //  delete m_force_w;
  //  delete m_fopr_dw;
  //  delete m_fopr_w;
}

//====================================================================
template<typename AFIELD>
void AForce_F_Domainwall_5din<AFIELD>::set_parameters(const Parameters& params)
{
#pragma omp barrier

  int ith = ThreadManager::get_thread_id();

  std::string vlevel;
  if (!params.fetch_string("verbose_level", vlevel)) {
    if(ith == 0) m_vl = vout.set_verbose_level(vlevel);
  }

  //- fetch and check input parameters
  double mq, M0, b, c;
  int Ns;
  std::vector<int> bc;
  std::string repr;

  int err = 0;
  err += params.fetch_double("quark_mass", mq);
  err += params.fetch_double("domain_wall_height", M0);
  err += params.fetch_int("extent_of_5th_dimension", Ns);
  err += params.fetch_int_vector("boundary_condition", bc);
  err += params.fetch_string("gamma_matrix_type", repr);

  if (err) {
    vout.crucial(m_vl, "Error at %s: input parameter not found.\n",
                 class_name.c_str());
    exit(EXIT_FAILURE);
  }

  if(ith == 0) m_repr = repr;

  int err2 = 0;
  err2 += params.fetch_double("coefficient_b", b);
  err2 += params.fetch_double("coefficient_c", c);

  if (err2) {
    vout.general(m_vl, "  coefficients b, c are not provided:"
                       " set to Shamir's form.\n");
    b = 1.0;
    c = 0.0;
  }

  set_parameters(real_t(mq), real_t(M0), Ns, bc, real_t(b), real_t(c));

  if( !m_fopr_dw) {
    m_fopr_dw.reset(new AFopr_Domainwall_5din<AFIELD>(params));
  } else {
    m_fopr_dw->set_parameters(params);
  }

#pragma omp barrier

}


//====================================================================
template<typename AFIELD>
void AForce_F_Domainwall_5din<AFIELD>::set_parameters(
                                          const real_t mq,
                                          const real_t M0,
                                          const int Ns,
                                          const std::vector<int> bc,
                                          const real_t b,
                                          const real_t c)
{
  const int Ndim = CommonParameters::Ndim();
  const int Nvol = CommonParameters::Nvol();

  //- range check
  int err = 0;
  err += ParameterCheck::non_zero(mq);
  err += ParameterCheck::non_zero(M0);
  err += ParameterCheck::non_zero(Ns);

  if (err) {
    vout.crucial(m_vl, "Error at %s: parameter range check failed.\n",
                 class_name.c_str());
    exit(EXIT_FAILURE);
  }

  assert(bc.size() == Ndim);

  int ith = ThreadManager::get_thread_id();
  if (ith == 0) {
    m_mq = mq;
    m_M0 = M0;
    m_Ns = Ns;
    m_NinF = m_Nvcd * m_Ns;

    bool need_reset_shift = false;
    if(m_boundary.size() != bc.size()) {
      need_reset_shift = true;
    } else {
      for(int mu=0; mu < m_boundary.size(); ++mu){
        need_reset_shift &= (bc[mu] != m_boundary[mu]);
      }
    }

    m_boundary.resize(Ndim);
    m_boundary = bc;

    if (m_b.size() != m_Ns) {
      m_b.resize(m_Ns);
      m_c.resize(m_Ns);
      m_v1.reset(m_NinF, Nvol, 1);
      m_v2.reset(m_NinF, Nvol, 1);
      m_h1.reset(m_NinF/2, Nvol, 1);
      m_h2.reset(m_NinF/2, Nvol, 1);
      m_z1.reset(m_NinF, Nvol, 1);
      m_z2.reset(m_NinF, Nvol, 1);
      need_reset_shift = true;
    }
    if (need_reset_shift) {
      int nin = m_NinF/2; // 2 spinor
      m_shift.reset(new ShiftAField_lex<AFIELD>(nin, m_boundary));
    }
    for (int is = 0; is < m_Ns; ++is) {
      m_b[is] = b;
      m_c[is] = c;
    }

  }
#pragma omp barrier

  //- print input parameters
  vout.general(m_vl, "%s: parameters:\n", class_name.c_str());
  vout.general(m_vl, "  gamma-matrix type = %s\n", m_repr.c_str());
  vout.general(m_vl, "  mq   = %12.8f\n", mq);
  vout.general(m_vl, "  M0   = %12.8f\n", M0);
  vout.general(m_vl, "  Ns   = %4d\n", Ns);
  for (int mu = 0; mu < Ndim; ++mu) {
    vout.general(m_vl, "  boundary[%d] = %2d\n", mu, bc[mu]);
  }
  vout.general(m_vl, "  coefficients:\n");
  for (int is = 0; is < m_Ns; ++is) {
    vout.general(m_vl, "  b[%2d] = %16.10f  c[%2d] = %16.10f\n",
                 is, m_b[is], is, m_c[is]);
  }


}

//====================================================================
template<typename AFIELD>
void AForce_F_Domainwall_5din<AFIELD>::get_parameters(Parameters& params) const
{
  std::vector<double> b(m_Ns);
  std::vector<double> c(m_Ns);
  for (int is = 0; is < m_Ns; ++is) {
    b[is] = double(m_b[is]);
    c[is] = double(m_c[is]);
  }

  params.set_double("quark_mass", double(m_mq));
  params.set_double("domain_wall_height", double(m_M0));
  params.set_int("extent_of_5th_dimension", m_Ns);
  params.set_int_vector("boundary_condition", m_boundary);
  params.set_double_vector("coefficient_b", b);
  params.set_double_vector("coefficient_c", c);
  params.set_string("gamma_matrix_type", m_repr);

  params.set_string("verbose_level", vout.get_verbose_level(m_vl));
}


//====================================================================
template<typename AFIELD>
void AForce_F_Domainwall_5din<AFIELD>::set_config(Field* U)
{
  int nth = ThreadManager::get_num_threads();
  if(nth > 1){
    set_config_impl(U);
  } else {
    set_config_omp(U);
  }
  m_fopr_dw->set_config(U);

  m_U = U;

}

//====================================================================
template<typename AFIELD>
void AForce_F_Domainwall_5din<AFIELD>::set_config_omp(Field* U)
{
#pragma omp parallel
 { set_config_impl(U); }
}

//====================================================================
template<typename AFIELD>
void AForce_F_Domainwall_5din<AFIELD>::set_config_impl(Field* U)
{
#pragma omp barrier

  AIndex_lex<real_t, AFIELD::IMPL> index;
  convert_gauge(index, m_Ucp, *U);

#pragma omp barrier
}

//====================================================================
template<typename AFIELD>
void AForce_F_Domainwall_5din<AFIELD>::set_mode(const std::string& mode)
{
#pragma omp barrier

  int ith = ThreadManager::get_thread_id();
  if (ith == 0) m_mode = mode;

#pragma omp barrier
}

//====================================================================
template<typename AFIELD>
void AForce_F_Domainwall_5din<AFIELD>::force_udiv(AFIELD& force,
                                             const AFIELD& eta)
{
#pragma omp barrier

  m_fopr_dw->set_mode("H");
  m_fopr_dw->mult(m_z1, eta);

  force_udiv1_H(m_force2, m_z1, eta);
  copy(force, m_force2);

  force_udiv1_Hdag(m_force2, eta, m_z1);
  axpy(force, real_t(1.0), m_force2);

#pragma omp barrier
}

//====================================================================
template<typename AFIELD>
void AForce_F_Domainwall_5din<AFIELD>::force_udiv1(AFIELD& force,
                                              const AFIELD& zeta,
                                              const AFIELD& eta)
{
  if(m_mode=="H"){
    force_udiv1_H(force, zeta, eta);
  } else if(m_mode=="Hdag"){
    force_udiv1_Hdag(force, zeta, eta);
  } else if(m_mode=="D"){
    force_udiv1_D(force, zeta, eta);
  } else if(m_mode=="Ddag"){
    force_udiv1_Ddag(force, zeta, eta);
  } else {
    vout.crucial(m_vl, "Error at %s: irrelevant mode: %s.\n",
                 class_name.c_str(), m_mode.c_str());
    exit(EXIT_FAILURE);
  }

}


template<typename AFIELD>
void AForce_F_Domainwall_5din<AFIELD>::force_udiv1_D(AFIELD& force,
                                              const AFIELD& zeta,
                                              const AFIELD& eta)
{
#pragma omp barrier

  int Nin4 = VLEN * NVCD;
  int Nin5 = Nin4 * m_Ns;
  int Nst  = eta.nvol();
  int Nstv = Nst/VLEN;
  real_t *yp = m_v1.ptr(0);  // working vector

  int ith, nth, site0, site1;
  set_threadtask(ith, nth, site0, site1, Nstv);
  svbool_t pg = set_predicate();

  for (int site = site0; site < site1; ++site) {
    real_t* __restrict wp2 = const_cast<real_t*>(eta.ptr(0)) + Nin5 * site;
    real_t* __restrict yp2 = &yp[Nin5 * site];

    for (int is = 0; is < m_Ns; ++is) {

      int    is_up = (is + 1) % m_Ns;
      real_t Fup   = -0.25*m_c[is];
      if (is == m_Ns - 1) Fup *= -m_mq;

      int    is_dn = (is - 1 + m_Ns) % m_Ns;
      real_t Fdn   = -0.25*m_c[is];
      if (is == 0) Fdn *= -m_mq;

      real_t b1 = -0.5*m_b[is];

      for (int ic = 0; ic < NC; ++ic) {
        svreal_t vt1r, vt1i, vt2r, vt2i, vt3r, vt3i, vt4r, vt4i;

        set_aPm5_dirac_vec(vt1r, vt1i, vt2r, vt2i,
                           vt3r, vt3i, vt4r, vt4i,
                           Fup, wp2, is_up, ic);

        add_aPp5_dirac_vec(vt1r, vt1i, vt2r, vt2i,
                           vt3r, vt3i, vt4r, vt4i,
                           Fdn, wp2, is_dn, ic);

        svreal_t wt1r, wt1i, wt2r, wt2i, wt3r, wt3i, wt4r, wt4i;
        int      offset = VLEN * (2 * ND * ic + NVCD * is);

        load_vec(pg, wt1r, wp2 + offset         );
        load_vec(pg, wt1i, wp2 + offset +   VLEN);
        load_vec(pg, wt2r, wp2 + offset + 2*VLEN);
        load_vec(pg, wt2i, wp2 + offset + 3*VLEN);
        load_vec(pg, wt3r, wp2 + offset + 4*VLEN);
        load_vec(pg, wt3i, wp2 + offset + 5*VLEN);
        load_vec(pg, wt4r, wp2 + offset + 6*VLEN);
        load_vec(pg, wt4i, wp2 + offset + 7*VLEN);

        axpy_vec(pg, vt1r, b1, wt1r);
        axpy_vec(pg, vt1i, b1, wt1i);
        axpy_vec(pg, vt2r, b1, wt2r);
        axpy_vec(pg, vt2i, b1, wt2i);
        axpy_vec(pg, vt3r, b1, wt3r);
        axpy_vec(pg, vt3i, b1, wt3i);
        axpy_vec(pg, vt4r, b1, wt4r);
        axpy_vec(pg, vt4i, b1, wt4i);

        save_vec(pg, yp2 + offset,          vt1r);
        save_vec(pg, yp2 + offset +   VLEN, vt1i);
        save_vec(pg, yp2 + offset + 2*VLEN, vt2r);
        save_vec(pg, yp2 + offset + 3*VLEN, vt2i);
        save_vec(pg, yp2 + offset + 4*VLEN, vt3r);
        save_vec(pg, yp2 + offset + 5*VLEN, vt3i);
        save_vec(pg, yp2 + offset + 6*VLEN, vt4r);
        save_vec(pg, yp2 + offset + 7*VLEN, vt4i);
      }
    }
  }
#pragma omp barrier
  { // xp
    constexpr int mu=0;
    for(int site=site0; site<site1; ++site){
      real_t* __restrict in0  = m_v1.ptr(0) + site * VLEN * NVCD * m_Ns;
      real_t* __restrict out0 = m_h1.ptr(0) + site * VLEN * 2*NC*ND2 * m_Ns;
      for(int is=0; is<m_Ns; ++is){
        real_t* __restrict in  = in0  + VLEN * NVCD * is;
        real_t* __restrict out = out0 + VLEN * 2*NC*ND2 * is;
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
    } // site
    m_shift->backward(m_h2, m_h1, mu);
    for(int site=site0; site<site1; ++site){
      real_t* __restrict in0  = m_h2.ptr(0)  + site * VLEN * 2*NC*ND2 * m_Ns;
      real_t* __restrict out0 = m_v2.ptr(0)  + site * VLEN * NVCD * m_Ns;
      for(int is=0; is<m_Ns; ++is){
        real_t* __restrict in  = in0  + VLEN * 2*NC*ND2 * is;
        real_t* __restrict out = out0 + VLEN * NVCD * is;
        for(int ic=0; ic<NC; ++ic){
          int ic2 = ND * 2 * ic;
          svreal_t sp1r, sp1i, sp2r, sp2i;
          svreal_t msp1i, msp2i;
          load_vec(pg, sp1r, in         );
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
    } // site
    QXS_Spinor::tensorProd_5din(force, mu, zeta, m_v2);
  }
  { // yp
    constexpr int mu=1;
    for(int site=site0; site<site1; ++site){
      real_t* __restrict in0  = m_v1.ptr(0) + site * VLEN * NVCD * m_Ns;
      real_t* __restrict out0 = m_h1.ptr(0) + site * VLEN * 2*NC*ND2 * m_Ns;
      for(int is=0; is<m_Ns; ++is){
        real_t* __restrict in  = in0  + VLEN * NVCD * is;
        real_t* __restrict out = out0 + VLEN * 2*NC*ND2 * is;
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
    } // site
    m_shift->backward(m_h2, m_h1, mu);
    for(int site=site0; site<site1; ++site){
      real_t* __restrict in0  = m_h2.ptr(0)  + site * VLEN * 2*NC*ND2 * m_Ns;
      real_t* __restrict out0 = m_v2.ptr(0)  + site * VLEN * NVCD * m_Ns;
      for(int is=0; is<m_Ns; ++is){
        real_t* __restrict in  = in0  + VLEN * 2*NC*ND2 * is;
        real_t* __restrict out = out0 + VLEN * NVCD * is;
        for(int ic=0; ic<NC; ++ic){
          int ic2 = ND * 2 * ic;
          svreal_t sp1r, sp1i, sp2r, sp2i;
          svreal_t msp1r, msp1i;
          load_vec(pg, sp1r, in         );
          load_vec(pg, sp1i, in +   VLEN);
          load_vec(pg, sp2r, in + 2*VLEN);
          load_vec(pg, sp2i, in + 3*VLEN);
          in  += 4*VLEN;

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
    }
    QXS_Spinor::tensorProd_5din(force, mu, zeta, m_v2);
  }
  { // zp
    constexpr int mu=2;
    for(int site=site0; site<site1; ++site){
      real_t* __restrict in0  = m_v1.ptr(0) + site * VLEN * NVCD * m_Ns;
      real_t* __restrict out0 = m_h1.ptr(0) + site * VLEN * 2*NC*ND2 * m_Ns;
      for(int is=0; is<m_Ns; ++is){
        real_t* __restrict in  = in0  + VLEN * NVCD * is;
        real_t* __restrict out = out0 + VLEN * 2*NC*ND2 * is;
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
    } // site
    m_shift->backward(m_h2, m_h1, mu);
    for(int site=site0; site<site1; ++site){
      real_t* __restrict in0  = m_h2.ptr(0)  + site * VLEN * 2*NC*ND2 * m_Ns;
      real_t* __restrict out0 = m_v2.ptr(0)  + site * VLEN * NVCD * m_Ns;
      for(int is=0; is<m_Ns; ++is){
        real_t* __restrict in  = in0  + VLEN * 2*NC*ND2 * is;
        real_t* __restrict out = out0 + VLEN * NVCD * is;
        for(int ic=0; ic<NC; ++ic){
          int ic2 = ND * 2 * ic;
          svreal_t sp1r, sp1i, sp2r, sp2i;
          svreal_t msp1i, msp2r;
          load_vec(pg, sp1r, in         );
          load_vec(pg, sp1i, in +   VLEN);
          load_vec(pg, sp2r, in + 2*VLEN);
          load_vec(pg, sp2i, in + 3*VLEN);
          in  += 4*VLEN;

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
    }
    QXS_Spinor::tensorProd_5din(force, mu, zeta, m_v2);
  }
  { // tp (Dirac rep.)
    constexpr int mu=3;
    for(int site=site0; site<site1; ++site){
      real_t* __restrict in0  = m_v1.ptr(0) + site * VLEN * NVCD * m_Ns;
      real_t* __restrict out0 = m_h1.ptr(0) + site * VLEN * 2*NC*ND2 * m_Ns;
      for(int is=0; is<m_Ns; ++is){
        real_t* __restrict in  = in0  + VLEN * NVCD * is;
        real_t* __restrict out = out0 + VLEN * 2*NC*ND2 * is;
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
    } // site
    m_shift->backward(m_h2, m_h1, mu);
    for(int site=site0; site<site1; ++site){
      real_t* __restrict in0  = m_h2.ptr(0)  + site * VLEN * 2*NC*ND2 * m_Ns;
      real_t* __restrict out0 = m_v2.ptr(0)  + site * VLEN * NVCD * m_Ns;
      for(int is=0; is<m_Ns; ++is){
        real_t* __restrict in  = in0  + VLEN * 2*NC*ND2 * is;
        real_t* __restrict out = out0 + VLEN * NVCD * is;
        for(int ic=0; ic<NC; ++ic){
          int ic2 = ND * 2 * ic;
          svreal_t sp1r, sp1i, sp2r, sp2i;
          svreal_t vzero;

          load_vec(pg, sp1r, in         );
          load_vec(pg, sp1i, in +   VLEN);
          load_vec(pg, sp2r, in + 2*VLEN);
          load_vec(pg, sp2i, in + 3*VLEN);
          in  += 4*VLEN;

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
    }
    QXS_Spinor::tensorProd_5din(force, mu, zeta, m_v2);
  }
#pragma omp barrier

}


template<typename AFIELD>
void AForce_F_Domainwall_5din<AFIELD>::force_udiv1_Ddag(AFIELD& force,
                                                        const AFIELD& zeta,
                                                        const AFIELD& eta)
{
#pragma omp barrier

  int Nin4 = VLEN * NVCD;
  int Nin5 = Nin4 * m_Ns;
  int Nst  = eta.nvol();
  int Nstv = Nst/VLEN;
  real_t *yp = m_v1.ptr(0);  // working vector

  int ith, nth, site0, site1;
  set_threadtask(ith, nth, site0, site1, Nstv);
  svbool_t pg = set_predicate();

  for (int site = site0; site < site1; ++site) { // Dirac rep.
    real_t* __restrict wp2 = const_cast<real_t*>(zeta.ptr(0)) + Nin5 * site;
    real_t* __restrict yp2 = &yp[Nin5 * site];

    for (int is = 0; is < m_Ns; ++is) {

      int    is_up = (is + 1) % m_Ns;
      real_t Fup   = -0.25*m_c[is];
      if (is == m_Ns - 1) Fup *= -m_mq;

      int    is_dn = (is - 1 + m_Ns) % m_Ns;
      real_t Fdn   = -0.25*m_c[is];
      if (is == 0) Fdn *= -m_mq;

      real_t b1 = -0.5*m_b[is];

      for (int ic = 0; ic < NC; ++ic) {
        svreal_t vt1r, vt1i, vt2r, vt2i, vt3r, vt3i, vt4r, vt4i;

        set_aPm5_dirac_vec(vt1r, vt1i, vt2r, vt2i,
                           vt3r, vt3i, vt4r, vt4i,
                           Fup, wp2, is_up, ic);

        add_aPp5_dirac_vec(vt1r, vt1i, vt2r, vt2i,
                           vt3r, vt3i, vt4r, vt4i,
                           Fdn, wp2, is_dn, ic);

        svreal_t wt1r, wt1i, wt2r, wt2i, wt3r, wt3i, wt4r, wt4i;
        int      offset = VLEN * (2 * ND * ic + NVCD * is);

        load_vec(pg, wt1r, wp2 + offset         );
        load_vec(pg, wt1i, wp2 + offset +   VLEN);
        load_vec(pg, wt2r, wp2 + offset + 2*VLEN);
        load_vec(pg, wt2i, wp2 + offset + 3*VLEN);
        load_vec(pg, wt3r, wp2 + offset + 4*VLEN);
        load_vec(pg, wt3i, wp2 + offset + 5*VLEN);
        load_vec(pg, wt4r, wp2 + offset + 6*VLEN);
        load_vec(pg, wt4i, wp2 + offset + 7*VLEN);

        axpy_vec(pg, vt1r, b1, wt1r);
        axpy_vec(pg, vt1i, b1, wt1i);
        axpy_vec(pg, vt2r, b1, wt2r);
        axpy_vec(pg, vt2i, b1, wt2i);
        axpy_vec(pg, vt3r, b1, wt3r);
        axpy_vec(pg, vt3i, b1, wt3i);
        axpy_vec(pg, vt4r, b1, wt4r);
        axpy_vec(pg, vt4i, b1, wt4i);

        save_vec(pg, yp2 + offset,          vt1r);
        save_vec(pg, yp2 + offset +   VLEN, vt1i);
        save_vec(pg, yp2 + offset + 2*VLEN, vt2r);
        save_vec(pg, yp2 + offset + 3*VLEN, vt2i);
        save_vec(pg, yp2 + offset + 4*VLEN, vt3r);
        save_vec(pg, yp2 + offset + 5*VLEN, vt3i);
        save_vec(pg, yp2 + offset + 6*VLEN, vt4r);
        save_vec(pg, yp2 + offset + 7*VLEN, vt4i);
      }
    }
  }
#pragma omp barrier
  { // xp
    constexpr int mu=0;
    for(int site=site0; site<site1; ++site){
      real_t* __restrict in0  = const_cast<real_t*>(eta.ptr(0)) + site * VLEN * NVCD * m_Ns;
      real_t* __restrict out0 = m_h1.ptr(0) + site * VLEN * 2*NC*ND2 * m_Ns;
      for(int is=0; is<m_Ns; ++is){
        real_t* __restrict in  = in0  + VLEN * NVCD * is;
        real_t* __restrict out = out0 + VLEN * 2*NC*ND2 * is;
        for(int ic=0; ic<NC; ++ic){
          svreal_t sp1r, sp1i, sp2r, sp2i;
          set_sp2_xm(pg, sp1r, sp1i, sp2r, sp2i, in, ic);
          save_vec(pg, out,          sp1r);
          save_vec(pg, out +   VLEN, sp1i);
          save_vec(pg, out + 2*VLEN, sp2r);
          save_vec(pg, out + 3*VLEN, sp2i);
          out += 4*VLEN;
        }
      }
    } // site
    m_shift->backward(m_h2, m_h1, mu);
    for(int site=site0; site<site1; ++site){
      real_t* __restrict in0  = m_h2.ptr(0)  + site * VLEN * 2*NC*ND2 * m_Ns;
      real_t* __restrict out0 = m_v2.ptr(0)  + site * VLEN * NVCD * m_Ns;
      for(int is=0; is<m_Ns; ++is){
        real_t* __restrict in  = in0  + VLEN * 2*NC*ND2 * is;
        real_t* __restrict out = out0 + VLEN * NVCD * is;
        for(int ic=0; ic<NC; ++ic){
          int ic2 = ND * 2 * ic;
          svreal_t sp1r, sp1i, sp2r, sp2i;
          svreal_t msp1r, msp2r;
          load_vec(pg, sp1r, in         );
          load_vec(pg, sp1i, in +   VLEN);
          load_vec(pg, sp2r, in + 2*VLEN);
          load_vec(pg, sp2i, in + 3*VLEN);
          in  += 4*VLEN;

          flip_sign(pg, msp1r, sp1r);
          flip_sign(pg, msp2r, sp2r);
          save_vec(pg, out + VLEN*(ic2   + ID1),  sp1r);  // sp1
          save_vec(pg, out + VLEN*(ic2+1 + ID1),  sp1i);
          save_vec(pg, out + VLEN*(ic2   + ID2),  sp2r);  // sp2
          save_vec(pg, out + VLEN*(ic2+1 + ID2),  sp2i);
          save_vec(pg, out + VLEN*(ic2   + ID3),  sp2i);  // i x sp2
          save_vec(pg, out + VLEN*(ic2+1 + ID3), msp2r);
          save_vec(pg, out + VLEN*(ic2   + ID4),  sp1i);  // i x sp1
          save_vec(pg, out + VLEN*(ic2+1 + ID4), msp1r);
        }
      }
    } // site
    QXS_Spinor::tensorProd_5din(force, mu, m_v1, m_v2);
  }
  { // yp
    constexpr int mu=1;
    for(int site=site0; site<site1; ++site){
      real_t* __restrict in0  = const_cast<real_t*>(eta.ptr(0)) + site * VLEN * NVCD * m_Ns;
      real_t* __restrict out0 = m_h1.ptr(0) + site * VLEN * 2*NC*ND2 * m_Ns;
      for(int is=0; is<m_Ns; ++is){
        real_t* __restrict in  = in0  + VLEN * NVCD * is;
        real_t* __restrict out = out0 + VLEN * 2*NC*ND2 * is;
        for(int ic=0; ic<NC; ++ic){
          svreal_t sp1r, sp1i, sp2r, sp2i;
          set_sp2_ym(pg, sp1r, sp1i, sp2r, sp2i, in, ic);
          save_vec(pg, out,          sp1r);
          save_vec(pg, out +   VLEN, sp1i);
          save_vec(pg, out + 2*VLEN, sp2r);
          save_vec(pg, out + 3*VLEN, sp2i);
          out += 4*VLEN;
        }
      }
    } // site
    m_shift->backward(m_h2, m_h1, mu);
    for(int site=site0; site<site1; ++site){
      real_t* __restrict in0  = m_h2.ptr(0)  + site * VLEN * 2*NC*ND2 * m_Ns;
      real_t* __restrict out0 = m_v2.ptr(0)  + site * VLEN * NVCD * m_Ns;
      for(int is=0; is<m_Ns; ++is){
        real_t* __restrict in  = in0  + VLEN * 2*NC*ND2 * is;
        real_t* __restrict out = out0 + VLEN * NVCD * is;
        for(int ic=0; ic<NC; ++ic){
          int ic2 = ND * 2 * ic;
          svreal_t sp1r, sp1i, sp2r, sp2i;
          svreal_t msp2r, msp2i;
          load_vec(pg, sp1r, in         );
          load_vec(pg, sp1i, in +   VLEN);
          load_vec(pg, sp2r, in + 2*VLEN);
          load_vec(pg, sp2i, in + 3*VLEN);
          in  += 4*VLEN;

          flip_sign(pg, msp2r, sp2r);
          flip_sign(pg, msp2i, sp2i);
          save_vec(pg, out + VLEN*(ic2   + ID1),  sp1r);  // sp1
          save_vec(pg, out + VLEN*(ic2+1 + ID1),  sp1i);
          save_vec(pg, out + VLEN*(ic2   + ID2),  sp2r);  // sp2
          save_vec(pg, out + VLEN*(ic2+1 + ID2),  sp2i);
          save_vec(pg, out + VLEN*(ic2   + ID3), msp2r);  // -sp2
          save_vec(pg, out + VLEN*(ic2+1 + ID3), msp2i);
          save_vec(pg, out + VLEN*(ic2   + ID4),  sp1r);  // sp1
          save_vec(pg, out + VLEN*(ic2+1 + ID4),  sp1i);
        }
      }
    }
    QXS_Spinor::tensorProd_5din(force, mu, m_v1, m_v2);
  }
  { // zp
    constexpr int mu=2;
    for(int site=site0; site<site1; ++site){
      real_t* __restrict in0  = const_cast<real_t*>(eta.ptr(0)) + site * VLEN * NVCD * m_Ns;
      real_t* __restrict out0 = m_h1.ptr(0) + site * VLEN * 2*NC*ND2 * m_Ns;
      for(int is=0; is<m_Ns; ++is){
        real_t* __restrict in  = in0  + VLEN * NVCD * is;
        real_t* __restrict out = out0 + VLEN * 2*NC*ND2 * is;
        for(int ic=0; ic<NC; ++ic){
          svreal_t sp1r, sp1i, sp2r, sp2i;
          set_sp2_zm(pg, sp1r, sp1i, sp2r, sp2i, in, ic);
          save_vec(pg, out,          sp1r);
          save_vec(pg, out +   VLEN, sp1i);
          save_vec(pg, out + 2*VLEN, sp2r);
          save_vec(pg, out + 3*VLEN, sp2i);
          out += 4*VLEN;
        }
      }
    } // site
    m_shift->backward(m_h2, m_h1, mu);
    for(int site=site0; site<site1; ++site){
      real_t* __restrict in0  = m_h2.ptr(0)  + site * VLEN * 2*NC*ND2 * m_Ns;
      real_t* __restrict out0 = m_v2.ptr(0)  + site * VLEN * NVCD * m_Ns;
      for(int is=0; is<m_Ns; ++is){
        real_t* __restrict in  = in0  + VLEN * 2*NC*ND2 * is;
        real_t* __restrict out = out0 + VLEN * NVCD * is;
        for(int ic=0; ic<NC; ++ic){
          int ic2 = ND * 2 * ic;
          svreal_t sp1r, sp1i, sp2r, sp2i;
          svreal_t msp1r, msp2i;
          load_vec(pg, sp1r, in         );
          load_vec(pg, sp1i, in +   VLEN);
          load_vec(pg, sp2r, in + 2*VLEN);
          load_vec(pg, sp2i, in + 3*VLEN);
          in  += 4*VLEN;

          flip_sign(pg, msp1r, sp1r);
          flip_sign(pg, msp2i, sp2i);
          save_vec(pg, out + VLEN*(ic2   + ID1),  sp1r);  // sp1
          save_vec(pg, out + VLEN*(ic2+1 + ID1),  sp1i);
          save_vec(pg, out + VLEN*(ic2   + ID2),  sp2r);  // sp2
          save_vec(pg, out + VLEN*(ic2+1 + ID2),  sp2i);
          save_vec(pg, out + VLEN*(ic2   + ID3),  sp1i);  // -i x sp1
          save_vec(pg, out + VLEN*(ic2+1 + ID3), msp1r);
          save_vec(pg, out + VLEN*(ic2   + ID4), msp2i);  //  i x sp2
          save_vec(pg, out + VLEN*(ic2+1 + ID4),  sp2r);
        }
      }
    }
    QXS_Spinor::tensorProd_5din(force, mu, m_v1, m_v2);
  }
  { // tp (Dirac rep.)
    constexpr int mu=3;
    for(int site=site0; site<site1; ++site){
      real_t* __restrict in0  = const_cast<real_t*>(eta.ptr(0)) + site * VLEN * NVCD * m_Ns;
      real_t* __restrict out0 = m_h1.ptr(0) + site * VLEN * 2*NC*ND2 * m_Ns;
      for(int is=0; is<m_Ns; ++is){
        real_t* __restrict in  = in0  + VLEN * NVCD * is;
        real_t* __restrict out = out0 + VLEN * 2*NC*ND2 * is;
        for(int ic=0; ic<NC; ++ic){
          svreal_t sp1r, sp1i, sp2r, sp2i;
          set_sp2_tm_dirac(pg, sp1r, sp1i, sp2r, sp2i, in, ic);
          save_vec(pg, out,          sp1r);
          save_vec(pg, out +   VLEN, sp1i);
          save_vec(pg, out + 2*VLEN, sp2r);
          save_vec(pg, out + 3*VLEN, sp2i);
          out += 4*VLEN;
        }
      }
    } // site
    m_shift->backward(m_h2, m_h1, mu);
    for(int site=site0; site<site1; ++site){
      real_t* __restrict in0  = m_h2.ptr(0)  + site * VLEN * 2*NC*ND2 * m_Ns;
      real_t* __restrict out0 = m_v2.ptr(0)  + site * VLEN * NVCD * m_Ns;
      for(int is=0; is<m_Ns; ++is){
        real_t* __restrict in  = in0  + VLEN * 2*NC*ND2 * is;
        real_t* __restrict out = out0 + VLEN * NVCD * is;
        for(int ic=0; ic<NC; ++ic){
          int ic2 = ND * 2 * ic;
          svreal_t sp1r, sp1i, sp2r, sp2i;
          svreal_t vzero;

          load_vec(pg, sp1r, in         );
          load_vec(pg, sp1i, in +   VLEN);
          load_vec(pg, sp2r, in + 2*VLEN);
          load_vec(pg, sp2i, in + 3*VLEN);
          in  += 4*VLEN;

          clear_vec(pg, vzero);

          save_vec(pg, out + VLEN*(ic2   + ID1),  sp1r);  // sp1
          save_vec(pg, out + VLEN*(ic2+1 + ID1),  sp1i);
          save_vec(pg, out + VLEN*(ic2   + ID2),  sp2r);  // sp2
          save_vec(pg, out + VLEN*(ic2+1 + ID2),  sp2i);
          save_vec(pg, out + VLEN*(ic2   + ID3), vzero);  // 0
          save_vec(pg, out + VLEN*(ic2+1 + ID3), vzero);
          save_vec(pg, out + VLEN*(ic2   + ID4), vzero);  // 0
          save_vec(pg, out + VLEN*(ic2+1 + ID4), vzero);
        }
      }
    }
    QXS_Spinor::tensorProd_5din(force, mu, m_v1, m_v2);
  }
#pragma omp barrier

}


//====================================================================
template<typename AFIELD>
void AForce_F_Domainwall_5din<AFIELD>::force_udiv1_H(AFIELD& force,
                                              const AFIELD& zeta,
                                              const AFIELD& eta)
{
  m_fopr_dw->mult_gm5R(m_z2, zeta);
  force_udiv1_D(force, m_z2, eta);

}



//====================================================================
template<typename AFIELD>
void AForce_F_Domainwall_5din<AFIELD>::force_udiv1_Hdag(AFIELD& force,
                                                   const AFIELD& zeta,
                                                   const AFIELD& eta)
{
  m_fopr_dw->mult_gm5R(m_z2, eta);
  force_udiv1_Ddag(force, zeta, m_z2);
}
//============================================================END=====
