/*!
        @file    aforce_F_Domainwall_5din-tmpl.h

        @brief

        @author  Issaku Kanamori (kanamori)
                 $LastChangedBy: kanamori $

        @date    $LastChangedDate: 2025-12-03 19:35:35 +0900 (2025年12月03日 (水)) $

        @version $LastChangedRevision: 2668 $
*/

#include "lib_alt_Accel/Force/Fermion/aforce_F_Domainwall_5din.h"
#include "lib/Force/Fermion/force_F_Domainwall.h"
#include "lib/Force/Fermion/tensorProd.h"
#include "lib_alt_Accel/Field/aindex_lex.h"

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
  //  int nth = ThreadManager::get_num_threads();
  //  if(nth > 1){
  //    set_config_impl(U);
  //  } else {
  //    set_config_omp(U);
  //  }
  m_fopr_dw->set_config(U);

  m_U = U;

}

//====================================================================
template<typename AFIELD>
void AForce_F_Domainwall_5din<AFIELD>::set_config_omp(Field* U)
{
  //#pragma omp parallel
  // { set_config_impl(U); }
}

//====================================================================
template<typename AFIELD>
void AForce_F_Domainwall_5din<AFIELD>::set_config_impl(Field* U)
{
  //#pragma omp barrier
  //
  //  AIndex_lex<real_t, AFIELD::IMPL> index;
  //  convert_gauge(index, m_Ucp, *U);
  //
  //#pragma omp barrier
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

  int ith = ThreadManager::get_thread_id();
  int nth = ThreadManager::get_num_threads();
  int ith_kernel = 0;
  //if(nth > 1) ith_kernel = 1;

  real_t *wp = const_cast<AFIELD *>(&eta)->ptr(0);

  real_t *yp = m_v1.ptr(0);    // working vector: 4-spinor
  real_t *vp = m_v2.ptr(0);    // working vector: 4-spinor
  real_t *h1p = m_h1.ptr(0);   // working vector: 2-spinor
  real_t *h2p = m_h2.ptr(0);   // working vector: 2-spinor

  int Nst=force.nvol();
  int Nst_pad = CEIL_NWP(Nst);
  int Ns = m_Ns;
  int Nin5 = NVCD * Ns;
  int size = Nin5 * Nst_pad;
  int size2 = size/2;  // for 2-spinor
  real_t *b = &m_b[0];
  real_t *c = &m_c[0];
  real_t mq = m_mq;
  //  real_t M0 = m_M0;
  int Nx = CommonParameters::Nx();
  int Ny = CommonParameters::Ny();
  int Nz = CommonParameters::Nz();
  int Nt = CommonParameters::Nt();

  if (ith == ith_kernel){

#pragma acc data present(yp[0:size], wp[0:size])  \
  copyin(Nst, Nst_pad, Ns, Nin5, mq, b[0:Ns], c[0:Ns])
    {
      #pragma acc parallel num_workers(NUM_WORKERS) vector_length(VECTOR_LENGTH)
      {
        #pragma acc loop gang worker vector
        for (int idx = 0; idx < NVC * Nst_pad; ++idx) {
          // idx = idx_in + NWP * (ivc + NVC*idx_out)
          int idx2_wp = idx / NWP;
          int idx_in  = idx % NWP;
          int idx_out = idx2_wp / NVC;
          int site = idx_in + NWP * idx_out;
          int ivc  = idx2_wp % NVC;  // color and r/i

          if(site < Nst){
            for (int is = 0; is < Ns; ++is) {

              real_t wt1, wt2, wt3, wt4;
              real_t vt1, vt2, vt3, vt4;
              real_t xt1, xt2, xt3, xt4;

              int is_up = (is+1) % Ns;
              real_t Fup   = 0.5;
              if (is == Ns-1) Fup = -0.5 * mq;

              int ivcs_up = ivc + NVCD * is_up;
              wt1 = wp[IDX2(Nin5, ID1 + ivcs_up, site)];
              wt2 = wp[IDX2(Nin5, ID2 + ivcs_up, site)];
              wt3 = wp[IDX2(Nin5, ID3 + ivcs_up, site)];
              wt4 = wp[IDX2(Nin5, ID4 + ivcs_up, site)];

              vt1 = Fup * (wt1 - wt3);
              vt2 = Fup * (wt2 - wt4);
              vt3 = Fup * (wt3 - wt1);
              vt4 = Fup * (wt4 - wt2);

              int is_dn = (is-1 + Ns) % Ns;
              real_t Fdn   = 0.5;
              if (is == 0) { Fdn = -0.5 * mq; }

              int ivcs_dn = ivc + NVCD * is_dn;
              wt1 = wp[IDX2(Nin5, ID1 + ivcs_dn, site)];
              wt2 = wp[IDX2(Nin5, ID2 + ivcs_dn, site)];
              wt3 = wp[IDX2(Nin5, ID3 + ivcs_dn, site)];
              wt4 = wp[IDX2(Nin5, ID4 + ivcs_dn, site)];

              vt1 += Fdn * (wt1 + wt3);
              vt2 += Fdn * (wt2 + wt4);
              vt3 += Fdn * (wt3 + wt1);
              vt4 += Fdn * (wt4 + wt2);

              int ivcs = ivc + NVCD * is;
              wt1 = wp[IDX2(Nin5, ID1 + ivcs, site)];
              wt2 = wp[IDX2(Nin5, ID2 + ivcs, site)];
              wt3 = wp[IDX2(Nin5, ID3 + ivcs, site)];
              wt4 = wp[IDX2(Nin5, ID4 + ivcs, site)];

              real_t B2 = -0.5 * b[is];
              real_t C2 = -0.5 * c[is];

              yp[IDX2(Nin5, ID1 + ivcs, site)] = B2 * wt1 + C2 * vt1;
              yp[IDX2(Nin5, ID2 + ivcs, site)] = B2 * wt2 + C2 * vt2;
              yp[IDX2(Nin5, ID3 + ivcs, site)] = B2 * wt3 + C2 * vt3;
              yp[IDX2(Nin5, ID4 + ivcs, site)] = B2 * wt4 + C2 * vt4;
            } // is
          } // ist<Nst
        } // idx

      } // acc parallel
    } // acc data
  } //ith==ith_kernel


#pragma omp barrier
  // indexing of fermion field
  // index for the field
  //   idx_in + NMP * (r/i + 2*ic + 6*id)
  //          + NWP * 24*is
  //          + NWP * 24*Ns*idx_out

  ///////////////////////////////////////////////
  int mu=0;   // xp
  ///////////////////////////////////////////////
  if (ith == ith_kernel){ //  make 2 spinor
#pragma acc data present(yp[0:size], h1p[0:size2]) \
                 copyin(Ns, Nst_pad)
    {
#pragma acc parallel \
            num_workers(NUM_WORKERS) vector_length(VECTOR_LENGTH)
      {
#pragma acc loop
        for (int idx = 0; idx < Ns * NVC * Nst_pad; ++idx) {

          int idx2_wp = idx / NWP;
          int idx_in  = idx % NWP;
          int idx_out = idx2_wp / Ns;
          int is      = idx2_wp % Ns;
          int site = idx_in + NWP * idx_out;

          if(site < Nst){
            for(int ic=0; ic<NC; ++ic) {
              real_t vt1r, vt1i, vt2r, vt2i;
              vt1r = yp[IDX2_SP_5D_R(ic, 0, is, Ns, site)] - yp[IDX2_SP_5D_I(ic, 3, is, Ns, site)];
              vt1i = yp[IDX2_SP_5D_I(ic, 0, is, Ns, site)] + yp[IDX2_SP_5D_R(ic, 3, is, Ns, site)];
              vt2r = yp[IDX2_SP_5D_R(ic, 1, is, Ns, site)] - yp[IDX2_SP_5D_I(ic, 2, is, Ns, site)];
              vt2i = yp[IDX2_SP_5D_I(ic, 1, is, Ns, site)] + yp[IDX2_SP_5D_R(ic, 2, is, Ns, site)];
              h1p[IDX2(NVCD2*Ns, (2*ic+0 + NVC*(0 + 2*is)), site)] = vt1r;
              h1p[IDX2(NVCD2*Ns, (2*ic+1 + NVC*(0 + 2*is)), site)] = vt1i;
              h1p[IDX2(NVCD2*Ns, (2*ic+0 + NVC*(1 + 2*is)), site)] = vt2r;
              h1p[IDX2(NVCD2*Ns, (2*ic+1 + NVC*(1 + 2*is)), site)] = vt2i;
            } // ic
          } // site < NST
        } // idx
      }
    }
  } // ith==kernel

  m_shift->backward(m_h2, m_h1, mu);

  if (ith == ith_kernel){ // accumlate to 4 spinor
#pragma acc data present(vp[0:size], h2p[0:size2]) \
  copyin(Ns, Nst_pad)
    {
#pragma acc parallel \
            num_workers(NUM_WORKERS) vector_length(VECTOR_LENGTH)
      {
#pragma acc loop
        for (int idx = 0; idx < Ns * Nst_pad; ++idx) {
          int idx2_wp = idx / NWP;
          int idx_in  = idx % NWP;
          int idx_out = idx2_wp / Ns;
          int is      = idx2_wp % Ns;
          int site = idx_in + NWP * idx_out;

          if(site < Nst){
            for(int ic=0; ic<NC; ++ic) {
              real_t vt1r, vt1i, vt2r, vt2i;
              vt1r = h2p[IDX2(NVCD2*Ns, (2*ic+0 + NVC*(0 + 2*is)), site)];
              vt1i = h2p[IDX2(NVCD2*Ns, (2*ic+1 + NVC*(0 + 2*is)), site)];
              vt2r = h2p[IDX2(NVCD2*Ns, (2*ic+0 + NVC*(1 + 2*is)), site)];
              vt2i = h2p[IDX2(NVCD2*Ns, (2*ic+1 + NVC*(1 + 2*is)), site)];

              vp[IDX2_SP_5D_R(ic, 0, is, Ns, site)] =  vt1r; // sp1
              vp[IDX2_SP_5D_I(ic, 0, is, Ns, site)] =  vt1i;
              vp[IDX2_SP_5D_R(ic, 1, is, Ns, site)] =  vt2r; // sp2
              vp[IDX2_SP_5D_I(ic, 1, is, Ns, site)] =  vt2i;
              vp[IDX2_SP_5D_R(ic, 2, is, Ns, site)] =  vt2i; // i x sp2
              vp[IDX2_SP_5D_I(ic, 2, is, Ns, site)] = -vt2r;
              vp[IDX2_SP_5D_R(ic, 3, is, Ns, site)] =  vt1i; // i x sp1
              vp[IDX2_SP_5D_I(ic, 3, is, Ns, site)] = -vt1r;
            } // ic
          } // site < NST
        } // idx
      }
    }
  } // ith==kernel
  Accel_Spinor::tensorProd_5din(force, mu, zeta, m_v2);

  ///////////////////////////////////////////////
  mu=1; // yp
  ///////////////////////////////////////////////
  if (ith == ith_kernel){ // make 2 spinor
#pragma acc data present(yp[0:size], h1p[0:size2]) \
                 copyin(Ns, Nst_pad)
    {
#pragma acc parallel \
            num_workers(NUM_WORKERS) vector_length(VECTOR_LENGTH)
      {
#pragma acc loop
        for (int idx = 0; idx < Ns * NVC * Nst_pad; ++idx) {
          // idx = idx_in + NWP* ( is + Ns*idx_out)

          int idx2_wp = idx / NWP;
          int idx_in  = idx % NWP;
          int idx_out = idx2_wp / Ns;
          int is      = idx2_wp % Ns;
          int site = idx_in + NWP * idx_out;

          if(site < Nst){
            for(int ic=0; ic<NC; ++ic) {
              real_t vt1r, vt1i, vt2r, vt2i;
              vt1r = yp[IDX2_SP_5D_R(ic, 0, is, Ns, site)] + yp[IDX2_SP_5D_R(ic, 3, is, Ns, site)];
              vt1i = yp[IDX2_SP_5D_I(ic, 0, is, Ns, site)] + yp[IDX2_SP_5D_I(ic, 3, is, Ns, site)];
              vt2r = yp[IDX2_SP_5D_R(ic, 1, is, Ns, site)] - yp[IDX2_SP_5D_R(ic, 2, is, Ns, site)];
              vt2i = yp[IDX2_SP_5D_I(ic, 1, is, Ns, site)] - yp[IDX2_SP_5D_I(ic, 2, is, Ns, site)];
              h1p[IDX2(NVCD2*Ns, (2*ic+0 + NVC*(0 + 2*is)), site)] = vt1r;
              h1p[IDX2(NVCD2*Ns, (2*ic+1 + NVC*(0 + 2*is)), site)] = vt1i;
              h1p[IDX2(NVCD2*Ns, (2*ic+0 + NVC*(1 + 2*is)), site)] = vt2r;
              h1p[IDX2(NVCD2*Ns, (2*ic+1 + NVC*(1 + 2*is)), site)] = vt2i;
            } // ic
          } // site < NST
        } // idx
      }
    }
  } // ith==kernel

  m_shift->backward(m_h2, m_h1, mu);

  if (ith == ith_kernel){ // accumlate to 4 spinor
#pragma acc data present(vp[0:size], h2p[0:size2]) \
  copyin(Ns, Nst_pad)
    {
#pragma acc parallel \
            num_workers(NUM_WORKERS) vector_length(VECTOR_LENGTH)
      {
#pragma acc loop
        for (int idx = 0; idx < Ns * Nst_pad; ++idx) {
          int idx2_wp = idx / NWP;
          int idx_in  = idx % NWP;
          int idx_out = idx2_wp / Ns;
          int is      = idx2_wp % Ns;
          int site = idx_in + NWP * idx_out;

          if(site < Nst){
            for(int ic=0; ic<NC; ++ic) {
              real_t vt1r, vt1i, vt2r, vt2i;
              vt1r = h2p[IDX2(NVCD2*Ns, (2*ic+0 + NVC*(0 + 2*is)), site)];
              vt1i = h2p[IDX2(NVCD2*Ns, (2*ic+1 + NVC*(0 + 2*is)), site)];
              vt2r = h2p[IDX2(NVCD2*Ns, (2*ic+0 + NVC*(1 + 2*is)), site)];
              vt2i = h2p[IDX2(NVCD2*Ns, (2*ic+1 + NVC*(1 + 2*is)), site)];

              vp[IDX2_SP_5D_R(ic, 0, is, Ns, site)] =  vt1r; // sp1
              vp[IDX2_SP_5D_I(ic, 0, is, Ns, site)] =  vt1i;
              vp[IDX2_SP_5D_R(ic, 1, is, Ns, site)] =  vt2r; // sp2
              vp[IDX2_SP_5D_I(ic, 1, is, Ns, site)] =  vt2i;
              vp[IDX2_SP_5D_R(ic, 2, is, Ns, site)] = -vt2r; // -sp2
              vp[IDX2_SP_5D_I(ic, 2, is, Ns, site)] = -vt2i;
              vp[IDX2_SP_5D_R(ic, 3, is, Ns, site)] =  vt1r; //  sp1
              vp[IDX2_SP_5D_I(ic, 3, is, Ns, site)] =  vt1i;
            } // ic
          } // site < NST
        } // idx
      }
    }
  } // ith==kernel
  Accel_Spinor::tensorProd_5din(force, mu, zeta, m_v2);

  ///////////////////////////////////////////////
  mu=2;        // zp
  ///////////////////////////////////////////////
  if (ith == ith_kernel){ //  make 2 spinor
#pragma acc data present(yp[0:size], h1p[0:size2]) \
                 copyin(Ns, Nst_pad)
    {
#pragma acc parallel \
            num_workers(NUM_WORKERS) vector_length(VECTOR_LENGTH)
      {
#pragma acc loop
        for (int idx = 0; idx < Ns * NVC * Nst_pad; ++idx) {

          int idx2_wp = idx / NWP;
          int idx_in  = idx % NWP;
          int idx_out = idx2_wp / Ns;
          int is      = idx2_wp % Ns;
          int site = idx_in + NWP * idx_out;

          if(site < Nst){
            for(int ic=0; ic<NC; ++ic) {
              real_t vt1r, vt1i, vt2r, vt2i;
              vt1r = yp[IDX2_SP_5D_R(ic, 0, is, Ns, site)] - yp[IDX2_SP_5D_I(ic, 2, is, Ns, site)];
              vt1i = yp[IDX2_SP_5D_I(ic, 0, is, Ns, site)] + yp[IDX2_SP_5D_R(ic, 2, is, Ns, site)];
              vt2r = yp[IDX2_SP_5D_R(ic, 1, is, Ns, site)] + yp[IDX2_SP_5D_I(ic, 3, is, Ns, site)];
              vt2i = yp[IDX2_SP_5D_I(ic, 1, is, Ns, site)] - yp[IDX2_SP_5D_R(ic, 3, is, Ns, site)];
              h1p[IDX2(NVCD2*Ns, (2*ic+0 + NVC*(0 + 2*is)), site)] = vt1r;
              h1p[IDX2(NVCD2*Ns, (2*ic+1 + NVC*(0 + 2*is)), site)] = vt1i;
              h1p[IDX2(NVCD2*Ns, (2*ic+0 + NVC*(1 + 2*is)), site)] = vt2r;
              h1p[IDX2(NVCD2*Ns, (2*ic+1 + NVC*(1 + 2*is)), site)] = vt2i;
            } // ic
          } // site < NST
        } // idx
      }
    }
  } // ith==kernel

  m_shift->backward(m_h2, m_h1, mu);

  if (ith == ith_kernel){ // accumlate to 4 spinor
#pragma acc data present(vp[0:size], h2p[0:size2]) \
  copyin(Ns, Nst_pad)
    {
#pragma acc parallel \
            num_workers(NUM_WORKERS) vector_length(VECTOR_LENGTH)
      {
#pragma acc loop
        for (int idx = 0; idx < Ns * Nst_pad; ++idx) {
          int idx2_wp = idx / NWP;
          int idx_in  = idx % NWP;
          int idx_out = idx2_wp / Ns;
          int is      = idx2_wp % Ns;
          int site = idx_in + NWP * idx_out;

          if(site < Nst){
            for(int ic=0; ic<NC; ++ic) {
              real_t vt1r, vt1i, vt2r, vt2i;
              vt1r = h2p[IDX2(NVCD2*Ns, (2*ic+0 + NVC*(0 + 2*is)), site)];
              vt1i = h2p[IDX2(NVCD2*Ns, (2*ic+1 + NVC*(0 + 2*is)), site)];
              vt2r = h2p[IDX2(NVCD2*Ns, (2*ic+0 + NVC*(1 + 2*is)), site)];
              vt2i = h2p[IDX2(NVCD2*Ns, (2*ic+1 + NVC*(1 + 2*is)), site)];

              vp[IDX2_SP_5D_R(ic, 0, is, Ns, site)] =  vt1r; // sp1
              vp[IDX2_SP_5D_I(ic, 0, is, Ns, site)] =  vt1i;
              vp[IDX2_SP_5D_R(ic, 1, is, Ns, site)] =  vt2r; // sp2
              vp[IDX2_SP_5D_I(ic, 1, is, Ns, site)] =  vt2i;
              vp[IDX2_SP_5D_R(ic, 2, is, Ns, site)] =  vt1i; // -i x sp1
              vp[IDX2_SP_5D_I(ic, 2, is, Ns, site)] = -vt1r;
              vp[IDX2_SP_5D_R(ic, 3, is, Ns, site)] = -vt2i; //  i x sp2
              vp[IDX2_SP_5D_I(ic, 3, is, Ns, site)] =  vt2r;
            } // ic
          } // site < NST
        } // idx
      }
    }
  } // ith==kernel
  Accel_Spinor::tensorProd_5din(force, mu, zeta, m_v2);


  ///////////////////////////////////////////////
  mu=3;   // tp
  ///////////////////////////////////////////////
  if (ith == ith_kernel){ // make 2 spinor, Dirac rep.
#pragma acc data present(yp[0:size], h1p[0:size2]) \
                 copyin(Ns, Nst_pad)
    {
#pragma acc parallel \
            num_workers(NUM_WORKERS) vector_length(VECTOR_LENGTH)
      {
#pragma acc loop
        for (int idx = 0; idx < Ns * NVC * Nst_pad; ++idx) {

          int idx2_wp = idx / NWP;
          int idx_in  = idx % NWP;
          int idx_out = idx2_wp / Ns;
          int is      = idx2_wp % Ns;
          int site = idx_in + NWP * idx_out;

          if(site < Nst){
            for(int ic=0; ic<NC; ++ic) {
              real_t vt1r, vt1i, vt2r, vt2i;
              vt1r = 2.0 * yp[IDX2_SP_5D_R(ic, 2, is, Ns, site)];
              vt1i = 2.0 * yp[IDX2_SP_5D_I(ic, 2, is, Ns, site)];
              vt2r = 2.0 * yp[IDX2_SP_5D_R(ic, 3, is, Ns, site)];
              vt2i = 2.0 * yp[IDX2_SP_5D_I(ic, 3, is, Ns, site)];
              h1p[IDX2(NVCD2*Ns, (2*ic+0 + NVC*(0 + 2*is)), site)] = vt1r;
              h1p[IDX2(NVCD2*Ns, (2*ic+1 + NVC*(0 + 2*is)), site)] = vt1i;
              h1p[IDX2(NVCD2*Ns, (2*ic+0 + NVC*(1 + 2*is)), site)] = vt2r;
              h1p[IDX2(NVCD2*Ns, (2*ic+1 + NVC*(1 + 2*is)), site)] = vt2i;
            } // ic
          } // site < NST
        } // idx
      }
    }
  } // ith==kernel

  m_shift->backward(m_h2, m_h1, mu);

  if (ith == ith_kernel){ // accumlate to 4 spinor
#pragma acc data present(vp[0:size], h2p[0:size2]) \
  copyin(Ns, Nst_pad)
    {
#pragma acc parallel \
            num_workers(NUM_WORKERS) vector_length(VECTOR_LENGTH)
      {
#pragma acc loop
        for (int idx = 0; idx < Ns * Nst_pad; ++idx) {
          int idx2_wp = idx / NWP;
          int idx_in  = idx % NWP;
          int idx_out = idx2_wp / Ns;
          int is      = idx2_wp % Ns;
          int site = idx_in + NWP * idx_out;

          if(site < Nst){
            for(int ic=0; ic<NC; ++ic) {
              real_t vt1r, vt1i, vt2r, vt2i;
              vt1r = h2p[IDX2(NVCD2*Ns, (2*ic+0 + NVC*(0 + 2*is)), site)];
              vt1i = h2p[IDX2(NVCD2*Ns, (2*ic+1 + NVC*(0 + 2*is)), site)];
              vt2r = h2p[IDX2(NVCD2*Ns, (2*ic+0 + NVC*(1 + 2*is)), site)];
              vt2i = h2p[IDX2(NVCD2*Ns, (2*ic+1 + NVC*(1 + 2*is)), site)];

              vp[IDX2_SP_5D_R(ic, 0, is, Ns, site)] =  0.0;  // 0
              vp[IDX2_SP_5D_I(ic, 0, is, Ns, site)] =  0.0;
              vp[IDX2_SP_5D_R(ic, 1, is, Ns, site)] =  0.0;  // 0
              vp[IDX2_SP_5D_I(ic, 1, is, Ns, site)] =  0.0;

              vp[IDX2_SP_5D_R(ic, 2, is, Ns, site)] =  vt1r; // sp1
              vp[IDX2_SP_5D_I(ic, 2, is, Ns, site)] =  vt1i;
              vp[IDX2_SP_5D_R(ic, 3, is, Ns, site)] =  vt2r; // sp2
              vp[IDX2_SP_5D_I(ic, 3, is, Ns, site)] =  vt2i;
            } // ic
          } // site < NST
        } // idx
      }
    }
  } // ith==kernel
  Accel_Spinor::tensorProd_5din(force, mu, zeta, m_v2);

#pragma omp barrier

}

//====================================================================
template<typename AFIELD>
void AForce_F_Domainwall_5din<AFIELD>::force_udiv1_Ddag(AFIELD& force,
                                                        const AFIELD& zeta,
                                                        const AFIELD& eta)
{
#pragma omp barrier


  int ith = ThreadManager::get_thread_id();
  int nth = ThreadManager::get_num_threads();
  int ith_kernel = 0;
  //if(nth > 1) ith_kernel = 1;

  real_t *wp = const_cast<AFIELD *>(&zeta)->ptr(0);
  real_t *etap = const_cast<AFIELD *>(&eta)->ptr(0);

  real_t *yp = m_v1.ptr(0);    // working vector: 4-spinor
  real_t *vp = m_v2.ptr(0);    // working vector: 4-spinor
  real_t *h1p = m_h1.ptr(0);   // working vector: 2-spinor
  real_t *h2p = m_h2.ptr(0);   // working vector: 2-spinor

  int Nst=force.nvol();
  int Nst_pad = CEIL_NWP(Nst);
  int Ns = m_Ns;
  int Nin5 = NVCD * Ns;
  int size = Nin5 * Nst_pad;
  int size2 = size/2;  // for 2-spinor
  real_t *b = &m_b[0];
  real_t *c = &m_c[0];
  real_t mq = m_mq;
  //  real_t M0 = m_M0;
  int Nx = CommonParameters::Nx();
  int Ny = CommonParameters::Ny();
  int Nz = CommonParameters::Nz();
  int Nt = CommonParameters::Nt();

  if (ith == ith_kernel){

#pragma acc data present(yp[0:size], wp[0:size])  \
  copyin(Nst, Nst_pad, Ns, Nin5, mq, b[0:Ns], c[0:Ns])
    {
      #pragma acc parallel num_workers(NUM_WORKERS) vector_length(VECTOR_LENGTH)
      {
        #pragma acc loop gang worker vector
        for (int idx = 0; idx < NVC * Nst_pad; ++idx) {
          // idx = idx_in + NWP * (ivc + NVC*idx_out)
          int idx2_wp = idx / NWP;
          int idx_in  = idx % NWP;
          int idx_out = idx2_wp / NVC;
          int site = idx_in + NWP * idx_out;
          int ivc  = idx2_wp % NVC;  // color and r/i

          if(site < Nst){
            for (int is = 0; is < Ns; ++is) {

              real_t wt1, wt2, wt3, wt4;
              real_t vt1, vt2, vt3, vt4;
              real_t xt1, xt2, xt3, xt4;

              int is_up = (is+1) % Ns;
              real_t Fup   = 0.5;
              if (is == Ns-1) Fup = -0.5 * mq;

              int ivcs_up = ivc + NVCD * is_up;
              wt1 = wp[IDX2(Nin5, ID1 + ivcs_up, site)];
              wt2 = wp[IDX2(Nin5, ID2 + ivcs_up, site)];
              wt3 = wp[IDX2(Nin5, ID3 + ivcs_up, site)];
              wt4 = wp[IDX2(Nin5, ID4 + ivcs_up, site)];

              vt1 = Fup * (wt1 - wt3);
              vt2 = Fup * (wt2 - wt4);
              vt3 = Fup * (wt3 - wt1);
              vt4 = Fup * (wt4 - wt2);

              int is_dn = (is-1 + Ns) % Ns;
              real_t Fdn   = 0.5;
              if (is == 0) { Fdn = -0.5 * mq; }

              int ivcs_dn = ivc + NVCD * is_dn;
              wt1 = wp[IDX2(Nin5, ID1 + ivcs_dn, site)];
              wt2 = wp[IDX2(Nin5, ID2 + ivcs_dn, site)];
              wt3 = wp[IDX2(Nin5, ID3 + ivcs_dn, site)];
              wt4 = wp[IDX2(Nin5, ID4 + ivcs_dn, site)];

              vt1 += Fdn * (wt1 + wt3);
              vt2 += Fdn * (wt2 + wt4);
              vt3 += Fdn * (wt3 + wt1);
              vt4 += Fdn * (wt4 + wt2);

              int ivcs = ivc + NVCD * is;
              wt1 = wp[IDX2(Nin5, ID1 + ivcs, site)];
              wt2 = wp[IDX2(Nin5, ID2 + ivcs, site)];
              wt3 = wp[IDX2(Nin5, ID3 + ivcs, site)];
              wt4 = wp[IDX2(Nin5, ID4 + ivcs, site)];

              real_t B2 = -0.5 * b[is];
              real_t C2 = -0.5 * c[is];

              yp[IDX2(Nin5, ID1 + ivcs, site)] = B2 * wt1 + C2 * vt1;
              yp[IDX2(Nin5, ID2 + ivcs, site)] = B2 * wt2 + C2 * vt2;
              yp[IDX2(Nin5, ID3 + ivcs, site)] = B2 * wt3 + C2 * vt3;
              yp[IDX2(Nin5, ID4 + ivcs, site)] = B2 * wt4 + C2 * vt4;
            } // is
          } // ist<Nst
        } // idx

      } // acc parallel
    } // acc data
  } //ith==ith_kernel


#pragma omp barrier
  // indexing of fermion field
  // index for the field
  //   idx_in + NMP * (r/i + 2*ic + 6*id)
  //          + NWP * 24*is
  //          + NWP * 24*Ns*idx_out

  ///////////////////////////////////////////////
  int mu=0;   // xp
  ///////////////////////////////////////////////
  if (ith == ith_kernel){ //  make 2 spinor
#pragma acc data present(etap[0:size], h1p[0:size2]) \
                 copyin(Ns, Nst_pad)
    {
#pragma acc parallel \
            num_workers(NUM_WORKERS) vector_length(VECTOR_LENGTH)
      {
#pragma acc loop
        for (int idx = 0; idx < Ns * NVC * Nst_pad; ++idx) {

          int idx2_wp = idx / NWP;
          int idx_in  = idx % NWP;
          int idx_out = idx2_wp / Ns;
          int is      = idx2_wp % Ns;
          int site = idx_in + NWP * idx_out;

          if(site < Nst){
            for(int ic=0; ic<NC; ++ic) {
              real_t vt1r, vt1i, vt2r, vt2i;
              vt1r = etap[IDX2_SP_5D_R(ic, 0, is, Ns, site)] + etap[IDX2_SP_5D_I(ic, 3, is, Ns, site)];
              vt1i = etap[IDX2_SP_5D_I(ic, 0, is, Ns, site)] - etap[IDX2_SP_5D_R(ic, 3, is, Ns, site)];
              vt2r = etap[IDX2_SP_5D_R(ic, 1, is, Ns, site)] + etap[IDX2_SP_5D_I(ic, 2, is, Ns, site)];
              vt2i = etap[IDX2_SP_5D_I(ic, 1, is, Ns, site)] - etap[IDX2_SP_5D_R(ic, 2, is, Ns, site)];
              h1p[IDX2(NVCD2*Ns, (2*ic+0 + NVC*(0 + 2*is)), site)] = vt1r;
              h1p[IDX2(NVCD2*Ns, (2*ic+1 + NVC*(0 + 2*is)), site)] = vt1i;
              h1p[IDX2(NVCD2*Ns, (2*ic+0 + NVC*(1 + 2*is)), site)] = vt2r;
              h1p[IDX2(NVCD2*Ns, (2*ic+1 + NVC*(1 + 2*is)), site)] = vt2i;
            } // ic
          } // site < NST
        } // idx
      }
    }
  } // ith==kernel

  m_shift->backward(m_h2, m_h1, mu);

  if (ith == ith_kernel){ // accumlate to 4 spinor
#pragma acc data present(vp[0:size], h2p[0:size2]) \
  copyin(Ns, Nst_pad)
    {
#pragma acc parallel \
            num_workers(NUM_WORKERS) vector_length(VECTOR_LENGTH)
      {
#pragma acc loop
        for (int idx = 0; idx < Ns * Nst_pad; ++idx) {
          int idx2_wp = idx / NWP;
          int idx_in  = idx % NWP;
          int idx_out = idx2_wp / Ns;
          int is      = idx2_wp % Ns;
          int site = idx_in + NWP * idx_out;

          if(site < Nst){
            for(int ic=0; ic<NC; ++ic) {
              real_t vt1r, vt1i, vt2r, vt2i;
              vt1r = h2p[IDX2(NVCD2*Ns, (2*ic+0 + NVC*(0 + 2*is)), site)];
              vt1i = h2p[IDX2(NVCD2*Ns, (2*ic+1 + NVC*(0 + 2*is)), site)];
              vt2r = h2p[IDX2(NVCD2*Ns, (2*ic+0 + NVC*(1 + 2*is)), site)];
              vt2i = h2p[IDX2(NVCD2*Ns, (2*ic+1 + NVC*(1 + 2*is)), site)];

              vp[IDX2_SP_5D_R(ic, 0, is, Ns, site)] =  vt1r; // sp1
              vp[IDX2_SP_5D_I(ic, 0, is, Ns, site)] =  vt1i;
              vp[IDX2_SP_5D_R(ic, 1, is, Ns, site)] =  vt2r; // sp2
              vp[IDX2_SP_5D_I(ic, 1, is, Ns, site)] =  vt2i;
              vp[IDX2_SP_5D_R(ic, 2, is, Ns, site)] = -vt2i; // -i x sp2
              vp[IDX2_SP_5D_I(ic, 2, is, Ns, site)] =  vt2r;
              vp[IDX2_SP_5D_R(ic, 3, is, Ns, site)] = -vt1i; // -i x sp1
              vp[IDX2_SP_5D_I(ic, 3, is, Ns, site)] =  vt1r;
            } // ic
          } // site < NST
        } // idx
      }
    }
  } // ith==kernel

  Accel_Spinor::tensorProd_5din(force, mu, m_v1, m_v2);

  ///////////////////////////////////////////////
  mu=1; // yp
  ///////////////////////////////////////////////
  if (ith == ith_kernel){ // make 2 spinor
#pragma acc data present(etap[0:size], h1p[0:size2]) \
                 copyin(Ns, Nst_pad)
    {
#pragma acc parallel \
            num_workers(NUM_WORKERS) vector_length(VECTOR_LENGTH)
      {
#pragma acc loop
        for (int idx = 0; idx < Ns * NVC * Nst_pad; ++idx) {
          // idx = idx_in + NWP* ( is + Ns*idx_out)

          int idx2_wp = idx / NWP;
          int idx_in  = idx % NWP;
          int idx_out = idx2_wp / Ns;
          int is      = idx2_wp % Ns;
          int site = idx_in + NWP * idx_out;

          if(site < Nst){
            for(int ic=0; ic<NC; ++ic) {
              real_t vt1r, vt1i, vt2r, vt2i;
              vt1r = etap[IDX2_SP_5D_R(ic, 0, is, Ns, site)] - etap[IDX2_SP_5D_R(ic, 3, is, Ns, site)];
              vt1i = etap[IDX2_SP_5D_I(ic, 0, is, Ns, site)] - etap[IDX2_SP_5D_I(ic, 3, is, Ns, site)];
              vt2r = etap[IDX2_SP_5D_R(ic, 1, is, Ns, site)] + etap[IDX2_SP_5D_R(ic, 2, is, Ns, site)];
              vt2i = etap[IDX2_SP_5D_I(ic, 1, is, Ns, site)] + etap[IDX2_SP_5D_I(ic, 2, is, Ns, site)];
              h1p[IDX2(NVCD2*Ns, (2*ic+0 + NVC*(0 + 2*is)), site)] = vt1r;
              h1p[IDX2(NVCD2*Ns, (2*ic+1 + NVC*(0 + 2*is)), site)] = vt1i;
              h1p[IDX2(NVCD2*Ns, (2*ic+0 + NVC*(1 + 2*is)), site)] = vt2r;
              h1p[IDX2(NVCD2*Ns, (2*ic+1 + NVC*(1 + 2*is)), site)] = vt2i;
            } // ic
          } // site < NST
        } // idx
      }
    }
  } // ith==kernel

  m_shift->backward(m_h2, m_h1, mu);

  if (ith == ith_kernel){ // accumlate to 4 spinor
#pragma acc data present(vp[0:size], h2p[0:size2]) \
  copyin(Ns, Nst_pad)
    {
#pragma acc parallel \
            num_workers(NUM_WORKERS) vector_length(VECTOR_LENGTH)
      {
#pragma acc loop
        for (int idx = 0; idx < Ns * Nst_pad; ++idx) {
          int idx2_wp = idx / NWP;
          int idx_in  = idx % NWP;
          int idx_out = idx2_wp / Ns;
          int is      = idx2_wp % Ns;
          int site = idx_in + NWP * idx_out;

          if(site < Nst){
            for(int ic=0; ic<NC; ++ic) {
              real_t vt1r, vt1i, vt2r, vt2i;
              vt1r = h2p[IDX2(NVCD2*Ns, (2*ic+0 + NVC*(0 + 2*is)), site)];
              vt1i = h2p[IDX2(NVCD2*Ns, (2*ic+1 + NVC*(0 + 2*is)), site)];
              vt2r = h2p[IDX2(NVCD2*Ns, (2*ic+0 + NVC*(1 + 2*is)), site)];
              vt2i = h2p[IDX2(NVCD2*Ns, (2*ic+1 + NVC*(1 + 2*is)), site)];

              vp[IDX2_SP_5D_R(ic, 0, is, Ns, site)] =  vt1r; // sp1
              vp[IDX2_SP_5D_I(ic, 0, is, Ns, site)] =  vt1i;
              vp[IDX2_SP_5D_R(ic, 1, is, Ns, site)] =  vt2r; // sp2
              vp[IDX2_SP_5D_I(ic, 1, is, Ns, site)] =  vt2i;
              vp[IDX2_SP_5D_R(ic, 2, is, Ns, site)] =  vt2r; //  sp2
              vp[IDX2_SP_5D_I(ic, 2, is, Ns, site)] =  vt2i;
              vp[IDX2_SP_5D_R(ic, 3, is, Ns, site)] = -vt1r; // -sp1
              vp[IDX2_SP_5D_I(ic, 3, is, Ns, site)] = -vt1i;
            } // ic
          } // site < NST
        } // idx
      }
    }
  } // ith==kernel

  Accel_Spinor::tensorProd_5din(force, mu, m_v1, m_v2);

  ///////////////////////////////////////////////
  mu=2; // zp
  ///////////////////////////////////////////////
  if (ith == ith_kernel){ // make 2 spinor
#pragma acc data present(etap[0:size], h1p[0:size2]) \
                 copyin(Ns, Nst_pad)
    {
#pragma acc parallel \
            num_workers(NUM_WORKERS) vector_length(VECTOR_LENGTH)
      {
#pragma acc loop
        for (int idx = 0; idx < Ns * NVC * Nst_pad; ++idx) {
          // idx = idx_in + NWP* ( is + Ns*idx_out)

          int idx2_wp = idx / NWP;
          int idx_in  = idx % NWP;
          int idx_out = idx2_wp / Ns;
          int is      = idx2_wp % Ns;
          int site = idx_in + NWP * idx_out;

          if(site < Nst){
            for(int ic=0; ic<NC; ++ic) {
              real_t vt1r, vt1i, vt2r, vt2i;
              vt1r = etap[IDX2_SP_5D_R(ic, 0, is, Ns, site)] + etap[IDX2_SP_5D_I(ic, 2, is, Ns, site)];
              vt1i = etap[IDX2_SP_5D_I(ic, 0, is, Ns, site)] - etap[IDX2_SP_5D_R(ic, 2, is, Ns, site)];
              vt2r = etap[IDX2_SP_5D_R(ic, 1, is, Ns, site)] - etap[IDX2_SP_5D_I(ic, 3, is, Ns, site)];
              vt2i = etap[IDX2_SP_5D_I(ic, 1, is, Ns, site)] + etap[IDX2_SP_5D_R(ic, 3, is, Ns, site)];
              h1p[IDX2(NVCD2*Ns, (2*ic+0 + NVC*(0 + 2*is)), site)] = vt1r;
              h1p[IDX2(NVCD2*Ns, (2*ic+1 + NVC*(0 + 2*is)), site)] = vt1i;
              h1p[IDX2(NVCD2*Ns, (2*ic+0 + NVC*(1 + 2*is)), site)] = vt2r;
              h1p[IDX2(NVCD2*Ns, (2*ic+1 + NVC*(1 + 2*is)), site)] = vt2i;
            } // ic
          } // site < NST
        } // idx
      }
    }
  } // ith==kernel

  m_shift->backward(m_h2, m_h1, mu);

  if (ith == ith_kernel){ // accumlate to 4 spinor
#pragma acc data present(vp[0:size], h2p[0:size2]) \
  copyin(Ns, Nst_pad)
    {
#pragma acc parallel \
            num_workers(NUM_WORKERS) vector_length(VECTOR_LENGTH)
      {
#pragma acc loop
        for (int idx = 0; idx < Ns * Nst_pad; ++idx) {
          int idx2_wp = idx / NWP;
          int idx_in  = idx % NWP;
          int idx_out = idx2_wp / Ns;
          int is      = idx2_wp % Ns;
          int site = idx_in + NWP * idx_out;

          if(site < Nst){
            for(int ic=0; ic<NC; ++ic) {
              real_t vt1r, vt1i, vt2r, vt2i;
              vt1r = h2p[IDX2(NVCD2*Ns, (2*ic+0 + NVC*(0 + 2*is)), site)];
              vt1i = h2p[IDX2(NVCD2*Ns, (2*ic+1 + NVC*(0 + 2*is)), site)];
              vt2r = h2p[IDX2(NVCD2*Ns, (2*ic+0 + NVC*(1 + 2*is)), site)];
              vt2i = h2p[IDX2(NVCD2*Ns, (2*ic+1 + NVC*(1 + 2*is)), site)];

              vp[IDX2_SP_5D_R(ic, 0, is, Ns, site)] =  vt1r; // sp1
              vp[IDX2_SP_5D_I(ic, 0, is, Ns, site)] =  vt1i;
              vp[IDX2_SP_5D_R(ic, 1, is, Ns, site)] =  vt2r; // sp2
              vp[IDX2_SP_5D_I(ic, 1, is, Ns, site)] =  vt2i;
              vp[IDX2_SP_5D_R(ic, 2, is, Ns, site)] = -vt1i; //  i x sp1
              vp[IDX2_SP_5D_I(ic, 2, is, Ns, site)] =  vt1r;
              vp[IDX2_SP_5D_R(ic, 3, is, Ns, site)] =  vt2i; // -i x sp2
              vp[IDX2_SP_5D_I(ic, 3, is, Ns, site)] = -vt2r;
            } // ic
          } // site < NST
        } // idx
      }
    }
  } // ith==kernel

  Accel_Spinor::tensorProd_5din(force, mu, m_v1, m_v2);

  ///////////////////////////////////////////////
  mu=3; // tp, Driac rep.
  ///////////////////////////////////////////////
  if (ith == ith_kernel){ // make 2 spinor
#pragma acc data present(etap[0:size], h1p[0:size2]) \
                 copyin(Ns, Nst_pad)
    {
#pragma acc parallel \
            num_workers(NUM_WORKERS) vector_length(VECTOR_LENGTH)
      {
#pragma acc loop
        for (int idx = 0; idx < Ns * NVC * Nst_pad; ++idx) {
          // idx = idx_in + NWP* ( is + Ns*idx_out)

          int idx2_wp = idx / NWP;
          int idx_in  = idx % NWP;
          int idx_out = idx2_wp / Ns;
          int is      = idx2_wp % Ns;
          int site = idx_in + NWP * idx_out;

          if(site < Nst){
            for(int ic=0; ic<NC; ++ic) {
              real_t vt1r, vt1i, vt2r, vt2i;
              vt1r = 2.0 * etap[IDX2_SP_5D_R(ic, 0, is, Ns, site)];
              vt1i = 2.0 * etap[IDX2_SP_5D_I(ic, 0, is, Ns, site)];
              vt2r = 2.0 * etap[IDX2_SP_5D_R(ic, 1, is, Ns, site)];
              vt2i = 2.0 * etap[IDX2_SP_5D_I(ic, 1, is, Ns, site)];
              h1p[IDX2(NVCD2*Ns, (2*ic+0 + NVC*(0 + 2*is)), site)] = vt1r;
              h1p[IDX2(NVCD2*Ns, (2*ic+1 + NVC*(0 + 2*is)), site)] = vt1i;
              h1p[IDX2(NVCD2*Ns, (2*ic+0 + NVC*(1 + 2*is)), site)] = vt2r;
              h1p[IDX2(NVCD2*Ns, (2*ic+1 + NVC*(1 + 2*is)), site)] = vt2i;
            } // ic
          } // site < NST
        } // idx
      }
    }
  } // ith==kernel

  m_shift->backward(m_h2, m_h1, mu);

  if (ith == ith_kernel){ // accumlate to 4 spinor
#pragma acc data present(vp[0:size], h2p[0:size2]) \
  copyin(Ns, Nst_pad)
    {
#pragma acc parallel \
            num_workers(NUM_WORKERS) vector_length(VECTOR_LENGTH)
      {
#pragma acc loop
        for (int idx = 0; idx < Ns * Nst_pad; ++idx) {
          int idx2_wp = idx / NWP;
          int idx_in  = idx % NWP;
          int idx_out = idx2_wp / Ns;
          int is      = idx2_wp % Ns;
          int site = idx_in + NWP * idx_out;

          if(site < Nst){
            for(int ic=0; ic<NC; ++ic) {
              real_t vt1r, vt1i, vt2r, vt2i;
              vt1r = h2p[IDX2(NVCD2*Ns, (2*ic+0 + NVC*(0 + 2*is)), site)];
              vt1i = h2p[IDX2(NVCD2*Ns, (2*ic+1 + NVC*(0 + 2*is)), site)];
              vt2r = h2p[IDX2(NVCD2*Ns, (2*ic+0 + NVC*(1 + 2*is)), site)];
              vt2i = h2p[IDX2(NVCD2*Ns, (2*ic+1 + NVC*(1 + 2*is)), site)];

              vp[IDX2_SP_5D_R(ic, 0, is, Ns, site)] =  vt1r; // sp1
              vp[IDX2_SP_5D_I(ic, 0, is, Ns, site)] =  vt1i;
              vp[IDX2_SP_5D_R(ic, 1, is, Ns, site)] =  vt2r; // sp2
              vp[IDX2_SP_5D_I(ic, 1, is, Ns, site)] =  vt2i;
              vp[IDX2_SP_5D_R(ic, 2, is, Ns, site)] =  0.0;  // 0
              vp[IDX2_SP_5D_I(ic, 2, is, Ns, site)] =  0.0;
              vp[IDX2_SP_5D_R(ic, 3, is, Ns, site)] =  0.0;  // 0
              vp[IDX2_SP_5D_I(ic, 3, is, Ns, site)] =  0.0;
            } // ic
          } // site < NST
        } // idx
      }
    }
  } // ith==kernel

  Accel_Spinor::tensorProd_5din(force, mu, m_v1, m_v2);

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
