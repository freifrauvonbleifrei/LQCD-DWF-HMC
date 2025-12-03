/*!
        @file    aforce_F_Domainwall_5din_eo-tmpl.h

        @brief

        @author  Issaku Kanamori (kanamori)
                 $LastChangedBy: kanamori $

        @date    $LastChangedDate: 2025-12-03 19:35:35 +0900 (2025年12月03日 (水)) $

        @version $LastChangedRevision: 2668 $
*/

// required by aforce_F_Domainwall_5din
#include "lib_alt_Accel/Fopr/afopr_Domainwall_5din.h"
#include "lib_alt_Accel/Force/Fermion/aforce_F_Domainwall_5din.h"
#include "lib_alt_Accel/Field/shiftAField_lex.h"

#include "lib_alt_Accel/Force/Fermion/aforce_F_Domainwall_5din_eo.h"
#include "lib_alt_Accel/Field/aindex_lex.h"
#include "lib_alt_Accel/Field/aindex_eo.h"
#include "lib_alt_Accel/Field/aindex_eo-inc.h"

template<typename AFIELD>
const std::string AForce_F_Domainwall_5din_eo<AFIELD>::class_name
                               = "AForce_F_Domainwall_5din_eo<AFIELD>";
//====================================================================
template<typename AFIELD>
void AForce_F_Domainwall_5din_eo<AFIELD>::init(const Parameters& params)
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
  int Nvol2 = Nvol/2;

  m_Nvcd = 2 * Nc * Nd;
  m_boundary.resize(Ndim);

  set_parameters(params);

  int err = 0;

  Parameters params_kernel = params;

  m_v1.reset(m_NinF, Nvol2, 1);
  m_v2.reset(m_NinF, Nvol2, 1);
  m_v3.reset(m_NinF, Nvol2, 1);
  m_z1.reset(m_NinF, Nvol2, 1);
  m_eta_lex.reset(m_NinF, Nvol, 1);
  m_zeta_lex.reset(m_NinF, Nvol, 1);

  m_force1.reset(NinG, Nvol, Ndim);

  vout.decrease_indent();

  vout.general(m_vl, "%s: construction finished.\n",
               class_name.c_str());

}

//====================================================================
template<typename AFIELD>
void AForce_F_Domainwall_5din_eo<AFIELD>::tidyup()
{
  //  delete m_force_w;
  //  delete m_fopr_dw;
  //  delete m_fopr_w;
}

//====================================================================
template<typename AFIELD>
void AForce_F_Domainwall_5din_eo<AFIELD>::set_parameters(const Parameters& params)
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

  if(ith == 0){
    if( !m_fopr_dw) {
      m_fopr_dw.reset(new AFopr_Domainwall_5din_eo<AFIELD>(params));
    } else {
      m_fopr_dw->set_parameters(params);
    }

    Parameters params_force_lex(params);
    params_force_lex.set_string("fermion_type", "Domainwall_5din");
    if( !m_force_lex) {
      m_force_lex.reset(new AForce_F_Domainwall_5din<AFIELD>(params_force_lex));
    } else {
      m_force_lex->set_parameters(params_force_lex);
    }
  }

#pragma omp barrier

}


//====================================================================
template<typename AFIELD>
void AForce_F_Domainwall_5din_eo<AFIELD>::set_parameters(
                                          const real_t mq,
                                          const real_t M0,
                                          const int Ns,
                                          const std::vector<int> bc,
                                          const real_t b,
                                          const real_t c)
{
  const int Ndim = CommonParameters::Ndim();
  const int Nvol = CommonParameters::Nvol();
  const int Nvol2 = Nvol/2;

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

    m_boundary.resize(Ndim);
    m_boundary = bc;

    if (m_b.size() != m_Ns) {
      m_b.resize(m_Ns);
      m_c.resize(m_Ns);
      m_v1.reset(m_NinF, Nvol2, 1);
      m_v2.reset(m_NinF, Nvol2, 1);
      m_v3.reset(m_NinF, Nvol2, 1);
      m_z1.reset(m_NinF, Nvol2, 1);
      m_eta_lex.reset(m_NinF, Nvol, 1);
      m_zeta_lex.reset(m_NinF, Nvol, 1);
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
void AForce_F_Domainwall_5din_eo<AFIELD>::get_parameters(Parameters& params) const
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
void AForce_F_Domainwall_5din_eo<AFIELD>::set_config(Field* U)
{
  int nth = ThreadManager::get_num_threads();
  if(nth > 1){
    set_config_impl(U);
  } else {
    set_config_omp(U);
  }
  m_fopr_dw->set_config(U);
  m_force_lex->set_config(U);

  m_U = U;

}

//====================================================================
template<typename AFIELD>
void AForce_F_Domainwall_5din_eo<AFIELD>::set_config_omp(Field* U)
{
#pragma omp parallel
 { set_config_impl(U); }
}

//====================================================================
template<typename AFIELD>
void AForce_F_Domainwall_5din_eo<AFIELD>::set_config_impl(Field* U)
{
#pragma omp barrier

  AIndex_lex<real_t, AFIELD::IMPL> index;
  convert_gauge(index, m_Ucp, *U);

#pragma omp barrier
}

//====================================================================
template<typename AFIELD>
void AForce_F_Domainwall_5din_eo<AFIELD>::set_mode(const std::string& mode)
{
#pragma omp barrier

  int ith = ThreadManager::get_thread_id();
  if (ith == 0) m_mode = mode;

#pragma omp barrier
}

//====================================================================
template<typename AFIELD>
void AForce_F_Domainwall_5din_eo<AFIELD>::force_udiv(AFIELD& force,
                                                     const AFIELD& eta)
{
#pragma omp barrier

  m_fopr_dw->set_mode("D");

  m_fopr_dw->set_config(m_U);
  m_fopr_dw->mult(m_z1, eta);

  force_udiv1_D(force, m_z1, eta);
  force_udiv1_Ddag(m_force1, eta, m_z1);
  axpy(force, real_t(1.0), m_force1);


#pragma omp barrier
}

//====================================================================
template<typename AFIELD>
void AForce_F_Domainwall_5din_eo<AFIELD>::force_udiv1(AFIELD& force,
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
void AForce_F_Domainwall_5din_eo<AFIELD>::force_udiv1_D(AFIELD& force,
                                                        const AFIELD& zeta,
                                                        const AFIELD& eta)
{
#pragma omp barrier

  m_fopr_dw->mult(m_v1, eta,  "Doe");
  m_fopr_dw->mult(m_v2, m_v1, "Doo_inv");
  m_index_eo.merge(m_eta_lex, eta, m_v2);

  m_fopr_dw->mult_dag(m_v3, zeta, "Dee_inv");
  m_fopr_dw->mult_dag(m_v1, m_v3, "Deo");
  m_fopr_dw->mult_dag(m_v2, m_v1, "Doo_inv");
  m_index_eo.merge(m_zeta_lex, m_v3, m_v2);

  m_force_lex->set_mode("D");
  m_force_lex->force_udiv1(force, m_zeta_lex, m_eta_lex);

  scal(force, real_t(-1.0));
#pragma omp barrier

}


template<typename AFIELD>
void AForce_F_Domainwall_5din_eo<AFIELD>::force_udiv1_Ddag(AFIELD& force,
                                                           const AFIELD& zeta,
                                                           const AFIELD& eta)
{
#pragma omp barrier
  m_fopr_dw->mult_dag(m_v3, eta, "Dee_inv");
  m_fopr_dw->mult_dag(m_v1,  m_v3, "Deo");
  m_fopr_dw->mult_dag(m_v2,  m_v1, "Doo_inv");
  m_index_eo.merge(m_eta_lex, m_v3, m_v2);

  m_fopr_dw->mult(m_v2, zeta, "Doe");
  m_fopr_dw->mult(m_v1, m_v2, "Doo_inv");
  m_index_eo.merge(m_zeta_lex, zeta, m_v1);

  m_force_lex->set_mode("Ddag");
  m_force_lex->force_udiv1(force, m_zeta_lex, m_eta_lex);

  scal(force, real_t(-1.0));

#pragma omp barrier

}


//====================================================================
template<typename AFIELD>
void AForce_F_Domainwall_5din_eo<AFIELD>::force_udiv1_H(AFIELD& force,
                                                        const AFIELD& zeta,
                                                        const AFIELD& eta)
{
#pragma omp barrier
  //  force.set(0.0);

  m_fopr_dw->mult(m_v1, eta,  "Doe");
  m_fopr_dw->mult(m_v2, m_v1, "Doo_inv");
  m_index_eo.merge(m_eta_lex, eta, m_v2);

  m_fopr_dw->mult_gm5R(m_v1, zeta);
  m_fopr_dw->mult_dag(m_v2, m_v1, "Dee_inv");
  m_fopr_dw->mult_gm5R(m_v3, m_v2);
  m_fopr_dw->mult_dag(m_v1, m_v2, "Deo");
  m_fopr_dw->mult_dag(m_v2, m_v1, "Doo_inv");
  m_fopr_dw->mult_gm5R(m_v1, m_v2);
  m_index_eo.merge(m_zeta_lex, m_v3, m_v1);

  m_force_lex->set_mode("H");
  m_force_lex->force_udiv1(force, m_zeta_lex, m_eta_lex);

  scal(force, real_t(-1.0));

#pragma omp barrier
}



//====================================================================
template<typename AFIELD>
void AForce_F_Domainwall_5din_eo<AFIELD>::force_udiv1_Hdag(AFIELD& force,
                                                           const AFIELD& zeta,
                                                           const AFIELD& eta)
{
#pragma omp barrier

  //  force.set(0.0);

  m_fopr_dw->mult_gm5R(m_v1, eta);
  m_fopr_dw->mult_dag(m_v2, m_v1, "Dee_inv");
  m_fopr_dw->mult_gm5R(m_v3, m_v2);
  m_fopr_dw->mult_dag(m_v1,  m_v2, "Deo");
  m_fopr_dw->mult_dag(m_v2,  m_v1, "Doo_inv");
  m_fopr_dw->mult_gm5R(m_v1, m_v2);


  m_index_eo.merge(m_eta_lex, m_v3, m_v1);

  m_fopr_dw->mult(m_v2, zeta, "Doe");
  m_fopr_dw->mult(m_v1, m_v2, "Doo_inv");
  m_index_eo.merge(m_zeta_lex, zeta, m_v1);

  m_force_lex->set_mode("Hdag");
  m_force_lex->force_udiv1(force, m_zeta_lex, m_eta_lex);

  scal(force, real_t(-1.0));

#pragma omp barrier
}
//============================================================END=====
