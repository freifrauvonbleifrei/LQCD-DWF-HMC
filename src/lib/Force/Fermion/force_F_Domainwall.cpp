/*!
        @file    force_F_Domainwall.cpp

        @brief

        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: kanamori $

        @date    $LastChangedDate:: 2025-12-03 19:35:35 #$

        @version $LastChangedRevision: 2668 $
*/

#include "force_F_Domainwall.h"

const std::string Force_F_Domainwall::class_name
                                       = "Force_F_Domainwall";

//====================================================================
void Force_F_Domainwall::init()
{
  m_fopr_w  = new Fopr_Wilson;
  m_fopr_dw = new Fopr_Domainwall(m_fopr_w);
  m_force_w = new Force_F_Wilson_Nf2;
  m_boundary.resize(CommonParameters::Ndim());
}

//====================================================================
void Force_F_Domainwall::tidyup()
{
  delete m_force_w;
  delete m_fopr_dw;
  delete m_fopr_w;
}

//====================================================================
void Force_F_Domainwall::set_parameters(const Parameters& params)
{
  const string str_vlevel = params.get_string("verbose_level");

  m_vl = vout.set_verbose_level(str_vlevel);

  //- fetch and check input parameters
  double mq, M0, b, c;
  int Ns;
  std::vector<int> bc;
  double           alpha;

  int err = 0;
  err += params.fetch_double("quark_mass", mq);
  err += params.fetch_double("domain_wall_height", M0);
  err += params.fetch_int("extent_of_5th_dimension", Ns);
  err += params.fetch_int_vector("boundary_condition", bc);
  err += params.fetch_double("coefficient_b", b);
  err += params.fetch_double("coefficient_c", c);

  if (err) {
    vout.crucial(m_vl, "Error at %s: input parameter not found.\n",
                 class_name.c_str());
    exit(EXIT_FAILURE);
  }

  std::vector<double> b_vec(Ns);
  std::vector<double> c_vec(Ns);
  for(int is = 0; is < Ns; ++is){
    b_vec[is] = b;
    c_vec[is] = c;
  }

  int err3 = params.fetch_double("parameter_alpha", alpha);
  if (err3) {
    vout.general(m_vl, "  parameter alpha is not provided: set to 1.0.\n");
    alpha = 1.0;
  }

  set_parameters(mq, M0, Ns, bc, b_vec, c_vec, alpha);

}

//====================================================================
void Force_F_Domainwall::get_parameters(Parameters& params) const
{
  params.set_double("quark_mass", m_mq);
  params.set_double("domain_wall_height", m_M0);
  params.set_int("extent_of_5th_dimension", m_Ns);
  params.set_int_vector("boundary_condition", m_boundary);
  params.set_double_vector("coefficient_b", m_b);
  params.set_double_vector("coefficient_c", m_c);
  params.set_double("parameter_alpha", m_alpha);
  params.set_string("gamma_matrix_type", m_repr);

  params.set_string("verbose_level", vout.get_verbose_level(m_vl));
}

//====================================================================
void Force_F_Domainwall::set_parameters(
                                        const double mq,
                                        const double M0,
                                        const int Ns,
                                        const std::vector<int> bc,
                                        const std::vector<double> b_vec,
                                        const std::vector<double> c_vec,
                                        const double alpha)
{
  int Ndim = CommonParameters::Ndim();

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
    m_alpha = alpha;

    m_boundary.resize(Ndim);
    for(int mu = 0; mu < Ndim; ++mu){
      m_boundary[mu] = bc[mu];
    }

    if (m_b.size() != m_Ns) {
      m_b.resize(m_Ns);
      m_c.resize(m_Ns);
    }
    for(int is = 0; is < m_Ns; ++is){
      m_b[is] = b_vec[is];
      m_c[is] = c_vec[is];
    }
  }
#pragma omp barrier

  //- print input parameters
  vout.general(m_vl, "%s:\n", class_name.c_str());
  vout.general(m_vl, "  mq   = %12.8f\n", m_mq);
  vout.general(m_vl, "  M0   = %12.8f\n", m_M0);
  vout.general(m_vl, "  Ns   = %4d\n", m_Ns);
  for (int mu = 0; mu < Ndim; ++mu) {
    vout.general(m_vl, "  boundary[%d] = %2d\n", mu, m_boundary[mu]);
  }
  vout.general(m_vl, "  coefficients:\n");
  for (int is = 0; is < m_Ns; ++is) {
    vout.general(m_vl, "  b[%2d] = %16.10f  c[%2d] = %16.10f\n",
                 is, m_b[is], is, m_c[is]);
  }
  vout.general(m_vl, "  alpha = %8.4f\n", m_alpha);

  //- Domain-wall operator
  m_fopr_dw->set_parameters(m_mq, m_M0, m_Ns, m_boundary,
                            m_b, m_c, alpha);

  //- Wilson opr/force
  const double kappa = 1.0 / (8.0 - 2.0 * m_M0);
  m_fopr_w->set_parameters(kappa, m_boundary);
  m_force_w->set_parameters(kappa, m_boundary);

}

//====================================================================
void Force_F_Domainwall::set_config(Field *U)
{
  m_U = (Field_G *)U;
  m_fopr_w->set_config(U);
  m_fopr_dw->set_config(U);
  m_force_w->set_config(U);
}

//====================================================================
void Force_F_Domainwall::set_mode(const std::string& mode)
{
  #pragma omp barrier

  int ith=ThreadManager::get_thread_id();
  if(ith == 0) {  m_mode = mode; }

  #pragma omp barrier
}

//====================================================================
void Force_F_Domainwall::force_core1(Field& force,
                                     const Field& zeta,
                                     const Field& eta)
{
  vout.crucial(m_vl, "Error at %s: unimplemented.\n",
               class_name.c_str());
  exit(EXIT_FAILURE);
}

//====================================================================
void Force_F_Domainwall::force_udiv(Field& force_, const Field& eta)
{
  const int Nvol = CommonParameters::Nvol();
  const int Ndim = CommonParameters::Ndim();

  Field_F zeta(Nvol, m_Ns);
  m_fopr_dw->set_mode("H");
  m_fopr_dw->mult(zeta, eta);

  Field_G force(Nvol, Ndim);
  set_mode("H");
  force_udiv1(force, zeta, eta);
  copy(force_, force); // force_ = force;

  set_mode("Hdag");
  force_udiv1(force, eta, zeta);
  axpy(force_, 1.0, force); // force_ += force;

}

//====================================================================
void Force_F_Domainwall::force_udiv1(Field& force_,
                                     const Field& zeta,
                                     const Field& eta)
{
  int Nc   = CommonParameters::Nc();
  int Nd   = CommonParameters::Nd();
  int Nvol = CommonParameters::Nvol();
  int Ndim = CommonParameters::Ndim();
  int NinF = 2 * Nc * Nd;
  int NinG = 2 * Nc * Nc;

  assert(eta.nin() == NinF);
  assert(eta.nvol() == Nvol);
  assert(eta.nex() == m_Ns);

  assert(zeta.nin() == NinF);
  assert(zeta.nvol() == Nvol);
  assert(zeta.nex() == m_Ns);

  Field force(NinG,Nvol,Ndim);

  if(m_mode=="H"){
    force_udiv1_H(force_, zeta, eta);
  }else if(m_mode=="Hdag"){
    force_udiv1_Hdag(force_, zeta, eta);
  }else{
    vout.crucial(m_vl, "Error at %s: irrelevant mult mode: %s.\n",
                 class_name.c_str(), m_mode.c_str());
    exit(EXIT_FAILURE);
  }

}

//====================================================================
void Force_F_Domainwall::force_udiv1_H(Field& force_,
                                       const Field& zeta,
                                       const Field& eta)
{
  int Nc   = CommonParameters::Nc();
  int Nd   = CommonParameters::Nd();
  int Nvol = CommonParameters::Nvol();
  int Ndim = CommonParameters::Ndim();
  int NinF = 2 * Nc * Nd;
  int NinG = 2 * Nc * Nc;

  Field eta2(NinF, Nvol, m_Ns);
  Field eta4(NinF, Nvol, 1), zeta4(NinF, Nvol, 1);

  Field force1(NinG, Nvol, Ndim);

  force_.set(0.0);

  force_Himpl(eta2, eta);

  for(int is = 0; is < m_Ns; ++is){

    copy(zeta4, 0, zeta, is);
    copy(eta4,  0, eta2, is);
    scal(eta4, 4.0-m_M0);

    // force += d_forcew->force_udiv1(zeta4,eta4);
    m_force_w->force_udiv1(force1, zeta4, eta4);
    axpy(force_, 1.0, force1);

  }

}

//====================================================================
void Force_F_Domainwall::force_udiv1_Hdag(Field& force_,
                                          const Field& zeta,
                                          const Field& eta)
{
  int Nc   = CommonParameters::Nc();
  int Nd   = CommonParameters::Nd();
  int Nvol = CommonParameters::Nvol();
  int Ndim = CommonParameters::Ndim();
  int NinF = 2 * Nc * Nd;
  int NinG = 2 * Nc * Nc;

  Field zeta2(NinF, Nvol, m_Ns);
  Field eta4(NinF, Nvol, 1), zeta4(NinF, Nvol, 1);

  Field force1(NinG, Nvol, Ndim);

  force_.set(0.0);

  force_Himpl(zeta2,zeta);

  for(int is = 0; is < m_Ns; ++is){

    copy(zeta4, 0, zeta2, is);
    copy(eta4,  0, eta, is);
    scal(eta4, 4.0-m_M0);

    m_force_w->force_udiv1(force1, zeta4, eta4);
    axpy(force_, 1.0, force1);

  }

}

//====================================================================
void Force_F_Domainwall::force_Himpl(Field& eta2, const Field& eta)
{

  int NinF = eta.nin();
  int Nvol = eta.nvol();

  Field w4(NinF,Nvol,1), w4_up(NinF,Nvol,1), w4_dn(NinF,Nvol,1);
  Field Lw4(NinF,Nvol,1);

  for(int is = 0; is < m_Ns; ++is){
    int isR = m_Ns - 1 - is;

    Lw4.set(0.0);
    int is_up = (is+1) % m_Ns;
    double Fup = -0.5;
    if(is == m_Ns-1) Fup = m_mq/2.0;
    copy(w4_up, 0, eta, is_up);
    m_fopr_w->mult_gm5(w4, w4_up);
    axpy(w4_up, -1.0, w4);

    axpy(Lw4, Fup, w4_up);

    int is_dn = (is-1+m_Ns) % m_Ns;
    double Fdn = -0.5;
    if(is == 0) Fdn = m_mq/2.0;
    copy(w4_dn, 0, eta, is_dn);
    m_fopr_w->mult_gm5(w4, w4_dn);
    axpy(w4_dn, 1.0, w4);
    axpy(Lw4, Fdn, w4_dn);

    copy(w4, 0, eta, is);
    scal(w4, m_b[is]);
    axpy(w4, -m_c[is], Lw4);

    copy(eta2, isR, w4, 0);

  }

}

//============================================================END=====
