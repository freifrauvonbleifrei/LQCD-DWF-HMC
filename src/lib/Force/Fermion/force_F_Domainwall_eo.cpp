/*!
        @file    force_F_Domainwall_eo.cpp

        @brief

        @author  Issaku Kanamori (kanamori)
                 $LastChangedBy: kanamori $

        @date    $LastChangedDate:: 2025-12-03 19:35:35 #$

        @version $LastChangedRevision: 2668 $
*/



#include "force_F_Domainwall_eo.h"
#include "lib/ResourceManager/threadManager.h"

const std::string Force_F_Domainwall_eo::class_name = "Force_F_Domainwall_eo";

//====================================================================
void Force_F_Domainwall_eo::set_parameters(const Parameters& params)
{
  std::string vlevel;
  if (!params.fetch_string("verbose_level", vlevel)) {
    m_vl = vout.set_verbose_level(vlevel);
  }

  //- fetch and check input parameters
  string           gmset_type;
  double           mq, M0;
  int              Ns;
  std::vector<int> bc;
  double           b, c;

  int err_optional = 0;
  err_optional += params.fetch_string("gamma_matrix_type", gmset_type);

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


  Parameters params_fopr(params);
  if (m_repr != ""){
    if(err_optional){
      params_fopr.set_string("gamma_matrix_type", m_repr);
    } else {
      vout.general("%s: Warning: gamma_matrix_type is overwritten\n", class_name.c_str());
    }
  }
  if(err_optional){ m_repr = gmset_type; }

  if(m_fopr){
    m_fopr->set_parameters(params_fopr);
  } else {
    m_fopr = new Fopr_Domainwall_eo(params_fopr);
    //m_fopr = new Fopr_Wilson_eo(params_fopr);
  }

}


//====================================================================
void Force_F_Domainwall_eo::get_parameters(Parameters& params) const
{
  if(m_fopr){
    m_fopr->get_parameters(params);
  }
  params.set_string("verbose_level", vout.get_verbose_level(m_vl));
}


/*
//====================================================================
void Force_F_Domainwall_eo::set_parameters(const real_t mq, const real_t M0,
                                           const int Ns,
                                           const std::vector<int> bc,
                                           const real_t b, const real_t c)
{
  const int Ndim = CommonParameters::Ndim();
  const int Nvol  = CommonParameters::Nvol();


  m_M0 = real_t(M0);
  m_mq = real_t(mq);
  m_Ns = Ns;

  m_boundary.resize(m_Ndim);
  assert(bc.size() == m_Ndim);
  for (int mu = 0; mu < m_Ndim; ++mu) {
    m_boundary[mu] = bc[mu];
  }

  if (m_b.size() != m_Ns) {
    m_b.resize(m_Ns);
    m_c.resize(m_Ns);
  }
  for (int is = 0; is < m_Ns; ++is) {
    m_b[is] = real_t(b);
    m_c[is] = real_t(c);
  }

  //- print input parameters
  vout.general(m_vl, "Parameters of %s:\n", class_name.c_str());
  vout.general(m_vl, "   mq   = %8.4f\n", m_mq);
  vout.general(m_vl, "   M0   = %8.4f\n", m_M0);
  vout.general(m_vl, "   Ns   = %4d\n", m_Ns);
  for (int mu = 0; mu < m_Ndim; ++mu) {
    vout.general(m_vl, "   boundary[%d] = %2d\n", mu, m_boundary[mu]);
  }
  //vout.general(m_vl, "  coefficients b = %16.10f  c = %16.10f\n",
  //                                                m_b[0], m_c[0]);
  vout.general(m_vl, "  coefficients:\n");
  for (int is = 0; is < m_Ns; ++is) {
    vout.general(m_vl, "  b[%2d] = %16.10f  c[%2d] = %16.10f\n",
                 is, m_b[is], is, m_c[is]);
  }


  //- post-process
  m_fopr->set_parameters(m_mq, m_M0, m_Ns, m_bc, m_b, m_c);

}
*/


//====================================================================
void Force_F_Domainwall_eo::set_mode(const std::string& mode)
{
  #pragma omp barrier

  int ith=ThreadManager::get_thread_id();
  if(ith == 0) {  m_mode = mode; }

  #pragma omp barrier
}

//====================================================================
void Force_F_Domainwall_eo::force_udiv(Field& force_, const Field& eta_e)
{
  const int Nvol  = CommonParameters::Nvol();
  const int Ndim  = CommonParameters::Ndim();
  const int Nvol2 = Nvol / 2;
  const int Nex   = eta_e.nex();

  // working vectors
  Field_G force(Nvol, Ndim);
  Field_F zeta_e(Nvol2, Nex);

  m_fopr->set_mode("H");                   // H = Gm5 (1-Mt_eo Mt_oe)
  m_fopr->mult(zeta_e, eta_e);
  force_udiv1_H(force_,   zeta_e, eta_e);  // <eta_e| Hdag H_div |eta_e>
  force_udiv1_Hdag(force , eta_e, zeta_e); // <eta_e| Hdag_div H |eta_e>
  axpy(force_, 1.0, force);

}


//====================================================================
void Force_F_Domainwall_eo::force_udiv1(Field& force_, const Field& zeta_e, const Field& eta_e)
{
  if(m_mode == "H"){
    // assumes eta_e=|eta>,  zeta_e = H|eta>
    //  calculates force_= <zeta_e| div(H) |eta_e>
    force_udiv1_H(force_, zeta_e, eta_e);
  } else if (m_mode == "Hdag") {
    // assumes zeta_e=|eta>,  eta_e = H|eta>
    //  calculates force_= <zeta_e| div(Hdag) |eta_e>
    force_udiv1_Hdag(force_, zeta_e, eta_e);
  } else {
    vout.crucial(m_vl, "%s: unknown mode %s.\n",
                 class_name.c_str(), m_mode.c_str());
  }

}

//====================================================================
void Force_F_Domainwall_eo::force_udiv1_H(Field& force_, const Field& zeta_e, const Field& eta_e)
{

  const int Nvol = CommonParameters::Nvol();
  const int Ndim = CommonParameters::Ndim();
  const int Nvol2 = Nvol / 2;
  const int Nex = eta_e.nex();

  Field_G force_eo(Nvol2, Ndim);
  Field_G U_eo(Nvol2, Ndim);

  //  Field_F gm5_zeta_eo(Nvol2, 1);
  Field_F work_eo(Nvol2,Nex);
  Field_F work_eo2(Nvol2,Nex);

  const Field_F &u_e = eta_e;
  Field_F u_o(Nvol2, Nex);
  Field_F v_e(Nvol2, Nex);
  Field_F v_o(Nvol2, Nex);


  m_fopr->mult(work_eo, u_e, "Doe");
  m_fopr->mult(u_o, work_eo, "Doo_inv");

  m_fopr->mult_gm5R(work_eo, zeta_e);       // work_eo = A eta
  m_fopr->mult_dag(v_e, work_eo, "Dee_inv");
  m_fopr->mult_dag(work_eo, v_e, "Deo");  // (Deo)^dag: odd <-- even
  m_fopr->mult_dag(v_o, work_eo, "Doo_inv");

  m_index.convertField(U_eo, *m_U, 0);
  force_eo.set(0.0);
  for (int mu = 0; mu < Ndim; ++mu) {
    work_eo.set(0.0);
#pragma omp barrier
    m_fopr->mult_up(mu, work_eo, u_o, "Deo");
    for(int s = 0; s < Nex; ++s) {
      mult_Field_Gd(work_eo2, s, U_eo, mu, work_eo, s);
    }
    tensorProd_Field_F(force_eo, mu, v_e, work_eo2);
  }
  m_index.reverseField(force_, force_eo, 0);


  m_index.convertField(U_eo, *m_U, 1);
  force_eo.set(0.0);
  for (int mu = 0; mu < Ndim; ++mu) {
    work_eo.set(0.0);
#pragma omp barrier
    m_fopr->mult_up(mu, work_eo, u_e, "Doe");
    for(int s = 0; s < Nex; ++s) {
      mult_Field_Gd(work_eo2, s, U_eo, mu, work_eo, s);
    }
    tensorProd_Field_F(force_eo, mu, v_o, work_eo2);
  }
  m_index.reverseField(force_, force_eo, 1);
  scal(force_, -1.0);

}


//====================================================================
void Force_F_Domainwall_eo::force_udiv1_Hdag(Field& force_, const Field& eta_e, const Field& zeta_e)
{
  // eta_e  =  |eta>
  // zeta_e = H|eta>
  const int Nvol = CommonParameters::Nvol();
  const int Ndim = CommonParameters::Ndim();
  const int Nvol2 = Nvol / 2;
  const int Nex = eta_e.nex();

  Field_G force_eo(Nvol2, Ndim);
  Field_G U_eo(Nvol2, Ndim);

  Field_F work_eo(Nvol2,Nex);
  Field_F work_eo2(Nvol2,Nex);

  Field_F u_e(Nvol2, Nex);
  Field_F u_o(Nvol2, Nex);
  const Field_F &v_e = eta_e;
  Field_F v_o(Nvol2, Nex);

  force_.set(0.0);

  // u_e = zeta_e
  // u_o = Doo_inv Doe zeta_e
  m_fopr->mult(work_eo, v_e, "Doe");
  m_fopr->mult(v_o, work_eo, "Doo_inv");

  // v_e = (Gm5 Dee_inv)^dag eta_e
  // v_o = (Gm5 Dee_inv Deo Doo_inv)^dag eta_e
  //     = (Deo Doo_inv)^dag v_e
  m_fopr->mult_gm5R(work_eo, zeta_e);
  m_fopr->mult_dag(u_e, work_eo, "Dee_inv");
  m_fopr->mult_dag(work_eo, u_e, "Deo");      // (Deo)^dag: odd <-- even
  m_fopr->mult_dag(u_o, work_eo, "Doo_inv");

  m_index.convertField(U_eo, *m_U, 0);
  force_eo.set(0.0);
  for (int mu = 0; mu < Ndim; ++mu) {
    work_eo.set(0.0);
#pragma omp barrier
    m_fopr->mult_up(mu, work_eo, u_o, "Ddag_oe");
    for(int s = 0; s < Nex; ++s) {
      mult_Field_Gd(work_eo2, s, U_eo, mu, work_eo, s);
    }
    tensorProd_Field_F(force_eo, mu, v_e, work_eo2);
  }
  m_index.reverseField(force_, force_eo, 0);

  m_index.convertField(U_eo, *m_U, 1);
  force_eo.set(0.0);

  for (int mu = 0; mu < Ndim; ++mu) {
    work_eo.set(0.0);
    double tmp0=work_eo.norm2();
#pragma omp barrier
    m_fopr->mult_up(mu, work_eo, u_e, "Ddag_eo");
    for(int s = 0; s < Nex; ++s) {
      mult_Field_Gd(work_eo2, s, U_eo, mu, work_eo, s);
    }
    tensorProd_Field_F(force_eo, mu, v_o, work_eo2);
  }

  m_index.reverseField(force_, force_eo, 1);
  scal(force_, -1.0);

}


//====================================================================
void Force_F_Domainwall_eo::force_udiv1_impl(Field_G& force, const Field_F& zeta, const Field_F& eta, const int ieo)
{
  vout.crucial("%s: force_udiv1_impl, no longer used\n", class_name.c_str());
  abort();
}


//====================================================================
//============================================================END=====
