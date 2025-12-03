/*!
        @file    force_F_Wilson_eo.cpp

        @brief

        @author  UEDA, Satoru  (sueda)
                 $LastChangedBy: kanamori $

        @date    $LastChangedDate:: 2025-12-03 19:35:35 #$

        @version $LastChangedRevision: 2668 $
*/



#include "force_F_Wilson_eo.h"

const std::string Force_F_Wilson_eo::class_name = "Force_F_Wilson_eo";

//====================================================================
void Force_F_Wilson_eo::set_parameters(const Parameters& params)
{
  std::string vlevel;
  if (!params.fetch_string("verbose_level", vlevel)) {
    m_vl = vout.set_verbose_level(vlevel);
  }

  //- fetch and check input parameters
  double           kappa;
  std::vector<int> bc;

  int err = 0;
  err += params.fetch_double("hopping_parameter", kappa);
  err += params.fetch_int_vector("boundary_condition", bc);

  if (err) {
    vout.crucial(m_vl, "Error at %s: input parameter not found.\n", class_name.c_str());
    exit(EXIT_FAILURE);
  }


  set_parameters(kappa, bc);
}


//====================================================================
void Force_F_Wilson_eo::get_parameters(Parameters& params) const
{
  params.set_double("hopping_parameter", m_kappa);
  params.set_int_vector("boundary_condition", m_boundary);

  params.set_string("verbose_level", vout.get_verbose_level(m_vl));
}


//====================================================================
void Force_F_Wilson_eo::set_parameters(const double kappa, const std::vector<int> bc)
{
  const int Ndim = CommonParameters::Ndim();
  const int Nvol  = CommonParameters::Nvol();

  //- print input parameters
  vout.general(m_vl, "%s:\n", class_name.c_str());
  vout.general(m_vl, "  kappa = %12.8f\n", kappa);
  for (int mu = 0; mu < Ndim; ++mu) {
    vout.general(m_vl, "  boundary[%d] = %2d\n", mu, bc[mu]);
  }

  //- range check
  // NB. kappa == 0 is allowed.
  assert(bc.size() == Ndim);

  //- store values
  m_kappa = kappa;

  m_boundary.resize(Ndim);
  m_boundary = bc;

  //- post-process
  m_fopr_w->set_parameters(m_kappa, m_boundary);

}


//====================================================================
void Force_F_Wilson_eo::force_udiv(Field& force_, const Field& eta_e)
{
  const int Nvol  = CommonParameters::Nvol();
  const int Ndim  = CommonParameters::Ndim();
  const int Nvol2 = Nvol / 2;

  // working vectors
  Field_G force(Nvol, Ndim);
  Field_F zeta_e(Nvol2, 1);

  m_fopr_w->set_mode("H");   // D = (1-Mt_eo Mt_oe)
  m_fopr_w->mult(zeta_e, eta_e);

  force_udiv1(force_,  zeta_e, eta_e); // zeta_e div(eta_e)
  force_udiv1(force ,  eta_e, zeta_e); // eta_d div(zeta_e)
  axpy(force_, 1.0, force);


  }


//====================================================================
void Force_F_Wilson_eo::force_udiv1(Field& force_, const Field& zeta_e, const Field& eta_e)
{
  const int Nvol = CommonParameters::Nvol();
  const int Ndim = CommonParameters::Ndim();
  const int Nvol2 = Nvol / 2;

  Field_G force_eo(Nvol2, Ndim);
  Field_G U_eo(Nvol2, Ndim);

  Field_F eta_o(Nvol2, 1);
  Field_F gm5_zeta_eo(Nvol2, 1);
  Field_F work_eo(Nvol2,1);
  Field_F work_eo2(Nvol2,1);

  m_index.convertField(U_eo, *m_U, 0);
  m_fopr_w->Meo(eta_o, eta_e, 1);
  m_fopr_w->mult_gm5(gm5_zeta_eo, zeta_e);
  force_eo.set(0.0);

  for (int mu = 0; mu < Ndim; ++mu) {
    work_eo.set(0.0);
#pragma omp barrier
    m_fopr_w->mult_up(mu, work_eo, eta_o, "Deo");
    mult_Field_Gd(work_eo2, 0, U_eo, mu, work_eo, 0);
    tensorProd_Field_F(force_eo, mu, gm5_zeta_eo, work_eo2);
  }
  m_index.reverseField(force_, force_eo, 0);

  m_fopr_w->Meo(work_eo, zeta_e, 1);
  m_fopr_w->mult_gm5(gm5_zeta_eo, work_eo);
  m_index.convertField(U_eo, *m_U, 1);
  force_eo.set(0.0);
  for (int mu = 0; mu < Ndim; ++mu) {
    work_eo.set(0.0);
#pragma omp barrier
    m_fopr_w->mult_up(mu, work_eo, eta_e, "Doe");
    mult_Field_Gd(work_eo2, 0, U_eo, mu, work_eo, 0);
    tensorProd_Field_F(force_eo, mu, gm5_zeta_eo, work_eo2);
  }
  m_index.reverseField(force_, force_eo, 1);

  //  scal(force_, -m_kappa); // force *= -m_kappa;
  //  scal(force_, -1.0);
  scal(force_, m_kappa);

}


//====================================================================
void Force_F_Wilson_eo::force_udiv1_impl(Field_G& force, const Field_F& zeta, const Field_F& eta, const int ieo)
{
  vout.crucial("%s: force_udiv1_impl, no longer used\n", class_name.c_str());
  abort();
}

//====================================================================
void Force_F_Wilson_eo::init_fopr(const std::string &repr){
    m_repr = repr;
    m_fopr_w      = new Fopr_Wilson_eo(repr);
}

//====================================================================
void Force_F_Wilson_eo::init_fopr(const Parameters &params){
    m_fopr_w      = new Fopr_Wilson_eo(params);
}

//====================================================================
//============================================================END=====
