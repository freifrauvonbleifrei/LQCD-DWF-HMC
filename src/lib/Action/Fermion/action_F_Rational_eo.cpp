/*!
        @file    action_F_Rational_eo.cpp

        @brief

        @author  Issaku Kanamori (kanamori)
                 $LastChangedBy: kanamori $

        @date    $LastChangedDate:: 2025-12-03 19:35:35 #$

        @version $LastChangedRevision: 2668 $
*/

#include "action_F_Rational_eo.h"
#include "lib/Field/index_eo.h"
//#include "lib/Fopr/fopr_Wilson_eo.h"

const std::string Action_F_Rational_eo::class_name = "Action_F_Rational_eo";

//====================================================================
void Action_F_Rational_eo::set_parameters(const Parameters& params)
{
  std::string vlevel;
  if (!params.fetch_string("verbose_level", vlevel)) {
    m_vl = vout.set_verbose_level(vlevel);
  }

  m_params_rational_H   = params.get_Parameters("Fopr_Rational_calcH");
  m_fopr_H.reset(new Fopr_Rational(m_fopr_eo, m_params_rational_H));
  m_params_rational_psf = m_params_rational_H;

  int exponent_numerator = m_params_rational_H.get_int("exponent_numerator");
  int exponent_denominator = m_params_rational_H.get_int("exponent_denominator");
  m_params_rational_psf.set_int("exponent_numerator", -exponent_numerator);
  m_params_rational_psf.set_int("exponent_denominator", 2*exponent_denominator);
  m_fopr_psf.reset(new Fopr_Rational(m_fopr_eo, m_params_rational_psf));
}


//====================================================================
void Action_F_Rational_eo::get_parameters(Parameters& params) const
{
  params.set_string("verbose_level", vout.get_verbose_level(m_vl));
  params.set_Parameters("Fopr_Rational_calcH", m_params_rational_H);
  //  params.set_Parameters("Fopr_Rational_psf",   m_params_rational_psf);
}



//====================================================================
double Action_F_Rational_eo::langevin(RandomNumbers *rand)
{
  const int Nvol = CommonParameters::Nvol();
  const int Ndim = CommonParameters::Ndim();

  const int NinF     = m_fopr_psf->field_nin();
  const int NvolF    = m_fopr_psf->field_nvol();
  const int NexF     = m_fopr_psf->field_nex();
  const int size_psf = NinF * NvolF * NexF * CommonParameters::NPE();

  assert(2*NvolF == Nvol);
  m_psf.reset(NinF, NvolF, NexF);

  vout.general(m_vl, "  %s: %s\n", class_name.c_str(), m_label.c_str());

  Field xi(NinF, NvolF, NexF);
  Field xi_lexical(NinF, Nvol, NexF);
  rand->gauss_lex_global(xi_lexical);

  Index_eo index;
  index.convertField(xi, xi_lexical, 0);
  m_fopr_psf->set_config(m_U);
  m_fopr_psf->mult(m_psf, xi);

  //  const double xi2   = xi.norm();
  //  const double H_psf = xi2 * xi2;
  const double H_psf = xi.norm2();

  vout.general(m_vl, "    H_Fstandard  = %18.8f\n", H_psf);
  vout.general(m_vl, "    H_F/dof      = %18.8f\n", H_psf / size_psf);

  return H_psf;
}


//====================================================================
double Action_F_Rational_eo::calcH()
{
  const int Nvol = CommonParameters::Nvol();
  const int Ndim = CommonParameters::Ndim();

  const int NinF     = m_fopr_H->field_nin();
  const int NvolF    = m_fopr_H->field_nvol();
  const int NexF     = m_fopr_H->field_nex();
  const int size_psf = NinF * NvolF * NexF * CommonParameters::NPE();

  vout.general(m_vl, "  %s: %s\n", class_name.c_str(), m_label.c_str());

  Field  v1(NinF, NvolF, NexF);
  int    Nconv;
  double diff;
  m_fopr_H->set_config(m_U);
  m_fopr_H->mult(v1, m_psf);
  const double H_psf = dot(v1, m_psf);

  vout.general(m_vl, "    H_Fstandard  = %18.8f\n", H_psf);
  vout.general(m_vl, "    H_F/dof      = %18.8f\n", H_psf / size_psf);

  return H_psf;
}


//====================================================================
void Action_F_Rational_eo::force(Field& force)
{
  const int Nin  = m_U->nin();
  const int Nvol = m_U->nvol();
  const int Nex  = m_U->nex();

  assert(force.nin() == Nin);
  assert(force.nvol() == Nvol);
  assert(force.nex() == Nex);

  vout.general(m_vl, "  %s: %s\n", class_name.c_str(), m_label.c_str());

  // solvers are in the force
  m_fopr_force_MD->set_config(m_U);
  m_fopr_force_MD->force_core(force, m_psf);

  double Fave, Fmax, Fdev;
  force.stat(Fave, Fmax, Fdev);
  vout.general(m_vl,
               "    Fstandard_ave = %12.6f  Fstandard_max = %12.6f  Fstandard_dev = %12.6f\n",
               Fave, Fmax, Fdev);
}


//====================================================================
//============================================================END=====
