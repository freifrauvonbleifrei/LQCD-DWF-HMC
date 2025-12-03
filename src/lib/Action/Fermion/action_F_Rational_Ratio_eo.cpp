/*!
        @file    action_F_Rational_Ratio_eo.cpp

        @brief

        @author  Issaku Kanamori (kanamori)
                 $LastChangedBy: kanamori $

        @date    $LastChangedDate:: 2025-12-03 19:35:35 #$

        @version $LastChangedRevision: 2668 $
*/

#include "action_F_Rational_Ratio_eo.h"
#include "lib/Field/index_eo.h"

const std::string Action_F_Rational_Ratio_eo::class_name = "Action_F_Rational_Ratio_eo";

//====================================================================
void Action_F_Rational_Ratio_eo::set_parameters(const Parameters& params)
{
  std::string vlevel;
  if (!params.fetch_string("verbose_level", vlevel)) {
    m_vl = vout.set_verbose_level(vlevel);
  }

  // operators for Hamiltonian
  m_params_rational_denom_H   = params.get_Parameters("Fopr_Rational_denominator_calcH");
  m_params_rational_numer_H   = params.get_Parameters("Fopr_Rational_numerator_calcH");
  vout.general(m_vl, "%s: building denominator operator for calcH.\n", class_name.c_str());
  m_fopr_denom_H.reset(new Fopr_Rational(m_fopr_denom, m_params_rational_denom_H));
  vout.general(m_vl, "%s: building numerator operator for calcH.\n", class_name.c_str());
  m_fopr_numer_H.reset(new Fopr_Rational(m_fopr_numer, m_params_rational_numer_H));

  // operators needed to generate pseudo fermion
  // the denominator
  m_params_rational_denom_psf = m_params_rational_denom_H;
  int denom_exponent_numerator = m_params_rational_denom_H.get_int("exponent_numerator");
  int denom_exponent_denominator = m_params_rational_denom_H.get_int("exponent_denominator");
  m_params_rational_denom_psf.set_int("exponent_numerator", -denom_exponent_numerator);
  m_params_rational_denom_psf.set_int("exponent_denominator", 2*denom_exponent_denominator);
  vout.general(m_vl, "%s: building denominator operator for generating pseudo fermion.\n", class_name.c_str());
  m_fopr_denom_psf.reset(new Fopr_Rational(m_fopr_denom, m_params_rational_denom_psf));

  // the numerator
  m_params_rational_numer_psf = m_params_rational_numer_H;
  int numer_exponent_numerator = m_params_rational_numer_H.get_int("exponent_numerator");
  int numer_exponent_denominator = m_params_rational_numer_H.get_int("exponent_denominator");
  m_params_rational_numer_psf.set_int("exponent_numerator", -numer_exponent_numerator);
  m_params_rational_numer_psf.set_int("exponent_denominator", numer_exponent_denominator);
  vout.general(m_vl, "%s: building numerator operator for generating pseudo fermion.\n", class_name.c_str());
  m_fopr_numer_psf.reset(new Fopr_Rational(m_fopr_numer, m_params_rational_numer_psf));

}


//====================================================================
void Action_F_Rational_Ratio_eo::get_parameters(Parameters& params) const
{
  params.set_string("verbose_level", vout.get_verbose_level(m_vl));

  Parameters params_rational_denom_H;
  Parameters params_rational_numer_H;
  m_fopr_denom_H->get_parameters(params_rational_denom_H);
  m_fopr_numer_H->get_parameters(params_rational_numer_H);
  params.set_Parameters("Fopr_Rational_denomminator_calcH", params_rational_denom_H);
  params.set_Parameters("Fopr_Rational_numerator_calcH", params_rational_numer_H);
}



//====================================================================
double Action_F_Rational_Ratio_eo::langevin(RandomNumbers *rand)
{
  const int Nvol = CommonParameters::Nvol();
  const int Ndim = CommonParameters::Ndim();

  const int NinF     = m_fopr_denom_psf->field_nin();
  const int NvolF    = m_fopr_denom_psf->field_nvol();
  const int NexF     = m_fopr_denom_psf->field_nex();
  const size_t size_psf = NinF * size_t(NvolF) * NexF * CommonParameters::NPE();

  assert(2*NvolF == Nvol);
  assert(m_fopr_numer_psf->field_nin() == NinF);
  assert(m_fopr_numer_psf->field_nvol() == NvolF);
  assert(m_fopr_numer_psf->field_nex() == NexF);

  m_psf.reset(NinF, NvolF, NexF);

  vout.general(m_vl, "  %s: %s\n", class_name.c_str(), m_label.c_str());

  Field xi(NinF, NvolF, NexF);
  Field xi_tmp(NinF, NvolF, NexF);
  Field xi_lexical(NinF, Nvol, NexF);

  // generate gaussian noise
  rand->gauss_lex_global(xi_lexical);
  Index_eo index;
  index.convertField(xi, xi_lexical, 0);

  // generate pseudo fermion:
  m_fopr_denom_psf->set_config(m_U);
  m_fopr_denom_psf->mult(xi_tmp, xi);    // (B^dag B)^{1/4}
  m_fopr_numer_psf->set_config(m_U);
  m_fopr_numer_psf->mult(m_psf, xi_tmp); // (A^dag A)^{-1/4}

  const double H_psf = xi.norm2();

  vout.general(m_vl, "    H_Fstandard  = %18.8f\n", H_psf);
  vout.general(m_vl, "    H_F/dof      = %18.8f\n", H_psf / size_psf);

  return H_psf;
}


//====================================================================
double Action_F_Rational_Ratio_eo::calcH()
{
  const int Nvol = CommonParameters::Nvol();
  const int Ndim = CommonParameters::Ndim();

  const int NinF     = m_fopr_denom_H->field_nin();
  const int NvolF    = m_fopr_denom_H->field_nvol();
  const int NexF     = m_fopr_denom_H->field_nex();
  const size_t size_psf = NinF * size_t(NvolF) * NexF * CommonParameters::NPE();

  vout.general(m_vl, "  %s: %s\n", class_name.c_str(), m_label.c_str());

  Field  v1(NinF, NvolF, NexF);
  Field  v2(NinF, NvolF, NexF);
  int    Nconv;
  double diff;
  m_fopr_numer_H->set_config(m_U);
  m_fopr_numer_H->mult(v1, m_psf); // (A^dag A)^{1/4}
  m_fopr_denom_H->set_config(m_U);
  m_fopr_denom_H->mult(v2, v1);    // (B^dag B)^{-1/2}
  const dcomplex H_psf_cmplx = dotc(v1, v2);
  const double H_psf = real(H_psf_cmplx);

  vout.general(m_vl, "    H_Fstandard  = %18.8f\n", H_psf);
  vout.general(m_vl, "    H_F/dof      = %18.8f\n", H_psf / size_psf);

  return H_psf;
}


//====================================================================
void Action_F_Rational_Ratio_eo::force(Field& force)
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
