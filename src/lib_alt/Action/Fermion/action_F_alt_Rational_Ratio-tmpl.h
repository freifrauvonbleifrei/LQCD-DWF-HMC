/*!
      @file    action_F_alt_Rational_Ratio-tmpl.h
      @brief
      @author  Hideo Matsufuru (matufuru)
               $LastChangedBy: kanamori $
      @date    $LastChangedDate:: 2025-12-03 19:35:35 #$
      @version $LastChangedRevision: 2668 $
*/

#include "lib/Tools/timer.h"
#include "lib/Field/index_eo.h"

template<typename AFIELD>
const std::string Action_F_alt_Rational_Ratio<AFIELD>::class_name
                                     = "Action_F_alt_Rational_Ratio";

//====================================================================
template<typename AFIELD>
void Action_F_alt_Rational_Ratio<AFIELD>::init()
{
  m_vl = CommonParameters::Vlevel();

  int Nin  = m_fopr1->field_nin();
  int Nvol = m_fopr1->field_nvol();
  int Nex  = m_fopr1->field_nex();

  m_v1.reset(Nin, Nvol, Nex);
  m_v2.reset(Nin, Nvol, Nex);
}

//====================================================================
template<typename AFIELD>
void Action_F_alt_Rational_Ratio<AFIELD>::tidyup()
{
  // do nothing.
}

//====================================================================
template<typename AFIELD>
void Action_F_alt_Rational_Ratio<AFIELD>::set_parameters(
                                            const Parameters& params)
{
  const string str_vlevel = params.get_string("verbose_level");
  m_vl = vout.set_verbose_level(str_vlevel);

  typedef AFopr_Rational<AFIELD> Fopr_rational;

  // operators for Hamiltonian
  m_params_rational1_H   = params.get_Parameters("Fopr_Rational_numerator_calcH");
  m_params_rational2_H   = params.get_Parameters("Fopr_Rational_denominator_calcH");
  vout.general(m_vl, "%s: building numerator operator (fopr1) for calcH.\n", class_name.c_str());
  m_fopr1_rational_H.reset(new Fopr_rational(m_fopr1, m_params_rational1_H));
  vout.general(m_vl, "%s: building denominator operator (fopr2) for calcH.\n", class_name.c_str());
  m_fopr2_rational_H.reset(new Fopr_rational(m_fopr2, m_params_rational2_H));
  // operators needed to generate pseudo fermion

  // the numerator
  m_params_rational1_langevin = m_params_rational1_H;
  int exponent1_numerator = m_params_rational1_H.get_int("exponent_numerator");
  int exponent1_denominator = m_params_rational1_H.get_int("exponent_denominator");
  m_params_rational1_langevin.set_int("exponent_numerator", -exponent1_numerator);
  m_params_rational1_langevin.set_int("exponent_denominator", exponent1_denominator);
  vout.general(m_vl, "%s: building numerator operator for generating pseudo fermion.\n", class_name.c_str());
  m_fopr1_rational_langevin.reset(new Fopr_rational(m_fopr1, m_params_rational1_langevin));

  // the denominator
  m_params_rational2_langevin = m_params_rational2_H;
  int exponent2_numerator = m_params_rational2_H.get_int("exponent_numerator");
  int exponent2_denominator = m_params_rational2_H.get_int("exponent_denominator");
  m_params_rational2_langevin.set_int("exponent_numerator", -exponent2_numerator);
  m_params_rational2_langevin.set_int("exponent_denominator", 2*exponent2_denominator);
  vout.general(m_vl, "%s: building denominator operator for generating pseudo fermion.\n", class_name.c_str());
  m_fopr2_rational_langevin.reset(new Fopr_rational(m_fopr2, m_params_rational2_langevin));


}

//====================================================================
template<typename AFIELD>
void Action_F_alt_Rational_Ratio<AFIELD>::get_parameters(
                                            Parameters& params) const
{
  params.set_string("verbose_level", vout.get_verbose_level(m_vl));

  Parameters params_rational1_H;
  Parameters params_rational2_H;
  m_fopr1_rational_H->get_parameters(params_rational1_H);
  m_fopr2_rational_H->get_parameters(params_rational2_H);
  params.set_Parameters("Fopr_Rational_numerator_calcH",    params_rational1_H);
  params.set_Parameters("Fopr_Rational_denomminator_calcH", params_rational2_H);
}


//====================================================================
template<typename AFIELD>
void Action_F_alt_Rational_Ratio<AFIELD>::set_config(Field *U)
{
  m_U = U;

  m_fopr1->set_config(m_U);
  m_fopr2->set_config(m_U);
  m_force_MD->set_config(m_U);

}

//====================================================================
template<typename AFIELD>
double Action_F_alt_Rational_Ratio<AFIELD>::langevin(RandomNumbers *rand)
{
  const int NinF     = m_fopr1_rational_langevin->field_nin();
  const int NvolF    = m_fopr1_rational_langevin->field_nvol();
  const int NexF     = m_fopr1_rational_langevin->field_nex();
  const int size_psf = NinF * NvolF * NexF * CommonParameters::NPE();

  Timer timer;
  timer.start();

  m_psf.reset(NinF, NvolF, NexF);

  vout.general(m_vl, "  %s: %s\n",
               class_name.c_str(), m_label.c_str());

  const int Nvol = CommonParameters::Nvol();
  const int Nc=CommonParameters::Nc();
  const int Nd=CommonParameters::Nd();
  int Nin=NinF;
  int Nex=NexF;
  if(2*Nc*Nd != NinF){ // temporal hack;
    Nin = 2*Nc*Nd;
    Nex = NinF/Nin;
  }
  Field xi(Nin, NvolF, Nex);

  AFIELD xiA(NinF, NvolF, NexF);

  if(NvolF == Nvol){
    vout.detailed(m_vl, "  Gaussian filed with lexical action.\n");
    rand->gauss_lex_global(xi);
  }else if(NvolF == Nvol/2){
    vout.detailed(m_vl, "  Gaussian field with even-odd action.\n");
    Field xi_lex(Nin, Nvol, Nex);
    rand->gauss_lex_global(xi_lex);
    Index_eo idx_eo;
    idx_eo.convertField(xi, xi_lex, 0);
  }else{
     vout.crucial(m_vl, "  Unsupported field volume.\n");
     exit(EXIT_FAILURE);
  }

  double H_psf;

  #pragma omp parallel
  {
    if(m_fopr1_rational_langevin->needs_convert()){
      m_fopr1_rational_langevin->convert(xiA, xi);
    }else{
      AIndex_lex<real_t,AFIELD::IMPL> index_alt;
      convert(index_alt, xiA, xi);
    }

    m_fopr2_rational_langevin->set_config(m_U);
    m_fopr2_rational_langevin->mult(m_v1, xiA);     // (B^dag B)^{1/4}

    m_fopr1_rational_langevin->set_config(m_U);
    m_fopr1_rational_langevin->mult(m_psf, m_v1);   // (A^dag A)^{-1/4}

    double xi2 = xi.norm2();

    int ith = ThreadManager::get_thread_id();
    if(ith == 0) H_psf = xi2;
  }

  timer.stop();
  double elapsed_time = timer.elapsed_sec();

  vout.general(m_vl, "    H_Frational  = %18.8f\n", H_psf);
  vout.general(m_vl, "    H_F/dof      = %18.8f\n", H_psf / size_psf);
  vout.general(m_vl, "    elapsed time: %12.6f [sec]\n", elapsed_time);

  // for confirmation
  // H_psf = calcH();

  return H_psf;

}

//====================================================================
template<typename AFIELD>
double Action_F_alt_Rational_Ratio<AFIELD>::calcH()
{
  const int NinF     = m_fopr1_rational_H->field_nin();
  const int NvolF    = m_fopr1_rational_H->field_nvol();
  const int NexF     = m_fopr1_rational_H->field_nex();
  const int size_psf = NinF * NvolF * NexF * CommonParameters::NPE();

  vout.general(m_vl, "  %s: %s\n", class_name.c_str(), m_label.c_str());

  Timer timer;
  timer.start();

  double H_psf;

#pragma omp parallel
 {
  m_fopr1_rational_H->set_config(m_U);
  m_fopr1_rational_H->mult(m_v1, m_psf);  // (A^dag A)^{1/4}

  m_fopr2_rational_H->set_config(m_U);
  m_fopr2_rational_H->mult(m_v2, m_v1);   // (B^dag B)^{-1/2}

  double H_psf1 = dot(m_v1, m_v2);

  int ith = ThreadManager::get_thread_id();
  if(ith == 0) H_psf = H_psf1;
 }

  timer.stop();
  double elapsed_time = timer.elapsed_sec();

  vout.general(m_vl, "    H_Frational  = %18.8f\n", H_psf);
  vout.general(m_vl, "    H_F/dof      = %18.8f\n", H_psf / size_psf);
  vout.general(m_vl, "    elapsed time: %12.6f [sec]\n", elapsed_time);

  return H_psf;

}

//====================================================================
template<typename AFIELD>
void Action_F_alt_Rational_Ratio<AFIELD>::force(Field& force)
{
  vout.general(m_vl, "  %s: %s\n", class_name.c_str(), m_label.c_str());

  Timer timer;
  timer.start();

  const int Nin  = m_U->nin();
  const int Nvol = m_U->nvol();
  const int Nex  = m_U->nex();

  assert(force.nin() == Nin);
  assert(force.nvol() == Nvol);
  assert(force.nex() == Nex);

  AFIELD force1(Nin, Nvol, Nex);

  // solvers are in the force
  #pragma omp parallel
  {

    m_force_MD->set_config(m_U);

    m_force_MD->force_core(force1, m_psf);
    AIndex_lex<real_t,AFIELD::IMPL> index_lex;
    reverse_gauge(index_lex, force, force1);

    double Fave, Fmax, Fdev;
    force.stat(Fave, Fmax, Fdev);
    vout.general(m_vl, "    Fave = %12.6f  Fmax = %12.6f  Fdev = %12.6f\n",
                 Fave, Fmax, Fdev);
  }

  timer.stop();
  double elapsed_time = timer.elapsed_sec();
  vout.general(m_vl, "    elapsed time: %12.6f [sec]\n", elapsed_time);

}

//============================================================END=====
