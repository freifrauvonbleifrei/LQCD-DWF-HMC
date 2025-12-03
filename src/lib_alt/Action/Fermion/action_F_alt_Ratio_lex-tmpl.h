/*!
      @file    action_F_alt_Ratio_lex-tmpl.h
      @brief
      @author  Hideo Matsufuru (matufuru)
               $LastChangedBy: kanamori $
      @date    $LastChangedDate:: 2025-12-03 19:35:35 #$
      @version $LastChangedRevision: 2668 $
*/

template<typename AFIELD>
const std::string Action_F_alt_Ratio_lex<AFIELD>::class_name
                                          = "Action_F_alt_Ratio_lex";

//====================================================================
template<typename AFIELD>
void Action_F_alt_Ratio_lex<AFIELD>::init()
{
  m_vl = CommonParameters::Vlevel();

}

//====================================================================
template<typename AFIELD>
void Action_F_alt_Ratio_lex<AFIELD>::tidyup()
{
  // do nothing.
}

//====================================================================
template<typename AFIELD>
void Action_F_alt_Ratio_lex<AFIELD>::set_parameters(const Parameters& params)
{
  const string str_vlevel = params.get_string("verbose_level");
  m_vl = vout.set_verbose_level(str_vlevel);
}

//====================================================================
template<typename AFIELD>
void Action_F_alt_Ratio_lex<AFIELD>::set_parameters()
{
}

//====================================================================
template<typename AFIELD>
void Action_F_alt_Ratio_lex<AFIELD>::get_parameters(
                                         Parameters& params) const
{
  params.set_string("verbose_level", vout.get_verbose_level(m_vl));
}

//====================================================================
template<typename AFIELD>
void Action_F_alt_Ratio_lex<AFIELD>::set_config(Field *U)
{
  m_U = U;

  m_fopr_prec->set_config(U);
  m_fopr_prec_force->set_config(U);

  m_fopr->set_config(U);
  m_fopr_force->set_config(U);

}

//====================================================================
template<typename AFIELD>
double Action_F_alt_Ratio_lex<AFIELD>::langevin(RandomNumbers *rand)
{
  const int Nvol = CommonParameters::Nvol();
  const int Ndim = CommonParameters::Ndim();

  const int NinF     = m_fopr_prec->field_nin();
  const int NvolF    = m_fopr_prec->field_nvol();
  const int NexF     = m_fopr_prec->field_nex();
  const int size_psf = NinF * NvolF * NexF * CommonParameters::NPE();

  assert(NvolF == Nvol);
  m_psf.reset(NinF, NvolF, NexF);

  Timer timer;
  timer.start();

  vout.general(m_vl, "  %s: %s\n", class_name.c_str(), m_label.c_str());

  const int Nc=CommonParameters::Nc();
  const int Nd=CommonParameters::Nd();
  int Nin=NinF;
  int Nex=NexF;
  if(2*Nc*Nd != NinF){ // temporal hack;
    Nin = 2*Nc*Nd;
    Nex = NinF/Nin;
  }
  Field xi(Nin, NvolF, Nex);
  rand->gauss_lex_global(xi);

  AFIELD xiA(NinF, NvolF, NexF);

  AFIELD v2(NinF, NvolF, NexF);
  AFIELD v1(NinF, NvolF, NexF);

#pragma omp parallel
 {
  if(m_fopr->needs_convert()){
    m_fopr->convert(xiA, xi);
  }else{
    AIndex_lex<real_t,AFIELD::IMPL> index_alt;
    convert(index_alt, xiA, xi);
  }

  m_fopr->set_config(m_U);

  m_fopr->set_mode("H");
  m_fopr->mult_dag(v2, xiA);
 }

  int    Nconv;
  double diff;
  m_fprop_H_prec->set_config(m_U);
  m_fprop_H_prec->set_mode("DdagD");
  m_fprop_H_prec->invert(v1, v2, Nconv, diff);
  vout.general(m_vl, "    Fprop_H_prec: Nconv = %6d  diff  = %12.6e\n",
               Nconv, diff);

  double H_psf;

#pragma omp parallel
 {
  m_fopr_prec->set_mode("H");
  m_fopr_prec->mult(m_psf, v1);

  double H_psf1 = xi.norm2();

  int ith = ThreadManager::get_thread_id();
  if(ith == 0) H_psf = H_psf1;
 }

  timer.stop();
  double elapsed_time = timer.elapsed_sec();

  vout.general(m_vl, "    H_Fratio     = %18.8f\n", H_psf);
  vout.general(m_vl, "    H_F/dof      = %18.8f\n", H_psf / size_psf);
  vout.general(m_vl, "    elapsed time: %12.6f [sec]\n", elapsed_time);

  return H_psf;

}

//====================================================================
template<typename AFIELD>
double Action_F_alt_Ratio_lex<AFIELD>::calcH()
{
  const int Nvol = CommonParameters::Nvol();
  const int Ndim = CommonParameters::Ndim();

  const int NinF     = m_fopr_prec->field_nin();
  const int NvolF    = m_fopr_prec->field_nvol();
  const int NexF     = m_fopr_prec->field_nex();
  const int size_psf = NinF * NvolF * NexF * CommonParameters::NPE();

  Timer timer;
  double elapsed_time;
  timer.start();

  vout.general(m_vl, "  %s: %s\n", class_name.c_str(), m_label.c_str());

  m_fopr_prec->set_config(m_U);
  m_fopr->set_config(m_U);

  AFIELD v1(NinF, NvolF, NexF);
  m_fopr_prec->set_mode("H");

#pragma omp parallel
 {
  m_fopr_prec->mult_dag(v1, m_psf);
 }

  AFIELD v2(NinF, NvolF, NexF);
  int    Nconv;
  double diff;

  m_fprop_H->set_config(m_U);
  m_fprop_H->reset_performance();
  m_fprop_H->set_mode("DdagD");
  m_fprop_H->invert(v2, v1, Nconv, diff);

  double flop_count, elapsed_time2, gflops;
  m_fprop_H->get_performance(flop_count, elapsed_time2);
  if(elapsed_time2 > 0.0){
    gflops = (flop_count/elapsed_time2) * 1.0e-9;
  }else{
    gflops = 0.0;
  }

  vout.general(m_vl, "    Fprop_H: %6d  %12.4e  %12.4e [GFlops]\n",
               Nconv, diff, gflops);

  double H_psf;
#pragma	omp parallel
 {
  double H_psf1 = dot(v1, v2);

  int ith = ThreadManager::get_thread_id();
  if(ith == 0) H_psf = H_psf1;
 }

  timer.stop();
  elapsed_time = timer.elapsed_sec();

  vout.general(m_vl, "    H_Fratio     = %18.8f\n", H_psf);
  vout.general(m_vl, "    H_F/dof      = %18.8f\n", H_psf / size_psf);
  vout.general(m_vl, "    elapsed time: %12.6f [sec]\n", elapsed_time);

  return H_psf;
}


//====================================================================
template<typename AFIELD>
void Action_F_alt_Ratio_lex<AFIELD>::force(Field& force)
{
  const int Nin  = m_U->nin();
  const int Nvol = m_U->nvol();
  const int Nex  = m_U->nex();
  const int Nc   = CommonParameters::Nc();
  const int Ndim = CommonParameters::Ndim();

  assert(force.nin() == Nin);
  assert(force.nvol() == Nvol);
  assert(force.nex() == Nex);

  AFIELD force1(Nin, Nvol, Nex);

  int NinF  = m_fopr_prec->field_nin();
  int NvolF = m_fopr_prec->field_nvol();
  int NexF  = m_fopr_prec->field_nex();

  Timer timer;
  timer.start();

  vout.general(m_vl, "  %s: %s\n", class_name.c_str(), m_label.c_str());

  m_fopr_prec->set_config(m_U);
  m_fopr->set_config(m_U);
  m_fopr_prec_force->set_config(m_U);
  m_fopr_force->set_config(m_U);

  AFIELD v1(NinF, NvolF, NexF);
  AFIELD v2(NinF, NvolF, NexF);

#pragma omp parallel
 {
  m_fopr_prec->set_mode("H");
  m_fopr_prec->mult_dag(v1, m_psf);
 }

  int    Nconv;
  double diff;
  m_fprop_MD->set_config(m_U);
  m_fprop_MD->set_mode("DdagD");
  m_fprop_MD->invert(v2, v1, Nconv, diff);
  vout.general(m_vl, "    Fprop_MD: %6d %18.15e\n", Nconv, diff);

  AFIELD force_tmp(Nin, Nvol, Nex);

#pragma omp parallel
{
  m_fopr_force->force_core(force1, v2);

  m_fopr_prec_force->set_mode("Hdag");
  m_fopr_prec_force->force_core1(force_tmp, v2, m_psf);

  axpy(force1, -1.0, force_tmp);

  m_fopr_prec_force->set_mode("H");
  m_fopr_prec_force->force_core1(force_tmp, m_psf, v2);

  axpy(force1, -1.0, force_tmp);

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
