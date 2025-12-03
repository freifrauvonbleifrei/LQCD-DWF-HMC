/*!
      @file    action_F_alt_Ratio_eo-tmpl.h
      @brief
      @author  Hideo Matsufuru (matufuru)
               $LastChangedBy: kanamori $
      @date    $LastChangedDate:: 2025-12-03 19:35:35 #$
      @version $LastChangedRevision: 2668 $
*/

template<typename AFIELD>
const std::string Action_F_alt_Ratio_eo<AFIELD>::class_name
                                          = "Action_F_alt_Ratio_eo";

#include "lib_alt/Measurements/Fermion/fprop_alt_Standard_eo.h"
//====================================================================
template<typename AFIELD>
void Action_F_alt_Ratio_eo<AFIELD>::init()
{
  m_vl = CommonParameters::Vlevel();

}

//====================================================================
template<typename AFIELD>
void Action_F_alt_Ratio_eo<AFIELD>::tidyup()
{
  // do nothing.
}

//====================================================================
template<typename AFIELD>
void Action_F_alt_Ratio_eo<AFIELD>::set_parameters(const Parameters& params)
{
  const string str_vlevel = params.get_string("verbose_level");
  m_vl = vout.set_verbose_level(str_vlevel);
}

//====================================================================
template<typename AFIELD>
void Action_F_alt_Ratio_eo<AFIELD>::set_parameters()
{
}

//====================================================================
template<typename AFIELD>
void Action_F_alt_Ratio_eo<AFIELD>::get_parameters(
                                         Parameters& params) const
{
  params.set_string("verbose_level", vout.get_verbose_level(m_vl));
}

//====================================================================
template<typename AFIELD>
void Action_F_alt_Ratio_eo<AFIELD>::set_config(Field *U)
{
  m_U = U;

  m_fopr1->set_config(U);
  m_fopr2->set_config(U);
  m_force_MD->set_config(U);

}

//====================================================================
template<typename AFIELD>
double Action_F_alt_Ratio_eo<AFIELD>::langevin(RandomNumbers *rand)
{
  const int Nx   = CommonParameters::Nx();
  const int Ny   = CommonParameters::Ny();
  const int Nz   = CommonParameters::Nz();
  const int Nt   = CommonParameters::Nt();
  const int Nvol = CommonParameters::Nvol();
  const int Ndim = CommonParameters::Ndim();
  const int Nx2 = Nx/2;

  const int NinF     = m_fopr1->field_nin();
  const int NvolF    = m_fopr1->field_nvol();
  const int NexF     = m_fopr1->field_nex();
  const int size_psf = NinF * NvolF * NexF * CommonParameters::NPE();

  assert((2 * NvolF) == Nvol);
  m_psf.reset(NinF, NvolF, NexF);

  Timer timer;
  timer.start();

  vout.general(m_vl, "  %s: %s\n", class_name.c_str(), m_label.c_str());

  Field xi_lex;
  const int Nc=CommonParameters::Nc();
  const int Nd=CommonParameters::Nd();
  int Nin=NinF;
  int Nex=NexF;
  if(2*Nc*Nd != NinF){ // temporal hack;
    Nin = 2*Nc*Nd;
    Nex = NinF/Nin;
  }
  xi_lex.reset(Nin, Nvol, Nex);
  rand->gauss_lex_global(xi_lex);

  Field xi(Nin, NvolF, Nex);

  Index_eo idx_eo;

  AFIELD xiA(NinF, NvolF, NexF);

  AFIELD v2(NinF, NvolF, NexF);
  AFIELD v1(NinF, NvolF, NexF);
  AForce_F<AFIELD> *m_forceF2;    //!< force of dynamical fermion

  double H_psf;

#pragma omp parallel
 {
  idx_eo.convertField(xi, xi_lex, 0);

  if(m_fopr1->needs_convert()){
    m_fopr1->convert(xiA, xi);
  }else{
    AIndex_lex<real_t,AFIELD::IMPL> index_alt(Nx2, Ny, Nz, Nt);
    convert(index_alt, xiA, xi);
  }

  m_fopr1->set_config(m_U);
  m_fopr2->set_config(m_U);

  m_fopr2->set_mode("H");
  m_fopr2->mult_dag(v2, xiA);

 }

  int    Nconv;
  double diff;
  m_fprop1_H->set_config(m_U);
  m_fprop1_H->set_mode("DdagD_even");
  Fprop_alt_Standard_eo<AFIELD> *fprop=dynamic_cast< Fprop_alt_Standard_eo<AFIELD>* >(m_fprop1_H);

  m_fprop1_H->invert(v1, v2, Nconv, diff);

  vout.general(m_vl, "    Fprop_H_prec: %6d %18.15e\n", Nconv, diff);

#pragma omp parallel
 {
   m_fopr1->set_mode("H");
   m_fopr1->mult(m_psf, v1);

  real_t H_psf1 = xiA.norm2();

  int ith = ThreadManager::get_thread_id();
  if(ith == 0) H_psf = H_psf1;

 }

  timer.stop();
  double elapsed_time = timer.elapsed_sec();

  vout.general(m_vl, "    H_Fratio     = %18.8f\n", H_psf);
  vout.general(m_vl, "    H_F/dof      = %18.8f\n", H_psf / size_psf);
  vout.general(m_vl, "    elapsed time: %12.6f [sec]\n", elapsed_time);

  double Halt = calcH();

  return H_psf;

}

//====================================================================
template<typename AFIELD>
double Action_F_alt_Ratio_eo<AFIELD>::calcH()
{
  const int Nvol = CommonParameters::Nvol();
  const int Ndim = CommonParameters::Ndim();

  const int NinF     = m_fopr1->field_nin();
  const int NvolF    = m_fopr1->field_nvol();
  const int NexF     = m_fopr1->field_nex();
  const int size_psf = NinF * NvolF * NexF * CommonParameters::NPE();

  Timer timer;
  double elapsed_time;
  timer.start();

  vout.general(m_vl, "  %s: %s\n", class_name.c_str(), m_label.c_str());

  AFIELD v1(NinF, NvolF, NexF);
  AFIELD v2(NinF, NvolF, NexF);

#pragma omp parallel
 {
  m_fopr1->set_config(m_U);

  m_fopr1->set_mode("H");
  m_fopr1->mult_dag(v1, m_psf);
 }

  int    Nconv;
  double diff;
  Timer timer_solver;

  m_fprop2_H->set_config(m_U);
  m_fprop2_H->reset_performance();
  m_fprop2_H->set_mode("DdagD_even");
  timer_solver.start();
  m_fprop2_H->invert(v2, v1, Nconv, diff);
  timer_solver.stop();
  double elapsed_time_solver = timer_solver.elapsed_sec();


  double H_psf;

#pragma omp parallel
 {
  real_t H_psf1 = dot(v1, v2);

  int ith = ThreadManager::get_thread_id();
  if(ith == 0) H_psf = H_psf1;
 }

 double flop_count, elapsed_time2, gflops;
  m_fprop2_H->get_performance(flop_count, elapsed_time2);
  if(elapsed_time2 > 0.0){
    gflops = (flop_count/elapsed_time2) * 1.0e-9;
  }else{
    gflops = 0.0;
  }

  vout.general(m_vl, "    Fprop_H: %6d  %12.4e  %12.4e [GFlops]\n",
               Nconv, diff, gflops);
  vout.general(m_vl, "      elapsed time: %12.6f [sec]\n", elapsed_time_solver);

  timer.stop();
  elapsed_time = timer.elapsed_sec();

  vout.general(m_vl, "    H_Fratio     = %18.8f\n", H_psf);
  vout.general(m_vl, "    H_F/dof      = %18.8f\n", H_psf / size_psf);
  vout.general(m_vl, "    elapsed time: %12.6f [sec]\n", elapsed_time);

  return H_psf;
}



//====================================================================
template<typename AFIELD>
void Action_F_alt_Ratio_eo<AFIELD>::force(Field& force)
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

  int NinF  = m_fopr1->field_nin();
  int NvolF = m_fopr1->field_nvol();
  int NexF  = m_fopr1->field_nex();

  AFIELD v1(NinF, NvolF, NexF);
  AFIELD v2(NinF, NvolF, NexF);

  m_fopr1->set_config(m_U);
  m_fopr1->set_mode("H");

#pragma omp parallel
  {
    m_fopr1->mult_dag(v1, m_psf);
  }

  {
    int    Nconv;
    double diff;

    Timer timer_solver;

    m_fprop2_MD->set_config(m_U);
    m_fprop2_MD->set_mode("DdagD_even");
    timer_solver.start();
    m_fprop2_MD->invert(v2, v1, Nconv, diff);
    timer_solver.stop();
    double elapsed_time = timer_solver.elapsed_sec();

    vout.general(m_vl, "    Fprop_MD: %6d %18.15e\n", Nconv, diff);
    vout.general(m_vl, "      elapsed time: %12.6f [sec]\n", elapsed_time);
  }


#pragma omp parallel
  {
    m_force_MD->set_config(m_U);

    m_force_MD->force_core1(force1, v2, m_psf);

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
