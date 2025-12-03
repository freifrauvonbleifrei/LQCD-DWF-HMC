/*!
        @file    aforce_F_Ratio_eo-tmpl.h

        @brief

        @author  Issaku Kanamori (kanamori)
                 $LastChangedBy: kanamori $

        @date    $LastChangedDate:: 2025-12-03 19:35:35 #$

        @version $LastChangedRevision: 2668 $
*/

#include "aforce_F_Ratio_eo.h"
#include "lib/Tools/timer.h"

template<typename AFIELD>
const std::string AForce_F_Ratio_eo<AFIELD>::class_name
                                             = "AForce_F_Ratio_eo";

//====================================================================
template<typename AFIELD>
void AForce_F_Ratio_eo<AFIELD>::init(const Parameters& params)
{
  ThreadManager::assert_single_thread(class_name);

  std::string vlevel;
  if (!params.fetch_string("verbose_level", vlevel)) {
    m_vl = vout.set_verbose_level(vlevel);
  } else {
    m_vl = CommonParameters::Vlevel();
  }
  init();
}

//====================================================================
template<typename AFIELD>
void AForce_F_Ratio_eo<AFIELD>::init()
{
  ThreadManager::assert_single_thread(class_name);
  vout.general(m_vl, "%s: construction\n", class_name.c_str());
  vout.increase_indent();

  //  set_parameters(params); // nothing to set

  int NinF  = m_fopr1->field_nin();
  int NvolF = m_fopr1->field_nvol();
  int NexF  = m_fopr1->field_nex();

  m_v1.reset(NinF, NvolF, NexF);
  m_v2.reset(NinF, NvolF, NexF);



  int Nc   = CommonParameters::Nc();
  int NinG = 2 * Nc * Nc;
  int Nvol = CommonParameters::Nvol();
  int Ndim = CommonParameters::Ndim();

  m_force.reset(NinG, Nvol, Ndim);

  vout.decrease_indent();
  vout.general(m_vl, "%s: construction finished.\n",
               class_name.c_str());

}

//====================================================================
template<typename AFIELD>
void AForce_F_Ratio_eo<AFIELD>::tidyup()
{
  //  delete m_solver;

}

//====================================================================
template<typename AFIELD>
void AForce_F_Ratio_eo<AFIELD>::set_parameters(const Parameters& params)
{
  std::string vlevel;
  if (!params.fetch_string("verbose_level", vlevel)) {
    m_vl = vout.set_verbose_level(vlevel);
  } else {
    m_vl = CommonParameters::Vlevel();
  }

}


//====================================================================
template<typename AFIELD>
void AForce_F_Ratio_eo<AFIELD>::get_parameters(Parameters& params) const
{
  params.set_string("verbose_level", vout.get_verbose_level(m_vl));
}

//====================================================================
template<typename AFIELD>
void AForce_F_Ratio_eo<AFIELD>::set_config(Field *U)
{
  m_U = (Field_G *)U;
  m_forceF1->set_config(U);
  m_forceF2->set_config(U);
}

//====================================================================
template<typename AFIELD>
void AForce_F_Ratio_eo<AFIELD>::force_udiv(AFIELD& force,
                                           const AFIELD& eta)
{
  force_udiv_impl(force, eta);
}

//====================================================================
template<typename AFIELD>
void AForce_F_Ratio_eo<AFIELD>::force_udiv1(AFIELD& force,
                                            const AFIELD& zeta,
                                            const AFIELD& eta)
{
  force_udiv1_impl(force, zeta, eta);
}


//====================================================================
template<typename AFIELD>
void AForce_F_Ratio_eo<AFIELD>::force_udiv1_impl(AFIELD& force,
                                                const AFIELD& zeta,
                                                const AFIELD& eta)
{
#pragma omp barrier

  const int Nin  = m_U->nin();
  const int Nvol = m_U->nvol();
  const int Nex  = m_U->nex();

  assert(force.nin() == Nin);
  assert(force.nvol() == Nvol);
  assert(force.nex() == Nex);

  m_forceF1->set_config(m_U);
  m_forceF2->set_config(m_U);

  m_forceF2->force_udiv(force, zeta);

  m_forceF1->set_mode("Hdag");
  m_forceF1->force_udiv1(m_force, zeta, eta);

  axpy(force, -1.0, m_force);
#pragma omp barrier

  m_forceF1->set_mode("H");
  m_forceF1->force_udiv1(m_force, eta, zeta);

  axpy(force, -1.0, m_force);
#pragma omp barrier

}


//====================================================================
template<typename AFIELD>
void AForce_F_Ratio_eo<AFIELD>::force_udiv_impl(AFIELD& force,
                                                const AFIELD& eta)
{
  vout.crucial(m_vl, "Error at %s: %s is not implemented.\n",
               class_name.c_str(), __func__);
  exit(EXIT_FAILURE);
}



//====================================================================
template<typename AFIELD>
void AForce_F_Ratio_eo<AFIELD>::force_core(
                             AFIELD&, const AFIELD&)
{
  vout.crucial(m_vl, "Error at %s: %s not implemented.\n",
               class_name.c_str(), __func__);
  exit(EXIT_FAILURE);
}


//===========================================================END======
