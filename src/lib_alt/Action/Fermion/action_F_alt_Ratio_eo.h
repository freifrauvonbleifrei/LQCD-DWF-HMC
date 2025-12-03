/*!
      @file    action_F_alt_Ratio_eo.h
      @brief
      @author  Hideo Matsufuru (matufuru)
               $LastChangedBy: kanamori $
      @date    $LastChangedDate:: 2025-12-03 19:35:35 #$
      @version $LastChangedRevision: 2668 $
*/

#ifndef ACTION_F_ALT_RATIO_EO_INCLUDED
#define ACTION_F_ALT_RATIO_EO_INCLUDED

#include "lib/Action/action.h"

#include "lib/Force/Fermion/force_F.h"
#include "lib/IO/bridgeIO.h"
using Bridge::vout;

#include "lib_alt/Measurements/Fermion/fprop_alt.h"

//! HMC action for Hasenbusch preconditioned fermions.

/*!
    Hasenbusch preconditioned fermion action in alternative
    implementation.
                                    [04 Feb 2019 H.Matsufuru]

    takes external force to avoid redundant force smearing
                                    [25 Sep 2023 I.Kanamori ]
*/

template<typename AFIELD>
class Action_F_alt_Ratio_eo : public Action
{
 public:
  //  typedef AField<double> AFIELD;
  typedef typename AFIELD::real_t real_t;
  static const std::string class_name;

 private:
  Bridge::VerboseLevel m_vl;      //!< verbose level

  AFopr<AFIELD>    *m_fopr1;      //!< preconditioner
  AFopr<AFIELD>    *m_fopr2;      //!< dynamical fermion
  AFIELD m_psf;                   //!< pseudofermion field
  std::string m_label;            //!< label of action

  Fprop_alt<AFIELD> *m_fprop1_H;
  Fprop_alt<AFIELD> *m_fprop2_MD;
  Fprop_alt<AFIELD> *m_fprop2_H;

  AForce_F<AFIELD> *m_force_MD;   //!< force of the ratio operator

  Field *m_U;

 public:
  Action_F_alt_Ratio_eo(
    AFopr<AFIELD> *fopr1, AFopr<AFIELD> *fopr2,
    Fprop_alt<AFIELD> *fprop1_H,
    Fprop_alt<AFIELD> *fprop2_MD, Fprop_alt<AFIELD> *fprop2_H,
    AForce_F<AFIELD> *force_MD )
    : Action(),
      m_fopr1(fopr1), m_fopr2(fopr2),
      m_fprop1_H(fprop1_H), m_fprop2_MD(fprop2_MD), m_fprop2_H(fprop2_H),
     m_force_MD(force_MD)
  {
    init();
  }



  ~Action_F_alt_Ratio_eo()
  { tidyup();}

  void set_parameters(const Parameters&);
  void set_parameters();

  void get_parameters(Parameters&) const;

  void set_label(const std::string label)
  {
    m_label = label;
    vout.detailed(m_vl, "  label: %s\n", m_label.c_str());
  }

  std::string get_label()
  { return m_label; }

  void set_config(Field *U);

  double langevin(RandomNumbers *);

  double calcH();

  void force(Field&);

 private:
  void init();
  void tidyup();

};
#endif
