/*!
        @file    aforce_F_Ratio_eo.h

        @brief

        @author  Issaku Kanamori (kanamori)
                 $LastChangedBy: kanamori $

        @date    $LastChangedDate:: 2025-12-03 19:35:35 #$

        @version $LastChangedRevision: 2668 $
*/

#ifndef AFORCE_F_RATIO_EO_INCLUDED
#define AFORCE_F_RATIO_EO_INCLUDED

#include "Force/Fermion/aforce_F.h"

#include "Field/field_F.h"
#include "Fopr/afopr.h"
#include "lib_alt/Measurements/Fermion/fprop_alt.h"

#include "IO/bridgeIO.h"
using Bridge::vout;

//! Force calculation for smeared fermion operators.

/*!
    This class determines the force of a ratio of
    given fermion operators,
       S = phi^dag A (B^dag B)^{-1} A^dag phi
       A: numerator (fopr1)
       A: denominator (fopr2)
                                     [25 Sep 2023 I.Kanamori]
 */

template<typename AFIELD>
class AForce_F_Ratio_eo : public AForce_F<AFIELD>
{
 public:
  static const std::string class_name;
  typedef typename AFIELD::real_t real_t;

 private:
  Bridge::VerboseLevel m_vl;

  Field_G *m_U;

  AFopr<AFIELD>    *m_fopr1;      //!< preconditioner  (numerator)
  AForce_F<AFIELD> *m_forceF1;    //!< force of preconditioner

  AForce_F<AFIELD> *m_forceF2;    //!< force of dynamical fermion

  Fprop_alt<AFIELD> *m_fprop2_MD;

  AFIELD m_v1, m_v2;
  AFIELD m_force;

public:

  AForce_F_Ratio_eo(
    AFopr<AFIELD> *fopr1,  // dummy, only to obtain field size
    AForce_F<AFIELD> *forceF1, AForce_F<AFIELD> *forceF2)
    : AForce_F<AFIELD>(),
      m_fopr1(fopr1), m_forceF1(forceF1),
      m_forceF2(forceF2)
  {
    m_vl = CommonParameters::Vlevel();
    init();
  }


  AForce_F_Ratio_eo(
    AFopr<AFIELD> *fopr1,
    AForce_F<AFIELD> *forceF1,AForce_F<AFIELD> *forceF2,
    const Parameters& params)
    : AForce_F<AFIELD>(),
      m_fopr1(fopr1),
      m_forceF1(forceF1),m_forceF2(forceF2)
  {
    m_vl = CommonParameters::Vlevel();
    init(params);
  }


  ~AForce_F_Ratio_eo() { tidyup(); }

  void set_parameters(const Parameters& params);
  void set_parameters();

  void get_parameters(Parameters& params) const;

  void set_config(Field *U);

  void force_udiv(AFIELD&, const AFIELD&);

  void force_core(AFIELD&, const AFIELD&);
  void force_udiv1(AFIELD&, const AFIELD&, const AFIELD&);

 private:

  void init();
  void init(const Parameters &params);

  void tidyup();

  void force_udiv_impl(AFIELD&, const AFIELD&);
  void force_udiv1_impl(AFIELD&, const AFIELD&, const AFIELD&);

#ifdef USE_FACTORY
 private:
  //static AForce_F<AFIELD> *create_object(AFopr<AFIELD> *fopr,
  //                                       AForce_F<AFIELD> *force_F)
  //{ return new AForce_F_Ratio_eo(fopr, force_F); }

 public:
  //static bool register_factory()
  //{
  //  bool init1 = AForce_F<AFIELD>::Factory_fopr_force::Register(
  //                                         "Ratio_eo", create_object);
  //  return init1;
  // }
#endif


};
#endif
