/*!
      @file    action_F_alt_Rational_Ratio.h
      @brief
      @author  Hideo Matsufuru (matufuru)
               $LastChangedBy: kanamori $
      @date    $LastChangedDate:: 2025-12-03 19:35:35 #$
      @version $LastChangedRevision: 2668 $
*/

#ifndef ACTION_F_ALT_RATIONAL_RATIO_INCLUDED
#define ACTION_F_ALT_RATIONAL_RATIO_INCLUDED

#include "lib/Action/action.h"

#include "lib/Fopr/afopr_Rational.h"
#include "lib_alt/Force/Fermion/aforce_F_Rational.h"

#include "lib/IO/bridgeIO.h"
using Bridge::vout;

// include files in alt-code
//#include "lib_alt/Field/afield.h"
//#include "lib_alt/Measurements/Fermion/afprop.h"

//! action class for RHMC, with externally constructed AFopr_Rational.

/*!
    Action class for RHMC that is an alternative to Action_F_Rational
    in the core library.
                                            [05 Feb 2019 H.Matsufuru]

  S = phi^dag (A^dag A)^{1/4} (B^dag B)^{-1/2} (A^dag A)^{1/4} \phi
    A: "Fopr1" (numerator)   the rational parameter must be exponet = -1/2
    B: "Fopr2" (denominator) the rational parameter must be exponent = 1/4

    The input Fopr must be A and B (not (A^dag A)^{1/4} etc. )
    This class requires the forces constructed externally specific
    to this action.
                                        [ 15 Aug 2023 I.Kanamori]
*/



template<typename AFIELD>
class Action_F_alt_Rational_Ratio : public Action
{
 public:
  //typedef AField<double> AFIELD;
  typedef typename AFIELD::real_t real_t;
  static const std::string class_name;

 private:
  std::string m_label;    // label of action


  AFopr<AFIELD>    *m_fopr1, *m_fopr2;
  AForce_F<AFIELD> *m_force_MD;

  unique_ptr<AFopr<AFIELD> > m_fopr1_rational_H;
  unique_ptr<AFopr<AFIELD> > m_fopr1_rational_langevin;
  unique_ptr<AFopr<AFIELD> > m_fopr2_rational_H;
  unique_ptr<AFopr<AFIELD> > m_fopr2_rational_langevin;


  Field *m_U;

  AFIELD m_psf;
  AFIELD m_v1, m_v2;

  Parameters m_params_rational1_H;
  Parameters m_params_rational1_langevin;
  Parameters m_params_rational2_H;
  Parameters m_params_rational2_langevin;


  Bridge::VerboseLevel m_vl;   //!< verbose level

 public:
  //! constructor.
  Action_F_alt_Rational_Ratio(AFopr<AFIELD> *fopr1_H,
                              AFopr<AFIELD> *fopr2_H,
                              AForce_F<AFIELD> *force_MD)
    : Action(),
      m_fopr1(fopr1_H), m_fopr2(fopr2_H), m_force_MD(force_MD)
      { init(); }

  Action_F_alt_Rational_Ratio(AFopr<AFIELD> *fopr1_H,
                              AFopr<AFIELD> *fopr2_H,
                              AForce_F<AFIELD> *force_MD,
                              Parameters &params)
    : Action(),
      m_fopr1(fopr1_H), m_fopr2(fopr2_H), m_force_MD(force_MD)
      {
        init();
        set_parameters(params);
      }

  //! destructor.
  ~Action_F_alt_Rational_Ratio() { tidyup(); }

  void set_parameters(const Parameters& params);

  void get_parameters(Parameters&) const;

  //! set the label of action.
  void set_label(const std::string label)
  {
    m_label = label;
    vout.detailed(m_vl, "  label: %s\n", m_label.c_str());
  }

  //! returns the label of action.
  std::string get_label()
  { return m_label; }

  //! setting gauge configuration.
  void set_config(Field *U);

  //! Langevin step called at the beginning of HMC.
  double langevin(RandomNumbers *);

  //! calculation of Hamiltonian.
  double calcH();

  //! returns the force for updating conjugate momentum.
  //const Field force();
  void force(Field&);

 private:
  void init();
  void tidyup();

};
#endif
