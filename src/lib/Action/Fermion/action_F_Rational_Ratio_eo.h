/*!
        @file    action_F_Rational_Ratio_eo.h

        @brief

        @author  Issaku Kanamori (kanamori)
                 $LastChangedBy: kanamori $

        @date    $LastChangedDate:: 2025-12-03 19:35:35 #$

        @version $LastChangedRevision: 2668 $
*/

#ifndef ACTION_F_RATIONAL_RATIO_EO_INCLUDED
#define ACTION_F_RATIONAL_RATIO_EO_INCLUDED

#include "Action/action.h"

#include "Solver/shiftsolver_CG.h"
#include "Tools/math_Rational.h"
#include "Force/Fermion/force_F_Rational.h"

#include "IO/bridgeIO.h"
using Bridge::vout;

//! action class for RHMC, ratio of two rational fanctions [even-odd version]

/*!
  S = phi^dag (A^dag A)^{1/4} (B^dag B)^{-1/2} (A^dag A)^{1/4} \phi
    A: numerator, the rational parameter must be exponet = -1/2
    B: denominator, the rational parameter must be exponent = 1/4

    This input Fopr must be A and B (not (A^dag A)^{1/4} etc. )
    This class requires the forces constructed externally specific
    to this action.
                                        [ 1 Jul 2023 I.Kanamori]
*/


class Action_F_Rational_Ratio_eo : public Action
{
 public:
  static const std::string class_name;

 private:
  Bridge::VerboseLevel m_vl;

  std::string m_label;    // label of action

  Fopr *m_fopr_denom, *m_fopr_numer;
  Force *m_fopr_force_MD;

  unique_ptr<Fopr_Rational> m_fopr_denom_H;
  unique_ptr<Fopr_Rational> m_fopr_denom_psf;
  unique_ptr<Fopr_Rational> m_fopr_numer_H;
  unique_ptr<Fopr_Rational> m_fopr_numer_psf;

  Field *m_U;

  Field m_psf;

  Parameters m_params_rational_denom_H;
  Parameters m_params_rational_denom_psf;
  Parameters m_params_rational_numer_H;
  Parameters m_params_rational_numer_psf;


public:
  //! constructor requires pointers to Fopr and Force instances.
  Action_F_Rational_Ratio_eo(Fopr *fopr_numer, Fopr *fopr_denom, Force *fopr_force_MD)
    : m_vl(CommonParameters::Vlevel()),
      m_fopr_numer(fopr_numer), m_fopr_denom(fopr_denom),
      m_fopr_force_MD(fopr_force_MD)
  {}

  Action_F_Rational_Ratio_eo(Fopr *fopr_numer, Fopr *fopr_denom, Force *fopr_force_MD,
                             const Parameters& params)
    : m_vl(CommonParameters::Vlevel()),
      m_fopr_numer(fopr_numer), m_fopr_denom(fopr_denom),
      m_fopr_force_MD(fopr_force_MD)
  {
    set_parameters(params);
  }

  //! destructor. constructed instances are deconstructed in tydyup().
  ~Action_F_Rational_Ratio_eo() {}

  void set_parameters(const Parameters& params);

  void get_parameters(Parameters& params) const;

  //! set the label of action.
  void set_label(const std::string label)
  {
    m_label = label;
    vout.detailed(m_vl, "  label: %s\n", m_label.c_str());
  }

  //! returns the label of action.
  std::string get_label() const
  {
    return m_label;
  }

  //! setting gauge configuration.
  void set_config(Field *U)
  {
    m_U = U;

    m_fopr_denom->set_config(m_U);
    m_fopr_numer->set_config(m_U);
    m_fopr_force_MD->set_config(m_U);
  }

  //! Langevin step called at the beginning of HMC.
  double langevin(RandomNumbers *);

  //! calculation of Hamiltonian.
  double calcH();

  //! returns the force for updating conjugate momentum.
  //const Field force();
  void force(Field&);

 private:
};
#endif
