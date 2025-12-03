/*!
        @file    action_F_Rational_eo.h

        @brief

        @author  Issaku Kanamori (kanamori)
                 $LastChangedBy: kanamori $

        @date    $LastChangedDate:: 2025-12-03 19:35:35 #$

        @version $LastChangedRevision: 2668 $
*/

#ifndef ACTION_F_RATIONAL_EO_INCLUDED
#define ACTION_F_RATIONAL_EO_INCLUDED

#include "Action/action.h"

#include "Solver/shiftsolver_CG.h"
#include "Tools/math_Rational.h"
#include "Force/Fermion/force_F_Rational.h"

#include "IO/bridgeIO.h"
using Bridge::vout;

//! action class for RHMC, with externally constructed Fopr_Rational.

/*!
    For the class, Fopr and Force objects are instantiated outside
    the class and specified at the construction.
    This class just provides the framework of rational actions.
                                        [28 Dec 2011 H.Matsufuru]
    unique_ptr is introduced to avoid memory leaks
                                        [21 Mar 2015 Y.Namekawa]
 */


class Action_F_Rational_eo : public Action
{
 public:
  static const std::string class_name;

 private:
  Bridge::VerboseLevel m_vl;

  std::string m_label;    // label of action

  Fopr *m_fopr_eo;
  Force *m_fopr_force_MD;

  unique_ptr<Fopr_Rational> m_fopr_H;
  unique_ptr<Fopr_Rational> m_fopr_psf;

  Field *m_U;

  Field m_psf;

  Parameters m_params_rational_H;
  Parameters m_params_rational_psf;

public:
  //! constructor requires pointers to Fopr and Force instances.
  Action_F_Rational_eo(
    Fopr *fopr_eo, Force *fopr_force_MD)
    : m_vl(CommonParameters::Vlevel()),
    m_fopr_eo(fopr_eo), m_fopr_force_MD(fopr_force_MD)
  {}

  Action_F_Rational_eo(Fopr *fopr_eo, Force *fopr_force_MD,
                    const Parameters& params)
    : m_vl(CommonParameters::Vlevel()),
    m_fopr_eo(fopr_eo), m_fopr_force_MD(fopr_force_MD)
  {
    set_parameters(params);
  }

  //! destructor. constructed instances are deconstructed in tydyup().
  ~Action_F_Rational_eo() {}

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
