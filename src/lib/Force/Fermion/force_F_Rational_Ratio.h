/*!
        @file    force_F_Rational_Ratio.h

        @brief

        @author  Hideo Matsufuru (matsufuru)
                 $LastChangedBy: kanamori $

        @date    $LastChangedDate:: 2025-12-03 19:35:35 #$

        @version $LastChangedRevision: 2668 $
*/

#ifndef FORCE_F_RATIONAL_RATIO_INCLUDED
#define FORCE_F_RATIONAL_RATIO_INCLUDED

#include "force_F.h"

#include "Field/field_F.h"
#include "lib/Solver/shiftsolver.h"

#include "IO/bridgeIO.h"
using Bridge::vout;

//! Force calculation for smeared fermion operators.

/*!
  Force for the following action:
  S = phi^dag (A^dag A)^{1/4} (B^dag B)^{-1/2} (A^dag A)^{1/4} \phi
    A: numerator, the rational parameter must be exponet = -1/2
    B: denominator, the rational parameter must be exponent = 1/4

                                    [ 1 Jul 2023 I.Kanamori]
 */


class Force_F_Rational_Ratio : public Force
{
 public:
  static const std::string class_name;

 private:
  Bridge::VerboseLevel m_vl;

  // for denominator operator
  int m_denom_Np;                  // number of poles in rational approx.
  int m_denom_n_exp, m_denom_d_exp;    // numerator and denominator of the exponent
  double m_denom_x_min, m_denom_x_max; // valid range of approximate sign function
  int m_denom_Niter;               // max iteration of shiftsolver
  double m_denom_Stop_cond;        // stopping condition of shift solver
  std::string m_denom_str_vlevel;  // verbose level of the shift solver
  Fopr *m_fopr_denom;
  Force *m_force_denom;

  // for numerator operator
  int m_numer_Np;                  // number of poles in rational approx.
  int m_numer_n_exp, m_numer_d_exp;    // numerator and denominator of the exponent
  double m_numer_x_min, m_numer_x_max; // valid range of approximate sign function
  int m_numer_Niter;               // max iteration of shiftsolver
  double m_numer_Stop_cond;        // stopping condition of shift solver
  std::string m_numer_str_vlevel;  // verbose level of the shift solver
  Fopr *m_fopr_numer;
  Force *m_force_numer;


  Field_G *m_U;

  // rational approx. coefficients
  double m_denom_a0;
  std::vector<double> m_denom_bl;
  std::vector<double> m_denom_cl;

  double m_numer_a0;
  std::vector<double> m_numer_bl;
  std::vector<double> m_numer_cl;

  // shiftsolvers: this class has the instances
  unique_ptr<Shiftsolver> m_shiftsolver_denom;
  unique_ptr<Shiftsolver> m_shiftsolver_numer;

public:
  Force_F_Rational_Ratio(Fopr *fopr_numer, Force *force_numer,
                         Fopr *fopr_denom, Force *force_denom)
    : m_vl(CommonParameters::Vlevel()),
      m_fopr_numer(fopr_numer), m_force_numer(force_numer),
      m_fopr_denom(fopr_denom), m_force_denom(force_denom)
  {}

  Force_F_Rational_Ratio(Fopr *fopr_numer, Force *force_numer,
                         Fopr *fopr_denom, Force *force_denom,
                         const Parameters& params)
    : m_vl(CommonParameters::Vlevel()),
      m_fopr_numer(fopr_numer), m_force_numer(force_numer),
      m_fopr_denom(fopr_denom), m_force_denom(force_denom)
  {
    set_parameters(params);
  }

  ~Force_F_Rational_Ratio() {}

  void set_parameters(const Parameters& params);
  void set_parameters(const int denom_Np, const int denom_n_exp, const int denom_d_exp,
                      const double denom_x_min, const double denom_x_max,
                      const int denom_Niter, const double denom_Stop_cond,
                      const std::string str_denom_vlevel,
                      const int numer_Np, const int numer_n_exp, const int numer_d_exp,
                      const double numer_x_min, const double numer_x_max,
                      const int numer_Niter, const double numer_Stop_cond,
                      const std::string str_numer_vlevel);

  void get_parameters(Parameters& params) const;

  void set_config(Field *U)
  {
    m_U = (Field_G *)U;
    m_fopr_denom->set_config(m_U);
    m_fopr_numer->set_config(m_U);
    m_force_denom->set_config(m_U);
    m_force_numer->set_config(m_U);
  }

  void force_udiv(Field&, const Field&);

  void force_core1(Field&, const Field&, const Field&);  // dummy entry
  void force_udiv1(Field&, const Field&, const Field&);  // dummy entry

 private:
  void force_udiv_impl(Field_G&, const Field_F&);

  void set_rational_parameters(
          double &a0, std::vector<double> &bl, std::vector<double> &cl,
          const int Np, const int n_exp, const int d_exp,
          const double x_min, const double x_max);

};
#endif
