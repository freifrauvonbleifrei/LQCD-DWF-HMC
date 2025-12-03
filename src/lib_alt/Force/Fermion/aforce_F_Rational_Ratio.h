/*!
        @file    aforce_F_Rational_Ratio.h

        @brief

        @author  Issaku Kanamori (kanamori)
                 $LastChangedBy: kanamori $

        @date    $LastChangedDate:: 2025-12-03 19:35:35 #$

        @version $LastChangedRevision: 2668 $
*/

#ifndef AFORCE_F_RATIONAL_RATIO_INCLUDED
#define AFORCE_F_RATIONAL_RATIO_INCLUDED

#include "Force/Fermion/aforce_F.h"

#include "Field/field_F.h"
#include "Fopr/afopr.h"
#include "Solver/ashiftsolver_CG.h"
#include "Tools/math_Rational.h"

#include "IO/bridgeIO.h"
using Bridge::vout;

//! Force calculation for smeared fermion operators.

/*!
  Force for the following action:
  S = phi^dag (A^dag A)^{1/4} (B^dag B)^{-1/2} (A^dag A)^{1/4} \phi
    A: numerator (fopr1),   assumes the rational parameter with exponet = -1/2
    B: denominator (fopr2), assumes the rational parameter with exponent = 1/4

                                    [ 15 Aug 2023 I.Kanamori]
 */


template<typename AFIELD>
class AForce_F_Rational_Ratio : public AForce_F<AFIELD>
{
 public:
  static const std::string class_name;
  typedef typename AFIELD::real_t real_t;

 private:
  // for fopr1 (numerator)
  int    m_Np1;               // number of poles in rational approx.
  int    m_n_exp1, m_d_exp1;  // numerator and denominator of the exponent
  real_t m_x_min1, m_x_max1;  // valid range of approximate sign function
  int    m_Niter1;            // max iteration of shiftsolver
  real_t m_Stop_cond1;        // stopping condition of shift solver
  AFopr<AFIELD>    *m_fopr1;  // kernel operataor A for the rational approx.
  AForce_F<AFIELD> *m_force1; // force of the kenel operator
  std::string m_str_vlevel1;  // verbose level of the shift solver

  // for fopr2 (denominator)
  int    m_Np2;               // number of poles in rational approx.
  int    m_n_exp2, m_d_exp2;  // numerator and denominator of the exponent
  real_t m_x_min2, m_x_max2;  // valid range of approximate sign function
  int    m_Niter2;            // max iteration of shiftsolver
  real_t m_Stop_cond2;        // stopping condition of shift solver
  AFopr<AFIELD>    *m_fopr2;  // kernel operataor B for the rational approx.
  AForce_F<AFIELD> *m_force2; // force of the kenel operator
  std::string m_str_vlevel2;  // verbose level of the shift solver

  Bridge::VerboseLevel m_vl;

  Field_G *m_U;

  // rational approx. coefficients
  real_t  m_a01;
  std::vector<real_t> m_bl1;
  std::vector<real_t> m_cl1;

  real_t  m_a02;
  std::vector<real_t> m_bl2;
  std::vector<real_t> m_cl2;

  std::vector<AFIELD> m_eta1;
  std::vector<AFIELD> m_eta2;
  std::vector<AFIELD> m_zeta1;
  AFIELD m_v1, m_v2;

  unique_ptr< AShiftsolver_CG<AFIELD, AFopr<AFIELD> > > m_solver1;
  unique_ptr< AShiftsolver_CG<AFIELD, AFopr<AFIELD> > > m_solver2;

  AFIELD m_force;

 public:

  AForce_F_Rational_Ratio(AFopr<AFIELD> *fopr1, AForce_F<AFIELD> *force1,
                          AFopr<AFIELD> *fopr2, AForce_F<AFIELD> *force2,
                          const Parameters& params)
    : AForce_F<AFIELD>(),
    m_fopr1(fopr1), m_force1(force1), m_fopr2(fopr2), m_force2(force2)
  { init(params); }

  ~AForce_F_Rational_Ratio() { tidyup(); }

  void set_parameters(const Parameters& params);
  void set_parameters(const int Np1, const int n_exp1, const int d_exp1,
                      const real_t x_min1, const real_t x_max1,
                      const int Niter1, const real_t Stop_cond1,
                      const std::string str_vlevel1,
                      const int Np2, const int n_exp2, const int d_exp2,
                      const real_t x_min2, const real_t x_max2,
                      const int Niter2, const real_t Stop_cond2,
                      const std::string str_vlevel2);


  void get_parameters(Parameters& params) const;

  void set_config(Field *U);

  void force_udiv(AFIELD&, const AFIELD&);

  void force_core1(AFIELD&, const AFIELD&, const AFIELD&);
  void force_udiv1(AFIELD&, const AFIELD&, const AFIELD&);

 private:

  void init(const Parameters& params);

  void tidyup();

  void force_udiv_impl(AFIELD&, const AFIELD&);
  void force_udiv1_impl(AFIELD&, const AFIELD&, const AFIELD&);

  void set_rational_parameters(
          double &a0, std::vector<double> &bl, std::vector<double> &cl,
          const int Np, const int n_exp, const int d_exp,
          const double x_min, const double x_max);

#ifdef USE_FACTORY
 private:
  //static AForce_F<AFIELD> *create_object(AFopr<AFIELD> *fopr1,
  //                                       AForce_F<AFIELD> *force_F1,
  //                                       AFopr<AFIELD> *fopr2,
  //                                       AForce_F<AFIELD> *force_F2)
  //{ return new AForce_F_Rational_Ratio(fopr1, force_F1, fopr2, force_F2); }

 public:
  //static bool register_factory()
  //{
  //  bool init1 = AForce_F<AFIELD>::Factory_fopr_force::Register(
  //                                         "Rational_Ratio", create_object);
  //  return init1;
  // }
#endif


};
#endif
