/*!
        @file    aforce_F_Domainwall_5din.h

        @brief

        @author  Issaku Kanamori (kanamori)
                 $LastChangedBy: kanamori $

        @date    $LastChangedDate: 2025-12-03 19:35:35 +0900 (2025年12月03日 (水)) $

        @version $LastChangedRevision: 2668 $
*/

#ifndef AFORCE_F_DOMAINWALL_5DIN_INCLUDED
#define AFORCE_F_DOMAINWALL_5DIN_INCLUDED

#include "lib/Force/Fermion/aforce_F.h"
#include "lib/Fopr/afopr.h"
#include "lib/IO/bridgeIO.h"
using Bridge::vout;


template<typename AFIELD> class ShiftAField_lex;

#include "lib_alt_Accel/Fopr/afopr_Domainwall_5din.h"



//! Force calculation for domain-wall fermions.



template<typename AFIELD>
class AForce_F_Domainwall_5din : public AForce_F<AFIELD>
{
 public:
  typedef typename AFIELD::real_t real_t;
  static const std::string class_name;

 private:
  // parameters common to overlap fermion
  real_t m_mq;                  //!< quark mass
  real_t m_M0;                  //!< domain-wall height
  int m_Ns;                     //!< size of fifth-dimension
  std::vector<int> m_boundary;  //!< boundary condition
  std::vector<real_t> m_b;      //!< coefficient b (array)
  std::vector<real_t> m_c;      //!< coefficient c (array)
  std::string m_repr;           //!< gamma matrix representation

  Bridge::VerboseLevel m_vl;

  std::string m_mode;

  using AForce_F<AFIELD>::m_Ucp;
  Field *m_U;

  unique_ptr<AFopr_Domainwall_5din<AFIELD> > m_fopr_dw;
  unique_ptr<ShiftAField_lex<AFIELD> > m_shift;

  AFIELD m_v1, m_v2;          //!< working vector
  AFIELD m_h1, m_h2;          //!< working vector (half spinor)
  AFIELD m_z1, m_z2;          //!< working vector
  AFIELD m_force1, m_force2;  //!< working gauge field

  int m_Nvcd;
  int m_NinF;

public:

  AForce_F_Domainwall_5din(const Parameters& params)
    : AForce_F<AFIELD>()
  { init(params); }

  ~AForce_F_Domainwall_5din() { tidyup(); }

  void set_parameters(const Parameters& params);

  void set_parameters(const real_t mq,
                      const real_t M0,
                      const int Ns,
                      const std::vector<int> bc,
                      const real_t b,
                      const real_t c);

  void get_parameters(Parameters& params) const;

  void set_config(Field *U);

  void set_mode(const std::string& mode);

  void force_udiv(AFIELD& force, const AFIELD& eta);

  void force_udiv1(AFIELD& force, const AFIELD& zeta, const AFIELD& eta);

 private:
  void init(const Parameters& params);

  void tidyup();

  void set_config_omp(Field *U);

  void set_config_impl(Field *U);

  void force_udiv1_H(AFIELD& force, const AFIELD& zeta, const AFIELD& eta);
  void force_udiv1_Hdag(AFIELD& force, const AFIELD& zeta, const AFIELD& eta);

  void force_udiv1_D(AFIELD& force, const AFIELD& zeta, const AFIELD& eta);
  void force_udiv1_Ddag(AFIELD& force, const AFIELD& zeta, const AFIELD& eta);

};

#endif
