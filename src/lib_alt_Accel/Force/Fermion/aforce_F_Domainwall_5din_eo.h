/*!
        @file    aforce_F_Domainwall_5din_eo.h

        @brief

        @author  Issaku Kanamori (kanamori)
                 $LastChangedBy: kanamori $

        @date    $LastChangedDate: 2025-12-03 19:35:35 +0900 (2025年12月03日 (水)) $

        @version $LastChangedRevision: 2668 $
*/

#ifndef AFORCE_F_DOMAINWALL_5DIN_EO_INCLUDED
#define AFORCE_F_DOMAINWALL_5DIN_EO_INCLUDED

#include "lib/Force/Fermion/aforce_F.h"
#include "lib/Fopr/afopr.h"
#include "lib/IO/bridgeIO.h"
using Bridge::vout;

template<typename REALTYPE, Impl IMPL> class AIndex_eo;

//! Force calculation for domain-wall fermions.

#include "lib_alt_Accel/Fopr/afopr_Domainwall_5din_eo.h"
#include "lib_alt_Accel/Force/Fermion/aforce_F_Domainwall_5din.h"


template<typename AFIELD>
class AForce_F_Domainwall_5din_eo : public AForce_F<AFIELD>
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

  unique_ptr<AFopr_Domainwall_5din_eo<AFIELD> >  m_fopr_dw;
  unique_ptr<AForce_F_Domainwall_5din<AFIELD> >  m_force_lex;

  AIndex_eo<real_t, AFIELD::IMPL>  m_index_eo;

  AFIELD m_v1, m_v2, m_v3;    //!< working vector
  AFIELD m_z1;                //!< working vector
  AFIELD m_force1;            //!< working gauge field
  AFIELD m_eta_lex, m_zeta_lex;  //!< lexical vector

  int m_Nvcd;
  int m_NinF;

public:

  AForce_F_Domainwall_5din_eo(const Parameters& params)
    : AForce_F<AFIELD>()
  { init(params); }

  ~AForce_F_Domainwall_5din_eo() { tidyup(); }

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
  AForce_F_Domainwall_5din_eo();

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
