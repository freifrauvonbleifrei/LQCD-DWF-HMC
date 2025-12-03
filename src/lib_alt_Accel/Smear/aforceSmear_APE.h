/*!
        @file    aforceSmear_APE.h

        @brief

        @author  Issaku Kanamori (kanamori)
                 $LastChangedBy: kanamori $

        @date    $LastChangedDate: 2025-12-03 19:35:35 +0900 (2025年12月03日 (水)) $

        @version $LastChangedRevision: 2668 $
*/

#ifndef ACCEL_AFORCESMEAR_APE_INCLUDED
#define ACCEL_AFORCESMEAR_APE_INCLUDED

#include "lib/Smear/aforceSmear.h"

#include "lib_alt_Accel/Field/shiftAField_lex.h"
#include "lib/Field/shiftField_lex.h"

#include "lib/Smear/aprojection.h"
#include "lib/Smear/smear_APE.h"

#include "lib/IO/bridgeIO.h"
using Bridge::vout;

//! Recursive calculation for APE smeared fermion force.

/*!
                                [08 Apr 2012 H.Matsufuru]
    (Coding history will be recovered from trac.)
    YAML is implemented.        [14 Nov 2012 Y.Namekawa]
    unique_ptr is introduced to avoid memory leaks
                                [21 Mar 2015 Y.Namekawa]
    Accel version, which is a copy of the QXS version
                                [21 Nov 2025 I.Kanamori]
*/

// todo: move to lib_alt/Smear for all alternative cases

template<typename AFIELD>
class AForceSmear_APE : public AForceSmear<AFIELD>
{
 public:
  static const std::string class_name;
  typedef typename AFIELD::real_t real_t;

 private:
  Bridge::VerboseLevel m_vl;

  int m_Ndim, m_Nvol;
  std::vector<double> m_rho;
  AProjection<AFIELD>* m_proj;
  ShiftField_lex* m_shift;
  std::vector<Field_G> m_U;
  std::vector<Field_G> m_iTheta;

  Field_G m_vt1, m_vt2, m_vt3;
  Field_G m_C, m_Ctmp, m_sigmap_tmp, m_Xi, m_sigma_tmp;

  std::vector<AFIELD> m_U_alt;
  std::vector<AFIELD> m_iTheta_alt;
  AFIELD m_C_alt;
  unique_ptr< ShiftAField_lex<AFIELD> > m_shift_alt;

  AFIELD m_av1;
  AFIELD m_vt1_alt, m_vt2_alt, m_vt3_alt;
  AFIELD m_sigma_tmp_alt, m_sigmap_tmp_alt;
  AFIELD m_Ctmp_alt;

public:
  AForceSmear_APE(AProjection<AFIELD> *proj)
    : m_vl(CommonParameters::Vlevel()), m_proj(proj)
  {
    init();
  }

  AForceSmear_APE(AProjection<AFIELD>* proj, const Parameters& params)
    : m_vl(CommonParameters::Vlevel()), m_proj(proj)
  {
    init();
    set_parameters(params);
  }

  ~AForceSmear_APE(){ tidyup(); }

  // Setting parameters with Parameters object.
  void set_parameters(const Parameters& params);

  // Setting parameters with uniform smearing parameter.
  void set_parameters(const double rho1);

  // Setting parameters with anisotropic smearing parameter.
  void set_parameters(const std::vector<double>& rho);

  // Getting parameters by Parameters object.
  void get_parameters(Parameters& params) const;

  // Force computation.
  void force_udiv(Field_G& Sigma, const Field_G& Sigma_p, const Field_G& U);

  // Force computation.
  void force_udiv(AFIELD& Sigma, const AFIELD& Sigma_p, const Field_G& U);

private:
  void init();

  void tidyup();

  double rho(const int mu, const int nu)
  {
    return m_rho[mu + nu * m_Ndim];
  }

  void force_each(Field_G&,
                  const Field_G&, const Field_G&,
                  const Field_G&, const Field_G&, const int mu, const int nu);


  void force_each(AFIELD&,
                  const AFIELD&, const AFIELD&,
                  const AFIELD&, const AFIELD&, const int mu, const int nu);

  void staple(Field_G&,
              const Field_G&, const Field_G&,
              const int mu, const int nu);

  void staple(AFIELD&,
              const AFIELD&, const AFIELD&,
              const int mu, const int nu);

#ifdef USE_FACTORY
 private:
  static AForceSmear<AFIELD> *create_object(AProjection<AFIELD> *proj)
  {
    return new AForceSmear_APE<AFIELD>(proj);
  }

  static AForceSmear<AFIELD> *create_object_with_params(
                     AProjection<AFIELD>* proj, const Parameters& params)
  {
    return new AForceSmear_APE<AFIELD>(proj, params);
  }

 public:
  static bool register_factory()
  {
    bool init = true;
    init &= AForceSmear<AFIELD>::Factory::Register("APE", create_object);
    init &= AForceSmear<AFIELD>::Factory_params::Register("APE",
                                              create_object_with_params);
    return init;
  }
#endif
};
#endif
