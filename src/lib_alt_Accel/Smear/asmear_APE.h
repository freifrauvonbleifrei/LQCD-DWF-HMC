/*!
        @file    asmear_APE.h

        @brief

        @author  Issaku Kanamori (kanamori)
                 $LastChangedBy: kanamori $

        @date    $LastChangedDate: 2025-12-03 19:35:35 +0900 (2025年12月03日 (水)) $

        @version $LastChangedRevision: 2668 $
*/

#ifndef ASMEAR_APE_INCLUDED
#define ASMEAR_APE_INCLUDED

#include "lib/Smear/asmear.h"

#include "lib/Measurements/Gauge/staple_lex.h"
#include "lib/IO/bridgeIO.h"
using Bridge::vout;

#include "lib_alt_Accel/Measurements/Gauge/astaple_lex.h"

//! APE type smearing of link variables.
/* taken from QXS version
    add real_t
   to do: move to lib_alt/Smear
 */


template<typename AFIELD>
class ASmear_APE : public ASmear<AFIELD>
{
 public:
  static const std::string class_name;
  typedef typename AFIELD::real_t real_t;

 private:
  Bridge::VerboseLevel m_vl;

  int m_Nvol, m_Ndim;
  std::valarray<double> m_rho;   //!< smearing parameter
  AProjection<AFIELD> *m_proj;   //!< projector to group element.
  AStaple_lex<AFIELD> *m_staple;

  AFIELD m_v1, m_v2, m_v3;
  AFIELD m_w1, m_w2;

 public:
  //! Constructor.
  ASmear_APE(AProjection<AFIELD> *proj, const Parameters& params)
    { init(proj, params); }

  //! Deconstructor
  ~ASmear_APE() { tidyup(); }

  //! Setting parameters with Parameters object.
  void set_parameters(const Parameters& params);

  //! Setting parameter with isotropic parameter.
  void set_parameters(const double rho1);

  //! Setting parameter with anisotropic parameter.
  void set_parameters(const std::vector<double>& rho);

  //! Getting parameters by Parameters object.
  void get_parameters(Parameters& params) const;

  //! Smearing of a given gauge field.
  void smear(Field_G& Usmear, const Field_G& U);

  //! Smearing of a given gauge field.
  void smear(AFIELD& Usmear, const AFIELD& U);

 private:

  //! Initial setup.
  void init(AProjection<AFIELD>* proj, const Parameters& params);

  //! Final clean-up.
  void tidyup();


#ifdef USE_FACTORY
 private:
  static ASmear<AFIELD> *create_object_with_params(
                                         AProjection<AFIELD> *proj,
                                         const Parameters& params)
  { return new ASmear_APE<AFIELD>(proj, params); }

 public:
  static bool register_factory()
  {
    bool init = true;
    init &= ASmear<AFIELD>::Factory_params::Register(
                                 "APE", create_object_with_params);
    return init;
  }
#endif
};
#endif
