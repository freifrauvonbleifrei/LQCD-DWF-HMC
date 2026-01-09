/*!
        @file    aforce_F.cpp
        @brief
        @author  Issaku Kanamori (kanamori)
                 $LastChangedBy: kanamori $
        @date    $LastChangedDate: 2026-01-09 16:11:38 +0900 (2026年01月09日 (金)) $
        @version $LastChangedRevision: 2687 $
*/

#include "lib/Force/Fermion/aforce_F.h"
#include "lib/ResourceManager/threadManager.h"


#ifdef USE_FACTORY

#ifdef USE_FACTORY_AUTOREGISTER
#else

//#include "lib_alt/Force/Fermion/aforce_F_Wilson_Nf2.h"
//#include "lib_alt_QXS/Force/Fermion/aforce_F_Clover_Nf2.h"
//#include "lib_alt/Force/Fermion/aforce_F_Smeared.h"
//#include "lib/Force/Fermion/aforce_F_Rational.h"

template<typename AFIELD>
bool AForce_F<AFIELD>::init_factory()
{
  bool result = true;
  // result &= AForce_F_Wilson_Nf2<AFIELD>::register_factory();
  //  result &= AForce_F_Clover_Nf2<AFIELD>::register_factory();
  //  result &= AForce_F_Smeared<AFIELD>::register_factory();
  //  result &= AForce_F_Rational<AFIELD>::register_factory();
  return result;
}

#endif /* USE_FACTORY_AUTOREGISTER */

#endif /* USE_FACTORY */



// include files in alt-code dorectories
#include "lib_alt_Accel/inline/define_params.h"

typedef double real_t;


// include files in alt-code dorectories
#include "lib_alt_Accel/Field/afield.h"
#include "lib_alt_Accel/Field/aindex_lex.h"
#include "lib_alt_Accel/Field/afield-inc.h"
#include "lib_alt_Accel/Field/afield_Gauge-inc.h"
namespace Alt_Gauge = Accel_Gauge;

typedef AField<double, ACCEL> AFIELD;

template<>
const std::string AForce_F<AField<double, ACCEL> >
::class_name = "AForce_F<AField<double, ACCEL> >";


namespace Accel_Gauge {
  void mult_generator(AFIELD& force, AFIELD &ut, AFIELD &Ucp)
{
  int Ndim = force.nex();

#pragma omp barrier
  for(int mu = 0; mu < Ndim; ++mu){
    copy(ut, 0, force, mu);
    mult_Gnn(force, mu, Ucp, mu, ut, 0);
    at_G(force, mu);
  }

  scal(force, -2.0);
#pragma omp barrier
}

}


//====================================================================
template<>
void AForce_F<AFIELD>::init()
{
  int Nc   = CommonParameters::Nc();
  int Nin  = 2 * Nc * Nc;
  int Nvol = CommonParameters::Nvol();
  int Ndim = CommonParameters::Ndim();

  m_Ucp.reset(Nin, Nvol, Ndim);
  m_ut.reset( Nin, Nvol, 1);
}

//====================================================================
template<>
void AForce_F<AFIELD>::tidyup()
{
  // do nothing.
}

//====================================================================
// Note that since mult_generator() is called in force_core() and
// force_core1(), its specialization must be placed before the
// latter functions.
template<>
void AForce_F<AFIELD>::mult_generator(AFIELD& force)
{
  Accel_Gauge::mult_generator(force, m_ut, m_Ucp);
}

template<>
void AForce_F<AFIELD>::mult_generator(Field_G& force)
{
  vout.crucial("%s: AForce_F<AFIELD>::mult_generator(Field_G& force) is not implemented\n", class_name.c_str());
#pragma omp barrier
  abort();


}


//====================================================================
template<>
void AForce_F<AFIELD>::force_core(AFIELD& force, const AFIELD& eta)
{
  force_udiv(force, eta);
  mult_generator(force);
  //  Accel_Gauge::mult_generator(force, m_ut, m_Ucp);
}

//====================================================================
template<>
void AForce_F<AFIELD>::force_core1(AFIELD& force,
                                   const AFIELD& zeta,
                                   const AFIELD& eta)
{
  force_udiv1(force, zeta, eta);
  mult_generator(force);
  //  Accel_Gauge::mult_generator(force, m_ut, m_Ucp);
}

//====================================================================

// explicit instanciation.

template class AForce_F<AField<double,ACCEL> >;
//template class AForce_F<AField<float,ACCEL> >;

//============================================================END=====
