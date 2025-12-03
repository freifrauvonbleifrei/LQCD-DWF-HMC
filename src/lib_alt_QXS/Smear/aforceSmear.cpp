/*!
        @file    aforceSmear.cpp
        @brief
        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: kanamori $
        @date    $LastChangedDate:: 2025-12-03 19:35:35 #$
        @version $LastChangedRevision: 2668 $
*/

#include "lib/Smear/aforceSmear.h"
#include "lib_alt_QXS/Field/afield.h"

#ifdef USE_FACTORY

#ifdef USE_FACTORY_AUTOREGISTER
#else
// setup factories for all subclasses

#include "lib_alt_QXS/Smear/aforceSmear_APE.h"

template<>
bool AForceSmear<AField<double,QXS> >::init_factory()
{
  typedef AField<double,QXS> AFIELD;

  bool result = true;

  result &= AForceSmear_APE<AFIELD>::register_factory();

  return result;
}


#endif /* USE_FACTORY_AUTOREGISTER */

#endif /* USE_FACTORY */

// explicit instanciation.
template class AForceSmear<AField<double,QXS> >;

//============================================================END=====
