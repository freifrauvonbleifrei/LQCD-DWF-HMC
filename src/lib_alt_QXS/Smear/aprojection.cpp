/*
        @file    aprojection.cpp
        @brief
        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: kanamori $
        @date    $LastChangedDate:: 2025-12-03 19:35:35 #$
        @version $LastChangedRevision: 2668 $
*/

#include "lib/Smear/aprojection.h"
#include "lib_alt_QXS/Field/afield.h"

#ifdef USE_FACTORY

#ifdef USE_FACTORY_AUTOREGISTER
#else

#include "lib_alt_QXS/Smear/aprojection_Stout_SU3.h"

template<>
bool AProjection<AField<double,QXS> >::init_factory()
{
  typedef AField<double,QXS> AFIELD;

  bool result = true;

  result &= AProjection_Stout_SU3<AFIELD>::register_factory();

  return result;
}

#endif /* USE_FACTORY_AUTOREGISTER */

#endif /* USE_FACTORY */


// explicit instanciation.
template class AProjection<AField<double,QXS> >;

//============================================================END=====
