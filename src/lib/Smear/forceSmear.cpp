/*!
        @file    forceSmear.cpp

        @brief

        @author  Hideo Matsufuru (matsufuru)
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate:: 2025-09-02 15:10:15 #$

        @version $LastChangedRevision: 2654 $
*/

#include "forceSmear.h"

#ifdef USE_FACTORY

#ifdef USE_FACTORY_AUTOREGISTER
#else
// setup factories for all subclasses

#include "forceSmear_APE.h"
#include "forceSmear_HYP.h"
#include "forceSmear_APE_SF.h"
#include "forceSmear_HYP_SF.h"

template<>
bool ForceSmear::init_factory()
{
  bool result = true;

  result &= ForceSmear_APE::register_factory();
  result &= ForceSmear_HYP::register_factory();
  result &= ForceSmear_APE_SF::register_factory();
  result &= ForceSmear_HYP_SF::register_factory();

  return result;
}


#endif /* USE_FACTORY_AUTOREGISTER */

#endif /* USE_FACTORY */
