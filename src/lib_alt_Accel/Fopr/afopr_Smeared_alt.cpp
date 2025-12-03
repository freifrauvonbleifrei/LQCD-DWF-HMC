/*!
      @file    afopr_Smeared_alt.cpp
      @brief
      @author  Issaku Kanamori (kanamori)
               $LastChangedBy: kanamori $
      @date    $LastChangedDate: 2025-12-03 19:35:35 +0900 (2025年12月03日 (水)) $
      @version $LastChangedRevision: 2668 $
*/

#include "lib/Fopr/afopr_Smeared_alt.h"
#include "lib/Fopr/afopr_Smeared_alt-tmpl.h"

#include "lib_alt_Accel/Field/afield.h"


template<>
const std::string AFopr_Smeared_alt<AField<double, ACCEL> >::class_name
  = "AFopr_Smeared_alt<AField<double,ACCEL> >";

template<>
const std::string AFopr_Smeared_alt<AField<float, ACCEL> >::class_name
  = "AFopr_Smeared_alt<AField<float,ACCEL> >";

#ifdef USE_FACTORY_AUTOREGISTER
namespace {
  bool init1 = AFopr<AField<double, ACCEL> >::Factory_params::Register(
    "Smeared_alt", create_object);

  bool init2 = AFopr<AField<float, ACCEL> >::Factory_params::Register(
    "Smeared_alt", create_object);
}
#endif

// explicit instanciation.
template class AFopr_Smeared_alt<AField<double, ACCEL> >;
template class AFopr_Smeared_alt<AField<float, ACCEL> >;

//============================================================END=====
