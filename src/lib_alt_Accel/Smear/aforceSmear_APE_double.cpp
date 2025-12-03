/*!
      @file    aforceSmear_APE_double.cpp
      @brief
      @author  Issaku Kanamori (kanamori)
               $LastChangedBy: kanamori $
      @date    $LastChangedDate: 2025-12-03 19:35:35 +0900 (2025年12月03日 (水)) $
      @version $LastChangedRevision: 2668 $
*/

#include "lib_alt_Accel/Smear/aforceSmear_APE.h"

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <vector>
using namespace std;


#include "lib_alt_Accel/inline/define_params.h"

typedef double real_t;

#include "lib_alt_Accel/Field/afield.h"
#include "lib_alt_Accel/Field/afield-inc.h"
#include "lib_alt_Accel/Field/aindex_lex.h"
#include "lib_alt_Accel/Field/afield_Gauge-inc.h"
namespace Alt_Gauge = Accel_Gauge;

#ifdef USE_FACTORY_AUTOREGISTER
namespace {
  bool init1 = AProjection<AField<double,ACCEL> >::Factory_params::Register(
                    "Projection_Stout_SU3", create_object_with_params);
}
#endif


template<>
const std::string AForceSmear_APE<AField<double,ACCEL> >::class_name
                       = "AforceSmear_APE<AField<double,ACCEL> >";


#include "lib_alt_Accel/Smear/aforceSmear_APE-tmpl.h"

// class instanciation.
template class AForceSmear_APE<AField<double,ACCEL> >;

//============================================================END=====
