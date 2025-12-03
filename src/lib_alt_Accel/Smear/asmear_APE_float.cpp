/*!
      @file    asmear_APE_float.cpp
      @brief
      @author  Issaku Kanamori (kanamori)
               $LastChangedBy: kanamori $
      @date    $LastChangedDate: 2025-12-03 19:35:35 +0900 (2025年12月03日 (水)) $
      @version $LastChangedRevision: 2668 $
*/

#include "lib_alt_Accel/Smear/asmear_APE.h"

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <vector>
using namespace std;


#include "lib_alt_Accel/inline/define_params.h"

typedef float real_t;

#include "lib_alt_Accel/Field/afield.h"
#include "lib_alt_Accel/Field/afield-inc.h"
#include "lib_alt_Accel/Field/aindex_lex.h"


#ifdef USE_FACTORY_AUTOREGISTER
namespace {
  bool init = ASmear<AField<float,ACCEL> >::Factory_params::Register(
                                   "APE", create_object_with_params);
}
#endif


template<>
const std::string ASmear_APE<AField<float,ACCEL> >::class_name
                       = "ASmear_APE<AField<float,ACCEL> >";


#include "lib_alt_Accel/Smear/asmear_APE-tmpl.h"

// class instanciation.
template class ASmear_APE<AField<float,ACCEL> >;

//============================================================END=====
