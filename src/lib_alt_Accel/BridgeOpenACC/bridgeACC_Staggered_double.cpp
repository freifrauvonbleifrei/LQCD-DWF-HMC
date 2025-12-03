/*!
      @file    bridgeACC_Staggered_double.cpp
      @brief
      @author  Hideo Matsufuru (matufuru)
      @date    $LastChangedDate: 2013-01-22 13:51:53 #$
      @version $LastChangedRevision: 2653 $
*/

#include "bridge_defs.h"

#include "lib_alt_Accel/inline/define_params.h"
#include "lib_alt_Accel/inline/define_index.h"

typedef double real_t;

namespace BridgeACC {

#include "src/mult_Staggered_openacc-inc.h"
#include "src/mult_Staggered_dir_openacc-inc.h"
#include "src/mult_Staggered_eo_openacc-inc.h"
#include "src/mult_Staggered_eo_dir_openacc-inc.h"

}

//============================================================END=====
