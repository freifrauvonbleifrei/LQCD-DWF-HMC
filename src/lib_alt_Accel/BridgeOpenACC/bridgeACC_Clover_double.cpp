/*!
      @file    bridgeACC_Clover_double.cpp
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

#define SU3_3RD_ROW_RECONST

#include "src/mult_Clover_openacc-inc.h"
#include "src/mult_CloverTerm_openacc-inc.h"

}

//============================================================END=====
