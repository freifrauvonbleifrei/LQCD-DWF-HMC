/*!
      @file    bridgeACC_AField_double.cpp
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

#include "src/afield_openacc-inc.h"
#include "src/afield_Gauge_openacc-inc.h"
#include "src/index_eo_alt_openacc-inc.h"
#include "src/shiftAField_lex_openacc-inc.h"

}

//============================================================END=====
