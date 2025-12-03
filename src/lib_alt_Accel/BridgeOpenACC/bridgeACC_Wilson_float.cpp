/*!
      @file    bridgeACC_Wilson_float.cpp
      @brief
      @author  Hideo Matsufuru (matufuru)
      @date    $LastChangedDate: 2013-01-22 13:51:53 #$
      @version $LastChangedRevision: 2653 $
*/

#include "bridge_defs.h"

#include "lib_alt_Accel/inline/define_params.h"
#include "lib_alt_Accel/inline/define_index.h"

typedef float real_t;

namespace BridgeACC {

#define SU3_3RD_ROW_RECONST

#include "src/mult_Wilson_openacc-inc.h"
#include "src/mult_Wilson_eo_openacc-inc.h"
#include "src/mult_Wilson_dir_openacc-inc.h"
#include "src/mult_Wilson_eo_dir_openacc-inc.h"

}

//============================================================END=====
