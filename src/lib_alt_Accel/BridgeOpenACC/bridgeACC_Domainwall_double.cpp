/*!
      @file    bridgeACC_Domainwall_double.cpp
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

#include "src/inc/mult_Wilson_inline_openacc-inc.h"

#include "src/mult_Domainwall_5din_openacc-inc.h"
#include "src/mult_Domainwall_5din_eo_openacc-inc.h"
#include "src/mult_Domainwall_5din_LUinv_openacc-inc.h"

#ifdef USE_DOMAINWALL_5DIN_4D_KERNEL
#include "src/mult_Domainwall_5din_4d_openacc-inc.h"
#include "src/mult_Domainwall_5din_eo_4d_openacc-inc.h"
#endif

#ifdef USE_DOMAINWALL_5DIN_EE_MATINV_KERNEL
#include "src/mult_Domainwall_5din_matinv_openacc-inc.h"
#endif

}

//============================================================END=====
