/*!
        @file    aforce_F_Domainwall.cpp
        @brief
        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: kanamori $
        @date    $LastChangedDate:: 2025-12-03 19:35:35 #$
        @version $LastChangedRevision: 2668 $
*/


#include<cassert>

// include files in core library
#include "lib/ResourceManager/threadManager.h"
#include "lib/IO/bridgeIO.h"
using Bridge::vout;

#include "lib_alt_QXS/inline/define_vlen.h"
#include "lib_alt_QXS/inline/define_params.h"

#define  VLEN     VLEND
#define  VLENX    VLENXD
#define  VLENY    VLENYD

typedef double real_t;

#include "lib_alt_QXS/inline/vsimd_double-inc.h"
#include "lib_alt_QXS/inline/vsimd_common_double-inc.h"
#include "lib_alt_QXS/inline/vsimd_Wilson_SU3_double-inc.h"
#include "lib_alt_QXS/inline/vsimd_Domainwall_SU3_double-inc.h"
#include "lib_alt_QXS/BridgeQXS/bridgeQXS_Domainwall.h"


#include "lib_alt_QXS/Field/afield.h"
#include "lib_alt_QXS/Field/afield-inc.h"
#include "lib_alt_QXS/Field/afield_Spinor-inc.h"
#include "lib_alt_QXS/Field/shiftAField_lex.h"
#include "lib_alt_QXS/Fopr/afopr_Domainwall_5din_eo.h"
#include "aforce_F_Domainwall_5din.h"
#include "aforce_F_Domainwall_5din_eo.h"


#include "aforce_F_Domainwall_5din_eo-tmpl.h"

// explicit instanciation for AField<double,QXS>.
template<>
const std::string AForce_F_Domainwall_5din_eo<AField<double,QXS> >::class_name
                         = "AForce_F_Domainwall_5din_eo<AField<double,QXS> >";


template class AForce_F_Domainwall_5din_eo<AField<double,QXS> >;

//============================================================END=====
