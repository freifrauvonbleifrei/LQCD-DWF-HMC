/*!
        @file    aforce_F_Domainwall.cpp
        @brief
        @author  Issaku Kanamori (kanamori)
                 $LastChangedBy: kanamori $
        @date    $LastChangedDate: 2025-12-03 19:35:35 +0900 (2025年12月03日 (水)) $
        @version $LastChangedRevision: 2668 $
*/

#include "lib_alt_Accel/inline/define_available.h"
#ifdef ACCEL_FORCE_DOMAINWALL_5DIN_AVAILABLE

#include <cassert>
#include <stdlib.h>
#include <assert.h>
#include <vector>

// include files in core library
#include "lib/ResourceManager/threadManager.h"
#include "lib/IO/bridgeIO.h"
using Bridge::vout;

#include "lib_alt_Accel/inline/define_params.h"

typedef double real_t;

#include "lib_alt_Accel/BridgeACC/bridgeACC_AField.h"
#include "lib_alt_Accel/BridgeACC/bridgeACC_Domainwall.h"


#include "lib_alt_Accel/Field/afield.h"
#include "lib_alt_Accel/Field/aindex_lex.h"
#include "lib_alt_Accel/Field/afield-inc.h"
#include "lib_alt_Accel/Field/afield_Spinor-inc.h"
#include "lib_alt_Accel/Field/shiftAField_lex.h"
#include "lib_alt_Accel/Fopr/afopr_Domainwall_5din.h"


#include "aforce_F_Domainwall_5din-tmpl.h"

// explicit instanciation for AField<double,ACCEL>.
template<>
const std::string AForce_F_Domainwall_5din<AField<double,ACCEL> >::class_name
                         = "AForce_F_Domainwall_5din<AField<double,ACCEL> >";


template class AForce_F_Domainwall_5din<AField<double,ACCEL> >;

#endif
//============================================================END=====
