/*!
        @file    aforce_F_Ratio_eo.cpp
        @brief
        @author  Issaku Kanamori (kanamori)
                 $LastChangedBy: kanamori $
        @date    $LastChangedDate: 2025-12-03 19:35:35 +0900 (2025年12月03日 (水)) $
        @version $LastChangedRevision: 2668 $
*/

#include "lib_alt/Force/Fermion/aforce_F_Ratio_eo.h"

#include<cassert>

// include files in core library
#include "lib/ResourceManager/threadManager.h"
#include "lib/IO/bridgeIO.h"
using Bridge::vout;

// include files in alt-code dorectories
#include "lib_alt_Accel/Field/afield.h"
#include "lib_alt_Accel/Field/aindex_lex.h"
#include "lib_alt_Accel/Field/afield-inc.h"

#include "lib_alt/Force/Fermion/aforce_F_Ratio_eo-tmpl.h"


// explicit instanciation for AField<double,ACCEL>.
template<>
const std::string AForce_F_Ratio_eo<AField<double,ACCEL> >::class_name
                           = "AForce_F_Ratio_eo<AField<double,ACCEL> >";


template class AForce_F_Ratio_eo<AField<double,ACCEL> >;

//============================================================END=====
