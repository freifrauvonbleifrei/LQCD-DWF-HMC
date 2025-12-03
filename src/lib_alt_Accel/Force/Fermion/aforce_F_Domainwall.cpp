/*!
        @file    aforce_F_Domainwall.cpp
        @brief
        @author  Issaku Kanamori (kanamori)
                 $LastChangedBy: kanamori $
        @date    $LastChangedDate: 2025-12-03 19:35:35 +0900 (2025年12月03日 (水)) $
        @version $LastChangedRevision: 2668 $
*/

#include "lib_alt/Force/Fermion/aforce_F_Domainwall.h"

#include<cassert>

// include files in core library
#include "lib/ResourceManager/threadManager.h"
#include "lib/IO/bridgeIO.h"
using Bridge::vout;

typedef double real_t;

#include "lib_alt_Accel/Field/afield.h"
#include "lib_alt_Accel/Field/aindex_lex.h"
#include "lib_alt_Accel/Field/afield-inc.h"
#include "lib_alt_Accel/Field/shiftAField_lex.h"
#include "lib_alt_Accel/Fopr/afopr_Wilson.h"


#include "lib_alt/Force/Fermion/aforce_F_Domainwall-tmpl.h"

// explicit instanciation for AField<double,QXS>.
template<>
const std::string AForce_F_Domainwall<AField<double,ACCEL> >::class_name
                         = "AForce_F_Domainwall<AField<double,ACCEL> >";


template class AForce_F_Domainwall<AField<double,ACCEL> >;

//============================================================END=====
