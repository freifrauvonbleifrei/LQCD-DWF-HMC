/*!
        @file    action_F_alt_Rational_Ratio.cpp
        @brief
        @author  Issaku Kanamori (kanamori)
                 $LastChangedBy: kanamori $
        @date    $LastChangedDate: 2025-12-03 19:35:35 +0900 (2025年12月03日 (水)) $
        @version $LastChangedRevision: 2668 $
*/

#include "lib_alt/Action/Fermion/action_F_alt_Rational_Ratio.h"

#include<cassert>

// include files in core library
#include "lib/ResourceManager/threadManager.h"
#include "lib/IO/bridgeIO.h"
using Bridge::vout;

#include "lib_alt_Accel/Field/afield.h"
#include "lib_alt_Accel/Field/aindex_lex.h"
#include "lib_alt_Accel/Field/afield-inc.h"

#include "lib_alt/Action/Fermion/action_F_alt_Rational_Ratio-tmpl.h"

// explicit instanciation for AField<double,ACCEL>.
template<>
const std::string Action_F_alt_Rational_Ratio<AField<double,ACCEL> >::
     class_name = "Action_F_alt_Rational_Ratio<AField<double,ACCEL> >";


template class Action_F_alt_Rational_Ratio<AField<double,ACCEL> >;

//============================================================END=====
