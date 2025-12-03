/*!
        @file    aforce_F_Rational_Ratio.cpp
        @brief
        @author  Issaku Kanamori (kanamori)
                 $LastChangedBy: kanamori $
        @date    $LastChangedDate:: 2025-12-03 19:35:35 #$
        @version $LastChangedRevision: 2668 $
*/

#include "lib_alt/Force/Fermion/aforce_F_Rational_Ratio.h"

#include<cassert>

// include files in core library
#include "lib/ResourceManager/threadManager.h"
#include "lib/IO/bridgeIO.h"
using Bridge::vout;

// include files in alt-code dorectories
#include "lib_alt_QXS/Field/afield.h"
#include "lib_alt_QXS/Field/aindex_lex.h"
#include "lib_alt_QXS/Field/afield-inc.h"

#include "lib_alt/Force/Fermion/aforce_F_Rational_Ratio-tmpl.h"


// explicit instanciation for AField<double,QXS>.
template<>
const std::string AForce_F_Rational_Ratio<AField<double,QXS> >::class_name
                           = "AForce_F_Rational_Ratio<AField<double,QXS> >";


template class AForce_F_Rational_Ratio<AField<double,QXS> >;

//============================================================END=====
