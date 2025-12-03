/*!
      @file    director_alt_Smear.cpp
      @brief
      @author  Issaku Kanamori
               $LastChangedBy: kanamori $
      @date    $LastChangedDate: 2025-12-03 19:35:35 +0900 (2025年12月03日 (水)) $
      @version $LastChangedRevision: 2668 $
*/

#include "lib_alt/Smear/director_alt_Smear.h"

#include "lib/ResourceManager/threadManager.h"
#include "lib/Fopr/afopr_Smeared.h"

#include "lib_alt_Accel/Field/afield.h"
#include "lib_alt_Accel/Field/afield-inc.h"
#include "lib_alt_Accel/Field/aindex_lex.h"

#include "lib_alt/Smear/director_alt_Smear-tmpl.h"

//====================================================================
// explicit instanciation for AField<double,ACCEL>.
template<>
const std::string Director_alt_Smear<AField<double,ACCEL> >::class_name
                         = "Director_alt_Smear<AField<double,ACCEL> >";

template class Director_alt_Smear<AField<double,ACCEL> >;

template<>
const std::string Director_alt_Smear<AField<float,ACCEL> >::class_name
                         = "Director_alt_Smear<AField<float,ACCEL> >";

template class Director_alt_Smear<AField<float,ACCEL> >;

//============================================================END=====
