/*!
      @file    astaple_lex_float.cpp
      @brief
      @author  Hideo Matsufuru (matufuru)
               $LastChangedBy: matufuru $
      @date    $LastChangedDate:: 2025-08-20 16:24:03 #$
      @version $LastChangedRevision: 2653 $
*/

#include "lib_alt_Accel/Measurements/Gauge/astaple_lex.h"

#include <string>

#include "lib/ResourceManager/threadManager.h"

#include "lib_alt_Accel/inline/define_params.h"
#include "lib_alt_Accel/Field/aindex_lex.h"
#include "lib_alt_Accel/Field/afield.h"

// function template files
#include "lib_alt_Accel/Field/afield-inc.h"
#include "lib_alt_Accel/Field/afield_Gauge-inc.h"

// class template
#include "lib_alt_Accel/Measurements/Gauge/astaple_lex-tmpl.h"

template<>
const std::string AStaple_lex<AField<float,ACCEL> >::class_name
                              = "AStaple_lex<AField<float,ACCEL> >";

// explicit instanciation.
template class AStaple_lex<AField<float,ACCEL> >;

//============================================================END=====
