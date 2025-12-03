/*!
      @file    afield_double.h
      @brief
      @author  Hideo Matsufuru (matufuru)
               $LastChangedBy: matufuru $
      @date    $LastChangedDate:: 2025-08-20 16:24:03 #$
      @version $LastChangedRevision: 2653 $
*/

#include "lib_alt_Accel/Field/afield.h"

#include <cassert>

#include "lib/ResourceManager/threadManager.h"

#include "lib_alt_Accel/Field/aindex_lex.h"

typedef double real_t;

#include "lib_alt_Accel/inline/afield_th-inc.h"

// library
#include "lib_alt_Accel/BridgeACC/bridgeACC.h"
#include "lib_alt_Accel/BridgeACC/bridgeACC_AField.h"
#include "lib_alt_Accel/BridgeACC/bridgeACC_AField_Gauge.h"

// class template file
#include "lib_alt_Accel/Field/afield-tmpl.h"

template<>
const std::string AField<double, ACCEL>::class_name = "AField<double, ACCEL>";

// explicit instanciation.
template class AField<double, ACCEL>;

//============================================================END=====
