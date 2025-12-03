/*!
      @file    fprop_alt_Standard_lex.cpp
      @brief
      @author  Hideo Matsufuru (matufuru)
               $LastChangedBy: kanamori $
      @date    $LastChangedDate:: 2025-12-03 19:35:35 #$
      @version $LastChangedRevision: 2668 $
*/

#define DIRECTOR_ALT_SMEARED_IMPLEMENTED

#include "lib_alt/Measurements/Fermion/fprop_alt_Standard_lex.h"

#include "lib/ResourceManager/threadManager.h"
#include "lib/Fopr/afopr_Smeared.h"
#include "lib/Fopr/afopr_Smeared_alt.h"

#include "lib_alt_QXS/Field/afield.h"
#include "lib_alt_QXS/Field/afield-inc.h"
#include "lib_alt_QXS/Field/aindex_lex.h"

#include "lib_alt/Measurements/Fermion/fprop_alt_Standard_lex-tmpl.h"

//====================================================================
// explicit instanciation for AField<double,QXS>.
template<>
const std::string Fprop_alt_Standard_lex<AField<double, QXS> >::class_name
  = "Fprop_alt_Standard_lex<Afield<double,QXS> >";


template class Fprop_alt_Standard_lex<AField<double, QXS> >;

//====================================================================
// explicit instanciation for AField<float,QXS>.
template<>
const std::string Fprop_alt_Standard_lex<AField<float, QXS> >::class_name
  = "Fprop_alt_Standard_lex<Afield<float,QXS> >";


template class Fprop_alt_Standard_lex<AField<float, QXS> >;

//============================================================END=====
