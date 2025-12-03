/*!
      @file    director_alt_Smear.cpp
      @brief
      @author  Hideo Matsufuru (matufuru)
               $LastChangedBy: kanamori $
      @date    $LastChangedDate:: 2025-12-03 19:35:35 #$
      @version $LastChangedRevision: 2668 $
*/

#include "lib_alt/Smear/director_alt_Smear.h"

#include "lib/ResourceManager/threadManager.h"
#include "lib/Fopr/afopr_Smeared.h"

#include "lib_alt_QXS/Field/afield.h"
#include "lib_alt_QXS/Field/afield-inc.h"
#include "lib_alt_QXS/Field/aindex_lex.h"

#include "lib_alt/Smear/director_alt_Smear-tmpl.h"

//====================================================================
// explicit instanciation for AField<double,QXS>.
template<>
const std::string Director_alt_Smear<AField<double,QXS> >::class_name
                         = "Director_alt_Smear<AField<double,QXS> >";

template class Director_alt_Smear<AField<double,QXS> >;

template<>
const std::string Director_alt_Smear<AField<float,QXS> >::class_name
                         = "Director_alt_Smear<AField<float,QXS> >";

template class Director_alt_Smear<AField<float,QXS> >;

//============================================================END=====
