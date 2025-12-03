/*!
        @file    fprop_alt.cpp
        @brief
        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $
        @date    $LastChangedDate:: 2025-08-20 16:24:03 #$
        @version $LastChangedRevision: 2653 $
*/

#include "lib_alt/Measurements/Fermion/fprop_alt.h"

#include "lib_alt_Accel/Field/afield.h"


// explicit instanciation

template class Fprop_alt<AField<double,ACCEL> >;

template class Fprop_alt<AField<float,ACCEL> >;

//============================================================END=====
