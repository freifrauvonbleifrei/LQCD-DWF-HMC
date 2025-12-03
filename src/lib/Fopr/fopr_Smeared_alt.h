/*!
        @file    fopr_Smeared_alt.h

        @brief

        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: kanamori $

        @date    $LastChangedDate:: 2025-12-03 19:35:35 #$

        @version $LastChangedRevision: 2668 $
*/

#ifndef FOPR_SMEARED_ALT_INCLUDED
#define FOPR_SMEARED_ALT_INCLUDED

#include "lib/Fopr/afopr_Smeared_alt.h"

#ifdef USE_ALT_CODE


//! Domain-wall fermion operator.

/*!
    Alternative smeared fermion operatior with AFIELD=Field.
                                    [29 Mar 2023 H.Matsufuru]
 */

typedef AFopr_Smeared_alt<Field> Fopr_Smeared_alt;

#endif

#endif
