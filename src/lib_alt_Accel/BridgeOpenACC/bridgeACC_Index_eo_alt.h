/*!
      @file    index_eo_alt_openacc.h
      @brief
      @author  Hideo Matsufuru (matufuru)
               $LastChangedBy: matufuru $
      @date    $LastChangedDate:: 2025-08-20 16:24:03 #$
      @version $LastChangedRevision: 2653 $
*/

#ifndef BRIDGEACC_INDEX_EO_ALT_INCLUDED
#define BRIDGEACC_INDEX_EO_ALT_INCLUDED

#include "bridge_defs.h"

namespace BridgeACC {

// real_t = double
void split(double* RESTRICT ve, double* RESTRICT vo,
           double* RESTRICT w, int ieo_origin,
           int nin, int* Nsize);

void merge(double* RESTRICT v,
           double* RESTRICT we, double* RESTRICT wo, int ieo_origin,
           int nin, int* Nsize);

// real_t = float
void split(float* RESTRICT ve, float* RESTRICT vo,
           float* RESTRICT w, int ieo_origin,
           int nin, int* Nsize);

void merge(float* RESTRICT v,
           float* RESTRICT we, float* RESTRICT wo, int ieo_origin,
           int nin, int* Nsize);

}

//============================================================END=====
#endif
