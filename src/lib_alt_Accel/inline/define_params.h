/*!
        @file    define_params_SU3.h
        @brief
        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $
        @date    $LastChangedDate:: 2025-08-20 16:24:03 #$
        @version $LastChangedRevision: 2653 $
*/

#ifndef ACCEL_DEFINE_PARAMS_INCLUDED
#define ACCEL_DEFINE_PARAMS_INCLUDED

#if defined USE_GROUP_SU3
//#include "lib_alt_Accel/inline/define_params_SU3.h"

//temporary setting
#define  NC      3
#define  NVC     6
#define  NDF    18
#define  NDF2    9
#define  ND      4
#define  ND2     2
#define  NCD    12
#define  NVCD   24
#define  NVCD2  12
#define  NDIM    4

#endif

// Warp length
//#define  NWP     1
#define  NWP     32

// OpenACC runtime paramters
#define NUM_WORKERS 1
#define VECTOR_LENGTH 32

// selection of implementation

#define USE_DOMAINWALL_5DIN_EE_MATINV_KERNEL
// set if implement Dee_inv by matrix inversion
 
#define USE_DOMAINWALL_5DIN_4D_KERNEL
// set if implement Domainwall kernel by 4D site loop

// function for padding
#define CEIL_NWP(nst) ((((nst) + NWP-1)/NWP) * NWP)
// this definition must be the same as ceil_nwp() function below

namespace {

  inline int ceil_nwp(const int nst)
  {
    int size = nst/NWP;
    if( (nst % NWP) > 0) ++size;
    return size * NWP;
  }

}

#endif
//============================================================END=====
