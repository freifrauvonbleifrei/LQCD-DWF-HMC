/*!
        @file    define_available.h
        @brief
        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: kanamori $
        @date    $LastChangedDate:: 2025-12-03 19:35:35 #$
        @version $LastChangedRevision: 2668 $
*/

#ifndef ACCEL_DEFINE_AVAILABLE_INCLUDED
#define ACCEL_DEFINE_AVAILABLE_INCLUDED

#ifdef USE_ACCEL_OPENACC

#define ACCEL_FOPR_WILSON_AVAILABLE
#define ACCEL_FOPR_WILSON_EO_AVAILABLE
#define ACCEL_FOPR_CLOVER_AVAILABLE
#define ACCEL_FOPR_CLOVER_EO_AVAILABLE
#define ACCEL_FOPR_STAGGERED_AVAILABLE
#define ACCEL_FOPR_STAGGERED_EO_AVAILABLE
#define ACCEL_FOPR_DOMAINWALL_AVAILABLE
#define ACCEL_FOPR_DOMAINWALL_EO_AVAILABLE
#define ACCEL_FOPR_DOMAINWALL_5DIN_AVAILABLE
#define ACCEL_FOPR_DOMAINWALL_5DIN_EO_AVAILABLE

#define ACCEL_FORCE_DOMAINWALL_5DIN_AVAILABLE
#define ACCEL_FORCE_DOMAINWALL_5DIN_EO_AVAILABLE

#endif  // USE_ACCEL_OPENACC


#endif  // ACCEL_DEFINE_AVAILABLE_INCLUDED
//============================================================END=====
