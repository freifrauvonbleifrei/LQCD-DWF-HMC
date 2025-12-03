/*!
        @file    bridge_init_factory_alt_Accel.h
        @brief
        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate:: 2025-08-20 16:24:03 #$

        @version $LastChangedRevision: 2653 $
*/

#ifndef BRIDGE_INIT_FACTORY_ALT_ACCEL_INCLUDED
#define BRIDGE_INIT_FACTORY_ALT_ACCEL_INCLUDED


#ifdef USE_FACTORY

#ifdef USE_FACTORY_AUTOREGISTER
#else
bool bridge_init_factory_alt_Accel();
#endif
  
void bridge_report_factory_alt_Accel();

#endif /* USE_FACTORY */

#endif
