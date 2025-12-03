/*!
        @file    afopr_common_th-inc.h
        @brief
        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $
        @date    $LastChangedDate:: 2025-08-20 16:24:03 #$
        @version $LastChangedRevision: 2653 $
*/

#ifndef ACCEL_AFOPR_COMMON_TH_INC_INCLUDED
#define ACCEL_AFOPR_COMMON_TH_INC_INCLUDED


namespace {

//====================================================================
// case (a): tasks are equally assgined to all threads

inline void set_threadtask_afopr(int& ith, int& nth, int& is, int& ns,
                           const int size)
{
  nth = ThreadManager::get_num_threads();
  ith = ThreadManager::get_thread_id();

  is = size * ith / nth;
  ns = size * (ith + 1) / nth;

}

} // nameless namespace end
#endif
//============================================================END=====
