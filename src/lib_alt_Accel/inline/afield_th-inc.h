/*!
        @file    afield_th-inc.h
        @brief
        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $
        @date    $LastChangedDate:: 2025-08-20 16:24:03 #$
        @version $LastChangedRevision: 2653 $
*/

#ifndef ACCEL_AFIELD_TH_INC_H
#define ACCEL_AFIELD_TH_INC_H

namespace {

//====================================================================
inline void set_thread(int& ith, int& nth)
{
  nth = ThreadManager::get_num_threads();
  ith = ThreadManager::get_thread_id();
}

//====================================================================
inline void set_threadtask(int& ith, int& nth, int& is, int& ns,
                           const int size)
{
  nth = ThreadManager::get_num_threads();
  ith = ThreadManager::get_thread_id();

  size_t is2 = size_t(size) * size_t(ith) / nth;
  size_t ns2 = size_t(size) * size_t(ith + 1) / nth;
  is = int(is2);
  ns = int(ns2);

}

//====================================================================

} // nameless namespace end
#endif
