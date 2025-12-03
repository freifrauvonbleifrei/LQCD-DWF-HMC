/*!
      @file    bridgeACC.h
      @brief
      @author  Hideo Matsufuru (matufuru)
               $LastChangedBy: matufuru $
      @date    $LastChangedDate:: 2025-08-20 16:24:03 #$
      @version $LastChangedRevision: 2653 $
*/

#ifndef BRIDGEACC_INCLUDED
#define BRIDGEACC_INCLUDED

#include <cstddef>

namespace BridgeACC {

// real_t = double
void enter_data_create(double* data, size_t size);
void exit_data_delete(double* data, size_t size);

void update_device(double* x, size_t offset, size_t size);
void update_host(double* x, size_t offset, size_t size);

// real_t = float
void enter_data_create(float* data, size_t size);
void exit_data_delete(float* data, size_t size);

void update_device(float* x, size_t offset, size_t size);
void update_host(float* x, size_t offset, size_t size);

}
#endif
