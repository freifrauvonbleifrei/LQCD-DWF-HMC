/*!
        @file    configure.h

        @brief

        @author  Tatsumi Aoyama  (aoym)
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate:: 2025-09-02 15:10:15 #$

        @version $LastChangedRevision: 2654 $
*/

#ifndef CONFIGURE_INCLUDED
#define CONFIGURE_INCLUDED

#include "version.h"

// some #define options are set as compiler options from Makefile

//#define USE_MPI
#define ENABLE_MULTI_INSTANCE

// smart pointer support
#include <memory>
using std::unique_ptr;

#define DEPRECATED    [[deprecated]]
// #define DEPRECATED __attribute__((deprecated))
// #define DEPRECATED

#endif /* CONFIGURE_INCLUDED */
