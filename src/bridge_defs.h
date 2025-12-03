/*!
        @file    bridge_defs.h

        @brief

        @author  Tatsumi Aoyama  (aoym)
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate:: 2025-09-02 15:10:15 #$

        @version $LastChangedRevision: 2654 $
*/

#ifndef DEFS_INCLUDED
#define DEFS_INCLUDED

#include <string>

// for debug
//#define LOG printf(">>> %s\n", __PRETTY_FUNCTION__)
#define LOG

// direction label
enum Direction
{
  XDIR = 0,
  YDIR = 1,
  ZDIR = 2,
  TDIR = 3,
  WDIR = 4
};

enum ForwardBackward
{
  Forward  =  1,  // +mu
  Backward = -1   // -mu
};

namespace Element_type
{
  enum type
  {
    REAL = 1, COMPLEX = 2
  };
}

// restrict keyword
#ifdef __GNUC__
#define RESTRICT  __restrict__
#else
#ifdef PC_NVIDIA
#define RESTRICT  restrict
#endif
#endif


#endif /* DEFS_INCLUDED */
