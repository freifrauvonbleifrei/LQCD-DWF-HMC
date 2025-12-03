#ifndef FAPP_MACROS_H_INCLUDED
#define FAPP_MACROS_H_INCLUDED

#ifdef USE_FAPP
#include "fj_tool/fapp.h"

#define START_FAPP(name, number, level)         \
  fapp_start( name, number, level )
#define STOP_FAPP(name, number, level)  \
  fapp_stop( name, number, level )

#else
#define START_FAPP(name, number, level)
#define STOP_FAPP(name, number, level)
#endif

#endif
