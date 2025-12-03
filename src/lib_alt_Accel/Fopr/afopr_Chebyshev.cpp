/*!                                                                                    @file    afopr_Chebyshev.cpp
        @brief
        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $
        @date    $LastChangedDate:: 2025-08-20 16:24:03 #$
        @version $LastChangedRevision: 2653 $
*/

#include "lib/Fopr/afopr_Chebyshev.h"

#include <cassert>

#include "lib/Field/field.h"

#include "lib_alt_Accel/Field/afield.h"
#include "lib_alt_Accel/Field/afield-inc.h"

template<>
const std::string AFopr_Chebyshev<AField<double,ACCEL> >::class_name
                            = "AFopr_Chebyshev<AField<double,ACCEL> >";

template<>
const std::string AFopr_Chebyshev<AField<float, ACCEL> >::class_name
                             = "AFopr_Chebyshev<AField<float,ACCEL> >";

#ifdef USE_FACTORY
namespace {
  AFopr<AField<double,ACCEL> > *create_object_1(
                                      AFopr<AField<double,ACCEL> > *fopr,
                                      const Parameters& params)
  {
    return new AFopr_Chebyshev<AField<double,ACCEL> >(fopr, params);
  }
  bool init1 = AFopr<AField<double,ACCEL> >::Factory_fopr_params::Register(
                        "Chebyshev", create_object_1);

  AFopr<AField<float,ACCEL> > *create_object_2(
                                      AFopr<AField<float,ACCEL> > *fopr,
                                      const Parameters& params)
  {
    return new AFopr_Chebyshev<AField<float, ACCEL> >(fopr, params);
  }
  bool init2 = AFopr<AField<float, ACCEL> >::Factory_fopr_params::Register(
                        "Chebyshev", create_object_2);
}
#endif

#include "lib/Fopr/afopr_Chebyshev-tmpl.h"

// class instanciation.
template class AFopr_Chebyshev<AField<double,ACCEL> >;
template class AFopr_Chebyshev<AField<float, ACCEL> >;

//============================================================END=====
