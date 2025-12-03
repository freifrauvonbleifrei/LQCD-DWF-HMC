/*!
        @file    aeigensolver_IRArnoldi.cpp
        @brief
        @author  Hideo Matsufuru (matsufuru)
                 $LastChangedBy: matufuru $
        @date    $LastChangedDate:: 2025-08-20 16:24:03 #$
        @version $LastChangedRevision: 2653 $
*/

#include "lib/Eigen/aeigensolver_IRArnoldi.h"

#include "lib/ResourceManager/threadManager.h"
#include "lib/Fopr/afopr.h"

#include "lib_alt_Accel/inline/afield_th-inc.h"
#include "lib_alt_Accel/Field/afield.h"
#include "lib_alt_Accel/Field/afield-inc.h"

#include "lib/Eigen/aeigensolver_IRArnoldi-tmpl.h"

// explicit instanciation for AField<double,ACCEL>.
template<>
const std::string
  AEigensolver_IRArnoldi<AField<double,ACCEL>,
                         AFopr<AField<double,ACCEL> > >::class_name
  = "AEigensolver_IRArnoldi<AField<double,ACCEL> >";

#ifdef USE_FACTORY_AUTOREGISTER
namespace {
  bool init2 = AEigensolver<AField<double,ACCEL>,
                            AFopr<AField<double, ACCEL> > >::
       Factory_params::Register("IRArnoldi", create_object);
}
#endif

template class AEigensolver_IRArnoldi<AField<double,ACCEL>,
                                      AFopr<AField<double,ACCEL> > >;

//============================================================END=====
