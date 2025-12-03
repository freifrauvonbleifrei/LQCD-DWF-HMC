/*!
        @file    complexTraits.h
        @brief
        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $
        @date    $LastChangedDate:: 2022-12-29 18:56:53 #$
        @version $LastChangedRevision: 2439 $
*/

#ifndef COMPLEX_TRAITS_INCLUDED
#define COMPLEX_TRAITS_INCLUDED

#include "bridge_complex.h"

template<typename REALTYPE>
class ComplexTraits;

template<>
class ComplexTraits<double> {
 public:
  typedef dcomplex complex_t;
};

template<>
class ComplexTraits<float> {
 public:
  typedef fcomplex complex_t;
};

#endif
