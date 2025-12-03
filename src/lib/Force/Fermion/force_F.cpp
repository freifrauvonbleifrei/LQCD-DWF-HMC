/*!
        @file    force_F.cpp

        @brief

        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate:: 2025-09-02 15:10:15 #$

        @version $LastChangedRevision: 2654 $
*/

#include "force_F.h"

template<>
const std::string Force::class_name = "Force";
//const std::string AForce_F<Field>::class_name = "Force_F";

//====================================================================
template<>
void Force::init()
{
  // do nothing.
}

//====================================================================
template<>
void Force::tidyup()
{
  // do nothing.
}

//====================================================================
template<>
void Force::mult_generator(Field& force)
{
  vout.crucial("%s: mult_generator(Field&) must not be called\n",
	       class_name.c_str());
  exit(EXIT_FAILURE);

}

//====================================================================
template<>
void Force::mult_generator(Field_G& force)
{
  const int Nvol = force.nvol();
  const int Ndim = force.nex();

  for (int mu = 0; mu < Ndim; ++mu) {
    for (int isite = 0; isite < Nvol; ++isite) {
      Mat_SU_N u = m_U->mat(isite, mu);

      u *= force.mat(isite, mu);
      u.at();
      u *= -2.0;

      force.set_mat(isite, mu, u);
    }
  }
}

//====================================================================
// default templates for core and core1
template<>
void Force::force_core(Field& force_, const Field& eta)
{
  const int Nvol = CommonParameters::Nvol();
  const int Ndim = CommonParameters::Ndim();

  Field_G force(Nvol, Ndim);

  force_udiv(force, eta);

  mult_generator(force);

  force_ = force;
}


//====================================================================
template<>
void Force::force_core1(Field& force_, const Field& zeta, const Field& eta)
{
  const int Nvol = CommonParameters::Nvol();
  const int Ndim = CommonParameters::Ndim();

  Field_G force(Nvol, Ndim);

  force_udiv1(force, zeta, eta);

  mult_generator(force);

  force_ = force;
}

//====================================================================


// explicit instanciation.
template class AForce_F<Field>;

//============================================================END=====
