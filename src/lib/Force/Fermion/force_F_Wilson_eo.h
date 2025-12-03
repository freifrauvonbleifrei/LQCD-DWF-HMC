/*!
        @file    force_F_Wilson_eo.h

        @brief

        @author  UEDA, Satoru (sueda)
                 $LastChangedBy: kanamori $

        @date    $LastChangedDate:: 2025-12-03 19:35:35 #$

        @version $LastChangedRevision: 2668 $
*/

#ifndef FORCE_F_WILSON_EO_INCLUDED
#define FORCE_F_WILSON_EO_INCLUDED

#include "force_F.h"
#include "tensorProd.h"

#include "Fopr/fopr_Wilson_eo.h"
#include "Fopr/fopr_Wilson.h"

#include "IO/bridgeIO.h"
using Bridge::vout;

//! Force for the Wilson fermion operator with even-odd precondition

/*!
    This class calculates the force of the standard Wilson fermion.
    The gamma matrix representation is given as control string
    "Dirac"(default) or "Chiral" at the construction, which is
    used to construct the Fopr_Wilson instance.
                                     [19 June 2012 S.UEDA]
    (Coding history will be recovered from trac.)
    YAML is implemented.             [14 Nov 2012 Y.Namekawa]

    Reivised and now it works with Action_F_Standared_eo.
    See also the implementation note by I.K.
                                     [21? Apr 2023 I.Kanamori ]

*/


class Force_F_Wilson_eo : public Force
{
 public:
  static const std::string class_name;

 private:
  Bridge::VerboseLevel m_vl;

  double m_kappa;
  std::vector<int> m_boundary;
  Fopr_Wilson_eo *m_fopr_w;

  Field_F m_psf;
  std::string m_repr;

  Index_eo m_index;

 public:
  DEPRECATED
  Force_F_Wilson_eo()
    : m_vl(CommonParameters::Vlevel())
  {
    std::string repr("Dirac");
    init_fopr(repr);
    m_boundary.resize(CommonParameters::Ndim());
  }

  DEPRECATED
  Force_F_Wilson_eo(const std::string repr)
    : m_vl(CommonParameters::Vlevel())
  {
    init_fopr(repr);
    m_boundary.resize(CommonParameters::Ndim());
  }

  Force_F_Wilson_eo(const Parameters& params)
    : m_vl(CommonParameters::Vlevel())
  {
    std::string repr = params.get_string("gamma_matrix_type");
    init_fopr(repr);
    m_boundary.resize(CommonParameters::Ndim());
    set_parameters(params);
  }

  ~Force_F_Wilson_eo()
  {
    delete m_fopr_w;
  }

  void set_parameters(const Parameters& params);
  void set_parameters(const double kappa, const std::vector<int> bc);

  void get_parameters(Parameters& params) const;

  void set_config(Field *U)
  {
    m_U = (Field_G *)U;
    m_fopr_w->set_config(U);
  }

  //! return the force field, differentiated by the gauge field U form the fermion field on the even site.
  //! force_udiv =
  //!   <eta_e| Ddag D_deriv | eta_e>  + <eta_e| Ddag_deriv D | eta_e>
  //!   <eta_e| H  H_deriv | eta_e>  + <eta_e| H_deriv H | eta_e>
  //! with Ddag = gm5 D gm5, H = gm5 D
  void force_udiv(Field& force, const Field& eta_e);

  // inuput must be |eta_e> and H|eta_e> 
  void force_udiv1(Field& force_eo, const Field& zeta_eo, const Field& eta_eo);

private:
  void force_udiv1_impl(Field_G& force_eo, const Field_F& zeta_eo, const Field_F& eta_eo, const int ieo); // not used

  void init_fopr(const std::string &repr);
  void init_fopr(const Parameters &params);

};
#endif
