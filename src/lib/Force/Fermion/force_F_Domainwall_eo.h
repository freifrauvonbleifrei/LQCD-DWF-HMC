/*!
        @file    force_F_Domainwall_eo.h

        @brief

        @author  Issaku Kanamori (kanamori)
                 $LastChangedBy: kanamori $

        @date    $LastChangedDate:: 2025-12-03 19:35:35 #$

        @version $LastChangedRevision: 2668 $
*/

#ifndef FORCE_F_DOMAINWALL_EO_INCLUDED
#define FORCE_F_DOMAINWALL_EO_INCLUDED

#include "force_F.h"
#include "tensorProd.h"

#include "Field/index_eo.h"
#include "Fopr/fopr_Domainwall_eo.h"
#include "Fopr/fopr_Wilson_eo.h"

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


class Force_F_Domainwall_eo : public Force
{
 public:
  static const std::string class_name;

 private:
  Bridge::VerboseLevel m_vl;

  // parameters for the fermion operator
  double m_mq;                 //!< quark mass
  double m_M0;                 //!< domain-wall height
  int m_Ns;                    //!< size of fifth-dimension
  std::vector<int> m_boundary; //!< boundary conditions
  std::vector<double> m_b;
  std::vector<double> m_c;

  std::string m_mode;
  Fopr_Domainwall_eo *m_fopr;
  //Fopr_Wilson_eo *m_fopr;

  Field_F m_psf;
  std::string m_repr;

  Index_eo m_index;

 public:
  DEPRECATED
  Force_F_Domainwall_eo()
    : m_vl(CommonParameters::Vlevel())
  {
    m_repr="";
    m_fopr=nullptr;
    m_boundary.resize(CommonParameters::Ndim());
  }

  DEPRECATED
  Force_F_Domainwall_eo(const std::string repr)
    : m_vl(CommonParameters::Vlevel())
  {
    m_repr=repr;
    m_fopr=nullptr;
    m_boundary.resize(CommonParameters::Ndim());
  }

  Force_F_Domainwall_eo(const Parameters& params)
    : m_vl(CommonParameters::Vlevel())
  {
    m_repr="";
    m_fopr=nullptr;
    m_boundary.resize(CommonParameters::Ndim());
    set_parameters(params);
  }

  ~Force_F_Domainwall_eo()
  {
    delete m_fopr;
  }

  void set_parameters(const Parameters& params);

  /*
  void set_parameters(const real_t mq, const real_t M0,
                      const int Ns,
                      const std::vector<int> bc,
                      const real_t b, const real_t c);
  */

  void get_parameters(Parameters& params) const;

  void set_config(Field *U)
  {
    m_U = (Field_G *)U;
    m_fopr->set_config(U);
  }

  void set_mode(const std::string& mode);

  //! return the force field, differentiated by the gauge field U form the fermion field on the even site.
  //! force_udiv =
  //!   <eta_e| Ddag D_deriv | eta_e>  + <eta_e| Ddag_deriv D | eta_e>
  //!   <eta_e| Hdag  Hdag_deriv | eta_e>  + <eta_e| Hdag_deriv Hdag | eta_e>

  //! with H = Gm5 D,  Hdag !=H  in general
  void force_udiv(Field& force, const Field& eta_e);

  //
  void force_udiv1(Field& force_eo, const Field& zeta_eo, const Field& eta_eo);

 private:
  void force_udiv1_H(Field& force_eo, const Field& zeta_eo, const Field& eta_eo);
  void force_udiv1_Hdag(Field& force_eo, const Field& zeta_eo, const Field& eta_eo);
  void force_udiv1_impl(Field_G& force_eo, const Field_F& zeta_eo, const Field_F& eta_eo, const int ieo);


};
#endif
