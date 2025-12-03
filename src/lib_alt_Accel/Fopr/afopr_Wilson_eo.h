/*!
      @file    afopr_Wilson_eo.h
      @brief
      @author  Hideo Matsufuru (matufuru)
      @date    $LastChangedDate: 2013-01-22 13:51:53 #$
      @version $LastChangedRevision: 2653 $
*/

#ifndef ACCEL_AFOPR_WILSON_EO_H
#define ACCEL_AFOPR_WILSON_EO_H

#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>

#include "lib/Communicator/communicator.h"
#include "lib/Communicator/communicator_impl.h"
#include "lib/Fopr/afopr.h"

#include "lib_alt_Accel/Field/aindex_eo.h"

class Field;

//! Wilson fermion operator (even-odd) for Accel branch.

/*!
    This class implements the standard Wilson fermion operator
    with the even-odd site oerrdering for the Accel branch.
                                        [14 Feb 2025 H.Matsufuru]
 */

template<typename AFIELD>
class AFopr_Wilson_eo : public AFopr<AFIELD>
{
 public:
  typedef typename AFIELD::real_t real_t;
  static const std::string class_name;

 protected:
  enum Gamma_repr {DIRAC, CHIRAL} m_repr;  //!< gamma matrix representation

  int m_Nc, m_Nd, m_Ncol, m_Nvc, m_Ndf;
  int m_Nx, m_Ny, m_Nz, m_Nt, m_Nst, m_Ndim;
  int m_Nx2, m_Nst2;

  real_t  m_CKs;               //!< hopping parameter.
  std::vector<int> m_boundary; //!< pointer to boundary condition
  Bridge::VerboseLevel m_vl;   //!< verbose level

  AIndex_eo<real_t,AFIELD::IMPL> m_index_eo;

  int m_Ieo_origin;  //!< e-o of local lattice origin

  Field *m_conf;     //!< pointer to original gauge configuration
  AFIELD m_Ueo;      //!< gauge config. with even-odd separation

  std::string m_mode; //!< mult mode

  AFIELD m_v1, m_v2;  //!< working vectors

  int do_comm[4];  //!< switchs of communication (4=Ndim): (0: n, 1: y).
  int do_comm_any; //!< switchs of communication (if any): (0: n, 1: y).

  int m_Nsize[4];  //!< lattice size for kernel
  int m_bc[4];     //!< boundary condition for boundary part
  int m_bc2[4];    //!< boundary condition for bulk part

  std::vector<int> m_Nbdsize;
  using Channel = Channel_impl<std::allocator<real_t> >;
  std::vector<Channel> chsend_up, chrecv_up, chsend_dn, chrecv_dn;
  ChannelSet chset_send, chset_recv;

public:

  //! constructor.
  AFopr_Wilson_eo(const Parameters& params){ init(params); }

  //! destructor.
  ~AFopr_Wilson_eo(){ tidyup(); }

  //! setting parameters by a Parameter object.
  void set_parameters(const Parameters& params);

  //! setting parameters by values.
  void set_parameters(real_t CKs, std::vector<int> bc);

  //! getting parameters as a Parameters object.
  void get_parameters(Parameters& params) const;

  //! setting gauge configuration.
  void set_config(Field* u);
    
  //! returns the pointer to gauge configuration.
  inline Field* get_conf() { return m_conf; }

  //! setting mult mode.
  void set_mode(std::string mode);

  //! return the current mult mode.
  std::string get_mode() const { return m_mode; }

  void mult(AFIELD&, const AFIELD&);
  void mult_dag(AFIELD&, const AFIELD&);
  void mult_gm5(AFIELD&, const AFIELD&);

  void mult(AFIELD&, const AFIELD&,
            const std::string mode);

  void DdagD(AFIELD&, const AFIELD&);
  void Ddag(AFIELD&, const AFIELD&);
  void D(AFIELD&, const AFIELD&);

  //! Fermion matrix with ieo = 0: even <-- odd, 1: odd <-- even.
  void Meo(AFIELD&, const AFIELD&, const int ieo);

  //! Fermion matrix with ieo = 0: even <-- odd, 1: odd <-- even.
  void Meo(AFIELD&, const AFIELD&, const AFIELD&,
           const int ieo, const int iflag);

  //! Fermion matrix with ieo = 0: even <-- odd, 1: odd <-- even.
  void Meo_alt(AFIELD&, const AFIELD&, const int ieo);

  //! Fermion matrix with ieo = 0: even <-- odd, 1: odd <-- even.
  void Meo_alt(AFIELD&, const AFIELD&, const AFIELD&,
               const int ieo, const int iflag);

  //! returns inner size parameter.
  int field_nin() { return 2 * m_Nc * m_Nd; }

  //! returns local volume size parameter.
  int field_nvol() { return m_Nst2; }

  //! returns external size parameter.
  int field_nex() { return 1; }

  //! returns floating operation counts.
  double flop_count();

  //! returns floating operation counts.
  double flop_count(const std::string mode);

 private:
  //! initial setup.
  void init(const Parameters& params);

  //! final tidy-up.
  void tidyup();

  //! setup channels for communication.                                        
  void setup_channels();

  void set_config_omp(Field* u);
  void set_config_impl(Field* u);

  void mult_xp(real_t*, real_t*, const int);
  void mult_xm(real_t*, real_t*, const int);
  void mult_yp(real_t*, real_t*, const int);
  void mult_ym(real_t*, real_t*, const int);
  void mult_zp(real_t*, real_t*, const int);
  void mult_zm(real_t*, real_t*, const int);
  void mult_tp(real_t*, real_t*, const int);
  void mult_tm(real_t*, real_t*, const int);


#ifdef USE_FACTORY
 private:
  static AFopr<AFIELD> *create_object_with_params(const Parameters& params)
  { return new AFopr_Wilson_eo(params); }

 public:
  static bool register_factory()
  {
    bool init1 = AFopr<AFIELD>::Factory_params::Register("Wilson_eo",
                                                create_object_with_params);
    return init1;
  }
#endif

};

#endif  // ACCEL_AFOPR_WILSON_EO_H
