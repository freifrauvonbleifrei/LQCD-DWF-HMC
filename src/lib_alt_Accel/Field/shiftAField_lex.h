/*!
        @file    shiftAField_lex.h
        @brief
        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $
        @date    $LastChangedDate:: 2025-08-20 16:24:03 #$
        @version $LastChangedRevision: 2653 $
*/

#ifndef ACCEL_SHIFTAFIELD_LEX_INCLUDED
#define ACCEL_SHIFTAFIELD_LEX_INCLUDED

#include <vector>

#include "lib/Communicator/communicator.h"
#include "lib/Communicator/communicator_impl.h"
#include "lib/IO/bridgeIO.h"
using Bridge::vout;

//! Shift of a field in the lexical site index.

/*!
   This class is an Accel version of the ShiftField in
   the Bridge core libirary.
                                   [08 Aug 2025 H.Matsufuru]
 */
template<typename AFIELD>
class ShiftAField_lex {

 public:
  typedef typename AFIELD::real_t real_t;
  static const std::string class_name;

 private:
  int m_Nin;             //!< internal degree of freedom.
  int m_Nvol;
  int m_Ndim;
  int m_Nx, m_Ny, m_Nz, m_Nt;
  int m_Nxv, m_Nstv;
  std::vector<int> m_boundary;
  Bridge::VerboseLevel m_vl;

  int do_comm[4];  // switchs of communication (4=Ndim): (0: n, 1: y).
  int do_comm_any; // switchs of communication (if any): (0: n, 1: y).

  int m_Nsize[4];
  int m_bc[4];
  int m_bc2[4];

  std::vector<int> m_Nbdsize;
  using Channel = Channel_impl<std::allocator<real_t> >;
  std::vector<Channel> chsend_up, chrecv_up, chsend_dn, chrecv_dn;
  ChannelSet chset_send, chset_recv;

 public:
  ShiftAField_lex(int nin) { init(nin); }

  ShiftAField_lex(int nin, std::vector<int>& bc)
  { init(nin, bc); }

  ~ShiftAField_lex()
  { tidyup(); }

 private:
  // non-copyable
  ShiftAField_lex(const ShiftAField_lex&);
  ShiftAField_lex& operator=(const ShiftAField_lex&);

 public:
  void forward(AFIELD&, const AFIELD&, const int mu);
  void backward(AFIELD&, const AFIELD&, const int mu);

  void forward( AFIELD&, const int, const AFIELD&, const int,
                                                const int mu);
  void backward(AFIELD&, const int, const AFIELD&, const int,
                                                const int mu);

 private:

  //! setup channels for communication.
  void setup_channels();

  void init(int Nin);
  void init(int Nin, std::vector<int>& bc);

  void tidyup();

  void up_x(real_t *, real_t *);
  void up_y(real_t *, real_t *);
  void up_z(real_t *, real_t *);
  void up_t(real_t *, real_t *);
  void dn_x(real_t *, real_t *);
  void dn_y(real_t *, real_t *);
  void dn_z(real_t *, real_t *);
  void dn_t(real_t *, real_t *);

};

#endif
