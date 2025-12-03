/*!
        @file    afopr_Domainwall_5din_eo-tmpl.h
        @brief
        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: kanamori $
        @date    $LastChangedDate:: 2025-12-03 19:35:35 #$
        @version $LastChangedRevision: 2668 $
*/

#include "lib_alt_Accel/Fopr/afopr_Domainwall_5din_eo.h"

#include "lib_alt_Accel/Fopr/afopr_common_th-inc.h"

template<typename AFIELD>
const std::string AFopr_Domainwall_5din_eo<AFIELD>::class_name
                                        = "AFopr_Domainwall_5din_eo";

#define AFOPR_TIMER

#ifdef AFOPR_TIMER
#include "lib/Tools/timer.h"
#define START_TIMER(var_timer) var_timer->start()
#define STOP_TIMER(var_timer)  var_timer->stop()
#else
#define START_TIMER(var_timer)
#define STOP_TIMER(var_timer)
#endif

//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_5din_eo<AFIELD>::init(const Parameters& params)
{
  ThreadManager::assert_single_thread(class_name);

  std::string vlevel;
  if (!params.fetch_string("verbose_level", vlevel)) {
    m_vl = vout.set_verbose_level(vlevel);
  } else {
    m_vl = CommonParameters::Vlevel();
  }

  vout.general(m_vl, "%s: construction\n", class_name.c_str());
  vout.increase_indent();

  m_repr = "Dirac";  // now only the Dirac repr is available.

  std::string repr;
  if (!params.fetch_string("gamma_matrix_type", repr)) {
    if (repr != "Dirac") {
      vout.crucial("Error at %s: unsupported gamma-matrix type: %s\n",
                   class_name.c_str(), repr.c_str());
      exit(EXIT_FAILURE);
    }
  }

  int Nc = CommonParameters::Nc();
  if (Nc != 3) {
    vout.crucial("%s: only applicable to Nc = 3\n",
                 class_name.c_str());
    exit(EXIT_FAILURE);
  }

  // reset the timers
#ifdef AFOPR_TIMER
  timer_mult_Deo.reset(       new Timer("afopr_Domainwall_5din_eo: Deo & dagger   "));
  timer_mult_Dee_inv.reset(   new Timer("afopr_Domainwall_5din_eo: Dee_inv & dag  "));
  timer_pack.reset(           new Timer("afopr_Domainwall_5din_eo: pack           "));
  timer_bulk.reset(           new Timer("afopr_Domainwall_5din_eo: bulk           "));
  timer_boundary.reset(       new Timer("afopr_Domainwall_5din_eo: boundary       "));
  timer_comm.reset(           new Timer("afopr_Domainwall_5din_eo: comm           "));
  timer_comm_recv_wait.reset( new Timer("afopr_Domainwall_5din_eo: comm_recv_wait "));
  timer_comm_send_wait.reset( new Timer("afopr_Domainwall_5din_eo: comm_send_wait "));
  timer_comm_recv_start.reset(new Timer("afopr_Domainwall_5din_eo: comm_recv_start"));
  timer_comm_send_start.reset(new Timer("afopr_Domainwall_5din_eo: comm_send_start"));
#pragma omp barrier
  vout.detailed(m_vl, "%s: detailed timer was initialized\n");
#endif

  m_Nd   = CommonParameters::Nd();
  m_Nd2  = m_Nd/2;
  m_Nvcd = 2 * Nc * m_Nd;
  m_Ndf  = 2 * Nc * Nc;

  m_Nx   = CommonParameters::Nx();
  m_Ny   = CommonParameters::Ny();
  m_Nz   = CommonParameters::Nz();
  m_Nt   = CommonParameters::Nt();
  m_Nvol = CommonParameters::Nvol();
  m_Ndim = CommonParameters::Ndim();

  m_Nx2  = m_Nx / 2;
  m_Nst2 = m_Nvol / 2;

  int ipe3 = Communicator::ipe(3);
  int ipe2 = Communicator::ipe(2);
  int ipe1 = Communicator::ipe(1);
  m_Ieo_origin = (ipe1 * m_Ny + ipe2 * m_Nz + ipe3 * m_Nt) % 2;

  // condition check
  if (m_Nx % 2 != 0) {
    vout.crucial(m_vl, "%s: Nx must be even.\n",
                 class_name.c_str());
    exit(EXIT_FAILURE);
  }

  m_Nsize[0] = m_Nx2;
  m_Nsize[1] = m_Ny;
  m_Nsize[2] = m_Nz;
  m_Nsize[3] = m_Nt;

  // switches for coomunication
  int req_comm = 0;  //  0: communication only if necessary
                     //  1: communication is enforced any time
  if (!params.fetch_int("require_communication", req_comm)) {
    vout.general(m_vl, "req_comm = %d (input)\n", req_comm);
  } else {
    vout.general(m_vl, "req_comm = %d (default)\n", req_comm);
  }

  do_comm_any = 0;
  for (int mu = 0; mu < m_Ndim; ++mu) {
    do_comm[mu] = 1;
    if ((req_comm == 0) && (Communicator::npe(mu) == 1)) do_comm[mu] = 0;
    do_comm_any += do_comm[mu];
    vout.general("do_comm[%d] = %d\n", mu, do_comm[mu]);
  }

  m_Ns = 0; // temporary set
  set_parameters(params);

  m_Nbdsize.resize(m_Ndim);
  int Nbdin = (m_Nvcd / 2) * m_Ns;
  m_Nbdsize[0] = Nbdin * ceil_nwp((m_Ny * m_Nz * m_Nt + 1)/2);
  m_Nbdsize[1] = Nbdin * ceil_nwp(m_Nx2 * m_Nz * m_Nt);
  m_Nbdsize[2] = Nbdin * ceil_nwp(m_Nx2 * m_Ny * m_Nt);
  m_Nbdsize[3] = Nbdin * ceil_nwp(m_Nx2 * m_Ny * m_Nz);

  setup_channels();

  // gauge configuration.
  int Nst_pad2 = 2 * ceil_nwp(m_Nst2);
  m_Ueo.reset(m_Ndf, Nst_pad2, m_Ndim);

  vout.decrease_indent();
  vout.detailed(m_vl, "%s: initalization finished.\n",
                class_name.c_str());
}

//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_5din_eo<AFIELD>::tidyup()
{
  ThreadManager::assert_single_thread(class_name);

  // openacc device memory clean up
  for(int mu = 0; mu < m_Ndim; ++mu){

#ifdef USE_MPI
    real_t* buf_dn1 = (real_t*)chsend_dn[mu].ptr();
    real_t* buf_dn2 = (real_t*)chrecv_dn[mu].ptr();
    real_t* buf_up1 = (real_t*)chrecv_up[mu].ptr();
    real_t* buf_up2 = (real_t*)chsend_up[mu].ptr();
    BridgeACC::afield_tidyup(buf_dn1, m_Nbdsize[mu]);
    BridgeACC::afield_tidyup(buf_dn2, m_Nbdsize[mu]);
    BridgeACC::afield_tidyup(buf_up1, m_Nbdsize[mu]);
    BridgeACC::afield_tidyup(buf_up2, m_Nbdsize[mu]);
#else
    real_t* buf_up = (real_t*)chsend_up[mu].ptr();
    real_t* buf_dn = (real_t*)chsend_dn[mu].ptr();
    BridgeACC::afield_tidyup(buf_up, m_Nbdsize[mu]);
    BridgeACC::afield_tidyup(buf_dn, m_Nbdsize[mu]);
#endif

  }

#ifdef AFOPR_TIMER
  timer_mult_Deo->report();
  timer_mult_Dee_inv->report();
  timer_pack->report();
  timer_bulk->report();
  timer_boundary->report();
  timer_comm->report();
  timer_comm_recv_wait->report();
  timer_comm_send_wait->report();
  timer_comm_recv_start->report();
  timer_comm_send_start->report();
#endif

}


//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_5din_eo<AFIELD>::setup_channels()
{
  ThreadManager::assert_single_thread(class_name);

  chsend_up.resize(m_Ndim);
  chrecv_up.resize(m_Ndim);
  chsend_dn.resize(m_Ndim);
  chrecv_dn.resize(m_Ndim);

  for (int mu = 0; mu < m_Ndim; ++mu) {
    size_t Nvsize = m_Nbdsize[mu] * sizeof(real_t);

    chsend_dn[mu].send_init(Nvsize, mu, -1);
    chsend_up[mu].send_init(Nvsize, mu, 1);
#ifdef USE_MPI
    chrecv_up[mu].recv_init(Nvsize, mu, 1);
    chrecv_dn[mu].recv_init(Nvsize, mu, -1);
#else
    void *buf_up = (void *)chsend_dn[mu].ptr();
    chrecv_up[mu].recv_init(Nvsize, mu, 1, buf_up);
    void *buf_dn = (void *)chsend_up[mu].ptr();
    chrecv_dn[mu].recv_init(Nvsize, mu, -1, buf_dn);
#endif

    if (do_comm[mu] == 1) {
      chset_send.append(chsend_up[mu]);
      chset_send.append(chsend_dn[mu]);
      chset_recv.append(chrecv_up[mu]);
      chset_recv.append(chrecv_dn[mu]);
    }

    // openacc device memory allocation
#ifdef USE_MPI
    real_t* buf_dn1 = (real_t*)chsend_dn[mu].ptr();
    real_t* buf_dn2 = (real_t*)chrecv_dn[mu].ptr();
    real_t* buf_up1 = (real_t*)chsend_up[mu].ptr();
    real_t* buf_up2 = (real_t*)chrecv_up[mu].ptr();
    BridgeACC::afield_init(buf_dn1, m_Nbdsize[mu]);
    BridgeACC::afield_init(buf_dn2, m_Nbdsize[mu]);
    BridgeACC::afield_init(buf_up1, m_Nbdsize[mu]);
    BridgeACC::afield_init(buf_up2, m_Nbdsize[mu]);
#else
    BridgeACC::afield_init((real_t*)buf_up, m_Nbdsize[mu]);
    BridgeACC::afield_init((real_t*)buf_dn, m_Nbdsize[mu]);
#endif

  }

}

//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_5din_eo<AFIELD>::set_parameters(
  const Parameters& params)
{
  std::string vlevel;
  if (!params.fetch_string("verbose_level", vlevel)) {
    m_vl = vout.set_verbose_level(vlevel);
  }
 
  string           gmset_type;
  double           mq, M0;
  int              Ns;
  std::vector<int> bc;
  double           b, c, alpha;

  int err_optional = 0;
  err_optional += params.fetch_string("gamma_matrix_type", m_repr);
  if(err_optional){
    vout.crucial(m_vl, "  gamma_matrix_type is not given\n");
  }

  err_optional = 0;
  err_optional += params.fetch_string("code_implementation", m_impl);
  if(err_optional){
    vout.crucial(m_vl, "  code_implementation is not given\n");
    m_impl = "5d";
  }

  int err = 0;
  err += params.fetch_double("quark_mass", mq);
  err += params.fetch_double("domain_wall_height", M0);
  err += params.fetch_int("extent_of_5th_dimension", Ns);
  err += params.fetch_int_vector("boundary_condition", bc);

  if (err) {
    vout.crucial(m_vl, "Error at %s: input parameter not found.\n",
                 class_name.c_str());
    exit(EXIT_FAILURE);
  }

  int err2 = 0;
  err2 += params.fetch_double("coefficient_b", b);
  err2 += params.fetch_double("coefficient_c", c);

  if (err2) {
    vout.general(m_vl, "  coefficients b, c are not provided:"
                       " set to Shamir's form.\n");
    b = 1.0;
    c = 0.0;
  }

  int err3 = 0;
  err3 += params.fetch_double("parameter_alpha", alpha);
  if (err3) {
    vout.general(m_vl, "  parameter alpha is not provided: set to 1.0.\n");
    alpha = 1.0;
  }

  set_parameters(real_t(mq), real_t(M0), Ns, bc,
                 real_t(b), real_t(c), real_t(alpha));

}

//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_5din_eo<AFIELD>::get_parameters(
                                            Parameters& params) const
{
  params.set_string("kernel_type", m_kernel_type);
  params.set_string("gamma_matrix_type", m_repr);
  params.set_string("code_implementation", m_impl);
  params.set_double("quark_mass", double(m_mq));
  params.set_double("domain_wall_height", double(m_M0));
  params.set_int("extent_of_5th_dimension", m_Ns);
  params.set_int_vector("boundary_condition", m_boundary);
  params.set_double("coefficient_b", double(m_b[0]));
  params.set_double("coefficient_c", double(m_c[0]));
  params.set_double("parameter_alpha", double(m_alpha));
  params.set_string("gamma_matrix_type", m_repr);

  params.set_string("verbose_level", vout.get_verbose_level(m_vl));
}


//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_5din_eo<AFIELD>::set_parameters(
                                           const real_t mq,
                                           const real_t M0,
                                           const int Ns,
                                           const std::vector<int> bc,
                                           const real_t b,
                                           const real_t c,
                                           const real_t alpha)
{
#pragma omp barrier

  int ith = ThreadManager::get_thread_id();

  if (ith == 0) {
    m_M0   = real_t(M0);
    m_mq   = real_t(mq);
    m_Ns   = Ns;
    m_NinF = m_Nvcd * m_Ns;
    m_alpha = alpha;

    assert(bc.size() == m_Ndim);
    if (m_boundary.size() != m_Ndim) m_boundary.resize(m_Ndim);

    for (int mu = 0; mu < m_Ndim; ++mu) {
      m_boundary[mu] = bc[mu];
      m_bc[mu]  = 1;
      if(do_comm[mu] > 0){ // do communication
        if(Communicator::ipe(mu) == 0) m_bc[mu] = m_boundary[mu];
        m_bc2[mu] = 0;
      }else{  // no communication
        m_bc[mu]  = 0;    // for boundar part (dummy)
        m_bc2[mu] = m_boundary[mu];  // for bulk part
      }
    }

    if (m_b.size() != m_Ns) {
      m_b.resize(m_Ns);
      m_c.resize(m_Ns);
    }
    for (int is = 0; is < m_Ns; ++is) {
      m_b[is] = real_t(b);
      m_c[is] = real_t(c);
    }
  }

#pragma omp barrier

  vout.general(m_vl, "%s: input parameters\n", class_name.c_str());
  vout.general(m_vl, "  gamma matrix repr.:  %s\n", m_repr.c_str());
  vout.general(m_vl, "  code implementation: %s\n", m_impl.c_str());
  vout.general(m_vl, "  mq   = %8.4f\n", m_mq);
  vout.general(m_vl, "  M0   = %8.4f\n", m_M0);
  vout.general(m_vl, "  Ns   = %4d\n", m_Ns);
  for (int mu = 0; mu < m_Ndim; ++mu) {
    vout.general(m_vl, "  boundary[%d] = %2d\n", mu, m_boundary[mu]);
  }
  vout.general(m_vl, "  coefficients:\n");
  for (int is = 0; is < m_Ns; ++is) {
    vout.general(m_vl, "    b[%2d] = %16.10f  c[%2d] = %16.10f\n",
                 is, m_b[is], is, m_c[is]);
  }
  vout.general(m_vl, "  alpha = %8.4f\n", m_alpha);

  // working 5d vectors.
  if (m_w1.nex() != Ns) {
    m_w1.reset(m_NinF, m_Nst2, 1);
    m_v1.reset(m_NinF, m_Nst2, 1);
    m_v2.reset(m_NinF, m_Nst2, 1);
  }

  set_precond_parameters();

#pragma omp barrier
}


//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_5din_eo<AFIELD>::set_coefficients(
                                     const std::vector<real_t> vec_b,
                                     const std::vector<real_t> vec_c)
{
#pragma omp barrier

  if ((vec_b.size() != m_Ns) || (vec_c.size() != m_Ns)) {
    vout.crucial(m_vl, "%s: size of coefficient vectors incorrect.\n",
                 class_name.c_str());
  }

  vout.general(m_vl, "%s: coefficient vectors are set:\n",
               class_name.c_str());

  int ith = ThreadManager::get_thread_id();
  if (ith == 0) {
    for (int is = 0; is < m_Ns; ++is) {
      m_b[is] = vec_b[is];
      m_c[is] = vec_c[is];
      vout.general(m_vl, "b[%2d] = %16.10f  c[%2d] = %16.10f\n",
                   is, m_b[is], is, m_c[is]);
    }
  }

  set_precond_parameters();

#pragma omp barrier
}

//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_5din_eo<AFIELD>::set_precond_parameters()
{
  int ith = ThreadManager::get_thread_id();
  if (ith == 0) {

    if (m_dp.size() != m_Ns) {
      m_dp.resize(m_Ns);
      m_dm.resize(m_Ns);
      m_dpinv.resize(m_Ns);
      m_e.resize(m_Ns - 1);
      m_f.resize(m_Ns - 1);
    }

    for (int is = 0; is < m_Ns; ++is) {
      m_dp[is] = m_alpha * (1.0 + m_b[is] * (4.0 - m_M0));
      m_dm[is] = m_alpha * (1.0 - m_c[is] * (4.0 - m_M0));
    }

    m_e[0] = m_mq * m_dm[m_Ns - 1] / m_dp[0];
    m_f[0] = m_mq * m_dm[0]/m_alpha;
    for (int is = 1; is < m_Ns - 1; ++is) {
      m_e[is] = m_e[is - 1] * m_dm[is - 1] / m_dp[is];
      m_f[is] = m_f[is - 1] * m_dm[is] / m_dp[is - 1];
    }

    m_g = m_e[m_Ns - 2] * m_dm[m_Ns - 2];

    for (int is = 0; is < m_Ns - 1; ++is) {
      m_dpinv[is] = 1.0 / m_dp[is];
    }
    m_dpinv[m_Ns - 1] = 1.0 / (m_dp[m_Ns - 1] + m_g);

  }

  // setting inverse matrix
  set_matrix5d_inverse();

}

//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_5din_eo<AFIELD>::set_matrix5d_inverse()
{
  int Nc = CommonParameters::Nc();

  int mat_size = m_Nd2 * m_Ns;
  m_mat_inv.resize(mat_size * mat_size);

  for(int isb = 0; isb < m_Ns; ++isb){
    for(int idb = 0; idb < m_Nd2; ++idb){
      m_w1.set(0.0);
      int inb  = 2 * Nc * (idb * 2 + m_Nd * isb);
      int site = 0;
      int idxb = m_index_eo.idxh(inb, m_NinF, site, 0);
      m_w1.set(idxb, 1.0);

      D_ee_inv(m_v2, m_w1, 0);

      for(int isa = 0; isa < m_Ns; ++isa){
        for(int ida = 0; ida < m_Nd2; ++ida){
          int ina   = 2 * Nc * (ida * 2 + m_Nd * isa);
          int idxa = m_index_eo.idxh(ina, m_NinF, site, 0);
          m_mat_inv[mat_index(ida, isa, idb, isb)] = m_v2.cmp(idxa);
	}
      }

    }
  }

}


//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_5din_eo<AFIELD>::set_config(Field *u)
{
  int nth = ThreadManager::get_num_threads();

  vout.detailed(m_vl, "%s: set_config is called: num_threads = %d\n",
                class_name.c_str(), nth);

  if (nth > 1) {
    set_config_impl(u);
  } else {
    set_config_omp(u);
  }

  vout.detailed(m_vl, "%s: set_config finished\n", class_name.c_str());
}

//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_5din_eo<AFIELD>::set_config_omp(Field *u)
{
  vout.detailed(m_vl, "  set_config_omp is called.\n");

#pragma omp parallel
  {
    set_config_impl(u);
  }
}

//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_5din_eo<AFIELD>::set_config_impl(Field *u)
{
#pragma omp barrier

  convert_gauge(m_index_eo, m_Ueo, *u);

#pragma omp barrier
}

//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_5din_eo<AFIELD>::convert(AFIELD& v, const Field& w)
{
#pragma omp barrier

  int Nst = v.nvol();
  if(w.nvol() != Nst){
    vout.crucial(m_vl, "%s: convert: field size irrelevant\n");
    exit(EXIT_FAILURE);
  }

  int ith, nth, isite, nsite;
  set_threadtask_afopr(ith, nth, isite, nsite, Nst);

  AIndex_lex<real_t, AFIELD::IMPL> index;

  for (int site = isite; site < nsite; ++site) {
    for (int is = 0; is < m_Ns; ++is) {
      for (int ivcd = 0; ivcd < NVCD; ++ivcd) {
         int in_alt = ivcd + NVCD * is;
         real_t vt = real_t(w.cmp(ivcd, site, is));
         v.set_host(index.idx(in_alt, m_NinF, site, 0), vt);
      }
    }
  } // site-llop

  v.update_device();

}


//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_5din_eo<AFIELD>::reverse(Field& v, const AFIELD& w)
{
#pragma omp barrier

  int Nst = v.nvol();
  if(w.nvol() != Nst){
    vout.crucial(m_vl, "%s: convert: field size irrelevant\n");
    exit(EXIT_FAILURE);
  }

  int ith, nth, isite, nsite;
  set_threadtask_afopr(ith, nth, isite, nsite, Nst);

  w.update_host();

  AIndex_lex<real_t, AFIELD::IMPL> index;

  for (int site = isite; site < nsite; ++site) {
    for (int is = 0; is < m_Ns; ++is) {
      for (int ivcd = 0; ivcd < NVCD; ++ivcd) {
        int in_alt = ivcd + NVCD * is;
        double vt = double(w.cmp_host(index.idx(in_alt, m_NinF, site, 0)));
        v.set(ivcd, site, is, vt);
      }
    }
  } // site-loop

#pragma omp barrier
}


//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_5din_eo<AFIELD>::set_mode(std::string mode)
{
#pragma omp barrier

  int ith = ThreadManager::get_thread_id();
  if (ith == 0) m_mode = mode;
  vout.paranoiac(m_vl, "  mode is set to %s\n", mode.c_str());

#pragma omp barrier
}


//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_5din_eo<AFIELD>::mult(AFIELD& v, const AFIELD& w)
{
  if (m_mode == "D") {
    D(v, w);
    //D_alt(v, w);
  } else if (m_mode == "Ddag") {
    Ddag(v, w);
    //Ddag_alt(v, w);
  } else if (m_mode == "DdagD") {
    DdagD(v, w);
    // DdagD_alt(v, w);
  //} else if (m_mode == "DDdag") { // this mode is not available
    // Ddag(m_w1, w);  // this fails because m_w1 is used in D.
    // D(v, m_w1);
  } else if (m_mode == "H") {
    H(v, w);
  }else if (m_mode == "Deo") {
    D_eo(v, w, 0);
  } else if (m_mode == "Doe") {
    D_eo(v, w, 1);
  } else if (m_mode == "Dee") {
    D_ee(v, w, 0);
  } else if (m_mode == "Doo") {
    D_ee(v, w, 1);
  } else if (m_mode == "Dee_inv") {
    D_ee_inv(v, w, 0);
    //D_ee_inv_alt(v, w, 0);
  } else if (m_mode == "Doo_inv") {
    D_ee_inv(v, w, 1);
    //D_ee_inv_alt(v, w, 1);
  } else {
    vout.crucial(m_vl, "mode undeifined in %s.\n", class_name.c_str());
    vout.crucial(m_vl, "in mult, mode=%s.\n", m_mode.c_str());
    exit(EXIT_FAILURE);
  }
}


//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_5din_eo<AFIELD>::mult_dag(AFIELD& v, const AFIELD& w)
{
  if (m_mode == "D") {
    Ddag(v, w);
    //Ddag_alt(v, w);
  } else if (m_mode == "Ddag") {
    D(v, w);
    //D_alt(v, w);
  } else if (m_mode == "DdagD") {
    DdagD(v, w);
    //DdagD_alt(v, w);
  } else if (m_mode == "H") {
    Hdag(v, w);
  } else if (m_mode == "Deo") {
    Ddag_eo(v, w, 1);
  } else if (m_mode == "Doe") {
    Ddag_eo(v, w, 0);
  } else if (m_mode == "Dee") {
    Ddag_ee(v, w, 0);
  } else if (m_mode == "Doo") {
    Ddag_ee(v, w, 1);
  } else if (m_mode == "Dee_inv") {
    Ddag_ee_inv(v, w, 0);
    //Ddag_ee_inv_alt(v, w, 0);
  } else if (m_mode == "Doo_inv") {
    Ddag_ee_inv(v, w, 1);
    //Ddag_ee_inv_alt(v, w, 1);
  } else {
    vout.crucial(m_vl, "mode undeifined in %s.\n", class_name.c_str());
    vout.crucial(m_vl, "in mult_dag, mode=%s.\n", m_mode.c_str());
    exit(EXIT_FAILURE);
  }
}


//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_5din_eo<AFIELD>::mult(AFIELD& v, const AFIELD& w,
                                       std::string mode)
{
  assert(w.check_size(m_NinF, m_Nst2, 1));
  assert(v.check_size(m_NinF, m_Nst2, 1));

  if (mode == "Deo") {
    D_eo(v, w, 0);
  } else if (mode == "Doe") {
    D_eo(v, w, 1);
  } else if (mode == "Dee") {
    D_ee(v, w, 0);
  } else if (mode == "Doo") {
    D_ee(v, w, 1);
  } else if (mode == "Dee_inv") {
    D_ee_inv(v, w, 0);
    //D_ee_inv_alt(v, w, 0);
  } else if (mode == "Doo_inv") {
    D_ee_inv(v, w, 1);
    //D_ee_inv_alt(v, w, 1);
  } else {
    vout.crucial(m_vl, "mode undeifined in %s.\n", class_name.c_str());
    vout.crucial(m_vl, "in mult, mode=%s.\n", m_mode.c_str());
    exit(EXIT_FAILURE);
  }

}

//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_5din_eo<AFIELD>::mult_dag(AFIELD& v,
                                                const AFIELD& w,
                                                std::string mode)
{
  assert(w.check_size(m_NinF, m_Nst2, 1));
  assert(v.check_size(m_NinF, m_Nst2, 1));

  if (mode == "Deo") {
    Ddag_eo(v, w, 1);
  } else if (mode == "Doe") {
    Ddag_eo(v, w, 0);
  } else if (mode == "Dee") {
    Ddag_ee(v, w, 0);
  } else if (mode == "Doo") {
    Ddag_ee(v, w, 1);
  } else if (mode == "Dee_inv") {
    Ddag_ee_inv(v, w, 0);
    //Ddag_ee_inv_alt(v, w, 0);
  } else if (mode == "Doo_inv") {
    Ddag_ee_inv(v, w, 1);
    //Ddag_ee_inv_alt(v, w, 1);
  } else {
    vout.crucial(m_vl, "mode undeifined in %s.\n", class_name.c_str());
    vout.crucial(m_vl, "in mult, mode=%s.\n", m_mode.c_str());
    exit(EXIT_FAILURE);
  }
}

//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_5din_eo<AFIELD>::DdagD(AFIELD& v, const AFIELD& w)
{
  D_eo(m_v1, w, 1);
  LU_inv(m_v2, m_v1);
  D_eo(m_v1, m_v2, 0);
  LU_inv(m_v2, m_v1);

  copy(m_v1, w);
  axpy(m_v1, real_t(-1.0), m_v2);

  LUdag_inv(v, m_v1);
  Ddag_eo(m_v2, v, 1);
  LUdag_inv(v, m_v2);
  Ddag_eo(m_v2, v, 0);

  copy(v, m_v1);
  axpy(v, real_t(-1.0), m_v2);
}


//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_5din_eo<AFIELD>::D(AFIELD& v, const AFIELD& w)
{
  D_eo(m_v1, w, 1);
  LU_inv(m_v2, m_v1);
  D_eo(m_v1, m_v2, 0);
  LU_inv(m_v2, m_v1);

  copy(v, w);
  axpy(v, real_t(-1.0), m_v2);
}



//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_5din_eo<AFIELD>::H(AFIELD& v, const AFIELD& w)
{
#pragma omp barrier

  D_eo(m_v1, w, 1);
  LU_inv(m_v2, m_v1);
  D_eo(m_v1, m_v2, 0);
  LU_inv(m_v2, m_v1);

  copy(m_v1, w);
  axpy(m_v1, real_t(-1.0), m_v2);
  mult_gm5R(v, m_v1);
#pragma omp barrier
}



//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_5din_eo<AFIELD>::Ddag(AFIELD& v, const AFIELD& w)
{
  assert(w.check_size(m_NinF, m_Nst2, 1));
  assert(v.check_size(m_NinF, m_Nst2, 1));

#pragma omp barrier

  LUdag_inv(m_v1, w);
  Ddag_eo(m_v2, m_v1, 1);
  LUdag_inv(m_v1, m_v2);
  Ddag_eo(m_v2, m_v1, 0);

  copy(v, w);
  axpy(v, real_t(-1.0), m_v2);

}

//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_5din_eo<AFIELD>::Hdag(AFIELD& v, const AFIELD& w)
{
#pragma omp barrier
  mult_gm5R(v, w);

  LUdag_inv(m_v1, v);
  Ddag_eo(m_v2, m_v1, 1);
  LUdag_inv(m_v1, m_v2);
  Ddag_eo(m_v2, m_v1, 0);

  axpy(v, real_t(-1.0), m_v2);
#pragma omp barrier
}


//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_5din_eo<AFIELD>::DdagD_alt(AFIELD& v, const AFIELD& w)
{
  D_eo(m_v1, w, 1);
  D_ee_inv_alt(m_v2, m_v1, 1);
  D_eo(m_v1, m_v2, 0);
  D_ee_inv_alt(m_v2, m_v1, 0);

  copy(m_v1, w);
  axpy(m_v1, real_t(-1.0), m_v2);

  Ddag_ee_inv_alt(v, m_v1, 0);
  Ddag_eo(m_v2, v, 1);
  Ddag_ee_inv_alt(v, m_v2, 1);
  Ddag_eo(m_v2, v, 0);

  copy(v, m_v1);
  axpy(v, real_t(-1.0), m_v2);
}

//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_5din_eo<AFIELD>::D_alt(AFIELD& v, const AFIELD& w)
{
  D_eo(m_v1, w, 1);
  D_ee_inv_alt(m_v2, m_v1, 1);
  D_eo(m_v1, m_v2, 0);
  D_ee_inv_alt(m_v2, m_v1, 0);

  copy(v, w);
  axpy(v, real_t(-1.0), m_v2);
}


//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_5din_eo<AFIELD>::Ddag_alt(AFIELD& v, const AFIELD& w)
{
  assert(w.check_size(m_NinF, m_Nst2, 1));
  assert(v.check_size(m_NinF, m_Nst2, 1));

#pragma omp barrier

  Ddag_ee_inv_alt(m_v1, w, 0);
  Ddag_eo(m_v2, m_v1, 1);
  Ddag_ee_inv_alt(m_v1, m_v2, 1);
  Ddag_eo(m_v2, m_v1, 0);

  copy(v, w);
  axpy(v, real_t(-1.0), m_v2);

}


//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_5din_eo<AFIELD>::mult_gm5(AFIELD& v,
                                                const AFIELD& w)
{
  real_t *vp = v.ptr(0);
  real_t *wp = const_cast<AFIELD *>(&w)->ptr(0);

#pragma omp barrier

  int ith = ThreadManager::get_thread_id();

  if (ith == 0)
    BridgeACC::mult_domainwall_5din_mult_gm5_dirac(
                                              vp, wp, m_Ns, m_Nsize);

#pragma omp barrier
}

//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_5din_eo<AFIELD>::mult_R(AFIELD& v,
                                              const AFIELD& w)
{
  assert(w.check_size(m_NinF, m_Nst2, 1));
  assert(v.check_size(m_NinF, m_Nst2, 1));

#pragma omp barrier

  real_t *vp = v.ptr(0);
  real_t *wp = const_cast<AFIELD *>(&w)->ptr(0);

  int ith = ThreadManager::get_thread_id();

  if (ith == 0)
    BridgeACC::mult_domainwall_5din_mult_R(vp, wp, m_Ns, m_Nsize);

#pragma omp barrier
}


//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_5din_eo<AFIELD>::mult_gm5R(AFIELD& v,
                                                 const AFIELD& w)
{
#pragma omp barrier

  real_t *vp = v.ptr(0);
  real_t *wp = const_cast<AFIELD *>(&w)->ptr(0);

  int ith = ThreadManager::get_thread_id();

  if (ith == 0)
    BridgeACC::mult_domainwall_5din_mult_gm5R_dirac(
                                              vp, wp, m_Ns, m_Nsize);

#pragma omp barrier
}


//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_5din_eo<AFIELD>::D_eo(AFIELD& v,
                                            const AFIELD& w,
                                            const int ieo)
{  // ieo = 0: even < odd, ieo = 1: odd <- even
#pragma omp barrier

  int ith = ThreadManager::get_thread_id();
  if(ith == 0) START_TIMER(timer_mult_Deo);

  int nth = ThreadManager::get_num_threads();
  int ith_kernel = 0;
  if(nth > 1) ith_kernel = 1;

  real_t *vp = v.ptr(0);
  real_t *wp = const_cast<AFIELD *>(&w)->ptr(0);
  real_t *yp = m_w1.ptr(0);  // working vector
  real_t *up = m_Ueo.ptr(0);
  int jeo    = (m_Ieo_origin + ieo) % 2;

  if(do_comm_any > 0 && ith == 0){
    START_TIMER(timer_comm_recv_start);
    chset_recv.start();
    STOP_TIMER(timer_comm_recv_start);
  }

  if (ith == ith_kernel){

    BridgeACC::mult_domainwall_5din_eo_5dir_dirac(
                              yp, wp, m_mq, m_M0, m_Ns,
                              &m_b[0], &m_c[0], m_alpha, m_Nsize);

    if (do_comm_any > 0) {
      START_TIMER(timer_pack);
      real_t *buf1_xp = (real_t *)chsend_dn[0].ptr();
      real_t *buf1_xm = (real_t *)chsend_up[0].ptr();
      real_t *buf1_yp = (real_t *)chsend_dn[1].ptr();
      real_t *buf1_ym = (real_t *)chsend_up[1].ptr();
      real_t *buf1_zp = (real_t *)chsend_dn[2].ptr();
      real_t *buf1_zm = (real_t *)chsend_up[2].ptr();
      real_t *buf1_tp = (real_t *)chsend_dn[3].ptr();
      real_t *buf1_tm = (real_t *)chsend_up[3].ptr();

      BridgeACC::mult_domainwall_5din_eo_hop1_dirac(
                       buf1_xp, buf1_xm, buf1_yp, buf1_ym,
                       buf1_zp, buf1_zm, buf1_tp, buf1_tm,
                       up, yp,
                       m_Ns, m_bc, m_Nsize, do_comm, ieo, jeo, 0);
      STOP_TIMER(timer_pack);
    }
  }

#pragma omp barrier

  if(do_comm_any > 0 && ith == 0){
    START_TIMER(timer_comm);
    START_TIMER(timer_comm_send_start);
    chset_send.start();
    STOP_TIMER(timer_comm_send_start);
  }

  if (ith == ith_kernel){
    START_TIMER(timer_bulk);
    if(m_impl == "5d"){
      BridgeACC::mult_domainwall_5din_eo_hopb_dirac_5d(
                         vp, up, yp, m_Ns, m_bc2,
                         m_Nsize, do_comm, ieo, jeo, 0);
    }else{
#ifdef USE_DOMAINWALL_5DIN_4D_KERNEL
      BridgeACC::mult_domainwall_5din_eo_hopb_dirac_4d(
                         vp, up, yp, m_Ns, m_bc2,
                         m_Nsize, do_comm, ieo, jeo, 0);
#else
      vout.crucial(m_vl, "%s: 4D kernel not compiled\n",
                   class_name.c_str());
      exit(EXIT_FAILURE);
#endif
    }
    STOP_TIMER(timer_bulk);
  }

  if(do_comm_any > 0 && ith == 0){
    START_TIMER(timer_comm_recv_wait);
    chset_recv.wait();
    STOP_TIMER(timer_comm_recv_wait);
  }

#pragma omp barrier

  if(do_comm_any > 0 && ith == 0){
    STOP_TIMER(timer_comm);
  }

  if(do_comm_any > 0 && ith == ith_kernel){
    START_TIMER(timer_boundary);

    real_t *buf2_xp = (real_t *)chrecv_up[0].ptr();
    real_t *buf2_xm = (real_t *)chrecv_dn[0].ptr();
    real_t *buf2_yp = (real_t *)chrecv_up[1].ptr();
    real_t *buf2_ym = (real_t *)chrecv_dn[1].ptr();
    real_t *buf2_zp = (real_t *)chrecv_up[2].ptr();
    real_t *buf2_zm = (real_t *)chrecv_dn[2].ptr();
    real_t *buf2_tp = (real_t *)chrecv_up[3].ptr();
    real_t *buf2_tm = (real_t *)chrecv_dn[3].ptr();

    BridgeACC::mult_domainwall_5din_eo_hop2_dirac(
                        vp, up, yp,
                        buf2_xp, buf2_xm, buf2_yp, buf2_ym,
                        buf2_zp, buf2_zm, buf2_tp, buf2_tm,
                        m_Ns, m_bc,
                        m_Nsize, do_comm, ieo, jeo);

    STOP_TIMER(timer_boundary);
  }
  if(do_comm_any > 0 && ith == 0){
    START_TIMER(timer_comm_send_wait);
    chset_send.wait();
    STOP_TIMER(timer_comm_send_wait);
  }

#pragma omp barrier

  if(ith == 0) STOP_TIMER(timer_mult_Deo);

}

//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_5din_eo<AFIELD>::Ddag_eo(AFIELD& v,
                                              const AFIELD& w,
                                              const int ieo)
  // ieo = 0: even < odd, ieo = 1: odd <- even
{
#pragma omp barrier

  int ith = ThreadManager::get_thread_id();
  if(ith == 0) START_TIMER(timer_mult_Deo);

  int nth = ThreadManager::get_num_threads();
  int ith_kernel = 0;
  if(nth > 1) ith_kernel = 1;

  real_t *vp = v.ptr(0);
  real_t *wp = const_cast<AFIELD *>(&w)->ptr(0);
  real_t *yp = m_w1.ptr(0);
  real_t *up = m_Ueo.ptr(0);
  int jeo    = (m_Ieo_origin + ieo) % 2;

  if(do_comm_any > 0 && ith == 0){
    START_TIMER(timer_comm_recv_start);
    chset_recv.start();
    STOP_TIMER(timer_comm_recv_start);
  }

  if (ith == ith_kernel){

    if (do_comm_any > 0) {
      START_TIMER(timer_pack);

      real_t *buf1_xp = (real_t *)chsend_dn[0].ptr();
      real_t *buf1_xm = (real_t *)chsend_up[0].ptr();
      real_t *buf1_yp = (real_t *)chsend_dn[1].ptr();
      real_t *buf1_ym = (real_t *)chsend_up[1].ptr();
      real_t *buf1_zp = (real_t *)chsend_dn[2].ptr();
      real_t *buf1_zm = (real_t *)chsend_up[2].ptr();
      real_t *buf1_tp = (real_t *)chsend_dn[3].ptr();
      real_t *buf1_tm = (real_t *)chsend_up[3].ptr();

      BridgeACC::mult_domainwall_5din_eo_hop1_dirac(
                         buf1_xp, buf1_xm, buf1_yp, buf1_ym,
                         buf1_zp, buf1_zm, buf1_tp, buf1_tm,
                         up, wp,
                         m_Ns, m_bc, m_Nsize, do_comm, ieo, jeo, 1);
      STOP_TIMER(timer_pack);
    }
  }
#pragma omp barrier

  if(do_comm_any > 0 && ith == 0){
    START_TIMER(timer_comm);
    START_TIMER(timer_comm_send_start);
    chset_send.start();
    STOP_TIMER(timer_comm_send_start);
  }

  if(ith == ith_kernel){
    START_TIMER(timer_bulk);
    if(m_impl == "5d"){
      BridgeACC::mult_domainwall_5din_eo_hopb_dirac_5d(
                           yp, up, wp, m_Ns, m_bc2,
                           m_Nsize, do_comm, ieo, jeo, 1);
    }else{
#ifdef USE_DOMAINWALL_5DIN_4D_KERNEL
      BridgeACC::mult_domainwall_5din_eo_hopb_dirac_4d(
                           yp, up, wp, m_Ns, m_bc2,
                           m_Nsize, do_comm, ieo, jeo, 1);
#else
      vout.crucial(m_vl, "%s: 4D kernel not compiled\n",
                   class_name.c_str());
      exit(EXIT_FAILURE);
#endif
    }
    STOP_TIMER(timer_bulk);
  }

  if(do_comm_any > 0 && ith == 0){
    START_TIMER(timer_comm_recv_wait);
    chset_recv.wait();
    STOP_TIMER(timer_comm_recv_wait);
  }

#pragma omp barrier

  if(do_comm_any > 0 && ith == 0){
    STOP_TIMER(timer_comm);
  }

  if (ith == ith_kernel) {
    if (do_comm_any > 0) {
      START_TIMER(timer_boundary);

      real_t *buf2_xp = (real_t *)chrecv_up[0].ptr();
      real_t *buf2_xm = (real_t *)chrecv_dn[0].ptr();
      real_t *buf2_yp = (real_t *)chrecv_up[1].ptr();
      real_t *buf2_ym = (real_t *)chrecv_dn[1].ptr();
      real_t *buf2_zp = (real_t *)chrecv_up[2].ptr();
      real_t *buf2_zm = (real_t *)chrecv_dn[2].ptr();
      real_t *buf2_tp = (real_t *)chrecv_up[3].ptr();
      real_t *buf2_tm = (real_t *)chrecv_dn[3].ptr();

      BridgeACC::mult_domainwall_5din_eo_hop2_dirac(
                            yp, up, vp,
                            buf2_xp, buf2_xm, buf2_yp, buf2_ym,
                            buf2_zp, buf2_zm, buf2_tp, buf2_tm,
                            m_Ns, m_bc,
                            m_Nsize, do_comm, ieo, jeo);
      STOP_TIMER(timer_boundary);

    }
    BridgeACC::mult_domainwall_5din_eo_5dirdag_dirac(
                            vp, yp, m_mq, m_M0, m_Ns,
                            &m_b[0], &m_c[0], m_alpha, m_Nsize);

  }
  if(do_comm_any > 0 && ith == 0){
    START_TIMER(timer_comm_send_wait);
    chset_send.wait();
    STOP_TIMER(timer_comm_send_wait);
  }

#pragma omp barrier

  if(ith == 0) STOP_TIMER(timer_mult_Deo);

}

//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_5din_eo<AFIELD>::D_ee(AFIELD& v, const AFIELD& w,
                                       const int ieo)
{
#pragma omp barrier

  int ith = ThreadManager::get_thread_id();
  if (ith == 0){
    real_t *vp = v.ptr(0);
    real_t *wp = const_cast<AFIELD *>(&w)->ptr(0);

    BridgeACC::mult_domainwall_5din_ee_5dir_dirac(
                                 vp, wp, m_mq, m_M0, m_Ns,
                                 &m_b[0], &m_c[0], m_alpha, m_Nsize);
  }

#pragma omp barrier

}


//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_5din_eo<AFIELD>::Ddag_ee(AFIELD& v, const AFIELD& w,
                                       const int ieo)
{
#pragma omp barrier

  int ith = ThreadManager::get_thread_id();
  if (ith == 0){
    real_t *vp = v.ptr(0);
    real_t *wp = const_cast<AFIELD *>(&w)->ptr(0);

    BridgeACC::mult_domainwall_5din_ee_5dirdag_dirac(
                                 vp, wp, m_mq, m_M0, m_Ns,
                                 &m_b[0], &m_c[0], m_alpha, m_Nsize);
  }

#pragma omp barrier


}

//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_5din_eo<AFIELD>::D_ee_inv_alt(AFIELD& v,
                                                    const AFIELD& w,
                                                    const int ieo)
{
#pragma omp barrier

  int ith = ThreadManager::get_thread_id();

  if(ith == 0) START_TIMER(timer_mult_Dee_inv);

#ifdef USE_DOMAINWALL_5DIN_EE_MATINV_KERNEL
  if (ith == 0){
    real_t *vp = v.ptr(0);
    real_t *wp = const_cast<AFIELD *>(&w)->ptr(0);

    if(m_impl == "5d"){
      BridgeACC::mult_domainwall_5din_ee_inv_dirac_5d(
                      vp, wp, 1, m_Ns, &m_mat_inv[0], m_Nsize);
    }else{
      BridgeACC::mult_domainwall_5din_ee_inv_dirac_4d(
                      vp, wp, 1, m_Ns, &m_mat_inv[0], m_Nsize);
    }
  }
#else
  vout.crucial(m_vl, "%s: Dee inverse matrix mult not compiled\n",
               class_name.c_str());
  exit(EXIT_FAILURE);
#endif

#pragma omp barrier

  if(ith == 0) STOP_TIMER(timer_mult_Dee_inv);

}

//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_5din_eo<AFIELD>::Ddag_ee_inv_alt(AFIELD& v,
                                                       const AFIELD& w,
                                                       const int ieo)
{
#pragma omp barrier

  int ith = ThreadManager::get_thread_id();

  if(ith == 0) START_TIMER(timer_mult_Dee_inv);

#ifdef USE_DOMAINWALL_5DIN_EE_MATINV_KERNEL
  if (ith == 0){
    real_t *vp = v.ptr(0);
    real_t *wp = const_cast<AFIELD *>(&w)->ptr(0);

    if(m_impl == "5d"){
      BridgeACC::mult_domainwall_5din_ee_inv_dirac_5d(
                      vp, wp, -1, m_Ns, &m_mat_inv[0], m_Nsize);
    }else{
      BridgeACC::mult_domainwall_5din_ee_inv_dirac_4d(
                      vp, wp, -1, m_Ns, &m_mat_inv[0], m_Nsize);
    }
  }
#else
  vout.crucial(m_vl, "%s: Dee inverse matrix mult not compiled\n",
               class_name.c_str());
  exit(EXIT_FAILURE);
#endif

#pragma omp barrier

  if(ith == 0) STOP_TIMER(timer_mult_Dee_inv);

}

//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_5din_eo<AFIELD>::D_ee_inv(AFIELD& v, const AFIELD& w,
                                       const int ieo)
{
#pragma omp barrier

  LU_inv(v, w);

#pragma omp barrier

}

//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_5din_eo<AFIELD>::Ddag_ee_inv(AFIELD& v, const AFIELD& w,
                                       const int ieo)
{
#pragma omp barrier

  LUdag_inv(v, w);

#pragma omp barrier

}

//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_5din_eo<AFIELD>::LU_inv(AFIELD& v, const AFIELD& w)
{

  int ith = ThreadManager::get_thread_id();

  if (ith == 0){

    START_TIMER(timer_mult_Dee_inv);

    real_t *vp = v.ptr(0);
    real_t *wp = const_cast<AFIELD *>(&w)->ptr(0);

    BridgeACC::mult_domainwall_5din_LUinv_dirac(
                   vp, wp, m_Ns, m_Nsize,
                   &m_e[0], &m_f[0], &m_dpinv[0], &m_dm[0], m_alpha);

  }
#pragma omp barrier

  if(ith == 0) STOP_TIMER(timer_mult_Dee_inv);

}


//====================================================================
template<typename AFIELD>
void AFopr_Domainwall_5din_eo<AFIELD>::LUdag_inv(AFIELD& v, const AFIELD& w)
{
  int ith = ThreadManager::get_thread_id();

  if (ith == 0){

    START_TIMER(timer_mult_Dee_inv);

    real_t *vp = v.ptr(0);
    real_t *wp = const_cast<AFIELD *>(&w)->ptr(0);

    BridgeACC::mult_domainwall_5din_LUdaginv_dirac(
                  vp, wp, m_Ns, m_Nsize,
                  &m_e[0], &m_f[0], &m_dpinv[0], &m_dm[0], m_alpha);

  }
#pragma omp barrier

  if(ith == 0) STOP_TIMER(timer_mult_Dee_inv);

}


//====================================================================
template<typename AFIELD>
double AFopr_Domainwall_5din_eo<AFIELD>::flop_count(std::string mode)
{
  // flop counting of this class is not confirmed.

  int    Lvol  = CommonParameters::Lvol();
  double vsite = static_cast<double>(Lvol);
  double vNs   = static_cast<double>(m_Ns);

  double flop_Wilson_hop;
  double flop_LU_inv;
  double flop_pre;
  int    Nc = CommonParameters::Nc();
  int    Nd = CommonParameters::Nd();
  if (m_repr == "Dirac") {
    flop_Wilson_hop = static_cast<double>(
      vNs * Nc * Nd * (6 * (4 * Nc + 2) + 2 * (4 * Nc + 1))) * vsite;
    flop_pre    = static_cast<double>(vNs * Nc * Nd * 14) * vsite;
    // chiral projectin:
    //   FLOP as (1 +- gm5), because the factor (1/2) can be absorbed
    //   into the other coefficients
    //
    // L_inv: chiral projection: Ns           [each: 2*Nc*Nd]
    //        axpy:              2(Ns-1)      [each: 4*Nc*Nd]
    // U_inv: chiral projection: Ns
    //        axpy:              2(Ns-1)
    //        scal:              Ns           [each: 2*Nc*Nd]
    flop_LU_inv = static_cast<double>(Nc * Nd * (22 * vNs -16)) * vsite;
  } else if (m_repr == "Chiral") {
    flop_Wilson_hop = static_cast<double>(
      vNs * Nc * Nd * (8 * (4 * Nc + 2))) * vsite;
    flop_pre    = static_cast<double>(vNs * Nc * Nd * 6) * vsite;
    // no FLOP is needed for chiral projection
    // L_inv: chiral projection: Ns
    //        axpy:              2(Ns-1)      [each: 2*Nc*Nd for chirality+/-]
    // U_inv: chiral projection: Ns
    //        axpy:              2(Ns-1)
    //        scal:              Ns           [each: 2*Nc*Nd]
    flop_LU_inv = static_cast<double>(Nc * Nd * (10 * vNs -8)) * vsite;
  }
  double flop_axpy = static_cast<double>(vNs * 2 * Nc * Nd) * vsite;

  //  double axpy1 = static_cast<double>(2 * m_NinF);
  //  double scal1 = static_cast<double>(1 * m_NinF);
  //  double flop_DW = vNs * (flop_Wilson + vsite*(6*axpy1 + 2*scal1));
  // In Ddag case, flop_Wilson + 7 axpy which equals flop_DW.
  double flop_Deo = flop_Wilson_hop + flop_pre;


  //  double flop_LU_inv = 2.0 * vsite *
  //            ((3.0*axpy1 + scal1)*(vNs-1.0) + axpy1 + 2.0*scal1);

  double flop = 0.0;
  if ((mode == "D") || (mode == "Ddag")) {
    flop = flop_Deo + flop_LU_inv + 0.5 * flop_axpy;  // axpy is for 1 - (Deo etc)
  } else if (mode == "DdagD") {
    flop = 2 * flop_Deo + 2 * flop_LU_inv + flop_axpy;
  } else if ((mode == "Dee_inv") || (mode == "Doo_inv")) {
    flop = 0.5 * flop_LU_inv; // 0.5 is for even-odd
  } else if ((mode == "Dee") || (mode == "Doo")) {
    flop = 0.5 * flop_pre;
  } else if ((mode == "Doe") || (mode == "Doe")) {
    flop = 0.5 * flop_Deo;  // 0.5 is for even-odd
  } else {
    vout.crucial(m_vl, "Error at %s: input mode is undefined: %s.\n",
                 class_name.c_str(), mode.c_str());
    exit(EXIT_FAILURE);
  }

  return flop;
}

//============================================================END=====
