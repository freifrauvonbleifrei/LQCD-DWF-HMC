/*!
        @file    aforceSmear_APE-tmpl.h

        @brief

        @author  Issaku Kanamori (kanamori)
                 $LastChangedBy: kanamori $

        @date    $LastChangedDate: 2025-12-03 19:35:35 +0900 (2025年12月03日 (水)) $

        @version $LastChangedRevision: 2668 $
*/

#include "lib_alt_Accel/Smear/aforceSmear_APE.h"

#ifdef USE_FACTORY_AUTOREGISTER
namespace {
  bool init = ForceSmear_APE::register_factory();
}
#endif

template<typename AFIELD>
const std::string AForceSmear_APE<AFIELD>::class_name
                                      = "AForceSmear_APE<AFIELD>";

//====================================================================
template<typename AFIELD>
void AForceSmear_APE<AFIELD>::init()
{
  m_Ndim = CommonParameters::Ndim();
  m_Nvol = CommonParameters::Nvol();

  m_rho.resize(m_Ndim * m_Ndim);
  m_U.resize(m_Ndim);
  m_iTheta.resize(m_Ndim);

  for (int mu = 0; mu < m_Ndim; ++mu) {
    for (int nu = 0; nu < m_Ndim; ++nu) {
      m_rho[mu + nu * m_Ndim] = 0.0;
    }
  }

  int Nc = CommonParameters::Nc();
  int Nin = 2 * Nc * Nc;
  m_shift = new ShiftField_lex(Nin);

  int Nvol = CommonParameters::Nvol();
  m_U_alt.resize(m_Ndim);
  m_iTheta_alt.resize(m_Ndim);
  for(int mu=0; mu<m_Ndim; ++mu){
    m_U_alt[mu].reset(Nin, Nvol, 1);
  }
  for(int mu=0; mu<m_Ndim; ++mu){
    m_iTheta_alt[mu].reset(Nin, Nvol, 1);
  }

  m_C_alt.reset(Nin, Nvol, 1);
  m_Ctmp_alt.reset(Nin, Nvol, 1);

  m_shift_alt.reset(new ShiftAField_lex<AFIELD>(Nin) );
  m_av1.reset(Nin, Nvol, 1);

  m_vt1_alt.reset(Nin, Nvol, 1);
  m_vt2_alt.reset(Nin, Nvol, 1);
  m_vt3_alt.reset(Nin, Nvol, 1);
  m_sigmap_tmp_alt.reset(Nin, Nvol, 1);
  m_sigma_tmp_alt.reset(Nin, Nvol, 1);
}


//====================================================================
template<typename AFIELD>
void AForceSmear_APE<AFIELD>::tidyup()
{
  delete m_shift;
}


//====================================================================
template<typename AFIELD>
void AForceSmear_APE<AFIELD>::set_parameters(const Parameters& params)
{
 std::string vlevel;
  if (!params.fetch_string("verbose_level", vlevel)) {
    m_vl = vout.set_verbose_level(vlevel);
  }

  //- fetch and check input parameters
  double rho1;

  int err = 0;
  err += params.fetch_double("rho_uniform", rho1);

  if (err) {
    vout.crucial(m_vl, "Error at %s: input parameter not found.\n",
                 class_name.c_str());
    exit(EXIT_FAILURE);
  }

  set_parameters(rho1);
}


//====================================================================
template<typename AFIELD>
void AForceSmear_APE<AFIELD>::get_parameters(Parameters& params) const
{
  params.set_double_vector("rho", m_rho);

  params.set_string("verbose_level", vout.get_verbose_level(m_vl));
}


//====================================================================
template<typename AFIELD>
void AForceSmear_APE<AFIELD>::set_parameters(const double rho1)
{
  //- print input parameters
  vout.general(m_vl, "%s:\n", class_name.c_str());
  vout.general(m_vl, "  rho = %8.4f\n", rho1);

  //- store values
  for (int mu = 0; mu < m_Ndim; ++mu) {
    for (int nu = 0; nu < m_Ndim; ++nu) {
      m_rho[mu + nu * m_Ndim] = rho1;
    }
  }
}


//====================================================================
template<typename AFIELD>
void AForceSmear_APE<AFIELD>::set_parameters(const std::vector<double>& rho)
{
  //- print input parameters
  vout.general(m_vl, "%s:\n", class_name.c_str());
  for (int mu = 0; mu < m_Ndim; ++mu) {
    vout.general(m_vl, "  rho[%d] = %8.4f\n", mu, rho[mu]);
  }

  //- range check
  // NB. rho == 0 is allowed.
  assert(rho.size() == m_Ndim * m_Ndim);

  //- store values
  for (int mu = 0; mu < m_Ndim; ++mu) {
    for (int nu = 0; nu < m_Ndim; ++nu) {
      m_rho[mu + nu * m_Ndim] = rho[mu + nu * m_Ndim];
    }
  }
}



//====================================================================
template<typename AFIELD>
void AForceSmear_APE<AFIELD>::force_udiv(Field_G& Sigma,
                                         const Field_G& Sigmap,
                                         const Field_G& U)
{
  const int Nc = CommonParameters::Nc();

  assert(Sigmap.nin() == (2 * Nc * Nc));
  assert(Sigmap.nvol() == m_Nvol);
  assert(Sigmap.nex() == m_Ndim);

  AIndex_lex<real_t, AFIELD::IMPL> index;

  for (int mu = 0; mu < m_Ndim; ++mu) {
    copy(m_U[mu], 0, U, mu);
  }
#pragma omp barrier

  for (int mu = 0; mu < m_Ndim; ++mu) {
    m_C.set(0.0);
#pragma omp barrier

    for (int nu = 0; nu < m_Ndim; ++nu) {
      if (nu == mu) continue;
      real_t rho = m_rho[mu + m_Ndim * nu];

      staple(m_Ctmp, m_U[mu], m_U[nu], mu, nu);

      axpy(m_C, 0, rho, m_Ctmp, 0);
#pragma omp barrier
    }

    copy(m_sigmap_tmp, 0, Sigmap, mu);
#pragma omp barrier

    real_t alpha = m_rho[mu + m_Ndim * mu];

    m_proj->force_recursive(m_Xi, m_iTheta[mu],
                            alpha, m_sigmap_tmp, m_C, m_U[mu]);
    copy(Sigma, mu, m_Xi, 0);
#pragma omp barrier
  }

  for (int mu = 0; mu < m_Ndim; ++mu) {
    for (int nu = 0; nu < m_Ndim; ++nu) {
      if (nu == mu) continue;
      real_t rho = m_rho[mu + m_Ndim * nu];

      force_each(m_sigma_tmp, m_U[mu], m_U[nu],
                 m_iTheta[mu], m_iTheta[nu], mu, nu);

      axpy(Sigma, mu, rho, m_sigma_tmp, 0);
#pragma omp barrier
    }
  }
}


//====================================================================
template<typename AFIELD>
void AForceSmear_APE<AFIELD>::force_udiv(AFIELD& Sigma,
                                         const AFIELD& Sigmap,
                                         const Field_G& U)
{
  // m_proj->force_recursive is not vectorized

  const int Nc = CommonParameters::Nc();

  assert(Sigmap.nin() == (2 * Nc * Nc));
  assert(Sigmap.nvol() == m_Nvol);
  assert(Sigmap.nex() == m_Ndim);

  AIndex_lex<real_t, AFIELD::IMPL> index;

  for (int mu = 0; mu < m_Ndim; ++mu) {
    copy(m_U[mu], 0, U, mu);
    convert_gauge(index, m_U_alt[mu], m_U[mu]);
  }
#pragma omp barrier

  for (int mu = 0; mu < m_Ndim; ++mu) {
    m_C_alt.set(0.0);
#pragma omp barrier

    for (int nu = 0; nu < m_Ndim; ++nu) {
      if (nu == mu) continue;
      real_t rho = m_rho[mu + m_Ndim * nu];

      staple(m_Ctmp_alt, m_U_alt[mu], m_U_alt[nu], mu, nu);
      axpy(m_C_alt, 0, rho, m_Ctmp_alt, 0);
#pragma omp barrier
    }

    copy(m_sigmap_tmp_alt, 0, Sigmap, mu);
#pragma omp barrier

    real_t alpha = m_rho[mu + m_Ndim * mu];

    reverse_gauge(index, m_sigmap_tmp, m_sigmap_tmp_alt);
    reverse_gauge(index, m_C, m_C_alt);
    m_proj->force_recursive(m_Xi, m_iTheta[mu],
                            alpha, m_sigmap_tmp, m_C, m_U[mu]);
    convert_gauge(index, m_av1, m_Xi);
    convert_gauge(index, m_iTheta_alt[mu], m_iTheta[mu]);
    copy(Sigma, mu, m_av1, 0);
#pragma omp barrier
  }

  for (int mu = 0; mu < m_Ndim; ++mu) {
    for (int nu = 0; nu < m_Ndim; ++nu) {
      if (nu == mu) continue;
      real_t rho = m_rho[mu + m_Ndim * nu];

      force_each(m_sigma_tmp_alt, m_U_alt[mu], m_U_alt[nu],
                 m_iTheta_alt[mu], m_iTheta_alt[nu], mu, nu);
      axpy(Sigma, mu, rho, m_sigma_tmp_alt, 0);
#pragma omp barrier
    }
  }
}


//====================================================================
template<typename AFIELD>
void AForceSmear_APE<AFIELD>::force_each(Field_G& Sigma_mu,
                                         const Field_G& V_mu,
                                         const Field_G& V_nu,
                                         const Field_G& iTheta_mu,
                                         const Field_G& iTheta_nu,
                                         const int mu, const int nu)
{
#pragma omp barrier

  Sigma_mu.set(0.0);
#pragma omp barrier

  m_shift->backward(m_vt1, V_nu, mu);
  m_shift->backward(m_vt2, V_mu, nu);

  mult_Field_Gnd(m_vt3, 0, m_vt1, 0, m_vt2, 0);
#pragma omp barrier

  multadd_Field_Gnd(Sigma_mu, 0, m_vt3, 0, iTheta_nu, 0, 1.0);
#pragma omp barrier

  mult_Field_Gdn(m_vt3, 0, iTheta_mu, 0, V_nu, 0);
#pragma omp barrier

  mult_Field_Gdn(m_vt2, 0, m_vt1, 0, m_vt3, 0);
#pragma omp barrier

  m_shift->forward(m_vt3, m_vt2, nu);
  axpy(Sigma_mu, 1.0, m_vt3);
#pragma omp barrier

  mult_Field_Gdn(m_vt3, 0, V_mu, 0, iTheta_nu, 0);
#pragma omp barrier

  mult_Field_Gdn(m_vt2, 0, m_vt1, 0, m_vt3, 0);
#pragma omp barrier

  m_shift->forward(m_vt3, m_vt2, nu);
  axpy(Sigma_mu, 1.0, m_vt3);
#pragma omp barrier

  m_shift->backward(m_vt1, iTheta_nu, mu);
  m_shift->backward(m_vt2, V_mu, nu);

  mult_Field_Gnd(m_vt3, 0, m_vt1, 0, m_vt2, 0);
#pragma omp barrier

  multadd_Field_Gnd(Sigma_mu, 0, m_vt3, 0, V_nu, 0, 1.0);
#pragma omp barrier

  mult_Field_Gdd(m_vt2, 0, m_vt1, 0, V_mu, 0);
#pragma omp barrier

  mult_Field_Gnn(m_vt3, 0, m_vt2, 0, V_nu, 0);
#pragma omp barrier

  m_shift->forward(m_vt2, m_vt3, nu);

  axpy(Sigma_mu, 1.0, m_vt2);
#pragma omp barrier

  m_shift->backward(m_vt1, V_nu, mu);
  m_shift->backward(m_vt2, iTheta_mu, nu);

  mult_Field_Gnd(m_vt3, 0, m_vt1, 0, m_vt2, 0);
#pragma omp barrier

  multadd_Field_Gnd(Sigma_mu, 0, m_vt3, 0, V_nu, 0, 1.0);
#pragma omp barrier
}

//====================================================================
template<typename AFIELD>
void AForceSmear_APE<AFIELD>::force_each(AFIELD& Sigma_mu,
                                         const AFIELD& V_mu,
                                         const AFIELD& V_nu,
                                         const AFIELD& iTheta_mu,
                                         const AFIELD& iTheta_nu,
                                         const int mu, const int nu)
{
#pragma omp barrier

  Sigma_mu.set(0.0);
#pragma omp barrier

  m_shift_alt->backward(m_vt1_alt, V_nu, mu);
  m_shift_alt->backward(m_vt2_alt, V_mu, nu);

  Alt_Gauge::mult_Gnd(m_vt3_alt, 0, m_vt1_alt, 0, m_vt2_alt, 0);
#pragma omp barrier

  Alt_Gauge::multadd_Gnd(Sigma_mu, 0, m_vt3_alt, 0, iTheta_nu, 0, real_t(1.0));
#pragma omp barrier

  Alt_Gauge::mult_Gdn(m_vt3_alt, 0, iTheta_mu, 0, V_nu, 0);
#pragma omp barrier

  Alt_Gauge::mult_Gdn(m_vt2_alt, 0, m_vt1_alt, 0, m_vt3_alt, 0);
#pragma omp barrier

  m_shift_alt->forward(m_vt3_alt, m_vt2_alt, nu);
  axpy(Sigma_mu, real_t(1.0), m_vt3_alt);
#pragma omp barrier

  Alt_Gauge::mult_Gdn(m_vt3_alt, 0, V_mu, 0, iTheta_nu, 0);
#pragma omp barrier

  Alt_Gauge::mult_Gdn(m_vt2_alt, 0, m_vt1_alt, 0, m_vt3_alt, 0);
#pragma omp barrier

  m_shift_alt->forward(m_vt3_alt, m_vt2_alt, nu);
  axpy(Sigma_mu, real_t(1.0), m_vt3_alt);
#pragma omp barrier

  m_shift_alt->backward(m_vt1_alt, iTheta_nu, mu);
  m_shift_alt->backward(m_vt2_alt, V_mu, nu);

  Alt_Gauge::mult_Gnd(m_vt3_alt, 0, m_vt1_alt, 0, m_vt2_alt, 0);
#pragma omp barrier

  Alt_Gauge::multadd_Gnd(Sigma_mu, 0, m_vt3_alt, 0, V_nu, 0, real_t(1.0));
#pragma omp barrier

  Alt_Gauge::mult_Gdd(m_vt2_alt, 0, m_vt1_alt, 0, V_mu, 0);
#pragma omp barrier

  Alt_Gauge::mult_Gnn(m_vt3_alt, 0, m_vt2_alt, 0, V_nu, 0);
#pragma omp barrier

  m_shift_alt->forward(m_vt2_alt, m_vt3_alt, nu);

  axpy(Sigma_mu, 1.0, m_vt2_alt);
#pragma omp barrier

  m_shift_alt->backward(m_vt1_alt, V_nu, mu);
  m_shift_alt->backward(m_vt2_alt, iTheta_mu, nu);

  Alt_Gauge::mult_Gnd(m_vt3_alt, 0, m_vt1_alt, 0, m_vt2_alt, 0);
#pragma omp barrier

  Alt_Gauge::multadd_Gnd(Sigma_mu, 0, m_vt3_alt, 0, V_nu, 0, real_t(1.0));
#pragma omp barrier
}


//====================================================================
template<typename AFIELD>
void AForceSmear_APE<AFIELD>::staple(Field_G& c,
                            const Field_G& u_mu, const Field_G& u_nu,
                            const int mu, const int nu)
{
#pragma omp barrier

  //- upper direction
  m_shift->backward(m_vt1, u_mu, nu);

  mult_Field_Gnn(m_vt2, 0, u_nu, 0, m_vt1, 0);
#pragma omp barrier

  m_shift->backward(m_vt1, u_nu, mu);

  mult_Field_Gnd(c, 0, m_vt2, 0, m_vt1, 0);
#pragma omp barrier

  //- lower direction
  m_shift->backward(m_vt2, u_nu, mu);

  mult_Field_Gnn(m_vt1, 0, u_mu, 0, m_vt2, 0);
#pragma omp barrier

  mult_Field_Gdn(m_vt2, 0, u_nu, 0, m_vt1, 0);
#pragma omp barrier

  m_shift->forward(m_vt1, m_vt2, nu);

  axpy(c, 0, real_t(1.0), m_vt1, 0);
#pragma omp barrier
}


//====================================================================
template<typename AFIELD>
void AForceSmear_APE<AFIELD>::staple(AFIELD& c,
                            const AFIELD& u_mu, const AFIELD& u_nu,
                            const int mu, const int nu)
{
#pragma omp barrier

  //- upper direction
  m_shift_alt->backward(m_vt1_alt, u_mu, nu);

  Alt_Gauge::mult_Gnn(m_vt2_alt, 0, u_nu, 0, m_vt1_alt, 0);
#pragma omp barrier

  m_shift_alt->backward(m_vt1_alt, u_nu, mu);

  Alt_Gauge::mult_Gnd(c, 0, m_vt2_alt, 0, m_vt1_alt, 0);
#pragma omp barrier

  //- lower direction
  m_shift_alt->backward(m_vt2_alt, u_nu, mu);

  Alt_Gauge::mult_Gnn(m_vt1_alt, 0, u_mu, 0, m_vt2_alt, 0);
#pragma omp barrier

  Alt_Gauge::mult_Gdn(m_vt2_alt, 0, u_nu, 0, m_vt1_alt, 0);
#pragma omp barrier

  m_shift_alt->forward(m_vt1_alt, m_vt2_alt, nu);

  axpy(c, 0, real_t(1.0), m_vt1_alt, 0);
#pragma omp barrier
}


//============================================================END=====
