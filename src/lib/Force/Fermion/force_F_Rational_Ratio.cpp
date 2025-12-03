/*!
        @file    force_F_Rational_Ratio.cpp

        @brief

        @author  Hideo Matsufuru (matsufuru)
                 $LastChangedBy: kanamori $

        @date    $LastChangedDate:: 2025-12-03 19:35:35 #$

        @version $LastChangedRevision: 2668 $
*/

#include "force_F_Rational_Ratio.h"
#include "lib/Solver/shiftsolver_CG.h"
#include "lib/Tools/math_Rational.h"
#include "lib/Tools/timer.h"

const std::string Force_F_Rational_Ratio::class_name = "Force_F_Rational_Ratio";

//====================================================================
void Force_F_Rational_Ratio::set_parameters(const Parameters& params)
{
  std::string vlevel;
  if (!params.fetch_string("verbose_level", vlevel)) {
    m_vl = vout.set_verbose_level(vlevel);
  }

  //- fetch and check input parameters
  int err = 0;
  Parameters params_denom;
  Parameters params_numer;
  std::string key_denom("Fopr_Rational_denominator_MD");
  std::string key_numer("Fopr_Rational_numerator_MD");
  if (!params.find_Parameters(key_denom)) {
    vout.crucial(m_vl, "Error at %s: input parameter not found (key = %s).\n", class_name.c_str(), key_denom.c_str());
    exit(EXIT_FAILURE);
  }
  if (!params.find_Parameters(key_numer)) {
    vout.crucial(m_vl, "Error at %s: input parameter not found (key = %s).\n", class_name.c_str(), key_numer.c_str());
    exit(EXIT_FAILURE);
  }
  params_numer = params.lookup(key_numer);
  params_denom = params.lookup(key_denom);

  // parameters for numerator: rational approximation of (A^dag A)^{1/4}
  int    numer_Np, numer_n_exp, numer_d_exp;
  double numer_x_min, numer_x_max;
  int    numer_Niter;
  double numer_Stop_cond;
  std::string str_numer_vlevel;

  err += params_numer.fetch_int("number_of_poles", numer_Np);
  err += params_numer.fetch_int("exponent_numerator", numer_n_exp);
  err += params_numer.fetch_int("exponent_denominator", numer_d_exp);
  err += params_numer.fetch_double("lower_bound", numer_x_min);
  err += params_numer.fetch_double("upper_bound", numer_x_max);
  err += params_numer.fetch_int("maximum_number_of_iteration", numer_Niter);
  err += params_numer.fetch_double("convergence_criterion_squared", numer_Stop_cond);
  err += params_numer.fetch_string("verbose_level", str_numer_vlevel);

  if (err) {
    vout.crucial(m_vl, "Error at %s: input parameter not found (in Fopr_Rational_numerminator_MD).\n", class_name.c_str());
    exit(EXIT_FAILURE);
  }

  // parameters for denominator: rational approximation of (B^dag B)^{-1/2}
  int    denom_Np, denom_n_exp, denom_d_exp;
  double denom_x_min, denom_x_max;
  int    denom_Niter;
  double denom_Stop_cond;
  std::string str_denom_vlevel;

  err += params_denom.fetch_int("number_of_poles", denom_Np);
  err += params_denom.fetch_int("exponent_numerator", denom_n_exp);
  err += params_denom.fetch_int("exponent_denominator", denom_d_exp);
  err += params_denom.fetch_double("lower_bound", denom_x_min);
  err += params_denom.fetch_double("upper_bound", denom_x_max);
  err += params_denom.fetch_int("maximum_number_of_iteration", denom_Niter);
  err += params_denom.fetch_double("convergence_criterion_squared", denom_Stop_cond);
  err += params_denom.fetch_string("verbose_level", str_denom_vlevel);
  if (err) {
    vout.crucial(m_vl, "Error at %s: input parameter not found (in Fopr_Rational_denominator_MD).\n", class_name.c_str());
    exit(EXIT_FAILURE);
  }


  set_parameters(numer_Np, numer_n_exp, numer_d_exp, numer_x_min, numer_x_max,
                 numer_Niter, numer_Stop_cond, str_numer_vlevel,
                 denom_Np, denom_n_exp, denom_d_exp, denom_x_min, denom_x_max,
                 denom_Niter, denom_Stop_cond, str_denom_vlevel);

}


//====================================================================
void Force_F_Rational_Ratio::get_parameters(Parameters& params) const
{

  params.set_string("verbose_level", vout.get_verbose_level(m_vl));

  Parameters params_numer;
  params_numer.set_int("number_of_poles", m_numer_Np);
  params_numer.set_int("exponent_numerator", m_numer_n_exp);
  params_numer.set_int("exponent_denominator", m_numer_d_exp);
  params_numer.set_double("lower_bound", m_numer_x_min);
  params_numer.set_double("upper_bound", m_numer_x_max);
  params_numer.set_int("maximum_number_of_iteration", m_numer_Niter);
  params_numer.set_double("convergence_criterion_squared", m_numer_Stop_cond);
  params_numer.set_string("verobse_level", m_numer_str_vlevel);
  params.set_Parameters("Fopr_Rational_numerator_MD", params_numer);

  Parameters params_denom;
  params_denom.set_int("number_of_poles", m_denom_Np);
  params_denom.set_int("exponent_numerator", m_denom_n_exp);
  params_denom.set_int("exponent_denominator", m_denom_d_exp);
  params_denom.set_double("lower_bound", m_denom_x_min);
  params_denom.set_double("upper_bound", m_denom_x_max);
  params_denom.set_int("maximum_number_of_iteration", m_denom_Niter);
  params_denom.set_double("convergence_criterion_squared", m_denom_Stop_cond);
  params_denom.set_string("verobse_level", m_denom_str_vlevel);
  params.set_Parameters("Fopr_Rational_denomminator_MD", params_denom);

}


//====================================================================
void Force_F_Rational_Ratio::set_parameters(
           const int numer_Np, const int numer_n_exp, const int numer_d_exp,
           const double numer_x_min, const double numer_x_max,
           const int numer_Niter, const double numer_Stop_cond,
           const std::string numer_str_vlevel,
           const int denom_Np, const int denom_n_exp, const int denom_d_exp,
           const double denom_x_min, const double denom_x_max,
           const int denom_Niter, const double denom_Stop_cond,
           const std::string denom_str_vlevel)
{
  //- print input parameters
  vout.general(m_vl, "%s:\n", class_name.c_str());
  vout.general(m_vl, "  numerator operator\n");
  vout.general(m_vl, "    Np        = %d\n", numer_Np);
  vout.general(m_vl, "    n_exp     = %d\n", numer_n_exp);
  vout.general(m_vl, "    d_exp     = %d\n", numer_d_exp);
  vout.general(m_vl, "    x_min     = %12.8f\n", numer_x_min);
  vout.general(m_vl, "    x_max     = %12.8f\n", numer_x_max);
  vout.general(m_vl, "    Niter     = %d\n", numer_Niter);
  vout.general(m_vl, "    Stop_cond = %8.2e\n", numer_Stop_cond);
  vout.general(m_vl, "  denominator operator\n");
  vout.general(m_vl, "    Np        = %d\n", denom_Np);
  vout.general(m_vl, "    n_exp     = %d\n", denom_n_exp);
  vout.general(m_vl, "    d_exp     = %d\n", denom_d_exp);
  vout.general(m_vl, "    x_min     = %12.8f\n", denom_x_min);
  vout.general(m_vl, "    x_max     = %12.8f\n", denom_x_max);
  vout.general(m_vl, "    Niter     = %d\n", denom_Niter);
  vout.general(m_vl, "    Stop_cond = %8.2e\n", denom_Stop_cond);

  //- range check
  int err = 0;
  err += ParameterCheck::non_zero(numer_Np);
  err += ParameterCheck::non_zero(numer_n_exp);
  err += ParameterCheck::non_zero(numer_d_exp);
  // NB. x_min,x_max=0 is allowed.
  err += ParameterCheck::non_zero(numer_Niter);
  err += ParameterCheck::square_non_zero(numer_Stop_cond);
  if (err) {
    vout.crucial(m_vl, "Error at %s: parameter range check failed (numerator).\n", class_name.c_str());
    exit(EXIT_FAILURE);
  }

  err += ParameterCheck::non_zero(denom_Np);
  err += ParameterCheck::non_zero(denom_n_exp);
  err += ParameterCheck::non_zero(denom_d_exp);
  // NB. x_min,x_max=0 is allowed.
  err += ParameterCheck::non_zero(denom_Niter);
  err += ParameterCheck::square_non_zero(denom_Stop_cond);
  if (err) {
    vout.crucial(m_vl, "Error at %s: parameter range check failed (denominator).\n", class_name.c_str());
    exit(EXIT_FAILURE);
  }


  //- store values
  m_numer_Np        = numer_Np;
  m_numer_n_exp     = numer_n_exp;
  m_numer_d_exp     = numer_d_exp;
  m_numer_x_min     = numer_x_min;
  m_numer_x_max     = numer_x_max;
  m_numer_Niter     = numer_Niter;
  m_numer_Stop_cond = numer_Stop_cond;
  m_numer_str_vlevel = numer_str_vlevel;

  m_denom_Np        = denom_Np;
  m_denom_n_exp     = denom_n_exp;
  m_denom_d_exp     = denom_d_exp;
  m_denom_x_min     = denom_x_min;
  m_denom_x_max     = denom_x_max;
  m_denom_Niter     = denom_Niter;
  m_denom_Stop_cond = denom_Stop_cond;
  m_denom_str_vlevel = denom_str_vlevel;

  //- Rational approximation
  vout.general(m_vl, " Rational parameters for the numerator\n");
  set_rational_parameters(m_numer_a0, m_numer_bl, m_numer_cl,
                          m_numer_Np, m_numer_n_exp, m_numer_d_exp,
                          m_numer_x_min, m_numer_x_max);

  vout.general(m_vl, " Rational parameters for the denominator\n");
  set_rational_parameters(m_denom_a0, m_denom_bl, m_denom_cl,
                          m_denom_Np, m_denom_n_exp, m_denom_d_exp,
                          m_denom_x_min, m_denom_x_max);

  // create shiftsolvers
  Parameters params_solver_numer;
  params_solver_numer.set_int("maximum_number_of_iteration", m_numer_Niter);
  params_solver_numer.set_double("convergence_criterion_squared", m_numer_Stop_cond);
  params_solver_numer.set_string("verbose_level", m_numer_str_vlevel);
  m_shiftsolver_numer.reset(new Shiftsolver_CG(m_fopr_numer));
  m_shiftsolver_numer->set_parameters(params_solver_numer);

  Parameters params_solver_denom;
  params_solver_denom.set_int("maximum_number_of_iteration", m_denom_Niter);
  params_solver_denom.set_double("convergence_criterion_squared", m_denom_Stop_cond);
  params_solver_denom.set_string("verbose_level", m_denom_str_vlevel);
  m_shiftsolver_denom.reset(new Shiftsolver_CG(m_fopr_denom));
  m_shiftsolver_denom->set_parameters(params_solver_denom);
}

//====================================================================
void Force_F_Rational_Ratio::set_rational_parameters(
          double &a0, std::vector<double> &bl, std::vector<double> &cl,
          const int Np, const int n_exp, const int d_exp,
          const double x_min, const double x_max)
{
    const double x_min2 = x_min * x_min;
    const double x_max2 = x_max * x_max;

    bl.resize(Np);
    cl.resize(Np);
    Math_Rational rational;
    rational.set_parameters(Np, n_exp, d_exp, x_min2, x_max2);
    rational.get_parameters(a0, bl, cl);

    vout.general(m_vl, "    a0 = %18.14f\n", a0);
    for (int i = 0; i < Np; i++) {
      vout.general(m_vl, "    bl[%d] = %18.14f  cl[%d] = %18.14f\n",
                   i, bl[i], i, cl[i]);
    }

}

//====================================================================
void Force_F_Rational_Ratio::force_udiv(Field& force_, const Field& phi_)
{
  const int Nvol = CommonParameters::Nvol();
  const int Ndim = CommonParameters::Ndim();

  Field_G force(Nvol, Ndim);
  //  Field_F eta(eta_);

  force_udiv_impl(force, phi_);

  copy(force_, force); // force_ = force;
}


//====================================================================
void Force_F_Rational_Ratio::force_udiv_impl(Field_G& force, const Field_F& phi)
{
  const int Nvol = CommonParameters::Nvol();
  const int Ndim = CommonParameters::Ndim();

  const int NinF  = phi.nin();
  const int NvolF = phi.nvol();
  const int NexF  = phi.nex();

  //- Shiftsolver
  const int Nshift_denom = m_denom_Np;
  const int Nshift_numer = m_numer_Np;

  std::vector<Field> eta1(Nshift_numer); // eta
  std::vector<Field> zeta(Nshift_numer); // zeta
  std::vector<Field> eta2(Nshift_denom); // eta_tilde

  for (int i = 0; i < Nshift_numer; ++i) {
    eta1[i].reset(NinF, NvolF, NexF);
  }
  for (int i = 0; i < Nshift_numer; ++i) {
    zeta[i].reset(NinF, NvolF, NexF);
  }
  for (int i = 0; i < Nshift_denom; ++i) {
    eta2[i].reset(NinF, NvolF, NexF);
  }

  // for shiftsolvers
  m_fopr_denom->set_mode("DdagD");
  m_fopr_numer->set_mode("DdagD");

  Field phi2(NinF, NvolF, NexF); // (A^dag A)^{1/4} phi
  Field xi(NinF, NvolF, NexF);   // (B^dag B)^{-1/2} phi2

  Field_G force1(Nvol, Ndim);
  Field_F H_eta1(NvolF, NexF);
  Field_F H_zeta(NvolF, NexF);

#pragma omp parallel
  {

    {
      vout.general(m_vl, "    Shift solver in force calculation (for numerator, 1st)\n");
      vout.general(m_vl, "      Number of shift values = %d\n", m_numer_cl.size());

      int            Nconv;
      double         diff;
      Timer timer("shiftsolver, 1st numerator");
      timer.start();
      m_shiftsolver_numer->solve(eta1, m_numer_cl, phi, Nconv, diff);
      timer.stop();
      vout.general(m_vl, "      Nconv     = %d  \n", Nconv);
      vout.general(m_vl, "      diff(max) = %22.15e  \n", diff);
      timer.report();

      copy(phi2, phi);
      scal(phi2, m_numer_a0);

      for (int i = 0; i < Nshift_numer; ++i) {
        axpy(phi2, m_numer_bl[i], eta1[i]);
      }
    } // phi2 and eta1 are ready
#pragma omp barrier

    {
      vout.general(m_vl, "    Shift solver in force calculation (for denominator)\n");
      vout.general(m_vl, "      Number of shift values = %d\n", m_denom_cl.size());

      int            Nconv;
      double         diff;
      Timer timer("shiftsolver, denominator");
      timer.start();
      m_shiftsolver_denom->solve(eta2, m_denom_cl, phi2, Nconv, diff);
      timer.stop();
      vout.general(m_vl, "      Nconv     = %d  \n", Nconv);
      vout.general(m_vl, "      diff(max) = %22.15e  \n", diff);
      timer.report();

      copy(xi, phi2);
      scal(xi, m_denom_a0);
      for (int i = 0; i < Nshift_denom; ++i) {
        axpy(xi, m_denom_bl[i], eta2[i]);
      }
    } // xi and eta2 are ready

    {
      vout.general(m_vl, "    Shift solver in force calculation (for numerator, 2nd)\n");
      vout.general(m_vl, "      Number of shift values = %d\n", m_numer_cl.size());

      int            Nconv;
      double         diff;
      Timer timer("shiftsolver, 2nd numerator");
      timer.start();
      m_shiftsolver_numer->solve(zeta, m_numer_cl, xi, Nconv, diff);
      timer.stop();
      vout.general(m_vl, "      Nconv     = %d  \n", Nconv);
      vout.general(m_vl, "      diff(max) = %22.15e  \n", diff);
      timer.report();
    } // zeta is ready

    force.set(0.0);

  } // omp parallel

  // N.B. force_udiv may not be thread safe
  for (int i = 0; i < Nshift_denom; ++i) {
    m_force_denom->force_udiv(force1, eta2[i]); // symmetric
    axpy(force, m_denom_bl[i], force1);
  }

  for (int i = 0; i < Nshift_numer; ++i) {
    m_fopr_numer->set_mode("H");
#pragma omp parallel
    {
      m_fopr_numer->mult(H_eta1, eta1[i]);
      m_fopr_numer->mult(H_zeta, zeta[i]);
    }
    double b=m_numer_bl[i];

    m_force_numer->set_mode("H");
    m_force_numer->force_udiv1(force1, H_eta1, zeta[i]); // <eta1_i| H^dag H_div |zeta_i>
    axpy(force, b, force1);
    m_force_numer->force_udiv1(force1, H_zeta, eta1[i]); // <zeta_i| H^dag H_div |eta1_i>
    axpy(force, b, force1);

    m_force_numer->set_mode("Hdag");
    m_force_numer->force_udiv1(force1, eta1[i], H_zeta); // <eta1_i| Hdag_div H |zeta_i>
    axpy(force, b, force1);
    m_force_numer->force_udiv1(force1, zeta[i], H_eta1); // <zeta_i| Hdag_div H |eta1_i>
    axpy(force, b, force1);

  }
}


//====================================================================
void Force_F_Rational_Ratio::force_core1(Field&, const Field&, const Field&)
{
  vout.crucial(m_vl, "Error at %s: not implemented.\n", __func__);
  exit(EXIT_FAILURE);
}


//====================================================================
void Force_F_Rational_Ratio::force_udiv1(Field&, const Field&, const Field&)
{
  vout.crucial(m_vl, "Error at %s: not implemented.\n", __func__);
  exit(EXIT_FAILURE);
}


//====================================================================
//===========================================================END======
