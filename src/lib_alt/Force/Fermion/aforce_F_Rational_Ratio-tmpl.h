/*!
        @file    aforce_F_Rational_Ratio-tmpl.h

        @brief

        @author  Issaku Kanamori (kanamori)
                 $LastChangedBy: kanamori $

        @date    $LastChangedDate:: 2025-12-03 19:35:35 #$

        @version $LastChangedRevision: 2668 $
*/

#include "aforce_F_Rational_Ratio.h"
#include "lib/Tools/timer.h"

template<typename AFIELD>
const std::string AForce_F_Rational_Ratio<AFIELD>::class_name
                                             = "AForce_F_Rational_Ratio";
//====================================================================
template<typename AFIELD>
void AForce_F_Rational_Ratio<AFIELD>::init(const Parameters& params)
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

  set_parameters(params);

  int NinF  = m_fopr1->field_nin();
  int NvolF = m_fopr1->field_nvol();
  int NexF  = m_fopr1->field_nex();
  assert(m_fopr2->field_nin() == NinF);
  assert(m_fopr2->field_nvol() == NvolF);
  assert(m_fopr2->field_nex() == NexF);

  m_v1.reset(NinF, NvolF, NexF);
  m_v2.reset(NinF, NvolF, NexF);



  int Nc   = CommonParameters::Nc();
  int NinG = 2 * Nc * Nc;
  int Nvol = CommonParameters::Nvol();
  int Ndim = CommonParameters::Ndim();

  m_force.reset(NinG, Nvol, Ndim);

  vout.decrease_indent();
  vout.general(m_vl, "%s: construction finished.\n",
               class_name.c_str());

}

//====================================================================
template<typename AFIELD>
void AForce_F_Rational_Ratio<AFIELD>::tidyup()
{
  //  delete m_solver;

}

//====================================================================
template<typename AFIELD>
void AForce_F_Rational_Ratio<AFIELD>::set_parameters(const Parameters& params)
{

  std::string vlevel;
  if (!params.fetch_string("verbose_level", vlevel)) {
    m_vl = vout.set_verbose_level(vlevel);
  } else {
    m_vl = CommonParameters::Vlevel();
  }

  int err=0;
  std::string key1("Fopr_Rational_numerator_MD");
  std::string key2("Fopr_Rational_denominator_MD");
  if (!params.find_Parameters(key1)) {
    vout.crucial(m_vl, "Error at %s: input parameter not found (key = %s).\n", class_name.c_str(), key1.c_str());
    exit(EXIT_FAILURE);
  }
  if (!params.find_Parameters(key2)) {
    vout.crucial(m_vl, "Error at %s: input parameter not found (key = %s).\n", class_name.c_str(), key2.c_str());
    exit(EXIT_FAILURE);
  }
  Parameters params1 = params.lookup(key1);
  Parameters params2 = params.lookup(key2);

  // parameters for numerator: rational approximation of (A^dag A)^{1/4}
  int    Np1, n_exp1, d_exp1;
  double x_min1, x_max1;
  int    Niter1;
  double Stop_cond1;
  std::string str_vlevel1;

  // for rational parameters
  err += params1.fetch_int(   "number_of_poles",      Np1);
  err += params1.fetch_int(   "exponent_numerator",   n_exp1);
  err += params1.fetch_int(   "exponent_denominator", d_exp1);
  err += params1.fetch_double("lower_bound",          x_min1);
  err += params1.fetch_double("upper_bound",          x_max1);
  // for shift solver
  err += params1.fetch_int(   "maximum_number_of_iteration",   Niter1);
  err += params1.fetch_double("convergence_criterion_squared", Stop_cond1);
  err += params1.fetch_string("verbose_level",                 str_vlevel1);

  if (err) {
    vout.crucial(m_vl, "Error at %s: input parameter not found (in %s).\n",
                 class_name.c_str(), key1.c_str());
    exit(EXIT_FAILURE);
  }

  // parameters for denominator: rational approximation of (B^dag B)^{-1/2}
  int    Np2, n_exp2, d_exp2;
  double x_min2, x_max2;
  int    Niter2;
  double Stop_cond2;
  std::string str_vlevel2;

  // for rational parameters
  err += params2.fetch_int(   "number_of_poles",      Np2);
  err += params2.fetch_int(   "exponent_numerator",   n_exp2);
  err += params2.fetch_int(   "exponent_denominator", d_exp2);
  err += params2.fetch_double("lower_bound",          x_min2);
  err += params2.fetch_double("upper_bound",          x_max2);
  // for shift solver
  err += params2.fetch_int(   "maximum_number_of_iteration",   Niter2);
  err += params2.fetch_double("convergence_criterion_squared", Stop_cond2);
  err += params2.fetch_string("verbose_level",                 str_vlevel2);
  if (err) {
    vout.crucial(m_vl, "Error at %s: input parameter not found (in %s).\n",
                 class_name.c_str(), key2.c_str());
    exit(EXIT_FAILURE);
  }


  set_parameters(Np1, n_exp1, d_exp1, x_min1, x_max1,
                 Niter1, Stop_cond1, str_vlevel1,
                 Np2, n_exp2, d_exp2, x_min2, x_max2,
                 Niter2, Stop_cond2, str_vlevel2);

}


//====================================================================
template<typename AFIELD>
void AForce_F_Rational_Ratio<AFIELD>::set_parameters(
                                  const int Np1, const int n_exp1, const int d_exp1,
                                  const real_t x_min1, const real_t x_max1,
                                  const int Niter1, const real_t Stop_cond1,
                                  const std::string str_vlevel1,
                                  const int Np2, const int n_exp2, const int d_exp2,
                                  const real_t x_min2, const real_t x_max2,
                                  const int Niter2, const real_t Stop_cond2,
                                  const std::string str_vlevel2)
{


  //- range check
  int err = 0;
  err += ParameterCheck::non_zero(Np1);
  err += ParameterCheck::non_zero(n_exp1);
  err += ParameterCheck::non_zero(d_exp1);
  // NB. x_min,x_max=0 is allowed.
  err += ParameterCheck::non_zero(Niter1);
  err += ParameterCheck::square_non_zero(Stop_cond1);
  if (err) {
    vout.crucial(m_vl, "Error at %s: parameter range check failed (numerator).\n", class_name.c_str());
    exit(EXIT_FAILURE);
  }

  err += ParameterCheck::non_zero(Np2);
  err += ParameterCheck::non_zero(n_exp2);
  err += ParameterCheck::non_zero(d_exp2);
  // NB. x_min,x_max=0 is allowed.
  err += ParameterCheck::non_zero(Niter2);
  err += ParameterCheck::square_non_zero(Stop_cond2);
  if (err) {
    vout.crucial(m_vl, "Error at %s: parameter range check failed (denominator).\n", class_name.c_str());
    exit(EXIT_FAILURE);
  }


  int ith = ThreadManager::get_thread_id();
  if (ith == 0) {
    m_Np1        = Np1;
    m_n_exp1     = n_exp1;
    m_d_exp1     = d_exp1;
    m_x_min1     = x_min1;
    m_x_max1     = x_max1;
    m_Niter1     = Niter1;
    m_Stop_cond1 = Stop_cond1;
    m_str_vlevel1 = str_vlevel1;

    m_Np2        = Np2;
    m_n_exp2     = n_exp2;
    m_d_exp2     = d_exp2;
    m_x_min2     = x_min2;
    m_x_max2     = x_max2;
    m_Niter2     = Niter2;
    m_Stop_cond2 = Stop_cond2;
    m_str_vlevel2 = str_vlevel2;

    //- Rational approximation
    vout.general(m_vl, " Rational parameters for the numerator (fopr1)\n");
    set_rational_parameters(m_a01, m_bl1, m_cl1,
                            m_Np1, m_n_exp1, m_d_exp1,
                            m_x_min1, m_x_max1);

    vout.general(m_vl, " Rational parameters for the denominator (fopr2)\n");
    set_rational_parameters(m_a02, m_bl2, m_cl2,
                            m_Np2, m_n_exp2, m_d_exp2,
                            m_x_min2, m_x_max2);


    // create shiftsolvers
    Parameters params_solver1; // numerator: (A^dag A)^{1/4}
    params_solver1.set_int("maximum_number_of_iteration", m_Niter1);
    params_solver1.set_double("convergence_criterion_squared", m_Stop_cond1);
    params_solver1.set_string("verbose_level", m_str_vlevel1);
    m_solver1.reset(new AShiftsolver_CG<AFIELD, AFopr<AFIELD> >(m_fopr1));
    m_solver1->set_parameters(params_solver1);

    Parameters params_solver2; // denominator: (B^dag B)^{-1/2}
    params_solver2.set_int("maximum_number_of_iteration", m_Niter2);
    params_solver2.set_double("convergence_criterion_squared", m_Stop_cond2);
    params_solver2.set_string("verbose_level", m_str_vlevel2);
    m_solver2.reset(new AShiftsolver_CG<AFIELD, AFopr<AFIELD> >(m_fopr2));
    m_solver2->set_parameters(params_solver2);


    m_eta1.resize(m_Np1);
    m_eta2.resize(m_Np2);
    m_zeta1.resize(m_Np1);
    const int NinF  = m_fopr1->field_nin();
    const int NvolF = m_fopr1->field_nvol();
    const int NexF  = m_fopr1->field_nex();
    for (int i = 0; i < m_Np1; ++i) {
      m_eta1[i].reset(NinF, NvolF, NexF);
    }
    for (int i = 0; i < m_Np1; ++i) {
      m_zeta1[i].reset(NinF, NvolF, NexF);
    }
    for (int i = 0; i < m_Np2; ++i) {
      m_eta2[i].reset(NinF, NvolF, NexF);
    }

  } // ith==0
}


//====================================================================
template<typename AFIELD>
void AForce_F_Rational_Ratio<AFIELD>::set_rational_parameters(
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
template<typename AFIELD>
void AForce_F_Rational_Ratio<AFIELD>::get_parameters(Parameters& params) const
{
  Parameters params1;
  params1.set_int("number_of_poles", m_Np1);
  params1.set_int("exponent_numerator", m_n_exp1);
  params1.set_int("exponent_denominator", m_d_exp1);
  params1.set_double("lower_bound", double(m_x_min1));
  params1.set_double("upper_bound", double(m_x_max1));
  params1.set_int("maximum_number_of_iteration", m_Niter1);
  params1.set_double("convergence_criterion_squared", double(m_Stop_cond1));
  params1.set_string("verobse_level", m_str_vlevel1);
  params.set_Parameters("Fopr_Rational_numerator_MD", params1);

  Parameters params2;
  params2.set_int("number_of_poles", m_Np2);
  params2.set_int("exponent_numerator", m_n_exp2);
  params2.set_int("exponent_denominator", m_d_exp2);
  params2.set_double("lower_bound", double(m_x_min2));
  params2.set_double("upper_bound", double(m_x_max2));
  params2.set_int("maximum_number_of_iteration", m_Niter2);
  params2.set_double("convergence_criterion_squared", double(m_Stop_cond2));
  params2.set_string("verobse_level", m_str_vlevel2);
  params.set_Parameters("Fopr_Rational_denomminator_MD", params2);

  params.set_string("verbose_level", vout.get_verbose_level(m_vl));
}

//====================================================================
template<typename AFIELD>
void AForce_F_Rational_Ratio<AFIELD>::set_config(Field *U)
{
  m_U = (Field_G *)U;
  m_fopr1->set_config(U);
  m_fopr2->set_config(U);
  m_force1->set_config(U);
  m_force2->set_config(U);
}

//====================================================================
template<typename AFIELD>
void AForce_F_Rational_Ratio<AFIELD>::force_udiv(AFIELD& force,
                                           const AFIELD& eta)
{
  force_udiv_impl(force, eta);
}


//====================================================================
template<typename AFIELD>
void AForce_F_Rational_Ratio<AFIELD>::force_udiv1(AFIELD& force,
                                            const AFIELD& zeta,
                                            const AFIELD& eta)
{
  force_udiv1_impl(force, zeta, eta);
}


//====================================================================
template<typename AFIELD>
void AForce_F_Rational_Ratio<AFIELD>::force_udiv_impl(AFIELD& force,
                                                const AFIELD& eta)
{
#pragma omp barrier

  int Nvol = CommonParameters::Nvol();
  int Ndim = CommonParameters::Ndim();

  int NinF  = eta.nin();
  int NvolF = eta.nvol();
  int NexF  = eta.nex();

  int Nc   = CommonParameters::Nc();
  int NinG = 2 * Nc * Nc;

  m_fopr1->set_mode("DdagD");
  m_fopr2->set_mode("DdagD");

  {
    vout.general(m_vl, "    Shift solver in force calculation (for numerator, 1st)\n");
    vout.general(m_vl, "      Number of shift values = %d\n", m_cl1.size());

    int            Nconv;
    double         diff;
    Timer timer("shiftsolver, 1st numerator");
    timer.start();
    m_solver1->solve(m_eta1, m_cl1, eta, Nconv, diff);
    timer.stop();
    vout.general(m_vl, "      Nconv     = %d  \n", Nconv);
    vout.general(m_vl, "      diff(max) = %22.15e  \n", diff);
    //timer.report();
    double elapsed_time = timer.elapsed_sec();
    vout.general(m_vl, "      elapsed time: %12.6f [sec]\n", elapsed_time);

    copy(m_v1, eta);
    scal(m_v1, m_a01);

    for (int i = 0; i < m_Np1; ++i) {
      axpy(m_v1, m_bl1[i], m_eta1[i]);
    }
  } // m_v1 = (A^dag A=^{1/4} eta and m_eta1 are ready
#pragma omp barrier

  {
    vout.general(m_vl, "    Shift solver in force calculation (for denominator)\n");
    vout.general(m_vl, "      Number of shift values = %d\n", m_cl2.size());

    int            Nconv;
    double         diff;
    Timer timer("shiftsolver, denominator");
    timer.start();
    m_solver2->solve(m_eta2, m_cl2, m_v1, Nconv, diff);
    timer.stop();
    vout.general(m_vl, "      Nconv     = %d  \n", Nconv);
    vout.general(m_vl, "      diff(max) = %22.15e  \n", diff);
    //    timer.report();
    double elapsed_time = timer.elapsed_sec();
    vout.general(m_vl, "      elapsed time: %12.6f [sec]\n", elapsed_time);

    copy(m_v2, m_v1);
    scal(m_v2, m_a02);
    for (int i = 0; i < m_Np2; ++i) {
      axpy(m_v2, m_bl2[i], m_eta2[i]);
    }
  } // m_v2 = (B^dag B)^{-1/2} m_v1  and m_eta2 are ready

#pragma omp barrier
  {
    vout.general(m_vl, "    Shift solver in force calculation (for numerator, 2nd)\n");
    vout.general(m_vl, "      Number of shift values = %d\n", m_cl1.size());

    int            Nconv;
    double         diff;
    Timer timer("shiftsolver, 2nd numerator");
    timer.start();
    m_solver1->solve(m_zeta1, m_cl1, m_v2, Nconv, diff);
    timer.stop();
    vout.general(m_vl, "      Nconv     = %d  \n", Nconv);
    vout.general(m_vl, "      diff(max) = %22.15e  \n", diff);
    //    timer.report();
    double elapsed_time = timer.elapsed_sec();
    vout.general(m_vl, "      elapsed time: %12.6f [sec]\n", elapsed_time);

  } // m_zeta is ready

  force.set(0.0);

  for (int i = 0; i < m_Np2; ++i) {
    m_force2->force_udiv(m_force, m_eta2[i]); // symmetric
    axpy(force, m_bl2[i], m_force);
  }

  m_fopr1->set_mode("H");

  for (int i = 0; i < m_Np1; ++i) {
    double b=m_bl1[i];
    m_fopr1->mult(m_v1, m_eta1[i]);
    m_fopr1->mult(m_v2, m_zeta1[i]);

    m_force1->set_mode("H");
    m_force1->force_udiv1(m_force, m_v1, m_zeta1[i]); // <eta1_i| H^dag H_div |zeta_i>
    axpy(force, b, m_force);
    m_force1->force_udiv1(m_force, m_v2, m_eta1[i]); // <zeta_i| H^dag H_div |eta1_i>
    axpy(force, b, m_force);

    m_force1->set_mode("Hdag");
    m_force1->force_udiv1(m_force, m_eta1[i], m_v2); // <eta1_i| Hdag_div H |zeta_i>
    axpy(force, b, m_force);
    m_force1->force_udiv1(m_force, m_zeta1[i], m_v1); // <zeta_i| Hdag_div H |eta1_i>
    axpy(force, b, m_force);

#pragma omp barrier
  }

}


//====================================================================
template<typename AFIELD>
void AForce_F_Rational_Ratio<AFIELD>::force_udiv1_impl(AFIELD& force,
                                                const AFIELD& zeta,
                                                const AFIELD& eta)
{
  vout.crucial(m_vl, "Error at %s: %s is not implemented.\n",
               class_name.c_str(), __func__);
  exit(EXIT_FAILURE);
}



//====================================================================
template<typename AFIELD>
void AForce_F_Rational_Ratio<AFIELD>::force_core1(
                             AFIELD&, const AFIELD&, const AFIELD&)
{
  vout.crucial(m_vl, "Error at %s: %s not implemented.\n",
               class_name.c_str(), __func__);
  exit(EXIT_FAILURE);
}


//===========================================================END======
