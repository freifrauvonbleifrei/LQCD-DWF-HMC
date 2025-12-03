/*!
      @file    fprop_alt_Standard_eo_Richardson.cpp
      @brief
      @author  Issaku Kanamori (kanamori)
               $LastChangedBy: kanamori $
      @date    $LastChangedDate:: 2025-12-03 19:35:35 #$
      @version $LastChangedRevision: 2668 $
*/

template<typename AFIELD_d, typename AFIELD_f>
const std::string Fprop_alt_Standard_eo_Richardson<AFIELD_d, AFIELD_f>::class_name
  = "Fprop_alt_Standard_eo_Richardson";
//====================================================================
template<typename AFIELD_d, typename AFIELD_f>
void Fprop_alt_Standard_eo_Richardson<AFIELD_d, AFIELD_f>::init(const Parameters& params_fopr,
                                                                const Parameters& params_solver)
{
  // this constructor assumes that the factories are available.
  vout.general(m_vl, "%s: smearing is not set\n", class_name.c_str());
  Director_alt_Smear<AFIELD_d> *dr_smear = nullptr;
  init(params_fopr, params_solver, dr_smear);
}


//====================================================================
template<typename AFIELD_d, typename AFIELD_f>
void Fprop_alt_Standard_eo_Richardson<AFIELD_d, AFIELD_f>::init(const Parameters& params_fopr,
                                                                const Parameters& params_solver,
                                                                Director_alt_Smear<AFIELD_d> *dr_smear)
{
  // this constructor assumes that the factories are available.
  vout.general(m_vl, "\n");
  vout.general(m_vl, "%s: being setup.\n", class_name.c_str());

  m_dr_smear = dr_smear;

  string fopr_type = params_fopr.get_string("fermion_type");
  if(fopr_type.substr(fopr_type.size()-3 ,3) != "_eo") {
    fopr_type += "_eo";
  }


  m_kernel_d = nullptr;
  m_kernel_f = nullptr;

  if(m_dr_smear){
    m_kernel_d = AFopr<AFIELD_d>::New(fopr_type, params_fopr);
    m_fopr_d   = new AFopr_Smeared_alt<AFIELD_d>(m_kernel_d, m_dr_smear);
    m_fopr_f   = AFopr<AFIELD_f>::New(fopr_type, params_fopr);
  } else {
    m_fopr_d   = AFopr<AFIELD_d>::New(fopr_type, params_fopr);
    m_fopr_f   = AFopr<AFIELD_f>::New(fopr_type, params_fopr);
  }

  // solver parameters for preconditioner solver
  Parameters params_solver2;
  int        max_iter    = params_solver.get_int("maximum_number_of_iteration");
  int        max_restart = params_solver.get_int("maximum_number_of_restart");
  params_solver2.set_int("maximum_number_of_iteration", max_iter * max_restart);
  params_solver2.set_int("maximum_number_of_restart", 1);
  params_solver2.set_double("convergence_criterion_squared", 1.0e-8);
  string str_vlevel = params_solver.get_string("verbose_level");
  params_solver2.set_string("verbose_level", str_vlevel);
  params_solver2.set_string("initial_guess_mode", "ZERO");

  string solver_type = params_solver.get_string("solver_type");
  m_solver_prec = ASolver<AFIELD_f>::New(solver_type, m_fopr_f);
  m_solver_prec->set_parameters(params_solver2);

  m_precond
    = new APrecond_Mixedprec<AFIELD_d, AFIELD_f>(m_solver_prec);

  m_solver = new
             ASolver_Richardson<AFIELD_d>(m_fopr_d, m_precond);
  m_solver->set_parameters(params_solver);

  m_use_DdagD = false;
  if (solver_type == "CGNR") {
    m_use_DdagD = true;
  }
  reset_performance();

  int nin   = m_fopr_d->field_nin();
  int nvol2 = m_fopr_d->field_nvol();
  int nex   = m_fopr_d->field_nex();
  int nvol  = 2 * nvol2;

  m_axq.reset(nin, nvol, nex);
  m_abq.reset(nin, nvol, nex);

  m_be.reset(nin, nvol2, nex);
  m_bo.reset(nin, nvol2, nex);
  m_xe.reset(nin, nvol2, nex);
  m_xo.reset(nin, nvol2, nex);
  //  m_y1.reset(nin, nvol2, nex);
  //  m_y2.reset(nin, nvol2, nex);
  //  m_y1f.reset(nin, nvol2, nex);
  //  m_y2f.reset(nin, nvol2, nex);

  vout.general(m_vl, "%s: setup finished.\n", class_name.c_str());

}

//====================================================================
template<typename AFIELD_d, typename AFIELD_f>
void Fprop_alt_Standard_eo_Richardson<AFIELD_d, AFIELD_f>::tidyup()
{
  delete m_solver;
  delete m_precond;
  delete m_solver_prec;
  delete m_fopr_f;
  delete m_fopr_d;
  if (m_kernel_d != 0) delete m_kernel_d;
  if (m_kernel_f != 0) delete m_kernel_f;
}


//====================================================================
template<typename AFIELD_d, typename AFIELD_f>
void Fprop_alt_Standard_eo_Richardson<AFIELD_d, AFIELD_f>::set_config(Field *U)
{
  m_fopr_d->set_config(U);

  if(m_dr_smear) {
    int Nsmear  = m_dr_smear->get_Nsmear();
    Field* Uptr = m_dr_smear->getptr_smearedConfig(Nsmear);
    m_fopr_f->set_config(Uptr);
  } else {
    m_fopr_f->set_config(U);
  }
}


//====================================================================
template<typename AFIELD_d, typename AFIELD_f>
void Fprop_alt_Standard_eo_Richardson<AFIELD_d, AFIELD_f>::invert(Field& xq, const Field& b,
                                                                  int& nconv, double& diff)
{
  vout.paranoiac(m_vl, "%s: invert is called.\n", class_name.c_str());
  vout.paranoiac(m_vl, "mode = %s.\n", m_mode.c_str());

  if (m_mode == "D") {
    invert_D(xq, b, nconv, diff);
  } else if (m_mode == "DdagD") {
    invert_DdagD(xq, b, nconv, diff);
  } else if (m_mode == "DdagD_even") {
    //    invert_DdagD_even(xq, b, nconv, diff);
    vout.crucial(m_vl, "%s: unsupported mode: %s (not yet ready for Field)\n",
                 class_name.c_str(), m_mode.c_str());
  } else {
    vout.crucial(m_vl, "%s: unsupported mode: %s\n",
                 class_name.c_str(), m_mode.c_str());
    exit(EXIT_FAILURE);
  }
}


//====================================================================
template<typename AFIELD_d, typename AFIELD_f>
void Fprop_alt_Standard_eo_Richardson<AFIELD_d, AFIELD_f>::invert(AFIELD_d& xq, const AFIELD_d& b,
                                                                  int& nconv, double& diff)
{
  vout.paranoiac(m_vl, "%s: invert is called.\n", class_name.c_str());
  vout.paranoiac(m_vl, "mode = %s.\n", m_mode.c_str());

  if (m_mode == "D") {
    invert_D(xq, b, nconv, diff);
  } else if (m_mode == "DdagD") {
    invert_DdagD(xq, b, nconv, diff);
  } else if (m_mode == "DdagD_even") {
    invert_DdagD_even(xq, b, nconv, diff);
  } else {
    vout.crucial(m_vl, "%s: unsupported mode: %s\n",
                 class_name.c_str(), m_mode.c_str());
    exit(EXIT_FAILURE);
  }
}


//====================================================================
template<typename AFIELD_d, typename AFIELD_f>
void Fprop_alt_Standard_eo_Richardson<AFIELD_d, AFIELD_f>::invert_D(
  Field& xq, const Field& b,
  int& nconv, double& diff)
{
  ThreadManager::assert_single_thread(class_name);
  vout.paranoiac(m_vl, "invert_D start.\n");
#pragma omp parallel
  {
    Timer timer;
    timer.reset();
    timer.start();

    AIndex_lex<real_t, AFIELD_d::IMPL> index_alt;

    if (m_fopr_d->needs_convert()) {
      vout.detailed(m_vl, "convert required.\n");
      m_fopr_d->convert(m_abq, b);
    } else {
      vout.detailed(m_vl, "convert not required.\n");
      convert(index_alt, m_abq, b);
    }

#pragma omp barrier
    timer.stop();
    invert_D(m_axq, m_abq, nconv, diff);
    timer.start();
#pragma omp barrier
    if (m_fopr_d->needs_convert()) {
      m_fopr_d->reverse(xq, m_axq);
    } else {
      reverse(index_alt, xq, m_axq);
    }
    timer.stop();
#pragma omp master
    {
      m_elapsed_time += timer.elapsed_sec();
    }
  } // omp parallel
}


//====================================================================
template<typename AFIELD_d, typename AFIELD_f>
void Fprop_alt_Standard_eo_Richardson<AFIELD_d, AFIELD_f>::invert_DdagD(
  Field& xq, const Field& b,
  int& nconv, double& diff)
{
    Timer timer;
    timer.reset();
    timer.start();

    AIndex_lex<real_t, AFIELD_d::IMPL> index_alt;
#pragma omp parallel
    {
      if (m_fopr_d->needs_convert()) {
        vout.detailed(m_vl, "convert required.\n");
        m_fopr_d->convert(m_abq, b);
      } else {
        vout.detailed(m_vl, "convert not required.\n");
        convert(index_alt, m_abq, b);
      }
    }
    timer.stop();
    invert_DdagD(m_axq, m_abq, nconv, diff);
    timer.start();
#pragma omp parallel
    {
      if (m_fopr_d->needs_convert()) {
        m_fopr_d->reverse(xq, m_axq);
      } else {
        reverse(index_alt, xq, m_axq);
      }
    }
    timer.stop();
    m_elapsed_time += timer.elapsed_sec();
}


//====================================================================
template<typename AFIELD_d, typename AFIELD_f>
void Fprop_alt_Standard_eo_Richardson<AFIELD_d, AFIELD_f>::invert_D(
  AFIELD_d& xq, const AFIELD_d& b,
  int& nconv, double& diff)
{
  int nth = ThreadManager::get_num_threads();
  vout.paranoiac("invert_D start: nth=%d.\n", nth);

  if (nth > 1 ) {
    invert_D_th(xq, b, nconv, diff);
  } else {
#pragma omp parallel
    {
    invert_D_th(xq, b, nconv, diff);
    }
  }
}

//====================================================================
template<typename AFIELD_d, typename AFIELD_f>
void Fprop_alt_Standard_eo_Richardson<AFIELD_d, AFIELD_f>::invert_D_th(
  AFIELD_d& xq, const AFIELD_d& b,
  int& nconv, double& diff)
{
  //  ThreadManager::assert_single_thread(class_name);
  vout.paranoiac("invert_D_th start.\n");

#pragma omp master
  {
    m_timer.reset();
    m_timer.start();
  }

  //AIndex_lex<real_t, AFIELD_d::IMPL> index_alt;
  AIndex_eo<real_t, AFIELD_d::IMPL>  index_eo;

  index_eo.split(m_be, m_bo, b);
#pragma omp barrier

  invert_De(m_xe, m_xo, m_be, m_bo, nconv, diff);

#pragma omp barrier
  index_eo.merge(xq, m_xe, m_xo);

#pragma omp master
  {
    m_timer.stop();
    m_elapsed_time += m_timer.elapsed_sec();
    m_flop_count   += m_solver->flop_count();
  }
}


//====================================================================
template<typename AFIELD_d, typename AFIELD_f>
void Fprop_alt_Standard_eo_Richardson<AFIELD_d, AFIELD_f>::invert_DdagD(
  AFIELD_d& xq, const AFIELD_d& b,
  int& nconv, double& diff)
{
  ThreadManager::assert_single_thread(class_name);
  vout.paranoiac(m_vl, "invert_DdagD start.\n");

  #pragma omp parallel
  {
#pragma omp master
    {
      m_timer.reset();
      m_timer.start();
    }

    AIndex_eo<real_t, AFIELD_f::IMPL> index_eo;
    index_eo.split(m_xe, m_xo, b);
#pragma omp barrier

    int    nconv1;
    double diff1;
    invert_De_dag(m_be, m_bo, m_xe, m_xo, nconv1, diff1);
#pragma omp barrier

    int nconv2;
    double diff2;
    invert_De(m_xe, m_xo, m_be, m_bo, nconv2, diff2);

#pragma omp barrier
    nconv1 += nconv2;
    diff1  += diff2;

    index_eo.merge(xq, m_xe, m_xo);

#pragma omp master
    {
      m_timer.stop();
      m_elapsed_time += m_timer.elapsed_sec();
      m_flop_count   += m_solver->flop_count();
      nconv = nconv1;
      diff  = diff1;
    }
  }
}


//====================================================================
template<typename AFIELD_d, typename AFIELD_f>
void Fprop_alt_Standard_eo_Richardson<AFIELD_d, AFIELD_f>::invert_De(
  AFIELD_d& xe, AFIELD_d& xo,
  AFIELD_d& be, AFIELD_d& bo, // be, bo are destroyed
  int& nconv, double& diff)
{
  // set even source vector.
#pragma omp barrier

  m_fopr_d->mult(xo, bo, "Doo_inv");
#pragma omp barrier

  m_fopr_d->mult(xe, xo, "Deo");
#pragma omp barrier

  aypx(real_t(-1.0), xe, be);  // xe = -Deo Doo_inv bo + be
#pragma omp barrier

  m_fopr_d->mult(be, xe, "Dee_inv");

  if (m_use_DdagD) {
    m_fopr_d->set_mode("Ddag");
    m_fopr_d->mult(xo, be);
    m_fopr_d->set_mode("DdagD");
    m_fopr_f->set_mode("DdagD");
    START_FAPP("sovler", 1, 0);
    m_solver->solve(xe, xo, nconv, diff);
    STOP_FAPP("sovler", 1, 0);
  } else {
    m_fopr_d->set_mode("D");
    m_fopr_f->set_mode("D");
    START_FAPP("sovler", 1, 0);
    m_solver->solve(xe, be, nconv, diff);
    STOP_FAPP("sovler", 1, 0);
  }
#pragma omp barrier

  m_fopr_d->mult(xo, xe, "Doe");
#pragma omp barrier

  axpy(bo, real_t(-1.0), xo); // bo := bo - Doe xe 
#pragma omp barrier

  m_fopr_d->mult(xo, bo, "Doo_inv");
#pragma omp barrier
}


//====================================================================
template<typename AFIELD_d, typename AFIELD_f>
void Fprop_alt_Standard_eo_Richardson<AFIELD_d, AFIELD_f>::invert_De_dag(
  AFIELD_d& xe, AFIELD_d& xo,
  AFIELD_d& be, AFIELD_d& bo, // be, bo are destroyed
  int& nconv, double& diff)
{
  // set even source vector.
#pragma omp barrier

  m_fopr_d->mult_dag(xo, bo, "Doo_inv");
#pragma omp barrier

  m_fopr_d->mult_dag(xe, xo, "Deo");
#pragma omp barrier

  axpy(be, real_t(-1.0), xe); // be = be - Deo Boo_inv bo
#pragma omp barrier

  if (m_use_DdagD) {
    m_fopr_d->set_mode("DdagD");
    m_fopr_f->set_mode("DdagD");
    m_solver->solve(xe, be, nconv, diff);
    m_fopr_d->set_mode("D");
    m_fopr_d->mult(xo, xe);
  } else {
    m_fopr_d->set_mode("Ddag");
    m_fopr_f->set_mode("Ddag");
    m_solver->solve(xo, be, nconv, diff);
  }
#pragma omp barrier
  m_fopr_d->mult_dag(xe, xo, "Dee_inv");
#pragma omp barrier

  m_fopr_d->mult_dag(xo, xe, "Doe");
#pragma omp barrier

  axpy(bo, real_t(-1.0), xo);  // bo := bo - Doe xe
#pragma omp barrier

  m_fopr_d->mult_dag(xo, bo, "Doo_inv"); // xo = Doo_inv (bo - Doe xe)
#pragma omp barrier
}




//====================================================================
template<typename AFIELD_d, typename AFIELD_f>
void Fprop_alt_Standard_eo_Richardson<AFIELD_d, AFIELD_f>::invert_DdagD_even(
  AFIELD_d& xq, const AFIELD_d& b,
  int& nconv, double& diff)
{
  ThreadManager::assert_single_thread(class_name);
  vout.paranoiac("invert_DdagD_even start.\n");

  m_timer.reset();
  m_timer.start();

  m_fopr_d->set_mode("DdagD");
  m_fopr_f->set_mode("DdagD");

#pragma omp parallel
  {
    m_solver->solve(xq, b, nconv, diff);
  }

  m_timer.stop();
  m_elapsed_time += m_timer.elapsed_sec();
  m_flop_count   += m_solver->flop_count();
}


//====================================================================
template<typename AFIELD_d, typename AFIELD_f>
double Fprop_alt_Standard_eo_Richardson<AFIELD_d, AFIELD_f>::flop_count()
{
  return m_solver->flop_count();
}


//====================================================================
template<typename AFIELD_d, typename AFIELD_f>
void Fprop_alt_Standard_eo_Richardson<AFIELD_d, AFIELD_f>::reset_performance()
{
  m_flop_count   = 0.0;
  m_elapsed_time = 0.0;
}


//====================================================================
template<typename AFIELD_d, typename AFIELD_f>
void Fprop_alt_Standard_eo_Richardson<AFIELD_d, AFIELD_f>::report_performance()
{
  double flops  = m_flop_count / m_elapsed_time;
  double gflops = flops * 1.0e-9;

  vout.general(m_vl, "\n");
  vout.general(m_vl, "%s: solver performance:\n", class_name.c_str());
  vout.general(m_vl, "  Elapsed time = %14.6f sec\n", m_elapsed_time);
  vout.general(m_vl, "  Flop(total)  = %18.0f\n", m_flop_count);
  vout.general(m_vl, "  Performance  = %11.3f GFlops\n", gflops);
}

//====================================================================
template<typename AFIELD_d, typename AFIELD_f>
void Fprop_alt_Standard_eo_Richardson<AFIELD_d, AFIELD_f>::get_performance(double& flop_count,
                                                               double& elapsed_time)
{
  //  flop_count   = m_flop_count;
  //  elapsed_time = m_elapsed_time;
}

//====================================================================
template<typename AFIELD_d, typename AFIELD_f>
void Fprop_alt_Standard_eo_Richardson<AFIELD_d, AFIELD_f>::mult_performance(
  const std::string mode,
  const int Nrepeat)
{
  int nin  = m_fopr_d->field_nin();
  int nvol = m_fopr_d->field_nvol();
  int nex  = m_fopr_d->field_nex();

  unique_ptr<Timer> timer(new Timer);

  std::string mode_prev_d = m_fopr_d->get_mode();
  std::string mode_prev_f = m_fopr_f->get_mode();

  m_fopr_d->set_mode(mode);
  m_fopr_f->set_mode(mode);

  {
    //  AFIELD_d axq(nin, nvol, nex), abq(nin, nvol, nex);
    AFIELD_d &y1 = m_be;
    AFIELD_d &y2 = m_bo;
#pragma omp parallel
    {
      y1.set(0.0);
    }
    y1.set(0, 1.0);

    timer->start();

#pragma omp parallel
    {
      for (int i = 0; i < Nrepeat; ++i) {
        m_fopr_d->mult(y2, y1);
        m_fopr_d->mult(y1, y2);
      }
    }

    timer->stop();
  }

  double flop_fopr  = m_fopr_d->flop_count();
  double flop_total = flop_fopr * double(2 * Nrepeat);

  double elapsed_time = timer->elapsed_sec();
  double flops        = flop_total / elapsed_time;
  double gflops       = flops * 1.0e-9;

  vout.general(m_vl, "\n");
  vout.general(m_vl, "%s: mult performance:\n", class_name.c_str());
  vout.general(m_vl, "  mult mode = %s\n", mode.c_str());
  vout.general(m_vl, "  Double precision:\n");
  vout.general(m_vl, "   Elapsed time = %14.6f sec\n", elapsed_time);
  vout.general(m_vl, "   Flop(Fopr)   = %18.0f\n", flop_fopr);
  vout.general(m_vl, "   Flop(total)  = %18.0f\n", flop_total);
  vout.general(m_vl, "   Performance  = %11.3f GFlops\n", gflops);

  {
    //  AFIELD_f axq(nin, nvol, nex), abq(nin, nvol, nex);
    AFIELD_f y1f;
    AFIELD_f y2f;
    int nin   = m_fopr_f->field_nin();
    int nvol2 = m_fopr_f->field_nvol();
    int nex   = m_fopr_f->field_nex();
    y1f.reset(nin, nvol2, nex);
    y2f.reset(nin, nvol2, nex);

#pragma omp parallel
    {
      y1f.set(0.0);
    }
    y2f.set(0, 1.0);

    timer->reset();
    timer->start();

#pragma omp parallel
    {
      for (int i = 0; i < Nrepeat; ++i) {
        m_fopr_f->mult(y2f, y1f);
        m_fopr_f->mult(y1f, y2f);
      }
    }

    timer->stop();
  }

  flop_fopr  = m_fopr_f->flop_count();
  flop_total = flop_fopr * double(2 * Nrepeat);

  elapsed_time = timer->elapsed_sec();
  flops        = flop_total / elapsed_time;
  gflops       = flops * 1.0e-9;

  vout.general(m_vl, "  Single precision:\n");
  vout.general(m_vl, "   Elapsed time = %14.6f sec\n", elapsed_time);
  vout.general(m_vl, "   Flop(Fopr)   = %18.0f\n", flop_fopr);
  vout.general(m_vl, "   Flop(total)  = %18.0f\n", flop_total);
  vout.general(m_vl, "   Performance  = %11.3f GFlops\n", gflops);

  m_fopr_d->set_mode(mode_prev_d);
  m_fopr_f->set_mode(mode_prev_f);
}


//============================================================END=====
