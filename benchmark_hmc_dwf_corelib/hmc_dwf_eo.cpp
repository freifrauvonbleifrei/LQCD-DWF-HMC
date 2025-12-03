/*!
        @file    $Id: hmc_dwf_eo.cpp #$

        @brief   HMC for DWF

        @author  Issaku Kanamori (kanamori)
                 $LastChangedBy: kanamori $

        @date    $LastChangedDate:: 2#$

        @version $LastChangedRevision: 2668 $
*/

#include "bridge.h"

#include "Fopr/fopr_Domainwall_eo.h"
#include "Force/Fermion/force_F_Domainwall_eo.h"
#include "Force/Fermion/force_F_Rational_Ratio.h"
#include "Action/Fermion/action_F_Rational_Ratio_eo.h"


using Bridge::vout;


namespace {
  const std::string parameter_file_main = "main.yaml";
  const std::string parameter_file = "HMC_DWF_Nf2p1_eo.yaml";
}


//====================================================================
  int update_Nf2p1_eo(const Parameters &params_all)
  {
    // #####  parameter setup  #####
    const int Nc   = CommonParameters::Nc();
    const int Nvol = CommonParameters::Nvol();
    const int Ndim = CommonParameters::Ndim();


    const Parameters params_test       = params_all.lookup("HMC_Domainwall");
    const Parameters params_action_G   = params_all.lookup("Action_G");
    const Parameters params_smear      = params_all.lookup("LinkSmearing");

    const Parameters params_pv        = params_all.lookup("Fopr_PauliVillars");
    const Parameters params_nf1       = params_all.lookup("Fopr_Nf1");
    const Parameters params_nf2       = params_all.lookup("Fopr_Nf2");
    const Parameters params_nf2_prec  = params_all.lookup("Fopr_Nf2_prec"); // Hasenbusch



    // Rational approximation in molecular dynamics
    const Parameters params_nf1_rational_MD    = params_all.lookup("Fopr_Nf1_Rational_MD");
    const Parameters params_nf1_pv_rational_MD = params_all.lookup("Fopr_Nf1_PauliVillars_Rational_MD");

    // Rational approximation in calculating Hamiltonian
    const Parameters params_nf1_rational_calcH    = params_all.lookup("Fopr_Nf1_Rational_calcH");
    const Parameters params_nf1_pv_rational_calcH = params_all.lookup("Fopr_Nf1_PauliVillars_Rational_calcH");

    const Parameters params_solver_MD  = params_all.lookup("Solver_MD");
    const Parameters params_solver_H   = params_all.lookup("Solver_H");

    const Parameters params_integrator = params_all.lookup("Builder_Integrator");
    const Parameters params_hmc        = params_all.lookup("HMC_General");

    const string        str_gconf_status = params_test.get_string("gauge_config_status");
    const string        str_gconf_read   = params_test.get_string("gauge_config_type_input");
    const string        readfile         = params_test.get_string("config_filename_input");
    const string        str_gconf_write  = params_test.get_string("gauge_config_type_output");
    const string        writefile        = params_test.get_string("config_filename_output");
    const string        str_rand_type    = params_test.get_string("random_number_type");
    const unsigned long seed             = params_test.get_unsigned_long("seed_for_random_number");
    int                 i_conf           = params_test.get_int("trajectory_number");
    const int           Ntraj            = params_test.get_int("trajectory_number_step");
    const int           i_save_conf      = params_test.get_int("save_config_interval");
    const string        str_vlevel       = params_test.get_string("verbose_level");
    //    const bool          rand_from_file   = params_test.get_bool("rand_from_file")
    const bool   do_check        = params_test.is_set("expected_result");
    const double expected_result = do_check ? params_test.get_double("expected_result") : 0.0;

    const string str_action_G_type = params_action_G.get_string("action_type");
    const string str_gmset_type    = params_nf1.get_string("gamma_matrix_type");
    const string str_solver_MD_type = params_solver_MD.get_string("solver_type");
    const string str_solver_H_type  = params_solver_H.get_string("solver_type");
    const int    Nlevels            = params_integrator.get_int("number_of_levels");
    const string str_proj_type      = params_smear.get_string("projection_type");
    const string str_smear_type     = params_smear.get_string("smear_type");

    const std::vector<int> level_action       = params_integrator.get_int_vector("level_of_actions");

    const Bridge::VerboseLevel vl = vout.set_verbose_level(str_vlevel);

    vout.crucial("parameter setup finished.\n");
    Communicator::sync();

    //- print input parameters
    vout.general(vl, "  gconf_status   = %s\n", str_gconf_status.c_str());
    vout.general(vl, "  gconf_read     = %s\n", str_gconf_read.c_str());
    vout.general(vl, "  readfile       = %s\n", readfile.c_str());
    vout.general(vl, "  gconf_write    = %s\n", str_gconf_write.c_str());
    vout.general(vl, "  writefile      = %s\n", writefile.c_str());
    vout.general(vl, "  rand_type      = %s\n", str_rand_type.c_str());
    vout.general(vl, "  seed           = %lu\n", seed);
    vout.general(vl, "  i_conf         = %d\n", i_conf);
    vout.general(vl, "  Ntraj          = %d\n", Ntraj);
    vout.general(vl, "  i_save_conf    = %d\n", i_save_conf);
    vout.general(vl, "  vlevel         = %s\n", str_vlevel.c_str());
    vout.general(vl, "  gmset_type     = %s\n", str_gmset_type.c_str());
    vout.general(vl, "  solver_MD_type = %s\n", str_solver_MD_type.c_str());
    vout.general(vl, "  solver_H_type  = %s\n", str_solver_H_type.c_str());
    vout.general(vl, "\n");

    vout.crucial("parameter output finished.\n");
    Communicator::sync();

    //- input parameter check
    int err = 0;
    err += ParameterCheck::non_NULL(str_gconf_status);
    err += ParameterCheck::non_negative(i_conf);
    err += ParameterCheck::non_negative(Ntraj);
    err += ParameterCheck::non_negative(i_save_conf);

    if (err) {
      vout.crucial(vl, "Error: input parameters have not been set\n");
      exit(EXIT_FAILURE);
    }

    RandomNumberManager::initialize(str_rand_type, seed);

    // #####  object setup  #####
    Field_G U(Nvol, Ndim);

    if (str_gconf_status == "Continue") {
      GaugeConfig(str_gconf_read).read(U, readfile);
    } else if (str_gconf_status == "Cold_start") {
      GaugeConfig("Unit").read(U);
    } else if (str_gconf_status == "Hot_start") {
      GaugeConfig("Random").read(U);
    } else {
      vout.crucial(vl, "Error: unsupported gconf status \"%s\"\n",
                   str_gconf_status.c_str());
      exit(EXIT_FAILURE);
    }

    GaugeConfig gconf_write(str_gconf_write);

    vout.crucial("gauge config. finished\n");
    Communicator::sync();

    // sanity check of the fermion action
    {
      const std::string fermion_type = params_nf1.get_string("fermion_type");
      vout.general(vl, "fermion type in params_nf1: %s\n", fermion_type.c_str());
      if(fermion_type != "Domainwall_5din_eo") {
        vout.crucial(vl, "Error: wrong fermion_type is given\n");
        exit(EXIT_FAILURE);
      };

      const std::string fermion_type_pv = params_pv.get_string("fermion_type");
      vout.general(vl, "fermion type in params_pv: %s\n", fermion_type_pv.c_str());
      if(fermion_type_pv != fermion_type){
        vout.crucial(vl, "Error: inconsistent fermion_type\n");
        exit(EXIT_FAILURE);
      }

      const std::string fermion_type_nf2 = params_nf2.get_string("fermion_type");
      vout.general(vl, "fermion type in params_nf2: %s\n", fermion_type_nf2.c_str());
      if(fermion_type_pv != fermion_type){
        vout.crucial(vl, "Error: inconsistent fermion_type\n");
        exit(EXIT_FAILURE);
      }

    }

    //- Gauge field action
    unique_ptr<Action> action_G(
                       Action::New(str_action_G_type, params_action_G));

    //- Smearing director
    // define smearing method (SA)
    unique_ptr<Projection> proj(Projection::New(str_proj_type, params_smear));
    unique_ptr<Smear> smear(Smear::New(str_smear_type, proj.get(), params_smear));

    // define force smearing method (SA)
    unique_ptr<ForceSmear> force_smear(ForceSmear::New(str_smear_type, proj.get(), params_smear));
    unique_ptr<Director_Smear> dr_smear(new Director_Smear(smear.get(), params_smear));

    ///////////////////////////////////////////////////////////
    // RHMC part
    //
    //- Domainwall part: "denominatoar"
    unique_ptr<Fopr_eo> fopr_Nf1(new Fopr_Domainwall_eo(params_nf1));
    unique_ptr<Force> forceF_Nf1_dw(new Force_F_Domainwall_eo(params_nf1));

    //- PauliVillar part: "numerator"
    unique_ptr<Fopr_eo> fopr_Nf1_pv(new Fopr_Domainwall_eo(params_pv));
    unique_ptr<Force> forceF_Nf1_pv(new Force_F_Domainwall_eo(params_pv));

    //- Ratio: force
    Parameters params_Nf1_rational_ratio_MD;
    params_Nf1_rational_ratio_MD.set_Parameters("Fopr_Rational_numerator_MD", params_nf1_pv_rational_MD);
    params_Nf1_rational_ratio_MD.set_Parameters("Fopr_Rational_denominator_MD", params_nf1_rational_MD);

    unique_ptr<Force>
      forceF_Nf1(
        new Force_F_Rational_Ratio(fopr_Nf1_pv.get(), forceF_Nf1_pv.get(),
                                   fopr_Nf1.get(), forceF_Nf1_dw.get(),
                                   params_Nf1_rational_ratio_MD)
                           );

    // smeared objects
    unique_ptr<Fopr> fopr_Nf1_smr(
         new Fopr_Smeared_eo(fopr_Nf1.get(), dr_smear.get()));

    unique_ptr<Fopr> fopr_Nf1_pv_smr(
         new Fopr_Smeared_eo(fopr_Nf1_pv.get(), dr_smear.get()));

    unique_ptr<Force> forceF_Nf1_smr(
         new Force_F_Smeared(forceF_Nf1.get(), force_smear.get(), dr_smear.get()));


    //- Domain-wall/Pauli-Villars action
    Parameters params_Nf1_rational_ratio_calcH;
    params_Nf1_rational_ratio_calcH.set_Parameters("Fopr_Rational_numerator_calcH", params_nf1_pv_rational_calcH);
    params_Nf1_rational_ratio_calcH.set_Parameters("Fopr_Rational_denominator_calcH", params_nf1_rational_calcH);

    unique_ptr<Action> action_F_Nf1(
         new Action_F_Rational_Ratio_eo(fopr_Nf1_pv_smr.get(), fopr_Nf1_smr.get(),
                                        forceF_Nf1_smr.get(),
                                        params_Nf1_rational_ratio_calcH));


    ///////////////////////////////////////////////////////////
    // Nf=2 part
    //
    // physical part
    unique_ptr<Fopr_eo > fopr_Nf2(new Fopr_Domainwall_eo(params_nf2));
    unique_ptr<Force> forceF_Nf2(new Force_F_Domainwall_eo(params_nf2));

    // Hasenbuche operator
    unique_ptr<Fopr_eo> fopr_Nf2_prec(new Fopr_Domainwall_eo(params_nf2_prec));
    unique_ptr<Force> forceF_Nf2_prec(new Force_F_Domainwall_eo(params_nf2_prec));

    // Pauli-Villars
    unique_ptr<Fopr_eo> fopr_Nf2_pv(new Fopr_Domainwall_eo(params_pv));
    unique_ptr<Force> forceF_Nf2_pv(new Force_F_Domainwall_eo(params_pv));

    // smeared objects
    unique_ptr<Fopr>  fopr_Nf2_smr(Fopr::New("Smeared_eo", fopr_Nf2.get(), dr_smear.get()));
    unique_ptr<Fopr>  fopr_Nf2_prec_smr(Fopr::New("Smeared_eo", fopr_Nf2_prec.get(), dr_smear.get()));
    unique_ptr<Fopr>  fopr_Nf2_pv_smr(Fopr::New("Smeared_eo", fopr_Nf2_pv.get(), dr_smear.get()));

    unique_ptr<Force> forceF_Nf2_smr(
       new Force_F_Smeared(forceF_Nf2.get(), force_smear.get(), dr_smear.get()));
    unique_ptr<Force> forceF_Nf2_prec_smr(
       new Force_F_Smeared(forceF_Nf2_prec.get(), force_smear.get(), dr_smear.get()));
    unique_ptr<Force> forceF_Nf2_pv_smr(
       new Force_F_Smeared(forceF_Nf2_pv.get(), force_smear.get(), dr_smear.get()));


    // solvers
    unique_ptr<Solver> solver_Nf2_MD(
       Solver::New(str_solver_MD_type, fopr_Nf2_smr.get(), params_solver_MD));
    unique_ptr<Fprop> fprop_Nf2_MD(
       new Fprop_Standard_eo(solver_Nf2_MD.get()));

    unique_ptr<Solver> solver_Nf2_prec_MD(
       Solver::New(str_solver_MD_type, fopr_Nf2_prec_smr.get(), params_solver_MD));
    unique_ptr<Fprop> fprop_Nf2_prec_MD(
       new Fprop_Standard_eo(solver_Nf2_prec_MD.get()));

    unique_ptr<Solver> solver_Nf2_H(
       Solver::New(str_solver_H_type, fopr_Nf2_smr.get(), params_solver_H));
    unique_ptr<Fprop> fprop_Nf2_H(
       new Fprop_Standard_eo(solver_Nf2_H.get()));

    unique_ptr<Solver> solver_Nf2_prec_H(
       Solver::New(str_solver_H_type, fopr_Nf2_prec_smr.get(), params_solver_H));
    unique_ptr<Fprop> fprop_Nf2_prec_H(
       new Fprop_Standard_eo(solver_Nf2_prec_H.get()));

    unique_ptr<Solver> solver_Nf2_pv_H(
       Solver::New(str_solver_H_type, fopr_Nf2_pv_smr.get(), params_solver_H));
    unique_ptr<Fprop> fprop_Nf2_pv_H(
       new Fprop_Standard_eo(solver_Nf2_pv_H.get()));

    //- (prec/phys)
    unique_ptr<Action> action_F_Nf2(new Action_F_Ratio_eo(
       fopr_Nf2_prec_smr.get(), forceF_Nf2_prec_smr.get(),
       fopr_Nf2_smr.get(),      forceF_Nf2_smr.get(),
       fprop_Nf2_prec_H.get(),
       fprop_Nf2_MD.get(), fprop_Nf2_H.get() ));


    //- (PV/prec)
    unique_ptr<Action> action_F_Nf2_prec(new Action_F_Ratio_eo(
       fopr_Nf2_pv_smr.get(),   forceF_Nf2_pv_smr.get(),
       fopr_Nf2_prec_smr.get(), forceF_Nf2_prec_smr.get(),
       fprop_Nf2_pv_H.get(),
       fprop_Nf2_prec_MD.get(), fprop_Nf2_prec_H.get() ));


    ///////////////////////////////////////////////////////////
    // register actions
    //
    if(level_action.size()<3){
      vout.crucial(vl, "level_of_actions are too small\n");
      exit(EXIT_FAILURE);
    }
    ActionList actions(Nlevels);
    actions.append(level_action[0], action_F_Nf1.get());
    actions.append(level_action[1], action_F_Nf2.get());
    actions.append(level_action[2], action_F_Nf2_prec.get());
    actions.append(level_action[3], action_G.get());

    std::vector<Director *> directors(1);
    directors[0] = static_cast<Director *>(dr_smear.get());

    unique_ptr<Builder_Integrator> builder(
      new Builder_Integrator(actions, directors, params_integrator));
    Integrator *integrator = builder->build();

    //    unique_ptr<RandomNumbers> rand(new RandomNumbers_Mseries(i_conf));
    //   unique_ptr<RandomNumbers> rand(new RandomNumbers_parallel<RandomNumbers_Mseries>(i_conf));
    RandomNumbers *rand = RandomNumberManager::getInstance();

    HMC_General hmc(actions, directors, integrator, rand, params_hmc);

    vout.crucial(vl, "Object construction finished.\n");
    Communicator::sync();


    Timer timer("update_Nf2p1_eo");


    // ####  Execution main part  ####
    timer.start();

    vout.general(vl, "HMC: Ntraj = %d\n", Ntraj);

    double result = 0.0;
    for (int traj = 0; traj < Ntraj; ++traj) {
      Timer timer_traj("each trajectory");
      timer_traj.start();
      vout.general(vl, "\n");
      vout.general(vl, "---------------------------------------------------\n");
      vout.general(vl, "traj = %d\n", traj);
      rand->reset(i_conf+traj);

      result = hmc.update(U);

      timer_traj.stop();
      timer_traj.report();

      if ((i_conf + traj + 1) % i_save_conf == 0) {
        std::string filename = FileUtils::generate_filename("%s-%06d", writefile.c_str(), (i_conf + traj + 1));
        gconf_write.write_file(U, filename);
      }
    }

    gconf_write.write_file(U, writefile);

    timer.report();

    RandomNumberManager::finalize();


    return 0;
  }

//====================================================================
int main(int argc, char *argv[])
{
  // ###  initial setup  ###
  bridge_initialize(&argc, &argv);
  Bridge::VerboseLevel vl = Bridge::GENERAL;

  Parameters params_main = ParameterManager::read(parameter_file_main);
  bridge_setup(params_main.lookup("Main"));

  Parameters params = ParameterManager::read(parameter_file);


  Timer timer("Main");
  timer.start();

  update_Nf2p1_eo(params);

  timer.stop();
  timer.report();

  bridge_finalize();

  return EXIT_SUCCESS;
}
