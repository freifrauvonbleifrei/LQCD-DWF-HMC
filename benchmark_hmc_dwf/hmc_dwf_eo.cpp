/*!
        @file    $Id: hmc_dwf_eo.cpp #$

        @brief   HMC for DWF

        @author  Issaku Kanamori (kanamori)
                 $LastChangedBy: kanamori $

        @date    $LastChangedDate:: 2#$

        @version $LastChangedRevision: 2668 $
*/

#include "lib/bridge_setup.h"
#include "lib/fapp_macros.h"

// corelib
#include "lib/Communicator/communicator.h"
#include "lib/ResourceManager/threadManager.h"
#include "lib/Parameters/commonParameters.h"
#include "lib/Parameters/parameters.h"
#include "lib/Parameters/parameterManager_YAML.h"
#include "lib/Tools/timer.h"
#include "lib/bridge_init_factory.h"
#include "lib/IO/bridgeIO.h"
#include "lib/IO/gaugeConfig.h"


#include "lib/Fopr/fopr_Smeared.h"
#include "lib/Fopr/fopr_Smeared_alt.h"

#include "lib/HMC/hmc_General.h"
#include "lib/HMC/builder_Integrator.h"

#include "lib/Tools/file_utils.h"
#include "lib/Tools/randomNumberManager.h"
#include "lib/Tools/randomNumbers_parallel.h"
#include "lib/Tools/randomNumbers_Mseries.h"


#ifdef USE_ALT_CODE
#define DIRECTOR_ALT_SMEARED_IMPLEMENTED
#include "lib_alt/bridge_alt_init.h"
#include "lib_alt/Measurements/Fermion/fprop_alt_Standard_eo_Richardson.h"
#include "lib_alt/Smear/director_alt_Smear.h"
#include "lib_alt/Action/Fermion/action_F_alt_Ratio_eo.h"
#include "lib_alt/Action/Fermion/action_F_alt_Rational_Ratio.h"
#include "lib_alt/Force/Fermion/aforce_F_Rational_Ratio.h"
#include "lib_alt/Force/Fermion/aforce_F_Ratio_eo.h"
#include "lib_alt/Force/Fermion/aforce_F_Smeared_alt.h"
#endif

#ifdef USE_ALT_QXS
#include "lib_alt_QXS/bridge_alt_qxs.h"
#include "lib_alt_QXS/Fopr/afopr_Domainwall_5din_eo.h"
#include "lib_alt_QXS/Force/Fermion/aforce_F_Domainwall_5din_eo.h"
#define IMPL QXS
#endif

#ifdef USE_ALT_ACCEL
#include "lib_alt_Accel/bridge_alt_accel.h"
#include "lib_alt_Accel/Fopr/afopr_Domainwall_5din_eo.h"
#include "lib_alt_Accel/Force/Fermion/aforce_F_Domainwall_5din_eo.h"
#define IMPL ACCEL
#endif


using Bridge::vout;


namespace {
  const std::string parameter_file_main = "main.yaml";
  const std::string parameter_file = "HMC_DWF_Nf2p1_eo.yaml";
}


typedef  AField<double,IMPL>  AFIELD;
typedef  AField<float,IMPL>  AFIELD_f;
typedef  Fprop_alt_Standard_eo_Richardson<AFIELD, AFIELD_f> Fprop_alt_eo;

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
    typedef Director_alt_Smear<AFIELD>  Director_alt;
    unique_ptr<Director_alt> dr_smear(new Director_alt(params_smear));



    ///////////////////////////////////////////////////////////
    // PauliVillar operator
    //
    unique_ptr<AFopr<AFIELD> > fopr_pv(
                           AFopr<AFIELD>::New("Domainwall_5din_eo", params_pv));

    unique_ptr<AForce_F<AFIELD> > forceF_pv(
                           new AForce_F_Domainwall_5din_eo<AFIELD>(params_pv));

    unique_ptr<AFopr<AFIELD> >  fopr_pv_smr(
       AFopr<AFIELD>::New("Smeared_alt", fopr_pv.get(), dr_smear.get()));


    ///////////////////////////////////////////////////////////
    // RHMC part
    //
    //- Domainwall part: "denominatoar"
    unique_ptr<AFopr<AFIELD> > fopr_Nf1(
                           AFopr<AFIELD>::New("Domainwall_5din_eo", params_nf1));
    unique_ptr<AForce_F<AFIELD> > forceF_Nf1_dw(
                           new AForce_F_Domainwall_5din_eo<AFIELD>(params_nf1));

    //- Ratio: force
    Parameters params_Nf1_rational_ratio_MD;
    params_Nf1_rational_ratio_MD.set_Parameters("Fopr_Rational_numerator_MD", params_nf1_pv_rational_MD);
    params_Nf1_rational_ratio_MD.set_Parameters("Fopr_Rational_denominator_MD", params_nf1_rational_MD);

    unique_ptr<AForce_F_Rational_Ratio<AFIELD> > forceF_Nf1(
            new AForce_F_Rational_Ratio<AFIELD>(
                             fopr_pv.get(), forceF_pv.get(),
                             fopr_Nf1.get(), forceF_Nf1_dw.get(),
                             params_Nf1_rational_ratio_MD));


    // smeared objects
    unique_ptr<AFopr<AFIELD> >  fopr_Nf1_smr(
       AFopr<AFIELD>::New("Smeared_alt", fopr_Nf1.get(), dr_smear.get()));

    unique_ptr<AForce_F<AFIELD> > forceF_Nf1_smr(new
            AForce_F_Smeared_alt<AFIELD>(forceF_Nf1.get(), dr_smear.get()));


    //- Domain-wall/Pauli-Villars action
    Parameters params_Nf1_rational_ratio_calcH;
    params_Nf1_rational_ratio_calcH.set_Parameters("Fopr_Rational_numerator_calcH", params_nf1_pv_rational_calcH);
    params_Nf1_rational_ratio_calcH.set_Parameters("Fopr_Rational_denominator_calcH", params_nf1_rational_calcH);

    unique_ptr<Action> action_F_Nf1(new Action_F_alt_Rational_Ratio<AFIELD>(
                             fopr_pv_smr.get(), fopr_Nf1_smr.get(),
                             forceF_Nf1_smr.get(),
                             params_Nf1_rational_ratio_calcH));


    ///////////////////////////////////////////////////////////
    // Nf=2 part
    //
    // physical part
    unique_ptr<AFopr<AFIELD> > fopr_Nf2(
                           AFopr<AFIELD>::New("Domainwall_5din_eo", params_nf2));

    unique_ptr<AForce_F<AFIELD> > forceF_Nf2(
                          new AForce_F_Domainwall_5din_eo<AFIELD>(params_nf2));

    // Hasenbuche operator
    unique_ptr<AFopr<AFIELD> > fopr_Nf2_prec(
                           AFopr<AFIELD>::New("Domainwall_5din_eo", params_nf2_prec));

    unique_ptr<AForce_F<AFIELD> > forceF_Nf2_prec(
                          new AForce_F_Domainwall_5din_eo<AFIELD>(params_nf2_prec));

    // smeared objects
    unique_ptr<AFopr<AFIELD> >  fopr_Nf2_smr(
       AFopr<AFIELD>::New("Smeared_alt", fopr_Nf2.get(), dr_smear.get()));

    unique_ptr<AFopr<AFIELD> >  fopr_Nf2_prec_smr(
       AFopr<AFIELD>::New("Smeared_alt", fopr_Nf2_prec.get(), dr_smear.get()));


    // solvers
    unique_ptr<Fprop_alt<AFIELD> > fprop_Nf2_MD(new
       Fprop_alt_eo(params_nf2, params_solver_MD, dr_smear.get()));

    unique_ptr<Fprop_alt<AFIELD> > fprop_Nf2_prec_MD(new
       Fprop_alt_eo(params_nf2_prec, params_solver_MD, dr_smear.get()));

    unique_ptr<Fprop_alt<AFIELD> > fprop_Nf2_H(new
       Fprop_alt_eo(params_nf2, params_solver_H, dr_smear.get()));

    unique_ptr<Fprop_alt<AFIELD> > fprop_Nf2_prec_H(new
       Fprop_alt_eo(params_nf2_prec, params_solver_H, dr_smear.get()));

    unique_ptr<Fprop_alt<AFIELD> > fprop_Nf2_pv_H(new
       Fprop_alt_eo(params_pv, params_solver_H, dr_smear.get()));


    //- (prec/phys)
    unique_ptr<AForce_F<AFIELD> > forceF_Nf2_ratio(
       new AForce_F_Ratio_eo<AFIELD>(fopr_Nf2_prec.get(),
                                     forceF_Nf2_prec.get(), forceF_Nf2.get()) );
    unique_ptr<AForce_F<AFIELD>> forceF_Nf2_ratio_smr(new
       AForce_F_Smeared_alt<AFIELD>(forceF_Nf2_ratio.get(), dr_smear.get()));

    unique_ptr<Action> action_F_Nf2(new Action_F_alt_Ratio_eo<AFIELD>(
                             fopr_Nf2_prec_smr.get(),  fopr_Nf2_smr.get(),
                             fprop_Nf2_prec_H.get(),
                             fprop_Nf2_MD.get(), fprop_Nf2_H.get(),
                             forceF_Nf2_ratio_smr.get() ));

    //- (PV/prec)
    unique_ptr<AForce_F<AFIELD> > forceF_Nf2_ratio_prec(
       new AForce_F_Ratio_eo<AFIELD>(fopr_pv.get(),
                                     forceF_pv.get(), forceF_Nf2_prec.get()) );
    unique_ptr<AForce_F<AFIELD>> forceF_Nf2_ratio_prec_smr(new
       AForce_F_Smeared_alt<AFIELD>(forceF_Nf2_ratio_prec.get(), dr_smear.get()));


    unique_ptr<Action> action_F_Nf2_prec(new Action_F_alt_Ratio_eo<AFIELD>(
                             fopr_pv_smr.get(), fopr_Nf2_prec_smr.get(),
                             fprop_Nf2_pv_H.get(),
                             fprop_Nf2_prec_MD.get(), fprop_Nf2_prec_H.get(),
                             forceF_Nf2_ratio_prec_smr.get() ));

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

      START_FAPP("update", 0, 0);
      result = hmc.update(U);
      STOP_FAPP("update", 0, 0);

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

#ifdef USE_ALT_CODE
  bridge_alt_init(params_main);
  vout.crucial(vl, "bridge_alt_init finished.\n");
#endif

  Parameters params = ParameterManager::read(parameter_file);


  Timer timer("Main");
  timer.start();

  update_Nf2p1_eo(params);

  timer.stop();
  timer.report();

  bridge_finalize();

  return EXIT_SUCCESS;
}
