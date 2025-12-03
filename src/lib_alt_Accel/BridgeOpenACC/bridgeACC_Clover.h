/*!
      @file    bridgeACC_Clover.h
      @brief
      @author  Hideo Matsufuru (matufuru)
      @date    $LastChangedDate: 2013-01-22 13:51:53 #$
      @version $LastChangedRevision: 2653 $
*/

#ifndef BRIDGEACC_CLOVER_INCLUDED
#define BRIDGEACC_CLOVER_INCLUDED

namespace BridgeACC {

  // real_t = double

  void mult_clover_D_dirac(double *RESTRICT v2, double *RESTRICT u,
                           double *RESTRICT ct, double *RESTRICT v1,
                           double kappa, int *Nsize, int *bc, int flag);

  void mult_clover_D_chiral(double *RESTRICT v2, double *RESTRICT u,
                            double *RESTRICT ct, double *RESTRICT v1,
                            double kappa, int *Nsize, int *bc, int flag);

  void mult_csw_dirac(double *RESTRICT v2, double *RESTRICT u,
                      double *RESTRICT v1, int *Nsize, int flag);

  void mult_csw_chiral(double *RESTRICT v2, double *RESTRICT u,
                       double *RESTRICT v1, int *Nsize, int flag);

  // real_t = float
  
  void mult_clover_D_dirac(float *RESTRICT v2, float *RESTRICT u,
                           float *RESTRICT ct, float *RESTRICT v1,
                           float kappa, int *Nsize, int *bc, int flag);

  void mult_clover_D_chiral(float *RESTRICT v2, float *RESTRICT u,
                            float *RESTRICT ct, float *RESTRICT v1,
                            float kappa, int *Nsize, int *bc, int flag);

  void mult_csw_dirac(float *RESTRICT v2, float *RESTRICT u,
                      float *RESTRICT v1, int *Nsize, int flag);

  void mult_csw_chiral(float *RESTRICT v2, float *RESTRICT u,
                       float *RESTRICT v1, int *Nsize, int flag);

}

#endif
//============================================================END=====
