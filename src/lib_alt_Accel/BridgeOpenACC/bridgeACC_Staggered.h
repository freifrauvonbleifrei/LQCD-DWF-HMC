/*!
      @file    bridgeACC_Staggered.h
      @brief
      @author  Hideo Matsufuru (matufuru)
      @date    $LastChangedDate: 2013-01-22 13:51:53 #$
      @version $LastChangedRevision: 2653 $
*/

#ifndef BRIDGEACC_STAGGERED_INCLUDED
#define BRIDGEACC_STAGGERED_INCLUDED

#include "bridge_defs.h"

namespace BridgeACC {

  // real_t = double

  void mult_staggered_D(double *RESTRICT v2, double *RESTRICT u_up,
                        double *RESTRICT u_dn, double *RESTRICT v1,
                        double mq, int *bc, int *Nsize, int jdag);

  void mult_staggered_1(double *RESTRICT buf_xp, double *RESTRICT buf_xm,
                        double *RESTRICT buf_yp, double *RESTRICT buf_ym,
                        double *RESTRICT buf_zp, double *RESTRICT buf_zm,
                        double *RESTRICT buf_tp, double *RESTRICT buf_tm,
                        double *RESTRICT u_dn, double *RESTRICT v1,
                        int *bc, int *Nsize, int *do_comm);

  void mult_staggered_2(double *RESTRICT v2, double *RESTRICT u_up,
                        double *RESTRICT buf_xp, double *RESTRICT buf_xm,
                        double *RESTRICT buf_yp, double *RESTRICT buf_ym,
                        double *RESTRICT buf_zp, double *RESTRICT buf_zm,
                        double *RESTRICT buf_tp, double *RESTRICT buf_tm,
                        double mq, int *bc,
                        int *Nsize, int *do_comm, int jdag);

  void mult_staggered_Meo(double *RESTRICT v2,
			  double *RESTRICT u_up, double *RESTRICT u_dn,
			  double *RESTRICT v1, double *RESTRICT x1, 
			  double mq2, int *bc, int *Nsize,
			  const int ieo, const int jeo, const int iflag);

  void mult_staggered_1eo(double *RESTRICT buf_xp, double *RESTRICT buf_xm,
                          double *RESTRICT buf_yp, double *RESTRICT buf_ym,
                          double *RESTRICT buf_zp, double *RESTRICT buf_zm,
                          double *RESTRICT buf_tp, double *RESTRICT buf_tm,
                          double *RESTRICT u, double *RESTRICT v1,
                          int *bc, int *Nsize, int *do_comm,
                          const int ieo, const int jeo);

  void mult_staggered_2eo(double *RESTRICT v2, double *RESTRICT u, 
                          double *RESTRICT buf_xp, double *RESTRICT buf_xm,
                          double *RESTRICT buf_yp, double *RESTRICT buf_ym,
                          double *RESTRICT buf_zp, double *RESTRICT buf_zm,
                          double *RESTRICT buf_tp, double *RESTRICT buf_tm,
                          double mq, int *bc,  int *Nsize, int *do_comm,
                          const int ieo, const int jeo, const int iflag);

  void mult_staggered_phase(double *RESTRICT u, double *RESTRICT ph,
                            int *Nsize, int Nc);

  void mult_staggered_gm5(double *RESTRICT v2, double *RESTRICT v1,
                          double *RESTRICT prty, int *Nsize, int Nc);

  void mult_staggered_gm5(double *RESTRICT v2,
                          double *RESTRICT prty, int *Nsize, int Nc);

  void mult_staggered_xp1(double *RESTRICT buf, double *RESTRICT v1,
                          int *Nsize, int *bc, int Nc);

  void mult_staggered_xp2(double *RESTRICT v2, double *RESTRICT u,
                          double *RESTRICT buf, 
                          int *Nsize, int *bc, int Nc);

  void mult_staggered_xpb(double *RESTRICT v2, double *RESTRICT u,
                          double *RESTRICT v1,
                          int *Nsize, int *bc, int Nc);

  void mult_staggered_xm1(double *RESTRICT buf, double *RESTRICT u,
                          double *RESTRICT v1,
                          int *Nsize, int *bc, int Nc);

  void mult_staggered_xm2(double *RESTRICT v2, double *RESTRICT buf, 
                          int *Nsize, int *bc, int Nc);

  void mult_staggered_xmb(double *RESTRICT v2, double *RESTRICT u,
                          double *RESTRICT v1,
                          int *Nsize, int *bc, int Nc);

  void mult_staggered_yp1(double *RESTRICT buf, double *RESTRICT v1,
                          int *Nsize, int *bc, int Nc);

  void mult_staggered_yp2(double *RESTRICT v2, double *RESTRICT u,
                          double *RESTRICT buf, 
                          int *Nsize, int *bc, int Nc);

  void mult_staggered_ypb(double *RESTRICT v2, double *RESTRICT u,
                          double *RESTRICT v1,
                          int *Nsize, int *bc, int Nc);

  void mult_staggered_ym1(double *RESTRICT buf, double *RESTRICT u,
                          double *RESTRICT v1,
                          int *Nsize, int *bc, int Nc);

  void mult_staggered_ym2(double *RESTRICT v2, double *RESTRICT buf, 
                          int *Nsize, int *bc, int Nc);

   void mult_staggered_ymb(double *RESTRICT v2, double *RESTRICT u,
                           double *RESTRICT v1,
                           int *Nsize, int *bc, int Nc);

  void mult_staggered_zp1(double *RESTRICT buf, double *RESTRICT v1,
                          int *Nsize, int *bc, int Nc);

  void mult_staggered_zp2(double *RESTRICT v2, double *RESTRICT u,
                          double *RESTRICT buf, 
                          int *Nsize, int *bc, int Nc);

  void mult_staggered_zpb(double *RESTRICT v2, double *RESTRICT u,
                          double *RESTRICT v1,
                          int *Nsize, int *bc, int Nc);

  void mult_staggered_zm1(double *RESTRICT buf, double *RESTRICT u,
                          double *RESTRICT v1,
                          int *Nsize, int *bc, int Nc);

  void mult_staggered_zm2(double *RESTRICT v2, double *RESTRICT buf, 
                          int *Nsize, int *bc, int Nc);

  void mult_staggered_zmb(double *RESTRICT v2, double *RESTRICT u,
                          double *RESTRICT v1,
                          int *Nsize, int *bc, int Nc);

  void mult_staggered_tp1(double *RESTRICT buf, double *RESTRICT v1,
                          int *Nsize, int *bc, int Nc);

  void mult_staggered_tp2(double *RESTRICT v2, double *RESTRICT u,
                          double *RESTRICT buf, 
                          int *Nsize, int *bc, int Nc);

  void mult_staggered_tpb(double *RESTRICT v2, double *RESTRICT u,
                          double *RESTRICT v1,
                          int *Nsize, int *bc, int Nc);

  void mult_staggered_tm1(double *RESTRICT buf, double *RESTRICT u,
                          double *RESTRICT v1,
                          int *Nsize, int *bc, int Nc);

  void mult_staggered_tm2(double *RESTRICT v2, double *RESTRICT buf, 
                          int *Nsize, int *bc, int Nc);

  void mult_staggered_tmb(double *RESTRICT v2, double *RESTRICT u,
                          double *RESTRICT v1,
                          int *Nsize, int *bc, int Nc);

  void mult_staggered_xp1_eo(double *RESTRICT buf, double *RESTRICT v1,
                             int *Nsize, int *bc, int ieo, int Nc);

  void mult_staggered_xp2_eo(double *RESTRICT v2, double *RESTRICT u,
                             double *RESTRICT buf, 
                             int *Nsize, int *bc, int ieo, int Nc);

  void mult_staggered_xpb_eo(double *RESTRICT v2, double *RESTRICT u,
                             double *RESTRICT v1,
                             int *Nsize, int *bc, int ieo, int Nc);

  void mult_staggered_xm1_eo(double *RESTRICT buf, double *RESTRICT u,
                             double *RESTRICT v1,
                             int *Nsize, int *bc, int ieo, int Nc);

  void mult_staggered_xm2_eo(double *RESTRICT v2, double *RESTRICT buf, 
                             int *Nsize, int *bc, int ieo, int Nc);

  void mult_staggered_xmb_eo(double *RESTRICT v2, double *RESTRICT u,
                             double *RESTRICT v1,
                             int *Nsize, int *bc, int ieo, int Nc);

  // real_t = float

  void mult_staggered_D(float *RESTRICT v2, float *RESTRICT u_up,
                        float *RESTRICT u_dn, float *RESTRICT v1,
                        float mq, int *bc, int *Nsize, int jdag);

  void mult_staggered_1(float *RESTRICT buf_xp, float *RESTRICT buf_xm,
                        float *RESTRICT buf_yp, float *RESTRICT buf_ym,
                        float *RESTRICT buf_zp, float *RESTRICT buf_zm,
                        float *RESTRICT buf_tp, float *RESTRICT buf_tm,
                        float *RESTRICT u_dn, float *RESTRICT v1,
                        int *bc, int *Nsize, int *do_comm);

  void mult_staggered_2(float *RESTRICT v2, float *RESTRICT u_up,
                        float *RESTRICT buf_xp, float *RESTRICT buf_xm,
                        float *RESTRICT buf_yp, float *RESTRICT buf_ym,
                        float *RESTRICT buf_zp, float *RESTRICT buf_zm,
                        float *RESTRICT buf_tp, float *RESTRICT buf_tm,
                        float mq, int *bc,
                        int *Nsize, int *do_comm, int jdag);

  void mult_staggered_Meo(float *RESTRICT v2,
			  float *RESTRICT u_up, float *RESTRICT u_dn,
			  float *RESTRICT v1, float *RESTRICT x1, 
			  float mq2, int *bc, int *Nsize,
			  const int ieo, const int jeo, const int iflag);

  void mult_staggered_1eo(float *RESTRICT buf_xp, float *RESTRICT buf_xm,
                          float *RESTRICT buf_yp, float *RESTRICT buf_ym,
                          float *RESTRICT buf_zp, float *RESTRICT buf_zm,
                          float *RESTRICT buf_tp, float *RESTRICT buf_tm,
                          float *RESTRICT u, float *RESTRICT v1,
                          int *bc, int *Nsize, int *do_comm,
                          const int ieo, const int jeo);

  void mult_staggered_2eo(float *RESTRICT v2, float *RESTRICT u, 
                          float *RESTRICT buf_xp, float *RESTRICT buf_xm,
                          float *RESTRICT buf_yp, float *RESTRICT buf_ym,
                          float *RESTRICT buf_zp, float *RESTRICT buf_zm,
                          float *RESTRICT buf_tp, float *RESTRICT buf_tm,
                          float mq, int *bc,  int *Nsize, int *do_comm,
                          const int ieo, const int jeo, const int iflag);

  void mult_staggered_phase(float *RESTRICT u, float *RESTRICT ph,
                            int *Nsize, int Nc);

  void mult_staggered_gm5(float *RESTRICT v2, float *RESTRICT v1,
                          float *RESTRICT prty, int *Nsize, int Nc);

  void mult_staggered_gm5(float *RESTRICT v2,
                          float *RESTRICT prty, int *Nsize, int Nc);

  void mult_staggered_xp1(float *RESTRICT buf, float *RESTRICT v1,
                          int *Nsize, int *bc, int Nc);

  void mult_staggered_xp2(float *RESTRICT v2, float *RESTRICT u,
                          float *RESTRICT buf, 
                          int *Nsize, int *bc, int Nc);

  void mult_staggered_xpb(float *RESTRICT v2, float *RESTRICT u,
                          float *RESTRICT v1,
                          int *Nsize, int *bc, int Nc);

  void mult_staggered_xm1(float *RESTRICT buf, float *RESTRICT u,
                          float *RESTRICT v1,
                          int *Nsize, int *bc, int Nc);

  void mult_staggered_xm2(float *RESTRICT v2, float *RESTRICT buf, 
                          int *Nsize, int *bc, int Nc);

  void mult_staggered_xmb(float *RESTRICT v2, float *RESTRICT u,
                          float *RESTRICT v1,
                          int *Nsize, int *bc, int Nc);

  void mult_staggered_yp1(float *RESTRICT buf, float *RESTRICT v1,
                          int *Nsize, int *bc, int Nc);

  void mult_staggered_yp2(float *RESTRICT v2, float *RESTRICT u,
                          float *RESTRICT buf, 
                          int *Nsize, int *bc, int Nc);

  void mult_staggered_ypb(float *RESTRICT v2, float *RESTRICT u,
                          float *RESTRICT v1,
                          int *Nsize, int *bc, int Nc);

  void mult_staggered_ym1(float *RESTRICT buf, float *RESTRICT u,
                          float *RESTRICT v1,
                          int *Nsize, int *bc, int Nc);

  void mult_staggered_ym2(float *RESTRICT v2, float *RESTRICT buf, 
                          int *Nsize, int *bc, int Nc);

   void mult_staggered_ymb(float *RESTRICT v2, float *RESTRICT u,
                           float *RESTRICT v1,
                           int *Nsize, int *bc, int Nc);

  void mult_staggered_zp1(float *RESTRICT buf, float *RESTRICT v1,
                          int *Nsize, int *bc, int Nc);

  void mult_staggered_zp2(float *RESTRICT v2, float *RESTRICT u,
                          float *RESTRICT buf, 
                          int *Nsize, int *bc, int Nc);

   void mult_staggered_zpb(float *RESTRICT v2, float *RESTRICT u,
                           float *RESTRICT v1,
                           int *Nsize, int *bc, int Nc);

  void mult_staggered_zm1(float *RESTRICT buf, float *RESTRICT u,
                          float *RESTRICT v1,
                          int *Nsize, int *bc, int Nc);

  void mult_staggered_zm2(float *RESTRICT v2, float *RESTRICT buf, 
                          int *Nsize, int *bc, int Nc);

  void mult_staggered_zmb(float *RESTRICT v2, float *RESTRICT u,
                          float *RESTRICT v1,
                          int *Nsize, int *bc, int Nc);

  void mult_staggered_tp1(float *RESTRICT buf, float *RESTRICT v1,
                          int *Nsize, int *bc, int Nc);

  void mult_staggered_tp2(float *RESTRICT v2, float *RESTRICT u,
                          float *RESTRICT buf, 
                          int *Nsize, int *bc, int Nc);

  void mult_staggered_tpb(float *RESTRICT v2, float *RESTRICT u,
                          float *RESTRICT v1,
                          int *Nsize, int *bc, int Nc);

  void mult_staggered_tm1(float *RESTRICT buf, float *RESTRICT u,
                          float *RESTRICT v1,
                          int *Nsize, int *bc, int Nc);

  void mult_staggered_tm2(float *RESTRICT v2, float *RESTRICT buf, 
                          int *Nsize, int *bc, int Nc);

  void mult_staggered_tmb(float *RESTRICT v2, float *RESTRICT u,
                          float *RESTRICT v1,
                          int *Nsize, int *bc, int Nc);

  void mult_staggered_xp1_eo(float *RESTRICT buf, float *RESTRICT v1,
                             int *Nsize, int *bc, int ieo, int Nc);

  void mult_staggered_xp2_eo(float *RESTRICT v2, float *RESTRICT u,
                             float *RESTRICT buf, 
                             int *Nsize, int *bc, int ieo, int Nc);

  void mult_staggered_xpb_eo(float *RESTRICT v2, float *RESTRICT u,
                             float *RESTRICT v1,
                             int *Nsize, int *bc, int ieo, int Nc);

  void mult_staggered_xm1_eo(float *RESTRICT buf, float *RESTRICT u,
                             float *RESTRICT v1,
                             int *Nsize, int *bc, int ieo, int Nc);

  void mult_staggered_xm2_eo(float *RESTRICT v2, float *RESTRICT buf, 
                             int *Nsize, int *bc, int ieo, int Nc);

  void mult_staggered_xmb_eo(float *RESTRICT v2, float *RESTRICT u,
                             float *RESTRICT v1,
                             int *Nsize, int *bc, int ieo, int Nc);

}

#endif
