/*!
        @file    field_F_imp_SU2-inc.h

        @brief

        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate:: 2024-12-31 19:58:26 #$

        @version $LastChangedRevision: 2624 $
*/

// This implementation only applies to SU(2) group and Nd=4 case.
#define NC      2
#define NC2     4
#define NDF     8
#define ND      4
#define NCD     8
#define NCD2    16

//====================================================================
namespace {
  void check_Nc()
  {
    vout.paranoiac(CommonParameters::Vlevel(),
                   "Field_F: implementation for SU(2).\n");
  }


  double mult_Gn_r(const double *g, const double *w, int Nc)
  {
    return g[0] * w[0] - g[1] * w[1]
           + g[2] * w[2] - g[3] * w[3];
  }


  double mult_Gn_i(const double *g, const double *w, int Nc)
  {
    return g[0] * w[1] + g[1] * w[0]
           + g[2] * w[3] + g[3] * w[2];
  }


  double mult_Gd_r(const double *g, const double *w, int Nc)
  {
    return g[0] * w[0] + g[1] * w[1]
           + g[4] * w[2] + g[5] * w[3];
  }


  double mult_Gd_i(const double *g, const double *w, int Nc)
  {
    return g[0] * w[1] - g[1] * w[0]
           + g[4] * w[3] - g[5] * w[2];
  }
} // end of nameless namespace
//====================================================================
//============================================================END=====
