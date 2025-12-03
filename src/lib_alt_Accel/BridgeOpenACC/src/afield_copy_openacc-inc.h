/*!
      @file    afield_copy_openacc-inc.h
      @brief
      @author  Hideo Matsufuru (matufuru)
               $LastChangedBy: matufuru $
      @date    $LastChangedDate:: 2025-08-20 16:24:03 #$
      @version $LastChangedRevision: 2653 $
*/

//====================================================================
void copy(double* RESTRICT v, int nv1, float* RESTRICT w, int nv2,
              int nin, int Nst)
{
  int Nst_pad = CEIL_NWP(Nst);
  int nv = nin * Nst_pad;

#pragma acc data present(v[nv1:nv], w[nv2:nv]) \
                 copyin(nin, Nst, nv1, nv2)

#pragma acc parallel num_workers(NUM_WORKERS) vector_length(VECTOR_LENGTH)
  {
#pragma acc loop gang worker vector
    for(int ist = 0; ist < Nst; ++ist){
      for(int in = 0; in < nin; ++in){
        v[nv1 + IDX2(nin, in, ist)] = w[nv2 + IDX2(nin, in, ist)];
      }
    }
  }

}

//====================================================================
void copy(float* RESTRICT v, int nv1, double* RESTRICT w, int nv2,
          int nin, int Nst)
{
  int Nst_pad = CEIL_NWP(Nst);
  int nv = nin * Nst_pad;

#pragma acc data present(v[nv1:nv], w[nv2:nv]) \
                 copyin(nin, Nst, nv1, nv2)

#pragma acc parallel num_workers(NUM_WORKERS) vector_length(VECTOR_LENGTH)
  {
#pragma acc loop gang worker vector
    for(int ist = 0; ist < Nst; ++ist){
      for(int in = 0; in < nin; ++in){
        v[nv1 + IDX2(nin, in, ist)] = w[nv2 + IDX2(nin, in, ist)];
      }
    }
  }

}

//====================================================================
