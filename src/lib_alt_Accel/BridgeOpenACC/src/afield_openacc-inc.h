/*!
      @file    afield_openacc-inc.h
      @brief
      @author  Hideo Matsufuru (matufuru)
               $LastChangedBy: matufuru $
      @date    $LastChangedDate:: 2025-08-20 16:24:03 #$
      @version $LastChangedRevision: 2653 $
*/

//#include "inline/define_params.h"
//#include "inline/define_index.h"

//====================================================================
void afield_init(real_t* RESTRICT data, const int size)
{
#pragma acc enter data create(data[0:size])

#pragma acc data present(data[0:size]) copyin(size)
#pragma acc parallel num_workers(NUM_WORKERS) vector_length(VECTOR_LENGTH)
  {
#pragma acc loop gang worker vector
    for(int i = 0; i < size; ++i){
      data[i] = 0.0;
    }
  }

}

//====================================================================
void afield_tidyup(real_t* RESTRICT data, const int size)
{
#pragma acc exit data delete(data[0:size])
}

//====================================================================
void afield_set(real_t* RESTRICT v, real_t a, const int nin, const int nst)
{
  int Nst_pad = CEIL_NWP(nst);
  int nv = nin * Nst_pad;
#pragma acc data present(v[0:nv]) copyin(a, nin, nst)
#pragma acc parallel num_workers(NUM_WORKERS) vector_length(VECTOR_LENGTH)
  {
#pragma acc loop gang worker vector
    for(int ist = 0; ist < nst; ++ist){

      for(int in = 0; in < nin; ++in){
        v[IDX2(nin, in, ist)] = a;
      }

    }
  }

}

//====================================================================
void copy_to_device(real_t* RESTRICT v, int nv)
{
#pragma acc update device(v[0:nv])
}

//====================================================================
void copy_to_device(real_t* RESTRICT v, int nv1, int nv)
{
#pragma acc update device(v[nv1:nv])
}

//====================================================================
void copy_from_device(real_t* RESTRICT v, int nv)
{
#pragma acc update host(v[0:nv])
}

//====================================================================
void copy_from_device(real_t* RESTRICT v, int nv1, int nv)
{

#pragma acc update host(v[nv1:nv])
}

//====================================================================
void convert(real_t* RESTRICT v, double* RESTRICT w,
                 int nin, int nvol, int nvol_pad)
{
  // w is assumed to be Field, thus its volume is not ceiled.
  int nv1 = nin * nvol;
  int nv2 = nin * nvol_pad;

#pragma acc data present(v[0:nv2]) copyin(w[0:nv1], nin, nvol, nvol_pad)
#pragma acc parallel num_workers(NUM_WORKERS) vector_length(VECTOR_LENGTH)
  {
#pragma acc loop gang worker vector
    for(int ist = 0; ist < nvol_pad; ++ist){
      if(ist < nvol){
        for(int in = 0; in < nin; ++in){
          v[IDX2(nin, in, ist)] = w[IDX_CORE(nin, in, ist)];
        }
      }else{
        for(int in = 0; in < nin; ++in){
          v[IDX2(nin, in, ist)] = 0.0;
        }
      }
    }
  }

}

//====================================================================
void reverse(double* RESTRICT v, real_t* RESTRICT w,
                 int nin, int nvol, int nvol_pad)
{
  // v is assumed to be Field, thus its volume is not ceiled.
  int nv1 = nin * nvol;
  int nv2 = nin * nvol_pad;

#pragma acc data present(w[0:nv2]) copyin(nin, nvol, nvol_pad) \
                 copyout(v[0:nv1])
#pragma acc parallel num_workers(NUM_WORKERS) vector_length(VECTOR_LENGTH)
  {
#pragma acc loop gang worker vector
    for(int ist = 0; ist < nvol; ++ist){
      for(int in = 0; in < nin; ++in){
        v[IDX_CORE(nin, in, ist)] = w[IDX2(nin, in, ist)];
      }
    }
  }

}

//====================================================================
void copy(real_t* RESTRICT v, real_t* RESTRICT w, int nin, int nvol)
{
  int Nvol_pad = CEIL_NWP(nvol);
  int nv = nin * Nvol_pad;

#pragma acc data present(v[0:nv], w[0:nv]) copyin(nin, nvol)
#pragma acc parallel num_workers(NUM_WORKERS) vector_length(VECTOR_LENGTH)
  {
#pragma acc loop gang worker vector
    for(int ist = 0; ist < nvol; ++ist){
      for(int in = 0; in < nin; ++in){
        v[IDX2(nin, in, ist)] = w[IDX2(nin, in, ist)];
      }
    }
  }

}

//====================================================================
void copy(real_t* RESTRICT v, int nv1, real_t* RESTRICT w, int nv2,
              int nin, int nvol)
{
  int Nvol_pad = CEIL_NWP(nvol);
  int nv = nin * Nvol_pad;

#pragma acc data present(v[nv1:nv], w[nv2:nv]) \
                 copyin(nin, nvol, nv1, nv2)

#pragma acc parallel num_workers(NUM_WORKERS) vector_length(VECTOR_LENGTH)
  {
#pragma acc loop gang worker vector
    for(int ist = 0; ist < nvol; ++ist){
      for(int in = 0; in < nin; ++in){
        v[nv1 + IDX2(nin, in, ist)] = w[nv2 + IDX2(nin, in, ist)];
      }
    }
  }

}

//====================================================================
void axpy(real_t *RESTRICT v, int nv1, real_t a,
              real_t *RESTRICT w, int nv2, int nin, int nvol)
{
  int Nvol_pad = CEIL_NWP(nvol);
  int nv = nin * Nvol_pad;

#pragma acc data present(v[nv1:nv],w[nv2:nv]) \
                 copyin(a, nin, nvol, nv1, nv2)
#pragma acc parallel num_workers(NUM_WORKERS) vector_length(VECTOR_LENGTH)
  {
#pragma acc loop gang worker vector
    for(int ist = 0; ist < nvol; ++ist){
      for(int in = 0; in < nin; ++in){
        v[nv1 + IDX2(nin, in, ist)] += a * w[nv2 + IDX2(nin, in, ist)];
      }
    }
  }

}

//====================================================================
void axpy(real_t *RESTRICT v, int nv1, real_t ar, real_t ai,
          real_t *RESTRICT w, int nv2, int nin, int nvol)
{
  int Nvol_pad = CEIL_NWP(nvol);
  int nv = nin * Nvol_pad;

#pragma acc data present(v[nv1:nv],w[nv2:nv]) \
                 copyin(ar, ai, nin, nvol, nv1, nv2)
#pragma acc parallel num_workers(NUM_WORKERS) vector_length(VECTOR_LENGTH)
  {
    int nin2 = nin/2;

#pragma acc loop gang worker vector
    for(int ist = 0; ist < nvol; ++ist){
      for(int in = 0; in < nin2; ++in){
        real_t wr = w[nv2 + IDX2(nin, 2*in,   ist)];
        real_t wi = w[nv2 + IDX2(nin, 2*in+1, ist)];
        v[nv1 + IDX2(nin, 2*in,   ist)] += ar * wr - ai * wi;
        v[nv1 + IDX2(nin, 2*in+1, ist)] += ai * wr + ar * wi;
      }
    }
  }

}

//====================================================================
void aypx(real_t a, real_t *RESTRICT v, int nv1,
          real_t *RESTRICT w, int nv2, int nin, int nvol)
{
  int Nvol_pad = CEIL_NWP(nvol);
  int nv = nin * Nvol_pad;

#pragma acc data present(v[nv1:nv],w[nv2:nv])
#pragma acc parallel num_workers(NUM_WORKERS) vector_length(VECTOR_LENGTH)
  {
#pragma acc loop gang worker vector
    for(int ist = 0; ist < nvol; ++ist){
      for(int in = 0; in < nin; ++in){
        v[nv1 + IDX2(nin, in, ist)] = a * v[nv1 + IDX2(nin, in, ist)]
                                        + w[nv2 + IDX2(nin, in, ist)];
      }
    }
  }

}

//====================================================================
void aypx(real_t ar, real_t ai, real_t *RESTRICT v, int nv1,
          real_t *RESTRICT w, int nv2, int nin, int nvol)
{
  int Nvol_pad = CEIL_NWP(nvol);
  int nv = nin * Nvol_pad;

#pragma acc data present(v[nv1:nv],w[nv2:nv]) \
                 copyin(ar, ai, nin, nvol)
#pragma acc parallel num_workers(NUM_WORKERS) vector_length(VECTOR_LENGTH)
  {
    int nin2 = nin/2;
#pragma acc loop gang worker vector
    for(int ist = 0; ist < nvol; ++ist){
      for(int in = 0; in < nin2; ++in){
        real_t vr = v[nv1 + IDX2(nin, 2*in,   ist)];
        real_t vi = v[nv1 + IDX2(nin, 2*in+1, ist)];
        v[nv1 + IDX2(nin, 2*in,   ist)] = ar * vr - ai * vi
                                  + w[nv2 + IDX2(nin, 2*in,   ist)];
        v[nv1 + IDX2(nin, 2*in+1, ist)] = ai * vr + ar * vi
                                  + w[nv2 + IDX2(nin, 2*in+1, ist)];
      }
    }

  }

}

//====================================================================
void scal(real_t* RESTRICT v, int nv1, real_t a, int nin, int nvol)
{
  int Nvol_pad = CEIL_NWP(nvol);
  int nv = nin * Nvol_pad;

#pragma acc data present(v[nv1:nv]) copyin(a, nin, nvol)
#pragma acc parallel num_workers(NUM_WORKERS) vector_length(VECTOR_LENGTH)
  {
#pragma acc loop gang worker vector
    for(int ist = 0; ist < nvol; ++ist){
      for(int in = 0; in < nin; ++in){
        v[nv1 + IDX2(nin, in, ist)] *= a;
      }
    }
  }

}

//====================================================================
void scal(real_t* RESTRICT v, int nv1, real_t ar, real_t ai,
              int nin, int nvol)
{
  int Nvol_pad = CEIL_NWP(nvol);
  int nv = nin * Nvol_pad;
  int nin2 = nin/2;

#pragma acc data present(v[nv1:nv]) copyin(ar, ai, nin, nvol)
#pragma acc parallel num_workers(NUM_WORKERS) vector_length(VECTOR_LENGTH)
  {
#pragma acc loop gang worker vector
    for(int ist = 0; ist < nvol; ++ist){
      for(int in = 0; in < nin2; ++in){
        real_t vr = v[nv1 + IDX2(nin, 2*in,   ist)];
        real_t vi = v[nv1 + IDX2(nin, 2*in+1, ist)];
        v[nv1 + IDX2(nin, 2*in,   ist)] = ar * vr - ai * vi;
        v[nv1 + IDX2(nin, 2*in+1, ist)] = ai * vr + ar * vi;
      }
    }

  }

}

//====================================================================
real_t dot(real_t* RESTRICT v1, real_t* RESTRICT v2,
               int nin, int nvol)
{
  int Nvol_pad = CEIL_NWP(nvol);
  int nv = nin * Nvol_pad;

  real_t a = 0.0;

#pragma acc data present(v1[0:nv], v2[0:nv]) copyin(nin, nvol)
#pragma acc parallel num_workers(NUM_WORKERS) vector_length(VECTOR_LENGTH)
  {
#pragma acc loop gang worker vector reduction(+:a)
    for(int ist = 0; ist < nvol; ++ist){
      for(int in = 0; in < nin; ++in){
        a += v1[IDX2(nin, in, ist)] * v2[IDX2(nin, in, ist)];
      }
    }
  }

  return a;

}

//====================================================================
real_t norm2(real_t* RESTRICT v1, int nin, int nvol){

  int Nvol_pad = CEIL_NWP(nvol);
  int nv = nin * Nvol_pad;

  double a = 0.0;

#pragma acc data present(v1[0:nv])
#pragma acc parallel num_workers(NUM_WORKERS) vector_length(VECTOR_LENGTH)
  {
#pragma acc loop gang worker vector reduction(+:a)
    for(int ist = 0; ist < nvol; ++ist){
      for(int in = 0; in < nin; ++in){
        a += double(v1[IDX2(nin, in, ist)] * v1[IDX2(nin, in, ist)]);
      }
    }
  }

  return real_t(a);

}

//====================================================================
void dotc(real_t* ar, real_t* ai,
          real_t *RESTRICT v1, real_t *RESTRICT v2,
          int nin, int nst)
{

  double ar2 = 0.0;
  double ai2 = 0.0;

  int Nst_pad = CEIL_NWP(nst);
  int nv = nin * Nst_pad;

#pragma acc data present(v1[0:nv], v2[0:nv]) copyin(nst, nin)
#pragma acc parallel num_workers(NUM_WORKERS) vector_length(VECTOR_LENGTH)
  {
    int nin2 = nin/2;

#pragma acc loop gang worker vector reduction(+:ar2,ai2)
    for(int ist = 0; ist < nst; ++ist){

    double ar1 = 0.0;
    double ai1 = 0.0;
    for(int in = 0; in < nin2; ++in){
      real_t v1r = v1[IDX2(nin, 2*in,   ist)];
      real_t v1i = v1[IDX2(nin, 2*in+1, ist)];
      real_t v2r = v2[IDX2(nin, 2*in,   ist)];
      real_t v2i = v2[IDX2(nin, 2*in+1, ist)];
      ar1 += double(v1r * v2r + v1i * v2i);
      ai1 += double(v1r * v2i - v1i * v2r);
    }
    ar2 += ar1;
    ai2 += ai1;
  }

  } // acc parallel

  *ar = real_t(ar2);
  *ai = real_t(ai2);

}

//====================================================================
void xI(real_t *RESTRICT v1, int nin, int nst){

  int Nst_pad = CEIL_NWP(nst);
  int nv = nin * Nst_pad;

#pragma acc data present(v1[0:nv]) copyin(nin, nst)
#pragma acc parallel num_workers(NUM_WORKERS) vector_length(VECTOR_LENGTH)
 {
    int nin2 = nin/2;

#pragma acc loop gang worker vector
  for(int ist = 0; ist < nst; ++ist){
    for(int in = 0; in < nin2; ++in){
      real_t vr = v1[IDX2(nin, 2*in,   ist)];
      real_t vi = v1[IDX2(nin, 2*in+1, ist)];
      v1[IDX2(nin, 2*in,   ist)] = -vi;
      v1[IDX2(nin, 2*in+1, ist)] =  vr;
    }
  }

 }

}

//====================================================================
void conjg(real_t *RESTRICT v1, int nin, int nst){

  int Nst_pad = CEIL_NWP(nst);
  int nv = nin * Nst_pad;

#pragma acc data present(v1[0:nv]) copyin(nin, nst)
#pragma acc parallel num_workers(NUM_WORKERS) vector_length(VECTOR_LENGTH)
 {
    int nin2 = nin/2;

#pragma acc loop gang worker vector
  for(int ist = 0; ist < nst; ++ist){
    for(int in = 0; in < nin2; ++in){
      real_t vi = v1[IDX2(nin, 2*in+1, ist)];
      v1[IDX2(nin, 2*in+1, ist)] = -vi;
    }
  }

 }

}

//====================================================================
