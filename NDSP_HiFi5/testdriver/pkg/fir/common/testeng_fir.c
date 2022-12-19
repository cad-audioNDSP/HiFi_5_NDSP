/* ------------------------------------------------------------------------ */
/* Copyright (c) 2021 by Cadence Design Systems, Inc. ALL RIGHTS RESERVED.  */
/* These coded instructions, statements, and computer programs ('Cadence    */
/* Libraries') are the copyrighted works of Cadence Design Systems Inc.     */
/* Cadence IP is licensed for use with Cadence processor cores only and     */
/* must not be used for any other processors and platforms. Your use of the */
/* Cadence Libraries is subject to the terms of the license agreement you   */
/* have entered into with Cadence Design Systems, or a sublicense granted   */
/* to you by a direct Cadence license.                                     */
/* ------------------------------------------------------------------------ */
/*  IntegrIT, Ltd.   www.integrIT.com, info@integrIT.com                    */
/*                                                                          */
/* NatureDSP_Baseband Library                                               */
/*                                                                          */
/* This library contains copyrighted materials, trade secrets and other     */
/* proprietary information of IntegrIT, Ltd. This software is licensed for  */
/* use with Cadence processor cores only and must not be used for any other */
/* processors and platforms. The license to use these sources was given to  */
/* Cadence, Inc. under Terms and Condition of a Software License Agreement  */
/* between Cadence, Inc. and IntegrIT, Ltd.                                 */
/* ------------------------------------------------------------------------ */
/*          Copyright (C) 2009-2021 IntegrIT, Limited.                      */
/*                      All Rights Reserved.                                */
/* ------------------------------------------------------------------------ */
/*
 * Test-engine add-on for older FIR API
 */

#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>

/* Cross-platform data type definitions. */
#include "types.h"
/* Fixed point arithmetics. */
#include "NatureDSP_Math.h"
/* Test environment configuration. */
#include "config.h"
/* DSP Library API. */
#include LIBRARY_HEADER(fir)
/* Test data vectors tools and SEQ-file reader. */
#include "vectools.h"
#include "testeng_fir.h"
/* Aligning memory allocator. */
#include "malloc16.h"

#define MIN(a,b)   ( (a)<(b) ? (a) : (b) )
#define MAX(a,b)   ( (a)>(b) ? (a) : (b) )

/* Allocate vectors and load the data set for autocorrelation functions:
 * vector X (in), vector Z (out) */
int te_loadFxn_autocorr( tTestEngContext * context )
{
  int  N;
  int res;
  ASSERT( context && context->seqFile );
  N = MAX( 0, context->args.dim[0] );
  if(N<0) { N=0; }
  memset( &context->dataSet, 0, sizeof(context->dataSet) );

  /* Allocate data vectors memory. */
  res = ( 4 == vecsAlloc( context->desc->isAligned, context->desc->fmt,
                         &context->dataSet.X  , N,
                         &context->dataSet.Z  , N,
                         &context->dataSet.Zlo, N,
                         &context->dataSet.Zhi, N, 0 ) );
  if ( res )
  {
    /* Load vectors data from the SEQ-file. */
    if ( !( res = seqFileReadVecs( context->seqFile,
                                  &context->dataSet.X,
                                  &context->dataSet.Zlo,
                                  &context->dataSet.Zhi, 0 ) ) )
    {
      printf( "te_loadFxn_autocorr(): failed to read vectors data\n");
    }
  }
  else
  {
    printf( "te_loadFxn_autocorr(): failed to allocate vectors\n");
  }

  /* Free vectors data if failed. */
  if ( !res )    te_freeVectors(context);
  return (res);

} /* te_loadFxn_autocorr() */

/* Allocate vectors and load the data set for crosscorrelation functions:
 * vector X (in), vector Z (out) */
int te_loadFxn_crosscorr( tTestEngContext * context )
{
  int M, N;
  int lenX=0,lenY=0,lenZ=0;
  int res,fmtX,fmtY,fmtZ;
  ASSERT( context && context->seqFile );
  memset( &context->dataSet, 0, sizeof(context->dataSet) );

  M = MAX( 0, context->args.dim[0] );
  N = MAX( 0, context->args.dim[1] );
  switch(context->desc->extraParam & TE_FIR_OTHER_API_MASK)
  {
      case TE_FIR_CROSSCORR_API:
      case TE_FIR_CONVOLVE_API:
            lenX=N;
            lenY=M;
            lenZ=N;
            break;
      case (TE_FIR_LINEAR_API|TE_FIR_CROSSCORR_API):
      case (TE_FIR_LINEAR_API|TE_FIR_CONVOLVE_API ):
            lenX=N;
            lenY=M;
            lenZ=N+M-1;
            break;
      default: /* should not happen */
            ASSERT(0);
            break;
  }

  fmtX=fmtY=fmtZ=context->desc->fmt;
  if ((context->desc->extraParam&TE_FIR_OTHER_TYPE_MASK) == TE_FIR_OTHER_32X16)
  {
      if (context->desc->fmt & FMT_CPLX)
      {
          fmtX=fmtZ=FMT_CPLX|FMT_FRACT32;
          fmtY=     FMT_CPLX|FMT_FRACT16;
      }
      else
      {
          fmtX=fmtZ=FMT_REAL|FMT_FRACT32;
          fmtY=     FMT_REAL|FMT_FRACT16;
      }
  }

  /* Allocate data vectors memory. */
  res = ( 1 == vecsAlloc( context->desc->isAligned, fmtX,
                         &context->dataSet.X  , lenX, 0 ) );
  res&= ( 1 == vecsAlloc( context->desc->isAligned, fmtY,
                         &context->dataSet.Y  , lenY, 0 ) );
  res&= ( 3 == vecsAlloc( context->desc->isAligned, fmtZ,
                         &context->dataSet.Z  , lenZ,
                         &context->dataSet.Zlo, lenZ,
                         &context->dataSet.Zhi, lenZ, 0 ) );
  if ( res )
  {
    /* Load vectors data from the SEQ-file. */
    if ( !( res = seqFileReadVecs( context->seqFile,
                                  &context->dataSet.X,
                                  &context->dataSet.Y,
                                  &context->dataSet.Zlo,
                                  &context->dataSet.Zhi, 0 ) ) )
    {
      printf( "te_loadFxn_crosscorr(): failed to read vectors data\n");
    }
  }
  else
  {
    printf( "te_loadFxn_crosscorr(): failed to allocate vectors\n");
  }

  /* Free vectors data if failed. */
  if ( !res )    te_freeVectors(context);
  return (res);

} /* te_loadFxn_crosscorr() */

/* Apply the target function to the test case data set (crosscorrelation):
 * vector X (in), vector Z (out) */
void te_processFxn_autocorr( tTestEngContext * context )
{
  typedef void tFxnautocorr_scratch( void* S, void * X, const void * Z,  int N );
  typedef void tFxnautocorr( void * X, const void * Z,  int N );
  void *X, *Z;
  int N;
  ASSERT( context && context->target.fut );
  X = vecGetElem( &context->dataSet.X, 0 );
  Z = vecGetElem( &context->dataSet.Z, 0 );
  N=context->args.dim[0];
  te_vReportStd(context);
  if(context->desc->isAligned)
  {
      ((tFxnautocorr*)context->target.fut)( Z, X, N );
  }
  else
  {
      /* allocate scratch and use another API */
      void * scratch;
      size_t sz=0;
      switch(context->desc->extraParam & (TE_FIR_OTHER_TYPE_MASK|TE_FIR_OTHER_API_MASK))
      {
          case TE_FIR_AUTOCORR_API:
              switch (context->desc->fmt & FMT_DTYPE_MASK)
              {
              case FMT_FRACT16:    sz=FIR_ACORRA16X16_SCRATCH_SIZE(N); break;
              case FMT_FRACT32:    sz=FIR_ACORRA32X32_SCRATCH_SIZE(N); break;
              case FMT_FLOAT32:    sz=FIR_ACORRAF_SCRATCH_SIZE(N); break;
              default: /* should not happen */
                  ASSERT(0);
              }
              break;
          case TE_FIR_LINEAR_API:
              switch (context->desc->fmt & FMT_DTYPE_MASK)
              {
              case FMT_FRACT16:    sz=FIR_LACORRA16X16_SCRATCH_SIZE(N); break;
              case FMT_FRACT32:    sz=FIR_LACORRA32X32_SCRATCH_SIZE(N); break;
              default: /* should not happen */
                  ASSERT(0);
              }
              break;
#if 0 // HiFi3/3z API
          case TE_FIR_AUTOCORR_API|TE_FIR_OTHER_24X24:
              sz=FIR_ACORRA24X24_SCRATCH_SIZE(N); break;
#endif
          case TE_FIR_AUTOCORR_API|TE_FIR_OTHER_EP:
              sz=FIR_ACORRA32X32EP_SCRATCH_SIZE(N); break;

          default: /* should not happen */
              ASSERT(0);
      }
      ASSERT(sz);
      scratch = mallocAlign(sz, 0);
      ((tFxnautocorr_scratch*)context->target.fut)( scratch, Z , X, N );
      freeAlign(scratch);
  }

} /* te_processFxn_autocorr() */ 


/* Apply the target function to the test case data set (crosscorrelation):
 * vector X (in), vector Z (out) */
void te_processFxn_crosscorr( tTestEngContext * context )
{
  typedef void tFxncrosscorr( void * Z, const void * X,  const void * Y,  int N , int M );
  typedef void tFxncrosscorr_scratch(void* S, void * Z, const void * X,  const void * Y,  int N , int M );
  void *X, *Y, *Z;
  int M,N;
  ASSERT( context && context->target.fut );
  X = vecGetElem( &context->dataSet.X, 0 );
  Y = vecGetElem( &context->dataSet.Y, 0 );
  Z = vecGetElem( &context->dataSet.Z, 0 );
  M = MAX( 0, context->args.dim[0] );
  N = MAX( 0, context->args.dim[1] );
  te_vReportStd(context);
  if(!context->desc->isAligned)
  {
      /* select right scratch size */
      size_t sz;
      sz=0;
      switch(context->desc->extraParam & (TE_FIR_OTHER_TYPE_MASK|TE_FIR_OTHER_API_MASK))
      {
          case TE_FIR_CROSSCORR_API:
              switch(context->desc->fmt)
              {
                  case FMT_REAL|FMT_FRACT16:    sz=FIR_XCORRA16X16_SCRATCH_SIZE(N,M); break;
                  case FMT_REAL|FMT_FRACT32:    sz=FIR_XCORRA32X32_SCRATCH_SIZE(N,M); break;
                  case FMT_REAL|FMT_FLOAT32:    sz=FIR_XCORRAF_SCRATCH_SIZE(N,M); break;
                  case FMT_CPLX|FMT_FRACT32:    sz=CXFIR_XCORRA32X32_SCRATCH_SIZE(N,M); break;
                  case FMT_CPLX|FMT_FLOAT32:    sz=CXFIR_XCORRAF_SCRATCH_SIZE(N,M); break;
                  default: ASSERT(0);
              }
              break;
          case (TE_FIR_CROSSCORR_API|TE_FIR_LINEAR_API):
              switch(context->desc->fmt)
              {
                  case FMT_REAL|FMT_FRACT16:    sz=FIR_LXCORRA16X16_SCRATCH_SIZE(N,M); break;
                  case FMT_REAL|FMT_FRACT32:    sz=FIR_LXCORRA32X32_SCRATCH_SIZE(N,M); break;
                  default: ASSERT(0);
              }
              break;
#if 0 // HiFi3/3z API
          case (TE_FIR_CROSSCORR_API|TE_FIR_OTHER_24X24):
              sz=FIR_XCORRA24X24_SCRATCH_SIZE(N,M); break;
#endif
          case (TE_FIR_CROSSCORR_API|TE_FIR_OTHER_32X16):
              sz=FIR_XCORRA32X16_SCRATCH_SIZE(N,M); break;
          case (TE_FIR_CROSSCORR_API|TE_FIR_OTHER_EP):
              sz=FIR_XCORRA32X32EP_SCRATCH_SIZE(N,M); break;
          case TE_FIR_CONVOLVE_API:
              switch(context->desc->fmt)
              {
                  case FMT_REAL|FMT_FRACT16:    sz=FIR_CONVOLA16X16_SCRATCH_SIZE(N,M); break;
                  case FMT_REAL|FMT_FRACT32:    sz=FIR_CONVOLA32X32_SCRATCH_SIZE(N,M); break;
                  case FMT_REAL|FMT_FLOAT32:    sz=FIR_CONVOLAF_SCRATCH_SIZE(N,M); break;
                  default: ASSERT(0);
              }
              break;
          case (TE_FIR_CONVOLVE_API|TE_FIR_LINEAR_API):
              switch(context->desc->fmt)
              {
              case FMT_REAL|FMT_FRACT16:    sz=FIR_LCONVOLA16X16_SCRATCH_SIZE(N,M); break;
              case FMT_REAL|FMT_FRACT32:    sz=FIR_LCONVOLA32X32_SCRATCH_SIZE(N,M); break;
              default: ASSERT(0);
              }
              break;
#if 0 // HiFi3/3z API
          case (TE_FIR_CONVOLVE_API|TE_FIR_OTHER_24X24):
              sz=FIR_CONVOLA24X24_SCRATCH_SIZE(N,M); break;
#endif
          case (TE_FIR_CONVOLVE_API|TE_FIR_OTHER_32X16):
              sz= (context->desc->fmt & FMT_CPLX) ?
              CXFIR_CONVOLA32X16_SCRATCH_SIZE(N,M) : FIR_CONVOLA32X16_SCRATCH_SIZE(N,M); break;
          case (TE_FIR_CONVOLVE_API|TE_FIR_OTHER_EP):
              sz=FIR_CONVOLA32X32EP_SCRATCH_SIZE(N,M); break;
          default:
              ASSERT(0);
      }
      if (0!=sz)
      {
          void * scratch;
          scratch = mallocAlign(sz,0);
          ((tFxncrosscorr_scratch*)context->target.fut)( scratch, Z, X, Y, N, M );
          freeAlign(scratch);
      }
      else
      {
          ((tFxncrosscorr_scratch*)context->target.fut)( NULL, Z, X, Y, N, M );
      }
  }
  else
  {
    ((tFxncrosscorr*)context->target.fut)( Z, X, Y, N, M );
  }
} /* te_processFxn_crosscorr() */ 

/* Allocate vectors and load the data set for GCC-PHAT functions:
 * vector X (in), vector Y (in), vector U (in), vector V (temp), vector Z (out) */
int te_loadFxn_gccphat( tTestEngContext * context )
{
    int N;
    int szXYUV, szZ;
    int fmtXYV, fmtU, fmtZ;
    int res=0;
    /* Retrieve test data dimensions and formats. */
    N = context->args.dim[0];
    szXYUV = MAX(0, N); szZ = 1;
    fmtXYV = context->desc->fmt;
    fmtU = (context->desc->fmt & FMT_DOMAIN_MASK) | FMT_FLOAT64;
    fmtZ = FMT_REAL|FMT_FLOAT64;
    /* Allocate vectors for input data, results and lower/upper result thresholds. */
    if (3!=vecsAlloc(context->desc->isAligned, fmtXYV,
                     &context->dataSet.X, szXYUV,
                     &context->dataSet.Y, szXYUV,
                     &context->dataSet.V, szXYUV, NULL) ||
        0==vecAlloc (&context->dataSet.U, szXYUV, TE_ALIGN_YES, fmtU, NULL) ||
        3!=vecsAlloc(TE_ALIGN_YES, fmtZ,
                     &context->dataSet.Z, szZ,
                     &context->dataSet.Zlo, szZ,
                     &context->dataSet.Zhi, szZ, NULL)) {
        printf("te_loadFxn_gccphat: failed to allocate vectors\n");
    /* Load input data and lower/upper result thresholds from the SEQ-file. */
    } else if (!seqFileReadVecs(context->seqFile, 
                                &context->dataSet.Y, &context->dataSet.X, &context->dataSet.U,
                                &context->dataSet.Zlo, &context->dataSet.Zhi, NULL)) {
        printf("te_loadFxn_gccphat: failed to load vectors from the SEQ-file\n");
    } else {
        *vecGetElem_fl64(&context->dataSet.Z, 0) = -INFINITY;
        res = 1;
    }
    return res;
} /* te_loadFxn_gccphat() */

/* Apply the target function to the test case data set (GCC-PHAT):
  * vector X (in), vector Y (in), vector U (in), vector V (temp), vector Z (out) */
void te_processFxn_gccphat( tTestEngContext * context )
{
    typedef int  tFxn_fxp(void * pScr, void * U, const void * X, const void * Y, int N);
    typedef void tFxn_flp(void * pScr, void * U, const void * X, const void * Y, int N);
    void *pX, *pY, *pV, *pScr;
    float64_t *pU, *pZ, *pVdp;
    float64_t p, q, sinad;
    tVec vecVdp;
    size_t szScr;
    int sft=0, logNup;
    int n, N;
    ASSERT( context && context->target.fut );
    /* Test data dimension */
    N = context->args.dim[0];
    /* Base-2 log of N rounded up to the next integer. */
    logNup = 31-S_exp0_l(N-1);
    /* Allocate vector for FUT results converted to DP FP. */
    if (!vecAlloc(&vecVdp, MAX(0, N), TE_ALIGN_YES, 
                  (context->desc->fmt & FMT_DOMAIN_MASK)|FMT_FLOAT64, NULL)) {
        printf("te_processFxn_gccphat: failed to allocate a temporary vector\n"); 
        return;
    }
    /* Retrieven pointers to in/out arrays. */
    pX = vecGetElem     ( &context->dataSet.X, 0 );
    pY = vecGetElem     ( &context->dataSet.Y, 0 );
    pU = vecGetElem_fl64( &context->dataSet.U, 0 );
    pZ = vecGetElem_fl64( &context->dataSet.Z, 0 );
    pV = vecGetElem     ( &context->dataSet.V, 0 );
    pVdp = vecGetElem_fl64(&vecVdp, 0);
    /* Allocate the scratch memory area, if required. */
    szScr = 0;
    if (context->desc->extraPtr) {
        const te_gccphat_api_t * api = (te_gccphat_api_t*)context->desc->extraPtr;
        szScr = api->getScratchSize(N);
    }
    pScr = szScr==0 ? NULL : mallocAlign(szScr, 0);
    /* Update the Test Coverage report. */
    te_vReportStd(context);
    /* Invoke the test target function */
    switch (context->desc->fmt & FMT_DTYPE_MASK) {
    case FMT_FRACT32: sft = logNup-((tFxn_fxp*)context->target.fut)(pScr, pV, pX, pY, N); break;
    case FMT_FLOAT32: ((tFxn_flp*)context->target.fut)(pScr, pV, pX, pY, N); break;
    }
    /* Convert the results to DP FP. */
    vecToFp64(pVdp, &context->dataSet.V, sft);
    /* Compare results with reference data from the SEQ-file and compute
     * the SINAD ratio: sinad = 10*log10(|reference|^2/|result-reference|^2)
     * Here |.| denotes the L2 norm. */
    NASSERT((context->desc->fmt & FMT_DOMAIN_MASK)==FMT_REAL); /* TBD for complex data! */
    p = q = 0;
    for ( n=0; n<N; n++ ) {
        float64_t e = pVdp[n]-pU[n];
        p += pU[n]*pU[n]; q += e*e;
    }
    pZ[0] = sinad = 10.*log10(p/(q+DBL_MIN));
#if 0
    printf(" N = %3d  SINAD = %.1f dB ", N, sinad);
#endif
    /* Cleanup */
    if (pScr) freeAlign(pScr);
    vecFree(&vecVdp);
} /* te_processFxn_gccphat() */

/* function reads impulse response from file and creates FIR structure. returns 0 if failed */
int te_create_lms(tTestEngContext * context)
{
    tLmsContext *lmsContext;
    lmsContext = (tLmsContext *)malloc(sizeof(*lmsContext));
    context->target.handle = (void*)lmsContext;
    if (lmsContext==NULL) return 0;
    memset(lmsContext, 0, sizeof(*lmsContext));
    return 1;
}

void lms_vecsFree(tTestEngContext * context)
{
    tLmsContext *lmsContext;
    lmsContext = (tLmsContext *)context->target.handle;
    if (lmsContext)
    {
        if(lmsContext->h.szBulk) vecFree(&lmsContext->h);
        if(lmsContext->normmu.szBulk) vecFree(&lmsContext->normmu);
    }
}
/* function destroys FIR structure, returns 0 if failed */
int te_destroy_lms(tTestEngContext * context)
{
    tLmsContext *lmsContext;
    lmsContext = (tLmsContext *)context->target.handle;
    if (lmsContext)
    {
        lms_vecsFree(context);
        free(lmsContext);
    }
    return 1;
}

/* Allocate vectors and load the data set for lms functions: */
int te_loadFxn_lms( tTestEngContext * context )
{
    int irfmt,errfmt,normmufmt;
    tLmsContext *lmsContext;
    lmsContext = (tLmsContext *)context->target.handle;
    int M, N, P;
    int res;
    ASSERT( context && context->seqFile );

    M = context->args.dim[0] ;
    N = context->args.dim[1] ;
    P = MAX( 0, M+N-1);
    M = MAX( 0, M );
    N = MAX( 0, N );
    memset( &context->dataSet, 0, sizeof(context->dataSet) );
    lms_vecsFree(context);
    /* use 32-bit IR data for 16x32-bit functions */
    irfmt    = ((context->desc->extraParam&TE_FIR_OTHER_TYPE_MASK) == TE_FIR_OTHER_32X16) ? (FMT_REAL|FMT_FRACT32) : context->desc->fmt;
    errfmt   = ((context->desc->extraParam&TE_FIR_OTHER_TYPE_MASK) == TE_FIR_OTHER_32X16) ? (FMT_REAL|FMT_FRACT32) : context->desc->fmt;
    normmufmt= ((context->desc->extraParam&TE_FIR_OTHER_TYPE_MASK) == TE_FIR_OTHER_32X16) ? (FMT_REAL|FMT_FRACT32) : context->desc->fmt & ~FMT_CPLX;

  /* Allocate data vectors memory. */
  res = ( 2 == vecsAlloc( context->desc->isAligned, context->desc->fmt,
                         &context->dataSet.X  , P,
                         &context->dataSet.Y  , N, /*reference (near end) data vector (r) */
                         0 )) ;
  res&= ( 3 == vecsAlloc( context->desc->isAligned, errfmt,
                         &context->dataSet.Z  , N, /*error vector                         */
                         &context->dataSet.Zlo, N, /*                                     */
                         &context->dataSet.Zhi, N, /*                                     */
                         0 )) ;
  res&= ( 1 == vecsAlloc( context->desc->isAligned, normmufmt,
                         &lmsContext->normmu,   2,
                         0 )) ;
  res&= ( 4 == vecsAlloc( context->desc->isAligned, irfmt,
                         &context->dataSet.W  , M, /* impulse response (h)                */
                         &context->dataSet.Wlo, M,
                         &context->dataSet.Whi, M,
                         &lmsContext->h,        M,
                         0 )) ;

  if ( res )
  {
    /* Load vectors data from the SEQ-file. */
    if ( !( res = seqFileReadVecs( context->seqFile,
                                  &lmsContext->h,
                                  &context->dataSet.Y,
                                  &context->dataSet.X,
                                  &lmsContext->normmu,
                                  &context->dataSet.Zlo,
                                  &context->dataSet.Zhi, 
                                  &context->dataSet.Wlo,
                                  &context->dataSet.Whi, 
                                  0 ) ) )
    {
      printf( "te_loadFxn_lms(): failed to read vectors data\n");
    }
  }
  else
  {
    printf( "te_loadFxn_lms(): failed to allocate vectors\n");
  }

  /* Free vectors data if failed. */
  if ( !res )    
  {
    te_freeVectors(context);
    lms_vecsFree(context);
  }
  return (res);

} /* te_loadFxn_lms() */

/* Apply the target function to the test case data set (LMS):
 * vector X (in), vector Z (out) */
void te_processFxn_lms( tTestEngContext * context )
{
    tLmsContext *lmsContext;
    typedef void tFxn_complexf(complex_float *e, complex_float *h, const complex_float *r, const complex_float *x, float32_t norm, float32_t mu, int N, int M);
    typedef void tFxn_float32x32(float32_t *e, float32_t *h, const float32_t *r, const float32_t *x, float32_t norm, float32_t mu, int N, int M);
    typedef void tFxn_fract32x32(fract32   *e, fract32   *h, const fract32   *r, const fract32   *x, fract32   norm, fract32   mu, int N, int M);
    typedef void tFxn_fract16x32(fract32   *e, fract32   *h, const fract16   *r, const fract16   *x, fract32   norm, fract16   mu, int N, int M);
    typedef void tFxn_fract16x16(fract16   *e, fract16   *h, const fract16   *r, const fract16   *x, fract16   norm, fract16   mu, int N, int M);
    int M, N;
    void *E, *Hout,*Hin,*R,*X;

    ASSERT( context && context->target.fut );
    lmsContext = (tLmsContext *)context->target.handle;
    N=context->args.dim[1]; M=context->args.dim[0];
    E    = vecGetElem( &context->dataSet.Z, 0 );
    Hout = vecGetElem( &context->dataSet.W, 0 );
    Hin  = vecGetElem( &lmsContext->h     , 0 );
    R    = vecGetElem( &context->dataSet.Y, 0 );
    X    = vecGetElem( &context->dataSet.X, 0 );
    memcpy(Hout,Hin,context->dataSet.W.nElem*context->dataSet.W.szElem);  /* copy input IR to the output */
    te_vReportStd(context);
    switch(context->desc->fmt & FMT_DTYPE_MASK)
    {
    case FMT_FLOAT32|FMT_CPLX:
    {
        const float32_t* _normmu=vecGetElem_fl32(&lmsContext->normmu,0);
        ((tFxn_complexf*)context->target.fut)((complex_float *)E,(complex_float *)Hout,(const complex_float *)R,(const complex_float *)X,_normmu[0],_normmu[1],N,M);
        break;
    }
    case FMT_FLOAT32:
    {
        const float32_t* _normmu=vecGetElem_fl32(&lmsContext->normmu,0);
        ((tFxn_float32x32*)context->target.fut)((float32_t *)E,(float32_t *)Hout,(const float32_t *)R,(const float32_t *)X,_normmu[0],_normmu[1],N,M);
        break;
    }
    case FMT_FRACT16:
    {
        if ((context->desc->extraParam&TE_FIR_OTHER_TYPE_MASK) == TE_FIR_OTHER_32X16)
        {
            const fract32  * _normmu=vecGetElem_fr32(&lmsContext->normmu,0);
            ((tFxn_fract16x32*)context->target.fut)((fract32 *)E,(fract32 *)Hout,(const fract16 *)R,(const fract16 *)X,_normmu[0],(int16_t)_normmu[1],N,M);
        }
        else
        {
            const fract16  * _normmu=vecGetElem_fr16(&lmsContext->normmu,0);
            ((tFxn_fract16x16*)context->target.fut)((fract16 *)E,(fract16 *)Hout,(const fract16 *)R,(const fract16 *)X,_normmu[0],(int16_t)_normmu[1],N,M);
        }
        break;
    }
    case FMT_FRACT32:
    {
        const fract32  * _normmu=vecGetElem_fr32(&lmsContext->normmu,0);
        ((tFxn_fract32x32*)context->target.fut)((fract32   *)E,(fract32   *)Hout,(const fract32   *)R,(const fract32   *)X,_normmu[0],_normmu[1],N,M);
        break;
    }
    default: 
        ASSERT(0); /* not supported yet */
        break;
    }
} /* te_processFxn_lms() */



/* Allocate vectors and load the data set for LMS convergence test: */
int te_loadFxn_lmsconv( tTestEngContext * context )
{
    int irfmt,errfmt,normmufmt;
    //tLmsContext *lmsContext;
    //lmsContext = (tLmsContext *)context->target.handle;
    int M, N, P,nAlloc;
    int res;
    ASSERT( context && context->seqFile );

    M = context->args.dim[0] ;
    N = context->args.dim[1] ;
    P = context->args.dim[2] ;
    memset( &context->dataSet, 0, sizeof(context->dataSet) );
    lms_vecsFree(context);
    /* use 32-bit IR data for 16x32-bit functions */
    irfmt    = ((context->desc->extraParam&TE_FIR_OTHER_TYPE_MASK) == TE_FIR_OTHER_32X16) ? (FMT_REAL|FMT_FRACT32) : context->desc->fmt;
    errfmt   = ((context->desc->extraParam&TE_FIR_OTHER_TYPE_MASK) == TE_FIR_OTHER_32X16) ? (FMT_REAL|FMT_FRACT32) : context->desc->fmt;
    normmufmt= ((context->desc->extraParam&TE_FIR_OTHER_TYPE_MASK) == TE_FIR_OTHER_32X16) ? (FMT_REAL|FMT_FRACT32) : context->desc->fmt & ~FMT_CPLX;

  /* Allocate data vectors memory. */
    nAlloc=0;
    nAlloc+=vecsAlloc( context->desc->isAligned, context->desc->fmt,
                      &context->dataSet.X  , M+N*P-1,
                      &context->dataSet.Y  , N*P, /*reference (near end) data vector (r) */
                      0 );
    nAlloc+=vecsAlloc( context->desc->isAligned, errfmt,
                      &context->dataSet.U  , N, /*error vector                         */
                      0 ) ;
    nAlloc+=vecsAlloc( context->desc->isAligned, normmufmt,
                       &context->dataSet.V,P,
                       &context->dataSet.auxVec[0],1,
                       0 ) ;
    nAlloc+=vecsAlloc( context->desc->isAligned, irfmt,
                       &context->dataSet.Z  , M, /* impulse response (h)                */
                       &context->dataSet.Zlo, M,
                       &context->dataSet.Zhi, M,
                       0 );
    res=nAlloc==8;
    if ( res )
    {
    /* Load vectors data from the SEQ-file. */
    if ( !( res = seqFileReadVecs( context->seqFile,
                                    &context->dataSet.auxVec[0],
                                    &context->dataSet.X,
                                    &context->dataSet.Y,
                                    &context->dataSet.V,
                                    &context->dataSet.Zlo,
                                    &context->dataSet.Zhi,
                                    0 ) ) )
    {
      printf( "te_loadFxn_lmsconv(): failed to read vectors data\n");
    }
  }
  else
  {
    printf( "te_loadFxn_lmsconv(): failed to allocate vectors\n");
  }

  /* Free vectors data if failed. */
  if ( !res )    
  {
    te_freeVectors(context);
  }
  return (res);
}

/* Apply the target function to the test case data set (LMS convergence test):
*/
void te_processFxn_lmsconv( tTestEngContext * context )
{
    typedef void tFxn_float32 (void *e, void *h, const void *r, const void *x, void *norm, float32_t mu, int N, int M,int P);
    typedef void tFxn_fract32 (void *e, void *h, const void *r, const void *x, void *norm, fract32   mu, int N, int M,int P);
    typedef void tFxn_fract16 (void *e, void *h, const void *r, const void *x, void *norm, fract16   mu, int N, int M,int P);
    int M, N,P;
    void *E, *Hout,*R,*X,*mu,*norm;

    ASSERT( context && context->target.fut );
    P=context->args.dim[2]; N=context->args.dim[1]; M=context->args.dim[0];
    E    = vecGetElem( &context->dataSet.U, 0 );
    Hout = vecGetElem( &context->dataSet.Z, 0 );
    R    = vecGetElem( &context->dataSet.Y, 0 );
    X    = vecGetElem( &context->dataSet.X, 0 );
    mu   = vecGetElem( &context->dataSet.auxVec[0], 0 );
    norm = vecGetElem( &context->dataSet.V, 0 );
        //float32_t* hlo= (float32_t* )vecGetElem( &context->dataSet.Zlo, 0 );
        //float32_t* hhi= (float32_t* )vecGetElem( &context->dataSet.Zhi, 0 );
    te_vReportStd(context);
    switch(context->desc->fmt & FMT_DTYPE_MASK)
    {
    case FMT_FLOAT32|FMT_CPLX:
    case FMT_FLOAT32:
    {
        ((tFxn_float32*)context->target.fut)(E,Hout,R,X,norm,*(float32_t*)mu,N,M,P);
        break;
    }
    case FMT_FRACT16:
    {
        if ((context->desc->extraParam&TE_FIR_OTHER_TYPE_MASK) == TE_FIR_OTHER_32X16)
        {
         ((tFxn_fract32*)context->target.fut)(E,Hout,R,X,norm,*(fract32*)mu,N,M,P);
        }
        else
        {
         ((tFxn_fract16*)context->target.fut)(E,Hout,R,X,norm,*(fract16*)mu,N,M,P);
        }
        break;
    }
    case FMT_FRACT32:
    {
        ((tFxn_fract32*)context->target.fut)(E,Hout,R,X,norm,*(fract32*)mu,N,M,P);
        break;
    }
    default: 
        ASSERT(0); /* not supported yet */
        break;
    }
} /* te_processFxn_lmsconv() */

