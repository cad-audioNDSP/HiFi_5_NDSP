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
* Test-engine add-on for vector mathematics
*/

#include <string.h>
#include <math.h>

/* Cross-platform data type definitions. */
#include "types.h"
/* Test engine API. */
#include "testeng_math.h"
#define MIN(a,b)   ( (a)<(b) ? (a) : (b) )
#define MAX(a,b)   ( (a)>(b) ? (a) : (b) )


/* processing function for vec_recip16x16 */
void te_math_processFxn_recip16x16( tTestEngContext * context )
{
  typedef void tFxn( int16_t *  frac, 
                  int16_t *exp, 
                  const int16_t *  x, 
                  int N);
  tVec frac,exponent;
  tFxn *fxn;
  const int16_t *X; 
  float32_t *Z;
  int16_t* pFrac,*pExp;
  int n,N;
  ASSERT( context && context->target.fut );
  X = vecGetElem_fr16( &context->dataSet.X, 0 );
  Z = vecGetElem_fl32( &context->dataSet.Z, 0 );
  N   = context->args.dim[0];
  vecsAlloc(context->desc->isAligned,FMT_FRACT16,&frac,N<0?0:N,&exponent,N<0?0:N,NULL);
  pFrac=vecGetElem_fr16(&frac,0);
  pExp =vecGetElem_fr16(&exponent,0);
  fxn = (tFxn*)context->target.fut;
  te_vReportStd(context);
  /* 
    test for overlapping X and pFrac for each 1st of 4 tests if overlapping 
    is enabled in extraParams
  */
  if (((context->args.caseNum & 3)==1) && (context->desc->extraParam & OVLP_XZ))
  {
      size_t szX,szfrac;
      szX=vecGetSize(&context->dataSet.X);
      szfrac=vecGetSize(&frac);
      (void)szfrac;
      ASSERT(szX==szfrac);
      memcpy(pFrac,X,szX);
      X=pFrac;
  }
  /* 
    test for overlapping X and pExp for each 2nd of 4 tests if overlapping 
    is enabled in extraParams
  */
  if (((context->args.caseNum & 3)==2) && (context->desc->extraParam & OVLP_XW))
  {
      size_t szX,szexp;
      szX=vecGetSize(&context->dataSet.X);
      szexp=vecGetSize(&exponent);
      (void)szexp;
      ASSERT(szX==szexp);
      memcpy(pExp,X,szX);
      X=pExp;
  }
  fxn( pFrac, pExp, X, N );
  /* convert results to single precision floating point */
  for(n=0; n<N; n++)
  {
      Z[n]=STDLIB_MATH(ldexpf)(pFrac[n],pExp[n]-15);
  }
  vecsFree(&frac,&exponent,NULL);
} /* te_math_processFxn_recip16x16() */

/* processing function for vec_divide16x16 */
void te_math_processFxn_divide16x16( tTestEngContext * context )
{
  typedef void tFxn( int16_t *  frac, 
                  int16_t *exp, 
                  const int16_t *  x, 
                  const int16_t *  y, 
                  int N);
  tVec frac,exponent;
  tFxn *fxn;
  const int16_t *X; 
  const int16_t *Y; 
  float32_t *Z;
  int16_t* pFrac,*pExp;
  int n,N;
  ASSERT( context && context->target.fut );
  X = vecGetElem_fr16( &context->dataSet.X, 0 );
  Y = vecGetElem_fr16( &context->dataSet.Y, 0 );
  Z = vecGetElem_fl32( &context->dataSet.Z, 0 );
  N   = context->args.dim[0];
  vecsAlloc(context->desc->isAligned,FMT_FRACT16,&frac,N<0?0:N,&exponent,N<0?0:N,NULL);
  pFrac=vecGetElem_fr16(&frac,0);
  pExp =vecGetElem_fr16(&exponent,0);
  fxn = (tFxn*)context->target.fut;
  /* 
    test for overlapping X and pFrac for each 1st of 8 tests if overlapping 
    is enabled in extraParams
  */
  if (((context->args.caseNum & 7)==1) && (context->desc->extraParam & OVLP_XZ))
  {
      size_t szX,szfrac;
      szX=vecGetSize(&context->dataSet.X);
      szfrac=vecGetSize(&frac);
      (void)szfrac;
      ASSERT(szX==szfrac);
      memcpy(pFrac,X,szX);
      X=pFrac;
  }
  /* 
    test for overlapping X and pExp for each 2nd of 8 tests if overlapping 
    is enabled in extraParams
  */
  if (((context->args.caseNum & 3)==2) && (context->desc->extraParam & OVLP_XW))
  {
      size_t szX,szexp;
      szX=vecGetSize(&context->dataSet.X);
      szexp=vecGetSize(&exponent);
      (void)szexp;
      ASSERT(szX==szexp);
      memcpy(pExp,X,szX);
      X=pExp;
  }
  /* 
    test for overlapping X and pFrac for each 3rd of 8 tests if overlapping 
    is enabled in extraParams
  */
  if (((context->args.caseNum & 7)==3) && (context->desc->extraParam & OVLP_XZ))
  {
      size_t szY,szfrac;
      szY=vecGetSize(&context->dataSet.Y);
      szfrac=vecGetSize(&frac);
      (void)szfrac;
      ASSERT(szY==szfrac);
      memcpy(pFrac,Y,szY);
      Y=pFrac;
  }
  /* 
    test for overlapping X and pExp for each 4th of 8 tests if overlapping 
    is enabled in extraParams
  */
  if (((context->args.caseNum & 7)==4) && (context->desc->extraParam & OVLP_XW))
  {
      size_t szY,szexp;
      szY=vecGetSize(&context->dataSet.Y);
      szexp=vecGetSize(&exponent);
      (void)szexp;
      ASSERT(szY==szexp);
      memcpy(pExp,Y,szY);
      Y=pExp;
  }
  te_vReportStd(context);
  fxn( pFrac, pExp, X,Y, N );
  /* convert results to single precision floating point */
  for(n=0; n<N; n++)
  {
      Z[n]=STDLIB_MATH(ldexpf)(pFrac[n],pExp[n]-15);
  }
  vecsFree(&frac,&exponent,NULL);
} /* te_math_processFxn_divide16x16() */

/* processing function for vec_divide32x32 */
void te_math_processFxn_divide32x32( tTestEngContext * context )
{
  typedef void tFxn( int32_t *  frac, 
                  int16_t *exp, 
                  const int32_t *  x, 
                  const int32_t *  y, 
                  int N);
  tVec frac,exponent;
  tFxn *fxn;
  const int32_t *X; 
  const int32_t *Y; 
  float64_t *Z;
  int32_t* pFrac;
  int16_t *pExp;
  int n,N;
  ASSERT( context && context->target.fut );
  X = vecGetElem_fr32( &context->dataSet.X, 0 );
  Y = vecGetElem_fr32( &context->dataSet.Y, 0 );
  Z = vecGetElem_fl64( &context->dataSet.Z, 0 );
  N   = context->args.dim[0];
  vecsAlloc(context->desc->isAligned,FMT_FRACT32,&frac,N<0?0:N,NULL);
  vecsAlloc(context->desc->isAligned,FMT_FRACT16,&exponent,N<0?0:N,NULL);
  pFrac=vecGetElem_fr32(&frac,0);
  pExp =vecGetElem_fr16(&exponent,0);
  fxn = (tFxn*)context->target.fut;
  /* 
    test for overlapping X and pFrac for each 1st of 8 tests if overlapping 
    is enabled in extraParams
  */
  if (((context->args.caseNum & 7)==1) && (context->desc->extraParam & OVLP_XZ))
  {
      size_t szX,szfrac;
      szX=vecGetSize(&context->dataSet.X);
      szfrac=vecGetSize(&frac);
      (void)szfrac;
      ASSERT(szX==szfrac);
      memcpy(pFrac,X,szX);
      X=pFrac;
  }
  /* 
    test for overlapping X and pFrac for each 3rd of 8 tests if overlapping 
    is enabled in extraParams
  */
  if (((context->args.caseNum & 7)==3) && (context->desc->extraParam & OVLP_XZ))
  {
      size_t szY,szfrac;
      szY=vecGetSize(&context->dataSet.Y);
      szfrac=vecGetSize(&frac);
      (void)szfrac;
      ASSERT(szY==szfrac);
      memcpy(pFrac,Y,szY);
      Y=pFrac;
  }

  te_vReportStd(context);
  fxn( pFrac, pExp, X,Y, N );
  /* convert results to single precision floating point */
  for(n=0; n<N; n++)
  {
      Z[n]=STDLIB_MATH(ldexp)((float64_t)pFrac[n],pExp[n]-31);
  }
  vecsFree(&frac,&exponent,NULL);
} /* te_math_processFxn_divide32x32() */

/* processing function for vec_recip32x32 vec_recip24x24 */
void te_math_processFxn_recip32x32( tTestEngContext * context )
{
  typedef void tFxn( int32_t *  frac, 
                  int16_t *exp, 
                  const int32_t * x, 
                  int N);
  tVec frac,exponent;
  tFxn *fxn;
  int32_t *X; 
  float64_t *Z;
  int32_t* pFrac;
  int16_t *pExp;
  int n,N;
  ASSERT( context && context->target.fut );
  X = vecGetElem_fr32( &context->dataSet.X, 0 );
  Z = vecGetElem_fl64( &context->dataSet.Z, 0 );
  N   = context->args.dim[0];
  vecsAlloc(context->desc->isAligned,FMT_FRACT32,&frac,N<0?0:N,NULL);
  vecsAlloc(context->desc->isAligned,FMT_FRACT16,&exponent,N<0?0:N,NULL);
  pFrac=vecGetElem_fr32(&frac,0);
  pExp =vecGetElem_fr16(&exponent,0);
  fxn = (tFxn*)context->target.fut;
  /* 
    test for overlapping X and pFrac for each 1st of 4 tests if overlapping 
    is enabled in extraParams
  */
  if (((context->args.caseNum & 3)==1) && (context->desc->extraParam & OVLP_XZ))
  {
      size_t szX,szfrac;
      szX=vecGetSize(&context->dataSet.X);
      szfrac=vecGetSize(&frac);
      (void)szfrac;
      ASSERT(szX==szfrac);
      memcpy(pFrac,X,szX);
      X=pFrac;
  }
  te_vReportStd(context);
  fxn( pFrac, pExp, X, N );
  /* convert results to single precision floating point */
  for(n=0; n<N; n++)
  {
      Z[n]=STDLIB_MATH(ldexp)((float64_t)pFrac[n],pExp[n]-31);
  }
  vecsFree(&frac,&exponent,NULL);
} /* te_math_processFxn_recip32x32() */

/* Apply a function to the test case data set:
 * scalar functions with single argument, e.g. cos() */
void te_math_processFxn_scl_vXvZ(tTestEngContext * context)
{
  typedef float16_t tFxn_fl16( float16_t);
  typedef float32_t tFxn_fl32( float32_t );
  typedef float64_t tFxn_fl64( float64_t );
  typedef fract16   tFxn_fr16( fract16   );
  typedef fract32   tFxn_fr32( fract32   );

  #define CALL_FXN( typeFxn, Fxn, typeXZ, X, Z, n )           \
      { te_errh_resetStates( context );                       \
        *(typeXZ*)(Z) = ( (typeFxn*)(Fxn) )( *(typeXZ*)(X) ); \
        te_errh_verifyStates( context, n );                   \
        (X) = (typeXZ*)(X) + 1; (Z) = (typeXZ*)(Z) + 1; }

  tTestEngTarget  fxn;
  void *X, *Z;
  int n, N;

  X = vecGetElem( &context->dataSet.X, 0 );
  Z = vecGetElem( &context->dataSet.Z, 0 );

  fxn = context->target.fut;
  N   = context->args.dim[0];
  te_vReportStd(context);

  switch ( context->desc->fmt & FMT_DTYPE_MASK )
  {
  case FMT_FLOAT16: for ( n=0; n<N; n++ ) CALL_FXN( tFxn_fl16, fxn, float16_t, X, Z, n ); break;
  case FMT_FLOAT32: for ( n=0; n<N; n++ ) CALL_FXN( tFxn_fl32, fxn, float32_t, X, Z, n ); break;
  case FMT_FLOAT64: for ( n=0; n<N; n++ ) CALL_FXN( tFxn_fl64, fxn, float64_t, X, Z, n ); break;
  case FMT_FRACT16: for ( n=0; n<N; n++ ) CALL_FXN( tFxn_fr16, fxn, fract16  , X, Z, n ); break;
  case FMT_FRACT32: for ( n=0; n<N; n++ ) CALL_FXN( tFxn_fr32, fxn, fract32  , X, Z, n ); break;
  default: ASSERT( 0 );
  }

  #undef CALL_FXN

} /* te_math_processFxn_scl_vXvZ() */

void te_math_processFxn_scl_vXvZ32( tTestEngContext * context )
{
  typedef int32_t tFxn_fl32( float32_t );
  typedef int32_t tFxn_fl64( float64_t );
  typedef int32_t tFxn_fr16( fract16   );
  typedef int32_t tFxn_fr32( fract32   );

  #define CALL_FXN( typeFxn, Fxn, typeXZ, X, Z ) \
      { *(int32_t*)(Z) = ( (typeFxn*)(Fxn) )( *(typeXZ*)(X) ); \
        (X) = (typeXZ*)(X) + 1; (Z) = (typeXZ*)(Z) + 1; }

  tTestEngTarget  fxn;
  void *X, *Z;
  int n, N;

  X = vecGetElem( &context->dataSet.X, 0 );
  Z = vecGetElem( &context->dataSet.Z, 0 );

  fxn = context->target.fut;
  N   = context->args.dim[0];
  te_vReportStd(context);

  switch ( context->desc->fmt & FMT_DTYPE_MASK )
  {
  case FMT_FLOAT32: for ( n=0; n<N; n++ ) CALL_FXN( tFxn_fl32, fxn, float32_t, X, Z ); break;
  case FMT_FLOAT64: for ( n=0; n<N; n++ ) CALL_FXN( tFxn_fl64, fxn, float64_t, X, Z ); break;
  case FMT_FRACT16: for ( n=0; n<N; n++ ) CALL_FXN( tFxn_fr16, fxn, fract16  , X, Z ); break;
  case FMT_FRACT32: for ( n=0; n<N; n++ ) CALL_FXN( tFxn_fr32, fxn, fract32  , X, Z ); break;
  default: ASSERT( 0 );
  }

  #undef CALL_FXN

} /* te_math_processFxn_scl_vXvZ32() */

void te_math_processFxn_scl_vXvZ16( tTestEngContext * context )
{
  typedef int16_t tFxn_fl32( float32_t );
  typedef int16_t tFxn_fl64( float64_t );
  typedef int16_t tFxn_fr16( fract16   );
  typedef int16_t tFxn_fr32( fract32   );

  #define CALL_FXN( typeFxn, Fxn, typeX, X, Z ) \
      { *(int16_t*)(Z) = ( (typeFxn*)(Fxn) )( *(typeX*)(X) ); \
        (X) = (typeX*)(X) + 1; (Z) = (int16_t*)(Z) + 1; }

  tTestEngTarget  fxn;
  void *X, *Z;
  int n, N;

  X = vecGetElem( &context->dataSet.X, 0 );
  Z = vecGetElem( &context->dataSet.Z, 0 );

  fxn = context->target.fut;
  N   = context->args.dim[0];
  te_vReportStd(context);

  switch ( context->desc->fmt & FMT_DTYPE_MASK )
  {
  case FMT_FLOAT32: for ( n=0; n<N; n++ ) CALL_FXN( tFxn_fl32, fxn, float32_t, X, Z ); break;
  case FMT_FLOAT64: for ( n=0; n<N; n++ ) CALL_FXN( tFxn_fl64, fxn, float64_t, X, Z ); break;
  case FMT_FRACT16: for ( n=0; n<N; n++ ) CALL_FXN( tFxn_fr16, fxn, fract16  , X, Z ); break;
  case FMT_FRACT32: for ( n=0; n<N; n++ ) CALL_FXN( tFxn_fr32, fxn, fract32  , X, Z ); break;
  default: ASSERT( 0 );
  }

  #undef CALL_FXN

} /* te_math_processFxn_scl_vXvZ16() */
void te_math_processFxn_scl_vX32sY32vZ( tTestEngContext * context )
{
  typedef float32_t tFxn_fl32( int32_t,int32_t );
  typedef float64_t tFxn_fl64( int32_t,int32_t );
  typedef fract16   tFxn_fr16( int32_t,int32_t );
  typedef fract32   tFxn_fr32( int32_t,int32_t );

  #define CALL_FXN( typeFxn, Fxn, typeX,typeY,typeZ, X,Y,Z ) \
      { *(typeZ*)(Z) = ( (typeFxn*)(Fxn) )( *(typeX*)(X),*(typeY*)(Y) ); \
        (X) = (typeX*)(X) + 1; (Y) = (typeY*)(Y) + 0; (Z) = (typeZ*)(Z) + 1; }

  tTestEngTarget  fxn;
  void *X, *Y, *Z;
  int n, N;

  X = vecGetElem( &context->dataSet.X, 0 );
  Y = vecGetElem( &context->dataSet.Y, 0 );
  Z = vecGetElem( &context->dataSet.Z, 0 );

  fxn = context->target.fut;
  N   = context->args.dim[0];
  te_vReportStd(context);

  switch ( context->desc->fmt & FMT_DTYPE_MASK )
  {
  case FMT_FLOAT32: for ( n=0; n<N; n++ ) CALL_FXN( tFxn_fl32, fxn, int32_t,int32_t,float32_t, X,Y,Z ); break;
  case FMT_FLOAT64: for ( n=0; n<N; n++ ) CALL_FXN( tFxn_fl64, fxn, int32_t,int32_t,float64_t, X,Y,Z ); break;
  case FMT_FRACT16: for ( n=0; n<N; n++ ) CALL_FXN( tFxn_fr16, fxn, int32_t,int32_t,fract16  , X,Y,Z ); break;
  case FMT_FRACT32: for ( n=0; n<N; n++ ) CALL_FXN( tFxn_fr32, fxn, int32_t,int32_t,fract32  , X,Y,Z ); break;
  default: ASSERT( 0 );
  }

  #undef CALL_FXN

} /* te_math_processFxn_scl_vX32sY32vZ() */

void te_math_processFxn_scl_vXsY32vZ32( tTestEngContext * context )
{
  typedef int32_t tFxn_fl32( float32_t,int32_t );
  typedef int32_t tFxn_fl64( float64_t,int32_t );
  typedef int32_t tFxn_fr16( fract16  ,int32_t );
  typedef int32_t tFxn_fr32( fract32  ,int32_t );

  #define CALL_FXN( typeFxn, Fxn, typeX,typeY,typeZ, X,Y,Z ) \
      { *(typeZ*)(Z) = ( (typeFxn*)(Fxn) )( *(typeX*)(X),*(typeY*)(Y) ); \
        (X) = (typeX*)(X) + 1; (Y) = (typeY*)(Y) + 0; (Z) = (typeZ*)(Z) + 1; }

  tTestEngTarget  fxn;
  void *X, *Y, *Z;
  int n, N;

  X = vecGetElem( &context->dataSet.X, 0 );
  Y = vecGetElem( &context->dataSet.Y, 0 );
  Z = vecGetElem( &context->dataSet.Z, 0 );

  fxn = context->target.fut;
  N   = context->args.dim[0];
  te_vReportStd(context);

  switch ( context->desc->fmt & FMT_DTYPE_MASK )
  {
  case FMT_FLOAT32: for ( n=0; n<N; n++ ) CALL_FXN( tFxn_fl32, fxn, float32_t,int32_t,int32_t, X,Y,Z ); break;
  case FMT_FLOAT64: for ( n=0; n<N; n++ ) CALL_FXN( tFxn_fl64, fxn, float64_t,int32_t,int32_t, X,Y,Z ); break;
  case FMT_FRACT16: for ( n=0; n<N; n++ ) CALL_FXN( tFxn_fr16, fxn, fract16  ,int32_t,int32_t, X,Y,Z ); break;
  case FMT_FRACT32: for ( n=0; n<N; n++ ) CALL_FXN( tFxn_fr32, fxn, fract32  ,int32_t,int32_t, X,Y,Z ); break;
  default: ASSERT( 0 );
  }

  #undef CALL_FXN

} /* te_math_processFxn_scl_vXsY32vZ32() */

/* Apply a function to the test case data set:
 * scalar functions with single argument, e.g. cos(). Input X is complex, output is real*/
void te_math_processFxn_scl_vXcvZ( tTestEngContext * context )
{
  typedef float32_t tFxn_fl32( complex_float  );
  typedef float64_t tFxn_fl64( complex_double );
  typedef fract16   tFxn_fr16( complex_fract16);
  typedef fract32   tFxn_fr32( complex_fract32);

  #define CALL_FXN( typeFxn, Fxn, typeX,typeZ, X, Z ) \
      { *(typeZ*)(Z) = ( (typeFxn*)(Fxn) )( *(typeX*)(X) ); \
        (X) = (typeX*)(X) + 1; (Z) = (typeZ*)(Z) + 1; }

  tTestEngTarget  fxn;
  void *X, *Z;
  int n, N;

  X = vecGetElem( &context->dataSet.X, 0 );
  Z = vecGetElem( &context->dataSet.Z, 0 );

  fxn = context->target.fut;
  N   = context->args.dim[0];
  te_vReportStd(context);
  switch ( context->desc->fmt & FMT_DTYPE_MASK )
  {
  case FMT_FLOAT32: for ( n=0; n<N; n++ ) CALL_FXN( tFxn_fl32, fxn, complex_float  ,float32_t, X, Z ); break;
  case FMT_FLOAT64: for ( n=0; n<N; n++ ) CALL_FXN( tFxn_fl64, fxn, complex_double ,float64_t, X, Z ); break;
  case FMT_FRACT16: for ( n=0; n<N; n++ ) CALL_FXN( tFxn_fr16, fxn, complex_fract16,fract16  , X, Z ); break;
  case FMT_FRACT32: for ( n=0; n<N; n++ ) CALL_FXN( tFxn_fr32, fxn, complex_fract32,fract32  , X, Z ); break;
  default: ASSERT( 0 );
  }

  #undef CALL_FXN

} /* te_math_processFxn_scl_vXcvZ() */


/* Apply a function to the test case data set:
 * scalar atan2(y,x) */
void te_math_processFxn_scl_atan2( tTestEngContext * context )
{
  typedef float32_t tFxn_fl32( float32_t, float32_t );
  typedef float64_t tFxn_fl64( float64_t, float64_t );
  typedef fract16   tFxn_fr16( fract16  , fract16   );
  typedef fract32   tFxn_fr32( fract32  , fract32   );

  #define CALL_FXN( typeFxn, Fxn, typeXYZ, Y, X, Z, n )                         \
      { te_errh_resetStates( context );                                         \
        *(typeXYZ*)(Z) = ( (typeFxn*)(Fxn) )( *(typeXYZ*)(Y), *(typeXYZ*)(X) ); \
        te_errh_verifyStates( context, n );                                     \
        (Y) = (typeXYZ*)(Y) + 1; (X) = (typeXYZ*)(X) + 1; (Z) = (typeXYZ*)(Z) + 1; }

  tTestEngTarget  fxn;
  void *X, *Y, *Z;
  int n, N;

  X = vecGetElem( &context->dataSet.X, 0 );
  Y = vecGetElem( &context->dataSet.Y, 0 );
  Z = vecGetElem( &context->dataSet.Z, 0 );

  fxn = context->target.fut;
  N   = context->args.dim[0];
  te_vReportStd(context);

  switch ( context->desc->fmt & FMT_DTYPE_MASK )
  {
  case FMT_FLOAT32: for ( n=0; n<N; n++ ) CALL_FXN( tFxn_fl32, fxn, float32_t, Y, X, Z, n ); break;
  case FMT_FLOAT64: for ( n=0; n<N; n++ ) CALL_FXN( tFxn_fl64, fxn, float64_t, Y, X, Z, n ); break;
  case FMT_FRACT16: for ( n=0; n<N; n++ ) CALL_FXN( tFxn_fr16, fxn, fract16  , Y, X, Z, n ); break;
  case FMT_FRACT32: for ( n=0; n<N; n++ ) CALL_FXN( tFxn_fr32, fxn, fract32  , Y, X, Z, n ); break;
  default: ASSERT( 0 );
  }

  #undef CALL_FXN

} /* te_math_processFxn_scl_atan2() */

void te_math_processFxn_scl_dividef( tTestEngContext * context )
{
  typedef float32_t tFxn_fl32( float32_t, float32_t );

  #define CALL_FXN( typeFxn, Fxn, typeXYZ, Y, X, Z, n )                         \
      { te_errh_resetStates( context );                                         \
        *(typeXYZ*)(Z) = ( (typeFxn*)(Fxn) )( *(typeXYZ*)(X), *(typeXYZ*)(Y) ); \
        te_errh_verifyStates( context, n );                                     \
        (Y) = (typeXYZ*)(Y) + 1; (X) = (typeXYZ*)(X) + 1; (Z) = (typeXYZ*)(Z) + 1; }

  tTestEngTarget  fxn;
  void *X, *Y, *Z;
  int n, N;

  X = vecGetElem( &context->dataSet.X, 0 );
  Y = vecGetElem( &context->dataSet.Y, 0 );
  Z = vecGetElem( &context->dataSet.Z, 0 );

  fxn = context->target.fut;
  N   = context->args.dim[0];
  te_vReportStd(context);

  switch ( context->desc->fmt & FMT_DTYPE_MASK )
  {
  case FMT_FLOAT32: for ( n=0; n<N; n++ ) CALL_FXN( tFxn_fl32, fxn, float32_t, Y, X, Z, n ); break;
  default: ASSERT( 0 );
  }

  #undef CALL_FXN

} /* te_math_processFxn_scl_dividef() */

/* processing function for scl_recip16x16 */
void te_math_processFxn_scl_recip16x16( tTestEngContext * context )
{
  typedef uint32_t tFxn( int16_t);
  tFxn *fxn;
  int16_t *X; 
  float32_t *Z;
  int n,N;
  ASSERT( context && context->target.fut );
  X = vecGetElem_fr16( &context->dataSet.X, 0 );
  Z = vecGetElem_fl32( &context->dataSet.Z, 0 );
  N   = context->args.dim[0];
  te_vReportStd(context);
  fxn = (tFxn*)context->target.fut;
  /* convert results to single precision floating point */
  for(n=0; n<N; n++)
  {
      uint32_t z;
      int16_t  frac,exponent;
      z=fxn(X[n]);
      frac=(int16_t)z;
      exponent=(int16_t)((int32_t)z>>16);
      Z[n]=STDLIB_MATH(ldexpf)(frac,exponent-15);
  }
} /* te_math_processFxn_scl_recip16x16() */

/* processing function for scl_recip32x32 scl_recip24x24 */
void te_math_processFxn_scl_recip32x32( tTestEngContext * context )
{
  typedef uint32_t tFxn( int32_t);
  tFxn *fxn;
  int32_t *X; 
  float64_t *Z;
  int n,N;
  ASSERT( context && context->target.fut );
  X = vecGetElem_fr32( &context->dataSet.X, 0 );
  Z = vecGetElem_fl64( &context->dataSet.Z, 0 );
  N   = context->args.dim[0];
  te_vReportStd(context);
  fxn = (tFxn*)context->target.fut;
  /* convert results to single precision floating point */
  for(n=0; n<N; n++)
  {
      uint32_t z;
      int32_t  frac;
      int16_t exponent;
      z=fxn(X[n]);
      frac=(int32_t)(z<<8);
      exponent=(int16_t)((int32_t)z>>24);
      Z[n]=STDLIB_MATH(ldexp)((float64_t)frac,exponent-31);
  }
} /* te_math_processFxn_scl_recip32x32() */

/* processing function for scl_divide32x32 scl_divide24x24 */
void te_math_processFxn_scl_divide32x32( tTestEngContext * context )
{
  typedef uint32_t tFxn( int32_t, int32_t);
  tFxn *fxn;
  int32_t *X; 
  int32_t *Y; 
  float64_t *Z;
  int n,N;
  ASSERT( context && context->target.fut );
  X = vecGetElem_fr32( &context->dataSet.X, 0 );
  Y = vecGetElem_fr32( &context->dataSet.Y, 0 );
  Z = vecGetElem_fl64( &context->dataSet.Z, 0 );
  N   = context->args.dim[0];
  te_vReportStd(context);
  fxn = (tFxn*)context->target.fut;
  /* convert results to single precision floating point */
  for(n=0; n<N; n++)
  {
      uint32_t z;
      int32_t  frac;
      int16_t exponent;
      z=fxn(X[n],Y[n]);
      frac=(int32_t)(z<<8);
      exponent=(int16_t)((int32_t)z>>24);
      Z[n]=STDLIB_MATH(ldexp)((float64_t)frac,exponent-31);
  }
} /* te_math_processFxn_scl_divide32x32() */

/* processing function for scl_divide16x16 */
void te_math_processFxn_scl_divide16x16( tTestEngContext * context )
{
  typedef uint32_t tFxn( int16_t, int16_t);
  tFxn *fxn;
  int16_t *X; 
  int16_t *Y; 
  float32_t *Z;
  int n,N;
  ASSERT( context && context->target.fut );
  X = vecGetElem_fr16( &context->dataSet.X, 0 );
  Y = vecGetElem_fr16( &context->dataSet.Y, 0 );
  Z = vecGetElem_fl32( &context->dataSet.Z, 0 );
  N   = context->args.dim[0];
  fxn = (tFxn*)context->target.fut;
  te_vReportStd(context);
  /* convert results to single precision floating point */
  for(n=0; n<N; n++)
  {
      uint32_t z;
      int16_t frac;
      int16_t exponent;
      z=fxn(X[n],Y[n]);
      frac=(int16_t)(z);
      exponent=(int16_t)((int32_t)z>>16);
      Z[n]=STDLIB_MATH(ldexpf)(frac,exponent-15);
  }
} /* te_math_processFxn_scl_divide16x16() */

/* processing function for scl_divide64x32 */
void te_math_processFxn_scl_divide64x32( tTestEngContext * context )
{
  typedef int32_t tFxn( int64_t, int32_t);
  tFxn *fxn;
  int64_t *X; 
  int32_t *Y; 
  int32_t *Z;
  int n,N;
  ASSERT( context && context->target.fut );
  X = vecGetElem_i64( &context->dataSet.X, 0 );
  Y = vecGetElem_fr32( &context->dataSet.Y, 0 );
  Z = vecGetElem_fr32( &context->dataSet.Z, 0 );
  N   = context->args.dim[0];
  fxn = (tFxn*)context->target.fut;
  te_vReportStd(context);
  /* convert results to single precision floating point */
  for(n=0; n<N; n++)
  {
      Z[n]=fxn(X[n],Y[n]);
  }
} /* te_math_processFxn_scl_divide64x32() */

/* load function for 64x64 reciprocal */
int te_math_loadFxn_recip64x64(tTestEngContext* context)
{
    int N;
    int res = 0;
    ASSERT( context && context->seqFile );
    memset( &context->dataSet, 0, sizeof(context->dataSet) );
    N = MAX( 0, context->args.dim[0] );
  /* Allocate input 
     X 
     reference output (Y,U) in pseudofloating point (64-bit mantissa/16-bit exponent 
     bounds Z,Zlo,Zhi (difference between real output and reference) in single precision
     */
    if (
        vecAlloc( &context->dataSet.X  , N, context->desc->isAligned, FMT_REAL|FMT_INT64, 0 )+
        vecAlloc( &context->dataSet.Y  , N, context->desc->isAligned, FMT_REAL|FMT_INT64, 0 )+
        vecAlloc( &context->dataSet.U  , N, context->desc->isAligned, FMT_REAL|FMT_INT16, 0 )+
        vecAlloc( &context->dataSet.Z  , N, context->desc->isAligned, FMT_REAL|FMT_FLOAT32, 0 )+
        vecAlloc( &context->dataSet.Zlo, N, context->desc->isAligned, FMT_REAL|FMT_FLOAT32, 0 )+
        vecAlloc( &context->dataSet.Zhi, N, context->desc->isAligned, FMT_REAL|FMT_FLOAT32, 0 )
        !=6)
    {
        printf( "te_math_loadFxn_recip64x64(): failed to allocate memory\n");
    }
    else
    {
        if ( !seqFileReadVecs( context->seqFile,
                                &context->dataSet.X,
                                &context->dataSet.Y,
                                &context->dataSet.U,
                                &context->dataSet.Zlo,
                                &context->dataSet.Zhi, 0 ) )
        {
            printf( "te_math_loadFxn_recip64x64(): failed to read vectors data\n");
        }
        else 
        {
            res=1;
        }
    }
    /* Free vectors data if failed. */
    if ( !res )
    {
        te_freeVectors(context);
    }
    return (res);

} /* te_math_loadFxn_recip64x64() */

/* processing function for vec_recip64x64 */
void te_math_processFxn_recip64x64( tTestEngContext * context, int isScalar )
{
  typedef void tFxnVec(int64_t *  frac, int16_t *exp, const int64_t * x, int N);
  typedef uint64_t tFxnScl(int64_t x);
  tVec vmant,vexp;
  const int64_t *X; 
  const int64_t *ref_mant; 
  const int16_t *ref_exp; 
        int64_t *mant; 
        int16_t *exp; 
  float32_t *Z;
  int n,N;
  ASSERT( context && context->target.fut );
  X = vecGetElem_i64( &context->dataSet.X, 0 );
  ref_mant= vecGetElem_i64( &context->dataSet.Y, 0 );
  ref_exp = vecGetElem_i16( &context->dataSet.U, 0 );
  Z = vecGetElem_fl32( &context->dataSet.Z, 0 );
  N   = context->args.dim[0];
  if(
        vecAlloc( &vmant , MAX(0,N), context->desc->isAligned, FMT_REAL|FMT_INT64, 0 )+
        vecAlloc( &vexp  , MAX(0,N), context->desc->isAligned, FMT_REAL|FMT_INT16, 0 )
        !=2)
  {
        printf( "te_math_processFxn_recip64x64(): failed to allocate memory\n");
  }
  mant= vecGetElem_i64(&vmant, 0 );
  exp = vecGetElem_i16(&vexp , 0 );
  
  te_vReportStd(context);

  if (isScalar)
  {
      tFxnScl *fxn;
      fxn =(tFxnScl*)context->target.fut;
      uint64_t res;
      NASSERT(N==1);
      res=(*fxn)(X[0]);
      mant[0]=res<<8;
      exp[0]=(int16_t)(res>>56);
      if (exp[0]==64)
      {
          mant[0]|=0xff;
      }
  }
  else
  {
      ((tFxnVec*)context->target.fut)(mant,exp,X, N );
  }
  /* convert results to single precision floating point */
  for(n=0; n<N; n++)
  {
      int e=exp[n],e_ref=ref_exp[n],d;
      int64_t m=mant[n];
      int64_t m_ref=ref_mant[n];
      int64_t diff;
      d=MAX(e,e_ref);
      m=m>>(d-e);
      m_ref=m_ref>>(d-e_ref);
      diff=m-m_ref;
      Z[n]=STDLIB_MATH(ldexpf)((float32_t)diff,d-63);
  }
  vecsFree(&vmant,&vexp,NULL);
} /* te_math_processFxn_recip16x16() */

void te_math_processFxn_vec_recip64x64( tTestEngContext * context )
{
    te_math_processFxn_recip64x64(context,0);
}

void te_math_processFxn_scl_recip64x64( tTestEngContext * context )
{
    te_math_processFxn_recip64x64(context,1);
}

/* load function for 64x64 division */
int te_math_loadFxn_divide64x64(tTestEngContext* context)
{
    int N;
    int res = 0;
    ASSERT( context && context->seqFile );
    memset( &context->dataSet, 0, sizeof(context->dataSet) );
    N = MAX( 0, context->args.dim[0] );
  /* Allocate input 
     X 
     reference output (U,V) in pseudofloating point (64-bit mantissa/16-bit exponent 
     bounds Z,Zlo,Zhi (difference between real output and reference) in single precision
     */
    if (
        vecAlloc( &context->dataSet.X  , N, context->desc->isAligned, FMT_REAL|FMT_INT64, 0 )+
        vecAlloc( &context->dataSet.Y  , N, context->desc->isAligned, FMT_REAL|FMT_INT64, 0 )+
        vecAlloc( &context->dataSet.U  , N, context->desc->isAligned, FMT_REAL|FMT_INT64, 0 )+
        vecAlloc( &context->dataSet.V  , N, context->desc->isAligned, FMT_REAL|FMT_INT16, 0 )+
        vecAlloc( &context->dataSet.Z  , N, context->desc->isAligned, FMT_REAL|FMT_FLOAT32, 0 )+
        vecAlloc( &context->dataSet.Zlo, N, context->desc->isAligned, FMT_REAL|FMT_FLOAT32, 0 )+
        vecAlloc( &context->dataSet.Zhi, N, context->desc->isAligned, FMT_REAL|FMT_FLOAT32, 0 )
        !=7)
    {
        printf( "te_math_loadFxn_divide64x64(): failed to allocate memory\n");
    }
    else
    {
        if ( !seqFileReadVecs( context->seqFile,
                                &context->dataSet.X,
                                &context->dataSet.Y,
                                &context->dataSet.U,
                                &context->dataSet.V,
                                &context->dataSet.Zlo,
                                &context->dataSet.Zhi, 0 ) )
        {
            printf( "te_math_loadFxn_divide64x64(): failed to read vectors data\n");
        }
        else 
        {
            res=1;
        }
    }
    /* Free vectors data if failed. */
    if ( !res )
    {
        te_freeVectors(context);
    }
    return (res);

} /* te_math_loadFxn_divide64x64() */

/* processing function for vec_recip64x64 */
void te_math_processFxn_divide64x64( tTestEngContext * context, int isScalar )
{
  typedef void tFxnVec(int64_t *  frac, int16_t *exp, const int64_t * x, const int64_t * y, int N);
  typedef uint64_t tFxnScl(int64_t x,int64_t y);
  tVec vmant,vexp;
  const int64_t *X; 
  const int64_t *Y; 
  const int64_t *ref_mant; 
  const int16_t *ref_exp; 
        int64_t *mant; 
        int16_t *exp; 
  float32_t *Z;
  int n,N;
  ASSERT( context && context->target.fut );
  X = vecGetElem_i64( &context->dataSet.X, 0 );
  Y = vecGetElem_i64( &context->dataSet.Y, 0 );
  ref_mant= vecGetElem_i64( &context->dataSet.U, 0 );
  ref_exp = vecGetElem_i16( &context->dataSet.V, 0 );
  Z = vecGetElem_fl32( &context->dataSet.Z, 0 );
  N   = context->args.dim[0];
  if(
        vecAlloc( &vmant , MAX(0,N), context->desc->isAligned, FMT_REAL|FMT_INT64, 0 )+
        vecAlloc( &vexp  , MAX(0,N), context->desc->isAligned, FMT_REAL|FMT_INT16, 0 )
        !=2)
  {
        printf( "te_math_processFxn_recip64x64(): failed to allocate memory\n");

  }
  mant= vecGetElem_i64(&vmant, 0 );
  exp = vecGetElem_i16(&vexp , 0 );
  
  te_vReportStd(context);

  if (isScalar)
  {
      tFxnScl *fxn;
      fxn =(tFxnScl*)context->target.fut;
      uint64_t res;
      NASSERT(N==1);
      res=(*fxn)(X[0],Y[0]);
      mant[0]=res<<8;
      exp[0]=(int16_t)((int64_t)res>>56);
      if (exp[0]==64)
      {
          if(X[0]>=0) mant[0]|=0xff;
      }
  }
  else
  {
      ((tFxnVec*)context->target.fut)(mant,exp,X, Y, N );
  }
  /* convert results to single precision floating point */
  for(n=0; n<N; n++)
  {
      int e=exp[n],e_ref=ref_exp[n],d;
      int64_t m=mant[n];
      int64_t m_ref=ref_mant[n];
      int64_t diff;
      d=MAX(e,e_ref);
      m=m>>(d-e);
      m_ref=m_ref>>(d-e_ref);
      diff=m-m_ref;
      Z[n]=STDLIB_MATH(ldexpf)((float32_t)diff,d-63);
  }
  vecsFree(&vmant,&vexp,NULL);
} /* te_math_processFxn_recip16x16() */

void te_math_processFxn_vec_divide64x64( tTestEngContext * context )
{
    te_math_processFxn_divide64x64(context,0);
}

void te_math_processFxn_scl_divide64x64( tTestEngContext * context )
{
    te_math_processFxn_divide64x64(context,1);
}

/* processing function for vec_pow_32x32 */
void te_math_processFxn_pow32x32( tTestEngContext * context )
{
  typedef void tFxn( int32_t *  frac, 
                  int16_t *exp, 
                  const int32_t *  x, 
                  const int32_t *  y, 
                  int N);
  tVec frac,exponent;
  tFxn *fxn;
  const int32_t *X; 
  const int32_t *Y; 
  float64_t *Z;
  int32_t* pFrac;
  int16_t *pExp;
  int n,N;
  ASSERT( context && context->target.fut );
  X = vecGetElem_fr32( &context->dataSet.X, 0 );
  Y = vecGetElem_fr32( &context->dataSet.Y, 0 );
  Z = vecGetElem_fl64( &context->dataSet.Z, 0 );
  N   = context->args.dim[0];
  vecsAlloc(context->desc->isAligned,FMT_FRACT32,&frac,N<0?0:N,NULL);
  vecsAlloc(context->desc->isAligned,FMT_FRACT16,&exponent,N<0?0:N,NULL);
  pFrac=vecGetElem_fr32(&frac,0);
  pExp =vecGetElem_fr16(&exponent,0);
  fxn = (tFxn*)context->target.fut;
  te_vReportStd(context);
#if 1
  /* 
    test for overlapping X and pFrac for each 1st of 4 tests if overlapping 
    is enabled in extraParams
  */
  if (((context->args.caseNum & 3)==1) && (context->desc->extraParam & OVLP_XZ))
  {
      size_t szX,szfrac;
      szX=vecGetSize(&context->dataSet.X);
      szfrac=vecGetSize(&frac);
      (void)szfrac;
      ASSERT(szX==szfrac);
      memcpy(pFrac,X,szX);
      X=pFrac;
  }
  /* 
    test for overlapping X and pExp for each 2nd of 4 tests if overlapping 
    is enabled in extraParams
  */
  if (((context->args.caseNum & 3)==2) && (context->desc->extraParam & OVLP_YZ))
  {
      size_t szY,szfrac;
      szY=vecGetSize(&context->dataSet.Y);
      szfrac=vecGetSize(&frac);
      (void)szfrac;
      ASSERT(szY==szfrac);
      memcpy(pFrac,Y,szY);
      Y=pFrac;
  }
#endif
  fxn( pFrac, pExp, X,Y, N );
  /* convert results to single precision floating point */
  for(n=0; n<N; n++)
  {
      Z[n]=STDLIB_MATH(ldexp)(pFrac[n],pExp[n]-31);
  }
  vecsFree(&frac,&exponent,NULL);
} /* te_math_processFxn_pow32x32() */
/* Apply a function to the test case data set:
 * scalar functions with two arguments, e.g. powf() */
void te_math_processFxn_scl_vXvYvZ(tTestEngContext * context)
{
  typedef float32_t tFxn_fl32( float32_t, float32_t );
  typedef int16_t   tFxn_fr16( int16_t, int16_t );
  typedef int32_t   tFxn_fr32( int32_t, int32_t );

  #define CALL_FXN( typeFxn, Fxn, typeXYZ, X, Y, Z, n )                         \
      { te_errh_resetStates( context );                                         \
        *(typeXYZ*)(Z) = ( (typeFxn*)(Fxn) )( *(typeXYZ*)(X), *(typeXYZ*)(Y) ); \
        te_errh_verifyStates( context, n );                                     \
        (X) = (typeXYZ*)(X) + 1; (Y) = (typeXYZ*)(Y) + 1; (Z) = (typeXYZ*)(Z) + 1; }

  tTestEngTarget  fxn;
  void *X, *Y, *Z;
  int n, N;

  te_vReportStd(context);
  X = vecGetElem( &context->dataSet.X, 0 );
  Y = vecGetElem( &context->dataSet.Y, 0 );
  Z = vecGetElem( &context->dataSet.Z, 0 );

  fxn = context->target.fut;
  N = context->args.dim[0];

  switch ( context->desc->fmt & FMT_DTYPE_MASK )
  {
  case FMT_FLOAT32: for ( n=0; n<N; n++ ) CALL_FXN( tFxn_fl32, fxn, float32_t, X, Y, Z, n ); break;
  case FMT_FRACT16: for ( n=0; n<N; n++ ) CALL_FXN( tFxn_fr16, fxn, fract16, X, Y, Z, n ); break;
  case FMT_FRACT32: for ( n=0; n<N; n++ ) CALL_FXN( tFxn_fr32, fxn, fract32, X, Y, Z, n ); break;
  default: ASSERT( 0 );
  }

  #undef CALL_FXN

} /* te_math_processFxn_s_vXvYvZ() */

