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
 * Test engine implementation
 * Collection of test case processing functions
 */

#include <string.h>

/* Cross-platform data type definitions. */
#include "types.h"
/* Test data vectors tools and SEQ-file reader. */
#include "vectools.h"
/* Test engine API. */
#include "testeng.h"

#define MIN(a,b)   ( (a)<(b) ? (a) : (b) )
#define MAX(a,b)   ( (a)>(b) ? (a) : (b) )

/* Apply the target function to the test case data set:
 * vector X (in), scalar Y (in), vector Z (out) */
void te_processFxn_vXsYvZ( tTestEngContext * context )
{
  typedef void tFxn_fr16  ( const int16_t         * x, int16_t         y, int16_t         * z, int N );
  typedef void tFxn_fr32  ( const int32_t         * x, int32_t         y, int32_t         * z, int N );
  typedef void tFxn_fl32  ( const float32_t       * x, float32_t       y, float32_t       * z, int N );
  typedef void tFxn_fl64  ( const float64_t       * x, float64_t       y, float64_t       * z, int N );
  typedef void tFxn_fr16c ( const complex_fract16 * x, complex_fract16 y, complex_fract16 * z, int N );
  typedef void tFxn_fr32c ( const complex_fract32 * x, complex_fract32 y, complex_fract32 * z, int N );
  typedef void tFxn_fl32c ( const complex_float   * x, complex_float   y, complex_float   * z, int N );
  typedef void tFxn_fl64c ( const complex_double  * x, complex_double  y, complex_double  * z, int N );

  tTestEngTarget   fxn;
  void *X, *Y, *Z;
  int N;

  ASSERT( context && context->target.fut );
  te_vReportStd(context);

  X = vecGetElem( &context->dataSet.X, 0 );
  Y = vecGetElem( &context->dataSet.Y, 0 );
  Z = vecGetElem( &context->dataSet.Z, 0 );

  N = context->args.dim[0];

  fxn = context->target.fut;
  /* 
    test for overlapping X and Z for each odd test if overlapping 
    is enabled in extraParams
  */
  if ((context->args.caseNum & 1) && (context->desc->extraParam & OVLP_XZ))
  {
      size_t szX,szZ;
      szX=vecGetSize(&context->dataSet.X);
      szZ=vecGetSize(&context->dataSet.Z);
      ASSERT(szX==szZ); (void)szZ;
      ASSERT(context->dataSet.X.fmt==context->dataSet.Z.fmt);
      memcpy(Z,X,szX);
      X=Z;
  }

  switch ( context->desc->fmt )
  {
  case FMT_REAL|FMT_FRACT16: ( (tFxn_fr16 *)fxn )( (const int16_t        *)X, *(int16_t        *)Y, (int16_t         *)Z, N ); break;
  case FMT_REAL|FMT_FRACT32: ( (tFxn_fr32 *)fxn )( (const int32_t        *)X, *(int32_t        *)Y, (int32_t         *)Z, N ); break;
  case FMT_REAL|FMT_FLOAT32: ( (tFxn_fl32 *)fxn )( (const float32_t      *)X, *(float32_t      *)Y, (float32_t       *)Z, N ); break;
  case FMT_REAL|FMT_FLOAT64: ( (tFxn_fl64 *)fxn )( (const float64_t      *)X, *(float64_t      *)Y, (float64_t       *)Z, N ); break;
  case FMT_CPLX|FMT_FRACT16: ( (tFxn_fr16c*)fxn )( (const complex_fract16*)X, *(complex_fract16*)Y, (complex_fract16 *)Z, N ); break;
  case FMT_CPLX|FMT_FRACT32: ( (tFxn_fr32c*)fxn )( (const complex_fract32*)X, *(complex_fract32*)Y, (complex_fract32 *)Z, N ); break;
  case FMT_CPLX|FMT_FLOAT32: ( (tFxn_fl32c*)fxn )( (const complex_float  *)X, *(complex_float  *)Y, (complex_float   *)Z, N ); break;
  case FMT_CPLX|FMT_FLOAT64: ( (tFxn_fl64c*)fxn )( (const complex_double *)X, *(complex_double *)Y, (complex_double  *)Z, N ); break;
  default: ASSERT( 0 );
  }

} /* te_processFxn_vXsYvZ() */

/* Apply the target function to the test case data set:
 * vector Z (out), vector X (in), scalar Y (in) */
void te_processFxn_vZvXsY( tTestEngContext * context )
{
  typedef void tFxn_fr16  ( fract16         * z, const fract16         * x, fract16         y, int N );
  typedef void tFxn_fr32  ( fract32         * z, const fract32         * x, fract32         y, int N );
  typedef void tFxn_fl32  ( float32_t       * z, const float32_t       * x, float32_t       y, int N );
  typedef void tFxn_fl64  ( float64_t       * z, const float64_t       * x, float64_t       y, int N );
  typedef void tFxn_fr16c ( complex_fract16 * z, const complex_fract16 * x, complex_fract16 y, int N );
  typedef void tFxn_fr32c ( complex_fract32 * z, const complex_fract32 * x, complex_fract32 y, int N );
  typedef void tFxn_fl32c ( complex_float   * z, const complex_float   * x, complex_float   y, int N );
  typedef void tFxn_fl64c ( complex_double  * z, const complex_double  * x, complex_double  y, int N );

  tTestEngTarget   fxn;
  void *X, *Y, *Z;
  int N;

  ASSERT( context && context->target.fut );

  X = vecGetElem( &context->dataSet.X, 0 );
  Y = vecGetElem( &context->dataSet.Y, 0 );
  Z = vecGetElem( &context->dataSet.Z, 0 );

  N   = context->args.dim[0];
  fxn = context->target.fut;
  /* 
    test for overlapping X and Z for each odd test if overlapping 
    is enabled in extraParams
  */
  if ((context->args.caseNum & 1) && (context->desc->extraParam & OVLP_XZ))
  {
      size_t szX,szZ;
      szX=vecGetSize(&context->dataSet.X);
      szZ=vecGetSize(&context->dataSet.Z);
      (void)szZ;
      ASSERT(szX==szZ);
      ASSERT(context->dataSet.X.fmt==context->dataSet.Z.fmt);
      memcpy(Z,X,szX);
      X=Z;
  }

  te_vReportStd(context);

  switch ( context->desc->fmt )
  {
  case FMT_REAL|FMT_FRACT16: ( (tFxn_fr16 *)fxn )((fract16         *)Z, (const fract16        *)X, *(fract16        *)Y,  N ); break;
  case FMT_REAL|FMT_FRACT32: ( (tFxn_fr32 *)fxn )((fract32         *)Z, (const fract32        *)X, *(fract32        *)Y,  N ); break;
  case FMT_REAL|FMT_FLOAT32: ( (tFxn_fl32 *)fxn )((float32_t       *)Z, (const float32_t      *)X, *(float32_t      *)Y,  N ); break;
  case FMT_REAL|FMT_FLOAT64: ( (tFxn_fl64 *)fxn )((float64_t       *)Z, (const float64_t      *)X, *(float64_t      *)Y,  N ); break;
  case FMT_CPLX|FMT_FRACT16: ( (tFxn_fr16c*)fxn )((complex_fract16 *)Z, (const complex_fract16*)X, *(complex_fract16*)Y,  N ); break;
  case FMT_CPLX|FMT_FRACT32: ( (tFxn_fr32c*)fxn )((complex_fract32 *)Z, (const complex_fract32*)X, *(complex_fract32*)Y,  N ); break;
  case FMT_CPLX|FMT_FLOAT32: ( (tFxn_fl32c*)fxn )((complex_float   *)Z, (const complex_float  *)X, *(complex_float  *)Y,  N ); break;
  case FMT_CPLX|FMT_FLOAT64: ( (tFxn_fl64c*)fxn )((complex_double  *)Z, (const complex_double *)X, *(complex_double *)Y,  N ); break;
  default: ASSERT( 0 );
  }

} /* te_processFxn_vZvXsY() */

/* Apply the target function to the test case data set:
 * vector Z (out), vector X (in), scalar int32 Y (in) */
void te_processFxn_vZvXsY32( tTestEngContext * context )
{
  typedef void tFxn_fr16  ( fract16         * z, const fract16         * x, int32_t y, int N );
  typedef void tFxn_fr32  ( fract32         * z, const fract32         * x, int32_t y, int N );
  typedef void tFxn_fl32  ( float32_t       * z, const float32_t       * x, int32_t y, int N );
  typedef void tFxn_fl64  ( float64_t       * z, const float64_t       * x, int32_t y, int N );
  typedef void tFxn_fr16c ( complex_fract16 * z, const complex_fract16 * x, int32_t y, int N );
  typedef void tFxn_fr32c ( complex_fract32 * z, const complex_fract32 * x, int32_t y, int N );
  typedef void tFxn_fl32c ( complex_float   * z, const complex_float   * x, int32_t y, int N );
  typedef void tFxn_fl64c ( complex_double  * z, const complex_double  * x, int32_t y, int N );

  tTestEngTarget   fxn;
  void *X, *Y, *Z;
  int N;

  ASSERT( context && context->target.fut );

  X = vecGetElem( &context->dataSet.X, 0 );
  Y = vecGetElem( &context->dataSet.Y, 0 );
  Z = vecGetElem( &context->dataSet.Z, 0 );

  N   = context->args.dim[0];
  fxn = context->target.fut;
  /* 
    test for overlapping X and Z for each odd test if overlapping 
    is enabled in extraParams
  */
  if ((context->args.caseNum & 1) && (context->desc->extraParam & OVLP_XZ))
  {
      size_t szX,szZ;
      szX=vecGetSize(&context->dataSet.X);
      szZ=vecGetSize(&context->dataSet.Z);
      (void)szZ;
      ASSERT(szX==szZ);
      ASSERT(context->dataSet.X.fmt==context->dataSet.Z.fmt);
      memcpy(Z,X,szX);
      X=Z;
  }
  te_vReportStd(context);

  switch ( context->desc->fmt )
  {
  case FMT_REAL|FMT_FRACT16: ( (tFxn_fr16 *)fxn )((fract16         *)Z, (const fract16        *)X, *(int32_t*)Y,  N ); break;
  case FMT_REAL|FMT_FRACT32: ( (tFxn_fr32 *)fxn )((fract32         *)Z, (const fract32        *)X, *(int32_t*)Y,  N ); break;
  case FMT_REAL|FMT_FLOAT32: ( (tFxn_fl32 *)fxn )((float32_t       *)Z, (const float32_t      *)X, *(int32_t*)Y,  N ); break;
  case FMT_REAL|FMT_FLOAT64: ( (tFxn_fl64 *)fxn )((float64_t       *)Z, (const float64_t      *)X, *(int32_t*)Y,  N ); break;
  case FMT_CPLX|FMT_FRACT16: ( (tFxn_fr16c*)fxn )((complex_fract16 *)Z, (const complex_fract16*)X, *(int32_t*)Y,  N ); break;
  case FMT_CPLX|FMT_FRACT32: ( (tFxn_fr32c*)fxn )((complex_fract32 *)Z, (const complex_fract32*)X, *(int32_t*)Y,  N ); break;
  case FMT_CPLX|FMT_FLOAT32: ( (tFxn_fl32c*)fxn )((complex_float   *)Z, (const complex_float  *)X, *(int32_t*)Y,  N ); break;
  case FMT_CPLX|FMT_FLOAT64: ( (tFxn_fl64c*)fxn )((complex_double  *)Z, (const complex_double *)X, *(int32_t*)Y,  N ); break;
  default: ASSERT( 0 );
  }

} /* te_processFxn_vZvXsY32() */

/* Apply the target function to the test case data set:
 * vector X (in), vector Y (in), vector Z (out) */
void te_processFxn_vXvYvZ( tTestEngContext * context )
{
  typedef void tFxn( const void * X, const void * Y, void * Z, int N );

  tFxn *fxn;
  void *X, *Y, *Z;
  int N;

  ASSERT( context && context->target.fut );
  te_vReportStd(context);

  X = vecGetElem( &context->dataSet.X, 0 );
  Y = vecGetElem( &context->dataSet.Y, 0 );
  Z = vecGetElem( &context->dataSet.Z, 0 );

  N = context->args.dim[0];

  fxn = (tFxn*)context->target.fut;
  /* 
    test for overlapping X and Z for each 1st of 4 test if overlapping 
    is enabled in extraParams
  */
  if (((context->args.caseNum & 3)==1) && (context->desc->extraParam & OVLP_XZ))
  {
      size_t szX,szZ;
      szX=vecGetSize(&context->dataSet.X);
      szZ=vecGetSize(&context->dataSet.Z);
      (void)szZ;
      ASSERT(szX==szZ);
      ASSERT(context->dataSet.X.fmt==context->dataSet.Z.fmt);
      memcpy(Z,X,szX);
      X=Z;
  }
  /* 
    test for overlapping Y and Z for each 2nd of 4 test if overlapping 
    is enabled in extraParams
  */
  if (((context->args.caseNum & 3)==2) && (context->desc->extraParam & OVLP_YZ))
  {
      size_t szY,szZ;
      szY=vecGetSize(&context->dataSet.Y);
      szZ=vecGetSize(&context->dataSet.Z);
      (void)szZ;
      ASSERT(szY==szZ);
      ASSERT(context->dataSet.Y.fmt==context->dataSet.Z.fmt);
      memcpy(Z,Y,szY);
      Y=Z;
  }
  te_errh_resetStates( context );
  fxn( X, Y, Z, N );
  te_errh_verifyStates( context, -1 );

} /* te_processFxn_vXvYvZ() */

/* Apply the target function to the test case data set:
 * vector Z (out), vector X (in), vector Y (in)  */
void te_processFxn_vZvXvY( tTestEngContext * context )
{
  typedef void tFxn( const void * Z, const void * X, void * Y, int N );

  tFxn *fxn;
  void *X, *Y, *Z;
  int N;

  ASSERT( context && context->target.fut );
  te_vReportStd(context);

  X = vecGetElem( &context->dataSet.X, 0 );
  Y = vecGetElem( &context->dataSet.Y, 0 );
  Z = vecGetElem( &context->dataSet.Z, 0 );

  N = context->args.dim[0];

  fxn = (tFxn*)context->target.fut;
  /* 
    test for overlapping X and Z for each 1st of 4 test if overlapping 
    is enabled in extraParams
  */
  if (((context->args.caseNum & 3)==1) && (context->desc->extraParam & OVLP_XZ))
  {
      size_t szX,szZ;
      szX=vecGetSize(&context->dataSet.X);
      szZ=vecGetSize(&context->dataSet.Z);
      ASSERT(szX==szZ); (void)szZ;
      ASSERT(context->dataSet.X.fmt==context->dataSet.Z.fmt);
      memcpy(Z,X,szX);
      X=Z;
  }
  /* 
    test for overlapping Y and Z for each 2nd of 4 test if overlapping 
    is enabled in extraParams
  */
  if (((context->args.caseNum & 3)==2) && (context->desc->extraParam & OVLP_YZ))
  {
      size_t szY,szZ;
      szY=vecGetSize(&context->dataSet.Y);
      szZ=vecGetSize(&context->dataSet.Z);
      ASSERT(szY==szZ); (void)szZ;
      ASSERT(context->dataSet.Y.fmt==context->dataSet.Z.fmt);
      memcpy(Z,Y,szY);
      Y=Z;
  }

  te_errh_resetStates( context );
  fxn( Z, X, Y, N );
  te_errh_verifyStates( context, -1 );

} /* te_processFxn_vZvXvY() */

/* Apply the target function to the test case data set:
 * vector Z (out), vector W (out), vector X (in) */
void te_processFxn_vZvWvX( tTestEngContext * context ) 
{
  typedef void tFxn( void * Z, void * W, const void * X, int N );

  tFxn *fxn;
  void *X, *Z, *W;
  int N;

  ASSERT( context && context->target.fut );
  te_vReportStd(context);

  X = vecGetElem( &context->dataSet.X, 0 );
  Z = vecGetElem( &context->dataSet.Z, 0 );
  W = vecGetElem( &context->dataSet.W, 0 );

  N = context->args.dim[0];

  fxn = (tFxn*)context->target.fut;

  fxn( Z, W, X, N );

} /* te_processFxn_vZvWvX() */

/* Apply the target function to the test case data set:
 * vector Z (out), vector W (out), vector X (in), vector Y (in) */
void te_processFxn_vZvWvXvY( tTestEngContext * context ) 
{
  typedef void tFxn( void * Z, void * W, const void * X, const void * Y, int N );

  tFxn *fxn;
  void *X, *Y, *Z, *W;
  int N;

  ASSERT( context && context->target.fut );
  te_vReportStd(context);

  X = vecGetElem( &context->dataSet.X, 0 );
  Y = vecGetElem( &context->dataSet.Y, 0 );
  Z = vecGetElem( &context->dataSet.Z, 0 );
  W = vecGetElem( &context->dataSet.W, 0 );

  N = context->args.dim[0];

  fxn = (tFxn*)context->target.fut;

  fxn( Z, W, X, Y, N );

} /* te_processFxn_vZvWvXvY() */

/* Apply the target function to the test case data set:
 * vector Y (in), vector X (in), vector Z (out) */
void te_processFxn_vYvXvZ( tTestEngContext * context )
{
  typedef void tFxn( const void * Y, const void * X, void * Z, int N );

  tFxn *fxn;
  void *X, *Y, *Z;
  int N;

  ASSERT( context && context->target.fut );
  te_vReportStd(context);

  X = vecGetElem( &context->dataSet.X, 0 );
  Y = vecGetElem( &context->dataSet.Y, 0 );
  Z = vecGetElem( &context->dataSet.Z, 0 );

  N = context->args.dim[0];

  fxn = (tFxn*)context->target.fut;
  /* 
    test for overlapping X and Z for each 1st of 4 test if overlapping 
    is enabled in extraParams
  */
  if (((context->args.caseNum & 3)==1) && (context->desc->extraParam & OVLP_XZ))
  {
      size_t szX,szZ;
      szX=vecGetSize(&context->dataSet.X);
      szZ=vecGetSize(&context->dataSet.Z);
      ASSERT(szX==szZ); (void)szZ;
      ASSERT(context->dataSet.X.fmt==context->dataSet.Z.fmt);
      memcpy(Z,X,szX);
      X=Z;
  }
  /* 
    test for overlapping Y and Z for each 2nd of 4 test if overlapping 
    is enabled in extraParams
  */
  if (((context->args.caseNum & 3)==2) && (context->desc->extraParam & OVLP_YZ))
  {
      size_t szY,szZ;
      szY=vecGetSize(&context->dataSet.Y);
      szZ=vecGetSize(&context->dataSet.Z);
      ASSERT(szY==szZ); (void)szZ;
      ASSERT(context->dataSet.Y.fmt==context->dataSet.Z.fmt);
      memcpy(Z,Y,szY);
      Y=Z;
  }
  te_errh_resetStates( context );
  fxn( Y, X, Z, N );
  te_errh_verifyStates( context, -1 );

} /* te_processFxn_vYvXvZ() */

/* Apply the target function to the test case data set:
 * vector X (in), vector Y (in), scalar Z (out) */
void te_processFxn_vXvYsZ( tTestEngContext * context )
{
  typedef int16_t         tFxn_fr16  ( const int16_t         * x, const int16_t         * y, int N );
  typedef int32_t         tFxn_fr32  ( const int32_t         * x, const int32_t         * y, int N );
  typedef float16_t       tFxn_fl16  ( const float16_t       * x, const float16_t       * y, int N );
  typedef float32_t       tFxn_fl32  ( const float32_t       * x, const float32_t       * y, int N );
  typedef float64_t       tFxn_fl64  ( const float64_t       * x, const float64_t       * y, int N );
  typedef complex_fract16 tFxn_fr16c ( const complex_fract16 * x, const complex_fract16 * y, int N );
  typedef complex_fract32 tFxn_fr32c ( const complex_fract32 * x, const complex_fract32 * y, int N );
  typedef complex_float16 tFxn_fl16c ( const complex_float16 * x, const complex_float16 * y, int N );
  typedef complex_float   tFxn_fl32c ( const complex_float   * x, const complex_float   * y, int N );
  typedef complex_double  tFxn_fl64c ( const complex_double  * x, const complex_double  * y, int N );

  tTestEngTarget   fxn;
  void *X, *Y, *Z;
  int N;

  ASSERT( context && context->target.fut );
  te_vReportStd(context);

  X = vecGetElem( &context->dataSet.X, 0 );
  Y = vecGetElem( &context->dataSet.Y, 0 );
  Z = vecGetElem( &context->dataSet.Z, 0 );

  N = context->args.dim[0];

  fxn = context->target.fut;

  switch ( context->desc->fmt )
  {
  case FMT_REAL|FMT_FRACT16: *(int16_t        *)Z = ( (tFxn_fr16 *)fxn )( ( const int16_t         *)X, ( const int16_t         *)Y, N ); break;
  case FMT_REAL|FMT_FRACT32: *(int32_t        *)Z = ( (tFxn_fr32 *)fxn )( ( const int32_t         *)X, ( const int32_t         *)Y, N ); break;
  case FMT_REAL|FMT_FLOAT16: *(float16_t      *)Z = ( (tFxn_fl16 *)fxn )( ( const float16_t       *)X, ( const float16_t       *)Y, N ); break;
  case FMT_REAL|FMT_FLOAT32: *(float32_t      *)Z = ( (tFxn_fl32 *)fxn )( ( const float32_t       *)X, ( const float32_t       *)Y, N ); break;
  case FMT_REAL|FMT_FLOAT64: *(float64_t      *)Z = ( (tFxn_fl64 *)fxn )( ( const float64_t       *)X, ( const float64_t       *)Y, N ); break;
  case FMT_CPLX|FMT_FRACT16: *(complex_fract16*)Z = ( (tFxn_fr16c*)fxn )( ( const complex_fract16 *)X, ( const complex_fract16 *)Y, N ); break;
  case FMT_CPLX|FMT_FRACT32: *(complex_fract32*)Z = ( (tFxn_fr32c*)fxn )( ( const complex_fract32 *)X, ( const complex_fract32 *)Y, N ); break;
  case FMT_CPLX|FMT_FLOAT16: *(complex_float16*)Z = ( (tFxn_fl16c*)fxn )( ( const complex_float16 *)X, ( const complex_float16 *)Y, N ); break;
  case FMT_CPLX|FMT_FLOAT32: *(complex_float  *)Z = ( (tFxn_fl32c*)fxn )( ( const complex_float   *)X, ( const complex_float   *)Y, N ); break;
  case FMT_CPLX|FMT_FLOAT64: *(complex_double *)Z = ( (tFxn_fl64c*)fxn )( ( const complex_double  *)X, ( const complex_double  *)Y, N ); break;
  default: ASSERT( 0 );
  }

} /* te_processFxn_vXvYsZ() */

/* Apply the target function to the test case data set:
 * vector Z (out), vector Y (in), vector X (in) */
void te_processFxn_vZvYvX( tTestEngContext * context )
{
  typedef void tFxn( const void * Y, const void * X, void * Z, int N );

  tFxn *fxn;
  void *X, *Y, *Z;
  int N;

  ASSERT( context && context->target.fut );

  X = vecGetElem( &context->dataSet.X, 0 );
  Y = vecGetElem( &context->dataSet.Y, 0 );
  Z = vecGetElem( &context->dataSet.Z, 0 );

  N   = context->args.dim[0];
  fxn = (tFxn*)context->target.fut;
  /* 
    test for overlapping X and Z for each 1st of 4 test if overlapping 
    is enabled in extraParams
  */
  if (((context->args.caseNum & 3)==1) && (context->desc->extraParam & OVLP_XZ))
  {
      size_t szX,szZ;
      szX=vecGetSize(&context->dataSet.X);
      szZ=vecGetSize(&context->dataSet.Z);
      (void)szZ;
      ASSERT(szX==szZ);
      ASSERT(context->dataSet.X.fmt==context->dataSet.Z.fmt);
      memcpy(Z,X,szX);
      X=Z;
  }
  /* 
    test for overlapping Y and Z for each 2nd of 4 test if overlapping 
    is enabled in extraParams
  */
  if (((context->args.caseNum & 3)==2) && (context->desc->extraParam & OVLP_YZ))
  {
      size_t szY,szZ;
      szY=vecGetSize(&context->dataSet.Y);
      szZ=vecGetSize(&context->dataSet.Z);
      (void)szZ;
      ASSERT(szY==szZ);
      ASSERT(context->dataSet.Y.fmt==context->dataSet.Z.fmt);
      memcpy(Z,Y,szY);
      Y=Z;
  }
  te_errh_resetStates( context );
  te_vReportStd(context);
  fxn( Z, Y, X, N );
  te_errh_verifyStates( context, -1 );

} /* te_processFxn_vZvYvX() */

/* Apply the target function to the test case data set:
 * vector X (in), vector Z (out) */
void te_processFxn_vXvZ( tTestEngContext * context )
{
  typedef void tFxn( const void * X, void * Z, int N );

  tFxn *fxn;
  void *X, *Z;
  int N;

  ASSERT( context && context->target.fut );
  te_vReportStd(context);

  X = vecGetElem( &context->dataSet.X, 0 );
  Z = vecGetElem( &context->dataSet.Z, 0 );

  N = context->args.dim[0];

  /* 
    test for overlapping X and Z for each 1st of 2 test if overlapping 
    is enabled in extraParams
   */
  if ((context->args.caseNum & 1) && (context->desc->extraParam & OVLP_XZ))
  {
      size_t szX,szZ;
      szX=vecGetSize(&context->dataSet.X);
      szZ=vecGetSize(&context->dataSet.Z);
      ASSERT(szX==szZ); (void)szZ;
      ASSERT(context->dataSet.X.fmt==context->dataSet.Z.fmt);
      memcpy(Z,X,szX);
      X=Z;
  }
  fxn = (tFxn*)context->target.fut;

  te_errh_resetStates( context );
  fxn( X, Z, N );
  te_errh_verifyStates( context, -1 );

} /* te_processFxn_vXvZ() */

/* Apply the target function to the test case data set:
 * vector Z (out), vector X (in) */
void te_processFxn_vZvX( tTestEngContext * context )
{
  typedef void tFxn( void * Z, const void * X, int N );

  tFxn *fxn;
  void *X, *Z;
  int N;

  ASSERT( context && context->target.fut );
  te_vReportStd(context);

  X = vecGetElem( &context->dataSet.X, 0 );
  Z = vecGetElem( &context->dataSet.Z, 0 );

  if ( context->desc->dimNum == TE_DIM_NUM_1 )
  {
    N = context->args.dim[0];
  }
  else
  {
    N = context->args.dim[1];
  }

  fxn = (tFxn*)context->target.fut;
  /* 
    test for overlapping X and Z for each 1st of 2 test if overlapping 
    is enabled in extraParams
  */
  if (((context->args.caseNum & 1)==1) && (context->desc->extraParam & OVLP_XZ))
  {
      size_t szX,szZ;
      szX=vecGetSize(&context->dataSet.X);
      szZ=vecGetSize(&context->dataSet.Z);
      ASSERT(szX==szZ); (void)szZ;
      ASSERT(context->dataSet.X.fmt==context->dataSet.Z.fmt);
      memcpy(Z,X,szX);
      X=Z;
  }

  te_errh_resetStates(context);
  fxn(Z, X, N);
  te_errh_verifyStates(context, -1);

} /* te_processFxn_vZvX() */

/* Apply the target function to the test case data set:
 * vector X (in), scalar Z (out) */
void te_processFxn_vXsZ( tTestEngContext * context )
{
  typedef int16_t         tFxn_fr16  ( const int16_t         * x, int N );
  typedef int32_t         tFxn_fr32  ( const int32_t         * x, int N );
  typedef float16_t       tFxn_fl16  ( const float16_t       * x, int N );
  typedef float32_t       tFxn_fl32  ( const float32_t       * x, int N );
  typedef float64_t       tFxn_fl64  ( const float64_t       * x, int N );
  typedef complex_fract16 tFxn_fr16c ( const complex_fract16 * x, int N );
  typedef complex_fract32 tFxn_fr32c ( const complex_fract32 * x, int N );
  typedef complex_float16 tFxn_fl16c ( const complex_float16 * x, int N );
  typedef complex_float   tFxn_fl32c ( const complex_float   * x, int N );
  typedef complex_double  tFxn_fl64c ( const complex_double  * x, int N );

  tTestEngTarget   fxn;
  void *X, *Z;
  int N;

  ASSERT( context && context->target.fut );
  te_vReportStd(context);

  X = vecGetElem( &context->dataSet.X, 0 );
  Z = vecGetElem( &context->dataSet.Z, 0 );

  N = context->args.dim[0];

  fxn = context->target.fut;

  switch ( context->desc->fmt )
  {
  case FMT_REAL|FMT_FRACT16: *(int16_t        *)Z = ( (tFxn_fr16 *)fxn )( (const int16_t         *)X, N ); break;
  case FMT_REAL|FMT_FRACT32: *(int32_t        *)Z = ( (tFxn_fr32 *)fxn )( (const int32_t         *)X, N ); break;
  case FMT_REAL|FMT_FLOAT16: *(float16_t      *)Z = ( (tFxn_fl16 *)fxn )( (const float16_t       *)X, N ); break;
  case FMT_REAL|FMT_FLOAT32: *(float32_t      *)Z = ( (tFxn_fl32 *)fxn )( (const float32_t       *)X, N ); break;
  case FMT_REAL|FMT_FLOAT64: *(float64_t      *)Z = ( (tFxn_fl64 *)fxn )( (const float64_t       *)X, N ); break;
  case FMT_CPLX|FMT_FRACT16: *(complex_fract16*)Z = ( (tFxn_fr16c*)fxn )( (const complex_fract16 *)X, N ); break;
  case FMT_CPLX|FMT_FRACT32: *(complex_fract32*)Z = ( (tFxn_fr32c*)fxn )( (const complex_fract32 *)X, N ); break;
  case FMT_CPLX|FMT_FLOAT16: *(complex_float16*)Z = ( (tFxn_fl16c*)fxn )( (const complex_float16 *)X, N ); break;
  case FMT_CPLX|FMT_FLOAT32: *(complex_float  *)Z = ( (tFxn_fl32c*)fxn )( (const complex_float   *)X, N ); break;
  case FMT_CPLX|FMT_FLOAT64: *(complex_double *)Z = ( (tFxn_fl64c*)fxn )( (const complex_double  *)X, N ); break;
  default: ASSERT( 0 );
  }

} /* te_processFxn_vXsZ() */

/* Apply the target function to the test case data set:
 * vector X (in), scalar int32 Z (out) */
void te_processFxn_vXsZ32( tTestEngContext * context )
{
  typedef int32_t tFxn( const void * X, int N );

  tFxn *fxn;
  void *X, *Z;
  int N;

  ASSERT( context && context->target.fut );
  te_vReportStd(context);

  X = vecGetElem( &context->dataSet.X, 0 );
  Z = vecGetElem( &context->dataSet.Z, 0 );

  N = context->args.dim[0];

  fxn = (tFxn*)context->target.fut;

  *(int32_t*)Z = fxn( X, N );

} /* te_processFxn_vXsZ32() */

/* Apply the target function to the test case data set:
 * vector X (in), vector Z (out), vector W (out) */
void te_processFxn_vXvZvW( tTestEngContext * context )
{
  typedef void tFxn( const void * X, void * Z, void * W, int N );

  tFxn *fxn;
  void *X, *Z, *W;
  int N;

  ASSERT( context && context->target.fut );
  te_vReportStd(context);

  X = vecGetElem( &context->dataSet.X, 0 );
  Z = vecGetElem( &context->dataSet.Z, 0 );
  W = vecGetElem( &context->dataSet.W, 0 );

  N = context->args.dim[0];

  fxn = (tFxn*)context->target.fut;

  fxn( X, Z, W, N );

} /* te_processFxn_vXvZvW() */

/* Apply the target function to the test case data set:
 * scalar int32 Z (out), vector X (in), vector Y (in), real scalar U (in) */
void te_processFxn_sZ32_vXvYsUr( tTestEngContext * context ) 
{
  void *X, *Y, *U, *Z;
  int N;

  ASSERT( context && context->target.fut );
  te_vReportStd(context);
  te_vReportStd(context);

  X = vecGetElem( &context->dataSet.X, 0 );
  Y = vecGetElem( &context->dataSet.Y, 0 );
  U = vecGetElem( &context->dataSet.U, 0 );
  Z = vecGetElem( &context->dataSet.Z, 0 );

  N = MAX( 0, context->args.dim[0] );

  if ( context->desc->fmt & FMT_CPLX )
  {
    //typedef void tFxn( void * Z, const void * X, const void * Y, int u, int N );
    //tFxn * fxn = (tFxn*)context->target.fut;
    //fxn( Z, X, Y, (int)*(int32_t*)U, N );
    //!!!!!!!!!!!!!!!!!!!!!!!!!!!
    typedef complex_fract32 tFxn(const void * X, const void * Y, int u, int N);
    tFxn * fxn = (tFxn*)context->target.fut;
    *(complex_fract32*)Z = fxn(X, Y, (int)*(int32_t*)U, N);
  }
  else
  {
    typedef int32_t tFxn( const void * X, const void * Y, int u, int N );
    tFxn * fxn = (tFxn*)context->target.fut;
    *(int32_t*)Z = fxn( X, Y, (int)*(int32_t*)U, N );
  }

} /* te_processFxn_sZ32_vXvYsUr() */

/* Apply the target function to the test case data set, streaming variant:
 * vector Z (out), vector X (in), vector Y (in), real scalar U (in) */
void te_processFxn_vZ_vXvYsUr( tTestEngContext * context )
{
  typedef void tFxn( void * Z, const void * X, const void * Y, int u, int N );

  tFxn * fxn;
  void *X, *Y, *U, *Z;
  int N;

  ASSERT( context && context->target.fut );
  te_vReportStd(context);

  X = vecGetElem( &context->dataSet.X, 0 );
  Y = vecGetElem( &context->dataSet.Y, 0 );
  U = vecGetElem( &context->dataSet.U, 0 );
  Z = vecGetElem( &context->dataSet.Z, 0 );

  N = MAX( 0, context->args.dim[0] );

  fxn = (tFxn*)context->target.fut;

  fxn( Z, X, Y, (int)*(int32_t*)U, N );

} /* te_processFxn_vZ_vXvYsUr() */

/* Apply the target function to the test case data set, streaming variant:
 * scalar int16 Z (out), vector X (in), real scalar Y (in) */
void te_processFxn_sZ16_vXsYr( tTestEngContext * context )
{
  void *X, *Y, *Z;
  int N;

  ASSERT( context && context->target.fut );
  te_vReportStd(context);

  X = vecGetElem( &context->dataSet.X, 0 );
  Y = vecGetElem( &context->dataSet.Y, 0 );
  Z = vecGetElem( &context->dataSet.Z, 0 );

  N = MAX( 0, context->args.dim[0] );

  if ( context->desc->fmt & FMT_CPLX )
  {
    typedef void tFxn( void * Z, const void * X, int y, int N );
    tFxn * fxn = (tFxn*)context->target.fut;
    fxn( Z, X, (int)*(int32_t*)Y, N );
  }
  else
  {
    typedef int16_t tFxn( const void * X, int y, int N );
    tFxn * fxn = (tFxn*)context->target.fut;
    *(int16_t*)Z = fxn( X, (int)*(int32_t*)Y, N );
  }

} /* te_processFxn_sZ16_vXsYr() */

/* Apply the target function to the test case data set, streaming variant:
 * scalar int32 Z (out), vector X (in), real scalar Y (in)  */
void te_processFxn_sZ32_vXsYr( tTestEngContext * context )
{
  void *X, *Y, *Z;
  int N;

  ASSERT( context && context->target.fut );
  te_vReportStd(context);

  X = vecGetElem( &context->dataSet.X, 0 );
  Y = vecGetElem( &context->dataSet.Y, 0 );
  Z = vecGetElem( &context->dataSet.Z, 0 );

  N = MAX( 0, context->args.dim[0] );

  if ( context->desc->fmt & FMT_CPLX )
  {
    typedef void tFxn( void * Z, const void * X, int y, int N );
    tFxn * fxn = (tFxn*)context->target.fut;
    fxn( Z, X, (int)*(int32_t*)Y, N );
  }
  else
  {
    typedef int32_t tFxn( const void * X, int y, int N );
    tFxn * fxn = (tFxn*)context->target.fut;
    *(int32_t*)Z = fxn( X, (int)*(int32_t*)Y, N );
  }

} /* te_processFxn_sZ32_vXsYr() */

/* Apply the target function to the test case data set, streaming variant:
 * vector Z (out), vector X (in), scalar Y (in) */
void te_processFxn_vZ_vXsY( tTestEngContext * context )
{
  typedef void tFxn_i16( void * Z, const void * X, int16_t   y, int N );
  typedef void tFxn_i32( void * Z, const void * X, int32_t   y, int N );
  typedef void tFxn_f32( void * Z, const void * X, float32_t y, int N );

  void *X, *Y, *Z;
  int N;

  ASSERT( context && context->target.fut );
  te_vReportStd(context);

  X = vecGetElem( &context->dataSet.X, 0 );
  Y = vecGetElem( &context->dataSet.Y, 0 );
  Z = vecGetElem( &context->dataSet.Z, 0 );

  N = MAX( 0, context->args.dim[0] );
  /* 
    test for overlapping X and Z for each odd test if overlapping 
    is enabled in extraParams
  */
  if ((context->args.caseNum & 1) && (context->desc->extraParam & OVLP_XZ))
  {
      size_t szX,szZ;
      szX=vecGetSize(&context->dataSet.X);
      szZ=vecGetSize(&context->dataSet.Z);
      ASSERT(szX==szZ); (void)szZ;
      ASSERT(context->dataSet.X.fmt==context->dataSet.Z.fmt);
      memcpy(Z,X,szX);
      X=Z;
  }
  switch(context->desc->fmt & FMT_DTYPE_MASK)
  {
  case FMT_FRACT16:
      ((tFxn_i16*)context->target.fut)( Z, X, *(int16_t*)Y, N );
      break;
  case FMT_FRACT32:
      ((tFxn_i32*)context->target.fut)( Z, X, ((int32_t*)Y)[0], N );
      break;
  case FMT_FLOAT32:
      ((tFxn_f32*)context->target.fut)( Z, X, *(float32_t*)Y, N );
      break;
  default:
      ASSERT(0);    // not implemented !
  }

} /* te_processFxn_vZ_vXsY() */

/* Apply the target function to the test case data set, streaming variant:
 * scalar int Z (out), vector X (in) */
void te_processFxn_sZi_vX( tTestEngContext * context )
{
  typedef int tFxn( const void * X, int N );

  tFxn * fxn;
  void *X, *Z;
  int N;

  ASSERT( context && context->target.fut );
  te_vReportStd(context);

  X = vecGetElem( &context->dataSet.X, 0 );
  Z = vecGetElem( &context->dataSet.Z, 0 );

  N = MAX( 0, context->args.dim[0] );

  fxn = (tFxn*)context->target.fut;

  switch ( context->desc->fmt & FMT_DTYPE_MASK )
  {
  case FMT_FRACT16: *(int16_t*)Z = (int16_t)fxn( X, N ); break;
  case FMT_FRACT32: *(int32_t*)Z = (int16_t)fxn( X, N ); break;
  default: ASSERT( 0 );
  }

} /* te_processFxn_sZi_vX() */

/* Apply the target function to the test case data set, scalar target functions:
 * vector Z (out), vector X (in) */
void te_processFxn_scl_vZvX( tTestEngContext * context )
{
  #define APPLY_FXN( Fxn, typeX, X, typeZ, Z, N )             \
      { int n;                                                \
        for ( n=0; n<N; n++ ) {                               \
          typedef typeZ typeFxn( typeX );                     \
          *(typeZ*)(Z) = ( (typeFxn*)(Fxn) )( *(typeX*)(X) ); \
          (X) = (typeX*)(X) + 1; (Z) = (typeZ*)(Z) + 1;       \
        } }

  tTestEngTarget fxn;
  void *X, *Z;
  int fmtX, fmtZ;
  int N;

  ASSERT( context && context->target.fut );
  te_vReportStd(context);

  fmtX = context->dataSet.X.fmt;
  fmtZ = context->dataSet.Z.fmt;

  ASSERT( 0 == ( (fmtX|fmtZ) & FMT_CPLX ) );

  X = vecGetElem( &context->dataSet.X, 0 );
  Z = vecGetElem( &context->dataSet.Z, 0 );

  N = context->args.dim[0];

  fxn = context->target.fut;

  switch ( ((fmtX & FMT_DTYPE_MASK)<<8)|(fmtZ & FMT_DTYPE_MASK) )
  {
  case (FMT_FRACT16<<8)|FMT_FRACT16: APPLY_FXN( fxn, int16_t, X, int16_t, Z, N ); break;
  case (FMT_FRACT16<<8)|FMT_FRACT32: APPLY_FXN( fxn, int16_t, X, int32_t, Z, N ); break;
  case (FMT_FRACT32<<8)|FMT_FRACT16: APPLY_FXN( fxn, int32_t, X, int16_t, Z, N ); break;
  case (FMT_FRACT32<<8)|FMT_FRACT32: APPLY_FXN( fxn, int32_t, X, int32_t, Z, N ); break;
  default: ASSERT(0);
  }

} /* te_processFxn_scl_vZvX() */

/* Apply the target function to the test case data set, streaming variant:
 * vector X (in), scalar int32 Y (in), scalar int64 Z (out) */
void te_processFxn_vXsY32sZ64( tTestEngContext * context )
{
  typedef int64_t tFxn( const void * X, int y, int N );

  tFxn *fxn;
  void *X, *Z;
  int32_t y;
  int N;

  ASSERT( context && context->target.fut );

  X = vecGetElem( &context->dataSet.X, 0 );
  Z = vecGetElem( &context->dataSet.Z, 0 );

  y = *vecGetElem_i32( &context->dataSet.Y, 0 );

  N   = context->args.dim[0];
  fxn = (tFxn*)context->target.fut;
  te_vReportStd(context);

  *(int64_t*)Z = fxn( X, y, N );

} /* te_processFxn_vXsY32sZ64() */

/* Apply the target function to the test case data set, streaming variant:
 * scalar int32 Z (out), vector X (in), vector Y (in) */
void te_processFxn_sZ32vXvY( tTestEngContext * context )
{
  typedef int32_t tFxn( const void * X, const void * Y, int N );

  tFxn *fxn;
  void *X, *Y, *Z;
  int N;

  ASSERT( context && context->target.fut );

  X = vecGetElem( &context->dataSet.X, 0 );
  Y = vecGetElem( &context->dataSet.Y, 0 );
  Z = vecGetElem( &context->dataSet.Z, 0 );

  N   = context->args.dim[0];
  fxn = (tFxn*)context->target.fut;
  te_vReportStd(context);

  *(int32_t*)Z = fxn( X, Y, N );

} /* te_processFxn_sZ32vXvY() */

/* Apply the target function to the test case data set, streaming variant:
 * scalar int64 Z (out), vector X (in), vector Y (in) */
void te_processFxn_sZ64vXvY( tTestEngContext * context )
{
  typedef int64_t tFxn( const void * X, const void * Y, int N );

  tFxn *fxn;
  void *X, *Y, *Z;
  int N;

  ASSERT( context && context->target.fut );

  X = vecGetElem( &context->dataSet.X, 0 );
  Y = vecGetElem( &context->dataSet.Y, 0 );
  Z = vecGetElem( &context->dataSet.Z, 0 );

  N   = context->args.dim[0];
  fxn = (tFxn*)context->target.fut;
  te_vReportStd(context);

  *(int64_t*)Z = fxn( X, Y, N );

} /* te_processFxn_sZ32vXvY() */

/* Apply the target function to the test case data set, streaming variant:
 * vector X (in), vector Y (in), vector Z (out) */
void te_processFxn_s_vXvYvZ( tTestEngContext * context )
{
  typedef void tFxn( const void * X, const void * Y, void * Z, int N, int L );

  tFxn *fxn;
  void *X, *Y, *Z;
  int N, L;

  ASSERT( context && context->target.fut );
  te_vReportStd(context);

  X = vecGetElem( &context->dataSet.X, 0 );
  Y = vecGetElem( &context->dataSet.Y, 0 );
  Z = vecGetElem( &context->dataSet.Z, 0 );

  N = context->args.dim[0];
  L = context->args.dim[1];

  fxn = (tFxn*)context->target.fut;

  fxn( X, Y, Z, N, L );

} /* te_processFxn_s_vXvYvZ() */

/* Apply the target function to the test case data set, streaming variant:
 * vector X (in), vector Z (out) */
void te_processFxn_s_vXvZ( tTestEngContext * context )
{
  typedef void tFxn( const void * X, void * Z, int N, int L );

  tFxn *fxn;
  void *X, *Z;
  int N, L;

  ASSERT( context && context->target.fut );
  te_vReportStd(context);

  X = vecGetElem( &context->dataSet.X, 0 );
  Z = vecGetElem( &context->dataSet.Z, 0 );

  N = context->args.dim[0];
  L = context->args.dim[1];

  fxn = (tFxn*)context->target.fut;

  fxn( X, Z, N, L );

} /* te_processFxn_s_vXvZ() */
