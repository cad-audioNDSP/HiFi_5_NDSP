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
 * Test procedures for arithmetic and logic functions on data vectors.
 */

#include <string.h>

/* Cross-platform data type definitions. */
#include "types.h"
/* Test environment configuration. */
#include "config.h"
#include "packages.h"
/* DSP Library API. */
#include LIBRARY_HEADER(fit)
/* Test engine API. */
#include "testeng.h"
/* Test data vectors tools and SEQ-file reader. */
#include "vectools.h"

#define MIN(a,b)   ( (a)<(b) ? (a) : (b) )
#define MAX(a,b)   ( (a)>(b) ? (a) : (b) )
#define MAX_FUNC_NUM   10
/* Initializer for a function pointer array, appends NULL to a sequence of pointers. */
#define FUNC_LIST(...) { __VA_ARGS__, NULL }
/* Initializer for a test description structure. */
#define TEST_DESC( fmt, dimNum, align, loadFxn, procFxn ) { (fmt),0,NULL,(dimNum),(align),NULL,NULL,(loadFxn),(procFxn) }

/* load data to the vec_polyxxx functions */
static int te_loadFxn_poly( tTestEngContext * context )
{
  int M, N;
  int nElemX, nElemY, nElemZ, res;

  ASSERT( context && context->seqFile );

  M = MAX( 0, context->args.dim[0] );
  N = MAX( 0, context->args.dim[1] );

  nElemX = M;
  nElemY = N;
  nElemZ = M;

  memset( &context->dataSet, 0, sizeof(context->dataSet) );

  /* Allocate data vectors memory. */
  res=1;
  res &= ( 1 == vecsAlloc( context->desc->isAligned, context->desc->fmt,
                         &context->dataSet.X  , nElemX, 0 ) );
  res &= ( 1 == vecsAlloc( context->desc->isAligned, context->desc->fmt,
                         &context->dataSet.Y  , nElemY, 0 ) );
  res &= ( 1 == vecsAlloc( context->desc->isAligned, FMT_INT32,
                         &context->dataSet.Wlo , 1, 0 ) );
  res &= ( 3 == vecsAlloc( context->desc->isAligned, context->desc->fmt,
                         &context->dataSet.Z  , nElemZ,
                         &context->dataSet.Zlo, nElemZ,
                         &context->dataSet.Zhi, nElemZ, 0 ) );
  if ( res )
  {
    /* Load vectors data from the SEQ-file. */
    if ( !( res = seqFileReadVecs( context->seqFile,
                                  &context->dataSet.X,
                                  &context->dataSet.Y,
                                  &context->dataSet.Wlo,
                                  &context->dataSet.Zlo,
                                  &context->dataSet.Zhi, 0 ) ) )
    {
      printf( "te_loadFxn_poly(): failed to read vectors data; \n" );
    }
  }
  else
  {
    printf( "te_loadFxn_poly(): failed to allocate vectors\n");
  }
  /* Free vectors data if failed. */
  if ( !res ) te_freeVectors(context);
  return (res);
} /* te_loadFxn_poly() */

/* process data by the vec_polyxxx functions  */
static void te_processFxn_poly( tTestEngContext * context )
{
  typedef void tFxn_fr32  ( fract32   * z, const fract32   * x, const fract32   *y, int lsh, int N );
  typedef void tFxn_fl32  ( float32_t * z, const float32_t * x, const float32_t *y, int N );

  tTestEngTarget   fxn;
  void *X, *Y, *Z, *W;
  int N;

  ASSERT( context && context->target.fut );
  te_vReportStd(context);
  X = vecGetElem( &context->dataSet.X, 0 );
  Y = vecGetElem( &context->dataSet.Y, 0 );
  W = vecGetElem( &context->dataSet.Wlo, 0 );
  Z = vecGetElem( &context->dataSet.Z, 0 );

  N   = context->args.dim[0];
  fxn = context->target.fut;

  switch ( context->desc->fmt )
  {
  case FMT_REAL|FMT_FRACT32: ( (tFxn_fr32 *)fxn )((fract32  *)Z, (const fract32   *)X,(const fract32   *)Y, *(fract32        *)W,  N ); break;
  case FMT_REAL|FMT_FLOAT32: ( (tFxn_fl32 *)fxn )((float32_t*)Z, (const float32_t *)X,(const float32_t *)Y,  N ); break;
  default: ASSERT( 0 );
  }
} /* te_processFxn_poly() */

/* vec API test definitions. */
static const struct 
{
  tTestEngTarget   funcList[MAX_FUNC_NUM];
  tTestEngDesc     testDesc;
}
testDefTbl[] =
{
  /*
    * Stage 2
    */
  {
    FUNC_LIST( (tTestEngTarget)&vec_poly4f,(tTestEngTarget)&vec_poly8f ),
    TEST_DESC( FMT_REAL|FMT_FLOAT32, TE_DIM_NUM_2, TE_ALIGN_NO, &te_loadFxn_poly, &te_processFxn_poly ) },
  { 
    FUNC_LIST( NULL ), TEST_DESC(  0, 0, 0, NULL, NULL ) } /* End of table */
};

/* Perform all tests for vec API functions. */
int func_pfit2(int isFull, int isVerbose, int breakOnError)
{
  int res = 1;

  #define DO_TEST(fxn, seqFile)                                                                  \
    if ( res || !breakOnError ) res &= ( 0 != te_Exec( testDefTbl,                               \
                                                       sizeof(testDefTbl)/sizeof(testDefTbl[0]), \
                                                       MAX_FUNC_NUM,                             \
                                                       (tTestEngTarget)(fxn), "fit2/" seqFile,    \
                                                       isFull, isVerbose, breakOnError ) )

  /*
   * Stage 2
   */
  DO_TEST( &vec_poly4f        , "vec_poly4f.seq"         );
  DO_TEST( &vec_poly8f        , "vec_poly8f.seq"         );
  return (res);
}
