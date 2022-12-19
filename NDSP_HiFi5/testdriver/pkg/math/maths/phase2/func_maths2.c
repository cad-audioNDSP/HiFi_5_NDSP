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
 * Test procedures for vector mathematics
 */

/* Cross-platform data type definitions. */
#include "types.h"
/* Test environemnt convfiguration. */
#include "config.h"
#include "packages.h"
/* DSP Library API: arithmetic and logic functions on data vectors. */
#include LIBRARY_HEADER(math)
/* Test engine API. */
#include "../../common/testeng_math.h"

/* Test executive function. Performs the specified test on a brief or full or sanity version
 * of the designated SEQ-file. Return the test result (non-zero if passed). */
static int testExec( tTestEngTarget targetFxn, const char * seqName,
                     int errhExtendedTest, int isFull, int isVerbose, int breakOnError, int testBitexactness );


static void processFxn_scl_vZvXvY( tTestEngContext * context )
{
  typedef float32_t tFxn_fl32( float32_t, float32_t );
  typedef float64_t tFxn_fl64( float64_t, float64_t );
  typedef fract16   tFxn_fr16( fract16  , fract16   );
  typedef fract32   tFxn_fr32( fract32  , fract32   );

  #define CALL_FXN( typeFxn, Fxn, typeXYZ, X, Y, Z, n )                         \
      { te_errh_resetStates( context );                                         \
        *(typeXYZ*)(Z) = ( (typeFxn*)(Fxn) )( *(typeXYZ*)(X), *(typeXYZ*)(Y) ); \
        te_errh_verifyStates( context, n );                                     \
        (X) = (typeXYZ*)(X) + 1; (Y) = (typeXYZ*)(Y) + 1; (Z) = (typeXYZ*)(Z) + 1; }

  tTestEngTarget  fxn;
  void *X, *Y, *Z;
  int n, N;

  X = vecGetElem( &context->dataSet.X, 0 );
  Y = vecGetElem( &context->dataSet.Y, 0 );
  Z = vecGetElem( &context->dataSet.Z, 0 );

  fxn = context->target.fut;
  N   = context->args.dim[0];

  te_vReportStd(context);

  switch ( context->desc->fmt & ~FMT_CPLX )
  {
  case FMT_FLOAT32: for ( n=0; n<N; n++ ) CALL_FXN( tFxn_fl32, fxn, float32_t, X, Y, Z, n ); break;
  case FMT_FLOAT64: for ( n=0; n<N; n++ ) CALL_FXN( tFxn_fl64, fxn, float64_t, X, Y, Z, n ); break;
  case FMT_FRACT16: for ( n=0; n<N; n++ ) CALL_FXN( tFxn_fr16, fxn, fract16  , X, Y, Z, n ); break;
  case FMT_FRACT32: for ( n=0; n<N; n++ ) CALL_FXN( tFxn_fr32, fxn, fract32  , X, Y, Z, n ); break;
  default: ASSERT( 0 );
  }
  #undef CALL_FXN
}

/* Perform all tests for math API functions. */
int func_maths2(int isFull, int isVerbose, int breakOnError)
{
  int res = 1;

  #define DO_TEST( fxn, seqFile, extraFlags )                                                         \
    if ( res || !breakOnError ) res &= ( 0 != testExec( (tTestEngTarget)(fxn), "maths2/" seqFile,       \
                                                        (extraFlags) | TE_ERRH_EXTENDED_TEST_DISABLE, \
                                                        isFull, isVerbose, breakOnError, 0 ) )
  /* Extended variant runs each SEQ-file twice, with extended Error Handling check in the second
   * invocation. */
  #define DO_TEST_EXT( fxn, seqFile, extraFlags )                                                     \
    { if ( res || !breakOnError ) res &= ( 0 != testExec(                                             \
                                                      (tTestEngTarget)(fxn),                          \
                                                      "maths2/" seqFile,                                      \
                                                      (extraFlags) | TE_ERRH_EXTENDED_TEST_DISABLE,   \
                                                      isFull, isVerbose, breakOnError, 0 ) );         \
      if ( res || !breakOnError ) res &= ( 0 != testExec(                                             \
                                                      (tTestEngTarget)(fxn),                          \
                                                      "maths2/" seqFile,                                      \
                                                      (extraFlags) | TE_ERRH_EXTENDED_TEST_ENABLE,    \
                                                      isFull, isVerbose, breakOnError, 0 ) ); }
  DO_TEST    ( &scl_int2float     , "scl_int2float.seq"       , 0 );
  DO_TEST    ( &scl_float2int     , "scl_float2int.seq"       , 0 );
  DO_TEST_EXT( &scl_sinef         , "scl_sinef.seq"           , 0 );
  DO_TEST_EXT( &scl_cosinef       , "scl_cosinef.seq"         , 0 );
  DO_TEST_EXT( &scl_tanf          , "scl_tanf.seq"            , 0 );
  DO_TEST_EXT( &scl_log2f         , "scl_log2f.seq"           , 0 );
  DO_TEST_EXT( &scl_lognf         , "scl_lognf.seq"           , 0 );
  DO_TEST_EXT( &scl_log10f        , "scl_log10f.seq"          , 0 );
  DO_TEST_EXT( &scl_antilog2f     , "scl_antilog2f.seq"       , 0 );
  DO_TEST_EXT( &scl_antilognf     , "scl_antilognf.seq"       , 0 );
  DO_TEST_EXT( &scl_antilog10f    , "scl_antilog10f.seq"      , 0 );
  DO_TEST_EXT( &scl_powf          , "scl_powf.seq"            , 0 );
  DO_TEST_EXT( &scl_atanf         , "scl_atanf.seq"           , 0 );
  DO_TEST_EXT( &scl_atan2f        , "scl_atan2f.seq"          , 0 );
  DO_TEST_EXT( &scl_tanhf         , "scl_tanhf.seq"           , 0 );
  DO_TEST    ( &scl_sigmoidf      , "scl_sigmoidf.seq"        , 0 );
  DO_TEST    ( &scl_reluf         , "scl_reluf.seq"           , 0 );
  return (res);
}

/* Test executive function. Performs the specified test on a brief or full version
 * of the designated SEQ-file. Return the test result (non-zero if passed). */
static int testExec( tTestEngTarget   targetFxn, const char * seqName, 
              int errhExtendedTest, int isFull, int isVerbose, int breakOnError, int testBitexactness )
{
  #define MAX_FUNC_NUM   16
  /* Initializer for a function pointer array, appends NULL to a sequence of pointers. */
  #define FUNC_LIST(...) { __VA_ARGS__, NULL }
  /* Initializer for a test description structure. */
  #define TEST_DESC( fmt, dimNum, align, loadFxn, procFxn ) { (fmt),0,NULL,(dimNum),(align),NULL,NULL,(loadFxn),(procFxn) }

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
      FUNC_LIST( (tTestEngTarget)&scl_int2float),
      TEST_DESC(FMT_FLOAT32, TE_DIM_NUM_2, TE_ALIGN_NO, &te_loadFxn_vX32sY32vZ, &te_math_processFxn_scl_vX32sY32vZ) },
    {
      FUNC_LIST( (tTestEngTarget)&scl_float2int),
      TEST_DESC(FMT_FLOAT32, TE_DIM_NUM_2, TE_ALIGN_NO, &te_loadFxn_vXsY32vZ32, &te_math_processFxn_scl_vXsY32vZ32) },
    {
      FUNC_LIST( (tTestEngTarget)&scl_sinef,(tTestEngTarget)&scl_cosinef,(tTestEngTarget)&scl_tanf,
                 (tTestEngTarget)&scl_log2f,(tTestEngTarget)&scl_lognf,(tTestEngTarget)&scl_log10f,
                 (tTestEngTarget)&scl_atanf,(tTestEngTarget)&scl_antilog2f,(tTestEngTarget)&scl_antilognf,
                 (tTestEngTarget)&scl_antilog10f),
      TEST_DESC( FMT_REAL|FMT_FLOAT32, TE_DIM_NUM_1, TE_ALIGN_NO, &te_loadFxn_vXvZ, &te_math_processFxn_scl_vXvZ ) },
    {
      FUNC_LIST( (tTestEngTarget)&scl_atan2f ),
      TEST_DESC(FMT_FLOAT32, TE_DIM_NUM_1, TE_ALIGN_NO, &te_loadFxn_vXvYvZ, &te_math_processFxn_scl_atan2) },
    {
      FUNC_LIST( (tTestEngTarget)&scl_tanhf, (tTestEngTarget)&scl_sigmoidf),
      TEST_DESC( FMT_FLOAT32, TE_DIM_NUM_1, TE_ALIGN_NO, &te_loadFxn_vXvZ, &te_math_processFxn_scl_vXvZ ) },
    {
      FUNC_LIST( (tTestEngTarget)&scl_reluf),
      TEST_DESC( FMT_FLOAT32, TE_DIM_NUM_1, TE_ALIGN_NO, &te_loadFxn_vXvYvZ, &processFxn_scl_vZvXvY ) },
    {
      FUNC_LIST( (tTestEngTarget)&scl_powf),
      TEST_DESC(FMT_REAL | FMT_FLOAT32, TE_DIM_NUM_1, TE_ALIGN_NO, &te_loadFxn_vXvYvZ, &te_math_processFxn_scl_vXvYvZ) },
    { 
      FUNC_LIST( NULL ), TEST_DESC(  0, 0, 0, NULL, NULL ) } /* End of table */
  };

  {
    int tblIx, funcIx;

    for ( tblIx=0; tblIx<(int)(sizeof(testDefTbl)/sizeof(testDefTbl[0])); tblIx++ )
    {
      for ( funcIx=0; funcIx<MAX_FUNC_NUM; funcIx++ )
      {
        if ( targetFxn == testDefTbl[tblIx].funcList[funcIx] )
        {
          tTestEngDesc testDesc = testDefTbl[tblIx].testDesc;
          testDesc.extraParam = (uint32_t)errhExtendedTest;

          return ( TestEngRun( targetFxn, &testDesc, 
                               seqName, isFull, 
                               isVerbose, breakOnError, testBitexactness ) );
        }
      }
    }

    ASSERT( !"Test not defined" );
    return (0);
  }
  return te_Exec(testDefTbl, sizeof(testDefTbl) / sizeof(testDefTbl[0]), MAX_FUNC_NUM, targetFxn, seqName, isFull, isVerbose, breakOnError);
}
