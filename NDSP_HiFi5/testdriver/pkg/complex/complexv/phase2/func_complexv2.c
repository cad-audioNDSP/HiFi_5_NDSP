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

#include <math.h>

/* Cross-platform data type definitions. */
#include "types.h"
/* Test environemnt convfiguration. */
#include "config.h"
#include "packages.h"
/* Library API */
#include LIBRARY_HEADER(complex)
/* Test engine API. */
#include "testeng.h"

/* Test executive function. Performs the specified test on a brief or full or sanity version
 * of the designated SEQ-file. Return the test result (non-zero if passed). */
static int testExec( tTestEngTarget targetFxn, const char * seqName,
                     int errhExtendedTest, int isFull, int isVerbose, int breakOnError, int testBitexactness );

/* Apply a function to the test case data set:
 * scalar functions with single argument, e.g. cos(). Input X is complex, output is real*/
static void processFxn_scl_vXcvZ( tTestEngContext * context )
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

} /* processFxn_scl_vXcvZ() */

/* table for checking the bitexactness */
static const tTestEngVecSclTbl bitexactnessTbl[]=
{
    { processFxn_scl_vXcvZ , (te_fun_ptr_t)&scl_complex2mag    , (te_fun_ptr_t)&vec_complex2mag    },
    { processFxn_scl_vXcvZ , (te_fun_ptr_t)&scl_complex2invmag , (te_fun_ptr_t)&vec_complex2invmag },
    {0, NULL, NULL},
};

int func_complexv2(int isFull, int isVerbose, int breakOnError)
{
  int res = 1;

#define DO_TEST_BITEXACTNESS(fxn, seqFile, extraFlags)                                                \
{                                                                                                     \
    if ( res || !breakOnError ) res &= ( 0 != testExec( (tTestEngTarget)(fxn), "complexv2/" seqFile,       \
                                                        (extraFlags) | TE_ERRH_EXTENDED_TEST_DISABLE, \
                                         isFull, isVerbose, breakOnError,0 ) );                       \
    if ( res || !breakOnError ) res &= ( 0 != testExec( (tTestEngTarget)(fxn), "complexv2/" seqFile,       \
                                                        (extraFlags) | TE_ERRH_EXTENDED_TEST_DISABLE, \
                                         isFull, isVerbose, breakOnError,1 ) ); }

  DO_TEST_BITEXACTNESS ( &vec_complex2mag   , "vec_complex2mag.seq"     , 0 );
  DO_TEST_BITEXACTNESS ( &vec_complex2invmag, "vec_complex2invmag.seq"  , 0 );

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
  #define TEST_DESC( fmt,extraParam,dimNum, align, loadFxn, procFxn ) { (fmt),extraParam,bitexactnessTbl,(dimNum),(align),NULL,NULL,(loadFxn),(procFxn) }

  /* vec API test definitions. */
  static const struct 
  {
    tTestEngTarget   funcList[MAX_FUNC_NUM];
    tTestEngDesc     testDesc;
  }
  testDefTbl[] =
  {
    {
      FUNC_LIST( (tTestEngTarget)&vec_complex2mag,(tTestEngTarget)&vec_complex2invmag ),
      TEST_DESC( FMT_REAL|FMT_FLOAT32, 0,TE_DIM_NUM_1, TE_ALIGN_NO, &te_loadFxn_vXcvZ, &te_processFxn_vZvX ) },
    { 
      FUNC_LIST( NULL ), TEST_DESC(  0, 0,0, 0, NULL, NULL ) } /* End of table */
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
