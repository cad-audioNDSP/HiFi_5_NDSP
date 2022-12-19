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
/* DSP Library API: arithmetic and logic functions on data vectors. */
#include LIBRARY_HEADER(math)

/* Test engine API. */
#include "../../common/testeng_math.h"
/* Test executive function. Performs the specified test on a brief or full or sanity version
 * of the designated SEQ-file. Return the test result (non-zero if passed). */
static int testExec( tTestEngTarget targetFxn, const char * seqName,
                     int errhExtendedTest, int isFull, int isVerbose, int breakOnError, int testBitexactness );


/* table for checking the bitexactness */
static const tTestEngVecSclTbl bitexactnessTbl[]=
{
  {0, NULL, NULL},
};

#define DO_TEST( fxn, seqFile, extraFlags )                                                           \
if ( res || !breakOnError ) res &= ( 0 != testExec( (tTestEngTarget)(fxn), "mathvf1/" seqFile,           \
                                                        (extraFlags) | TE_ERRH_EXTENDED_TEST_DISABLE, \
                                                        isFull, isVerbose, breakOnError, 0 ) )

/* Perform all tests for math API functions. */
int func_mathvf1(int isFull, int isVerbose, int breakOnError)
{
  int res = 1;
  DO_TEST( &vec_divide16x16_fast   , "vec_divide16x16_fast.seq" , 0 );
  DO_TEST( &vec_divide32x32_fast   , "vec_divide32x32_fast.seq" , 0 );
#if 0 //HiFi3/3z API
  DO_TEST( &vec_divide24x24_fast   , "vec_divide24x24_fast.seq" , 0 );
#endif
  DO_TEST( &vec_sqrt32x32_fast     , "vec_sqrt32x32_fast.seq"   , 0 );
#if 0 //HiFi3/3z API
  DO_TEST( &vec_sqrt24x24_fast     , "vec_sqrt24x24_fast.seq"   , 0 );
#endif
  DO_TEST( &vec_sine32x32_fast     , "vec_sine32x32_fast.seq"   , 0 );
  DO_TEST( &vec_cosine32x32_fast   , "vec_cosine32x32_fast.seq" , 0 );
#if 0 //HiFi3/3z API
  DO_TEST( &vec_sine24x24_fast     , "vec_sine24x24_fast.seq"   , 0 );
  DO_TEST( &vec_cosine24x24_fast   , "vec_cosine24x24_fast.seq" , 0 );
#endif
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
      FUNC_LIST( (tTestEngTarget)&vec_sine32x32_fast,(tTestEngTarget)&vec_cosine32x32_fast),
      TEST_DESC( FMT_REAL|FMT_FRACT32, OVLP_XZ, TE_DIM_NUM_1, TE_ALIGN_YES, &te_loadFxn_vXvZ, &te_processFxn_vZvX ) },
    {
      FUNC_LIST( (tTestEngTarget)&vec_sqrt32x32_fast),
      TEST_DESC( FMT_REAL|FMT_FRACT32, OVLP_XZ, TE_DIM_NUM_1, TE_ALIGN_YES, &te_loadFxn_vXvZ, &te_processFxn_vZvX ) },
    {
      FUNC_LIST( (tTestEngTarget)&vec_divide16x16_fast),
      TEST_DESC(FMT_REAL | FMT_FRACT16, OVLP_XZ | OVLP_XW | OVLP_YZ | OVLP_YW,TE_DIM_NUM_1, TE_ALIGN_YES, &te_loadFxn_vXvYvZf, &te_math_processFxn_divide16x16) },
    {
      FUNC_LIST( (tTestEngTarget)&vec_divide32x32_fast),
      TEST_DESC(FMT_REAL | FMT_FRACT32, 0,TE_DIM_NUM_1, TE_ALIGN_YES, &te_loadFxn_vXvYvZd, &te_math_processFxn_divide32x32) },
    { 
      FUNC_LIST( NULL ), TEST_DESC(  0, 0, 0, 0, NULL, NULL ) } /* End of table */
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
          testDesc.extraParam |= (uint32_t)errhExtendedTest;

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
