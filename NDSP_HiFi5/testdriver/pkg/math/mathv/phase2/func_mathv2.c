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
  { te_math_processFxn_scl_vX32sY32vZ  , (te_fun_ptr_t)&scl_int2float       , (te_fun_ptr_t)&vec_int2float       },
  { te_math_processFxn_scl_vXsY32vZ32  , (te_fun_ptr_t)&scl_float2int       , (te_fun_ptr_t)&vec_float2int       },
  { te_math_processFxn_scl_vXvZ        , (te_fun_ptr_t)&scl_sinef           , (te_fun_ptr_t)&vec_sinef           },
  { te_math_processFxn_scl_vXvZ        , (te_fun_ptr_t)&scl_cosinef         , (te_fun_ptr_t)&vec_cosinef         },
  { te_math_processFxn_scl_vXvZ        , (te_fun_ptr_t)&scl_tanf            , (te_fun_ptr_t)&vec_tanf            },
  { te_math_processFxn_scl_vXvZ        , (te_fun_ptr_t)&scl_log2f           , (te_fun_ptr_t)&vec_log2f           },
  { te_math_processFxn_scl_vXvZ        , (te_fun_ptr_t)&scl_lognf           , (te_fun_ptr_t)&vec_lognf           },
  { te_math_processFxn_scl_vXvZ        , (te_fun_ptr_t)&scl_log10f          , (te_fun_ptr_t)&vec_log10f          },
  { te_math_processFxn_scl_vXvZ        , (te_fun_ptr_t)&scl_antilog2f       , (te_fun_ptr_t)&vec_antilog2f       },
  { te_math_processFxn_scl_vXvZ        , (te_fun_ptr_t)&scl_antilognf       , (te_fun_ptr_t)&vec_antilognf       },
  { te_math_processFxn_scl_vXvZ        , (te_fun_ptr_t)&scl_antilog10f      , (te_fun_ptr_t)&vec_antilog10f      },
  { te_math_processFxn_scl_vXvZ        , (te_fun_ptr_t)&scl_atanf           , (te_fun_ptr_t)&vec_atanf           },
  { te_math_processFxn_scl_atan2       , (te_fun_ptr_t)&scl_atan2f          , (te_fun_ptr_t)&vec_atan2f          },
  { te_math_processFxn_scl_vXvZ        , (te_fun_ptr_t)&scl_tanhf           , (te_fun_ptr_t)&vec_tanhf           },
  { te_math_processFxn_scl_vXvZ        , (te_fun_ptr_t)&scl_sigmoidf        , (te_fun_ptr_t)&vec_sigmoidf        },
  { te_math_processFxn_scl_vXvYvZ      , (te_fun_ptr_t)&scl_powf            , (te_fun_ptr_t)&vec_powf            },
  { te_math_processFxn_scl_vXvZ        , (te_fun_ptr_t)&scl_tanh_fp16       , (te_fun_ptr_t)&vec_tanh_fp16       },
  { te_math_processFxn_scl_vXvZ        , (te_fun_ptr_t)&scl_sigmoid_fp16    , (te_fun_ptr_t)&vec_sigmoid_fp16    },
	{0, NULL, NULL},
};

#define DO_TEST( fxn, seqFile, extraFlags )                                                           \
if ( res || !breakOnError ) res &= ( 0 != testExec( (tTestEngTarget)(fxn), "mathv2/" seqFile,           \
                                                        (extraFlags) | TE_ERRH_EXTENDED_TEST_DISABLE, \
                                                        isFull, isVerbose, breakOnError, 0 ) )

#define DO_TEST_BITEXACTNESS(fxn, seqFile, extraFlags )                                               \
{                                                                                                     \
    if ( res || !breakOnError ) res &= ( 0 != testExec( (tTestEngTarget)(fxn), "mathv2/" seqFile,       \
                                            (extraFlags) | TE_ERRH_EXTENDED_TEST_DISABLE,             \
                                         isFull, isVerbose, breakOnError,0 ) );                       \
    if ( res || !breakOnError ) res &= ( 0 != testExec( (tTestEngTarget)(fxn), "mathv2/" seqFile,       \
                                            (extraFlags) | TE_ERRH_EXTENDED_TEST_DISABLE,             \
                                         isFull, isVerbose, breakOnError,1 ) ); }

/* Perform all tests for math API functions. */
int func_mathv2(int isFull, int isVerbose, int breakOnError )
{
	int res = 1;
  DO_TEST_BITEXACTNESS( &vec_int2float     , "vec_int2float.seq"       , 0 );
  DO_TEST_BITEXACTNESS( &vec_float2int     , "vec_float2int.seq"       , 0 );
  DO_TEST_BITEXACTNESS( &vec_sinef         , "vec_sinef.seq"           , 0 );
  DO_TEST_BITEXACTNESS( &vec_cosinef       , "vec_cosinef.seq"         , 0 );
  DO_TEST_BITEXACTNESS( &vec_tanf          , "vec_tanf.seq"            , 0 );
  DO_TEST_BITEXACTNESS( &vec_log2f         , "vec_log2f.seq"           , 0 );
  DO_TEST_BITEXACTNESS( &vec_lognf         , "vec_lognf.seq"           , 0 );
  DO_TEST_BITEXACTNESS( &vec_log10f        , "vec_log10f.seq"          , 0 );
  DO_TEST_BITEXACTNESS( &vec_antilog2f     , "vec_antilog2f.seq"       , 0 );
  DO_TEST_BITEXACTNESS( &vec_antilognf     , "vec_antilognf.seq"       , 0 );
  DO_TEST_BITEXACTNESS( &vec_antilog10f    , "vec_antilog10f.seq"      , 0 );
  DO_TEST_BITEXACTNESS( &vec_powf          , "vec_powf.seq"            , 0 );
  DO_TEST_BITEXACTNESS( &vec_atanf         , "vec_atanf.seq"           , 0 );
  DO_TEST_BITEXACTNESS( &vec_atan2f        , "vec_atan2f.seq"          , 0 );
  DO_TEST_BITEXACTNESS( &vec_tanhf         , "vec_tanhf.seq"           , 0 );
  DO_TEST_BITEXACTNESS( &vec_sigmoidf      , "vec_sigmoidf.seq"        , 0 );
  DO_TEST             ( &vec_softmaxf      , "vec_softmaxf.seq"        , 0 );
  DO_TEST             ( &vec_reluf         , "vec_reluf.seq"           , 0 );
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
        FUNC_LIST((tTestEngTarget)&vec_powf),
        TEST_DESC(FMT_REAL | FMT_FLOAT32,OVLP_XZ|OVLP_YZ, TE_DIM_NUM_1, TE_ALIGN_NO, &te_loadFxn_vXvYvZ, &te_processFxn_vZvXvY) },
    { 
      FUNC_LIST( (tTestEngTarget)&vec_int2float ),
      TEST_DESC( FMT_REAL|FMT_FLOAT32, 0,TE_DIM_NUM_2, TE_ALIGN_NO, &te_loadFxn_vX32sY32vZ, &te_processFxn_vZvXsY32 ) },
    { 
      FUNC_LIST( (tTestEngTarget)&vec_float2int ),
      TEST_DESC( FMT_REAL|FMT_FLOAT32, 0,TE_DIM_NUM_2, TE_ALIGN_NO, &te_loadFxn_vXsY32vZ32, &te_processFxn_vZvXsY32 ) },
    {
      FUNC_LIST( (tTestEngTarget)&vec_sinef,(tTestEngTarget)&vec_cosinef,(tTestEngTarget)&vec_tanf,
                 (tTestEngTarget)&vec_log2f,(tTestEngTarget)&vec_lognf,(tTestEngTarget)&vec_log10f,
                 (tTestEngTarget)&vec_atanf,(tTestEngTarget)&vec_antilog2f,(tTestEngTarget)&vec_antilognf,
                 (tTestEngTarget)&vec_antilog10f),
      TEST_DESC( FMT_REAL|FMT_FLOAT32, OVLP_XZ, TE_DIM_NUM_1, TE_ALIGN_NO, &te_loadFxn_vXvZ, &te_processFxn_vZvX ) },
    {
      FUNC_LIST( (tTestEngTarget)&vec_atan2f ),
      TEST_DESC( FMT_FLOAT32, OVLP_XZ|OVLP_YZ,TE_DIM_NUM_1, TE_ALIGN_NO, &te_loadFxn_vXvYvZ, &te_processFxn_vZvYvX ) },
    {
      FUNC_LIST( (tTestEngTarget)&vec_tanhf, (tTestEngTarget)&vec_sigmoidf, (tTestEngTarget)&vec_softmaxf),
      TEST_DESC( FMT_FLOAT32, OVLP_XZ,TE_DIM_NUM_1, TE_ALIGN_NO, &te_loadFxn_vXvZ, &te_processFxn_vZvX ) },
    {
      FUNC_LIST( (tTestEngTarget)&vec_reluf ),
      TEST_DESC( FMT_FLOAT32, OVLP_XZ,TE_DIM_NUM_1, TE_ALIGN_NO, &te_loadFxn_vXsYvZ, &te_processFxn_vZ_vXsY ) },
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
