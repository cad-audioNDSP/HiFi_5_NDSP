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
 * Test procedures for FIR
 */

/* Cross-platform data type definitions. */
#include "types.h"
/* Test environment configuration. */
#include "config.h"
#include "packages.h"
/* DSP Library API */
#include LIBRARY_HEADER(fir)
/* Test data vectors tools and SEQ-file reader. */
#include "vectools.h"
#include "testeng.h"
#include "malloc16.h"
#include <string.h>

#define MAX(a,b)   ( (a)>(b) ? (a) : (b) )

#define MAX_FUNC_NUM   2
/* Initializer for a function pointer array, appends NULL to a sequence of pointers. */
#define FUNC_LIST(...) { __VA_ARGS__, NULL }
/* Initializer for a test description structure. */
#define TEST_DESC( fmt, argNum, extraPtr,align, loadFxn, procFxn ) { (fmt),0,extraPtr,(argNum),(align),NULL,NULL,(loadFxn),(procFxn) }

typedef struct
{
    size_t (*getScratchSize)(int P,int Q);
    uint32_t fmtx;
    uint32_t fmty;
}
tConv2dApi;

static const tConv2dApi apiconv2d_gen_3x3_fp16 = {conv2d_gen_3x3_fp16_getScratchSize  ,FMT_REAL|FMT_FLOAT16,FMT_REAL|FMT_FLOAT16 };
static const tConv2dApi apiconv2d_gen_5x5_fp16 = {conv2d_gen_5x5_fp16_getScratchSize  ,FMT_REAL|FMT_FLOAT16,FMT_REAL|FMT_FLOAT16 };
static const tConv2dApi apiconv2d_gen_11x7_fp16= {conv2d_gen_11x7_fp16_getScratchSize ,FMT_REAL|FMT_FLOAT16,FMT_REAL|FMT_FLOAT16 };

static const tConv2dApi apiconv2d_3x3_fp16 = {conv2d_3x3_fp16_getScratchSize  ,FMT_REAL|FMT_FLOAT16,FMT_REAL|FMT_FLOAT16 };
static const tConv2dApi apiconv2d_5x5_fp16 = {conv2d_5x5_fp16_getScratchSize  ,FMT_REAL|FMT_FLOAT16,FMT_REAL|FMT_FLOAT16 };
static const tConv2dApi apiconv2d_11x7_fp16= {conv2d_11x7_fp16_getScratchSize ,FMT_REAL|FMT_FLOAT16,FMT_REAL|FMT_FLOAT16 };

static int te_loadFxn_conv2d(tTestEngContext * context)
{
  const tConv2dApi* api=(const tConv2dApi*)context->desc->extraPtr;
  int M, N, P, Q;
  int nElemX,nElemY,nElemZ, res;

  ASSERT( context && context->seqFile );

  M = MAX( 0, context->args.dim[0] );
  N = MAX( 0, context->args.dim[1] );
  P = MAX( 0, context->args.dim[2] );
  Q = MAX( 0, context->args.dim[3] );

  nElemX = M*N;
  nElemY = P*Q;
  nElemZ = (nElemX==0 || nElemY==0) ? 0: (M+P-1)*(N+Q-1);

  memset( &context->dataSet, 0, sizeof(context->dataSet) );

  /* Allocate data vectors memory. */
  res  =  vecsAlloc( context->desc->isAligned, api->fmtx,
                         &context->dataSet.X  , nElemX, 0 ) ;
  res +=  vecsAlloc( context->desc->isAligned, api->fmty,
                         &context->dataSet.Y  , nElemY,
                         &context->dataSet.Z  , nElemZ,
                         &context->dataSet.Zlo, nElemZ,
                         &context->dataSet.Zhi, nElemZ, 0 ) ;
  res=res==5;
  if (res)
  {
    /* Load vectors data from the SEQ-file. */
    if ( !( res = seqFileReadVecs( context->seqFile,
                                  &context->dataSet.X,
                                  &context->dataSet.Y,
                                  &context->dataSet.Zlo,
                                  &context->dataSet.Zhi, 0 ) ) )
    {
      printf( "te_loadFxn_conv2d(): failed to read vectors data\n");
    }
  }
  else
  {
    printf( "te_loadFxn_conv2d(): failed to allocate vectors\n");
  }

  /* Free vectors data if failed. */
  if ( !res ) te_freeVectors(context);

  return (res);
}

static void te_processFxn_conv2d(tTestEngContext * context)
{
  typedef void tFxn(void* pScr,void * Z, const void * X, const void * Y, int rsh ,int P, int Q );
  typedef void tFxnFloat(void* pScr,void * Z, const void * X, const void * Y,  int P, int Q );
  tFxn *fxn;
  void *X, *Y, *Z,*pScr;
  int P,Q,rsh;
  const tConv2dApi* api=(const tConv2dApi*)context->desc->extraPtr;

  ASSERT( context && context->target.fut );

  X = vecGetElem( &context->dataSet.X, 0 );
  Y = vecGetElem( &context->dataSet.Y, 0 );
  Z = vecGetElem( &context->dataSet.Z, 0 );

  P = context->args.dim[2] ;
  Q = context->args.dim[3] ;
  rsh = context->args.dim[4] ;
  fxn = (tFxn*)context->target.fut;

  te_vReportStd(context);

  pScr=mallocAlign(api->getScratchSize(P,Q),0);
  if(api->fmtx==(FMT_REAL|FMT_FLOAT32) || api->fmtx == (FMT_REAL | FMT_FLOAT16))
  {
      ((tFxnFloat*)fxn)(pScr, Z, X, Y, P, Q);
  }
  else
  {
      fxn(pScr, Z, X, Y, rsh, P, Q);
  }
  freeAlign(pScr);
}

/* vec API test definitions. */
static const struct 
{
tTestEngTarget   funcList[MAX_FUNC_NUM];
tTestEngDesc     testDesc;
}
testDefTbl[] =
{
    {   FUNC_LIST((tTestEngTarget)&conv2d_3x3_fp16  ),TEST_DESC( FMT_REAL|FMT_FLOAT16 , TE_DIM_NUM_5, &apiconv2d_3x3_fp16 , TE_ALIGN_YES, &te_loadFxn_conv2d, &te_processFxn_conv2d ) },
    {   FUNC_LIST((tTestEngTarget)&conv2d_5x5_fp16  ),TEST_DESC( FMT_REAL|FMT_FLOAT16 , TE_DIM_NUM_5, &apiconv2d_5x5_fp16 , TE_ALIGN_YES, &te_loadFxn_conv2d, &te_processFxn_conv2d ) },
    {   FUNC_LIST((tTestEngTarget)&conv2d_11x7_fp16 ),TEST_DESC( FMT_REAL|FMT_FLOAT16 , TE_DIM_NUM_5, &apiconv2d_11x7_fp16, TE_ALIGN_YES, &te_loadFxn_conv2d, &te_processFxn_conv2d ) },
    {   FUNC_LIST((tTestEngTarget)&conv2d_gen_3x3_fp16)  ,TEST_DESC(FMT_REAL | FMT_FLOAT16 , TE_DIM_NUM_5, &apiconv2d_gen_3x3_fp16 , TE_ALIGN_NO, &te_loadFxn_conv2d, &te_processFxn_conv2d) },
    {   FUNC_LIST((tTestEngTarget)&conv2d_gen_5x5_fp16)  ,TEST_DESC(FMT_REAL | FMT_FLOAT16 , TE_DIM_NUM_5, &apiconv2d_gen_5x5_fp16 , TE_ALIGN_NO, &te_loadFxn_conv2d, &te_processFxn_conv2d) },
    {   FUNC_LIST((tTestEngTarget)&conv2d_gen_11x7_fp16) ,TEST_DESC(FMT_REAL | FMT_FLOAT16 , TE_DIM_NUM_5, &apiconv2d_gen_11x7_fp16, TE_ALIGN_NO, &te_loadFxn_conv2d, &te_processFxn_conv2d) },
    { {NULL}  } /* End of table */
};

/* Test executive function. Performs the specified test on a brief or full version
 * of the designated SEQ-file. Return the test result (non-zero if passed). */
static int testExec( tTestEngTarget  targetFxn, const char * seqName,
                     int isFull, int isVerbose, int breakOnError )
{
    return te_Exec(testDefTbl,sizeof(testDefTbl)/sizeof(testDefTbl[0]),MAX_FUNC_NUM,targetFxn, seqName,isFull, isVerbose, breakOnError);
}

/* Perform all tests for conv2d API functions. */
int func_conv2d3(int isFull, int isVerbose, int breakOnError)
{
    int res = 1;

    #define DO_TEST(fxn, seqFile) if ( res || !breakOnError ) res &= ( 0 != testExec((tTestEngTarget)(fxn), "conv2d3/" seqFile, isFull, isVerbose, breakOnError ) )

    DO_TEST(&conv2d_3x3_fp16    , "conv2d_3x3_fp16.seq"     );
    DO_TEST(&conv2d_gen_3x3_fp16, "conv2d_gen_3x3_fp16.seq" );
    DO_TEST(&conv2d_5x5_fp16    , "conv2d_5x5_fp16.seq"     );
    DO_TEST(&conv2d_gen_5x5_fp16, "conv2d_gen_5x5_fp16.seq" );
    DO_TEST(&conv2d_11x7_fp16   , "conv2d_11x7_fp16.seq"    );
    DO_TEST(&conv2d_gen_11x7_fp16,"conv2d_gen_11x7_fp16.seq");
  
    return (res);
}
