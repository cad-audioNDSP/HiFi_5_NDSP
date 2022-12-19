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
/* Library API. */
#include LIBRARY_HEADER(vector)
/* Test Engine add-on for Vector Mathematics functions */
#include "../common/testeng_vector.h"

#define MIN(a,b)   ( (a)<(b) ? (a) : (b) )
#define MAX(a,b)   ( (a)>(b) ? (a) : (b) )

#define MAX_FUNC_NUM   10
/* Initializer for a function pointer array, appends NULL to a sequence of pointers. */
#define FUNC_LIST(...) { __VA_ARGS__, NULL }
/* Initializer for a test description structure. */
#define TEST_DESC( fmt, dimNum, align, loadFxn, procFxn ) { (fmt),0,NULL,(dimNum),(align),NULL,NULL,(loadFxn),(procFxn) }
#define TEST_DESC_EX( fmt, dimNum, align, loadFxn, procFxn,api ) { (fmt),0,(void*)api,(dimNum),(align),NULL,NULL,(loadFxn),(procFxn) }

typedef struct
{
	int fmtX;
	int fmtY;
	int fmtZ;
	int sizeofX;
	int sizeofY;
	int sizeofZ;
}
vec_dot_batch_api;

static const vec_dot_batch_api vec_dot_batch_fp16_api={FMT_REAL|FMT_FLOAT16,FMT_REAL|FMT_FLOAT16,FMT_REAL|FMT_FLOAT16,sizeof(float16_t),sizeof(float16_t),sizeof(float16_t)};

/* Allocate vectors and load the data set for vec_dot_batch*/
static int te_loadFxn_vec_dot_batch( tTestEngContext * context )
{
    const vec_dot_batch_api* pApi=(const vec_dot_batch_api*)context->desc->extraPtr;
    int res=0;
    int M,N,P,rsh;
    int nAlloc;
    int fmtPtr;
    ASSERT( context && context->seqFile );

    M   = MAX( 0, context->args.dim[0] );
    N   = MAX( 0, context->args.dim[1] );
    P   = MAX( 0, context->args.dim[2] );
    rsh = context->args.dim[3];
    (void)M,(void)N,(void)P,(void)rsh;
    fmtPtr = sizeof(uintptr_t)==8 ? FMT_REAL|FMT_INT64: FMT_REAL|FMT_INT32;

    /* Allocate output/reference data vectors memory. */
    nAlloc =vecsAlloc( context->desc->isAligned, pApi->fmtX,&context->dataSet.X  , N,NULL);
    nAlloc+=vecsAlloc( context->desc->isAligned, pApi->fmtY,&context->dataSet.Y  , P,NULL);
    nAlloc+=vecsAlloc( context->desc->isAligned, fmtPtr    ,&context->dataSet.U  , M,NULL);
    nAlloc+=vecsAlloc( context->desc->isAligned, pApi->fmtZ,
                        &context->dataSet.Z  , M,
                        &context->dataSet.Zlo, M,
                        &context->dataSet.Zhi, M,
                        0 );
    if ( 6 !=  nAlloc)
    {
        printf( "te_loadFxn_vec_dot_batch(): failed to allocate data\n");
    }
    /* Load vectors data from the SEQ-file. */
    else if ( !seqFileReadVecs( context->seqFile,
                                &context->dataSet.U,
                                &context->dataSet.X,
                                &context->dataSet.Y,
                                &context->dataSet.Zlo,
                                &context->dataSet.Zhi, 
                                0 ) )
    {
        printf( "te_loadFxn_vec_dot_batch(): failed to read data\n");
        te_freeVectors(context);
    }
    else
    {
        res = 1;
    }
    return (res);
} 

/* test interleave functions */
static void te_processFxn_vec_dot_batch( tTestEngContext * context )
{
    const vec_dot_batch_api* pApi=(const vec_dot_batch_api*)context->desc->extraPtr;
    typedef void tFxn(void   *z, const void* x,const void* y,int rsh, int N, int M);
    typedef void tFxn_float32(void   *z, const void* x,const void* y,int N, int M);
    tTestEngTarget   fxn;
    int M,N,P,rsh,m;
    uintptr_t* y;
    uintptr_t y0;

    M   = MAX( 0, context->args.dim[0] );
    N   = MAX( 0, context->args.dim[1] );
    P   = MAX( 0, context->args.dim[2] );
    rsh = context->args.dim[3];
    (void)M,(void)N,(void)P,(void)rsh;

    ASSERT( context && context->target.fut );
    te_vReportStd(context);

    fxn = context->target.fut;
    y=(uintptr_t*)vecGetElem( &context->dataSet.U, 0 );
    y0=(uintptr_t)(uintptr_t)vecGetElem( &context->dataSet.Y, 0 );
    for (m=0; m<M; m++) 
    {
        y[m]=(y[m]*pApi->sizeofY)+y0;
    }

    if(pApi->fmtX==(FMT_REAL|FMT_FLOAT32) || pApi->fmtX == (FMT_REAL | FMT_FLOAT16))
    {
        ( (tFxn_float32*)fxn )( vecGetElem( &context->dataSet.Z, 0 ), 
                         vecGetElem( &context->dataSet.X, 0 ), (void*)y,N,M);
    }
    else
    {
        ( (tFxn *)fxn )( vecGetElem( &context->dataSet.Z, 0 ), 
                         vecGetElem( &context->dataSet.X, 0 ), 
                         (void*)y,rsh,N,M);
    }
} /* te_processFxn_vec_dot_batch() */

/* vec API test definitions. */
static const struct 
{
  tTestEngTarget   funcList[MAX_FUNC_NUM];
  tTestEngDesc     testDesc;
}
testDefTbl[] =
{
  {
    FUNC_LIST((tTestEngTarget)&vec_dot_batch_fp16),
    TEST_DESC_EX(FMT_REAL | FMT_FLOAT16, TE_DIM_NUM_4, TE_ALIGN_NO, &te_loadFxn_vec_dot_batch, &te_processFxn_vec_dot_batch,&vec_dot_batch_fp16_api) },
  {
    FUNC_LIST((tTestEngTarget)&vec_dot_batch_fp16_fast),
	TEST_DESC_EX(FMT_REAL | FMT_FLOAT16, TE_DIM_NUM_4, TE_ALIGN_YES, &te_loadFxn_vec_dot_batch, &te_processFxn_vec_dot_batch, &vec_dot_batch_fp16_api) },
  {
	FUNC_LIST(NULL), TEST_DESC(0, 0, 0, NULL, NULL) } /* End of table */
};

/* Perform all functional tests for Vector Mathematics API functions. */
int func_vector3(int isFull, int isVerbose, int breakOnError)
{
	int res = 1;

#define DO_TEST(fxn, seqFile)                                                                   \
    if ( res || !breakOnError ) res &= ( 0 != te_Exec( testDefTbl,                              \
                                                       sizeof(testDefTbl)/sizeof(testDefTbl[0]),\
                                                       MAX_FUNC_NUM,                            \
                                                       (tTestEngTarget)(fxn), "vec3/" seqFile,   \
                                                       isFull, isVerbose, breakOnError ) )
	DO_TEST(&vec_dot_batch_fp16		, "vec_dot_batch_fp16.seq");
	DO_TEST(&vec_dot_batch_fp16_fast, "vec_dot_batch_fp16_fast.seq");
	return (res);
}

