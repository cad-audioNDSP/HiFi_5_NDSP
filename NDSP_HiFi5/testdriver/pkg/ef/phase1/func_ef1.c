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
 * Test procedures for emulated float
 */
#include <string.h> 

/* Cross-platform data type definitions. */
#include "types.h"
/* Test environment configuration. */
#include "config.h"
#include "packages.h"
/* DSP Library API: arithmetic and logic functions on data vectors. */
#include LIBRARY_HEADER(vector)
/* Test engine API. */
#include "../../vector/common/testeng_vector.h"

#define MIN(a,b)   ( (a)<(b) ? (a) : (b) )
#define MAX(a,b)   ( (a)>(b) ? (a) : (b) )

#define MAX_FUNC_NUM   10
/* Initializer for a function pointer array, appends NULL to a sequence of pointers. */
#define FUNC_LIST(...) { __VA_ARGS__, NULL }
/* Initializer for a test description structure. */
#define TEST_DESC( fmt, dimNum, align, loadFxn, procFxn ) { (fmt),0,NULL,(dimNum),(align),NULL,NULL,(loadFxn),(procFxn) }

tTestEngProcessFxn te_processFxn_addmulef;     /* routine for emulated float add/mul */
tTestEngProcessFxn te_processFxn_macef;        /* routine for emulated float mac */
tTestEngProcessFxn te_processFxn_dotef;        /* routine for emulated float dot product */
tTestEngLoadFxn te_loadFxn_macef;

typedef struct 
{
    int32_t *mant;
    int16_t *exp;
} tEfPtr;

// get emulated float pointers from original uint64_t data
static tEfPtr getEfPtr(uint64_t* x,int N)
{
    tEfPtr ef;
    NASSERT(x);
    ef.mant=(int32_t*)x;
    ef.exp =(int16_t*)(ef.mant+N);
    return ef;
}

/* in-place splitting ef data */
tEfPtr inplace_uint64_ef(uint64_t *x,int N)
{
    tEfPtr ef={0,0};
    tVec mant,exp;
    if (N<=0) return ef;
    if (vecAlloc( &mant, N, 1, FMT_REAL|FMT_FRACT32, NULL) &&
        vecAlloc( &exp , N, 1, FMT_REAL|FMT_FRACT16, NULL))
    {
        int n;
        for(n=0; n<N; n++)
        {
            vecGetElem_fr32(&mant,n)[0]=(uint32_t)x[n];
            vecGetElem_fr16(&exp ,n)[0]=(uint16_t)(x[n]>>32);
        }
        ef=getEfPtr(x,N);
        memcpy(ef.mant,vecGetElem_fr32(&mant,0),N*sizeof(int32_t));
        memcpy(ef.exp ,vecGetElem_fr16(&exp ,0),N*sizeof(int16_t));
        vecsFree(&mant,&exp,NULL);
    }
    // failed 
    ASSERT("inplace_uint64_ef failed");
    return ef;
}

// combine uint64_t from emulated float
static void inplace_ef_uint64(uint64_t *x,int N)
{
    tEfPtr ef=getEfPtr(x,N);
    tVec temp;
    if (N<=0) return;
    if (vecAlloc( &temp, N, 1, FMT_REAL|FMT_INT64, NULL))
    {
        int n;
        for(n=0; n<N; n++)
        {
            vecGetElem_i64(&temp,n)[0]=(((uint64_t)ef.mant[n])&0xffffffff)|(((uint64_t)ef.exp[n])<<32);
        }
        memcpy(x,vecGetElem_i64(&temp,0),N*sizeof(int64_t));
        vecsFree(&temp,NULL);
    }
    // failed 
    ASSERT("inplace_ef_uint64 failed"!=NULL);
}

// processing function for vector add/mul operation
void te_processFxn_addmulef(tTestEngContext * context)
{
    tEfPtr efX,efY,efZ;
    typedef void fxn (      int32_t  *  zmant,       int16_t  *  zexp, 
                      const int32_t  *  xmant, const int16_t  *  xexp, 
                      const int32_t  *  ymant, const int16_t  *  yexp, 
                      int N);
    int N;
    uint64_t *X;
    uint64_t *Y;
    uint64_t *Z;
    ASSERT(context && context->target.fut);
    te_vReportStd(context);
    X = (uint64_t*)vecGetElem(&context->dataSet.X, 0);
    Y = (uint64_t*)vecGetElem(&context->dataSet.Y, 0);
    Z = (uint64_t*)vecGetElem(&context->dataSet.Z, 0);
    N = context->args.dim[0];
    efX=inplace_uint64_ef(X,N);
    efY=inplace_uint64_ef(Y,N);
    efZ=getEfPtr(Z,N);
    ((fxn*)context->target.fut)(efZ.mant,efZ.exp,efX.mant,efX.exp,efY.mant,efY.exp,N);
    inplace_ef_uint64(Z,N);
}

// processing function for scalar add/mul operation
void te_processFxn_scladdmulef(tTestEngContext * context)
{
    tEfPtr efX,efY,efZ;
    typedef void fxn (int32_t  *  zmant, int16_t  *  zexp, 
                      int32_t     xmant, int16_t     xexp, 
                      int32_t     ymant, int16_t     yexp);
    int n,N;
    uint64_t *X;
    uint64_t *Y;
    uint64_t *Z;
    ASSERT(context && context->target.fut);
    te_vReportStd(context);
    X = (uint64_t*)vecGetElem(&context->dataSet.X, 0);
    Y = (uint64_t*)vecGetElem(&context->dataSet.Y, 0);
    Z = (uint64_t*)vecGetElem(&context->dataSet.Z, 0);
    N = context->args.dim[0];
    efX=inplace_uint64_ef(X,N);
    efY=inplace_uint64_ef(Y,N);
    efZ=getEfPtr(Z,N);
    for (n=0; n<N; n++)
    {
        ((fxn*)context->target.fut)(efZ.mant+n,efZ.exp+n,efX.mant[n],efX.exp[n],efY.mant[n],efY.exp[n]);
    }
    inplace_ef_uint64(Z,N);
}

// processing function for vector mac operation
void te_processFxn_macef(tTestEngContext * context)
{
    tEfPtr efX,efY,efZ;
    typedef void fxn (      int32_t  *  zmant,       int16_t  *  zexp, 
                      const int32_t  *  xmant, const int16_t  *  xexp, 
                      const int32_t     ymant, const int16_t     yexp, 
                      int N);
    int N;
    uint64_t *X;
    uint64_t *Y;
    uint64_t *Z;
    ASSERT(context && context->target.fut);
    te_vReportStd(context);
    X = (uint64_t*)vecGetElem(&context->dataSet.X, 0);
    Y = (uint64_t*)vecGetElem(&context->dataSet.Y, 0);
    Z = (uint64_t*)vecGetElem(&context->dataSet.Z, 0);
    N = context->args.dim[0];
    efX=inplace_uint64_ef(X,N);
    efY=inplace_uint64_ef(Y,1);
    efZ=inplace_uint64_ef(Z,N);
    ((fxn*)context->target.fut)(efZ.mant,efZ.exp,efX.mant,efX.exp,efY.mant[0],efY.exp[0],N);
    inplace_ef_uint64(Z,N);
}

// processing function for scalar mac operation
void te_processFxn_sclmacef(tTestEngContext * context)
{
    tEfPtr efX,efY,efZ;
    typedef void fxn ( int32_t  *  zmant, int16_t  *  zexp, 
                       int32_t     xmant, int16_t     xexp, 
                       int32_t     ymant, int16_t     yexp);
    int n,N;
    uint64_t *X;
    uint64_t *Y;
    uint64_t *Z;
    ASSERT(context && context->target.fut);
    te_vReportStd(context);
    X = (uint64_t*)vecGetElem(&context->dataSet.X, 0);
    Y = (uint64_t*)vecGetElem(&context->dataSet.Y, 0);
    Z = (uint64_t*)vecGetElem(&context->dataSet.Z, 0);
    N = context->args.dim[0];
    efX=inplace_uint64_ef(X,N);
    efY=inplace_uint64_ef(Y,1);
    efZ=inplace_uint64_ef(Z,N);
    for (n=0; n<N; n++)
    {
        ((fxn*)context->target.fut)(efZ.mant+n,efZ.exp+n,efX.mant[n],efX.exp[n],efY.mant[0],efY.exp[0]);
    }
    inplace_ef_uint64(Z,N);
}

// processing function for vector dot product
void te_processFxn_dotef(tTestEngContext * context)
{
    tEfPtr efX,efY,efZ;
    typedef void fxn (      int32_t  *  zmant,       int16_t  *  zexp, 
                      const int32_t  *  xmant, const int16_t  *  xexp, 
                      const int32_t  *  ymant, const int16_t  *  yexp, 
                      int N);
    int N;
    uint64_t *X;
    uint64_t *Y;
    uint64_t *Z;
    ASSERT(context && context->target.fut);
    te_vReportStd(context);
    X = (uint64_t*)vecGetElem(&context->dataSet.X, 0);
    Y = (uint64_t*)vecGetElem(&context->dataSet.Y, 0);
    Z = (uint64_t*)vecGetElem(&context->dataSet.Z, 0);
    N = context->args.dim[0];
    efX=inplace_uint64_ef(X,N);
    efY=inplace_uint64_ef(Y,N);
    efZ=getEfPtr(Z,1);
    ((fxn*)context->target.fut)(efZ.mant,efZ.exp,efX.mant,efX.exp,efY.mant,efY.exp,N);
    inplace_ef_uint64(Z,1);
}

// load data for mac operation
int te_loadFxn_macef(tTestEngContext * context)
{
  int M, N, L;
  int fmt, nElemXZ;
  int res = 0;

  ASSERT( context && context->seqFile );
  N = MAX( 0, context->args.dim[0] );
  M = MAX( 0, context->args.dim[1] );
  L = MAX( 0, context->args.dim[2] );
  nElemXZ = M*N*L;
  fmt = context->desc->fmt;
  ASSERT( context && context->seqFile );
  memset( &context->dataSet, 0, sizeof(context->dataSet) );

  if (5 != vecsAlloc( context->desc->isAligned, fmt,
                           &context->dataSet.X  , nElemXZ,
                           &context->dataSet.Y  , 1,
                           &context->dataSet.Z  , nElemXZ,
                           &context->dataSet.Zlo, nElemXZ,
                           &context->dataSet.Zhi, nElemXZ, 0 ))
  {
    printf( "te_loadFxn_macef(): failed to allocate data\n ");
  }
  /* Load vectors data from the SEQ-file. */
  else if ( !seqFileReadVecs( context->seqFile,
                             &context->dataSet.X,
                             &context->dataSet.Y,
                             &context->dataSet.Z,
                             &context->dataSet.Zlo,
                             &context->dataSet.Zhi, 0 ) )
  {
    printf( "te_loadFxn_macef(): failed to read vectors data;\n");
  }
  else
  {
    res = 1;
  }

  /* Free vectors data if failed. */
  if ( !res ) te_freeVectors(context);
  return (res);
}

/* vec API test definitions. */
static const struct 
{
  tTestEngTarget   funcList[MAX_FUNC_NUM];
  tTestEngDesc     testDesc;
}
testDefTbl[] =
{
  /*
   * Stage 1
   */
  {
    FUNC_LIST( (tTestEngTarget)&vec_add_32x16ef,(tTestEngTarget)&vec_mul_32x16ef ),
    TEST_DESC( FMT_REAL|FMT_EF32X16, TE_DIM_NUM_2, TE_ALIGN_NO, &te_loadFxn_vXvYvZ, &te_processFxn_addmulef ) },
  {
    FUNC_LIST( (tTestEngTarget)&vec_mac_32x16ef),
    TEST_DESC( FMT_REAL|FMT_EF32X16, TE_DIM_NUM_1, TE_ALIGN_NO, &te_loadFxn_macef, &te_processFxn_macef ) },
  {
    FUNC_LIST( (tTestEngTarget)&vec_dot_32x16ef),
    TEST_DESC(FMT_REAL | FMT_EF32X16, TE_DIM_NUM_1, TE_ALIGN_NO, &te_loadFxn_vXvYsZ, &te_processFxn_dotef) },
  {
    FUNC_LIST( (tTestEngTarget)&scl_add_32x16ef,(tTestEngTarget)&scl_mul_32x16ef ),
    TEST_DESC( FMT_REAL|FMT_EF32X16, TE_DIM_NUM_2, TE_ALIGN_NO, &te_loadFxn_vXvYvZ, &te_processFxn_scladdmulef ) },
  {
    FUNC_LIST( (tTestEngTarget)&scl_mac_32x16ef),
    TEST_DESC( FMT_REAL|FMT_EF32X16, TE_DIM_NUM_1, TE_ALIGN_NO, &te_loadFxn_macef, &te_processFxn_sclmacef ) },

  { 
    FUNC_LIST( NULL ), TEST_DESC(  0, 0, 0, NULL, NULL ) } /* End of table */
};

/* Perform all tests for emulated float API functions. */
int func_ef1(int isFull, int isVerbose, int breakOnError)
{
	int res = 1;

#define DO_TEST(fxn, seqFile)                                                                \
    if ( res || !breakOnError ) res &= ( 0 != te_Exec( testDefTbl,                               \
                                                       sizeof(testDefTbl)/sizeof(testDefTbl[0]), \
                                                       MAX_FUNC_NUM,                             \
                                                       (tTestEngTarget)(fxn), "ef1/" seqFile,     \
                                                       isFull, isVerbose, breakOnError ) )

	DO_TEST(&vec_add_32x16ef, "vec_add_32x16ef.seq");
	DO_TEST(&vec_mul_32x16ef, "vec_mul_32x16ef.seq");
	DO_TEST(&vec_mac_32x16ef, "vec_mac_32x16ef.seq");
	DO_TEST(&vec_dot_32x16ef, "vec_dot_32x16ef.seq");
	DO_TEST(&scl_add_32x16ef, "scl_add_32x16ef.seq");
	DO_TEST(&scl_mul_32x16ef, "scl_mul_32x16ef.seq");
	DO_TEST(&scl_mac_32x16ef, "scl_mac_32x16ef.seq");
	
	return (res);

}
