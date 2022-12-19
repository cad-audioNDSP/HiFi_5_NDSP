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

static const vec_dot_batch_api vec_dot_batch16x16_api={FMT_REAL|FMT_INT16  ,FMT_REAL|FMT_INT16  ,FMT_REAL|FMT_INT32  ,sizeof(int16_t)  ,sizeof(int16_t)  ,sizeof(int32_t)}  ;
static const vec_dot_batch_api vec_dot_batch8x16_api ={FMT_REAL|FMT_INT8   ,FMT_REAL|FMT_INT16  ,FMT_REAL|FMT_INT16  ,sizeof(int8_t )  ,sizeof(int16_t)  ,sizeof(int16_t)}  ;
static const vec_dot_batch_api vec_dot_batch8x8_api  ={FMT_REAL|FMT_INT8   ,FMT_REAL|FMT_INT8   ,FMT_REAL|FMT_INT16  ,sizeof(int8_t )  ,sizeof(int8_t )  ,sizeof(int16_t)}  ;
static const vec_dot_batch_api vec_dot_batchf_api    ={FMT_REAL|FMT_FLOAT32,FMT_REAL|FMT_FLOAT32,FMT_REAL|FMT_FLOAT32,sizeof(float32_t),sizeof(float32_t),sizeof(float32_t)};

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
  /*
  * Stage 1
  */
#if 0// HiFi3/3z API
  {
    FUNC_LIST( (tTestEngTarget)&vec_dot24x24 ),
    TEST_DESC( FMT_REAL|FMT_FRACT32, TE_DIM_NUM_2, TE_ALIGN_NO, &te_loadFxn_vXvYsZ64, &te_processFxn_sZ64vXvY ) },
  {
    FUNC_LIST( (tTestEngTarget)&vec_dot24x24_fast ),
    TEST_DESC( FMT_REAL|FMT_FRACT32, TE_DIM_NUM_2, TE_ALIGN_YES, &te_loadFxn_vXvYsZ64, &te_processFxn_sZ64vXvY ) },
  {
    FUNC_LIST( (tTestEngTarget)&vec_add24x24 ),
    TEST_DESC( FMT_REAL|FMT_FRACT32, TE_DIM_NUM_2, TE_ALIGN_NO, &te_loadFxn_vXvYvZ, &te_processFxn_vZvXvY ) },
  {
    FUNC_LIST( (tTestEngTarget)&vec_add24x24_fast ),
    TEST_DESC( FMT_REAL|FMT_FRACT32, TE_DIM_NUM_2, TE_ALIGN_YES, &te_loadFxn_vXvYvZ, &te_processFxn_vZvXvY ) },
  {
    FUNC_LIST( (tTestEngTarget)&vec_power24x24 ),
    TEST_DESC( FMT_REAL|FMT_FRACT32, TE_DIM_NUM_2, TE_ALIGN_NO, &te_loadFxn_vXsY32sZ64, &te_processFxn_vXsY32sZ64 ) },
  {
    FUNC_LIST( (tTestEngTarget)&vec_power24x24_fast ),
    TEST_DESC( FMT_REAL|FMT_FRACT32, TE_DIM_NUM_2, TE_ALIGN_YES, &te_loadFxn_vXsY32sZ64, &te_processFxn_vXsY32sZ64 ) },
  { 
    FUNC_LIST((tTestEngTarget)&vec_min24x24,
              (tTestEngTarget)&vec_max24x24),
    TEST_DESC( FMT_REAL|FMT_FRACT32, TE_DIM_NUM_2, TE_ALIGN_NO, &te_loadFxn_vXsZ  , &te_processFxn_vXsZ ) },
  { 
    FUNC_LIST( (tTestEngTarget)&vec_min24x24_fast,
               (tTestEngTarget)&vec_max24x24_fast),
    TEST_DESC( FMT_REAL|FMT_FRACT32, TE_DIM_NUM_2, TE_ALIGN_YES, &te_loadFxn_vXsZ  , &te_processFxn_vXsZ ) },
  { 
    FUNC_LIST((tTestEngTarget)&(tTestEngTarget)&vec_scale24x24,
               (tTestEngTarget)&vec_shift24x24),
    TEST_DESC( FMT_REAL|FMT_FRACT32, TE_DIM_NUM_2, TE_ALIGN_NO, &te_loadFxn_vXsYvZ, &te_processFxn_vZvXsY ) },
  { 
    FUNC_LIST( (tTestEngTarget)&vec_scale32x24_fast, (tTestEngTarget)&vec_scale24x24_fast,
                ((tTestEngTarget)&vec_shift24x24_fast),
    TEST_DESC( FMT_REAL|FMT_FRACT32, TE_DIM_NUM_2, TE_ALIGN_YES, &te_loadFxn_vXsYvZ, &te_processFxn_vZvXsY ) },
  { 
    FUNC_LIST( (tTestEngTarget)&vec_bexp24 ),
    TEST_DESC( FMT_REAL|FMT_FRACT32, TE_DIM_NUM_2, TE_ALIGN_NO, &te_loadFxn_vXsZ  , &te_processFxn_vXsZ ) },
  { 
    FUNC_LIST( (tTestEngTarget)&vec_bexp24_fast ),
    TEST_DESC( FMT_REAL|FMT_FRACT32, TE_DIM_NUM_2, TE_ALIGN_YES, &te_loadFxn_vXsZ  , &te_processFxn_vXsZ ) },
  {
    FUNC_LIST( (tTestEngTarget)&scl_bexp24 ),
    TEST_DESC( FMT_FRACT32, TE_DIM_NUM_2, TE_ALIGN_NO, &te_loadFxn_vXsZ, &processFxn_scl_vXvZ32 ) },
  { 
    FUNC_LIST( (tTestEngTarget)&vec_scale32x24_fast),
    TEST_DESC( FMT_REAL|FMT_FRACT32, TE_DIM_NUM_2, TE_ALIGN_YES, &te_loadFxn_vXsYvZ, &te_processFxn_vZvXsY ) },
#endif
  {
    FUNC_LIST( (tTestEngTarget)&vec_dot_batch8x8 ),
    TEST_DESC_EX( FMT_REAL|FMT_FRACT16, TE_DIM_NUM_4, TE_ALIGN_NO, &te_loadFxn_vec_dot_batch, &te_processFxn_vec_dot_batch,&vec_dot_batch8x8_api ) },
  {
    FUNC_LIST( (tTestEngTarget)&vec_dot_batch8x8_fast ),
    TEST_DESC_EX( FMT_REAL|FMT_FRACT16, TE_DIM_NUM_4, TE_ALIGN_YES, &te_loadFxn_vec_dot_batch, &te_processFxn_vec_dot_batch,&vec_dot_batch8x8_api ) },
  {
    FUNC_LIST( (tTestEngTarget)&vec_dot_batch8x16 ),
    TEST_DESC_EX( FMT_REAL|FMT_FRACT16, TE_DIM_NUM_4, TE_ALIGN_NO, &te_loadFxn_vec_dot_batch, &te_processFxn_vec_dot_batch,&vec_dot_batch8x16_api ) },
  {
    FUNC_LIST( (tTestEngTarget)&vec_dot_batch8x16_fast ),
    TEST_DESC_EX( FMT_REAL|FMT_FRACT16, TE_DIM_NUM_4, TE_ALIGN_YES, &te_loadFxn_vec_dot_batch, &te_processFxn_vec_dot_batch,&vec_dot_batch8x16_api ) },
  {
    FUNC_LIST( (tTestEngTarget)&vec_dot_batch16x16 ),
    TEST_DESC_EX( FMT_REAL|FMT_FRACT16, TE_DIM_NUM_4, TE_ALIGN_NO, &te_loadFxn_vec_dot_batch, &te_processFxn_vec_dot_batch,&vec_dot_batch16x16_api ) },
  {
    FUNC_LIST( (tTestEngTarget)&vec_dot_batch16x16_fast ),
    TEST_DESC_EX( FMT_REAL|FMT_FRACT16, TE_DIM_NUM_4, TE_ALIGN_YES, &te_loadFxn_vec_dot_batch, &te_processFxn_vec_dot_batch,&vec_dot_batch16x16_api ) },
  {
    FUNC_LIST( (tTestEngTarget)&vec_dot_batchf ),
    TEST_DESC_EX( FMT_REAL|FMT_FLOAT32, TE_DIM_NUM_4, TE_ALIGN_NO, &te_loadFxn_vec_dot_batch, &te_processFxn_vec_dot_batch,&vec_dot_batchf_api ) },
  {
    FUNC_LIST( (tTestEngTarget)&vec_dot_batchf_fast ),
    TEST_DESC_EX( FMT_REAL|FMT_FLOAT32, TE_DIM_NUM_4, TE_ALIGN_YES, &te_loadFxn_vec_dot_batch, &te_processFxn_vec_dot_batch,&vec_dot_batchf_api ) },
  {
    FUNC_LIST( (tTestEngTarget)&vec_dot16x16 ),
    TEST_DESC( FMT_REAL|FMT_FRACT16, TE_DIM_NUM_2, TE_ALIGN_NO, &te_loadFxn_vXvYsZ64, &te_processFxn_sZ64vXvY ) },
  {
    FUNC_LIST( (tTestEngTarget)&vec_dot32x16 ),
    TEST_DESC( FMT_REAL|FMT_FRACT32, TE_DIM_NUM_2, TE_ALIGN_NO, &te_loadFxn_vXvY16sZ64, &te_processFxn_sZ64vXvY ) },
  {
    FUNC_LIST( (tTestEngTarget)&vec_dot32x32 ),
    TEST_DESC(FMT_REAL | FMT_FRACT32, TE_DIM_NUM_2, TE_ALIGN_NO, &te_loadFxn_vXvYsZ64, &te_processFxn_sZ64vXvY) },
  {
    FUNC_LIST( (tTestEngTarget)&vec_dot32x32_fast ),
    TEST_DESC(FMT_REAL | FMT_FRACT32, TE_DIM_NUM_2, TE_ALIGN_YES, &te_loadFxn_vXvYsZ64, &te_processFxn_sZ64vXvY) },
  {
    FUNC_LIST( (tTestEngTarget)&vec_dot32x16_fast ),
    TEST_DESC( FMT_REAL|FMT_FRACT32, TE_DIM_NUM_2, TE_ALIGN_YES, &te_loadFxn_vXvY16sZ64, &te_processFxn_sZ64vXvY ) },
  {
    FUNC_LIST( (tTestEngTarget)&vec_dot16x16_fast ),
    TEST_DESC( FMT_REAL|FMT_FRACT16, TE_DIM_NUM_2, TE_ALIGN_YES, &te_loadFxn_vXvYsZ32, &te_processFxn_sZ32vXvY ) },
  {
    FUNC_LIST( (tTestEngTarget)&vec_dot64x32),
    TEST_DESC( FMT_REAL|FMT_INT64, TE_DIM_NUM_2, TE_ALIGN_NO, &te_loadFxn_vXvY32sZ64, &te_processFxn_sZ64vXvY ) },
  {
    FUNC_LIST( (tTestEngTarget)&vec_dot64x32_fast),
    TEST_DESC( FMT_REAL|FMT_INT64, TE_DIM_NUM_2, TE_ALIGN_YES, &te_loadFxn_vXvY32sZ64, &te_processFxn_sZ64vXvY ) },
  {
    FUNC_LIST( (tTestEngTarget)&vec_dot64x64, 
                (tTestEngTarget)&vec_dot64x64i),
    TEST_DESC( FMT_REAL|FMT_INT64, TE_DIM_NUM_2, TE_ALIGN_NO, &te_loadFxn_vXvYsZ64, &te_processFxn_sZ64vXvY ) },
   {
    FUNC_LIST(  (tTestEngTarget)&vec_dot64x64_fast,
                (tTestEngTarget)&vec_dot64x64i_fast),
    TEST_DESC( FMT_REAL|FMT_INT64, TE_DIM_NUM_2, TE_ALIGN_YES, &te_loadFxn_vXvYsZ64, &te_processFxn_sZ64vXvY ) },
  {
    FUNC_LIST( (tTestEngTarget)&vec_add32x32 ),
    TEST_DESC( FMT_REAL|FMT_FRACT32, TE_DIM_NUM_2, TE_ALIGN_NO, &te_loadFxn_vXvYvZ, &te_processFxn_vZvXvY ) },
  {
    FUNC_LIST( (tTestEngTarget)&vec_add32x32_fast),
    TEST_DESC( FMT_REAL|FMT_FRACT32, TE_DIM_NUM_2, TE_ALIGN_YES, &te_loadFxn_vXvYvZ, &te_processFxn_vZvXvY ) },
  {
    FUNC_LIST( (tTestEngTarget)&vec_add16x16 ),
    TEST_DESC( FMT_REAL|FMT_FRACT16, TE_DIM_NUM_2, TE_ALIGN_NO, &te_loadFxn_vXvYvZ, &te_processFxn_vZvXvY ) },
  {
    FUNC_LIST( (tTestEngTarget)&vec_add16x16_fast ),
    TEST_DESC( FMT_REAL|FMT_FRACT16, TE_DIM_NUM_2, TE_ALIGN_YES, &te_loadFxn_vXvYvZ, &te_processFxn_vZvXvY ) },
  {
    FUNC_LIST( (tTestEngTarget)&vec_power16x16 ),
    TEST_DESC( FMT_REAL|FMT_FRACT16, TE_DIM_NUM_2, TE_ALIGN_NO, &te_loadFxn_vXsY32sZ64, &te_processFxn_vXsY32sZ64 ) },
  {
    FUNC_LIST( (tTestEngTarget)&vec_power16x16_fast ),
    TEST_DESC( FMT_REAL|FMT_FRACT16, TE_DIM_NUM_2, TE_ALIGN_YES, &te_loadFxn_vXsY32sZ64, &te_processFxn_vXsY32sZ64 ) },
  {
    FUNC_LIST( (tTestEngTarget)&vec_power32x32 ),
    TEST_DESC( FMT_REAL|FMT_FRACT32, TE_DIM_NUM_2, TE_ALIGN_NO, &te_loadFxn_vXsY32sZ64, &te_processFxn_vXsY32sZ64 ) },
  {
    FUNC_LIST( (tTestEngTarget)&vec_power32x32_fast ),
    TEST_DESC( FMT_REAL|FMT_FRACT32, TE_DIM_NUM_2, TE_ALIGN_YES, &te_loadFxn_vXsY32sZ64, &te_processFxn_vXsY32sZ64 ) },
  { 
    FUNC_LIST( (tTestEngTarget)&vec_min16x16, (tTestEngTarget)&vec_max16x16),
    TEST_DESC( FMT_REAL|FMT_FRACT16, TE_DIM_NUM_2, TE_ALIGN_NO, &te_loadFxn_vXsZ  , &te_processFxn_vXsZ ) },
  { 
    FUNC_LIST( (tTestEngTarget)&vec_min16x16_fast, (tTestEngTarget)&vec_max16x16_fast),
    TEST_DESC( FMT_REAL|FMT_FRACT16, TE_DIM_NUM_2, TE_ALIGN_YES, &te_loadFxn_vXsZ  , &te_processFxn_vXsZ ) },
  { 
    FUNC_LIST( (tTestEngTarget)&vec_min32x32,(tTestEngTarget)&vec_max32x32),
    TEST_DESC( FMT_REAL|FMT_FRACT32, TE_DIM_NUM_2, TE_ALIGN_NO, &te_loadFxn_vXsZ  , &te_processFxn_vXsZ ) },
  { 
    FUNC_LIST( (tTestEngTarget)&vec_min32x32_fast,(tTestEngTarget)&vec_max32x32_fast),
    TEST_DESC( FMT_REAL|FMT_FRACT32, TE_DIM_NUM_2, TE_ALIGN_YES, &te_loadFxn_vXsZ  , &te_processFxn_vXsZ ) },
  { 
    FUNC_LIST( (tTestEngTarget)&vec_scale32x32, 
                (tTestEngTarget)&vec_shift32x32),
    TEST_DESC( FMT_REAL|FMT_FRACT32, TE_DIM_NUM_2, TE_ALIGN_NO, &te_loadFxn_vXsYvZ, &te_processFxn_vZvXsY ) },
  { 
    FUNC_LIST( (tTestEngTarget)&vec_scale32x32_fast,
                (tTestEngTarget)&vec_shift32x32_fast),
    TEST_DESC( FMT_REAL|FMT_FRACT32, TE_DIM_NUM_2, TE_ALIGN_YES, &te_loadFxn_vXsYvZ, &te_processFxn_vZvXsY ) },
  { 
    FUNC_LIST( (tTestEngTarget)&vec_scale16x16), 
    TEST_DESC( FMT_REAL|FMT_FRACT16, TE_DIM_NUM_2, TE_ALIGN_NO, &te_loadFxn_vXsYvZ, &te_processFxn_vZvXsY ) },
  { 
    FUNC_LIST( (tTestEngTarget)&vec_scale16x16_fast), 
    TEST_DESC( FMT_REAL|FMT_FRACT16, TE_DIM_NUM_2, TE_ALIGN_YES, &te_loadFxn_vXsYvZ, &te_processFxn_vZvXsY ) },
  { 
    FUNC_LIST( (tTestEngTarget)&vec_shift16x16), 
    TEST_DESC( FMT_REAL|FMT_FRACT16, TE_DIM_NUM_2, TE_ALIGN_NO, &te_loadFxn_vXsY32vZ, &te_processFxn_vZvXsY32 ) },
  { 
    FUNC_LIST( (tTestEngTarget)&vec_shift16x16_fast), 
    TEST_DESC( FMT_REAL|FMT_FRACT16, TE_DIM_NUM_2, TE_ALIGN_YES, &te_loadFxn_vXsY32vZ, &te_processFxn_vZvXsY32 ) },
  { 
    FUNC_LIST( (tTestEngTarget)&vec_bexp32),
    TEST_DESC( FMT_REAL|FMT_FRACT32, TE_DIM_NUM_2, TE_ALIGN_NO, &te_loadFxn_vXsZ  , &te_processFxn_vXsZ ) },
  { 
    FUNC_LIST( (tTestEngTarget)&vec_bexp32_fast ),
    TEST_DESC( FMT_REAL|FMT_FRACT32, TE_DIM_NUM_2, TE_ALIGN_YES, &te_loadFxn_vXsZ  , &te_processFxn_vXsZ ) },
  { 
    FUNC_LIST( (tTestEngTarget)&vec_bexp16),
    TEST_DESC( FMT_REAL|FMT_FRACT16, TE_DIM_NUM_2, TE_ALIGN_NO, &te_loadFxn_vXsZ  , &te_processFxn_vXsZ ) },
  { 
    FUNC_LIST( (tTestEngTarget)&vec_bexp16_fast),
    TEST_DESC( FMT_REAL|FMT_FRACT16, TE_DIM_NUM_2, TE_ALIGN_YES, &te_loadFxn_vXsZ  , &te_processFxn_vXsZ ) },
  {
    FUNC_LIST( (tTestEngTarget)&scl_bexp16),
    TEST_DESC( FMT_FRACT16, TE_DIM_NUM_2, TE_ALIGN_NO, &te_loadFxn_vXsZ32, &processFxn_scl_vXvZ32 ) },
  {
    FUNC_LIST((tTestEngTarget)&scl_bexp32 ),
    TEST_DESC( FMT_FRACT32, TE_DIM_NUM_2, TE_ALIGN_NO, &te_loadFxn_vXsZ, &processFxn_scl_vXvZ32 ) },
  { 
    FUNC_LIST( NULL ), TEST_DESC(  0, 0, 0, NULL, NULL ) } /* End of table */
};

/* Perform all functional tests for Vector Mathematics API functions. */
int func_vector1(int isFull, int isVerbose, int breakOnError)
{
	int res = 1;

#define DO_TEST(fxn, seqFile)                                                                   \
    if ( res || !breakOnError ) res &= ( 0 != te_Exec( testDefTbl,                              \
                                                       sizeof(testDefTbl)/sizeof(testDefTbl[0]),\
                                                       MAX_FUNC_NUM,                            \
                                                       (tTestEngTarget)(fxn), "vec1/" seqFile,   \
                                                       isFull, isVerbose, breakOnError ) )
	
  DO_TEST( &vec_dot16x16        , "vec_dot16x16.seq"        );
#if 0//HiFi3/3z API
  DO_TEST( &vec_dot24x24        , "vec_dot24x24.seq"        );
#endif
  DO_TEST( &vec_dot32x16        , "vec_dot32x16.seq"        );
  DO_TEST( &vec_dot32x32        , "vec_dot32x32.seq"        );
  DO_TEST( &vec_dot64x32        , "vec_dot64x32.seq"        );
  DO_TEST( &vec_dot64x64        , "vec_dot64x64.seq"        );
  DO_TEST( &vec_dot64x64i       , "vec_dot64x64i.seq"       );
  DO_TEST( &vec_dot16x16_fast   , "vec_dot16x16_fast.seq"   );
#if 0//HiFi3/3z API
  DO_TEST( &vec_dot24x24_fast   , "vec_dot24x24_fast.seq"   );
#endif
  DO_TEST( &vec_dot32x16_fast   , "vec_dot32x16_fast.seq"   );
  DO_TEST( &vec_dot32x32_fast   , "vec_dot32x32_fast.seq"   );
  DO_TEST( &vec_dot64x32_fast   , "vec_dot64x32_fast.seq"   );
  DO_TEST( &vec_dot64x64_fast   , "vec_dot64x64_fast.seq"   );
  DO_TEST( &vec_dot64x64i_fast  , "vec_dot64x64i_fast.seq"  );
  DO_TEST( &vec_dot_batch8x8          , "vec_dot_batch8x8.seq"          );
  DO_TEST( &vec_dot_batch8x8_fast     , "vec_dot_batch8x8_fast.seq"     );
  DO_TEST( &vec_dot_batch8x16         , "vec_dot_batch8x16.seq"         );
  DO_TEST( &vec_dot_batch8x16_fast    , "vec_dot_batch8x16_fast.seq"    );
  DO_TEST( &vec_dot_batch16x16        , "vec_dot_batch16x16.seq"        );
  DO_TEST( &vec_dot_batch16x16_fast   , "vec_dot_batch16x16_fast.seq"   );
  DO_TEST( &vec_add32x32        , "vec_add32x32.seq"        );
#if 0 //HiFi3/3z API
  DO_TEST( &vec_add24x24        , "vec_add24x24.seq"        );
#endif
  DO_TEST( &vec_add16x16        , "vec_add16x16.seq"        );
  DO_TEST( &vec_add32x32_fast   , "vec_add32x32_fast.seq"   );
#if 0 //HiFi3/3z API
  DO_TEST( &vec_add24x24_fast   , "vec_add24x24_fast.seq"   );
#endif
  DO_TEST( &vec_add16x16_fast   , "vec_add16x16_fast.seq"   );
  DO_TEST( &vec_power32x32      , "vec_power32x32.seq"      );
#if 0 //HiFi3/3z API
  DO_TEST( &vec_power24x24      , "vec_power24x24.seq"      );
#endif
  DO_TEST( &vec_power16x16      , "vec_power16x16.seq"      );
  DO_TEST( &vec_power32x32_fast , "vec_power32x32_fast.seq" );
#if 0//HiFi3/3z API
  DO_TEST( &vec_power24x24_fast , "vec_power24x24_fast.seq" );
#endif
  DO_TEST( &vec_power16x16_fast , "vec_power16x16_fast.seq" );
  DO_TEST( &vec_shift32x32      , "vec_shift32x32.seq"      );
#if 0 //HiFi3/3z API
  DO_TEST( &vec_shift24x24      , "vec_shift24x24.seq"      );
#endif
  DO_TEST( &vec_shift16x16      , "vec_shift16x16.seq"      );
#if 0 //HiFi3/3z API
  DO_TEST( &vec_scale32x24      , "vec_scale32x24.seq"      );
  DO_TEST( &vec_scale24x24      , "vec_scale24x24.seq"      );
#endif
  DO_TEST( &vec_scale16x16      , "vec_scale16x16.seq"      );
  DO_TEST( &vec_scale32x32      , "vec_scale32x32.seq"      );
  DO_TEST( &vec_shift32x32_fast , "vec_shift32x32_fast.seq" );
#if 0//HiFi3/3z API
  DO_TEST( &vec_shift24x24_fast , "vec_shift24x24_fast.seq" );
#endif
  DO_TEST( &vec_shift16x16_fast , "vec_shift16x16_fast.seq" );
#if 0//HiFi3/3z API
  DO_TEST( &vec_scale32x24_fast , "vec_scale32x24_fast.seq" );
  DO_TEST( &vec_scale24x24_fast , "vec_scale24x24_fast.seq" );
#endif
  DO_TEST( &vec_scale16x16_fast , "vec_scale16x16_fast.seq" );
  DO_TEST( &vec_scale32x32_fast , "vec_scale32x32_fast.seq" );
  DO_TEST( &vec_min32x32        , "vec_min32x32.seq"        );
#if 0 //HiFi3/3z API
  DO_TEST( &vec_min24x24        , "vec_min24x24.seq"        );
#endif
  DO_TEST( &vec_min16x16        , "vec_min16x16.seq"        );
  DO_TEST( &vec_max32x32        , "vec_max32x32.seq"        );
#if 0 //HiFi3/3z API
  DO_TEST( &vec_max24x24        , "vec_max24x24.seq"        );
#endif
  DO_TEST( &vec_max16x16        , "vec_max16x16.seq"        );
  DO_TEST( &vec_min32x32_fast   , "vec_min32x32_fast.seq"   );
#if 0 //HiFi3/3z API
  DO_TEST( &vec_min24x24_fast   , "vec_min24x24_fast.seq"   );
#endif
  DO_TEST( &vec_min16x16_fast   , "vec_min16x16_fast.seq"   );
  DO_TEST( &vec_max32x32_fast   , "vec_max32x32_fast.seq"   );
#if 0 //HiFi3/3z API
  DO_TEST( &vec_max24x24_fast   , "vec_max24x24_fast.seq"   );
#endif
  DO_TEST( &vec_max16x16_fast   , "vec_max16x16_fast.seq"   );
  DO_TEST( &scl_bexp32             , "scl_bexp32.seq"       );
#if 0 //HiFi3/3z API
  DO_TEST( &scl_bexp24             , "scl_bexp24.seq"       );
#endif
  DO_TEST( &scl_bexp16             , "scl_bexp16.seq"       );
  DO_TEST( &vec_bexp32             , "vec_bexp32.seq"       );
#if 0 //HiFi3/3z API
  DO_TEST( &vec_bexp24             , "vec_bexp24.seq"       );
#endif
  DO_TEST( &vec_bexp16             , "vec_bexp16.seq"       );
  DO_TEST( &vec_bexp32_fast        , "vec_bexp32_fast.seq"  );
#if 0 //HiFi3/3z API
  DO_TEST( &vec_bexp24_fast        , "vec_bexp24_fast.seq"  );
#endif
  DO_TEST( &vec_bexp16_fast        , "vec_bexp16_fast.seq"  );

	return (res);
}

