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
 * Test procedures for Cholesky decomposition (block data)
 */
#ifndef __TEST_CHOL_H__
#define __TEST_CHOL_H__
#include "NatureDSP_types.h"
#include "NatureDSP_Math.h"
#include "NatureDSP_Signal_matinv.h"
#include "malloc16.h"
#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
/* Test Engine API. */
#include "testeng.h"

#define MIN(a,b)   ( (a)<(b) ? (a) : (b) )
#define MAX(a,b)   ( (a)>(b) ? (a) : (b) )

/* chol API flags */
#define CHOLAPI_FIX_SZ        0x0001 /* flag for fixed-size API like chol8x8n            */
#define CHOLAPI_FLT_SZ        0x0002 /* flag for floating-size API like cholmxnn         */
#define CHOLAPI_SZTYPE_MASK   0x0003

#define CHOLAPI_PACKED        0x0004 /* flag for block-ordered functions like cholmxnn  */
#define CHOLAPI_STREAM        0x0008 /* flag for stream-ordered functions like cholmxns */
#define CHOLAPI_LAYOUT_MASK   0x000C
#define CHOLAPI_USESIGMA2     0x0010 /* flag for functions that use sigma2 regularization term */

void cholmxnn_32b_unpackD(void * dst, const void * src, int N, int L);
void cholmxnn_32b_packD(void * dst, const void * src, int N, int L);

/* API definition structure. */
typedef struct 
{
    size_t (*cholScratchSz)(); 
    size_t (*fwdScratchSz)(); 
    size_t (*bkwScratchSz)();        
    size_t (*mmseScratchSz)(); 
    size_t (*preprocessScratchSz)(); 
    size_t (*pinvScratchSz)(); 

    /* Helping function, convert reciprocal of diagonal elements from the format used in SEQ-file into the format used in function */
    void (*te_cholFxn_unpackD)(void * dst, const void * src, int N, int L);
    /* Helping function, convert reciprocal of diagonal elements from the format used in function into the format used in SEQ-file */
    void (*te_cholFxn_packD)(void * dst, const void * src, int N, int L);

} tCholApi;

static const tCholApi cmatchol4x4_32x32_Api   =    { cmatcholdecomp4x4_32x32_getScratchSize,   cmatcholfwdsubst4x4_32x32_getScratchSize  ,cmatcholbkwsubst4x4_32x32_getScratchSize   ,cmatcholmmsesolver4x4_32x32_getScratchSize   ,cmatcholpreprocess4x4_32x32_getScratchSize  ,cmatcholpseudoinv4x4_32x32_getScratchSize  , cholmxnn_32b_unpackD, cholmxnn_32b_packD };
static const tCholApi cmatchol6x6_32x32_Api   =    { cmatcholdecomp6x6_32x32_getScratchSize,   cmatcholfwdsubst6x6_32x32_getScratchSize  ,cmatcholbkwsubst6x6_32x32_getScratchSize   ,cmatcholmmsesolver6x6_32x32_getScratchSize   ,cmatcholpreprocess6x6_32x32_getScratchSize  ,cmatcholpseudoinv6x6_32x32_getScratchSize  , cholmxnn_32b_unpackD, cholmxnn_32b_packD };
static const tCholApi cmatchol8x8_32x32_Api   =    { cmatcholdecomp8x8_32x32_getScratchSize,   cmatcholfwdsubst8x8_32x32_getScratchSize  ,cmatcholbkwsubst8x8_32x32_getScratchSize   ,cmatcholmmsesolver8x8_32x32_getScratchSize   ,cmatcholpreprocess8x8_32x32_getScratchSize  ,cmatcholpseudoinv8x8_32x32_getScratchSize  , cholmxnn_32b_unpackD, cholmxnn_32b_packD };
static const tCholApi cmatchol10x10_32x32_Api =    { cmatcholdecomp10x10_32x32_getScratchSize, cmatcholfwdsubst10x10_32x32_getScratchSize,cmatcholbkwsubst10x10_32x32_getScratchSize ,cmatcholmmsesolver10x10_32x32_getScratchSize ,cmatcholpreprocess10x10_32x32_getScratchSize,cmatcholpseudoinv10x10_32x32_getScratchSize, cholmxnn_32b_unpackD, cholmxnn_32b_packD };

static const tCholApi  matchol4x4_32x32_Api   =    {  matcholdecomp4x4_32x32_getScratchSize,   matcholfwdsubst4x4_32x32_getScratchSize   ,matcholbkwsubst4x4_32x32_getScratchSize    ,matcholmmsesolver4x4_32x32_getScratchSize    ,cmatcholpreprocess4x4_32x32_getScratchSize  , matcholpseudoinv4x4_32x32_getScratchSize  , cholmxnn_32b_unpackD, cholmxnn_32b_packD };
static const tCholApi  matchol6x6_32x32_Api   =    {  matcholdecomp6x6_32x32_getScratchSize,   matcholfwdsubst6x6_32x32_getScratchSize   ,matcholbkwsubst6x6_32x32_getScratchSize    ,matcholmmsesolver6x6_32x32_getScratchSize    ,cmatcholpreprocess6x6_32x32_getScratchSize  , matcholpseudoinv6x6_32x32_getScratchSize  , cholmxnn_32b_unpackD, cholmxnn_32b_packD };
static const tCholApi  matchol8x8_32x32_Api   =    {  matcholdecomp8x8_32x32_getScratchSize,   matcholfwdsubst8x8_32x32_getScratchSize   ,matcholbkwsubst8x8_32x32_getScratchSize    ,matcholmmsesolver8x8_32x32_getScratchSize    ,cmatcholpreprocess8x8_32x32_getScratchSize  , matcholpseudoinv8x8_32x32_getScratchSize  , cholmxnn_32b_unpackD, cholmxnn_32b_packD };
static const tCholApi  matchol10x10_32x32_Api =    {  matcholdecomp10x10_32x32_getScratchSize, matcholfwdsubst10x10_32x32_getScratchSize ,matcholbkwsubst10x10_32x32_getScratchSize  ,matcholmmsesolver10x10_32x32_getScratchSize  ,cmatcholpreprocess10x10_32x32_getScratchSize, matcholpseudoinv10x10_32x32_getScratchSize, cholmxnn_32b_unpackD, cholmxnn_32b_packD };

static const tCholApi chol4x4x1nfApi =    { cmatcholdecomp4x4f_getScratchSize,   cmatcholfwdsubst4x4f_getScratchSize,   cmatcholbkwsubst4x4f_getScratchSize,   cmatcholmmsesolver4x4f_getScratchSize,   cmatcholpreprocess4x4f_getScratchSize,   cmatcholpseudoinv4x4f_getScratchSize };
static const tCholApi chol6x6x1nfApi =    { cmatcholdecomp6x6f_getScratchSize,   cmatcholfwdsubst6x6f_getScratchSize,   cmatcholbkwsubst6x6f_getScratchSize,   cmatcholmmsesolver6x6f_getScratchSize,   cmatcholpreprocess6x6f_getScratchSize,   cmatcholpseudoinv6x6f_getScratchSize };
static const tCholApi chol8x8x1nfApi =    { cmatcholdecomp8x8f_getScratchSize,   cmatcholfwdsubst8x8f_getScratchSize,   cmatcholbkwsubst8x8f_getScratchSize,   cmatcholmmsesolver8x8f_getScratchSize,   cmatcholpreprocess8x8f_getScratchSize,   cmatcholpseudoinv8x8f_getScratchSize };
static const tCholApi chol10x10x1nfApi =  { cmatcholdecomp10x10f_getScratchSize, cmatcholfwdsubst10x10f_getScratchSize, cmatcholbkwsubst10x10f_getScratchSize, cmatcholmmsesolver10x10f_getScratchSize, cmatcholpreprocess10x10f_getScratchSize, cmatcholpseudoinv10x10f_getScratchSize };
static const tCholApi rchol4x4x1nfApi =   { matcholdecomp4x4f_getScratchSize,    matcholfwdsubst4x4f_getScratchSize,    matcholbkwsubst4x4f_getScratchSize,    matcholmmsesolver4x4f_getScratchSize,    matcholpreprocess4x4f_getScratchSize,    matcholpseudoinv4x4f_getScratchSize };
static const tCholApi rchol6x6x1nfApi =   { matcholdecomp6x6f_getScratchSize,    matcholfwdsubst6x6f_getScratchSize,    matcholbkwsubst6x6f_getScratchSize,    matcholmmsesolver6x6f_getScratchSize,    matcholpreprocess6x6f_getScratchSize,    matcholpseudoinv6x6f_getScratchSize };
static const tCholApi rchol8x8x1nfApi =   { matcholdecomp8x8f_getScratchSize,    matcholfwdsubst8x8f_getScratchSize,    matcholbkwsubst8x8f_getScratchSize,    matcholmmsesolver8x8f_getScratchSize,    matcholpreprocess8x8f_getScratchSize,    matcholpseudoinv8x8f_getScratchSize };
static const tCholApi rchol10x10x1nfApi = { matcholdecomp10x10f_getScratchSize,  matcholfwdsubst10x10f_getScratchSize,  matcholbkwsubst10x10f_getScratchSize,  matcholmmsesolver10x10f_getScratchSize,  matcholpreprocess10x10f_getScratchSize,  matcholpseudoinv10x10f_getScratchSize };

typedef struct
{
    int NA;
    int NR;
    int ND;
    int NY;
    int NB;
    int NX;
    int NQ;
    int NSigma2;
} tCholNelems;

/* Allocate vectors and load the data set for Cholesky decomposition functions */
int te_loadFxn_chol( tTestEngContext * context );

/* Allocate vectors and load the data set for cholnf functions */
int te_loadFxn_cholf( tTestEngContext * context );

/* processing function for Cholesky decomposition functions */
void te_processFxn_chol( tTestEngContext * context );

/* Allocate vectors and load the data set for Cholesky backward recursion functions */
int te_loadFxn_cholbkw( tTestEngContext * context );

/* processing function for Cholesky backward recursion functions */
void te_processFxn_cholbkw( tTestEngContext * context );

/* Allocate vectors and load the data set for Cholesky forward substitution functions */
int te_loadFxn_cholfwd( tTestEngContext * context );

/* processing function for Cholesky forward substitution functions */
void te_processFxn_cholfwd( tTestEngContext * context );

/* Allocate vectors and load the data set for Cholesky MMSE functions */
int te_loadFxn_cholmmse( tTestEngContext * context );

/* processing function for Cholesky MMSE functions */
void te_processFxn_cholmmse(tTestEngContext * context);

/* Allocate vectors and load the data set for cholfwdf functions */
int te_loadFxn_cholfwdf( tTestEngContext * context );

/* Allocate vectors and load the data set for cholbkwf functions */
int te_loadFxn_cholbkwf( tTestEngContext * context );

/* Allocate vectors and load the data set for cholmmsef functions */
int te_loadFxn_cholmmsef( tTestEngContext * context );

/* Allocate vectors and load the data set for cholpreprocessf functions */
int te_loadFxn_cholpreprocessf(tTestEngContext * context);

/* Allocate vectors and load the data set for cholpinvf functions */
int te_loadFxn_cholpinvf(tTestEngContext * context);

/* processing function for cholf */
void te_processFxn_cholf( tTestEngContext * context );

/* processing function for cholfwdf */
void te_processFxn_cholfwdf( tTestEngContext * context );

/* processing function for cholbkwf */
void te_processFxn_cholbkwf( tTestEngContext * context );

/* processing function for cholmmsef */
void te_processFxn_cholmmsef( tTestEngContext * context );

/* processing function for cholpreprocessf */
void te_processFxn_cholpreprocessf(tTestEngContext * context);

/* processing function for cholpinvf */
void te_processFxn_cholpinvf(tTestEngContext * context);

#define MAX_FUNC_NUM 4
/* Initializer for a function pointer array, appends NULL to a sequence of pointers. */
#define FUNC_LIST(...) { __VA_ARGS__, NULL }
/* Initializer for a test description structure. */
#define TEST_DESC( fmt, extraParam, extraPtr,argNum, align, loadFxn, procFxn ) { (fmt),extraParam,extraPtr,(argNum),(align),NULL,NULL,(loadFxn),(procFxn) }

#endif  /* __TEST_CHOL_H__ */
