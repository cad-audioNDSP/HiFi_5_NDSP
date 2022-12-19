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
 * Test-engine add-on for older FIR API
 */
#ifndef TESTENG__FIR_OLD_H__
#define TESTENG__FIR_OLD_H__

/* Cross-platform data type definitions. */
#include "types.h"
/* Test engine API. */
#include "testeng.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef size_t tFirOldDecimaFxnAlloc(int D,int M);
typedef void * tFirOldDecimaFxnInit (void * objmem, int D,int M, const void* h);
typedef size_t tFirOldextIrFxnAlloc(int M, int extIR);
typedef void * tFirOldextIrFxnInit (void * objmem, int M, int extIR, const void* h);
typedef size_t tFirOldFxnAlloc(int M);
typedef void * tFirOldFxnInit (void * objmem, int M, const void* h);
typedef void   tFirOldFxnProcess (void *handle, void * y, const void * x , int N);
typedef void   tFirOldPolyphaseFxnProcess (void *handle, void * y, const void * x);
typedef void   tFirOldPolyphaseFxnProcess_lsh (void *handle, void * y, const void * x, int lsh);
typedef size_t tFirOldextIrFxnAllocExtIr(int M);
typedef void   tFirOldextIrFxnCopyExtIr(void* pExtIr, const void* ir, int M);

typedef struct
{
    tFirOldDecimaFxnAlloc   *alloc;
    tFirOldDecimaFxnInit    *init;
    tFirOldFxnProcess       *process;
}
tFirOldDecimaDescr;

typedef struct
{
    tFirOldextIrFxnAlloc      *alloc;
    tFirOldextIrFxnInit       *init;
    tFirOldFxnProcess         *process;
    tFirOldextIrFxnAllocExtIr *allocExtIr;
    tFirOldextIrFxnCopyExtIr  *copyExtIr;
}
tFirOldextIrDescr;

typedef struct
{
    tFirOldFxnAlloc   *alloc;
    tFirOldFxnInit    *init;
    tFirOldFxnProcess *process;
    tFirOldextIrFxnAllocExtIr *allocExtIr;
    tFirOldextIrFxnCopyExtIr  *copyExtIr;
}
tFirOldGenericDescr;

typedef union
{
    tFirOldGenericDescr g;
    tFirOldextIrDescr   e;
    tFirOldDecimaDescr  d;
}
tFirOldDescr;

/* function reads IR from file and creates FIR structure. returns 0 if failed */
int te_create_fir_old(tTestEngContext * context);

/* function destroys FIR structure, returns 0 if failed */
int te_destroy_fir_old(tTestEngContext * context);

/* 
   Allocate vectors and load the data set for FIR:
*  vector X (in), vector Z (out) */
int te_loadFxn_fir_old(tTestEngContext * context);

/* Apply FIR function to the test case data set.
*  vector X (in), vector Z (out) */
void te_processFxn_fir_old(tTestEngContext * context);

typedef size_t tFirStereoDecimaFxnAlloc(int D,int M);
typedef void * tFirStereoDecimaFxnInit (void * objmem, int D,int M, const void* hl, const void* hr);
typedef size_t tFirStereoextIrFxnAlloc(int M, int extIR);
typedef void * tFirStereoextIrFxnInit (void * objmem, int M, int extIR, const void* hl, const void* hr);
typedef size_t tFirStereoFxnAlloc(int M);
typedef void * tFirStereoFxnInit (void * objmem, int M, const void* hl, const void* hr);
typedef void   tFirStereoFxnProcess (void *handle, void * y, const void * x , int N);
typedef size_t tFirStereoextIrFxnAllocExtIr(int M);
typedef void   tFirStereoextIrFxnCopyExtIr(void* pExtIrl, void* pExtIrr, const void* irl, const void* irr, int M);

typedef struct
{
    tFirStereoDecimaFxnAlloc   *alloc;
    tFirStereoDecimaFxnInit    *init;
    tFirStereoFxnProcess       *process;
}
tFirStereoDecimaDescr;

typedef struct
{
    tFirStereoextIrFxnAlloc      *alloc;
    tFirStereoextIrFxnInit       *init;
    tFirStereoFxnProcess         *process;
    tFirStereoextIrFxnAllocExtIr *allocExtIr;
    tFirStereoextIrFxnCopyExtIr  *copyExtIr;
}
tFirStereoextIrDescr;

typedef struct
{
    tFirStereoFxnAlloc   *alloc;
    tFirStereoFxnInit    *init;
    tFirStereoFxnProcess *process;
    tFirStereoextIrFxnAllocExtIr *allocExtIr;
    tFirStereoextIrFxnCopyExtIr  *copyExtIr;
}
tFirStereoGenericDescr;

typedef union
{
    tFirStereoGenericDescr g;
    tFirStereoextIrDescr   e;
    tFirStereoDecimaDescr  d;
}
tFirStereoDescr;

/* function reads IR from file and creates FIR structure. returns 0 if failed */
int te_create_fir_stereo(tTestEngContext * context);

/* function destroys FIR structure, returns 0 if failed */
int te_destroy_fir_stereo(tTestEngContext * context);

/* 
   Allocate vectors and load the data set for FIR:
*  vector X (in), vector Z (out) */
int te_loadFxn_fir_stereo(tTestEngContext * context);

/* Apply FIR function to the test case data set.
*  vector X (in), vector Z (out) */
void te_processFxn_fir_stereo(tTestEngContext * context);

/* management of external IR */
size_t bkfir16x16_allocExtIr(int M);
size_t bkfir32x16_allocExtIr(int M);
size_t bkfir32x32_allocExtIr(int M);
size_t bkfir24x24_allocExtIr(int M);
size_t bkfir24x24p_allocExtIr(int M);
size_t cxfir16x16_allocExtIr(int M);
size_t cxfir24x24_allocExtIr(int M);
size_t cxfir32x16_allocExtIr(int M);
size_t cxfir32x32_allocExtIr(int M);
size_t cxfir32x32ep_allocExtIr(int M);
size_t bkfirf_allocExtIr(int M);
size_t cxfirf_allocExtIr(int M);
size_t stereo_bkfir16x16_allocExtIr(int M);
size_t stereo_bkfir32x32_allocExtIr(int M);
size_t stereo_bkfirf_allocExtIr(int M);

/* left padding 2 bytes, right padding 6 bytes, inverted order */
void   bkfir16x16_copyExtIr(int16_t* pExtIr, const int16_t* ir, int M);
/* left padding 2 or 10 bytes, right padding 6 bytes, inverted order */
void   bkfir32x16_copyExtIr(int16_t* pExtIr, const int16_t* ir, int M);
/* left padding 4 bytes, right padding 12 bytes, inverted order */
void   bkfir32x32_copyExtIr(int32_t* pExtIr, const int32_t* ir, int M);
/* left padding 3 bytes, right padding 9 bytes, inverted order */
void bkfir24x24_copyExtIr(f24* pExtIr, const f24* ir, int M);
/* left padding ((-M&4)+1)*3 bytes, right padding 7, inverted order */
void   bkfir24x24p_copyExtIr(uint8_t* pExtIr, const uint8_t* ir, int M);
/* left padding 4 bytes, inverted order, complex coefficients */
void   cxfir16x16_copyExtIr(int16_t* pExtIr, const int16_t* ir, int M);
/* no padding, inverted order, complex coefficients */
void   cxfir24x24_copyExtIr(f24* pExtIr, const f24* ir, int M);
/* left padding 4 bytes, right padding 4 bytes, inverted order, complex coefficients */
void   cxfir32x16_copyExtIr(int16_t* pExtIr, const int16_t* ir, int M);
/* no padding, inverted order, complex conjugated coefficients */
void   cxfir32x32_copyExtIr(int32_t* pExtIr, const int32_t* ir, int M);
/* no padding, inverted order*/
void   cxfir32x32ep_copyExtIr(int32_t* pExtIr, const int32_t* ir, int M);
/* no padding, direct order */
void   bkfirf_copyExtIr(float32_t* pExtIr, const float32_t* ir, int M);
/* no padding, direct order */
void   cxfirf_copyExtIr(complex_float* pExtIr, const complex_float* ir, int M);
/* left padding 2 bytes, right padding 6 bytes, inverted order, first left then right */
void   stereo_bkfir16x16_copyExtIr(int16_t* pExtIrl, int16_t* pExtIrr, const int16_t* irl, const int16_t* irr, int M);
/* left padding 4 bytes, right padding 12 bytes, inverted order, first left then right */
void   stereo_bkfir32x32_copyExtIr(int32_t* pExtIrl, int32_t* pExtIrr, const int32_t* irl, const int32_t* irr, int M);
/* no padding, direct order, first left then right */
void   stereo_bkfirf_copyExtIr(float32_t* pExtIrl, float32_t* pExtIrr, const float32_t* irl, const float32_t* irr, int M);

#ifdef __cplusplus
};
#endif

#endif

