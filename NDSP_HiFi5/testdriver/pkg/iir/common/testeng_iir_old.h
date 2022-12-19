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
 * Test-engine add-on for old IIR categories 
 */
#ifndef TESTENG_IIR_OLD_H__
#define TESTENG_IIR_OLD_H__

/* Cross-platform data type definitions. */
#include "types.h"
/* Test engine API. */
#include "testeng.h"
#include "testeng_iir.h"

typedef size_t tIirOldFxnAlloc (int);
typedef void * tIirOldFxnInit (void * objmem, int M, const void * coef_sos,const int16_t * coef_g,int16_t gain);
typedef void * tIirOldFxnInitFloat (void * objmem, int M, const void * coef_sos,int16_t gain);
typedef size_t tIirOldFxnGetScratch (int,int);
typedef void   tIirOldFxnProcess (void *, void * ,void * ,const void *, int );
typedef void   tIirOldFxnProcessFloat (void *, void * ,const void *, int );

typedef size_t tIirFxnAlloc(int);
typedef void * tIirFxnInit(void * objmem, int M, const void * coef_sos, const int16_t * coef_g, int16_t gain);
typedef void * tIirFxnInitFloat(void * objmem, int M, const void * coef_sos, int16_t gain);
typedef size_t tIirFxnGetScratch(int, int);
typedef void   tIirFxnProcess(void *, void *, void *, const void *, int);
typedef void   tIirFxnProcessFloat(void *, void *, const void *, int);
typedef size_t tIirFxnDelay(void *);

typedef struct
{
    tIirOldFxnAlloc      *alloc;
    tIirOldFxnInit       *init;
    tIirOldFxnGetScratch *getScratch;
    tIirOldFxnProcess    *process;
}
tIirOldDescr;

typedef struct
{
    tIirFxnAlloc      *alloc;
    tIirFxnInit       *init;
    tIirFxnGetScratch *getScratch;
    tIirFxnProcess    *process;
    tIirFxnDelay      *delay;
}
tIirDescr;


/* function reads coefficients from file and creates IIR structure. returns 0 if failed */
int te_create_iir_old(tTestEngContext * context);

/* function destroys IIR structure, returns 0 if failed */
int te_destroy_iir_old(tTestEngContext * context);

/* 
   Allocate vectors and load the data set for IIR:
*  vector X (in), vector Z (out) */
int te_loadFxn_iir_old(tTestEngContext * context);

/* Apply IIR function to the test case data set.
*  vector X (in), vector Z (out) */
void te_processFxn_iir_old(tTestEngContext * context);

/* Functions with delay imitation */

/* function reads coefficients from file and creates IIR structure. returns 0 if failed */
int te_create_iir(tTestEngContext * context);

/* function destroys IIR structure, returns 0 if failed */
int te_destroy_iir(tTestEngContext * context);

/*
Allocate vectors and load the data set for IIR:
*  vector X (in), vector Z (out) */
int te_loadFxn_iir(tTestEngContext * context);

/* Apply IIR function to the test case data set.
*  vector X (in), vector Z (out) */
void te_processFxn_iir(tTestEngContext * context);



typedef size_t tIirOldStereoFxnAlloc (int M);
typedef void * tIirOldStereoFxnInit(void * objmem, int M, const void * coef_sosl, const int16_t * coef_gl, int16_t gainl, const void * coef_sosr, const int16_t * coef_gr, int16_t gainr);
typedef void * tIirOldStereoFxnInitFloat(void * objmem, int M, const void * coef_sosl, int16_t gainl, const void * coef_sosr, int16_t gainr);
typedef size_t tIirOldStereoFxnGetScratch(int N, int M);
typedef void   tIirOldStereoFxnProcess(void * bqriir, void * s, void * r, const void * x, int N);
typedef void   tIirOldStereoFxnProcessFloat(void * bqriir, void * r, const void * x, int N);

typedef size_t tIirStereoFxnAlloc(int M);
typedef void * tIirStereoFxnInit(void * objmem, int M, const void * coef_sosl, const int16_t * coef_gl, int16_t gainl, const void * coef_sosr, const int16_t * coef_gr, int16_t gainr);
typedef void * tIirStereoFxnInitFloat(void * objmem, int M, const void * coef_sosl, int16_t gainl, const void * coef_sosr, int16_t gainr);
typedef size_t tIirStereoFxnGetScratch(int N, int M);
typedef void   tIirStereoFxnProcess(void * bqriir, void * s, void * r, const void * x, int N);
typedef void   tIirStereoFxnProcessFloat(void * bqriir, void * r, const void * x, int N);
typedef size_t tIirStereoFxnDelay(void *);

typedef struct
{
    tIirOldStereoFxnAlloc      *alloc;
    tIirOldStereoFxnInit       *init;
    tIirOldStereoFxnGetScratch *getScratch;
    tIirOldStereoFxnProcess    *process;
}
tIirOldStereoDescr;

typedef struct
{
    tIirStereoFxnAlloc      *alloc;
    tIirStereoFxnInit       *init;
    tIirStereoFxnGetScratch *getScratch;
    tIirStereoFxnProcess    *process;
    tIirStereoFxnDelay      *delay;
}
tIirStereoDescr;

/* function reads coefficients from file and creates IIR structure. returns 0 if failed */
int te_create_iir_stereo_old(tTestEngContext * context);

/* function destroys IIR structure, returns 0 if failed */
int te_destroy_iir_stereo_old(tTestEngContext * context);

/* 
   Allocate vectors and load the data set for IIR:
*  vector X (in), vector Z (out) */
int te_loadFxn_iir_stereo_old(tTestEngContext * context);

/* Apply IIR function to the test case data set.
*  vector X (in), vector Z (out) */
void te_processFxn_iir_stereo_old(tTestEngContext * context);

/* Functions with delay imitation */

/* function reads coefficients from file and creates IIR structure. returns 0 if failed */
int te_create_iir_stereo(tTestEngContext * context);

/* function destroys IIR structure, returns 0 if failed */
int te_destroy_iir_stereo(tTestEngContext * context);

/*
Allocate vectors and load the data set for IIR:
*  vector X (in), vector Z (out) */
int te_loadFxn_iir_stereo(tTestEngContext * context);

/* Apply IIR function to the test case data set.
*  vector X (in), vector Z (out) */
void te_processFxn_iir_stereo(tTestEngContext * context);



#endif

