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
 * Test-engine add-on for FIR categories 
 */
#ifndef TESTENG__FIR_H__
#define TESTENG__FIR_H__

/* Cross-platform data type definitions. */
#include "types.h"
/* Test engine API. */
#include "testeng.h"

/* FIR filters types (bits 1:0 from extraParam): */
#define TE_FIR_FIR 0x0
#define TE_FIR_DN  0x1 /* downsampler */
#define TE_FIR_UP  0x2 /* upsampler */
#define TE_FIR_PP  0x3 /* polyphase filter */
#define TE_FIR_FILTER_TYPE_MASK 0x3

#define TE_FIR_FILTER_32X16 0x4 /* 32x16 filters */

/* FIR filters API (bits 5:4 from extraParam): */
#define TE_FIR_OLDGENERIC   0x00 /* legacy API (generic filters)*/
#define TE_FIR_OLDEXTIR     0x10 /* legacy API (filters with external IR) */
#define TE_FIR_OLDDECIMATOR 0x20 /* legacy API (decimator/interpolator) */
#define TE_FIR_POLYPHASE    0x30 /* Polyphase filter API */
#define TE_FIR_FILTER_API_MASK 0x70

/* set extIR==1 for filters with external IR */
#define TE_FIR_EXTIR        0x80

/* Fixed-point FIR with Left Shift amount specified for each input block. */
#define TE_FIR_LSH          0x100

/*
 * parameters in extraParam (bits 6:0) for correlation/convolution descriptions
 */
#define TE_FIR_AUTOCORR_API  0x0
#define TE_FIR_CROSSCORR_API 0x1
#define TE_FIR_CONVOLVE_API  0x2
#define TE_FIR_LINEAR_API    0x4 /* linear correlation/convolution */
#define TE_FIR_OTHER_API_MASK 0x7

#define TE_FIR_OTHER_32X16    0x08 /* 32x16 functions */
#define TE_FIR_OTHER_24X24    0x10 /* 24x24 functions */
#define TE_FIR_OTHER_EP       0x20 /* 32x32 extended precision functions */
#define TE_FIR_OTHER_TYPE_MASK 0x38

/* function reads IR from file and creates FIR structure. returns 0 if failed */
int te_create_fir(tTestEngContext * context);

/* function destroys FIR structure, returns 0 if failed */
int te_destroy_fir(tTestEngContext * context);

/* Allocate vectors and load the data set for FIR:
 * vector X (in), vector Z (out) */
int te_loadFxn_fir(tTestEngContext * context);
/* Apply FIR function to the test case data set.
 * vector X (in), vector Z (out) */
void te_processFxn_fir(tTestEngContext * context);

/* Allocate vectors and load the data set for autocorrelation functions:
 * vector X (in), vector Z (out) */
int te_loadFxn_autocorr( tTestEngContext * context );
/* Allocate vectors and load the data set for crosscorrelation functions:
 * vector X (in), vector Z (out) */
int te_loadFxn_crosscorr( tTestEngContext * context );
/* Apply the target function to the test case data set (crosscorrelation):
 * vector X (in), vector Z (out) */
void te_processFxn_autocorr( tTestEngContext * context );
/* Apply the target function to the test case data set (crosscorrelation):
 * vector X (in), vector Z (out) */
void te_processFxn_crosscorr( tTestEngContext * context );

typedef struct {
    size_t (*getScratchSize)(int N);
} te_gccphat_api_t;

/* Allocate vectors and load the data set for GCC-PHAT functions:
 * vector X (in), vector Y (in), vector U (in), vector V (temp), vector Z (out) */
int te_loadFxn_gccphat( tTestEngContext * context );
/* Apply the target function to the test case data set (GCC-PHAT):
  * vector X (in), vector Y (in), vector U (in), vector V (temp), vector Z (out) */
void te_processFxn_gccphat( tTestEngContext * context );

typedef struct
{
    tVec h;          /* impulse response */
    tVec normmu;     /* delay line */
}
tLmsContext;

/* function reads impulse response from file and creates FIR structure. returns 0 if failed */
int te_create_lms(tTestEngContext * context);
void lms_vecsFree(tTestEngContext * context);
/* function destroys FIR structure, returns 0 if failed */
int te_destroy_lms(tTestEngContext * context);
/* Allocate vectors and load the data set for lms functions: */
int te_loadFxn_lms( tTestEngContext * context );
/* Apply the target function to the test case data set (LMS):
 * vector X (in), vector Z (out) */
void te_processFxn_lms( tTestEngContext * context );

/* LMS convergence tess: */
int te_loadFxn_lmsconv( tTestEngContext * context );
void te_processFxn_lmsconv( tTestEngContext * context );
#endif

