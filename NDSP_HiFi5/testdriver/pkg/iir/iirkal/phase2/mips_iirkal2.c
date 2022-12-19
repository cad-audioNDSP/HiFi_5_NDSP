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
* Test module for testing cycle performance (Kalman Filter functions)
*/

/* Cross-platform data type definitions. */
#include "types.h"
/* Test environment configuration. */
#include "config.h"
#include "packages.h"
/* MIPS benchmarks instrumentation. */
#include "mips.h"
/* DSP Library API. */
#include LIBRARY_HEADER(iir)

#define REAL_f      f32

#define PROFILE_kalmanupd1f(cond,verb,N) {              \
    size_t szScr = kalmanupd1f_getScratchSize(N);       \
    (void)szScr; NASSERT(szScr<=sizeof(mips.scratch0)); \
    PROFILE_NORMALIZED(cond,verb,kalmanupd1f,(mips.scratch0.i8,mips.out0.f32,mips.inp1.f32,mips.inp0.f32,mips.inp2.f32,N),fout,"N="#N,prf_cyclespts,(N)); }

void mips_iirkal2(int isFull, int isVerbose, FILE * fout)
{
    PROFILE_kalmanupd1f(isFull, isVerbose, 32);
    PROFILE_kalmanupd1f(isFull, isVerbose, 96);
    PROFILE_kalmanupd1f(     1, isVerbose, 160);
    PROFILE_kalmanupd1f(isFull, isVerbose, 320);
}
