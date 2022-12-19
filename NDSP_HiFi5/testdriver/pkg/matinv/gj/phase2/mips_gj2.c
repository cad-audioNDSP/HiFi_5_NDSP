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
* Test module for testing cycle performance (Matrix Decomposition 
* and Inversion)
*/

/* Cross-platform data type definitions. */
#include "types.h"
/* Test environment configuration. */
#include "config.h"
#include "packages.h"
/* DSP Library API. */
#include LIBRARY_HEADER(matinv)
/* MIPS measurement means. */
#include "mips.h"
/* Utility functions and macros. */
#include "utils.h"

#define PROFILE_MATINV(cond,verb,fun,N,suffix)          PROFILE_NORMALIZED(cond,verb,fun,(mips.scratch0.suffix,mips.inp0.suffix   ),fout,"",prf_cyclesmtx,1);
#define PROFILE_GJELIM(cond,verb,fun,N,suffix)          PROFILE_NORMALIZED(cond,verb,fun,(mips.scratch0.suffix,mips.out0.suffix,mips.inp0.suffix,mips.inp1.suffix),fout,"",prf_cyclesmtx,1);

void mips_gj2(int isFull, int isVerbose, FILE * fout)
{
    int i;
    /* fill with random floating point data */
    Rand_reset(12,273);
    for (i=0; i<200; i++)
    {
        mips.inp0.f32[i]=((int16_t)Rand())*(1.f/32768.f);
    }
    PROFILE_MATINV(     1, isVerbose, cmtx_inv2x2f  , 2,cf32);
    PROFILE_MATINV(isFull, isVerbose, cmtx_inv3x3f  , 3,cf32);
    PROFILE_MATINV(     1, isVerbose, cmtx_inv4x4f  , 4,cf32);
    PROFILE_MATINV(isFull, isVerbose, cmtx_inv6x6f  , 6,cf32);
    PROFILE_MATINV(     1, isVerbose, cmtx_inv8x8f  , 8,cf32);
    PROFILE_MATINV(isFull, isVerbose, cmtx_inv10x10f,10,cf32);

    PROFILE_MATINV(     1, isVerbose, mtx_inv2x2f  , 2,f32);
    PROFILE_MATINV(isFull, isVerbose, mtx_inv3x3f  , 3,f32);
    PROFILE_MATINV(     1, isVerbose, mtx_inv4x4f  , 4,f32);
    PROFILE_MATINV(isFull, isVerbose, mtx_inv6x6f  , 6,f32);
    PROFILE_MATINV(     1, isVerbose, mtx_inv8x8f  , 8,f32);
    PROFILE_MATINV(isFull, isVerbose, mtx_inv10x10f,10,f32);

    PROFILE_GJELIM(     1, isVerbose, cmtx_gjelim2x2f  , 2,cf32);
    PROFILE_GJELIM(isFull, isVerbose, cmtx_gjelim3x3f  , 3,cf32);
    PROFILE_GJELIM(     1, isVerbose, cmtx_gjelim4x4f  , 4,cf32);
    PROFILE_GJELIM(isFull, isVerbose, cmtx_gjelim6x6f  , 6,cf32);
    PROFILE_GJELIM(     1, isVerbose, cmtx_gjelim8x8f  , 8,cf32);
    PROFILE_GJELIM(isFull, isVerbose, cmtx_gjelim10x10f,10,cf32);
    PROFILE_GJELIM(     1, isVerbose, mtx_gjelim2x2f   , 2,f32);
    PROFILE_GJELIM(isFull, isVerbose, mtx_gjelim3x3f   , 3,f32);
    PROFILE_GJELIM(     1, isVerbose, mtx_gjelim4x4f   , 4,f32);
    PROFILE_GJELIM(isFull, isVerbose, mtx_gjelim6x6f   , 6,f32);
    PROFILE_GJELIM(     1, isVerbose, mtx_gjelim8x8f   , 8,f32);
    PROFILE_GJELIM(isFull, isVerbose, mtx_gjelim10x10f ,10,f32);
}
