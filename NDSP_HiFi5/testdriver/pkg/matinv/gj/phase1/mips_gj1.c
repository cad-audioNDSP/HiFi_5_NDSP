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

#define PROFILE_MATINV_32X32(cond,verb,fun,N,suffix)    PROFILE_NORMALIZED(cond,verb,fun,(mips.scratch0.suffix,mips.inp0.suffix,31),fout,"",prf_cyclesmtx,1);
#define PROFILE_GJELIM_32X32(cond,verb,fun,N,suffix)    PROFILE_NORMALIZED(cond,verb,fun,(mips.scratch0.suffix,mips.out0.suffix,mips.inp0.suffix,mips.inp1.suffix,31,31),fout,"",prf_cyclesmtx,1);

void mips_gj1(int isFull, int isVerbose, FILE * fout)
{
    int i;
    /* fill with random floating point data */
    Rand_reset(12,273);
    for (i=0; i<200; i++)
    {
        mips.inp0.i32[i]=(int32_t)Rand();
        mips.inp1.i32[i]=(int32_t)Rand();
    }
    PROFILE_MATINV_32X32(     1, isVerbose, cmtx_inv2x2_32x32  , 2,ci32);
    PROFILE_MATINV_32X32(isFull, isVerbose, cmtx_inv3x3_32x32  , 3,ci32);
    PROFILE_MATINV_32X32(     1, isVerbose, cmtx_inv4x4_32x32  , 4,ci32);
    PROFILE_MATINV_32X32(isFull, isVerbose, cmtx_inv6x6_32x32  , 6,ci32);
    PROFILE_MATINV_32X32(     1, isVerbose, cmtx_inv8x8_32x32  , 8,ci32);
    PROFILE_MATINV_32X32(isFull, isVerbose, cmtx_inv10x10_32x32,10,ci32);
    PROFILE_MATINV_32X32(     1, isVerbose, mtx_inv2x2_32x32   , 2,i32);
    PROFILE_MATINV_32X32(isFull, isVerbose, mtx_inv3x3_32x32   , 3,i32);
    PROFILE_MATINV_32X32(     1, isVerbose, mtx_inv4x4_32x32   , 4,i32);
    PROFILE_MATINV_32X32(isFull, isVerbose, mtx_inv6x6_32x32   , 6,i32);
    PROFILE_MATINV_32X32(     1, isVerbose, mtx_inv8x8_32x32   , 8,i32);
    PROFILE_MATINV_32X32(isFull, isVerbose, mtx_inv10x10_32x32 ,10,i32);

    PROFILE_GJELIM_32X32(     1, isVerbose, cmtx_gjelim2x2_32x32  , 2,ci32);
    PROFILE_GJELIM_32X32(isFull, isVerbose, cmtx_gjelim3x3_32x32  , 3,ci32);
    PROFILE_GJELIM_32X32(     1, isVerbose, cmtx_gjelim4x4_32x32  , 4,ci32);
    PROFILE_GJELIM_32X32(isFull, isVerbose, cmtx_gjelim6x6_32x32  , 6,ci32);
    PROFILE_GJELIM_32X32(     1, isVerbose, cmtx_gjelim8x8_32x32  , 8,ci32);
    PROFILE_GJELIM_32X32(isFull, isVerbose, cmtx_gjelim10x10_32x32,10,ci32);
    PROFILE_GJELIM_32X32(     1, isVerbose, mtx_gjelim2x2_32x32   , 2,i32);
    PROFILE_GJELIM_32X32(isFull, isVerbose, mtx_gjelim3x3_32x32   , 3,i32);
    PROFILE_GJELIM_32X32(     1, isVerbose, mtx_gjelim4x4_32x32   , 4,i32);
    PROFILE_GJELIM_32X32(isFull, isVerbose, mtx_gjelim6x6_32x32   , 6,i32);
    PROFILE_GJELIM_32X32(     1, isVerbose, mtx_gjelim8x8_32x32   , 8,i32);
    PROFILE_GJELIM_32X32(isFull, isVerbose, mtx_gjelim10x10_32x32 ,10,i32);
}
