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
* Test module for testing cycle performance (Filtering)
*/
#include "config.h"
#include "packages.h"
#include LIBRARY_HEADER(fir)
#include "mips.h"

#define PROFILE_CONV2D(isFull,isVerbose,fun,M,N,P,Q,suffix1,suffix2)           \
{                                                                              \
    size_t inSz0=sizeof(mips.inp0.suffix2[0])*M*N;                             \
    size_t inSz1=sizeof(mips.inp1.suffix1[0])*P*Q;                             \
    size_t outSz0=sizeof(mips.out0.suffix1[0])*(M+P-1)*(N+Q-1);                \
    (void)inSz0,(void)inSz1,(void)outSz0;                                      \
    ASSERT(sizeof(mips.inp0)>=inSz0);                                          \
    ASSERT(sizeof(mips.inp1)>=inSz1);                                          \
    ASSERT(sizeof(mips.out0)>=outSz0);                                         \
    PROFILE_INVERTED(isFull,isVerbose,fun,                                     \
        ((void*)mips.scratch0.i8,mips.out0.suffix1,mips.inp0.suffix2,mips.inp1.suffix1,0,P,Q),\
        fout,"M=" #M ",N=" #N ",P=" #P ",Q=" #Q,prf_maccycle,(M*N*P*Q));       \
}

#define PROFILE_CONV2DF(isFull,isVerbose,fun,M,N,P,Q,suffix1,suffix2)         \
{                                                                             \
    size_t inSz0=sizeof(mips.inp0.suffix2[0])*M*N;                            \
    size_t inSz1=sizeof(mips.inp1.suffix1[0])*P*Q;                            \
    size_t outSz0=sizeof(mips.out0.suffix1[0])*(M+P-1)*(N+Q-1);               \
    (void)inSz0,(void)inSz1,(void)outSz0;                                     \
    ASSERT(sizeof(mips.inp0)>=inSz0);                                         \
    ASSERT(sizeof(mips.inp1)>=inSz1);                                         \
    ASSERT(sizeof(mips.out0)>=outSz0);                                        \
    PROFILE_INVERTED(isFull,isVerbose,fun,                                    \
        ((void*)mips.scratch0.i8,mips.out0.suffix1,mips.inp0.suffix2,mips.inp1.suffix1,P,Q), \
        fout,"M=" #M ",N=" #N ",P=" #P ",Q=" #Q,prf_maccycle,(M*N*P*Q));      \
}

void mips_conv2d2(int isFull, int isVerbose, FILE* fout)
{
    PROFILE_CONV2DF(isFull,isVerbose,conv2d_3x3f  ,3 ,3,256,256,f32,f32);
    PROFILE_CONV2DF(isFull,isVerbose,conv2d_5x5f  ,5 ,5,256,256,f32,f32);
    PROFILE_CONV2DF(isFull,isVerbose,conv2d_11x7f ,11,7,256,256,f32,f32);
    PROFILE_CONV2DF(1     ,isVerbose,conv2d_3x3f  ,3 ,3,128,256,f32,f32);
    PROFILE_CONV2DF(1     ,isVerbose,conv2d_5x5f  ,5 ,5,128,256,f32,f32);
    PROFILE_CONV2DF(1     ,isVerbose,conv2d_11x7f ,11,7,128,256,f32,f32);
    PROFILE_CONV2DF(isFull,isVerbose,conv2d_3x3f  ,3 ,3, 64, 64,f32,f32);
    PROFILE_CONV2DF(isFull,isVerbose,conv2d_5x5f  ,5 ,5, 64, 64,f32,f32);
    PROFILE_CONV2DF(isFull,isVerbose,conv2d_11x7f ,11,7, 64, 64,f32,f32);

    PROFILE_CONV2DF(isFull, isVerbose, conv2d_gen_3x3f, 3, 3, 256, 256, f32, f32);
    PROFILE_CONV2DF(isFull, isVerbose, conv2d_gen_5x5f, 5, 5, 256, 256, f32, f32);
    PROFILE_CONV2DF(isFull, isVerbose, conv2d_gen_11x7f, 11, 7, 256, 256, f32, f32);
    PROFILE_CONV2DF(1, isVerbose, conv2d_gen_3x3f, 3, 3, 128, 256, f32, f32);
    PROFILE_CONV2DF(1, isVerbose, conv2d_gen_5x5f, 5, 5, 128, 256, f32, f32);
    PROFILE_CONV2DF(1, isVerbose, conv2d_gen_11x7f, 11, 7, 128, 256, f32, f32);
    PROFILE_CONV2DF(isFull, isVerbose, conv2d_gen_3x3f, 3, 3, 64, 64, f32, f32);
    PROFILE_CONV2DF(isFull, isVerbose, conv2d_gen_5x5f, 5, 5, 64, 64, f32, f32);
    PROFILE_CONV2DF(isFull, isVerbose, conv2d_gen_11x7f, 11, 7, 64, 64, f32, f32);
}
