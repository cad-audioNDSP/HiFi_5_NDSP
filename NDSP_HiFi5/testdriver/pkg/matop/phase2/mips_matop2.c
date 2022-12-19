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
* Test module for testing cycle performance (Matrix Operations)
*/
#include "config.h"
#include "packages.h"
#include LIBRARY_HEADER(matop)
#include "mips.h"

#define PROFILE_mpyf(cond,verb,M,N,P)          PROFILE_INVERTED(cond,verb,mtx_mpyf,(pScr,mips.out0.f32, mips.inp0.f32, mips.inp1.f32, M, N, P),fout,#M"x"#N" x "#N"x"#P,prf_maccycle, (N*M*P));
#define PROFILE_mpyf_fast(cond,verb,M,N,P)     PROFILE_INVERTED(cond,verb,mtx_mpyf_fast,(pScr,mips.out0.f32, mips.inp0.f32, mips.inp1.f32, M, N, P),fout,#M"x"#N" x "#N"x"#P,prf_maccycle, (N*M*P));

#define PROFILE_mpytf(cond,verb,M,N,P)             PROFILE_INVERTED(cond,verb,mtx_mpytf,(pScr,mips.out0.f32, mips.inp0.f32, mips.inp1.f32, M, N, P),fout,#M"x"#N" x "#N"x"#P,prf_maccycle, (N*M*P));
#define PROFILE_mpytf_fast(cond,verb,M,N,P)        PROFILE_INVERTED(cond,verb,mtx_mpytf_fast,(pScr,mips.out0.f32, mips.inp0.f32, mips.inp1.f32, M, N, P),fout,#M"x"#N" x "#N"x"#P,prf_maccycle, (N*M*P));
#define PROFILE_mtx_vecmpyf(cond,verb,fun,M,N)     PROFILE_INVERTED(cond,verb,fun,(mips.out0.f32, mips.inp0.f32, mips.inp1.f32, M, N),fout,#M"x"#N" x "#N"x1",prf_maccycle, (N*M));

#define PROFILE_cmtx_mpytf_fast(isFull,isVerbose,fun,M,N,P)             \
{                                                                       \
    ASSERT(M*N*sizeof(mips.inp0.cf32[0])<=sizeof(mips.inp0));           \
    ASSERT(P*N*sizeof(mips.inp1.cf32[0])<=sizeof(mips.out1));           \
    ASSERT(M*P*sizeof(mips.out0.cf32[0])<=sizeof(mips.out1));           \
    ASSERT(fun##_getScratchSize(M,N,P)<=sizeof(mips.scratch0));         \
    PROFILE_INVERTED(isFull,isVerbose,fun,(&mips.scratch0,mips.out0.cf32,mips.inp1.cf32,mips.inp0.cf32,M,N,P),fout,#M"x"#N" x "#N"x"#P,prf_maccycle,(4*M*N*P)); \
}

#define PROFILE_cmtx_lrmpyf_fast(isFull,isVerbose,fun,N)                \
{                                                                       \
    ASSERT(  N*sizeof(mips.inp0.cf32[0])<=sizeof(mips.inp0));           \
    ASSERT(N*N*sizeof(mips.out1.cf32[0])<=sizeof(mips.out1));           \
    ASSERT(fun##_getScratchSize(N)<=sizeof(mips.scratch0));             \
    PROFILE_INVERTED(isFull,isVerbose,fun,(&mips.scratch0,mips.out1.cf32,mips.inp1.cf32,mips.inp0.cf32,N),fout,#N"x"#N" x "#N"x"#N" x "#N"x"#N,prf_maccycle,(2*4*N*N*N)); \
}

#define PROFILE_cmtx_vecmpytf_fast(isFull,isVerbose,fun,N)              \
{                                                                       \
    ASSERT(  N*sizeof(mips.inp0.cf32[0])<=sizeof(mips.inp0));           \
    ASSERT(N*N*sizeof(mips.out1.cf32[0])<=sizeof(mips.out1));           \
    PROFILE_INVERTED(isFull,isVerbose,fun,(mips.out1.cf32,mips.inp0.cf32,N),fout,#N"x1 x 1x"#N,prf_maccycle,(4*N*N)); \
}

#define PROFILE_TRANSP(isFull,isVerbose,fun,coef,M,N,suffix)            \
{                                                                       \
    ASSERT(M*N*sizeof(mips.inp0.suffix[0])<=sizeof(mips.inp0));         \
    ASSERT(M*N*sizeof(mips.out0.suffix[0])<=sizeof(mips.out0));         \
    PROFILE_INVERTED(isFull,isVerbose,fun,(mips.out0.suffix,mips.inp0.suffix,M,N),fout,"M=" #M ",N=" #N ,prf_ptscycle2,(1*M*N)); \
}

void mips_matop2(int isFull, int isVerbose, FILE * fout)
{
  void* pScr = (void*)mips.scratch0.i32;

  PROFILE_mpyf(     1, isVerbose, 16,16,16);
  PROFILE_mpyf(     1, isVerbose, 32,32,32);
  PROFILE_mpyf(isFull, isVerbose, 40,80,8);
  PROFILE_mpyf(isFull, isVerbose, 40,81,8);
  PROFILE_mpyf(isFull, isVerbose, 40,82,8);
  PROFILE_mpyf(isFull, isVerbose, 40,83,8);
  PROFILE_mpyf(isFull, isVerbose, 2,100,8);
  PROFILE_mpyf(isFull, isVerbose, 8,80,2);
  PROFILE_mpyf(isFull, isVerbose, 8,4,2);
  PROFILE_mpyf(isFull, isVerbose, 8,16,2);
  PROFILE_mpyf(isFull, isVerbose, 8,32,2);
  PROFILE_mpyf_fast(     1, isVerbose, 16,16,16);
  PROFILE_mpyf_fast(     1, isVerbose, 32,32,32);
  PROFILE_mpyf_fast(isFull, isVerbose,8,80,4);
  PROFILE_mpyf_fast(isFull, isVerbose,8,84,4);
  PROFILE_mpyf_fast(isFull, isVerbose,8,4 ,4);
  PROFILE_mpyf_fast(     1, isVerbose,8,16,4);
  PROFILE_mpyf_fast(isFull, isVerbose,8,32,4);

  PROFILE_mpytf(     1, isVerbose, 16,16,16);
  PROFILE_mpytf(     1, isVerbose, 32,32,32);
  PROFILE_mpytf(isFull, isVerbose, 40,80,8);
  PROFILE_mpytf(isFull, isVerbose, 40,81,8);
  PROFILE_mpytf(isFull, isVerbose, 40,82,8);
  PROFILE_mpytf(isFull, isVerbose, 40,83,8);
  PROFILE_mpytf(isFull, isVerbose, 2,100,8);
  PROFILE_mpytf(isFull, isVerbose, 8,80,2);
  PROFILE_mpytf(isFull, isVerbose, 8,4,2);
  PROFILE_mpytf(isFull, isVerbose, 8,16,2);
  PROFILE_mpytf(isFull, isVerbose, 8,32,2);
  PROFILE_mpytf_fast(     1, isVerbose, 16,16,16);
  PROFILE_mpytf_fast(     1, isVerbose, 32,32,32);
  PROFILE_mpytf_fast(isFull, isVerbose,8,80,4);
  PROFILE_mpytf_fast(isFull, isVerbose,8,84,4);
  PROFILE_mpytf_fast(isFull, isVerbose,8,4 ,4);
  PROFILE_mpytf_fast(     1, isVerbose,8,16,4);
  PROFILE_mpytf_fast(isFull, isVerbose,8,32,4);

  PROFILE_mtx_vecmpyf(isFull, isVerbose, mtx_vecmpyf,     16,100);
  PROFILE_mtx_vecmpyf(isFull, isVerbose, mtx_vecmpyf,     16,101);
  PROFILE_mtx_vecmpyf(isFull, isVerbose, mtx_vecmpyf,     16,102);
  PROFILE_mtx_vecmpyf(isFull, isVerbose, mtx_vecmpyf,     16,103);
  PROFILE_mtx_vecmpyf(isFull, isVerbose, mtx_vecmpyf,     16,104);
  PROFILE_mtx_vecmpyf(isFull, isVerbose, mtx_vecmpyf,     40,40);
  PROFILE_mtx_vecmpyf(isFull, isVerbose, mtx_vecmpyf_fast,16,100);
  PROFILE_mtx_vecmpyf(     1, isVerbose, mtx_vecmpyf_fast,16,104);
  PROFILE_mtx_vecmpyf(isFull, isVerbose, mtx_vecmpyf_fast,40,40);

  PROFILE_cmtx_mpytf_fast(isFull, isVerbose, cmtx_mpytf_fast,  64,  64,  64);
  PROFILE_cmtx_mpytf_fast(     1, isVerbose, cmtx_mpytf_fast,  96,  96,  96);
  PROFILE_cmtx_mpytf_fast(isFull, isVerbose, cmtx_mpytf_fast, 128, 128, 128);
  PROFILE_cmtx_mpytf_fast(isFull, isVerbose, cmtx_mpytf_fast, 160, 160, 160);

  PROFILE_cmtx_lrmpyf_fast(isFull, isVerbose, cmtx_lrmpyf_fast,  64);
  PROFILE_cmtx_lrmpyf_fast(     1, isVerbose, cmtx_lrmpyf_fast,  96);
  PROFILE_cmtx_lrmpyf_fast(isFull, isVerbose, cmtx_lrmpyf_fast, 128);
  PROFILE_cmtx_lrmpyf_fast(isFull, isVerbose, cmtx_lrmpyf_fast, 160);

  PROFILE_cmtx_vecmpytf_fast(isFull, isVerbose, cmtx_vecmpytf_fast,  64);
  PROFILE_cmtx_vecmpytf_fast(     1, isVerbose, cmtx_vecmpytf_fast,  96);
  PROFILE_cmtx_vecmpytf_fast(isFull, isVerbose, cmtx_vecmpytf_fast, 128);
  PROFILE_cmtx_vecmpytf_fast(isFull, isVerbose, cmtx_vecmpytf_fast, 160);

  PROFILE_TRANSP(isFull,isVerbose,mtx_transposef ,sizeof(float32_t),16,16,f32);
  PROFILE_TRANSP(isFull,isVerbose,mtx_transposef ,sizeof(float32_t),27,27,f32);
  PROFILE_TRANSP(1     ,isVerbose,mtx_transposef ,sizeof(float32_t),32,32,f32);
  PROFILE_TRANSP(isFull,isVerbose,mtx_transposef ,sizeof(float32_t),39,39,f32);
  PROFILE_TRANSP(isFull,isVerbose,mtx_transposef ,sizeof(float32_t),48,48,f32);

  PROFILE_TRANSP(isFull,isVerbose,mtx_transposef_fast ,sizeof(float32_t), 8, 8,f32);
  PROFILE_TRANSP(isFull,isVerbose,mtx_transposef_fast ,sizeof(float32_t),16,16,f32);
  PROFILE_TRANSP(1     ,isVerbose,mtx_transposef_fast ,sizeof(float32_t),32,32,f32);
  PROFILE_TRANSP(isFull,isVerbose,mtx_transposef_fast ,sizeof(float32_t),48,48,f32);
}
