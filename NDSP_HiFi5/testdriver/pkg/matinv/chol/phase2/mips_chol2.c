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
 * test module for testing cycle performance (matrix inversion and related functions)
 */
 
/* Cross-platform data type definitions. */
#include "NatureDSP_types.h"
/* Matrix Inversion functions API. */
#include "NatureDSP_Signal_matinv.h"
/* measurement utilities */
#include "mips.h"
#include "vectools.h"
#include "packages.h"

#define TEST_REAL       1
#define TEST_COMPLEX    1

#define TEST_GEN        1
#define TEST_FWD        1
#define TEST_BKW        1
#define TEST_MMSE       1
#define TEST_PREPROCESS 1
#define TEST_PINV       1

#define PROFILE_NORMALIZED_NEW(_cond, _verbose, _f_name, _f_arg, _file, _info, _fmt, _norm) \
{                                                                \
  if (_cond)                                                     \
    {                                                            \
      float32_t __norm;                                          \
      if (_verbose)                                              \
	        {                                                    \
           extern const char ANNOTATE_FUN_REF(_f_name)[];        \
           perf_info(_file, "%-12s\t", #_f_name);                \
           perf_info(_file,"%s\t", ANNOTATE_FUN_REF(_f_name));   \
           perf_info(_file, "%-16s\t",_info);                    \
	        } else                                               \
	        {                                                    \
           perf_info(_file, "%-12s\t%-16s\t", #_f_name,_info);   \
	        }                                                    \
      if (IS_PRESENT(_f_name))                                   \
	        {                                                    \
           PROFILE(_f_name _f_arg);                              \
           __norm=profiler_clk/(float32_t)(_norm);               \
           perf_info(_file, _fmt, profiler_clk,__norm);          \
           perf_info(_file, "\n");                               \
	        }                                                    \
	        else                                                 \
	        {                                                    \
           perf_info(_file, "NOT TESTED\n");                     \
	        }                                                    \
    }                                                            \
}

#define PROFILE_cmatcholdecompf(_cond,_verbose,M,N)                                                                                                     \
{                                                                                                                                                       \
	NASSERT(cmatcholdecomp##M##x##N##f_getScratchSize() <= sizeof(mips.scratch0)); 																		\
    PROFILE_NORMALIZED_NEW(_cond, _verbose, cmatcholdecomp##M##x##N##f, (mips.scratch0.cf32, mips.out0.cf32, mips.out1.cf32, mips.inp0.cf32, 1e-8f), fout, #M "x" #N, prf_cyclesmtx, 1);\
}
#define PROFILE_cmatcholfwdsubstf(_cond,_verbose,M,N,P)                                                                                                                          \
{                                                                                                                                                                              \
	NASSERT(cmatcholfwdsubst##M##x##N##f_getScratchSize() <= sizeof(mips.scratch0));\
    PROFILE_NORMALIZED_NEW(_cond, _verbose, cmatcholfwdsubst##M##x##N##f, (mips.scratch0.cf32, mips.out0.cf32, mips.inp0.cf32, mips.inp1.cf32, mips.inp0.cf32, mips.inp1.cf32), fout, #M "x" #N "x" #P, prf_cyclesmtx, 1); \
}
#define PROFILE_cmatcholbkwsubstf(_cond,_verbose,N,P)                                                                                                       \
{                                                                                                                                                         \
	NASSERT(cmatcholbkwsubst##N##x##N##f_getScratchSize() <= sizeof(mips.scratch0));\
    PROFILE_NORMALIZED_NEW(_cond, _verbose, cmatcholbkwsubst##N##x##N##f, (mips.scratch0.cf32, mips.out0.cf32, mips.inp0.cf32, mips.inp1.cf32, mips.inp0.cf32), fout, #N "x" #P, prf_cyclesmtx, 1); \
}
#define PROFILE_cmatcholmmsesolverf(_cond,_verbose,M,N,P)                                                                                                              \
{                                                                                                                                                                   \
	NASSERT(cmatcholmmsesolver##M##x##N##f_getScratchSize() <= sizeof(mips.scratch0));\
    PROFILE_NORMALIZED_NEW(_cond, _verbose, cmatcholmmsesolver##M##x##N##f, (mips.scratch0.cf32, mips.out0.cf32, mips.inp0.cf32, mips.inp1.cf32, 1e-8f), fout, #M "x" #N "x" #P, prf_cyclesmtx, 1); \
}
#define PROFILE_cmatcholpreprocessf(_cond,_verbose,M,N)                                                                                                    \
{                                                                                                                                                            \
	NASSERT(cmatcholpreprocess##M##x##N##f_getScratchSize() <= sizeof(mips.scratch0));																			 \
	PROFILE_NORMALIZED_NEW(_cond, _verbose, cmatcholpreprocess##M##x##N##f, (mips.scratch0.cf32, mips.out0.cf32, mips.inp0.cf32, 1e-8f), fout, #M "x" #N, prf_cyclesmtx, 1); \
}
#define PROFILE_cmatcholpseudoinvf(_cond,_verbose,M,N)                                                                                                      \
{                                                                                                                                                     \
	NASSERT(cmatcholpseudoinv##M##x##N##f_getScratchSize() <= sizeof(mips.scratch0));\
    PROFILE_NORMALIZED_NEW(_cond, _verbose, cmatcholpseudoinv##M##x##N##f, (mips.scratch0.cf32, mips.out0.cf32, mips.inp0.cf32, 1e-8f), fout, #M "x" #N, prf_cyclesmtx, 1); \
}

#define PROFILE_matcholdecompf(_cond,_verbose,M,N)                                                                                                   \
{                                                                                                                                                 \
	NASSERT(matcholdecomp##M##x##N##f_getScratchSize() <= sizeof(mips.scratch0));\
    PROFILE_NORMALIZED_NEW(_cond, _verbose, matcholdecomp##M##x##N##f, (mips.scratch0.f32, mips.out0.f32, mips.out1.f32, mips.inp0.f32, 1e-8f), fout, #M "x" #N, prf_cyclesmtx, 1); \
}
#define PROFILE_matcholfwdsubstf(_cond,_verbose,M,N,P)                                                                                                                     \
{                                                                                                                                                                        \
	NASSERT(matcholfwdsubst##M##x##N##f_getScratchSize() <= sizeof(mips.scratch0));\
    PROFILE_NORMALIZED_NEW(_cond, _verbose, matcholfwdsubst##M##x##N##f, (mips.scratch0.f32, mips.out0.f32, mips.inp0.f32, mips.inp1.f32, mips.inp0.f32, mips.inp1.f32), fout, #M "x" #N "x" #P, prf_cyclesmtx, 1); \
}
#define PROFILE_matcholbkwsubstf(_cond,_verbose,N,P)                                                                                                   \
{                                                                                                                                                    \
	NASSERT(matcholbkwsubst##N##x##N##f_getScratchSize() <= sizeof(mips.scratch0));\
    PROFILE_NORMALIZED_NEW(_cond, _verbose, matcholbkwsubst##N##x##N##f, (mips.scratch0.f32, mips.out0.f32, mips.inp0.f32, mips.inp1.f32, mips.inp0.f32), fout, #N "x" #P, prf_cyclesmtx, 1); \
}
#define PROFILE_matcholmmsesolverf(_cond,_verbose,M,N,P)                                                                                                           \
{                                                                                                                                                               \
	NASSERT(matcholmmsesolver##M##x##N##f_getScratchSize() <= sizeof(mips.scratch0));\
    PROFILE_NORMALIZED_NEW(_cond, _verbose, matcholmmsesolver##M##x##N##f, (mips.scratch0.f32, mips.out0.f32, mips.inp0.f32, mips.inp1.f32, 1e-8f), fout, #M "x" #N "x" #P, prf_cyclesmtx, 1); \
}
#define PROFILE_matcholpreprocessf(_cond,_verbose,M,N)                                                                                                      \
{                                                                                                                                                     \
	NASSERT(matcholpreprocess##M##x##N##f_getScratchSize() <= sizeof(mips.scratch0));\
    PROFILE_NORMALIZED_NEW(_cond, _verbose, matcholpreprocess##M##x##N##f, (mips.scratch0.f32, mips.out0.f32, mips.inp0.f32, 1e-8f), fout, #M "x" #N, prf_cyclesmtx, 1); \
}
#define PROFILE_matcholpseudoinvf(_cond,_verbose,M,N)                                                                                                      \
{                                                                                                                                                     \
	NASSERT(matcholpseudoinv##M##x##N##f_getScratchSize() <= sizeof(mips.scratch0));\
    PROFILE_NORMALIZED_NEW(_cond, _verbose, matcholpseudoinv##M##x##N##f, (mips.scratch0.f32, mips.out0.f32, mips.inp0.f32, 1e-8f), fout, #M "x" #N, prf_cyclesmtx, 1); \
}

static void mips_cchol2(int isFull, int isVerbose, FILE* fout )
{
#if TEST_COMPLEX
#if TEST_GEN
	PROFILE_cmatcholdecompf(isFull, isVerbose, 4, 4);
	PROFILE_cmatcholdecompf(isFull, isVerbose, 6, 6);
	PROFILE_cmatcholdecompf(isFull, isVerbose, 8, 8);
	PROFILE_cmatcholdecompf(1, isVerbose, 10, 10);
#endif

#if TEST_FWD
	PROFILE_cmatcholfwdsubstf(isFull, isVerbose, 4, 4, 1);
	PROFILE_cmatcholfwdsubstf(isFull, isVerbose, 6, 6, 1);
	PROFILE_cmatcholfwdsubstf(isFull, isVerbose, 8, 8, 1);
	PROFILE_cmatcholfwdsubstf(1, isVerbose, 10, 10, 1);
#endif

#if TEST_BKW
	PROFILE_cmatcholbkwsubstf(isFull, isVerbose, 4, 1);
	PROFILE_cmatcholbkwsubstf(isFull, isVerbose, 6, 1);
	PROFILE_cmatcholbkwsubstf(isFull, isVerbose, 8, 1);
	PROFILE_cmatcholbkwsubstf(1, isVerbose, 10, 1);
#endif

#if TEST_MMSE
	PROFILE_cmatcholmmsesolverf(isFull, isVerbose, 4, 4, 1);
	PROFILE_cmatcholmmsesolverf(isFull, isVerbose, 6, 6, 1);
	PROFILE_cmatcholmmsesolverf(isFull, isVerbose, 8, 8, 1);
	PROFILE_cmatcholmmsesolverf(1, isVerbose, 10, 10, 1);
#endif

#if TEST_PREPROCESS
	PROFILE_cmatcholpreprocessf(isFull, isVerbose, 4, 4);
	PROFILE_cmatcholpreprocessf(isFull, isVerbose, 6, 6);
	PROFILE_cmatcholpreprocessf(isFull, isVerbose, 8, 8);
	PROFILE_cmatcholpreprocessf(1, isVerbose, 10, 10);
#endif

#if TEST_PINV
	PROFILE_cmatcholpseudoinvf(isFull, isVerbose, 4, 4);
	PROFILE_cmatcholpseudoinvf(isFull, isVerbose, 6, 6);
	PROFILE_cmatcholpseudoinvf(isFull, isVerbose, 8, 8);
	PROFILE_cmatcholpseudoinvf(1, isVerbose, 10, 10);
#endif
#endif
}

void mips_rchol2(int isFull, int isVerbose, FILE* fout )
{
#if TEST_REAL
#if TEST_GEN
	PROFILE_matcholdecompf(isFull, isVerbose, 4, 4);
	PROFILE_matcholdecompf(isFull, isVerbose, 6, 6);
	PROFILE_matcholdecompf(isFull, isVerbose, 8, 8);
	PROFILE_matcholdecompf(1, isVerbose, 10, 10);
#endif

#if TEST_FWD
	PROFILE_matcholfwdsubstf(isFull, isVerbose, 4, 4, 1);
	PROFILE_matcholfwdsubstf(isFull, isVerbose, 6, 6, 1);
	PROFILE_matcholfwdsubstf(isFull, isVerbose, 8, 8, 1);
	PROFILE_matcholfwdsubstf(1, isVerbose, 10, 10, 1);
#endif

#if TEST_BKW
	PROFILE_matcholbkwsubstf(isFull, isVerbose, 4, 1);
	PROFILE_matcholbkwsubstf(isFull, isVerbose, 6, 1);
	PROFILE_matcholbkwsubstf(isFull, isVerbose, 8, 1);
	PROFILE_matcholbkwsubstf(1, isVerbose, 10, 1);
#endif

#if TEST_MMSE
	PROFILE_matcholmmsesolverf(isFull, isVerbose, 4, 4, 1);
	PROFILE_matcholmmsesolverf(isFull, isVerbose, 6, 6, 1);
	PROFILE_matcholmmsesolverf(isFull, isVerbose, 8, 8, 1);
	PROFILE_matcholmmsesolverf(1, isVerbose, 10, 10, 1);
#endif

#if TEST_PREPROCESS
	PROFILE_matcholpreprocessf(isFull, isVerbose, 4, 4);
	PROFILE_matcholpreprocessf(isFull, isVerbose, 6, 6);
	PROFILE_matcholpreprocessf(isFull, isVerbose, 8, 8);
	PROFILE_matcholpreprocessf(1, isVerbose, 10, 10);
#endif

#if TEST_PINV
	PROFILE_matcholpseudoinvf(isFull, isVerbose, 4, 4);
	PROFILE_matcholpseudoinvf(isFull, isVerbose, 6, 6);
	PROFILE_matcholpseudoinvf(isFull, isVerbose, 8, 8);
	PROFILE_matcholpseudoinvf(1, isVerbose, 10, 10);
#endif
#endif
}

void mips_chol2(int isFull, int isVerbose, FILE* fout)
{
    mips_cchol2(isFull, isVerbose, fout );
    mips_rchol2(isFull, isVerbose, fout );
}


