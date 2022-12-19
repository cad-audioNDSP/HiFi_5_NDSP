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

#define PROFILE_cmatcholdecomp_32x32(_cond,_verbose,M,N)                          \
{                                                                              \
    NASSERT(cmatcholdecomp##M##x##N##_32x32_getScratchSize() <= sizeof(mips.scratch0));\
    PROFILE_NORMALIZED_NEW(_cond, _verbose, cmatcholdecomp##M##x##N##_32x32, (mips.scratch0.u8, mips.out0.ci32, mips.out1.ci32, mips.inp0.ci32, 0,0), fout, #M "x" #N, prf_cyclesmtx, 1);\
}
#define PROFILE_cmatcholfwdsubst_32x32(_cond,_verbose,M,N,P)                        \
{                                                                                 \
    NASSERT(cmatcholfwdsubst##M##x##N##_32x32_getScratchSize() <= sizeof(mips.scratch0));\
    PROFILE_NORMALIZED_NEW(_cond, _verbose, cmatcholfwdsubst##M##x##N##_32x32, (mips.scratch0.u8, mips.out0.ci32, mips.inp0.ci32, mips.inp1.ci32, mips.inp0.ci32, mips.inp1.ci32,0), fout, #M "x" #N "x" #P, prf_cyclesmtx, 1); \
}
#define PROFILE_cmatcholbkwsubst_32x32(_cond,_verbose,N,P)                          \
{                                                                                 \
    NASSERT(cmatcholbkwsubst##N##x##N##_32x32_getScratchSize() <= sizeof(mips.scratch0));\
    PROFILE_NORMALIZED_NEW(_cond, _verbose, cmatcholbkwsubst##N##x##N##_32x32, (mips.scratch0.u8, mips.out0.ci32, mips.inp0.ci32, mips.inp1.ci32, mips.inp0.ci32,0), fout, #N "x" #P, prf_cyclesmtx, 1); \
}
#define PROFILE_cmatcholmmsesolver_32x32(_cond,_verbose,M,N,P)                        \
{                                                                                  \
    NASSERT(cmatcholmmsesolver##M##x##N##_32x32_getScratchSize() <= sizeof(mips.scratch0));\
    PROFILE_NORMALIZED_NEW(_cond, _verbose, cmatcholmmsesolver##M##x##N##_32x32, (mips.scratch0.ci32, mips.out0.ci32, mips.inp0.ci32, mips.inp1.ci32, 0,0,0,0), fout, #M "x" #N "x" #P, prf_cyclesmtx, 1); \
}
#define PROFILE_cmatcholpreprocess_32x32(_cond,_verbose,M,N)                            \
{                                                                                          \
    NASSERT(cmatcholpreprocess##M##x##N##_32x32_getScratchSize() <= sizeof(mips.scratch0));  \
    PROFILE_NORMALIZED_NEW(_cond, _verbose, cmatcholpreprocess##M##x##N##_32x32, (mips.scratch0.u8, mips.out0.ci64, mips.inp0.ci32, 0,0), fout, #M "x" #N, prf_cyclesmtx, 1); \
}
#define PROFILE_cmatcholpseudoinv_32x32(_cond,_verbose,M,N)                          \
{                                                                                  \
    NASSERT(cmatcholpseudoinv##M##x##N##_32x32_getScratchSize() <= sizeof(mips.scratch0));\
    PROFILE_NORMALIZED_NEW(_cond, _verbose, cmatcholpseudoinv##M##x##N##_32x32, (mips.scratch0.u8, mips.out0.ci32, mips.inp0.ci32, 0,0,0,0), fout, #M "x" #N, prf_cyclesmtx, 1); \
}
#define PROFILE_matcholdecomp_32x32(_cond,_verbose,M,N)                           \
{                                                                              \
    NASSERT(matcholdecomp##M##x##N##_32x32_getScratchSize() <= sizeof(mips.scratch0)); \
    PROFILE_NORMALIZED_NEW(_cond, _verbose, matcholdecomp##M##x##N##_32x32, (mips.scratch0.u8, mips.out0.i32, mips.out1.ci32, mips.inp0.i32, 0,0), fout, #M "x" #N, prf_cyclesmtx, 1);\
}
#define PROFILE_matcholfwdsubst_32x32(_cond,_verbose,M,N,P)                        \
{                                                                                \
    NASSERT(matcholfwdsubst##M##x##N##_32x32_getScratchSize() <= sizeof(mips.scratch0));\
    PROFILE_NORMALIZED_NEW(_cond, _verbose, matcholfwdsubst##M##x##N##_32x32, (mips.scratch0.u8, mips.out0.i32, mips.inp0.i32, mips.inp1.ci32, mips.inp0.i32, mips.inp1.i32,0), fout, #M "x" #N "x" #P, prf_cyclesmtx, 1); \
}
#define PROFILE_matcholbkwsubst_32x32(_cond,_verbose,N,P)                          \
{                                                                                \
    NASSERT(matcholbkwsubst##N##x##N##_32x32_getScratchSize() <= sizeof(mips.scratch0));\
    PROFILE_NORMALIZED_NEW(_cond, _verbose, matcholbkwsubst##N##x##N##_32x32, (mips.scratch0.u8, mips.out0.i32, mips.inp0.i32, mips.inp1.ci32, mips.inp0.i32,0), fout, #N "x" #P, prf_cyclesmtx, 1); \
}
#define PROFILE_matcholmmsesolver_32x32(_cond,_verbose,M,N,P)                        \
{                                                                                 \
    NASSERT(matcholmmsesolver##M##x##N##_32x32_getScratchSize() <= sizeof(mips.scratch0));\
    PROFILE_NORMALIZED_NEW(_cond, _verbose, matcholmmsesolver##M##x##N##_32x32, (mips.scratch0.u8, mips.out0.i32, mips.inp0.i32, mips.inp1.i32, 0,0,0,0), fout, #M "x" #N "x" #P, prf_cyclesmtx, 1); \
}
#define PROFILE_matcholpreprocess_32x32(_cond,_verbose,M,N)                            \
{                                                                                         \
    NASSERT(matcholpreprocess##M##x##N##_32x32_getScratchSize() <= sizeof(mips.scratch0));  \
    PROFILE_NORMALIZED_NEW(_cond, _verbose, matcholpreprocess##M##x##N##_32x32, (mips.scratch0.u8, mips.out0.i64, mips.inp0.i32, 0,0), fout, #M "x" #N, prf_cyclesmtx, 1); \
}
#define PROFILE_matcholpseudoinv_32x32(_cond,_verbose,M,N)                          \
{                                                                                 \
    NASSERT(matcholpseudoinv##M##x##N##_32x32_getScratchSize() <= sizeof(mips.scratch0));\
    PROFILE_NORMALIZED_NEW(_cond, _verbose, matcholpseudoinv##M##x##N##_32x32, (mips.scratch0.u8, mips.out0.i32, mips.inp0.i32, 0,0,0,0), fout, #M "x" #N, prf_cyclesmtx, 1); \
}

static void mips_cchol2(int isFull, int isVerbose, FILE* fout )
{
    PROFILE_cmatcholdecomp_32x32(isFull, isVerbose,  4,  4);
    PROFILE_cmatcholdecomp_32x32(isFull, isVerbose,  6,  6);
    PROFILE_cmatcholdecomp_32x32(isFull, isVerbose,  8,  8);
    PROFILE_cmatcholdecomp_32x32(1     , isVerbose, 10, 10);

    PROFILE_cmatcholfwdsubst_32x32(isFull, isVerbose, 4,  4, 1);
    PROFILE_cmatcholfwdsubst_32x32(isFull, isVerbose, 6,  6, 1);
    PROFILE_cmatcholfwdsubst_32x32(isFull, isVerbose, 8,  8, 1);
    PROFILE_cmatcholfwdsubst_32x32(1     , isVerbose,10, 10, 1);

    PROFILE_cmatcholbkwsubst_32x32(isFull, isVerbose,  4, 1);
    PROFILE_cmatcholbkwsubst_32x32(isFull, isVerbose,  6, 1);
    PROFILE_cmatcholbkwsubst_32x32(isFull, isVerbose,  8, 1);
    PROFILE_cmatcholbkwsubst_32x32(1     , isVerbose, 10, 1);

    PROFILE_cmatcholmmsesolver_32x32(isFull, isVerbose,  4,  4, 1);
    PROFILE_cmatcholmmsesolver_32x32(isFull, isVerbose,  6,  6, 1);
    PROFILE_cmatcholmmsesolver_32x32(isFull, isVerbose,  8,  8, 1);
    PROFILE_cmatcholmmsesolver_32x32(1     , isVerbose, 10, 10, 1);

    PROFILE_cmatcholpreprocess_32x32(isFull, isVerbose,  4,  4);
    PROFILE_cmatcholpreprocess_32x32(isFull, isVerbose,  6,  6);
    PROFILE_cmatcholpreprocess_32x32(isFull, isVerbose,  8,  8);
    PROFILE_cmatcholpreprocess_32x32(1     , isVerbose, 10, 10);

    PROFILE_cmatcholpseudoinv_32x32(isFull, isVerbose,  4,  4);
    PROFILE_cmatcholpseudoinv_32x32(isFull, isVerbose,  6,  6);
    PROFILE_cmatcholpseudoinv_32x32(isFull, isVerbose,  8,  8);
    PROFILE_cmatcholpseudoinv_32x32(1     , isVerbose, 10, 10);

    PROFILE_matcholdecomp_32x32(isFull, isVerbose,  4,  4);
    PROFILE_matcholdecomp_32x32(isFull, isVerbose,  6,  6);
    PROFILE_matcholdecomp_32x32(isFull, isVerbose,  8,  8);
    PROFILE_matcholdecomp_32x32(1     , isVerbose, 10, 10);

    PROFILE_matcholfwdsubst_32x32(isFull, isVerbose, 4,  4, 1);
    PROFILE_matcholfwdsubst_32x32(isFull, isVerbose, 6,  6, 1);
    PROFILE_matcholfwdsubst_32x32(isFull, isVerbose, 8,  8, 1);
    PROFILE_matcholfwdsubst_32x32(1     , isVerbose,10, 10, 1);

    PROFILE_matcholbkwsubst_32x32(isFull, isVerbose,  4, 1);
    PROFILE_matcholbkwsubst_32x32(isFull, isVerbose,  6, 1);
    PROFILE_matcholbkwsubst_32x32(isFull, isVerbose,  8, 1);
    PROFILE_matcholbkwsubst_32x32(1     , isVerbose, 10, 1);

    PROFILE_matcholmmsesolver_32x32(isFull, isVerbose,  4,  4, 1);
    PROFILE_matcholmmsesolver_32x32(isFull, isVerbose,  6,  6, 1);
    PROFILE_matcholmmsesolver_32x32(isFull, isVerbose,  8,  8, 1);
    PROFILE_matcholmmsesolver_32x32(1     , isVerbose, 10, 10, 1);

    PROFILE_matcholpreprocess_32x32(isFull, isVerbose,  4,  4);
    PROFILE_matcholpreprocess_32x32(isFull, isVerbose,  6,  6);
    PROFILE_matcholpreprocess_32x32(isFull, isVerbose,  8,  8);
    PROFILE_matcholpreprocess_32x32(1     , isVerbose, 10, 10);

    PROFILE_matcholpseudoinv_32x32(isFull, isVerbose,  4,  4);
    PROFILE_matcholpseudoinv_32x32(isFull, isVerbose,  6,  6);
    PROFILE_matcholpseudoinv_32x32(isFull, isVerbose,  8,  8);
    PROFILE_matcholpseudoinv_32x32(1     , isVerbose, 10, 10);
}

void mips_chol1(int isFull, int isVerbose, FILE* fout)
{
    mips_cchol2(isFull, isVerbose, fout );
}
