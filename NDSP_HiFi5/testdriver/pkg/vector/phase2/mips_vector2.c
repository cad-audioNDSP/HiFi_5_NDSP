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
* Test module for testing cycle performance (Vector Operations)
*/
#include "types.h"
#include "config.h"
#include "packages.h"
#include LIBRARY_HEADER(vector)
#include "mips.h"

typedef void tFxn(void   *z, const void* x, const void* y, int rsh, int N, int M);
#define PROFILE_VEC_DOT_BATCH(isFull,isVerbose,fn,zfmt,xfmt,yptr,fout,N,M)      \
{																				\
    uintptr_t y[20];															\
    int m;																		\
    NASSERT(M<=sizeof(y)/sizeof(y[0]));											\
    for (m=0; m<M; m++) y[m]=((uintptr_t)mips.inp1.i16)+16*m;					\
    PROFILE_NORMALIZED(isFull, isVerbose, fn,( zfmt,  xfmt,						\
                       (const yptr*)y,0,N,M),									\
                       fout,"N=" #N ", M=" #M,prf_cyclespts,M*N);				\
}

#define PROFILE_VEC_DOT_BATCHF(isFull,isVerbose,fn,zfmt,xfmt,yptr,fout,N,M)     \
{																				\
    uintptr_t y[20];															\
    int m;																		\
    NASSERT(M<=sizeof(y)/sizeof(y[0]));											\
    for (m=0; m<M; m++) y[m]=((uintptr_t)mips.inp1.i16)+16*m;					\
    PROFILE_NORMALIZED(isFull, isVerbose, fn,( zfmt, xfmt,						\
                       (const yptr*)y,N,M),										\
                       fout,"N=" #N ", M=" #M,prf_cyclespts,M*N);				\
}

void mips_vector2(int isFull, int isVerbose, FILE * fout)
{
    PROFILE_NORMALIZED(    1, isVerbose, vec_dotf,       (mips.inp0.f32,mips.inp1.f32,200           ),fout,"N=200",prf_cyclespts,200);
    PROFILE_VEC_DOT_BATCHF(1,isVerbose,vec_dot_batchf       ,mips.out0.f32,mips.inp0.f32,cfloat32ptr_t,fout,200,16);
    PROFILE_VEC_DOT_BATCHF(1,isVerbose,vec_dot_batchf_fast  ,mips.out0.f32,mips.inp0.f32,cfloat32ptr_t,fout,200,16);
    PROFILE_NORMALIZED(    1, isVerbose, vec_addf,       (mips.out0.f32,mips.inp0.f32,mips.inp1.f32,200  ),fout,"N=200",prf_cyclespts,200);
    PROFILE_NORMALIZED(    1, isVerbose, vec_powerf,     (mips.inp0.f32,200                    ),fout,"N=200",prf_cyclespts,200);
    PROFILE_NORMALIZED(    1, isVerbose, vec_shiftf,     (mips.out0.f32, mips.inp0.f32, 1, 200      ),fout,"N=200",prf_cyclespts,200);
    PROFILE_NORMALIZED(    1, isVerbose, vec_scalef,     (mips.out0.f32, mips.inp0.f32, 1.   , 200  ),fout,"N=200",prf_cyclespts,200);
    PROFILE_NORMALIZED(    1, isVerbose, vec_scale_sf,   (mips.out0.f32, mips.inp0.f32, 32767 ,- 1000, 1000, 200), fout,"N=200",prf_cyclespts,200);
    PROFILE_NORMALIZED(    1, isVerbose, vec_minf,       (mips.inp0.f32,200                    ),fout,"N=200",prf_cyclespts,200);
    PROFILE_NORMALIZED(    1, isVerbose, vec_maxf,       (mips.inp0.f32, 200                   ),fout,"N=200",prf_cyclespts,200);
    PROFILE_NORMALIZED(    1, isVerbose, vec_bexpf,      (mips.inp0.f32, 200                   ),fout,"N=200",prf_cyclespts,200);
    PROFILE_SIMPLE(     1, isVerbose, scl_bexpf,(9621325.f),      fout,""     ,prf_cycle);
}

