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

void mips_vector1(int isFull, int isVerbose, FILE * fout)
{
    PROFILE_NORMALIZED(isFull, isVerbose, vec_dot16x16,       (mips.inp0.i16,mips.inp1.i16,200           ),fout,"N=200",prf_cyclespts,200);
#if 0// HiFi3/3z API
    PROFILE_NORMALIZED(isFull, isVerbose, vec_dot24x24,       (mips.inp0.i32, mips.inp1.i32 , 200        ),fout,"N=200",prf_cyclespts,200);
#endif
    PROFILE_NORMALIZED(isFull, isVerbose, vec_dot32x16,       (mips.inp0.i32,mips.inp1.i16  ,200         ),fout,"x aligned, N=200",prf_cyclespts,200);
    PROFILE_NORMALIZED(isFull, isVerbose, vec_dot32x16,       (mips.inp0.i32,mips.inp1.i16+1,200         ),fout,"x unaligned, N=200",prf_cyclespts,200);
    PROFILE_NORMALIZED(isFull, isVerbose, vec_dot32x32,       (mips.inp0.i32, mips.inp1.i32, 200         ),fout,"x aligned, N=200",prf_cyclespts,200);
    PROFILE_NORMALIZED(isFull, isVerbose, vec_dot32x32,       (mips.inp0.i32+1, mips.inp1.i32, 200       ),fout,"x unaligned, N=200",prf_cyclespts,200);
    PROFILE_NORMALIZED(isFull, isVerbose, vec_dot64x32,       (mips.inp0.i64, mips.inp1.i32, 200         ),fout,"N=200",prf_cyclespts,200);
    PROFILE_NORMALIZED(isFull, isVerbose, vec_dot64x64,       (mips.inp0.i64, mips.inp1.i64, 200         ),fout,"N=200",prf_cyclespts,200);
    PROFILE_NORMALIZED(isFull, isVerbose, vec_dot64x64i,      (mips.inp0.i64, mips.inp1.i64, 200         ),fout,"N=200",prf_cyclespts,200);
    PROFILE_NORMALIZED(     1, isVerbose, vec_dot16x16_fast , (mips.inp0.i16,mips.inp1.i16,200           ),fout,"N=200",prf_cyclespts,200);
#if 0// HiFi3/3z API
    PROFILE_NORMALIZED(     1, isVerbose, vec_dot24x24_fast,  (mips.inp0.i32, mips.inp1.i32, 200         ),fout,"N=200",prf_cyclespts,200);
#endif
    PROFILE_NORMALIZED(     1, isVerbose, vec_dot32x16_fast , (mips.inp0.i32,mips.inp1.i16,200           ),fout,"N=200",prf_cyclespts,200);
    PROFILE_NORMALIZED(     1, isVerbose, vec_dot32x32_fast,  (mips.inp0.i32, mips.inp1.i32, 200         ),fout,"N=200",prf_cyclespts,200);
    PROFILE_NORMALIZED(     1, isVerbose, vec_dot64x32_fast,  (mips.inp0.i64, mips.inp1.i32, 200         ),fout,"N=200",prf_cyclespts,200);
    PROFILE_NORMALIZED(     1, isVerbose, vec_dot64x64_fast,  (mips.inp0.i64, mips.inp1.i64, 200         ),fout,"N=200",prf_cyclespts,200);
    PROFILE_NORMALIZED(     1, isVerbose, vec_dot64x64i_fast, (mips.inp0.i64, mips.inp1.i64, 200         ),fout,"N=200",prf_cyclespts,200);

    PROFILE_VEC_DOT_BATCH (1,isVerbose,vec_dot_batch8x8  ,mips.out0.i16,mips.inp0.i8 ,cint8ptr_t ,fout,200,16);
    PROFILE_VEC_DOT_BATCH (1,isVerbose,vec_dot_batch8x16 ,mips.out0.i16,mips.inp0.i8 ,cint16ptr_t,fout,200,16);
    PROFILE_VEC_DOT_BATCH (1,isVerbose,vec_dot_batch16x16,mips.out0.i32,mips.inp0.i16,cint16ptr_t,fout,200,16);
    PROFILE_VEC_DOT_BATCH (1,isVerbose,vec_dot_batch8x8_fast  ,mips.out0.i16,mips.inp0.i8 ,cint8ptr_t ,fout,200,16);
    PROFILE_VEC_DOT_BATCH (1,isVerbose,vec_dot_batch8x16_fast ,mips.out0.i16,mips.inp0.i8 ,cint16ptr_t,fout,200,16);
    PROFILE_VEC_DOT_BATCH (1,isVerbose,vec_dot_batch16x16_fast,mips.out0.i32,mips.inp0.i16,cint16ptr_t,fout,200,16);

    PROFILE_NORMALIZED(isFull, isVerbose, vec_add16x16     ,  (mips.out0.i16,mips.inp0.i16,mips.inp1.i16,200  ),fout,"x aligned, N=200",prf_cyclespts,200);
    PROFILE_NORMALIZED(isFull, isVerbose, vec_add16x16     ,  (mips.out0.i16,mips.inp0.i16+1,mips.inp1.i16,200  ),fout,"x unaligned, N=200",prf_cyclespts,200);
#if 0// HiFi3/3z API
    PROFILE_NORMALIZED(isFull, isVerbose, vec_add24x24,       (mips.out0.i32,mips.inp0.i32,mips.inp1.i32,200  ),fout,"N=200",prf_cyclespts,200);
#endif
    PROFILE_NORMALIZED(isFull, isVerbose, vec_add32x32     ,  (mips.out0.i32,mips.inp0.i32,mips.inp1.i32,200  ),fout,"x aligned, N=200",prf_cyclespts,200);
    PROFILE_NORMALIZED(isFull, isVerbose, vec_add32x32     ,  (mips.out0.i32,mips.inp0.i32+1,mips.inp1.i32,200  ),fout,"x unaligned, N=200",prf_cyclespts,200);
    PROFILE_NORMALIZED(     1, isVerbose, vec_add16x16_fast,  (mips.out0.i16,mips.inp0.i16,mips.inp1.i16,200  ),fout,"N=200",prf_cyclespts,200);
#if 0// HiFi3/3z API
    PROFILE_NORMALIZED(     1, isVerbose, vec_add24x24_fast,  (mips.out0.i32,mips.inp0.i32,mips.inp1.i32, 200 ),fout,"N=200",prf_cyclespts,200);
#endif
    PROFILE_NORMALIZED(     1, isVerbose, vec_add32x32_fast,  (mips.out0.i32,mips.inp0.i32,mips.inp1.i32,200  ),fout,"N=200",prf_cyclespts,200);
    PROFILE_NORMALIZED(isFull, isVerbose, vec_power16x16,     (mips.inp0.i16,24,200                 ),fout,"x aligned, N=200",prf_cyclespts,200);
    PROFILE_NORMALIZED(isFull, isVerbose, vec_power16x16,     (mips.inp0.i16+1,24,200               ),fout,"x unaligned, N=200",prf_cyclespts,200);
#if 0// HiFi3/3z API
    PROFILE_NORMALIZED(isFull, isVerbose, vec_power24x24,     (mips.inp0.i32, 44, 200               ),fout,"N=200",prf_cyclespts,200);
#endif
    PROFILE_NORMALIZED(isFull, isVerbose, vec_power32x32,     (mips.inp0.i32,44,200                 ),fout,"x aligned, N=200",prf_cyclespts,200);
    PROFILE_NORMALIZED(isFull, isVerbose, vec_power32x32,     (mips.inp0.i32+1,44,200               ),fout,"x unaligned, N=200",prf_cyclespts,200);
    PROFILE_NORMALIZED(     1, isVerbose, vec_power16x16_fast,(mips.inp0.i16,24,200                 ),fout,"N=200",prf_cyclespts,200);
#if 0// HiFi3/3z API
    PROFILE_NORMALIZED(     1, isVerbose, vec_power24x24_fast,(mips.inp0.i32, 44, 200               ),fout,"N=200",prf_cyclespts,200);
#endif
    PROFILE_NORMALIZED(     1, isVerbose, vec_power32x32_fast,(mips.inp0.i32,44,200                 ),fout,"N=200",prf_cyclespts,200);
    PROFILE_NORMALIZED(isFull, isVerbose, vec_shift16x16,     (mips.out0.i16,mips.inp0.i16, 1,200        ),fout,"shift>0, x aligned, N=200",prf_cyclespts,200);
    PROFILE_NORMALIZED(isFull, isVerbose, vec_shift16x16,     (mips.out0.i16+1,mips.inp0.i16+1, 1,200    ),fout,"shift>0, x unaligned, N=200",prf_cyclespts,200);
    PROFILE_NORMALIZED(isFull, isVerbose, vec_shift16x16,     (mips.out0.i16,mips.inp0.i16,-1,200        ),fout,"shift<0, x aligned, N=200",prf_cyclespts,200);
    PROFILE_NORMALIZED(isFull, isVerbose, vec_shift16x16,     (mips.out0.i16+1,mips.inp0.i16+1, -1,200    ),fout,"shift<0, x unaligned, N=200",prf_cyclespts,200);
#if 0// HiFi3/3z API
    PROFILE_NORMALIZED(isFull, isVerbose, vec_shift24x24,     (mips.out0.i32, mips.inp0.i32, 1, 200      ),fout,"N=200",prf_cyclespts,200);
#endif
    PROFILE_NORMALIZED(isFull, isVerbose, vec_shift32x32,     (mips.out0.i32, mips.inp0.i32, 1, 200      ),fout,"x aligned, N=200",prf_cyclespts,200);
    PROFILE_NORMALIZED(isFull, isVerbose, vec_shift32x32,     (mips.out0.i32+1, mips.inp0.i32+1, 1, 200      ),fout,"x unaligned, N=200",prf_cyclespts,200);
    PROFILE_NORMALIZED(isFull, isVerbose, vec_scale16x16,     (mips.out0.i16,mips.inp0.i16,32767,200     ),fout,"x aligned, N=200",prf_cyclespts,200);
    PROFILE_NORMALIZED(isFull, isVerbose, vec_scale16x16,     (mips.out0.i16+1,mips.inp0.i16+1,32767,200     ),fout,"x unaligned, N=200",prf_cyclespts,200);
    PROFILE_NORMALIZED(isFull, isVerbose, vec_scale32x32,     (mips.out0.i32, mips.inp0.i32, 32767, 200  ),fout,"x aligned, N=200",prf_cyclespts,200);
    PROFILE_NORMALIZED(isFull, isVerbose, vec_scale32x32,     (mips.out0.i32, mips.inp0.i32+1, 32767, 200),fout,"x unaligned, N=200",prf_cyclespts,200);
    PROFILE_NORMALIZED(     1, isVerbose, vec_shift16x16_fast,(mips.out0.i16,mips.inp0.i16, 1,200        ),fout,"shift>0, N=200",prf_cyclespts,200);
    PROFILE_NORMALIZED(     1, isVerbose, vec_shift16x16_fast,(mips.out0.i16,mips.inp0.i16,-1,200        ),fout,"shift<0, N=200",prf_cyclespts,200);
    PROFILE_NORMALIZED(     1, isVerbose, vec_shift32x32_fast,(mips.out0.i32, mips.inp0.i32, 1, 200      ),fout,"N=200",prf_cyclespts,200);
    PROFILE_NORMALIZED(     1, isVerbose, vec_scale16x16_fast,(mips.out0.i16,mips.inp0.i16,32767,200     ),fout,"N=200",prf_cyclespts,200);
#if 0// HiFi3/3z API
    PROFILE_NORMALIZED(     1, isVerbose, vec_scale24x24_fast,(mips.out0.i32, mips.inp0.i32, 32767, 200  ),fout,"N=200",prf_cyclespts,200);
    PROFILE_NORMALIZED(     1, isVerbose, vec_scale32x24_fast,(mips.out0.i32, mips.inp0.i32, 32767, 200  ),fout,"N=200",prf_cyclespts,200);
#endif
    PROFILE_NORMALIZED(     1, isVerbose, vec_scale32x32_fast,(mips.out0.i32, mips.inp0.i32, 32767, 200  ),fout,"N=200",prf_cyclespts,200);
#if 0// HiFi3/3z API
    PROFILE_NORMALIZED(     1, isVerbose, vec_shift24x24_fast,(mips.out0.i32, mips.inp0.i32, 1, 200      ),fout,"N=200",prf_cyclespts,200);
#endif
    PROFILE_NORMALIZED(isFull, isVerbose, vec_max16x16,       (mips.inp0.i16,200                    ),fout,"x aligned, N=200",prf_cyclespts,200);
    PROFILE_NORMALIZED(isFull, isVerbose, vec_max16x16,       (mips.inp0.i16+1,200                    ),fout,"x unaligned, N=200",prf_cyclespts,200);
    PROFILE_NORMALIZED(isFull, isVerbose, vec_min16x16,       (mips.inp0.i16,200                    ),fout,"x aligned, N=200",prf_cyclespts,200);
    PROFILE_NORMALIZED(isFull, isVerbose, vec_min16x16,       (mips.inp0.i16+1,200                    ),fout,"x unaligned, N=200",prf_cyclespts,200);
#if 0// HiFi3/3z API
    PROFILE_NORMALIZED(isFull, isVerbose, vec_max24x24,       (mips.inp0.i32, 200                   ),fout,"N=200",prf_cyclespts,200);
    PROFILE_NORMALIZED(isFull, isVerbose, vec_min24x24,       (mips.inp0.i32, 200                   ),fout,"N=200",prf_cyclespts,200);
#endif
    PROFILE_NORMALIZED(isFull, isVerbose, vec_max32x32,       (mips.inp0.i32, 200                   ),fout,"N=200",prf_cyclespts,200);
    PROFILE_NORMALIZED(isFull, isVerbose, vec_min32x32,       (mips.inp0.i32, 200                   ),fout,"N=200",prf_cyclespts,200);
    PROFILE_NORMALIZED(     1, isVerbose, vec_max16x16_fast,  (mips.inp0.i16, 200                   ),fout,"N=200",prf_cyclespts,200);
    PROFILE_NORMALIZED(     1, isVerbose, vec_min16x16_fast,  (mips.inp0.i16, 200                   ),fout,"N=200",prf_cyclespts,200);
#if 0// HiFi3/3z API
    PROFILE_NORMALIZED(     1, isVerbose, vec_max24x24_fast,  (mips.inp0.i32, 200                   ),fout,"N=200",prf_cyclespts,200);
    PROFILE_NORMALIZED(     1, isVerbose, vec_min24x24_fast,  (mips.inp0.i32, 200                   ),fout,"N=200",prf_cyclespts,200);
#endif
    PROFILE_NORMALIZED(     1, isVerbose, vec_max32x32_fast,  (mips.inp0.i32, 200                   ),fout,"N=200",prf_cyclespts,200);
    PROFILE_NORMALIZED(     1, isVerbose, vec_min32x32_fast,  (mips.inp0.i32, 200                   ),fout,"N=200",prf_cyclespts,200);
    PROFILE_NORMALIZED(     1, isVerbose, vec_bexp16,  (mips.inp0.i16, 200                          ),fout,"N=200",prf_cyclespts,200);
    PROFILE_NORMALIZED(     1, isVerbose, vec_bexp32,  (mips.inp0.i32, 200                          ),fout,"N=200",prf_cyclespts,200);
    PROFILE_NORMALIZED(     1, isVerbose, vec_bexp16_fast,  (mips.inp0.i16, 200                     ),fout,"N=200",prf_cyclespts,200);
    PROFILE_NORMALIZED(     1, isVerbose, vec_bexp32_fast,  (mips.inp0.i32, 200                     ),fout,"N=200",prf_cyclespts,200);
    PROFILE_SIMPLE(     1, isVerbose, scl_bexp16,(3917),                           fout,"",prf_cycle);
    PROFILE_SIMPLE(     1, isVerbose, scl_bexp32,(9621325),                        fout,"",prf_cycle);
}
