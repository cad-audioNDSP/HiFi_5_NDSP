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
  NatureDSP Signal Processing Library. Matrix Operations
  Matrix multiply
  * Optimized code for HiFi4
  IntegrIT, 2006-2019
  */

/* Code optimized for HiFi5 core */

/* Cross-platform data type definitions. */
#include "NatureDSP_types.h"
/* Signal Processing Library API. */
#include "NatureDSP_Signal_matop.h"
/* Common helper macros. */
#include "common.h"
#include "common_fpu.h"

#define SZ_F32 (sizeof(float32_t))

#if (HAVE_FPU==0 && HAVE_VFPU==0)
DISCARD_FUN(void,mtx_mpytf_fast,( void* pScr, float32_t * z, const float32_t * x,  const float32_t * yt, int M, int N, int P ))
#elif (HAVE_VFPU)
/*-------------------------------------------------------------------------
  Matrix Multiply
  These functions compute the expression z = 2^lsh * x * y for the matrices 
  x and y. The columnar dimension of x must match the row dimension of y. 
  The resulting matrix has the same number of rows as x and the same number 
  of columns as y.
  Transposing API allows to interpret input yt as transposed matrix y.

  NOTE: lsh factor is not relevant for floating point routines.

  Functions require scratch memory for storing intermediate data. This 
  scratch memory area should be aligned on 16 byte boundary and its size is 
  calculated by dedicated scratch allocation functions.

  Two versions of functions available: regular version (mtx_mpy[t]32x32, 
  mtx_mpy[t]16x16, mtx_mpy[t]8x16, mtx_mpy[t]8x8, mtx[t]_mpyf) with 
  arbitrary arguments and faster version (mtx_mpy[t]32x32_fast, 
  mtx_mpy[t]16x16_fast, mtx_mpy[t]8x16_fast, mtx_mpy[t]8x8_fast, 
  mtx_mpy[t]f_fast, cntx_mpyt32x32_fast, cntx_mpytf_fast) that apply 
  some restrictions

  Precision:
  32x32 32-bit inputs, 32-bit output
  16x16 16-bit inputs, 16-bit output
  8x8   8-bit inputs, 8-bit output
  8x16  8/16-bit inputs, 16-bit output
  f     floating point

  Input:
  x[M*N]      input matrix x, Q7, Q15, Q31 or floating point
  y[N*P]      input matrix y, Q7, Q15, Q31 or floating point
  yt[P*N]     transposed input matrix y. Q31,Q15, Q7 floating point. (for 
              transposing API only)
  M           number of rows in matrix x and z
  N           number of columns in matrix x and number of rows in matrix y
  P           number of columns in matrices y and z
  lsh         left shift applied to the result (applied to the fixed-
              point functions only) 
  Output:
  z[M*P]      output matrix z, Q7, Q15, Q31 or floating point 
  Scratch:
  pScr        size in bytes defined by corresponding scratch allocation 
              functions

  Restrictions:
  For regular routines mpy[t]32x32, mtx_mpy[t]16x16, mtx_mpy[t]8x8, 
  mtx_mpy[t]8x16, mtx_mpy[t]f):
  pScr    aligned on 16-byte boundary
  x,y,z   should not overlap

  For faster routines (mtx_mpy[t]32x32_fast, mtx_mpy[t]16x16_fast, 
  mtx_mpy[t]8x8_fast, mtx_mpy[t]8x16_fast, 
  mtx_mpy[t]f_fast):
  x,y,z       should not overlap
  x,y,z,pScr  aligned on 16-byte boundary
  M,N,P       multiplies of 4 for mtx_mpy[t]32x32_fast, mtx_mpy[t]16x16_fast, 
              mtx_mpy[t]8x8_fast, mtx_mpy[t]8x16_fast, mtx_mpy[t]f_fast
              multiplies of 32 for cntx_mpyt32x32_fast, cntx_mpytf_fast
  lsh         should be in range:
              -31...31 for mtx_mpy32x32, mtx_mpy32x32_fast, cntx_mpyt32x32_fast, 
                       cntx_mpytf_fast
              -15...15 for mtx_mpy16x16, mtx_mpy16x16_fast, mtx_mpy[t]8x8, 
                       mtx_mpy[t]8x8_fast, mtx_mpy[t]8x16, 
                       mtx_mpy[t]8x16_fast 

-------------------------------------------------------------------------*/
void mtx_mpytf_fast(void* pScr, float32_t * z, const float32_t * x, const float32_t * yt, int M, int N, int P)
{
    const xtfloatx4 *restrict pX0;
    const xtfloatx4 *restrict pY0;
          xtfloatx4 *restrict pZ;

    xtfloatx2 Y00, Y01, Y10, Y11;
    xtfloatx2 Y20, Y21, Y30, Y31;
    xtfloatx2 X00, X01, X10, X11;
    xtfloatx2 X20, X21, X30, X31;
    xtfloatx2 acc00, acc01, acc02, acc03;
    xtfloatx2 acc10, acc11, acc12, acc13;
    xtfloatx2 acc20, acc21, acc22, acc23;
    xtfloatx2 acc30, acc31, acc32, acc33;
    xtfloatx2 z00, z01, z02, z03;
    xtfloatx2 z10, z11, z12, z13;
    int m, n, p;

    NASSERT(x);
    NASSERT(yt);
    NASSERT(z);
    NASSERT((z != x) && (z != yt));
    NASSERT(N % 4 == 0);
    NASSERT(M % 4 == 0);
    NASSERT(P % 4 == 0);
    NASSERT_ALIGN(x, HIFI_SIMD_WIDTH);
    NASSERT_ALIGN(yt, HIFI_SIMD_WIDTH);
    NASSERT_ALIGN(z, HIFI_SIMD_WIDTH);

    if ((M <= 0) || (P <= 0)) return;
    if (N <= 0)
    {
        z00 = (xtfloatx2)0.0f;
        pZ = (xtfloatx4 *)(z);
        for (m = 0; m < ((M * P) >> 2); m++)
        {
            AE_SSX2X2_IP(z00, z00, pZ, SZ_F32 * 4);
        }
        return;
    }
    __Pragma("loop_count min=1");
    for (p = 0; p < (P>>2); p++)
    {
        pZ  = (xtfloatx4 *)(z);
        pX0 = (const xtfloatx4 *)(x);

        __Pragma("loop_count min=1");
        for (m = 0; m < (M>>2); m++)
        {
            pY0 = (const xtfloatx4 *)yt;

            acc00 = acc01 = acc02 = acc03 = 0.0f;
            acc10 = acc11 = acc12 = acc13 = 0.0f;
            acc20 = acc21 = acc22 = acc23 = 0.0f;
            acc30 = acc31 = acc32 = acc33 = 0.0f;

            __Pragma("loop_count min=1");
            __Pragma("no_unroll");
            for (n = 0; n < (N >> 2); n++)
            {
                AE_LSX2X2_XP(Y00, Y01, pY0, N*SZ_F32);
                AE_LSX2X2_XP(Y10, Y11, pY0, N*SZ_F32);
                AE_LSX2X2_XP(Y20, Y21, pY0, N*SZ_F32);
                AE_LSX2X2_XP(Y30, Y31, pY0, (4-3*N)*SZ_F32);

                AE_LSX2X2_XP(X00, X01, pX0, N*SZ_F32);
                AE_LSX2X2_XP(X10, X11, pX0, N*SZ_F32);
                AE_LSX2X2_XP(X20, X21, pX0, N*SZ_F32);
                AE_LSX2X2_XP(X30, X31, pX0, (4-3*N)*SZ_F32);

                MADDQ_S(acc00, acc01, Y00, Y10, X00);
                MADDQ_S(acc00, acc01, Y01, Y11, X01);
                MADDQ_S(acc02, acc03, Y20, Y30, X00);
                MADDQ_S(acc02, acc03, Y21, Y31, X01);

                MADDQ_S(acc10, acc11, Y00, Y10, X10);
                MADDQ_S(acc10, acc11, Y01, Y11, X11);
                MADDQ_S(acc12, acc13, Y20, Y30, X10);
                MADDQ_S(acc12, acc13, Y21, Y31, X11);

                MADDQ_S(acc20, acc21, Y00, Y10, X20);
                MADDQ_S(acc20, acc21, Y01, Y11, X21);
                MADDQ_S(acc22, acc23, Y20, Y30, X20);
                MADDQ_S(acc22, acc23, Y21, Y31, X21);

                MADDQ_S(acc30, acc31, Y00, Y10, X30);
                MADDQ_S(acc30, acc31, Y01, Y11, X31);
                MADDQ_S(acc32, acc33, Y20, Y30, X30);
                MADDQ_S(acc32, acc33, Y21, Y31, X31);
            }
            z00 = XT_SEL32_HL_SX2(acc00, acc01);
            z01 = XT_SEL32_LH_SX2(acc00, acc01);
            z00 = z00 + z01;
            z10 = XT_SEL32_HL_SX2(acc10, acc11);
            z11 = XT_SEL32_LH_SX2(acc10, acc11);
            z10 = z10 + z11;
            z02 = XT_SEL32_HL_SX2(acc02, acc03);
            z12 = XT_SEL32_LH_SX2(acc02, acc03);
            z01 = z02 + z12;
            z03 = XT_SEL32_HL_SX2(acc12, acc13);
            z13 = XT_SEL32_LH_SX2(acc12, acc13);
            z11 = z03 + z13;

            AE_SSX2X2_XP(z00, z01, pZ, P*SZ_F32);
            AE_SSX2X2_XP(z10, z11, pZ, P*SZ_F32);

            z00 = XT_SEL32_HL_SX2(acc20, acc21);
            z01 = XT_SEL32_LH_SX2(acc20, acc21);
            z00 = z00 + z01;
            z10 = XT_SEL32_HL_SX2(acc30, acc31);
            z11 = XT_SEL32_LH_SX2(acc30, acc31);
            z10 = z10 + z11;
            z02 = XT_SEL32_HL_SX2(acc22, acc23);
            z12 = XT_SEL32_LH_SX2(acc22, acc23);
            z01 = z02 + z12;
            z03 = XT_SEL32_HL_SX2(acc32, acc33);
            z13 = XT_SEL32_LH_SX2(acc32, acc33);
            z11 = z03 + z13;

            AE_SSX2X2_XP(z00, z01, pZ, P*SZ_F32);
            AE_SSX2X2_XP(z10, z11, pZ, P*SZ_F32);

            pX0 = (const xtfloatx4 *)XT_ADDX4(3*N, (uintptr_t)pX0);
        }
        yt = (const float32_t *)XT_ADDX4(4*N, (uintptr_t)yt);
        z  += 4;
    }
} /* mtx_mpytf_fast() */
#elif (HAVE_FPU)
void mtx_mpytf_fast( void* pScr, float32_t * z, const float32_t * x,  const float32_t * yt, int M, int N, int P )
{
    xtfloat * restrict pZ;
    const xtfloat * restrict pX;
    const xtfloat * restrict pY0;
    const xtfloat * restrict pY1;
    const xtfloat * restrict pY2;
    const xtfloat * restrict pY3;
    int m,n,p;
    NASSERT(x);
    NASSERT(yt);
    NASSERT(z);
    NASSERT((z != x) && (z != yt));
    NASSERT(N % 4 == 0);
    NASSERT(M % 4 == 0);
    NASSERT(P % 4 == 0);
    NASSERT_ALIGN(x, 8);
    NASSERT_ALIGN(yt, 8);
    NASSERT_ALIGN(z, 8);

    if ((M <= 0) || (P <= 0)) return;
    if (N <= 0)
    {
        pZ = (xtfloat *)(z);
        for (m = 0; m < (M * P); m++)
        {
            XT_SSIP(XT_CONST_S(0), pZ, sizeof(xtfloat));
        }
        return;
    }
    pZ = (xtfloat  *)(z);
    for (m = 0; m<M; m++)
    {
        for (p = 0; p<P; p+=4)
        {
            xtfloat A0,A1,A2,A3,x0,y0,y1,y2,y3;
            pX = (const xtfloat  *)(x+m*N);
            pY0 = (const xtfloat  *)(yt+p*N);
            pY1 = (const xtfloat  *)(pY0+N);
            pY2 = (const xtfloat  *)(pY1+N);
            pY3 = (const xtfloat  *)(pY2+N);
            A0=A1=A2=A3=XT_CONST_S(0);
            for ( n=0; n<N; n++ )
            {
                XT_LSIP(x0, pX, sizeof(xtfloat));
                XT_LSIP(y0, pY0, sizeof(xtfloat));
                XT_LSIP(y1, pY1, sizeof(xtfloat));
                XT_LSIP(y2, pY2, sizeof(xtfloat));
                XT_LSIP(y3, pY3, sizeof(xtfloat));
                XT_MADD_S(A0, x0, y0);
                XT_MADD_S(A1, x0, y1);
                XT_MADD_S(A2, x0, y2);
                XT_MADD_S(A3, x0, y3);
            }
            XT_SSIP(A0, pZ, sizeof(xtfloat));
            XT_SSIP(A1, pZ, sizeof(xtfloat));
            XT_SSIP(A2, pZ, sizeof(xtfloat));
            XT_SSIP(A3, pZ, sizeof(xtfloat));
        }
    }
}

#endif

size_t mtx_mpytf_fast_getScratchSize(int M, int N, int P)
{
    return 0;
} /* mtx_mpytf_fast_getScratchSize() */
