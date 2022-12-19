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

#if (HAVE_FPU==0 && HAVE_VFPU==0)
DISCARD_FUN(void,mtx_mpyf_fast,(void* pScr, float32_t * z, const float32_t * x,  const float32_t * y, int M, int N, int P ))
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
#define SZ_F32 (sizeof(float32_t))

void mtx_mpyf_fast(void* pScr, float32_t * z, const float32_t * x, const float32_t * y, int M, int N, int P)
{
    const xtfloatx4 *restrict pX0;
    const xtfloatx4 *restrict pX1;
    const xtfloatx4 *restrict pX2;
    const xtfloatx4 *restrict pX3;
    const xtfloatx4 *restrict pY;
          xtfloatx4 *restrict pZ;

    xtfloatx2 Y0, Y1, Y2, Y3, Y4, Y5, Y6, Y7;
    xtfloatx2 X00, X01, X10, X11;
    xtfloatx2 X20, X21, X30, X31;
    xtfloatx2 acc00, acc01, acc02, acc03, acc10, acc11, acc12, acc13;
    xtfloatx2 acc20, acc21, acc22, acc23, acc30, acc31, acc32, acc33;
    xtfloatx2 z00, z10, z01, z11, z02, z12, z03, z13;
    ae_int16x4 C0, C1, C2, C3, ind;
    int m, n, p;
    static const short ALIGN(16) dsel_ind[4] = { 1797, 1540, 769, 512 };

    (void)pScr;
    NASSERT(x);
    NASSERT(y);
    NASSERT(z);
    NASSERT((z != x) && (z != y));
    NASSERT(N % 4 == 0);
    NASSERT(M % 4 == 0);
    NASSERT(P % 4 == 0);
    NASSERT_ALIGN(x, HIFI_SIMD_WIDTH);
    NASSERT_ALIGN(y, HIFI_SIMD_WIDTH);
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

    ind = AE_L16X4_I((ae_int16x4*)(dsel_ind), 0);

    __Pragma("loop_count min=1");
    for (p = 0; p < P; p += 4)
    {
        pZ  = (xtfloatx4 *)(z);
        pX0 = (const xtfloatx4 *)(x);

        __Pragma("loop_count min=1");
        for (m = 0; m < (M>>2); m++)
        {
            pX1 = (const xtfloatx4 *)((float32_t *)pX0 + N);
            pX2 = (const xtfloatx4 *)((float32_t *)pX1 + N);
            pX3 = (const xtfloatx4 *)((float32_t *)pX2 + N);
            pY = (const xtfloatx4 *)(y);

            acc00 = acc01 = acc02 = acc03 = 0.0f;
            acc10 = acc11 = acc12 = acc13 = 0.0f;
            acc20 = acc21 = acc22 = acc23 = 0.0f;
            acc30 = acc31 = acc32 = acc33 = 0.0f;
            __Pragma("loop_count min=1");
            for (n = 0; n < (N >> 2); n++)
            {
                AE_LSX2X2_IP(X00, X01, pX0, 4*SZ_F32);
                AE_LSX2X2_IP(X10, X11, pX1, 4*SZ_F32);
                AE_LSX2X2_IP(X20, X21, pX2, 4*SZ_F32);
                AE_LSX2X2_IP(X30, X31, pX3, 4*SZ_F32);

                AE_LSX2X2_XP(Y0, Y1, pY, P*SZ_F32);
                AE_LSX2X2_XP(Y2, Y3, pY, P*SZ_F32);
                AE_LSX2X2_XP(Y4, Y5, pY, P*SZ_F32);
                AE_LSX2X2_XP(Y6, Y7, pY, P*SZ_F32);

                AE_DSEL16X4(C0, C1, AE_MOVINT16X4_FROMINT32X2(AE_MOVINT32X2_FROMXTFLOATX2(Y0)), AE_MOVINT16X4_FROMINT32X2(AE_MOVINT32X2_FROMXTFLOATX2(Y2)), ind);
                AE_DSEL16X4(C2, C3, AE_MOVINT16X4_FROMINT32X2(AE_MOVINT32X2_FROMXTFLOATX2(Y1)), AE_MOVINT16X4_FROMINT32X2(AE_MOVINT32X2_FROMXTFLOATX2(Y3)), ind);
                Y0 = AE_MOVXTFLOATX2_FROMINT32X2(AE_MOVINT32X2_FROMINT16X4(C0));
                Y2 = AE_MOVXTFLOATX2_FROMINT32X2(AE_MOVINT32X2_FROMINT16X4(C1));
                Y1 = Y2;
                Y2 = AE_MOVXTFLOATX2_FROMINT32X2(AE_MOVINT32X2_FROMINT16X4(C2));
                Y3 = AE_MOVXTFLOATX2_FROMINT32X2(AE_MOVINT32X2_FROMINT16X4(C3));

                AE_DSEL16X4(C0, C1, AE_MOVINT16X4_FROMINT32X2(AE_MOVINT32X2_FROMXTFLOATX2(Y4)), AE_MOVINT16X4_FROMINT32X2(AE_MOVINT32X2_FROMXTFLOATX2(Y6)), ind);
                AE_DSEL16X4(C2, C3, AE_MOVINT16X4_FROMINT32X2(AE_MOVINT32X2_FROMXTFLOATX2(Y5)), AE_MOVINT16X4_FROMINT32X2(AE_MOVINT32X2_FROMXTFLOATX2(Y7)), ind);
                Y4 = AE_MOVXTFLOATX2_FROMINT32X2(AE_MOVINT32X2_FROMINT16X4(C0));
                Y6 = AE_MOVXTFLOATX2_FROMINT32X2(AE_MOVINT32X2_FROMINT16X4(C1));
                Y5 = Y6;
                Y6 = AE_MOVXTFLOATX2_FROMINT32X2(AE_MOVINT32X2_FROMINT16X4(C2));
                Y7 = AE_MOVXTFLOATX2_FROMINT32X2(AE_MOVINT32X2_FROMINT16X4(C3));

                MADDQ_S(acc00, acc10, X00, X10, Y0);
                MADDQ_S(acc00, acc10, X01, X11, Y4);
                MADDQ_S(acc01, acc11, X00, X10, Y1);
                MADDQ_S(acc01, acc11, X01, X11, Y5);
                MADDQ_S(acc02, acc12, X00, X10, Y2);
                MADDQ_S(acc02, acc12, X01, X11, Y6);
                MADDQ_S(acc03, acc13, X00, X10, Y3);
                MADDQ_S(acc03, acc13, X01, X11, Y7);

                MADDQ_S(acc20, acc30, X20, X30, Y0);
                MADDQ_S(acc20, acc30, X21, X31, Y4);
                MADDQ_S(acc21, acc31, X20, X30, Y1);
                MADDQ_S(acc21, acc31, X21, X31, Y5);
                MADDQ_S(acc22, acc32, X20, X30, Y2);
                MADDQ_S(acc22, acc32, X21, X31, Y6);
                MADDQ_S(acc23, acc33, X20, X30, Y3);
                MADDQ_S(acc23, acc33, X21, X31, Y7);
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

            pX0 = (const xtfloatx4 *)pX3;
        }
        y += 4;
        z += 4;
    }
} /* mtx_mpyf_fast() */
#elif (HAVE_FPU)
void mtx_mpyf_fast(void* pScr, float32_t * z, const float32_t * x,  const float32_t * y, int M, int N, int P )
{
    xtfloat * restrict pZ;
    const xtfloat * restrict pX;
    const xtfloat * restrict pY;
    int m,n,p;
    (void)pScr;
    NASSERT(x);
    NASSERT(y);
    NASSERT(z);
    NASSERT((z != x) && (z != y));
    NASSERT(N % 4 == 0);
    NASSERT(M % 4 == 0);
    NASSERT(P % 4 == 0);
    NASSERT_ALIGN(x, 8);
    NASSERT_ALIGN(y, 8);
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
            pY = (const xtfloat  *)(y+p);
            A0=A1=A2=A3=XT_CONST_S(0);
            for ( n=0; n<N; n++ )
            {
                XT_LSIP(x0, pX, sizeof(xtfloat));
                y1=XT_LSI(pY, 1*sizeof(xtfloat));
                y2=XT_LSI(pY, 2*sizeof(xtfloat));
                y3=XT_LSI(pY, 3*sizeof(xtfloat));
                XT_LSXP(y0, pY, P*sizeof(xtfloat));
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


size_t mtx_mpyf_fast_getScratchSize(int M, int N, int P)
{
    (void)M; (void)N; (void)P;
    return 0;
}
