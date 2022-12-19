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
  NatureDSP Signal Processing Library. FIR part
    Real data circular convolution/correlation, floating point, helper file
    C code optimized for HiFi5
  IntegrIT, 2006-2019
*/

/* Portable data types. */
#include "NatureDSP_types.h"
/* Common utility and macros declarations. */
#include "common.h"
#include "raw_corrf.h"
#include "common_fpu.h"

#if (HAVE_VFPU)

#define AE_DSELSX2(out0, out1, in0, in1, sel_mask)                       \
{                                                                        \
    ae_int16x4 tmp0, tmp1;                                               \
    tmp0 = AE_MOVINT16X4_FROMINT32X2(AE_MOVINT32X2_FROMXTFLOATX2(in0));  \
    tmp1 = AE_MOVINT16X4_FROMINT32X2(AE_MOVINT32X2_FROMXTFLOATX2(in1));  \
    AE_DSEL16X4(tmp0, tmp1, tmp0, tmp1, sel_mask);                       \
    out0 = AE_MOVXTFLOATX2_FROMINT32X2(AE_MOVINT32X2_FROMINT16X4(tmp0)); \
    out1 = AE_MOVXTFLOATX2_FROMINT32X2(AE_MOVINT32X2_FROMINT16X4(tmp1)); \
}

/*
    raw linear correlation:
    Restrictions:
    x[N+M-1] - aligned on 16 byte boundary
    y[M]     - aligned on 16 byte boundary
    actual length of vector y must be a multiple of 4
    with zeroed excess values;
    N>0 && M>0
    N>=M-1
*/
void raw_corrf( float32_t * restrict r,
          const float32_t * restrict x,
          const float32_t * restrict y,
          int N, int M )
{
    //
    // Circular cross-correlation algorithm:
    //
    //   r[n] = sum( x[mod(n+m,N)]*y[m] )
    //        m=0..M-1
    //
    //   where n = 0..N-1
    //
    const xtfloatx4 * restrict pX0;
    const xtfloatx4 * restrict pX1;
    const xtfloatx4 * restrict pX2;
    const xtfloatx4 * restrict pX3;
    const xtfloatx4 * restrict pX4;
    const xtfloatx4 * restrict pY;
          xtfloatx4 * restrict pR;
    ae_valignx2 aR, aX3, aX4;
    xtfloatx2 x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, xa;
    xtfloatx2 y0, y1, y2, y3;
    xtfloatx2 r0, r1, r2, r3;
    int n, m;

    static const ALIGN(16) int16_t Sel[4] = { 0x0705, 0x0604, 0x0301, 0x0200 };
    ae_int16x4 sel;
    sel = AE_L16X4_I((ae_int16x4*)&Sel, 0);

    NASSERT(r);
    NASSERT(x);
    NASSERT(y);
    NASSERT_ALIGN(x, 16);
    NASSERT_ALIGN(y, 16);
    NASSERT(N > 0 && M > 0);
    NASSERT(N >= M - 1);

    M = (M + 3)&~3;

    pR = (xtfloatx4 *)r;
    aR = AE_ZALIGN128();
    for (n = 0; n < (N&~7); n += 8)
    {
        xtfloatx2 r00, r10, r20, r30, r40, r50, r60, r70;
        xtfloatx2 r01, r11, r21, r31, r41, r51, r61, r71;

        pY = (const xtfloatx4 *)(y);
        pX0 = (const xtfloatx4 *)(x + n);
        pX1 = (const xtfloatx4 *)(x + n + 4);
        pX2 = (const xtfloatx4 *)(x + n + 8);
        pX3 = (const xtfloatx4 *)(x + n + 3);
        pX4 = (const xtfloatx4 *)(x + n + 7);
        aX3 = AE_LA128_PP(pX3);
        aX4 = AE_LA128_PP(pX4);

        CONST_SX2X2(r00, r01, 0); CONST_SX2X2(r10, r11, 0);
        CONST_SX2X2(r20, r21, 0); CONST_SX2X2(r30, r31, 0);
        CONST_SX2X2(r40, r41, 0); CONST_SX2X2(r50, r51, 0);
        CONST_SX2X2(r60, r61, 0); CONST_SX2X2(r70, r71, 0);

        /* preload data from x */
        AE_LSX2X2_IP(x0, x2, pX0, 4 * sizeof(float32_t));
        x1 = XT_SEL32_LH_SX2(x0, x2);

        __Pragma("loop_count min=1");
        for (m = 0; m < (M >> 2); m++)
        {
            /* load data from x */
            AE_LSX2X2_IP(x4, x6, pX1, 4 * sizeof(float32_t));
            AE_LSX2IP(x8, castxcc(xtfloatx2, pX2), 4 * sizeof(float32_t));
            AE_LASX2X2_IP(x3, x5, aX3, pX3);
            AE_LASX2X2_IP(x7, x9, aX4, pX4);
            /* load data from y */
            AE_LSX2X2_IP(y0, y1, pY, 4 * sizeof(float32_t));
            /* compute correlation of 8 values */
            MADD_SX2X2(r00, r10, x0, x1, y0, y0);
            MADD_SX2X2(r20, r30, x2, x3, y0, y0);
            MADD_SX2X2(r40, r50, x4, x5, y0, y0);
            MADD_SX2X2(r60, r70, x6, x7, y0, y0);
            MADD_SX2X2(r01, r11, x2, x3, y1, y1);
            MADD_SX2X2(r21, r31, x4, x5, y1, y1);
            MADD_SX2X2(r41, r51, x6, x7, y1, y1);
            MADD_SX2X2(r61, r71, x8, x9, y1, y1);
            /* reload data and shift input line for the next iteration */
            AE_LSX2X2_IP(x0, x2, pX0, 4 * sizeof(float32_t));
            x1 = x5;
        }
        /* reduction add */
        ADD_SX2X2(r00, r10, r00, r10, r01, r11);
        ADD_SX2X2(r20, r30, r20, r30, r21, r31);
        AE_DSELSX2(r0, r1, r00, r10, sel);
        AE_DSELSX2(r2, r3, r20, r30, sel);
        ADD_SX2X2(r0, r1, r0, r2, r1, r3);
        /* save computed samples */
        AE_SASX2X2_IP(r0, r1, aR, pR);

        /* reduction add */
        ADD_SX2X2(r40, r50, r40, r50, r41, r51);
        ADD_SX2X2(r60, r70, r60, r70, r61, r71);
        AE_DSELSX2(r0, r1, r40, r50, sel);
        AE_DSELSX2(r2, r3, r60, r70, sel);
        ADD_SX2X2(r0, r1, r0, r2, r1, r3);
        /* save computed samples */
        AE_SASX2X2_IP(r0, r1, aR, pR);
    }
    if (N & 4)
    {
        xtfloatx2 r00, r10, r20, r30;
        xtfloatx2 r01, r11, r21, r31;
        xtfloatx2 r02, r12, r22, r32;
        xtfloatx2 r03, r13, r23, r33;

        pY = (const xtfloatx4 *)(y);
        pX0 = (const xtfloatx4 *)(x + n);
        pX3 = (const xtfloatx4 *)(x + n + 3);
        aX3 = AE_LA128_PP(pX3);

        CONST_SX2X2(r00, r01, 0); CONST_SX2X2(r02, r03, 0);
        CONST_SX2X2(r10, r11, 0); CONST_SX2X2(r12, r13, 0);
        CONST_SX2X2(r20, r21, 0); CONST_SX2X2(r22, r23, 0);
        CONST_SX2X2(r30, r31, 0); CONST_SX2X2(r32, r33, 0);

        /* preload data from x */
        AE_LSX2X2_IP(x0, x2, pX0, 4 * sizeof(float32_t));
        x1 = XT_SEL32_LH_SX2(x0, x2);

        for (m = 0; m < (M >> 3); m++)
        {
            /* load data from x */
            AE_LSX2X2_IP(x4, x6, pX0, 4 * sizeof(float32_t));
            AE_LSX2X2_IP(x8, xa, pX0, 4 * sizeof(float32_t));
            AE_LASX2X2_IP(x3, x5, aX3, pX3);
            AE_LASX2X2_IP(x7, x9, aX3, pX3);
            /* load data from y */
            AE_LSX2X2_IP(y0, y1, pY, 4 * sizeof(float32_t));
            AE_LSX2X2_IP(y2, y3, pY, 4 * sizeof(float32_t));
            /* compute correlation of 4 values */
            MADD_SX2X2(r00, r10, x0, x1, y0, y0);
            MADD_SX2X2(r20, r30, x2, x3, y0, y0);
            MADD_SX2X2(r01, r11, x2, x3, y1, y1);
            MADD_SX2X2(r21, r31, x4, x5, y1, y1);
            MADD_SX2X2(r02, r12, x4, x5, y2, y2);
            MADD_SX2X2(r22, r32, x6, x7, y2, y2);
            MADD_SX2X2(r03, r13, x6, x7, y3, y3);
            MADD_SX2X2(r23, r33, x8, x9, y3, y3);
            /* shift input line for the next iteration */
            x0 = x8;
            x1 = x9;
            x2 = xa;
        }
        if (M & 4)
        {
            /* load data from x */
            AE_LSX2X2_IP(x4, x6, pX0, 4 * sizeof(float32_t));
            AE_LASX2X2_IP(x3, x5, aX3, pX3);
            /* load data from y */
            AE_LSX2X2_IP(y0, y1, pY, 4 * sizeof(float32_t));
            /* compute correlation of 4 values */
            MADD_SX2X2(r00, r10, x0, x1, y0, y0);
            MADD_SX2X2(r20, r30, x2, x3, y0, y0);
            MADD_SX2X2(r01, r11, x2, x3, y1, y1);
            MADD_SX2X2(r21, r31, x4, x5, y1, y1);
        }
        /* reduction add */
        ADD_SX2X2(r00, r10, r00, r10, r01, r11);
        ADD_SX2X2(r20, r30, r20, r30, r21, r31);
        ADD_SX2X2(r02, r12, r02, r12, r03, r13);
        ADD_SX2X2(r22, r32, r22, r32, r23, r33);
        ADD_SX2X2(r00, r10, r00, r10, r02, r12);
        ADD_SX2X2(r20, r30, r20, r30, r22, r32);
        AE_DSELSX2(r0, r1, r00, r10, sel);
        AE_DSELSX2(r2, r3, r20, r30, sel);
        ADD_SX2X2(r0, r1, r0, r2, r1, r3);
        /* save computed samples */
        AE_SASX2X2_IP(r0, r1, aR, pR);
        n += 4;
    }
    AE_SA128POS_FP(aR, pR);

    N &= 3;
    if (N)
    {
        xtfloatx2 r00, r10, r20, r30;
        xtfloatx2 r01, r11, r21, r31;
        xtfloatx2 r02, r12, r22, r32;
        xtfloatx2 r03, r13, r23, r33;

        pY = (const xtfloatx4 *)(y);
        pX0 = (const xtfloatx4 *)(x + n);
        pX3 = (const xtfloatx4 *)(x + n + 3);
        aX3 = AE_LA128_PP(pX3);

        CONST_SX2X2(r00, r01, 0); CONST_SX2X2(r02, r03, 0);
        CONST_SX2X2(r10, r11, 0); CONST_SX2X2(r12, r13, 0);
        CONST_SX2X2(r20, r21, 0); CONST_SX2X2(r22, r23, 0);
        CONST_SX2X2(r30, r31, 0); CONST_SX2X2(r32, r33, 0);

        /* preload data from x */
        AE_LSX2X2_IP(x0, x2, pX0, 4 * sizeof(float32_t));
        x1 = XT_SEL32_LH_SX2(x0, x2);

        for (m = 0; m < (M >> 3); m++)
        {
            /* load data from x */
            AE_LSX2X2_IP(x4, x6, pX0, 4 * sizeof(float32_t));
            AE_LSX2X2_IP(x8, xa, pX0, 4 * sizeof(float32_t));
            AE_LASX2X2_IP(x3, x5, aX3, pX3);
            AE_LASX2X2_IP(x7, x9, aX3, pX3);
            /* load data from y */
            AE_LSX2X2_IP(y0, y1, pY, 4 * sizeof(float32_t));
            AE_LSX2X2_IP(y2, y3, pY, 4 * sizeof(float32_t));
            /* compute correlation of 4 values */
            MADD_SX2X2(r00, r10, x0, x1, y0, y0);
            MADD_SX2X2(r20, r30, x2, x3, y0, y0);
            MADD_SX2X2(r01, r11, x2, x3, y1, y1);
            MADD_SX2X2(r21, r31, x4, x5, y1, y1);
            MADD_SX2X2(r02, r12, x4, x5, y2, y2);
            MADD_SX2X2(r22, r32, x6, x7, y2, y2);
            MADD_SX2X2(r03, r13, x6, x7, y3, y3);
            MADD_SX2X2(r23, r33, x8, x9, y3, y3);
            /* shift input line for the next iteration */
            x0 = x8;
            x1 = x9;
            x2 = xa;
        }
        if (M & 4)
        {
            /* load data from x */
            AE_LSX2X2_IP(x4, x6, pX0, 4 * sizeof(float32_t));
            AE_LASX2X2_IP(x3, x5, aX3, pX3);
            /* load data from y */
            AE_LSX2X2_IP(y0, y1, pY, 4 * sizeof(float32_t));
            /* compute correlation of 4 values */
            MADD_SX2X2(r00, r10, x0, x1, y0, y0);
            MADD_SX2X2(r20, r30, x2, x3, y0, y0);
            MADD_SX2X2(r01, r11, x2, x3, y1, y1);
            MADD_SX2X2(r21, r31, x4, x5, y1, y1);
        }
        /* reduction add */
        ADD_SX2X2(r00, r10, r00, r10, r01, r11);
        ADD_SX2X2(r20, r30, r20, r30, r21, r31);
        ADD_SX2X2(r02, r12, r02, r12, r03, r13);
        ADD_SX2X2(r22, r32, r22, r32, r23, r33);
        ADD_SX2X2(r00, r10, r00, r10, r02, r12);
        ADD_SX2X2(r20, r30, r20, r30, r22, r32);
        AE_DSELSX2(r0, r1, r00, r10, sel);
        AE_DSELSX2(r2, r3, r20, r30, sel);
        ADD_SX2X2(r0, r1, r0, r2, r1, r3);
        /* save computed samples */
        AE_S32_H_I(AE_MOVINT32X2_FROMXTFLOATX2(r0), (ae_int32 *)pR, 0 * sizeof(float32_t));
        if (N > 1) AE_S32_L_I(AE_MOVINT32X2_FROMXTFLOATX2(r0), (ae_int32 *)pR, 1 * sizeof(float32_t));
        if (N > 2) AE_S32_H_I(AE_MOVINT32X2_FROMXTFLOATX2(r1), (ae_int32 *)pR, 2 * sizeof(float32_t));
    }
}/* raw_corrf() */
#elif HAVE_FPU
// for scalar FPU
/*
    raw linear correlation:
    Restrictions:
    x - aligned on 8 byte boundary
    N>0 && M>0
    N>=M-1
*/
void raw_corrf(
                      float32_t  * restrict r,
                const float32_t  * restrict x,
                const float32_t  * restrict y,
                int N, int M )
{
  /*
  * Circular cross-correlation algorithm:
  *
  *   r[n] = sum( x[mod(n+m,N)]*y[m] )
  *        m=0..M-1
  *
  *   where n = 0..N-1
  */

  xtfloat A0,A1,A2,A3,X0,X1,X2,X3,X4,Y0,Y1;
  const xtfloat * restrict pX;
  const xtfloat* restrict pY;
  xtfloat * restrict pR;
  int n, m;

  NASSERT(r);
  NASSERT(x);
  NASSERT(y);
  NASSERT_ALIGN(x, 8);
  NASSERT(N>0 && M>0);
  NASSERT(N >= M - 1);

  pR=(      xtfloat*)r;
  /* set circular buffer boundaries */
  WUR_AE_CBEGIN0( (uintptr_t)( x + 0 ) );
  WUR_AE_CEND0  ( (uintptr_t)( x + N ) );
  for (n = 0; n<(N&~3); n += 4)
  {
    A0=A1=A2=A3=XT_CONST_S(0);
    pX = (const xtfloat*)(x + n);
    pY = (const xtfloat*)y;
    XT_LSXC(X0, pX, 4);
    XT_LSXC(X1, pX, 4);
    XT_LSXC(X2, pX, 4);
    for ( m=0; m<(M>>1); m++ )
    {
      XT_LSXC(X3, pX, 4);
      XT_LSXC(X4, pX, 4);

      XT_LSIP(Y0, pY, 4);
      XT_LSIP(Y1, pY, 4);
      XT_MADD_S(A0, X0, Y0);
      XT_MADD_S(A1, X1, Y0);
      XT_MADD_S(A2, X2, Y0);
      XT_MADD_S(A3, X3, Y0);
      XT_MADD_S(A0, X1, Y1);
      XT_MADD_S(A1, X2, Y1);
      XT_MADD_S(A2, X3, Y1);
      XT_MADD_S(A3, X4, Y1);
      X0 = X2;
      X1 = X3;
      X2 = X4;
    }
    if (M & 1)
    {
      XT_LSXC(X3, pX, 4);
      XT_LSXC(X4, pX, 4);

      XT_LSIP(Y0, pY, 4);
      XT_MADD_S(A0, X0, Y0);
      XT_MADD_S(A1, X1, Y0);
      XT_MADD_S(A2, X2, Y0);
      XT_MADD_S(A3, X3, Y0);
    }
    XT_SSIP(A0,pR,4);
    XT_SSIP(A1,pR,4);
    XT_SSIP(A2,pR,4);
    XT_SSIP(A3,pR,4);
  }
  /* last 1...3 iterations */
  N&=3;
  if (N)
  {
    A0 = A1 = A2 = XT_CONST_S(0);
    pX = (const xtfloat*)(x+n);
    pY = (const xtfloat*)y;
    XT_LSXC(X0, pX, 4);
    XT_LSXC(X1, pX, 4);
    for (m = 0; m<(M >> 1); m++)
    {
      XT_LSXC(X2, pX, 4);
      XT_LSXC(X3, pX, 4);

      XT_LSIP(Y0, pY, 4);
      XT_LSIP(Y1, pY, 4);
      XT_MADD_S(A0, X0, Y0);
      XT_MADD_S(A1, X1, Y0);
      XT_MADD_S(A2, X2, Y0);
      XT_MADD_S(A0, X1, Y1);
      XT_MADD_S(A1, X2, Y1);
      XT_MADD_S(A2, X3, Y1);
      X0 = X2;
      X1 = X3;
    }
    if (M & 1)
    {
      XT_LSXC(X2, pX, 4);
      XT_LSIP(Y0, pY, 4);
      XT_MADD_S(A0, X0, Y0);
      XT_MADD_S(A1, X1, Y0);
      XT_MADD_S(A2, X2, Y0);
    }
    if (N>2) XT_SSI(A2, pR, 2 * 4);
    if (N>1) XT_SSI(A1, pR, 1 * 4);
    if (N>0) XT_SSI(A0, pR, 0 * 4);
  }
} /* raw_corrf() */

#endif
