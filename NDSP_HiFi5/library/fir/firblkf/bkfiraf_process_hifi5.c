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
    Real block FIR filter, floating point
    C code optimized for HiFi5
  IntegrIT, 2006-2019
*/

/* Portable data types. */
#include "NatureDSP_types.h"
/* Signal Processing Library API. */
#include "NatureDSP_Signal_fir.h"
/* Common utility and macros declarations. */
#include "common.h"
#include "common_fpu.h"
#include "bkfiraf.h"

/*-------------------------------------------------------------------------
  Real FIR filter.
  Computes a real FIR filter (direct-form) using IR stored in vector h. The 
  real data input is stored in vector x. The filter output result is stored 
  in vector y. The filter calculates N output samples using M coefficients 
  and requires last M-1 samples in the delay line.
  These functions implement FIR filter with no limitation on size of data
  block, alignment and length of impulse response at the cost of increased
  processing complexity.
  NOTE: 
  User application is not responsible for management of delay lines.

  Precision: 
  16x16    16-bit data, 16-bit coefficients, 16-bit outputs
  32x16    32-bit data, 16-bit coefficients, 32-bit outputs
  32x32    32-bit data, 32-bit coefficients, 32-bit outputs
  32x32ep  the same as above but using 72-bit accumulator for intermediate 
           computations
  f        floating point
  Input:
  x[N]     input samples, Q15, Q31, floating point
  h[M]     filter coefficients in normal order, Q15, Q31, floating point
  N        length of sample block
  M        length of filter
  Output:
  y[N]     input samples, Q15, Q31, floating point 

  Restrictions:
  x,y      should not be overlapping
-------------------------------------------------------------------------*/

#if (HAVE_FPU==0 && HAVE_VFPU==0)
DISCARD_FUN(void,bkfiraf_process,( bkfiraf_handle_t _bkfir,
                         float32_t * restrict  y,
                   const float32_t * restrict  x, int N ))
#elif HAVE_VFPU

#define SMALLER_CODESIZE 1

#define AE_DSELSX2(out0, out1, in0, in1, sel_mask)                       \
{                                                                        \
    ae_int16x4 tmp0, tmp1;                                               \
    tmp0 = AE_MOVINT16X4_FROMINT32X2(AE_MOVINT32X2_FROMXTFLOATX2(in0));  \
    tmp1 = AE_MOVINT16X4_FROMINT32X2(AE_MOVINT32X2_FROMXTFLOATX2(in1));  \
    AE_DSEL16X4(tmp0, tmp1, tmp0, tmp1, sel_mask);                       \
    out0 = AE_MOVXTFLOATX2_FROMINT32X2(AE_MOVINT32X2_FROMINT16X4(tmp0)); \
    out1 = AE_MOVXTFLOATX2_FROMINT32X2(AE_MOVINT32X2_FROMINT16X4(tmp1)); \
}

/* process block of samples */
void bkfiraf_process( bkfiraf_handle_t handle,
                         float32_t * restrict  y,
                   const float32_t * restrict  x, int N )
{
    bkfiraf_t * bkfir = (bkfiraf_ptr_t)handle;

    const xtfloatx4 * restrict pX;
          xtfloatx4 * restrict pDw;
    const xtfloatx4 * restrict pX0;
    const xtfloatx4 * restrict pX1;
    const xtfloatx4 * restrict pX2;
    const xtfloatx4 * restrict pH0;
    const xtfloatx4 * restrict pH1;
          xtfloatx4 * restrict pY;
    ae_valignx2 aX, aY, aX1, aX2, aH1;
    xtfloatx2 x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, xa;
    xtfloatx2 y0, y1, y2, y3;
    xtfloatx2 r0, r1, r2, r3;

    int M;
    int n, m;
    int ua_cnt;

    static const ALIGN(16) int16_t Sel[4] = { 0x0705, 0x0604, 0x0301, 0x0200 };
    ae_int16x4 sel;
    sel = AE_L16X4_I((ae_int16x4*)&Sel, 0);

    M = bkfir->M;
    NASSERT(bkfir && bkfir->magic == BKFIRAF_MAGIC && y && x);
    NASSERT_ALIGN((bkfir->coef), 16);
    NASSERT_ALIGN((bkfir->delayLine), 16);
    NASSERT_ALIGN((bkfir->delayLine + bkfir->delayLen), 16);
    if (N <= 0) return;

    // Setup pointers and circular delay line buffer.
    pX = (const xtfloatx4 *)x;
    pY = (      xtfloatx4 *)y;
    pDw= (      xtfloatx4 *)bkfir->delayPos;
    WUR_AE_CBEGIN0((uintptr_t)(bkfir->delayLine));
    WUR_AE_CEND0  ((uintptr_t)(bkfir->delayLine + bkfir->delayLen));

    ua_cnt = ((-(intptr_t)bkfir->delayPos) & 15) >> 2;
    ua_cnt = XT_MIN(ua_cnt, N);
    //process first 0..3 samples (until pDw is aligned to 16 bytes)
    if (ua_cnt)
    {
        xtfloatx2 r00, r10, r20, r30;
        xtfloatx2 r01, r11, r21, r31;
        xtfloatx2 r02, r12, r22, r32;
        xtfloatx2 r03, r13, r23, r33;

        __Pragma("loop_count min=1, max=3");
        for (n = 0; n < ua_cnt; n++)
        {
            xtfloat d0;
            AE_LSIP(d0, castxcc(xtfloat, pX), sizeof(float32_t));
            AE_SSXC(d0, castxcc(xtfloat, pDw), sizeof(float32_t));
        }

        pH0 = (const xtfloatx4 *)bkfir->coef;
        pX0 = pDw;
        AE_ADDCIRC_XC(castxcc(ae_int64, pX0), (8 - ua_cnt) * sizeof(float32_t));
        pX1 = pDw;
        AE_ADDCIRC_XC(castxcc(ae_int64, pX1), (11 - ua_cnt) * sizeof(float32_t));
        AE_LASX2X2POS_PC(aX, pX0);
        AE_LASX2X2POS_PC(aX1, pX1);

        CONST_SX2X2(r00, r01, 0); CONST_SX2X2(r02, r03, 0);
        CONST_SX2X2(r10, r11, 0); CONST_SX2X2(r12, r13, 0);
        CONST_SX2X2(r20, r21, 0); CONST_SX2X2(r22, r23, 0);
        CONST_SX2X2(r30, r31, 0); CONST_SX2X2(r32, r33, 0);

        /* preload data from x */
        AE_LASX2X2_IC(x0, x2, aX, pX0);
        x1 = XT_SEL32_LH_SX2(x0, x2);

        for (m = 0; m < ((M + 4) >> 3); m++)
        {
            /* load data from x */
            AE_LASX2X2_IC(x4, x6, aX, pX0);
            AE_LASX2X2_IC(x8, xa, aX, pX0);
            AE_LASX2X2_IC(x3, x5, aX1, pX1);
            AE_LASX2X2_IC(x7, x9, aX1, pX1);
            /* load data from y */
            AE_LSX2X2_IP(y0, y1, pH0, 4 * sizeof(float32_t));
            AE_LSX2X2_IP(y2, y3, pH0, 4 * sizeof(float32_t));
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
        if ((M + 4) & 4)
        {
            /* load data from x */
            AE_LASX2X2_IC(x4, x6, aX, pX0);
            AE_LASX2X2_IC(x3, x5, aX1, pX1);
            /* load data from y */
            AE_LSX2X2_IP(y0, y1, pH0, 4 * sizeof(float32_t));
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
        AE_S32_H_IP(AE_MOVINT32X2_FROMXTFLOATX2(r0), castxcc(ae_int32, pY), sizeof(float32_t));
        if (ua_cnt > 1) AE_S32_L_IP(AE_MOVINT32X2_FROMXTFLOATX2(r0), castxcc(ae_int32, pY), sizeof(float32_t));
        if (ua_cnt > 2) AE_S32_H_IP(AE_MOVINT32X2_FROMXTFLOATX2(r1), castxcc(ae_int32, pY), sizeof(float32_t));
        N -= ua_cnt;
    }

    aX = AE_LA128_PP(pX);
    aY = AE_ZALIGN128();
    for (n = 0; n < (N >> 3); n++)
    {
        xtfloatx2 r00, r10, r20, r30, r40, r50, r60, r70;
        xtfloatx2 r01, r11, r21, r31, r41, r51, r61, r71;

        AE_LASX2X2_IP(x0, x1, aX, pX);
        AE_LASX2X2_IP(x2, x3, aX, pX);
        AE_SSX2X2_XC(x0, x1, pDw, 4 * sizeof(float32_t));
        AE_SSX2X2_XC(x2, x3, pDw, 4 * sizeof(float32_t));

        pH0 = (const xtfloatx4 *)bkfir->coef;
        pH1 = (const xtfloatx4 *)XT_ADDX4(3, (uintptr_t)pH0);
        aH1 = AE_LA128_PP(pH1);

        /* preload data from x */
        pX0 = pDw;
        AE_LSX2XC(x0, castxcc(xtfloatx2, pX0), 2 * sizeof(float32_t));
        pX1 = pX0;
        pX2 = pX0;
        AE_ADDCIRC_XC(castxcc(ae_int64, pX2), 4 * sizeof(float32_t));
        AE_LASX2X2POS_PC(aX1, pX1);
        AE_LASX2X2POS_PC(aX2, pX2);
        AE_ADDCIRC_XC(castxcc(ae_int64, pX0), 2 * sizeof(float32_t));
        AE_LASX2X2_IC(x1, x2, aX1, pX1);

        /* load data from x */
        AE_LASX2X2_IC(x3, x4, aX2, pX2);
        /* load data from y */
        AE_LSX2X2_IP(y2, y3, pH0, 4 * sizeof(float32_t));
        y1 = XT_SEL32_LH_SX2(y2, y3);
        /* compute correlation of 8 values */
        CONST_SX2X2(r10, r30, 0);
        CONST_SX2X2(r50, r70, 0);
        MUL_SX2X2(r11, r31, x1, x2, y1, y1);
        MUL_SX2X2(r51, r71, x3, x4, y1, y1);
        MUL_SX2X2(r00, r01, x0, x1, y2, y3);
        MUL_SX2X2(r20, r21, x1, x2, y2, y3);
        MUL_SX2X2(r40, r41, x2, x3, y2, y3);
        MUL_SX2X2(r60, r61, x3, x4, y2, y3);
        /* reload data for the next iteration */
        AE_LSX2XC(x0, castxcc(xtfloatx2, pX0), 4 * sizeof(float32_t));
        AE_LASX2X2_IC(x1, x2, aX1, pX1);

        __Pragma("loop_count min=1");
        for (m = 0; m < (M >> 2); m++)
        {
            /* load data from x */
            AE_LASX2X2_IC(x3, x4, aX2, pX2);
            /* load data from y */
            AE_LASX2X2_IP(y0, y1, aH1, pH1);
            AE_LSX2X2_IP(y2, y3, pH0, 4 * sizeof(float32_t));
            /* compute correlation of 8 values */
            MADD_SX2X2(r10, r11, x0, x1, y0, y1);
            MADD_SX2X2(r30, r31, x1, x2, y0, y1);
            MADD_SX2X2(r50, r51, x2, x3, y0, y1);
            MADD_SX2X2(r70, r71, x3, x4, y0, y1);
            MADD_SX2X2(r00, r01, x0, x1, y2, y3);
            MADD_SX2X2(r20, r21, x1, x2, y2, y3);
            MADD_SX2X2(r40, r41, x2, x3, y2, y3);
            MADD_SX2X2(r60, r61, x3, x4, y2, y3);
            /* reload data for the next iteration */
            AE_LSX2XC(x0, castxcc(xtfloatx2, pX0), 4 * sizeof(float32_t));
            AE_LASX2X2_IC(x1, x2, aX1, pX1);
        }
        /* reduction add */
        ADD_SX2X2(r00, r10, r00, r10, r01, r11);
        ADD_SX2X2(r20, r30, r20, r30, r21, r31);
        AE_DSELSX2(r0, r1, r00, r10, sel);
        AE_DSELSX2(r2, r3, r20, r30, sel);
        ADD_SX2X2(r0, r1, r0, r2, r1, r3);
        /* save computed samples */
        AE_SASX2X2_IP(r0, r1, aY, pY);

        /* reduction add */
        ADD_SX2X2(r40, r50, r40, r50, r41, r51);
        ADD_SX2X2(r60, r70, r60, r70, r61, r71);
        AE_DSELSX2(r0, r1, r40, r50, sel);
        AE_DSELSX2(r2, r3, r60, r70, sel);
        ADD_SX2X2(r0, r1, r0, r2, r1, r3);
        /* save computed samples */
        AE_SASX2X2_IP(r0, r1, aY, pY);
    }
#if SMALLER_CODESIZE
    AE_SA128POS_FP(aY, pY);
    {
        const xtfloat * pH;
              xtfloat * pD;
        for (n = (8 * n); n < N; n++)
        {
            xtfloat acc0;
            xtfloat x0, h0;

            XT_LSIP(x0, castxcc(xtfloat, pX), sizeof(float32_t));
            XT_SSXC(x0, castxcc(xtfloat, pDw), sizeof(float32_t));
            pH = (xtfloat *)bkfir->coef + 1;
            pD = (xtfloat *)pDw;
            AE_ADDCIRC_XC(castxcc(ae_int64, pD), 8 * sizeof(float32_t));
            acc0 = XT_CONST_S(0);
            __Pragma("concurrent");
            __Pragma("no_unroll");
            __Pragma("loop_count min=4, factor=4");
            for (m = 0; m < M; m++)
            {
                XT_LSXC(x0, pD, sizeof(float32_t));
                XT_LSIP(h0, pH, sizeof(float32_t));
                XT_MADD_S(acc0, x0, h0);
            }
            XT_SSIP(acc0, castxcc(xtfloat, pY), 4);
        }
    }
#else
    if (N & 4)
    {
        xtfloatx2 r00, r10, r20, r30;
        xtfloatx2 r01, r11, r21, r31;
        xtfloatx2 r02, r12, r22, r32;
        xtfloatx2 r03, r13, r23, r33;

        AE_LASX2X2_IP(x0, x1, aX, pX);
        AE_SSX2X2_XC(x0, x1, pDw, 4 * sizeof(float32_t));

        pH0 = (const xtfloatx4 *)bkfir->coef;
        pX0 = pDw;
        AE_ADDCIRC_XC(castxcc(ae_int64, pX0), 4 * sizeof(float32_t));
        pX1 = (const xtfloatx4 *)XT_ADDX4(3, (uintptr_t)pX0);
        AE_LASX2X2POS_PC(aX1, pX1);

        CONST_SX2X2(r00, r01, 0); CONST_SX2X2(r02, r03, 0);
        CONST_SX2X2(r10, r11, 0); CONST_SX2X2(r12, r13, 0);
        CONST_SX2X2(r20, r21, 0); CONST_SX2X2(r22, r23, 0);
        CONST_SX2X2(r30, r31, 0); CONST_SX2X2(r32, r33, 0);

        /* preload data from x */
        AE_LSX2X2_XC(x0, x2, pX0, 4 * sizeof(float32_t));
        x1 = XT_SEL32_LH_SX2(x0, x2);

        for (m = 0; m < ((M + 4) >> 3); m++)
        {
            /* load data from x */
            AE_LSX2X2_XC(x4, x6, pX0, 4 * sizeof(float32_t));
            AE_LSX2X2_XC(x8, xa, pX0, 4 * sizeof(float32_t));
            AE_LASX2X2_IC(x3, x5, aX1, pX1);
            AE_LASX2X2_IC(x7, x9, aX1, pX1);
            /* load data from y */
            AE_LSX2X2_IP(y0, y1, pH0, 4 * sizeof(float32_t));
            AE_LSX2X2_IP(y2, y3, pH0, 4 * sizeof(float32_t));
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
        if ((M + 4) & 4)
        {
            /* load data from x */
            AE_LSX2X2_XC(x4, x6, pX0, 4 * sizeof(float32_t));
            AE_LASX2X2_IC(x3, x5, aX1, pX1);
            /* load data from y */
            AE_LSX2X2_IP(y0, y1, pH0, 4 * sizeof(float32_t));
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
        AE_SASX2X2_IP(r0, r1, aY, pY);
    }
    AE_SA128POS_FP(aY, pY);

    N &= 3;
    if (N)
    {
        xtfloatx2 r00, r10, r20, r30;
        xtfloatx2 r01, r11, r21, r31;
        xtfloatx2 r02, r12, r22, r32;
        xtfloatx2 r03, r13, r23, r33;

        __Pragma("loop_count min=1, max=3");
        for (n = 0; n < N; n++)
        {
            xtfloat d0;
            AE_LSIP(d0, castxcc(xtfloat, pX), sizeof(float32_t));
            AE_SSXC(d0, castxcc(xtfloat, pDw), sizeof(float32_t));
        }

        pH0 = (const xtfloatx4 *)bkfir->coef;
        pX0 = pDw;
        AE_ADDCIRC_XC(castxcc(ae_int64, pX0), (8 - N) * sizeof(float32_t));
        pX1 = (const xtfloatx4 *)XT_ADDX4(3, (uintptr_t)pX0);
        AE_LASX2X2POS_PC(aX1, pX1);

        CONST_SX2X2(r00, r01, 0); CONST_SX2X2(r02, r03, 0);
        CONST_SX2X2(r10, r11, 0); CONST_SX2X2(r12, r13, 0);
        CONST_SX2X2(r20, r21, 0); CONST_SX2X2(r22, r23, 0);
        CONST_SX2X2(r30, r31, 0); CONST_SX2X2(r32, r33, 0);

        /* preload data from x */
        AE_LSX2X2_XC(x0, x2, pX0, 4 * sizeof(float32_t));
        x1 = XT_SEL32_LH_SX2(x0, x2);

        for (m = 0; m < ((M + 4) >> 3); m++)
        {
            /* load data from x */
            AE_LSX2X2_XC(x4, x6, pX0, 4 * sizeof(float32_t));
            AE_LSX2X2_XC(x8, xa, pX0, 4 * sizeof(float32_t));
            AE_LASX2X2_IC(x3, x5, aX1, pX1);
            AE_LASX2X2_IC(x7, x9, aX1, pX1);
            /* load data from y */
            AE_LSX2X2_IP(y0, y1, pH0, 4 * sizeof(float32_t));
            AE_LSX2X2_IP(y2, y3, pH0, 4 * sizeof(float32_t));
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
        if ((M + 4) & 4)
        {
            /* load data from x */
            AE_LSX2X2_XC(x4, x6, pX0, 4 * sizeof(float32_t));
            AE_LASX2X2_IC(x3, x5, aX1, pX1);
            /* load data from y */
            AE_LSX2X2_IP(y0, y1, pH0, 4 * sizeof(float32_t));
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
        AE_S32_H_IP(AE_MOVINT32X2_FROMXTFLOATX2(r0), castxcc(ae_int32, pY), sizeof(float32_t));
        if (N > 1) AE_S32_L_IP(AE_MOVINT32X2_FROMXTFLOATX2(r0), castxcc(ae_int32, pY), sizeof(float32_t));
        if (N > 2) AE_S32_H_IP(AE_MOVINT32X2_FROMXTFLOATX2(r1), castxcc(ae_int32, pY), sizeof(float32_t));
    }
#endif
    bkfir->delayPos = (float32_t*)pDw;
}/* bkfiraf_process() */
#else
// for SFPU
void bkfiraf_process( bkfiraf_handle_t _bkfir,
                         float32_t * restrict  y,
                   const float32_t * restrict  x, int N )
{
    bkfiraf_t *state;
    int n, m, M;
    xtfloat h0;
    xtfloat x0, x1, x2, x3;
    xtfloat s0, s1, s2, s3;
    xtfloat acc0, acc1, acc2, acc3;
    const xtfloat* restrict pX = (const xtfloat*)x;
    xtfloat* restrict pD;
    const xtfloat*   restrict pH;
    xtfloat*          pZ = (      xtfloat*)y;
    NASSERT(_bkfir);
    state = (bkfiraf_t*)_bkfir;
    NASSERT(state->coef);
    NASSERT(state->delayLine);
    NASSERT(state->delayPos);
    NASSERT_ALIGN(state->coef,8);
    NASSERT_ALIGN(state->delayLine,8);
    NASSERT((state->M%4)==0);
    NASSERT(x);
    NASSERT(y);
    if(N<=0) return;
    pD = ( xtfloat*)state->delayPos;
    pH = (const xtfloat*  )state->coef;

    M=state->M;
    NASSERT(N>0);
    NASSERT(M>0);
    WUR_AE_CBEGIN0((uintptr_t)(state->delayLine));
    WUR_AE_CEND0((uintptr_t)(state->delayLine + M));
    for (n = 0; n<(N>>2); n ++)
    {
        XT_LSIP(x0, castxcc(xtfloat, pX), 4);
        XT_LSIP(x1, castxcc(xtfloat, pX), 4);
        XT_LSIP(x2, castxcc(xtfloat, pX), 4);
        XT_LSIP(x3, castxcc(xtfloat, pX), 4);
        acc0 = XT_CONST_S(0);
        acc1 = XT_CONST_S(0);
        acc2 = XT_CONST_S(0);
        acc3 = XT_CONST_S(0);
        s0 = x0;
        s1 = x1;
        s2 = x2;
        s3 = x3;
        { xtfloat dummy;  XT_LSXC(dummy, pD, -4); }
        for (m = 0; m<M; m ++)
        {
            h0 = pH[M-m];
            XT_MADD_S(acc0, x0, h0);
            XT_MADD_S(acc1, x1, h0);
            XT_MADD_S(acc2, x2, h0);
            XT_MADD_S(acc3, x3, h0);
            x3 = x2;
            x2 = x1;
            x1 = x0;
            XT_LSXC(x0, pD, -4);
        }
        { xtfloat dummy;  XT_LSXC(dummy, pD, 4); }
        XT_SSXC(s0, pD, 4);
        XT_SSXC(s1, pD, 4);
        XT_SSXC(s2, pD, 4);
        XT_SSXC(s3, pD, 4);

        XT_SSIP(acc0, pZ, 4);
        XT_SSIP(acc1, pZ, 4);
        XT_SSIP(acc2, pZ, 4);
        XT_SSIP(acc3, pZ, 4);
    }
    for (n = (4*n); n<N; n++)
    {
        XT_LSIP(x0, castxcc(xtfloat, pX), 4);
        acc0 = XT_CONST_S(0);
        s0 = x0;
        { xtfloat dummy;  XT_LSXC(dummy, pD, -4); }
        for (m = 0; m<M; m ++)
        {
            h0 = pH[M-m];
            XT_MADD_S(acc0, x0, h0);
            XT_LSXC(x0, pD, -4);
        }
        { xtfloat dummy;  XT_LSXC(dummy, pD, 4); }
        XT_SSXC(s0, pD, 4);
        XT_SSIP(acc0, pZ, 4);
    }
    state->delayPos = (float32_t*)pD;
}
#endif
