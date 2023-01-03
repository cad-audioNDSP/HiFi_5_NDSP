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
DISCARD_FUN(void,mtx_mpytf,( void* pScr, float32_t * z, const float32_t * x,  const float32_t * yt, int M, int N, int P ))
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

#define DSEL_F32X2(out0, out1, in0, in1, sel) \
{ \
    ae_int16x4 _out0, _out1; \
    AE_DSEL16X4(_out0, _out1, AE_MOVINT16X4_FROMINT32X2(AE_MOVINT32X2_FROMXTFLOATX2(in0)), \
                              AE_MOVINT16X4_FROMINT32X2(AE_MOVINT32X2_FROMXTFLOATX2(in1)), sel); \
    out0 = AE_MOVXTFLOATX2_FROMINT32X2(AE_MOVINT32X2_FROMINT16X4(_out0)); \
    out1 = AE_MOVXTFLOATX2_FROMINT32X2(AE_MOVINT32X2_FROMINT16X4(_out1)); \
} 

void mtx_mpytf(void* pScr, float32_t * z, const float32_t * x, const float32_t * yt, int M, int N, int P)
{
    const float32_t *restrict px;
    const xtfloatx4 *restrict pX0;
    const xtfloatx4 *restrict pX1;
    const xtfloatx4 *restrict pX2;
    const xtfloatx4 *restrict pX3;
    const xtfloatx4 *restrict pY0;
    const xtfloatx4 *restrict pY1;
    const xtfloatx4 *restrict pY2;
    const xtfloatx4 *restrict pY3;
          xtfloatx4 *restrict pZ0;
          xtfloatx4 *restrict pZ1;
          xtfloatx4 *restrict pZ2;
          xtfloatx4 *restrict pZ3;
          xtfloatx2 *restrict pZ0_;
          xtfloatx2 *restrict pZ1_;
          xtfloatx2 *restrict pZ2_;
          xtfloatx2 *restrict pZ3_;
          xtfloatx4 *restrict pS0;


    ae_valignx2 v_Y0, v_Y1, v_Y2, v_Y3;
    ae_valignx2 v_X0, v_X1, v_X2, v_X3;
    ae_valignx2 v_Z0, v_Z1, v_Z2, v_Z3;
    xtfloatx2 acc00, acc01, acc02, acc03;
    xtfloatx2 acc10, acc11, acc12, acc13;
    xtfloatx2 acc20, acc21, acc22, acc23;
    xtfloatx2 acc30, acc31, acc32, acc33;
    xtfloatx2 x00, x01, x10, x11;
    xtfloatx2 x20, x21, x30, x31;   
    xtfloatx2 y00, y01, y10, y11;
    xtfloatx2 y20, y21, y30, y31;
    xtfloatx2 z00, z01, z10, z11;
    xtfloatx2 z20, z21, z30, z31; 
    ae_int32x2 zero_v;
    ae_int16x4 ind;
    xtbool2 bmask0, bmask1;
    int m, n, p;

    static const short ALIGN(16) dsel_ind[4] = { 1797, 1540, 769, 512 };

    NASSERT(x);
    NASSERT(yt);
    NASSERT(z);
    NASSERT((z != x) && (z != yt));

    if (M <= 0 || P <= 0) return;
    /* If N<=0 then clear output matrix and return */
    if (N <= 0)
    {
        xtfloatx2 * pZ0;
        ae_valign v_z0;
        int MP = M*P;
        pZ0 = (xtfloatx2 *)z;
        v_z0 = AE_ZALIGN64();
        z00 = (xtfloatx2)0.0f;
        for (n = 0; n < (MP >> 1); n++)
        {
            XT_SASX2IP(z00, v_z0, pZ0);
        }
        XT_SASX2POSFP(v_z0, pZ0);
        if (MP & 1) XT_SSI(z00, (xtfloat *)pZ0, 0);
        return;
    }

    ind = AE_L16X4_I((ae_int16x4*)(dsel_ind), 0);
    bmask0 = AE_MOVBA2((N & 2) + ((N & 2) >> 1) + 2 * (int)((N & 3) == 1));  // 2*(N%4) if N%4<2, else 3 
    bmask1 = AE_MOVBA2(((int)((N & 3) == 3)) << 1);  // 2 if (N%4)=3, else 0

    zero_v = AE_MOVINT32X2_FROMXTFLOATX2(0.f);

    if (P >= 4)
    {
        __Pragma("loop_count min=1");
        for (p = 0; p < (P >> 2); p++)
        {
            pS0 = (xtfloatx4 *)pScr;
            pY0 = (const xtfloatx4 *)(yt + (p << 2)*N);
            pY1 = (const xtfloatx4 *)((float32_t *)pY0 + N);
            pY2 = (const xtfloatx4 *)((float32_t *)pY1 + N);
            pY3 = (const xtfloatx4 *)((float32_t *)pY2 + N);

            v_Y0 = AE_LA128_PP(pY0);
            v_Y1 = AE_LA128_PP(pY1);
            v_Y2 = AE_LA128_PP(pY2);
            v_Y3 = AE_LA128_PP(pY3);
            for (n = 0; n < (N >> 2); n++)
            {
                AE_LASX2X2_IP(y00, y01, v_Y0, pY0);
                AE_LASX2X2_IP(y10, y11, v_Y1, pY1);
                AE_LASX2X2_IP(y20, y21, v_Y2, pY2);
                AE_LASX2X2_IP(y30, y31, v_Y3, pY3);
                AE_SSX2X2_IP(y00, y01, pS0, 4 * SZ_F32);
                AE_SSX2X2_IP(y10, y11, pS0, 4 * SZ_F32);
                AE_SSX2X2_IP(y20, y21, pS0, 4 * SZ_F32);
                AE_SSX2X2_IP(y30, y31, pS0, 4 * SZ_F32);
            }

            AE_LASX2X2_IP(y00, y01, v_Y0, pY0);
            AE_LASX2X2_IP(y10, y11, v_Y1, pY1);
            AE_LASX2X2_IP(y20, y21, v_Y2, pY2);
            AE_LASX2X2_IP(y30, y31, v_Y3, pY3);

            MOVF_SX2(y00, zero_v, bmask0);
            MOVF_SX2(y01, zero_v, bmask1);
            MOVF_SX2(y10, zero_v, bmask0);
            MOVF_SX2(y11, zero_v, bmask1);
            MOVF_SX2(y20, zero_v, bmask0);
            MOVF_SX2(y21, zero_v, bmask1);
            MOVF_SX2(y30, zero_v, bmask0);
            MOVF_SX2(y31, zero_v, bmask1);

            AE_SSX2X2_IP(y00, y01, pS0, 4 * SZ_F32);
            AE_SSX2X2_IP(y10, y11, pS0, 4 * SZ_F32);
            AE_SSX2X2_IP(y20, y21, pS0, 4 * SZ_F32);
            AE_SSX2X2_IP(y30, y31, pS0, 4 * SZ_F32);

            px  = x;
            pZ0 = (xtfloatx4 *)(z);

            for (m = 0; m < (M >> 2); m++)
            {
                pX0 = (const xtfloatx4 *)(px);
                pX1 = (const xtfloatx4 *)XT_ADDX4(N,(uintptr_t)pX0);
                pX2 = (const xtfloatx4 *)XT_ADDX8(N,(uintptr_t)pX0);
                pX3 = (const xtfloatx4 *)XT_ADDX4(N,(uintptr_t)pX2);
                px  = (const float32_t *)XT_ADDX8(N,(uintptr_t)pX2);

                pY0 = (const xtfloatx4 *)pScr;

                pZ1 = (xtfloatx4 *)XT_ADDX4(P, (uintptr_t)pZ0);
                pZ2 = (xtfloatx4 *)XT_ADDX4(P, (uintptr_t)pZ1);
                pZ3 = (xtfloatx4 *)XT_ADDX4(P, (uintptr_t)pZ2);  

                acc00 = acc01 = acc02 = acc03 = 0.0f;
                acc10 = acc11 = acc12 = acc13 = 0.0f;
                acc20 = acc21 = acc22 = acc23 = 0.0f;
                acc30 = acc31 = acc32 = acc33 = 0.0f;

                v_X0 = AE_LA128_PP(pX0);
                v_X1 = AE_LA128_PP(pX1);
                v_X2 = AE_LA128_PP(pX2);
                v_X3 = AE_LA128_PP(pX3);

                for (n = 0; n < (N >> 2); n++)
                {
                    AE_LASX2X2_IP(x00, x01, v_X0, pX0);
                    AE_LASX2X2_IP(x10, x11, v_X1, pX1);
                    AE_LASX2X2_IP(x20, x21, v_X2, pX2);
                    AE_LASX2X2_IP(x30, x31, v_X3, pX3);

                    AE_LSX2X2_IP(y00, y01, pY0, 4 * SZ_F32);
                    AE_LSX2X2_IP(y10, y11, pY0, 4 * SZ_F32);
                    AE_LSX2X2_IP(y20, y21, pY0, 4 * SZ_F32);
                    AE_LSX2X2_IP(y30, y31, pY0, 4 * SZ_F32);

                    MADDQ_S(acc00, acc01, y00, y10, x00);
                    MADDQ_S(acc00, acc01, y01, y11, x01);
                    MADDQ_S(acc02, acc03, y20, y30, x00);
                    MADDQ_S(acc02, acc03, y21, y31, x01);

                    MADDQ_S(acc10, acc11, y00, y10, x10);
                    MADDQ_S(acc10, acc11, y01, y11, x11);
                    MADDQ_S(acc12, acc13, y20, y30, x10);
                    MADDQ_S(acc12, acc13, y21, y31, x11);

                    MADDQ_S(acc20, acc21, y00, y10, x20);
                    MADDQ_S(acc20, acc21, y01, y11, x21);
                    MADDQ_S(acc22, acc23, y20, y30, x20);
                    MADDQ_S(acc22, acc23, y21, y31, x21);

                    MADDQ_S(acc30, acc31, y00, y10, x30);
                    MADDQ_S(acc30, acc31, y01, y11, x31);
                    MADDQ_S(acc32, acc33, y20, y30, x30);
                    MADDQ_S(acc32, acc33, y21, y31, x31);
                }

                AE_LASX2X2_IP(x00, x01, v_X0, pX0);
                AE_LASX2X2_IP(x10, x11, v_X1, pX1);
                AE_LASX2X2_IP(x20, x21, v_X2, pX2);
                AE_LASX2X2_IP(x30, x31, v_X3, pX3);

                AE_LSX2X2_IP(y00, y01, pY0, 4 * SZ_F32);
                AE_LSX2X2_IP(y10, y11, pY0, 4 * SZ_F32);
                AE_LSX2X2_IP(y20, y21, pY0, 4 * SZ_F32);
                AE_LSX2X2_IP(y30, y31, pY0, 4 * SZ_F32);

                MOVF_SX2(x00, zero_v, bmask0);
                MOVF_SX2(x01, zero_v, bmask1);
                MOVF_SX2(x10, zero_v, bmask0);
                MOVF_SX2(x11, zero_v, bmask1);
                MOVF_SX2(x20, zero_v, bmask0);
                MOVF_SX2(x21, zero_v, bmask1);
                MOVF_SX2(x30, zero_v, bmask0);
                MOVF_SX2(x31, zero_v, bmask1);

                MADDQ_S(acc00, acc01, y00, y10, x00);
                MADDQ_S(acc00, acc01, y01, y11, x01);
                MADDQ_S(acc02, acc03, y20, y30, x00);
                MADDQ_S(acc02, acc03, y21, y31, x01);

                MADDQ_S(acc10, acc11, y00, y10, x10);
                MADDQ_S(acc10, acc11, y01, y11, x11);
                MADDQ_S(acc12, acc13, y20, y30, x10);
                MADDQ_S(acc12, acc13, y21, y31, x11);

                MADDQ_S(acc20, acc21, y00, y10, x20);
                MADDQ_S(acc20, acc21, y01, y11, x21);
                MADDQ_S(acc22, acc23, y20, y30, x20);
                MADDQ_S(acc22, acc23, y21, y31, x21);

                MADDQ_S(acc30, acc31, y00, y10, x30);
                MADDQ_S(acc30, acc31, y01, y11, x31);
                MADDQ_S(acc32, acc33, y20, y30, x30);
                MADDQ_S(acc32, acc33, y21, y31, x31);

#if 1
                DSEL_F32X2(acc00, acc01,acc00, acc01,ind);
                DSEL_F32X2(acc02, acc03,acc02, acc03,ind);
                DSEL_F32X2(acc10, acc11,acc10, acc11,ind);
                DSEL_F32X2(acc12, acc13,acc12, acc13,ind);
                z00 = acc00 + acc01;
                z01 = acc02 + acc03;
                z10 = acc10 + acc11;
                z11 = acc12 + acc13;

                DSEL_F32X2(acc20, acc21,acc20, acc21,ind);
                DSEL_F32X2(acc22, acc23,acc22, acc23,ind);
                DSEL_F32X2(acc30, acc31,acc30, acc31,ind);
                DSEL_F32X2(acc32, acc33,acc32, acc33,ind);
                z20 = acc20 + acc21;
                z21 = acc22 + acc23;
                z30 = acc30 + acc31;
                z31 = acc32 + acc33;
#else
                z00 = XT_SEL32_HL_SX2(acc00, acc01) + XT_SEL32_LH_SX2(acc00, acc01);
                z01 = XT_SEL32_HL_SX2(acc02, acc03) + XT_SEL32_LH_SX2(acc02, acc03);
                z10 = XT_SEL32_HL_SX2(acc10, acc11) + XT_SEL32_LH_SX2(acc10, acc11);
                z11 = XT_SEL32_HL_SX2(acc12, acc13) + XT_SEL32_LH_SX2(acc12, acc13);

                z20 = XT_SEL32_HL_SX2(acc20, acc21) + XT_SEL32_LH_SX2(acc20, acc21);
                z21 = XT_SEL32_HL_SX2(acc22, acc23) + XT_SEL32_LH_SX2(acc22, acc23);
                z30 = XT_SEL32_HL_SX2(acc30, acc31) + XT_SEL32_LH_SX2(acc30, acc31);
                z31 = XT_SEL32_HL_SX2(acc32, acc33) + XT_SEL32_LH_SX2(acc32, acc33);
#endif

                v_Z0 = AE_ZALIGN128();
                v_Z1 = AE_ZALIGN128();
                v_Z2 = AE_ZALIGN128();
                v_Z3 = AE_ZALIGN128();

                AE_SASX2X2_IP(z00, z01, v_Z0, pZ0); AE_SA128POS_FP(v_Z0, pZ0);
                AE_SASX2X2_IP(z10, z11, v_Z1, pZ1); AE_SA128POS_FP(v_Z1, pZ1);
                AE_SASX2X2_IP(z20, z21, v_Z2, pZ2); AE_SA128POS_FP(v_Z2, pZ2);
                AE_SASX2X2_IP(z30, z31, v_Z3, pZ3); AE_SA128POS_FP(v_Z3, pZ3);

                pZ0 = (xtfloatx4 *) XT_ADDX4(P - 4, (uintptr_t)pZ3);
            }

            if (M & 2)
            {
                pX0 = (const xtfloatx4 *)px;
                pX1 = (const xtfloatx4 *)XT_ADDX4(N,(uintptr_t)pX0);
                px  = (const float32_t *)XT_ADDX8(N,(uintptr_t)pX0);

                pY0 = (const xtfloatx4 *)pScr;

                pZ1 = (xtfloatx4 *)XT_ADDX4(P, (uintptr_t)pZ0);

                acc00 = acc01 = acc02 = acc03 = 0.0f;
                acc10 = acc11 = acc12 = acc13 = 0.0f;
                acc20 = acc21 = acc22 = acc23 = 0.0f;
                acc30 = acc31 = acc32 = acc33 = 0.0f;

                v_X0 = AE_LA128_PP(pX0);
                v_X1 = AE_LA128_PP(pX1);

                for (n = 0; n < (N >> 2); n++)
                {
                    AE_LASX2X2_IP(x00, x01, v_X0, pX0);
                    AE_LASX2X2_IP(x10, x11, v_X1, pX1);

                    AE_LSX2X2_IP(y00, y01, pY0, 4 * SZ_F32);
                    AE_LSX2X2_IP(y10, y11, pY0, 4 * SZ_F32);
                    AE_LSX2X2_IP(y20, y21, pY0, 4 * SZ_F32);
                    AE_LSX2X2_IP(y30, y31, pY0, 4 * SZ_F32);

                    MADDQ_S(acc00, acc01, y00, y10, x00);
                    MADDQ_S(acc20, acc21, y01, y11, x01);
                    MADDQ_S(acc02, acc03, y20, y30, x00);
                    MADDQ_S(acc22, acc23, y21, y31, x01);

                    MADDQ_S(acc10, acc11, y00, y10, x10);
                    MADDQ_S(acc30, acc31, y01, y11, x11);
                    MADDQ_S(acc12, acc13, y20, y30, x10);
                    MADDQ_S(acc32, acc33, y21, y31, x11);
                }

                AE_LASX2X2_IP(x00, x01, v_X0, pX0);
                AE_LASX2X2_IP(x10, x11, v_X1, pX1);

                AE_LSX2X2_IP(y00, y01, pY0, 4 * SZ_F32);
                AE_LSX2X2_IP(y10, y11, pY0, 4 * SZ_F32);
                AE_LSX2X2_IP(y20, y21, pY0, 4 * SZ_F32);
                AE_LSX2X2_IP(y30, y31, pY0, 4 * SZ_F32);

                MOVF_SX2(x00, zero_v, bmask0);
                MOVF_SX2(x01, zero_v, bmask1);
                MOVF_SX2(x10, zero_v, bmask0);
                MOVF_SX2(x11, zero_v, bmask1);

                MADDQ_S(acc00, acc01, y00, y10, x00);
                MADDQ_S(acc20, acc21, y01, y11, x01);
                MADDQ_S(acc02, acc03, y20, y30, x00);
                MADDQ_S(acc22, acc23, y21, y31, x01);

                MADDQ_S(acc10, acc11, y00, y10, x10);
                MADDQ_S(acc30, acc31, y01, y11, x11);
                MADDQ_S(acc12, acc13, y20, y30, x10);
                MADDQ_S(acc32, acc33, y21, y31, x11);

                acc00 += acc20; acc01 += acc21;
                acc02 += acc22; acc03 += acc23;
                acc10 += acc30; acc11 += acc31;
                acc12 += acc32; acc13 += acc33;
#if 1
                DSEL_F32X2(acc00, acc01,acc00, acc01,ind);
                DSEL_F32X2(acc02, acc03,acc02, acc03,ind);
                DSEL_F32X2(acc10, acc11,acc10, acc11,ind);
                DSEL_F32X2(acc12, acc13,acc12, acc13,ind);

                z00 = acc00 + acc01;
                z01 = acc02 + acc03;
                z10 = acc10 + acc11;
                z11 = acc12 + acc13;
#else
                z00 = XT_SEL32_HL_SX2(acc00, acc01) + XT_SEL32_LH_SX2(acc00, acc01);
                z01 = XT_SEL32_HL_SX2(acc02, acc03) + XT_SEL32_LH_SX2(acc02, acc03);
                z10 = XT_SEL32_HL_SX2(acc10, acc11) + XT_SEL32_LH_SX2(acc10, acc11);
                z11 = XT_SEL32_HL_SX2(acc12, acc13) + XT_SEL32_LH_SX2(acc12, acc13);
#endif

                v_Z0 = AE_ZALIGN128();
                v_Z1 = AE_ZALIGN128();

                AE_SASX2X2_IP(z00, z01, v_Z0, pZ0); AE_SA128POS_FP(v_Z0, pZ0);
                AE_SASX2X2_IP(z10, z11, v_Z1, pZ1); AE_SA128POS_FP(v_Z1, pZ1);

                pZ0 = (xtfloatx4 *) XT_ADDX4(P - 4, (uintptr_t)pZ1);
            }

            if (M & 1)
            {
                pX0 = (const xtfloatx4 *)px;
                px  = (const float32_t *)XT_ADDX8(N,(uintptr_t)pX0);

                pY0 = (const xtfloatx4 *)pScr;

                acc00 = acc01 = acc10 = acc11 = 0.0f;
                acc20 = acc21 = acc30 = acc31 = 0.0f;

                v_X0 = AE_LA128_PP(pX0);

                for (n = 0; n < (N >> 2); n++)
                {
                    AE_LASX2X2_IP(x00, x01, v_X0, pX0);

                    AE_LSX2X2_IP(y00, y01, pY0, 4*SZ_F32);
                    AE_LSX2X2_IP(y10, y11, pY0, 4*SZ_F32);
                    AE_LSX2X2_IP(y20, y21, pY0, 4*SZ_F32);
                    AE_LSX2X2_IP(y30, y31, pY0, 4*SZ_F32);

                    MADD_SX2X2(acc00, acc01, y00, y01, x00, x01);
                    MADD_SX2X2(acc10, acc11, y10, y11, x00, x01);
                    MADD_SX2X2(acc20, acc21, y20, y21, x00, x01);
                    MADD_SX2X2(acc30, acc31, y30, y31, x00, x01);
                }

                AE_LASX2X2_IP(x00, x01, v_X0, pX0);

                AE_LSX2X2_IP(y00, y01, pY0, 4*SZ_F32);
                AE_LSX2X2_IP(y10, y11, pY0, 4*SZ_F32);
                AE_LSX2X2_IP(y20, y21, pY0, 4*SZ_F32);
                AE_LSX2X2_IP(y30, y31, pY0, 4*SZ_F32);

                MOVF_SX2(x00, zero_v, bmask0);
                MOVF_SX2(x01, zero_v, bmask1);

                MADD_SX2X2(acc00, acc01, y00, y01, x00, x01);
                MADD_SX2X2(acc10, acc11, y10, y11, x00, x01);
                MADD_SX2X2(acc20, acc21, y20, y21, x00, x01);
                MADD_SX2X2(acc30, acc31, y30, y31, x00, x01);

                acc00 += acc01; acc10 += acc11;
                acc20 += acc21; acc30 += acc31;
#if 1
                DSEL_F32X2(acc00, acc10,acc00, acc10,ind);
                DSEL_F32X2(acc20, acc30,acc20, acc30,ind);
                z00 = acc00 + acc10;
                z01 = acc20 + acc30;
#else
                z00 = XT_SEL32_HL_SX2(acc00, acc10) + XT_SEL32_LH_SX2(acc00, acc10);
                z01 = XT_SEL32_HL_SX2(acc20, acc30) + XT_SEL32_LH_SX2(acc20, acc30);
#endif
                v_Z0 = AE_ZALIGN128();
                AE_SASX2X2_IP(z00, z01, v_Z0, pZ0);
                AE_SA128POS_FP(v_Z0, pZ0);
            }
            z += 4;
        }
    }


    if (P & 2)
    {
        ae_valign v_Z0_, v_Z1_, v_Z2_, v_Z3_;

        pS0 = (xtfloatx4 *)pScr;
        pY0 = (const xtfloatx4 *)(yt + (P&~3)*N);
        pY1 = (const xtfloatx4 *)((float32_t *)pY0 + N);

        v_Y0 = AE_LA128_PP(pY0);
        v_Y1 = AE_LA128_PP(pY1);
        for (n = 0; n < (N >> 2); n++)
        {
            AE_LASX2X2_IP(y00, y01, v_Y0, pY0);
            AE_LASX2X2_IP(y10, y11, v_Y1, pY1);
            AE_SSX2X2_IP(y00, y01, pS0, 4 * SZ_F32);
            AE_SSX2X2_IP(y10, y11, pS0, 4 * SZ_F32);
        }

        AE_LASX2X2_IP(y00, y01, v_Y0, pY0);
        AE_LASX2X2_IP(y10, y11, v_Y1, pY1);

        MOVF_SX2(y00, zero_v, bmask0);
        MOVF_SX2(y01, zero_v, bmask1);
        MOVF_SX2(y10, zero_v, bmask0);
        MOVF_SX2(y11, zero_v, bmask1);

        AE_SSX2X2_IP(y00, y01, pS0, 4*SZ_F32);
        AE_SSX2X2_IP(y10, y11, pS0, 4*SZ_F32);

        px = x;
        pZ0_ = (xtfloatx2 *)(z);

        for (m = 0; m < (M >> 2); m++)
        {
            pX0 = (const xtfloatx4 *)(px);
            pX1 = (const xtfloatx4 *)XT_ADDX4(N,(uintptr_t)pX0);
            pX2 = (const xtfloatx4 *)XT_ADDX8(N,(uintptr_t)pX0);
            pX3 = (const xtfloatx4 *)XT_ADDX4(N,(uintptr_t)pX2);
            px  = (const float32_t *)XT_ADDX8(N,(uintptr_t)pX2);

            pY0 = (const xtfloatx4 *)pScr;
            pY1 = (const xtfloatx4 *)((const float32_t *)pScr + 4);

            pZ1_ = (xtfloatx2 *)XT_ADDX4(P, (uintptr_t)pZ0_);
            pZ2_ = (xtfloatx2 *)XT_ADDX4(P, (uintptr_t)pZ1_);
            pZ3_ = (xtfloatx2 *)XT_ADDX4(P, (uintptr_t)pZ2_);

            acc00 = acc01 = acc02 = acc03 = 0.0f;
            acc10 = acc11 = acc12 = acc13 = 0.0f;
            acc20 = acc21 = acc22 = acc23 = 0.0f;
            acc30 = acc31 = acc32 = acc33 = 0.0f;

            v_X0 = AE_LA128_PP(pX0);
            v_X1 = AE_LA128_PP(pX1);
            v_X2 = AE_LA128_PP(pX2);
            v_X3 = AE_LA128_PP(pX3);

            //__Pragma("no_unroll");
            __Pragma("ymemory(pY0)")
            __Pragma("ymemory(pY1)")
            for (n = 0; n < (N >> 2); n++)
            {
                AE_LSX2X2_IP(y00, y01, pY0, 2*4*SZ_F32);
                AE_LSX2X2_IP(y10, y11, pY1, 2*4*SZ_F32);

                AE_LASX2X2_IP(x00, x01, v_X0, pX0);
                AE_LASX2X2_IP(x10, x11, v_X1, pX1);
                AE_LASX2X2_IP(x20, x21, v_X2, pX2);
                AE_LASX2X2_IP(x30, x31, v_X3, pX3);

                MADDQ_S(acc00, acc01, y00, y10, x00);
                MADDQ_S(acc02, acc03, y01, y11, x01);

                MADDQ_S(acc10, acc11, y00, y10, x10);
                MADDQ_S(acc12, acc13, y01, y11, x11);

                MADDQ_S(acc20, acc21, y00, y10, x20);
                MADDQ_S(acc22, acc23, y01, y11, x21);

                MADDQ_S(acc30, acc31, y00, y10, x30);
                MADDQ_S(acc32, acc33, y01, y11, x31);
            }

            AE_LSX2X2_IP(y00, y01, pY0, 2*4*SZ_F32);
            AE_LSX2X2_IP(y10, y11, pY1, 2*4*SZ_F32);

            AE_LASX2X2_IP(x00, x01, v_X0, pX0);
            AE_LASX2X2_IP(x10, x11, v_X1, pX1);
            AE_LASX2X2_IP(x20, x21, v_X2, pX2);
            AE_LASX2X2_IP(x30, x31, v_X3, pX3);
            
            MOVF_SX2(x00, zero_v, bmask0);
            MOVF_SX2(x01, zero_v, bmask1);
            MOVF_SX2(x10, zero_v, bmask0);
            MOVF_SX2(x11, zero_v, bmask1);
            MOVF_SX2(x20, zero_v, bmask0);
            MOVF_SX2(x21, zero_v, bmask1);
            MOVF_SX2(x30, zero_v, bmask0);
            MOVF_SX2(x31, zero_v, bmask1);

            MADDQ_S(acc00, acc01, y00, y10, x00);
            MADDQ_S(acc02, acc03, y01, y11, x01);

            MADDQ_S(acc10, acc11, y00, y10, x10);
            MADDQ_S(acc12, acc13, y01, y11, x11);

            MADDQ_S(acc20, acc21, y00, y10, x20);
            MADDQ_S(acc22, acc23, y01, y11, x21);

            MADDQ_S(acc30, acc31, y00, y10, x30);
            MADDQ_S(acc32, acc33, y01, y11, x31);

            acc00 += acc02; acc01 += acc03;
            acc10 += acc12; acc11 += acc13;
            acc20 += acc22; acc21 += acc23;
            acc30 += acc32; acc31 += acc33;
#if 1
            DSEL_F32X2(acc00, acc01,acc00, acc01,ind);
            DSEL_F32X2(acc10, acc11,acc10, acc11,ind);

            DSEL_F32X2(acc20, acc21,acc20, acc21,ind);
            DSEL_F32X2(acc30, acc31,acc30, acc31,ind);
            z00 = acc00 + acc01;
            z10 = acc10 + acc11;
            
            z20 = acc20 + acc21;
            z30 = acc30 + acc31;
#else
            z00 = XT_SEL32_HL_SX2(acc00, acc01) + XT_SEL32_LH_SX2(acc00, acc01);
            z10 = XT_SEL32_HL_SX2(acc10, acc11) + XT_SEL32_LH_SX2(acc10, acc11);

            z20 = XT_SEL32_HL_SX2(acc20, acc21) + XT_SEL32_LH_SX2(acc20, acc21);
            z30 = XT_SEL32_HL_SX2(acc30, acc31) + XT_SEL32_LH_SX2(acc30, acc31);
#endif

            v_Z0_ = AE_ZALIGN64();
            v_Z1_ = AE_ZALIGN64();
            v_Z2_ = AE_ZALIGN64();
            v_Z3_ = AE_ZALIGN64();
            AE_SASX2IP(z00, v_Z0_, pZ0_); AE_SA64POS_FP(v_Z0_, pZ0_);
            AE_SASX2IP(z10, v_Z1_, pZ1_); AE_SA64POS_FP(v_Z1_, pZ1_);
            AE_SASX2IP(z20, v_Z2_, pZ2_); AE_SA64POS_FP(v_Z2_, pZ2_);
            AE_SASX2IP(z30, v_Z3_, pZ3_); AE_SA64POS_FP(v_Z3_, pZ3_);

            pZ0_ = (xtfloatx2 *) XT_ADDX4(P - 2, (uintptr_t)pZ3_);
        }

        if (M & 2)
        {
            pX0 = (const xtfloatx4 *)(px);
            pX1 = (const xtfloatx4 *)XT_ADDX4(N,(uintptr_t)pX0);
            px  = (const float32_t *)XT_ADDX8(N,(uintptr_t)pX0);

            pY0 = (const xtfloatx4 *)pScr;

            pZ1_ = (xtfloatx2 *)XT_ADDX4(P, (uintptr_t)pZ0_);
            pZ2_ = (xtfloatx2 *)XT_ADDX4(P, (uintptr_t)pZ1_);
            pZ3_ = (xtfloatx2 *)XT_ADDX4(P, (uintptr_t)pZ2_);  

            acc00 = acc01 = acc02 = acc03 = 0.0f;
            acc10 = acc11 = acc12 = acc13 = 0.0f;

            v_X0 = AE_LA128_PP(pX0);
            v_X1 = AE_LA128_PP(pX1);

            //__Pragma("no_unroll");
            for (n = 0; n < (N >> 2); n++)
            {
                AE_LSX2X2_IP(y00, y01, pY0, 4*SZ_F32);
                AE_LSX2X2_IP(y10, y11, pY0, 4*SZ_F32);

                AE_LASX2X2_IP(x00, x01, v_X0, pX0);
                AE_LASX2X2_IP(x10, x11, v_X1, pX1);

                MADDQ_S(acc00, acc01, y00, y10, x00);
                MADDQ_S(acc02, acc03, y01, y11, x01);

                MADDQ_S(acc10, acc11, y00, y10, x10);
                MADDQ_S(acc12, acc13, y01, y11, x11);
            }

            AE_LSX2X2_IP(y00, y01, pY0, 4*SZ_F32);
            AE_LSX2X2_IP(y10, y11, pY0, 4*SZ_F32);

            AE_LASX2X2_IP(x00, x01, v_X0, pX0);
            AE_LASX2X2_IP(x10, x11, v_X1, pX1);
            
            MOVF_SX2(x00, zero_v, bmask0);
            MOVF_SX2(x01, zero_v, bmask1);
            MOVF_SX2(x10, zero_v, bmask0);
            MOVF_SX2(x11, zero_v, bmask1);

            MADDQ_S(acc00, acc01, y00, y10, x00);
            MADDQ_S(acc02, acc03, y01, y11, x01);

            MADDQ_S(acc10, acc11, y00, y10, x10);
            MADDQ_S(acc12, acc13, y01, y11, x11);

            acc00 += acc02; acc01 += acc03;
            acc10 += acc12; acc11 += acc13;

#if 1
            DSEL_F32X2(acc00, acc01,acc00, acc01,ind);
            DSEL_F32X2(acc10, acc11,acc10, acc11,ind);
            z00 = acc00 + acc01;
            z10 = acc10 + acc11;
#else
            z00 = XT_SEL32_HL_SX2(acc00, acc01) + XT_SEL32_LH_SX2(acc00, acc01);
            z10 = XT_SEL32_HL_SX2(acc10, acc11) + XT_SEL32_LH_SX2(acc10, acc11);
#endif

            v_Z0_ = AE_ZALIGN64();
            v_Z1_ = AE_ZALIGN64();
            v_Z2_ = AE_ZALIGN64();
            v_Z3_ = AE_ZALIGN64();
            AE_SASX2IP(z00, v_Z0_, pZ0_); AE_SA64POS_FP(v_Z0_, pZ0_);
            AE_SASX2IP(z10, v_Z1_, pZ1_); AE_SA64POS_FP(v_Z1_, pZ1_);

            pZ0_ = (xtfloatx2 *) XT_ADDX4(P - 2, (uintptr_t)pZ1_);
        }

        if ( M&1 )
        {
            pX0 = (const xtfloatx4 *)(px);
            pX1 = (const xtfloatx4 *)XT_ADDX4(N,(uintptr_t)pX0);

            pY0 = (const xtfloatx4 *)pScr;

            pZ1_ = (xtfloatx2 *)XT_ADDX4(P, (uintptr_t)pZ0_);
            pZ2_ = (xtfloatx2 *)XT_ADDX4(P, (uintptr_t)pZ1_);
            pZ3_ = (xtfloatx2 *)XT_ADDX4(P, (uintptr_t)pZ2_);  

            v_X0 = AE_LA128_PP(pX0);

            acc00 = acc01 = acc10 = acc11 = 0.0f;

            for (n = 0; n < (N >> 2); n++)
            {
                AE_LASX2X2_IP(x00, x01, v_X0, pX0);

                AE_LSX2X2_IP(y00, y01, pY0, 4*SZ_F32);
                AE_LSX2X2_IP(y10, y11, pY0, 4*SZ_F32);

                MADD_SX2X2(acc00, acc01, y00, y01, x00, x01);
                MADD_SX2X2(acc10, acc11, y10, y11, x00, x01);
            }

            AE_LASX2X2_IP(x00, x01, v_X0, pX0);

            AE_LSX2X2_IP(y00, y01, pY0, 4*SZ_F32);
            AE_LSX2X2_IP(y10, y11, pY0, 4*SZ_F32);

            MOVF_SX2(x00, zero_v, bmask0);
            MOVF_SX2(x01, zero_v, bmask1);;

            MADD_SX2X2(acc00, acc01, y00, y01, x00, x01);
            MADD_SX2X2(acc10, acc11, y10, y11, x00, x01);

            acc00 += acc01; acc10 += acc11;

#if 1
            DSEL_F32X2(acc00, acc10,acc00, acc10,ind);
            z00 = acc00 + acc10;
#else
            z00 = XT_SEL32_HL_SX2(acc00, acc10) + XT_SEL32_LH_SX2(acc00, acc10);
#endif

            v_Z0_ = AE_ZALIGN64();
            AE_SASX2IP(z00, v_Z0_, pZ0_);
            AE_SA64POS_FP(v_Z0_, pZ0_);
        }
        z += 2;
    }

    if (P & 1){
        pS0 = (xtfloatx4 *)pScr;
        pY0 = (const xtfloatx4 *)(yt + (P-1)*N);

        v_Y0 = AE_LA128_PP(pY0);
        for (n = 0; n < (N >> 2); n++)
        {
            AE_LASX2X2_IP(y00, y01, v_Y0, pY0);
            AE_SSX2X2_IP(y00, y01, pS0, 4 * SZ_F32);
        }
        AE_LASX2X2_IP(y00, y01, v_Y0, pY0);

        MOVF_SX2(y00, zero_v, bmask0);
        MOVF_SX2(y01, zero_v, bmask1);

        AE_SSX2X2_IP(y00, y01, pS0, 4*SZ_F32);

        px = x;

        for (m = 0; m < (M >> 2); m++)
        {
            pX0 = (const xtfloatx4 *)(px);
            pX1 = (const xtfloatx4 *)XT_ADDX4(N,(uintptr_t)pX0);
            pX2 = (const xtfloatx4 *)XT_ADDX8(N,(uintptr_t)pX0);
            pX3 = (const xtfloatx4 *)XT_ADDX4(N,(uintptr_t)pX2);
            px  = (const float32_t *)XT_ADDX8(N,(uintptr_t)pX2);

            pY0 = (const xtfloatx4 *)pScr;

            acc00 = acc01 = 0.0f;
            acc10 = acc11 = 0.0f;
            acc20 = acc21 = 0.0f;
            acc30 = acc31 = 0.0f;

            v_X0 = AE_LA128_PP(pX0);
            v_X1 = AE_LA128_PP(pX1);
            v_X2 = AE_LA128_PP(pX2);
            v_X3 = AE_LA128_PP(pX3);

            for (n = 0; n < (N >> 2); n++)
            {
                AE_LSX2X2_IP(y00, y01, pY0, 4*SZ_F32);

                AE_LASX2X2_IP(x00, x01, v_X0, pX0);
                AE_LASX2X2_IP(x10, x11, v_X1, pX1);
                AE_LASX2X2_IP(x20, x21, v_X2, pX2);
                AE_LASX2X2_IP(x30, x31, v_X3, pX3);

                MADD_SX2X2(acc00, acc01, y00, y01, x00, x01);
                MADD_SX2X2(acc10, acc11, y00, y01, x10, x11);
                MADD_SX2X2(acc20, acc21, y00, y01, x20, x21);
                MADD_SX2X2(acc30, acc31, y00, y01, x30, x31);
            }

            AE_LSX2X2_IP(y00, y01, pY0, 4*SZ_F32);

            AE_LASX2X2_IP(x00, x01, v_X0, pX0);
            AE_LASX2X2_IP(x10, x11, v_X1, pX1);
            AE_LASX2X2_IP(x20, x21, v_X2, pX2);
            AE_LASX2X2_IP(x30, x31, v_X3, pX3);
            
            MOVF_SX2(x00, zero_v, bmask0);
            MOVF_SX2(x01, zero_v, bmask1);
            MOVF_SX2(x10, zero_v, bmask0);
            MOVF_SX2(x11, zero_v, bmask1);
            MOVF_SX2(x20, zero_v, bmask0);
            MOVF_SX2(x21, zero_v, bmask1);
            MOVF_SX2(x30, zero_v, bmask0);
            MOVF_SX2(x31, zero_v, bmask1);
            
            MADD_SX2X2(acc00, acc01, y00, y01, x00, x01);
            MADD_SX2X2(acc10, acc11, y00, y01, x10, x11);
            MADD_SX2X2(acc20, acc21, y00, y01, x20, x21);
            MADD_SX2X2(acc30, acc31, y00, y01, x30, x31);

            acc00 += acc01; 
            acc10 += acc11;
            acc20 += acc21;
            acc30 += acc31;

            AE_SSXP(RADD_SX2(acc00), z, P*SZ_F32);
            AE_SSXP(RADD_SX2(acc10), z, P*SZ_F32);
            AE_SSXP(RADD_SX2(acc20), z, P*SZ_F32);
            AE_SSXP(RADD_SX2(acc30), z, P*SZ_F32);
        }

        if (M & 2)
        {
            pX0 = (const xtfloatx4 *)(px);
            pX1 = (const xtfloatx4 *)XT_ADDX4(N,(uintptr_t)pX0);
            px  = (const float32_t *)XT_ADDX8(N,(uintptr_t)pX0);

            pY0 = (const xtfloatx4 *)pScr;

            acc00 = acc01 = 0.0f;
            acc10 = acc11 = 0.0f;
            acc20 = acc21 = 0.0f;
            acc30 = acc31 = 0.0f;

            v_X0 = AE_LA128_PP(pX0);
            v_X1 = AE_LA128_PP(pX1);

            for (n = 0; n < (N >> 2); n++)
            {
                AE_LSX2X2_IP(y00, y01, pY0, 4*SZ_F32);

                AE_LASX2X2_IP(x00, x01, v_X0, pX0);
                AE_LASX2X2_IP(x10, x11, v_X1, pX1);

                MADD_SX2X2(acc00, acc01, y00, y01, x00, x01);
                MADD_SX2X2(acc10, acc11, y00, y01, x10, x11);
            }

            AE_LSX2X2_IP(y00, y01, pY0, 4*SZ_F32);

            AE_LASX2X2_IP(x00, x01, v_X0, pX0);
            AE_LASX2X2_IP(x10, x11, v_X1, pX1);
            
            MOVF_SX2(x00, zero_v, bmask0);
            MOVF_SX2(x01, zero_v, bmask1);
            MOVF_SX2(x10, zero_v, bmask0);
            MOVF_SX2(x11, zero_v, bmask1);
           
            MADD_SX2X2(acc00, acc01, y00, y01, x00, x01);
            MADD_SX2X2(acc10, acc11, y00, y01, x10, x11);

            acc00 += acc01; 
            acc10 += acc11;

            AE_SSXP(RADD_SX2(acc00), z, P*SZ_F32);
            AE_SSXP(RADD_SX2(acc10), z, P*SZ_F32);
        }

        if (M & 1)
        {
            pX0 = (const xtfloatx4 *)px;
            pY0 = (const xtfloatx4 *)pScr;

            v_X0 = AE_LA128_PP(pX0);

            acc00 = acc01 = acc10 = acc11 = 0.0f;

            for (n = 0; n < (N >> 2); n++)
            {
                AE_LASX2X2_IP(x00, x01, v_X0, pX0);
                AE_LSX2X2_IP(y00, y01, pY0, 4*SZ_F32);
                MADD_SX2X2(acc00, acc01, y00, y01, x00, x01);
            }

            AE_LASX2X2_IP(x00, x01, v_X0, pX0);
            AE_LSX2X2_IP(y00, y01, pY0, 4*SZ_F32);

            MOVF_SX2(x00, zero_v, bmask0);
            MOVF_SX2(x01, zero_v, bmask1);

            MADD_SX2X2(acc00, acc01, y00, y01, x00, x01);
            acc00 += acc01;
            AE_SSXP(RADD_SX2(acc00), z, P*SZ_F32);
        }
    }
} /* mtx_mpytf() */
#elif (HAVE_FPU)
void mtx_mpytf( void* pScr, float32_t * z, const float32_t * x,  const float32_t * yt, int M, int N, int P )
{
    int m, n, p, P0=(P&~3);
    xtfloat * restrict pZ;
    const xtfloat * restrict pX;
    const xtfloat * restrict pY;
    const xtfloat * restrict pY0;
    const xtfloat * restrict pY1;
    const xtfloat * restrict pY2;
    const xtfloat * restrict pY3;

    NASSERT(x);
    NASSERT(yt);
    NASSERT(z);
    NASSERT((z != x) && (z != yt));
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
    if (P0>0)
    {
        for (m = 0; m<M; m++)
        {
            pZ = (xtfloat  *)(z+m*P);
            for (p = 0; p<P0; p+=4)
            {
                xtfloat A0,A1,A2,A3,x0,y0,y1,y2,y3;
                pX = (const xtfloat  *)(x+m*N);
                pY = (const xtfloat  *)(yt+p*N);
                pY0 = (const xtfloat  *)pY;
                pY1 = (const xtfloat  *)(pY0+N);
                pY2 = (const xtfloat  *)(pY1+N);
                pY3 = (const xtfloat  *)(pY2+N);
                A0=A1=A2=A3=XT_CONST_S(0);
                for ( n=0; n<N; n++ )
                {
                    XT_LSIP(x0, pX, sizeof(xtfloat));
                    /*
                    y1=XT_LSI(pY, 1*sizeof(xtfloat));
                    y2=XT_LSI(pY, 2*sizeof(xtfloat));
                    y3=XT_LSI(pY, 3*sizeof(xtfloat));
                    XT_LSXP(y0, pY, P*sizeof(xtfloat));
                    */
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
    if (P&2)
    {
        for (m = 0; m<M; m++)
        {
            xtfloat A0,A1,x0,y0,y1;
            p = (P&~3);
            pZ = (xtfloat  *)(z+m*P+p);
            pX = (const xtfloat  *)(x+m*N);
            pY = (const xtfloat  *)(yt+p*N);
            pY0 = (const xtfloat  *)pY;
            pY1 = (const xtfloat  *)(pY0+N);
            A0=A1=XT_CONST_S(0);
            for ( n=0; n<N; n++ )
            {
                XT_LSIP(x0, pX, sizeof(xtfloat));
                /*
                y1=XT_LSI(pY, 1*sizeof(xtfloat));
                XT_LSXP(y0, pY, P*sizeof(xtfloat));
                */
                XT_LSIP(y0, pY0, sizeof(xtfloat));
                XT_LSIP(y1, pY1, sizeof(xtfloat));
                XT_MADD_S(A0, x0, y0);
                XT_MADD_S(A1, x0, y1);
            }
            XT_SSIP(A0, pZ, sizeof(xtfloat));
            XT_SSIP(A1, pZ, sizeof(xtfloat));
        }
    }
    if (P&1)
    {
        for (m = 0; m<M; m++)
        {
            xtfloat A0,x0,y0;
            p = (P&~1);
            pZ = (xtfloat  *)(z+m*P+p);
            pX = (const xtfloat  *)(x+m*N);
            pY = (const xtfloat  *)(yt+p*N);
            pY0 = (const xtfloat  *)pY;
            A0=XT_CONST_S(0);
            for ( n=0; n<N; n++ )
            {
                XT_LSIP(x0, pX, sizeof(xtfloat));
                //              XT_LSXP(y0, pY, P*sizeof(xtfloat));
                XT_LSIP(y0, pY0, sizeof(xtfloat));
                XT_MADD_S(A0, x0, y0);
            }
            XT_SSIP(A0, pZ, sizeof(xtfloat));
        }
    }

} /* mtx_mpytf() */

#endif

size_t mtx_mpytf_getScratchSize(int M, int N, int P)
{
    (void)M; (void)P;
    return N <= 0 ? 0 : (((N)+4) * 4 * sizeof(float32_t));
} /* mtx_mpytf_getScratchSize() */
