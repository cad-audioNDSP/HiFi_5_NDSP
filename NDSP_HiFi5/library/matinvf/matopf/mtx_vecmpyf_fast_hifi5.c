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
#include "common_fpu.h"

#if (HAVE_FPU==0 && HAVE_VFPU==0)
DISCARD_FUN(void ,mtx_vecmpyf_fast,( float32_t * z, const float32_t * x,  const float32_t * y, int M, int N ))
#elif (HAVE_VFPU)

#define SZ_F32 (sizeof(float32_t))

/*-------------------------------------------------------------------------
  Matrix by Vector Multiply
  These functions compute the expression z = 2^lsh * x * y for the matrices 
  x and vector y. 
  NOTE: lsh factor is not relevant for floating point routines.

  Two versions of functions available: regular version (mtx_vecmpy32x32,  
  mtx_vecmpy16x16, mtx_vecmpy8x8, mtx_vecmpy8x16, mtx_vecmpyf) with arbitrary 
  arguments and faster version (mtx_vecmpy32x32_fast, mtx_vecmpy16x16_fast, 
  mtx_vecmpy8x8_fast, mtx_vecmpy8x16_fast,  mtx_vecmpyf_fast) that apply 
  some restrictions

  Precision: 
  32x32 32-bit input, 32-bit output
  16x16 16-bit input, 16-bit output
  8x8   8-bit inputs, 8-bit output
  8x16  8/16-bit inputs, 16-bit output
  f     floating point

  Input:
  x[M*N] input matrix,Q31,Q15 or floating point
  y[N]   input vector,Q31,Q15 or floating point
  M      number of rows in matrix x
  N      number of columns in matrix x
  lsh    additional left shift(applied to the fixed-
         point functions only) 
  Output:
  z[M]   output vector,Q31,Q15 or floating point

  Restriction:
  For regular routines (mtx_vecmpy32x32, mtx_vecmpy16x16, mtx_vecmpy8x8,
  mtx_vecmpy8x16,  mtx_vecmpyf)
  x,y,z should not overlap

  For faster routines  (mtx_vecmpy32x32_fast, mtx_vecmpy16x16_fast, 
  mtx_vecmpy8x8_fast, mtx_vecmpy8x16_fast, mtx_vecmpyf_fast)
  x,y,z   should not overlap
  x,y     aligned on 16-byte boundary
  N, M    multiples of 4
  lsh     should be in range:
          -31...31 for mtx_vecmpy32x32, mtx_vecmpy32x32_fast
          -15...15 for mtx_vecmpy16x16, mtx_vecmpy16x16_fast, 
                   mtx_vecmpy8x8_fast, mtx_vecmpy8x16_fast   
-------------------------------------------------------------------------*/
void mtx_vecmpyf_fast(float32_t * z, const float32_t * x, const float32_t * y, int M, int N)
{
    const xtfloatx4 *restrict px0;
    const xtfloatx4 *restrict px1;
    const xtfloatx4 *restrict py;
    xtfloatx4 *restrict pz;
    xtfloatx2 y0, y1, y2, y3;
    xtfloatx2 z0, z1, z2, z3;
    xtfloatx2 x00, x01, x02, x03, x10, x11, x12, x13;
    xtfloatx2 x20, x21, x22, x23, x30, x31, x32, x33;
    xtfloatx2 x40, x41, x42, x43, x50, x51, x52, x53;
    xtfloatx2 x60, x61, x62, x63, x70, x71, x72, x73;
    xtfloatx2 acc00, acc01, acc10, acc11;
    xtfloatx2 acc20, acc21, acc30, acc31;
    xtfloatx2 acc40, acc41, acc50, acc51;
    xtfloatx2 acc60, acc61, acc70, acc71;
    ae_valignx2 v_z;
    int m, n;

    NASSERT(x);
    NASSERT(y);
    NASSERT(z);
    NASSERT((z != x) && (z != y));
    NASSERT_ALIGN(x, HIFI_SIMD_WIDTH);
    NASSERT_ALIGN(y, HIFI_SIMD_WIDTH);
    NASSERT(M % 4 == 0);
    NASSERT(N % 4 == 0);

    if (M < 4) return;

    pz = (xtfloatx4 *)z;
    /* If N<4 then clear output veector and return */
    if (N < 4)
    {
        z0 = z1 = (xtfloatx2)0.0f;
        __Pragma("loop_count min=1")
        for (m = 0; m < (M >> 2); m++)
        {
            AE_SSX2X2_IP(z0, z1, pz, 4 * SZ_F32);
        }
        return;
    }
    
    v_z = AE_ZALIGN128();
    px0 = (const xtfloatx4 *)(x);

    if (N == 4){
        px1 = (const xtfloatx4 *)((float32_t *)px0 + N);
        for (m = 0; m < (M>>3); m++)
        {
            py  = (const xtfloatx4 *)(y);
     
            AE_LSX2X2_XP(x00, x01, px0, 2*N*SZ_F32);
            AE_LSX2X2_XP(x10, x11, px1, 2*N*SZ_F32);
            AE_LSX2X2_XP(x20, x21, px0, 2*N*SZ_F32);
            AE_LSX2X2_XP(x30, x31, px1, 2*N*SZ_F32);

            AE_LSX2X2_XP(x40, x41, px0, 2*N*SZ_F32);
            AE_LSX2X2_XP(x50, x51, px1, 2*N*SZ_F32);
            AE_LSX2X2_XP(x60, x61, px0, 2*N*SZ_F32);
            AE_LSX2X2_XP(x70, x71, px1, 2*N*SZ_F32);

            AE_LSX2X2_IP(y0, y1, py, 4*SZ_F32);

            MULQ_S(acc00, acc10, x00, x10, y0);
            MULQ_S(acc01, acc11, x01, x11, y1);
            MULQ_S(acc20, acc30, x20, x30, y0);
            MULQ_S(acc21, acc31, x21, x31, y1);

            MULQ_S(acc40, acc50, x40, x50, y0);
            MULQ_S(acc41, acc51, x41, x51, y1);
            MULQ_S(acc60, acc70, x60, x70, y0);
            MULQ_S(acc61, acc71, x61, x71, y1); 

            acc00 = acc00 + acc01;
            acc10 = acc10 + acc11;
            y0 = XT_SEL32_HL_SX2(acc00, acc10);
            y1 = XT_SEL32_LH_SX2(acc00, acc10);
            z0 = y0 + y1;

            acc20 = acc20 + acc21;
            acc30 = acc30 + acc31;
            y0 = XT_SEL32_HL_SX2(acc20, acc30);
            y1 = XT_SEL32_LH_SX2(acc20, acc30);
            z1 = y0 + y1;

            acc40 = acc40 + acc41;
            acc50 = acc50 + acc51;
            y0 = XT_SEL32_HL_SX2(acc40, acc50);
            y1 = XT_SEL32_LH_SX2(acc40, acc50);
            z2 = y0 + y1;

            acc60 = acc60 + acc61;
            acc70 = acc70 + acc71;
            y0 = XT_SEL32_HL_SX2(acc60, acc70);
            y1 = XT_SEL32_LH_SX2(acc60, acc70);
            z3 = y0 + y1;

            AE_SASX2X2_IP(z0, z1, v_z, pz);
            AE_SASX2X2_IP(z2, z3, v_z, pz);
        }
        if ((M>>2)&1) {
            py  = (const xtfloatx4 *)(y);

            AE_LSX2X2_XP(x00, x01, px0, 2*N*SZ_F32);
            AE_LSX2X2_XP(x10, x11, px1, 2*N*SZ_F32);
            AE_LSX2X2_XP(x20, x21, px0, (4-2*N)*SZ_F32);
            AE_LSX2X2_XP(x30, x31, px1, (4-2*N)*SZ_F32);

            AE_LSX2X2_IP(y0, y1, py, 4*SZ_F32);

            MULQ_S(acc00, acc10, x00, x10, y0);
            MULQ_S(acc01, acc11, x01, x11, y1);
            MULQ_S(acc20, acc30, x20, x30, y0);
            MULQ_S(acc21, acc31, x21, x31, y1);

            acc00 = acc00 + acc01;
            acc10 = acc10 + acc11;
            y0 = XT_SEL32_HL_SX2(acc00, acc10);
            y1 = XT_SEL32_LH_SX2(acc00, acc10);
            z0 = y0 + y1;

            acc20 = acc20 + acc21;
            acc30 = acc30 + acc31;
            y0 = XT_SEL32_HL_SX2(acc20, acc30);
            y1 = XT_SEL32_LH_SX2(acc20, acc30);
            z1 = y0 + y1;

            AE_SASX2X2_IP(z0, z1, v_z, pz);
        }
    }  else if ( (N>>2)&1 ) { // not a multiple of 8
        px1 = (const xtfloatx4 *)((float32_t *)px0 + N);

        for (m = 0; m < (M>>3); m++)
        {            
            py = (const xtfloatx4 *)(y);

            acc00 = acc01 = acc10 = acc11 = 0;
            acc20 = acc21 = acc30 = acc31 = 0;
            acc40 = acc41 = acc50 = acc51 = 0;
            acc60 = acc61 = acc70 = acc71 = 0;

            __Pragma("loop_count min=1")
            for (n = 0; n < (N >> 3); n++)
            {
                AE_LSX2X2_XP(x00, x01, px0, 2*N*SZ_F32);
                AE_LSX2X2_XP(x10, x11, px1, 2*N*SZ_F32);
                AE_LSX2X2_XP(x20, x21, px0, 2*N*SZ_F32);
                AE_LSX2X2_XP(x30, x31, px1, 2*N*SZ_F32);

                AE_LSX2X2_XP(x40, x41, px0, 2*N*SZ_F32);
                AE_LSX2X2_XP(x50, x51, px1, 2*N*SZ_F32);
                AE_LSX2X2_XP(x60, x61, px0, (4-6*N)*SZ_F32);
                AE_LSX2X2_XP(x70, x71, px1, (4-6*N)*SZ_F32);

                AE_LSX2X2_IP(y0, y1, py, 4*SZ_F32);

                MADDQ_S(acc00, acc10, x00, x10, y0);
                MADDQ_S(acc01, acc11, x01, x11, y1);
                MADDQ_S(acc20, acc30, x20, x30, y0);
                MADDQ_S(acc21, acc31, x21, x31, y1);

                MADDQ_S(acc40, acc50, x40, x50, y0);
                MADDQ_S(acc41, acc51, x41, x51, y1);
                MADDQ_S(acc60, acc70, x60, x70, y0);
                MADDQ_S(acc61, acc71, x61, x71, y1); 

                AE_LSX2X2_XP(x00, x01, px0, 2*N*SZ_F32);
                AE_LSX2X2_XP(x10, x11, px1, 2*N*SZ_F32);
                AE_LSX2X2_XP(x20, x21, px0, 2*N*SZ_F32);
                AE_LSX2X2_XP(x30, x31, px1, 2*N*SZ_F32);

                AE_LSX2X2_XP(x40, x41, px0, 2*N*SZ_F32);
                AE_LSX2X2_XP(x50, x51, px1, 2*N*SZ_F32);
                AE_LSX2X2_XP(x60, x61, px0, (4-6*N)*SZ_F32);
                AE_LSX2X2_XP(x70, x71, px1, (4-6*N)*SZ_F32);

                AE_LSX2X2_IP(y0, y1, py, 4*SZ_F32);

                MADDQ_S(acc00, acc10, x00, x10, y0);
                MADDQ_S(acc01, acc11, x01, x11, y1);
                MADDQ_S(acc20, acc30, x20, x30, y0);
                MADDQ_S(acc21, acc31, x21, x31, y1);

                MADDQ_S(acc40, acc50, x40, x50, y0);
                MADDQ_S(acc41, acc51, x41, x51, y1);
                MADDQ_S(acc60, acc70, x60, x70, y0);
                MADDQ_S(acc61, acc71, x61, x71, y1); 
            }

            AE_LSX2X2_XP(x00, x01, px0, 2*N*SZ_F32);
            AE_LSX2X2_XP(x10, x11, px1, 2*N*SZ_F32);
            AE_LSX2X2_XP(x20, x21, px0, 2*N*SZ_F32);
            AE_LSX2X2_XP(x30, x31, px1, 2*N*SZ_F32);

            AE_LSX2X2_XP(x40, x41, px0, 2*N*SZ_F32);
            AE_LSX2X2_XP(x50, x51, px1, 2*N*SZ_F32);
            AE_LSX2X2_XP(x60, x61, px0, (4+N)*SZ_F32);
            AE_LSX2X2_XP(x70, x71, px1, (4+N)*SZ_F32);

            AE_LSX2X2_IP(y0, y1, py, 4*SZ_F32);

            MADDQ_S(acc00, acc10, x00, x10, y0);
            MADDQ_S(acc01, acc11, x01, x11, y1);
            MADDQ_S(acc20, acc30, x20, x30, y0);
            MADDQ_S(acc21, acc31, x21, x31, y1);

            MADDQ_S(acc40, acc50, x40, x50, y0);
            MADDQ_S(acc41, acc51, x41, x51, y1);
            MADDQ_S(acc60, acc70, x60, x70, y0);
            MADDQ_S(acc61, acc71, x61, x71, y1); 

            acc00 = acc00 + acc01;
            acc10 = acc10 + acc11;
            y0 = XT_SEL32_HL_SX2(acc00, acc10);
            y1 = XT_SEL32_LH_SX2(acc00, acc10);
            z0 = y0 + y1;

            acc20 = acc20 + acc21;
            acc30 = acc30 + acc31;
            y0 = XT_SEL32_HL_SX2(acc20, acc30);
            y1 = XT_SEL32_LH_SX2(acc20, acc30);
            z1 = y0 + y1;

            acc40 = acc40 + acc41;
            acc50 = acc50 + acc51;
            y0 = XT_SEL32_HL_SX2(acc40, acc50);
            y1 = XT_SEL32_LH_SX2(acc40, acc50);
            z2 = y0 + y1;

            acc60 = acc60 + acc61;
            acc70 = acc70 + acc71;
            y0 = XT_SEL32_HL_SX2(acc60, acc70);
            y1 = XT_SEL32_LH_SX2(acc60, acc70);
            z3 = y0 + y1;

            AE_SASX2X2_IP(z0, z1, v_z, pz);
            AE_SASX2X2_IP(z2, z3, v_z, pz);
        }
        if ((M>>2)&1) {
            py  = (const xtfloatx4 *)(y);

            acc00 = acc01 = acc10 = acc11 =
            acc20 = acc21 = acc30 = acc31 = (xtfloatx2)0.0f;

            __Pragma("loop_count min=2")
            for (n = 0; n < (N >> 2); n++)
            {                
                AE_LSX2X2_XP(x00, x01, px0, 2*N*SZ_F32);
                AE_LSX2X2_XP(x10, x11, px1, 2*N*SZ_F32);
                AE_LSX2X2_XP(x20, x21, px0, (4-2*N)*SZ_F32);
                AE_LSX2X2_XP(x30, x31, px1, (4-2*N)*SZ_F32);

                AE_LSX2X2_IP(y0, y1, py, 4*SZ_F32);

                MADDQ_S(acc00, acc10, x00, x10, y0);
                MADDQ_S(acc01, acc11, x01, x11, y1);
                MADDQ_S(acc20, acc30, x20, x30, y0);
                MADDQ_S(acc21, acc31, x21, x31, y1);
            }

            acc00 = acc00 + acc01;
            acc10 = acc10 + acc11;
            y0 = XT_SEL32_HL_SX2(acc00, acc10);
            y1 = XT_SEL32_LH_SX2(acc00, acc10);
            z0 = y0 + y1;

            acc20 = acc20 + acc21;
            acc30 = acc30 + acc31;
            y0 = XT_SEL32_HL_SX2(acc20, acc30);
            y1 = XT_SEL32_LH_SX2(acc20, acc30);
            z1 = y0 + y1;

            AE_SASX2X2_IP(z0, z1, v_z, pz);
        }
    } else { // multiple of 8
        for (m = 0; m < (M>>3); m++)
        {
            px1 = (const xtfloatx4 *)((float32_t *)px0 + 4);
            py  = (const xtfloatx4 *)(y);

            acc00 = acc01 = acc10 = acc11 = 0;
            acc20 = acc21 = acc30 = acc31 = 0;
            acc40 = acc41 = acc50 = acc51 = 0;
            acc60 = acc61 = acc70 = acc71 = 0;

            __Pragma("loop_count min=1")
            for (n = 0; n < (N >> 3); n++)
            {
                AE_LSX2X2_XP(x00, x01, px0, N*SZ_F32);
                AE_LSX2X2_XP(x02, x03, px1, N*SZ_F32);

                AE_LSX2X2_XP(x10, x11, px0, N*SZ_F32);
                AE_LSX2X2_XP(x12, x13, px1, N*SZ_F32);

                AE_LSX2X2_XP(x20, x21, px0, N*SZ_F32);
                AE_LSX2X2_XP(x22, x23, px1, N*SZ_F32);

                AE_LSX2X2_XP(x30, x31, px0, N*SZ_F32);
                AE_LSX2X2_XP(x32, x33, px1, N*SZ_F32);

                AE_LSX2X2_XP(x40, x41, px0, N*SZ_F32);
                AE_LSX2X2_XP(x42, x43, px1, N*SZ_F32);

                AE_LSX2X2_XP(x50, x51, px0, N*SZ_F32);
                AE_LSX2X2_XP(x52, x53, px1, N*SZ_F32);

                AE_LSX2X2_XP(x60, x61, px0, N*SZ_F32);
                AE_LSX2X2_XP(x62, x63, px1, N*SZ_F32);

                AE_LSX2X2_XP(x70, x71, px0, (8 - 7*N)*SZ_F32);
                AE_LSX2X2_XP(x72, x73, px1, (8 - 7*N)*SZ_F32);

                AE_LSX2X2_IP(y0, y1, py, 4*SZ_F32);
                AE_LSX2X2_IP(y2, y3, py, 4*SZ_F32);

                MADDQ_S(acc00, acc10, x00, x10, y0);
                MADDQ_S(acc00, acc10, x02, x12, y2);

                MADDQ_S(acc01, acc11, x01, x11, y1);
                MADDQ_S(acc01, acc11, x03, x13, y3);

                MADDQ_S(acc20, acc30, x20, x30, y0);
                MADDQ_S(acc20, acc30, x22, x32, y2);

                MADDQ_S(acc21, acc31, x21, x31, y1);
                MADDQ_S(acc21, acc31, x23, x33, y3);

                MADDQ_S(acc40, acc50, x40, x50, y0);
                MADDQ_S(acc40, acc50, x42, x52, y2);

                MADDQ_S(acc41, acc51, x41, x51, y1);
                MADDQ_S(acc41, acc51, x43, x53, y3);

                MADDQ_S(acc60, acc70, x60, x70, y0);
                MADDQ_S(acc60, acc70, x62, x72, y2);

                MADDQ_S(acc61, acc71, x61, x71, y1);
                MADDQ_S(acc61, acc71, x63, x73, y3);
            }

            acc00 = acc00 + acc01;
            acc10 = acc10 + acc11;
            y0 = XT_SEL32_HL_SX2(acc00, acc10);
            y1 = XT_SEL32_LH_SX2(acc00, acc10);
            z0 = y0 + y1;

            acc20 = acc20 + acc21;
            acc30 = acc30 + acc31;
            y0 = XT_SEL32_HL_SX2(acc20, acc30);
            y1 = XT_SEL32_LH_SX2(acc20, acc30);
            z1 = y0 + y1;

            acc40 = acc40 + acc41;
            acc50 = acc50 + acc51;
            y0 = XT_SEL32_HL_SX2(acc40, acc50);
            y1 = XT_SEL32_LH_SX2(acc40, acc50);
            z2 = y0 + y1;

            acc60 = acc60 + acc61;
            acc70 = acc70 + acc71;
            y0 = XT_SEL32_HL_SX2(acc60, acc70);
            y1 = XT_SEL32_LH_SX2(acc60, acc70);
            z3 = y0 + y1;

            AE_SASX2X2_IP(z0, z1, v_z, pz);
            AE_SASX2X2_IP(z2, z3, v_z, pz);
            /* Go to the next 8 rows of the input matrix */
            px0 = (const xtfloatx4 *)((float32_t *)px0 + 7*N);
        }
        if ((M >> 2) & 1) {
            px1 = (const xtfloatx4 *)((float32_t *)px0 + 4);
            py  = (const xtfloatx4 *)(y);

            acc00 = acc01 = acc10 = acc11 =
            acc20 = acc21 = acc30 = acc31 = (xtfloatx2)0.0f;

            __Pragma("loop_count min=1")
            for (n = 0; n < (N >> 3); n++)
            {
                AE_LSX2X2_XP(x00, x01, px0, N*SZ_F32);
                AE_LSX2X2_XP(x02, x03, px1, N*SZ_F32);

                AE_LSX2X2_XP(x10, x11, px0, N*SZ_F32);
                AE_LSX2X2_XP(x12, x13, px1, N*SZ_F32);

                AE_LSX2X2_XP(x20, x21, px0, N*SZ_F32);
                AE_LSX2X2_XP(x22, x23, px1, N*SZ_F32);

                AE_LSX2X2_XP(x30, x31, px0,(8-3*N)*SZ_F32);
                AE_LSX2X2_XP(x32, x33, px1,(8-3*N)*SZ_F32);

                AE_LSX2X2_IP(y0, y1, py, 4*SZ_F32);
                AE_LSX2X2_IP(y2, y3, py, 4*SZ_F32);

                MADDQ_S(acc00, acc10, x00, x10, y0);
                MADDQ_S(acc00, acc10, x02, x12, y2);

                MADDQ_S(acc01, acc11, x01, x11, y1);
                MADDQ_S(acc01, acc11, x03, x13, y3);

                MADDQ_S(acc20, acc30, x20, x30, y0);
                MADDQ_S(acc20, acc30, x22, x32, y2);

                MADDQ_S(acc21, acc31, x21, x31, y1);
                MADDQ_S(acc21, acc31, x23, x33, y3);
            }

            acc00 = acc00 + acc01;
            acc10 = acc10 + acc11;
            y0 = XT_SEL32_HL_SX2(acc00, acc10);
            y1 = XT_SEL32_LH_SX2(acc00, acc10);
            z0 = y0 + y1;

            acc20 = acc20 + acc21;
            acc30 = acc30 + acc31;
            y0 = XT_SEL32_HL_SX2(acc20, acc30);
            y1 = XT_SEL32_LH_SX2(acc20, acc30);
            z1 = y0 + y1;

            AE_SASX2X2_IP(z0, z1, v_z, pz);
        }
    }
    AE_SA128POS_FP(v_z, pz);
} /* mtx_vecmpyf() */
#elif (HAVE_FPU)
void mtx_vecmpyf_fast( float32_t * z, const float32_t * x,  const float32_t * y, int M, int N )
{
    int m, n;
    xtfloat * restrict pZ;
    const xtfloat * restrict pX;
    const xtfloat * restrict pY;

    NASSERT(x);
    NASSERT(y);
    NASSERT(z);
    NASSERT((z != x) && (z != y));
    NASSERT_ALIGN(x,8);
    NASSERT_ALIGN(y,8);
    NASSERT(M%4==0);
    NASSERT(N%4==0);
    pZ=(xtfloat*)z;
    if (N<0)
    {
        for ( m=0; m<M; m++) XT_SSIP(XT_CONST_S(0),pZ,sizeof(xtfloat));
        return;
    }
    for ( m=0; m<M; m+=4)
    {
        xtfloat a0,a1,a2,a3,x0,x1,x2,x3,y0;
        a0=a1=a2=a3=XT_CONST_S(0);
        pX=(const xtfloat*)(x+m*N);
        pY=(const xtfloat*)(y);
        for ( n=0; n<N; n++ )
        {
            XT_LSIP(y0,pY,sizeof(xtfloat));
            x1=XT_LSX(pX,1*N*sizeof(xtfloat));
            x2=XT_LSX(pX,2*N*sizeof(xtfloat));
            x3=XT_LSX(pX,3*N*sizeof(xtfloat));
            XT_LSIP(x0,pX,sizeof(xtfloat));
            XT_MADD_S(a0,x0,y0);
            XT_MADD_S(a1,x1,y0);
            XT_MADD_S(a2,x2,y0);
            XT_MADD_S(a3,x3,y0);
        }
        XT_SSIP(a0,pZ,sizeof(xtfloat));
        XT_SSIP(a1,pZ,sizeof(xtfloat));
        XT_SSIP(a2,pZ,sizeof(xtfloat));
        XT_SSIP(a3,pZ,sizeof(xtfloat));
    }
}

#endif
