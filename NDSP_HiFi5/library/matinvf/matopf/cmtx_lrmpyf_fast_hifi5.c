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
  NatureDSP Signal Processing Library. Matrix operations part
    Three Matrices Product
    Floating-point complex data variant
    C code optimized for HiFi5 core
  IntegrIT, 2006-2019
*/

/* Portable data types. */
#include "NatureDSP_types.h"
/* Common utility and macros declarations. */
#include "common.h"
#include "common_fpu.h"
/* Signal Processing Library API. */
#include "NatureDSP_Signal_matop.h"

#if (HAVE_FPU==0 && HAVE_VFPU==0)
DISCARD_FUN(void,cmtx_lrmpyf_fast,( void* pScr,
                        complex_float* z,
                  const complex_float* x, 
                  const complex_float* y,
                  int N ))
#elif (HAVE_VFPU)
/*-------------------------------------------------------------------------
  Three Matrices Product
  These functions compute the expression z = 2^lsh*((2^-rsh*x*y)*xt) for 
  square input matrices x and y, where xt is the Hermitian transpose of the
  matrix x.

  NOTE: 2^lsh and 2^-rsh factors are not relevant for floating point routines.

  Functions require scratch memory for storing intermediate data. This scratch
  memory area should be aligned on 16-byte boundary; its size is calculated
  by the corresponding scratch allocation function.

  Precision:
  32x32   32-bit input, 32-bit output
  f       floating point

  Input:
  x[N*N]  input matrix x, Q31 or floating point
  y[N*N]  input matrix y, Q31 or floating point
  N       number of rows and columns in matrices x, y and z
  rsh     right shift applied to intermediate results to avoid overflow (fixed 
          point functions only)
  lsh     bidirectional left shift applied to the result (fixed point functions
          only)
  Output:
  z[N*N]  output matrix z, Q31 or floating point
  Scratch:
  pScr    scratch memory area with size in bytes defined by the corresponding
          scratch allocation function

  Restrictions:
  x,y,z   should not overlap
  x,y,z   aligned on 16-byte boundary
  N       multiple of 32
  rsh     >=0
  lsh     -31...31
-------------------------------------------------------------------------*/
void cmtx_lrmpyf_fast ( void* pScr,
                        complex_float* z,
                  const complex_float* x, 
                  const complex_float* y,
                  int N )
{
    complex_float * u;
    const xtfloatx4 * restrict px0;
    const xtfloatx2 * restrict py0;
    xtfloatx4 * restrict pz0;
    xtfloatx4 * restrict pz1;
    xtfloatx4 * restrict pz2;
    xtfloatx4 * restrict pz3;

    int m, n, p;

    NASSERT_ALIGN(pScr, HIFI_SIMD_WIDTH);
    NASSERT_ALIGN(z, HIFI_SIMD_WIDTH);
    NASSERT_ALIGN(x, HIFI_SIMD_WIDTH);
    NASSERT_ALIGN(y, HIFI_SIMD_WIDTH);
    NASSERT(0==(N%32));

    if (N <= 0) return;
    u = (complex_float *)pScr;

    pz3 = (xtfloatx4 *)(u);
    __Pragma("loop_count min=1");
    for (p = 0; p < N; p+=4)
    {
        xtfloatx2 x00, x10, x20, x30, x01, x11, x21, x31;
        xtfloatx2 y00, y01, y02, y03, y10, y11, y12, y13;
        xtfloatx2 z00, z01, z02, z03, z10, z11, z12, z13;
        xtfloatx2 z20, z21, z22, z23, z30, z31, z32, z33;

        pz0 = pz3;
        pz1 = (xtfloatx4 *)((complex_float *)pz0 + N);
        pz2 = (xtfloatx4 *)((complex_float *)pz1 + N);
        pz3 = (xtfloatx4 *)((complex_float *)pz2 + N);
        __Pragma("loop_count min=1");
        for (m = 0; m < N; m+=4)
        {
            px0 = (const xtfloatx4 *)(y + m*N);
            py0 = (const xtfloatx2 *)(x + p*N);

            z00 = z01 = z02 = z03 = CONST_SX2(0);
            z10 = z11 = z12 = z13 = CONST_SX2(0);
            z20 = z21 = z22 = z23 = CONST_SX2(0);
            z30 = z31 = z32 = z33 = CONST_SX2(0);

            __Pragma("loop_count min=2, factor=2");
            __Pragma("ymemory(py0)");
            __Pragma("ymemory(px0)");
            for (n = 0; n < (N >> 1); n++)
            {
                /* load data */
                AE_LSX2X2_XP(x00, x01, px0, N*sizeof(complex_float));
                AE_LSX2X2_XP(x10, x11, px0, N*sizeof(complex_float));
                AE_LSX2X2_XP(x20, x21, px0, N*sizeof(complex_float));
                AE_LSX2X2_XP(x30, x31, px0, (-3*N+2)*sizeof(complex_float));

                AE_LSX2XP(y00, py0, N*sizeof(complex_float));
                AE_LSX2XP(y01, py0, N*sizeof(complex_float));
                AE_LSX2XP(y02, py0, N*sizeof(complex_float));
                AE_LSX2XP(y03, py0, (-3*N+1)*(int)sizeof(complex_float));
                AE_LSX2XP(y10, py0, N*sizeof(complex_float));
                AE_LSX2XP(y11, py0, N*sizeof(complex_float));
                AE_LSX2XP(y12, py0, N*sizeof(complex_float));
                AE_LSX2XP(y13, py0, (-3*N+1)*(int)sizeof(complex_float));

                /* perform multiplications */
                MADDCCONJ_SX2(z00, z01, x00, x00, y00, y01);
                MADDCCONJ_SX2(z02, z03, x00, x00, y02, y03);
                MADDCCONJ_SX2(z10, z11, x10, x10, y00, y01);
                MADDCCONJ_SX2(z12, z13, x10, x10, y02, y03);
                MADDCCONJ_SX2(z20, z21, x20, x20, y00, y01);
                MADDCCONJ_SX2(z22, z23, x20, x20, y02, y03);
                MADDCCONJ_SX2(z30, z31, x30, x30, y00, y01);
                MADDCCONJ_SX2(z32, z33, x30, x30, y02, y03);

                MADDCCONJ_SX2(z00, z01, x01, x01, y10, y11);
                MADDCCONJ_SX2(z02, z03, x01, x01, y12, y13);
                MADDCCONJ_SX2(z10, z11, x11, x11, y10, y11);
                MADDCCONJ_SX2(z12, z13, x11, x11, y12, y13);
                MADDCCONJ_SX2(z20, z21, x21, x21, y10, y11);
                MADDCCONJ_SX2(z22, z23, x21, x21, y12, y13);
                MADDCCONJ_SX2(z30, z31, x31, x31, y10, y11);
                MADDCCONJ_SX2(z32, z33, x31, x31, y12, y13);
            }
            /* save values */
            AE_SSX2X2_IP(z00, z10, pz0, 2*sizeof(complex_float));
            AE_SSX2X2_IP(z20, z30, pz0, 2*sizeof(complex_float));
            AE_SSX2X2_IP(z01, z11, pz1, 2*sizeof(complex_float));
            AE_SSX2X2_IP(z21, z31, pz1, 2*sizeof(complex_float));
            AE_SSX2X2_IP(z02, z12, pz2, 2*sizeof(complex_float));
            AE_SSX2X2_IP(z22, z32, pz2, 2*sizeof(complex_float));
            AE_SSX2X2_IP(z03, z13, pz3, 2*sizeof(complex_float));
            AE_SSX2X2_IP(z23, z33, pz3, 2*sizeof(complex_float));
        }
    }

    pz3 = (xtfloatx4 *)(z);
    __Pragma("loop_count min=1");
    for (m = 0; m < N; m+=4)
    {
        xtfloatx2 x00, x10, x20, x30, x01, x11, x21, x31;
        xtfloatx2 y00, y01, y02, y03, y10, y11, y12, y13;
        xtfloatx2 z00, z01, z02, z03, z10, z11, z12, z13;
        xtfloatx2 z20, z21, z22, z23, z30, z31, z32, z33;

        pz0 = pz3;
        pz1 = (xtfloatx4 *)((complex_float *)pz0 + N);
        pz2 = (xtfloatx4 *)((complex_float *)pz1 + N);
        pz3 = (xtfloatx4 *)((complex_float *)pz2 + N);
        __Pragma("loop_count min=1");
        for (p = 0; p < N; p+=4)
        {
            px0 = (const xtfloatx4 *)(x + m*N);
            py0 = (const xtfloatx2 *)(u + p*N);

            z00 = z01 = z02 = z03 = CONST_SX2(0);
            z10 = z11 = z12 = z13 = CONST_SX2(0);
            z20 = z21 = z22 = z23 = CONST_SX2(0);
            z30 = z31 = z32 = z33 = CONST_SX2(0);

            __Pragma("loop_count min=2, factor=2");
            __Pragma("ymemory(py0)");
            __Pragma("ymemory(px0)");
            for (n = 0; n < (N >> 1); n++)
            {
                /* load data */
                AE_LSX2X2_XP(x00, x01, px0, N*sizeof(complex_float));
                AE_LSX2X2_XP(x10, x11, px0, N*sizeof(complex_float));
                AE_LSX2X2_XP(x20, x21, px0, N*sizeof(complex_float));
                AE_LSX2X2_XP(x30, x31, px0, (-3*N+2)*sizeof(complex_float));

                AE_LSX2XP(y00, py0, N*sizeof(complex_float));
                AE_LSX2XP(y01, py0, N*sizeof(complex_float));
                AE_LSX2XP(y02, py0, N*sizeof(complex_float));
                AE_LSX2XP(y03, py0, (-3*N+1)*(int)sizeof(complex_float));
                AE_LSX2XP(y10, py0, N*sizeof(complex_float));
                AE_LSX2XP(y11, py0, N*sizeof(complex_float));
                AE_LSX2XP(y12, py0, N*sizeof(complex_float));
                AE_LSX2XP(y13, py0, (-3*N+1)*(int)sizeof(complex_float));

                /* perform multiplications */
                MADDC_SX2(z00, z01, x00, x00, y00, y01);
                MADDC_SX2(z02, z03, x00, x00, y02, y03);
                MADDC_SX2(z10, z11, x10, x10, y00, y01);
                MADDC_SX2(z12, z13, x10, x10, y02, y03);
                MADDC_SX2(z20, z21, x20, x20, y00, y01);
                MADDC_SX2(z22, z23, x20, x20, y02, y03);
                MADDC_SX2(z30, z31, x30, x30, y00, y01);
                MADDC_SX2(z32, z33, x30, x30, y02, y03);

                MADDC_SX2(z00, z01, x01, x01, y10, y11);
                MADDC_SX2(z02, z03, x01, x01, y12, y13);
                MADDC_SX2(z10, z11, x11, x11, y10, y11);
                MADDC_SX2(z12, z13, x11, x11, y12, y13);
                MADDC_SX2(z20, z21, x21, x21, y10, y11);
                MADDC_SX2(z22, z23, x21, x21, y12, y13);
                MADDC_SX2(z30, z31, x31, x31, y10, y11);
                MADDC_SX2(z32, z33, x31, x31, y12, y13);
            }
            /* save values */
            AE_SSX2X2_IP(z00, z01, pz0, 2*sizeof(complex_float));
            AE_SSX2X2_IP(z02, z03, pz0, 2*sizeof(complex_float));
            AE_SSX2X2_IP(z10, z11, pz1, 2*sizeof(complex_float));
            AE_SSX2X2_IP(z12, z13, pz1, 2*sizeof(complex_float));
            AE_SSX2X2_IP(z20, z21, pz2, 2*sizeof(complex_float));
            AE_SSX2X2_IP(z22, z23, pz2, 2*sizeof(complex_float));
            AE_SSX2X2_IP(z30, z31, pz3, 2*sizeof(complex_float));
            AE_SSX2X2_IP(z32, z33, pz3, 2*sizeof(complex_float));
        }
    }
} /* cmtx_lrmpyf_fast() */
#else
void cmtx_lrmpyf_fast ( void* pScr,
                        complex_float* z,
                  const complex_float* x, 
                  const complex_float* y,
                  int N )
{
    complex_float * u;
    const float32_t * restrict px0;
    const float32_t * restrict px1;
    const float32_t * restrict py0;
    const float32_t * restrict py1;
    float32_t * restrict pz0;
    float32_t * restrict pz1;
    float32_t x00re, x00im;
    float32_t x10re, x10im;
    float32_t y00re, y00im, y01re, y01im;
    float32_t z00re, z00im, z01re, z01im;
    float32_t z10re, z10im, z11re, z11im;
    int m, n, p;

    NASSERT_ALIGN(pScr, HIFI_SIMD_WIDTH);
    NASSERT_ALIGN(z, HIFI_SIMD_WIDTH);
    NASSERT_ALIGN(x, HIFI_SIMD_WIDTH);
    NASSERT_ALIGN(y, HIFI_SIMD_WIDTH);
    NASSERT(0==(N%32));

    if (N<=0) return;
    u = (complex_float *)pScr;

    pz1 = (float32_t *)(u);
    __Pragma("loop_count min=1");
    for ( p=0; p<N; p+=2 )
    {
        pz0 = pz1;
        pz1 = pz0 + N*2;

        __Pragma("loop_count min=1");
        for ( m=0; m<N; m+=2 )
        {
            px0 = (const float32_t *)(y + m*N);
            px1 = px0 + N*2;
            py0 = (const float32_t *)(x + p*N);
            py1 = py0 + N*2;

            z00re = z00im = z01re = z01im = z10re = z10im = z11re = z11im = 0.0f;
            __Pragma("loop_count min=1");
            for ( n=0; n<N; n++ )
            {
                XT_LSIP(x00re, px0, sizeof(float32_t));  XT_LSIP(x00im, px0, sizeof(float32_t));
                XT_LSIP(x10re, px1, sizeof(float32_t));  XT_LSIP(x10im, px1, sizeof(float32_t));

                XT_LSIP(y00re, py0, sizeof(float32_t));  XT_LSIP(y00im, py0, sizeof(float32_t));
                XT_LSIP(y01re, py1, sizeof(float32_t));  XT_LSIP(y01im, py1, sizeof(float32_t));

                z00re += x00re*y00re + x00im*y00im;
                z00im += x00im*y00re - x00re*y00im;
                z01re += x00re*y01re + x00im*y01im;
                z01im += x00im*y01re - x00re*y01im;
                z10re += x10re*y00re + x10im*y00im;
                z10im += x10im*y00re - x10re*y00im;
                z11re += x10re*y01re + x10im*y01im;
                z11im += x10im*y01re - x10re*y01im;
            }
            XT_SSIP(z00re, pz0, sizeof(float32_t));  XT_SSIP(z00im, pz0, sizeof(float32_t));
            XT_SSIP(z10re, pz0, sizeof(float32_t));  XT_SSIP(z10im, pz0, sizeof(float32_t));
            XT_SSIP(z01re, pz1, sizeof(float32_t));  XT_SSIP(z01im, pz1, sizeof(float32_t));
            XT_SSIP(z11re, pz1, sizeof(float32_t));  XT_SSIP(z11im, pz1, sizeof(float32_t));
        }
    }

    pz1 = (float32_t *)(z);
    __Pragma("loop_count min=1");
    for ( m=0; m<N; m+=2 )
    {
        pz0 = pz1;
        pz1 = pz0 + N*2;

        __Pragma("loop_count min=1");
        for ( p=0; p<N; p+=2 )
        {
            px0 = (const float32_t *)(x + m*N);
            px1 = px0 + N*2;
            py0 = (const float32_t *)(u + p*N);
            py1 = py0 + N*2;

            z00re = z00im = z01re = z01im = z10re = z10im = z11re = z11im = 0.0f;
            __Pragma("loop_count min=1");
            for ( n=0; n<N; n++ )
            {
                XT_LSIP(x00re, px0, sizeof(float32_t));  XT_LSIP(x00im, px0, sizeof(float32_t));
                XT_LSIP(x10re, px1, sizeof(float32_t));  XT_LSIP(x10im, px1, sizeof(float32_t));

                XT_LSIP(y00re, py0, sizeof(float32_t));  XT_LSIP(y00im, py0, sizeof(float32_t));
                XT_LSIP(y01re, py1, sizeof(float32_t));  XT_LSIP(y01im, py1, sizeof(float32_t));

                z00re += x00re*y00re - x00im*y00im;
                z00im += x00im*y00re + x00re*y00im;
                z01re += x00re*y01re - x00im*y01im;
                z01im += x00im*y01re + x00re*y01im;
                z10re += x10re*y00re - x10im*y00im;
                z10im += x10im*y00re + x10re*y00im;
                z11re += x10re*y01re - x10im*y01im;
                z11im += x10im*y01re + x10re*y01im;
            }
            XT_SSIP(z00re, pz0, sizeof(float32_t));  XT_SSIP(z00im, pz0, sizeof(float32_t));
            XT_SSIP(z01re, pz0, sizeof(float32_t));  XT_SSIP(z01im, pz0, sizeof(float32_t));
            XT_SSIP(z10re, pz1, sizeof(float32_t));  XT_SSIP(z10im, pz1, sizeof(float32_t));
            XT_SSIP(z11re, pz1, sizeof(float32_t));  XT_SSIP(z11im, pz1, sizeof(float32_t));
        }
    }
}
#endif

size_t cmtx_lrmpyf_fast_getScratchSize(int N)
{
    NASSERT(0==(N%32));
    return N>0 ? N*N*sizeof(complex_float) : 0;
}
