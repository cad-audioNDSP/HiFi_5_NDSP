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
    Vector by Vector Multiply 
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
DISCARD_FUN(void,cmtx_vecmpytf_fast,( complex_float* y,
                    const complex_float* x,
                    int N ))
#elif (HAVE_VFPU)
/*-------------------------------------------------------------------------
  Vector by Vector Multiply 
  These functions compute the expression y = 2^lsh*x*xt for the input column
  vector x and its Hermitian transpose xt.

  NOTE: lsh factor is not relevant for floating point routines.

  Precision:
  32x32   32-bit input, 32-bit output
  f       floating point

  Input:
  x[N]    input vector, Q31 or floating point
  N       size of vector x
  lsh     bidirectional left shift applied to the result (fixed point 
          functions only).
  Output:
  y[N*N]  output matrix, Q31 or floating point

  Restrictions:
  x,y     should not overlap
  x,y     aligned on 16-byte boundary
  N       multiple of 32
  lsh     -31...31
-------------------------------------------------------------------------*/
void cmtx_vecmpytf_fast ( complex_float* y,
                    const complex_float* x,
                    int N )
{
    const xtfloatx4 * restrict px0;
    const xtfloatx4 * restrict py0;
    xtfloatx4 * restrict pz0;
    xtfloatx4 * restrict pz1;
    xtfloatx4 * restrict pz2;
    xtfloatx4 * restrict pz3;
    int m, n;

    NASSERT_ALIGN(x, HIFI_SIMD_WIDTH);
    NASSERT_ALIGN(y, HIFI_SIMD_WIDTH);
    NASSERT(0==(N%32));

    px0 = (const xtfloatx4 *)(x);
    for (m = 0; m < N; m+=4)
    {
        xtfloatx2 x0, x1, x2, x3;
        xtfloatx2 y0, y1, y2, y3;
        xtfloatx2 z00, z01, z02, z03, z10, z11, z12, z13;
        xtfloatx2 z20, z21, z22, z23, z30, z31, z32, z33;

        py0 = px0;
        pz0 = (xtfloatx4 *)(y + m*N + m);
        pz1 = pz0 + 1;
        pz2 = (xtfloatx4 *)((complex_float *)pz0 + N*4);
        pz3 = pz2 + 1;

        AE_LSX2X2_IP(x0, x1, px0, 2*sizeof(complex_float));
        AE_LSX2X2_IP(x2, x3, px0, 2*sizeof(complex_float));

        {
            py0 += 2;

            /* perform multiplications */
            MULCCONJ_SX2(z01, z02, x0, x0, x1, x2);
            MULCCONJ_SX2(z03, z12, x0, x1, x3, x2);
            MULCCONJ_SX2(z13, z23, x1, x2, x3, x3);
            MULCCONJ_SX2(z00, z11, x0, x1, x0, x1);
            MULCCONJ_SX2(z22, z33, x2, x3, x2, x3);
            CONJC_SX2X2(z10, z20, z01, z02);
            CONJC_SX2X2(z21, z30, z12, z03);
            CONJC_SX2X2(z31, z32, z13, z23);
            /* set imaginary part of diagonal elements to zero */
            MULQ_S(z00, z11, z00, z11, FLOAT_SX2(AE_MOVDA32X2(1,0),0));
            MULQ_S(z22, z33, z22, z33, FLOAT_SX2(AE_MOVDA32X2(1,0),0));

            /* save values */
            AE_SSX2X2_XP(z00, z01, pz0, N*sizeof(complex_float));
            AE_SSX2X2_XP(z10, z11, pz0, N*sizeof(complex_float));
            AE_SSX2X2_XP(z20, z21, pz0, N*sizeof(complex_float));
            AE_SSX2X2_XP(z30, z31, pz0, (-3*N+4)*(int)sizeof(complex_float));
            AE_SSX2X2_XP(z02, z03, pz1, N*sizeof(complex_float));
            AE_SSX2X2_XP(z12, z13, pz1, N*sizeof(complex_float));
            AE_SSX2X2_XP(z22, z23, pz1, N*sizeof(complex_float));
            AE_SSX2X2_XP(z32, z33, pz1, (-3*N+4)*(int)sizeof(complex_float));
        }
        for (n = 0; n < ((N-m-4) >> 2); n++)
        {
            /* load data */
            AE_LSX2X2_IP(y0, y1, py0, 2*sizeof(complex_float));
            AE_LSX2X2_IP(y2, y3, py0, 2*sizeof(complex_float));

            /* perform multiplications */
            MULMUXQ_S(z00, z10, x0, x1, y0, 1);  MADDMUXQ_S(z00, z10, x0, x1, y0, 5);
            MULMUXQ_S(z20, z30, x2, x3, y0, 1);  MADDMUXQ_S(z20, z30, x2, x3, y0, 5);
            MULMUXQ_S(z01, z11, x0, x1, y1, 1);  MADDMUXQ_S(z01, z11, x0, x1, y1, 5);
            MULMUXQ_S(z21, z31, x2, x3, y1, 1);  MADDMUXQ_S(z21, z31, x2, x3, y1, 5);
            MULMUXQ_S(z02, z12, x0, x1, y2, 1);  MADDMUXQ_S(z02, z12, x0, x1, y2, 5);
            MULMUXQ_S(z22, z32, x2, x3, y2, 1);  MADDMUXQ_S(z22, z32, x2, x3, y2, 5);
            MULMUXQ_S(z03, z13, x0, x1, y3, 1);  MADDMUXQ_S(z03, z13, x0, x1, y3, 5);
            MULMUXQ_S(z23, z33, x2, x3, y3, 1);  MADDMUXQ_S(z23, z33, x2, x3, y3, 5);

            /* save values */
            AE_SSX2X2_XP(z00, z01, pz0, N*sizeof(complex_float));
            AE_SSX2X2_XP(z10, z11, pz0, N*sizeof(complex_float));
            AE_SSX2X2_XP(z20, z21, pz0, N*sizeof(complex_float));
            AE_SSX2X2_XP(z30, z31, pz0, (-3*N+4)*(int)sizeof(complex_float));
            AE_SSX2X2_XP(z02, z03, pz1, N*sizeof(complex_float));
            AE_SSX2X2_XP(z12, z13, pz1, N*sizeof(complex_float));
            AE_SSX2X2_XP(z22, z23, pz1, N*sizeof(complex_float));
            AE_SSX2X2_XP(z32, z33, pz1, (-3*N+4)*(int)sizeof(complex_float));
            CONJC_SX2X2(z00, z01, z00, z01);
            CONJC_SX2X2(z10, z11, z10, z11);
            CONJC_SX2X2(z20, z21, z20, z21);
            CONJC_SX2X2(z30, z31, z30, z31);
            CONJC_SX2X2(z02, z03, z02, z03);
            CONJC_SX2X2(z12, z13, z12, z13);
            CONJC_SX2X2(z22, z23, z22, z23);
            CONJC_SX2X2(z32, z33, z32, z33);
            AE_SSX2X2_XP(z00, z10, pz2, N*sizeof(complex_float));
            AE_SSX2X2_XP(z01, z11, pz2, N*sizeof(complex_float));
            AE_SSX2X2_XP(z02, z12, pz2, N*sizeof(complex_float));
            AE_SSX2X2_XP(z03, z13, pz2, N*sizeof(complex_float));
            AE_SSX2X2_XP(z20, z30, pz3, N*sizeof(complex_float));
            AE_SSX2X2_XP(z21, z31, pz3, N*sizeof(complex_float));
            AE_SSX2X2_XP(z22, z32, pz3, N*sizeof(complex_float));
            AE_SSX2X2_XP(z23, z33, pz3, N*sizeof(complex_float));
        }
    }
} /* cmtx_vecmpytf_fast() */
#else
void cmtx_vecmpytf_fast ( complex_float* y,
                    const complex_float* x,
                    int N )
{
    const float32_t * restrict px0;
    const float32_t * restrict py0;
    float32_t * restrict pz0;
    float32_t * restrict pz1;
    float32_t * restrict pz2;
    float32_t x0re, x0im, x1re, x1im;
    float32_t y0re, y0im, y1re, y1im;
    float32_t z00re, z00im, z01re, z01im;
    float32_t z10re, z10im, z11re, z11im;
    int m, n;

    NASSERT_ALIGN(x, HIFI_SIMD_WIDTH);
    NASSERT_ALIGN(y, HIFI_SIMD_WIDTH);
    NASSERT(0==(N%32));

    px0 = (const float32_t *)(x);
    for ( m=0; m<N; m+=2 )
    {
        pz0 = (float32_t *)(y + m*N + m);
        pz1 = pz0 + N*2;
        pz2 = pz1 + N*2;

        XT_LSIP(x0re, px0, sizeof(float32_t));  XT_LSIP(x0im, px0, sizeof(float32_t));
        XT_LSIP(x1re, px0, sizeof(float32_t));  XT_LSIP(x1im, px0, sizeof(float32_t));

        z00re = x0re*x0re + x0im*x0im;
        z00im = (float32_t)XT_CONST_S(0)*z00re;/* multiply by the real part for correct propagation of NaN */
        z01re = x0re*x1re + x0im*x1im;
        z01im = x0im*x1re - x0re*x1im;
        z10re =  z01re;
        z10im = -z01im;
        z11re = x1re*x1re + x1im*x1im;
        z11im = (float32_t)XT_CONST_S(0)*z11re;/* multiply by the real part for correct propagation of NaN */

        XT_SSIP(z00re, pz0, sizeof(float32_t));  XT_SSIP(z00im, pz0, sizeof(float32_t));
        XT_SSIP(z01re, pz0, sizeof(float32_t));  XT_SSIP(z01im, pz0, sizeof(float32_t));
        XT_SSIP(z10re, pz1, sizeof(float32_t));  XT_SSIP(z10im, pz1, sizeof(float32_t));
        XT_SSIP(z11re, pz1, sizeof(float32_t));  XT_SSIP(z11im, pz1, sizeof(float32_t));

        py0 = px0;
        for ( n=0; n<((N-m-2)>>1); n++ )
        {
            XT_LSIP(y0re, py0, sizeof(float32_t));  XT_LSIP(y0im, py0, sizeof(float32_t));
            XT_LSIP(y1re, py0, sizeof(float32_t));  XT_LSIP(y1im, py0, sizeof(float32_t));

            z00re = x0re*y0re + x0im*y0im;
            z00im = x0im*y0re - x0re*y0im;
            z01re = x0re*y1re + x0im*y1im;
            z01im = x0im*y1re - x0re*y1im;
            z10re = x1re*y0re + x1im*y0im;
            z10im = x1im*y0re - x1re*y0im;
            z11re = x1re*y1re + x1im*y1im;
            z11im = x1im*y1re - x1re*y1im;

            XT_SSIP(z00re, pz0, sizeof(float32_t));  XT_SSIP(z00im, pz0, sizeof(float32_t));
            XT_SSIP(z01re, pz0, sizeof(float32_t));  XT_SSIP(z01im, pz0, sizeof(float32_t));
            XT_SSIP(z10re, pz1, sizeof(float32_t));  XT_SSIP(z10im, pz1, sizeof(float32_t));
            XT_SSIP(z11re, pz1, sizeof(float32_t));  XT_SSIP(z11im, pz1, sizeof(float32_t));
            z00im = -z00im;
            z01im = -z01im;
            z10im = -z10im;
            z11im = -z11im;
            XT_SSIP(z00re, pz2, sizeof(float32_t));  XT_SSIP(z00im, pz2, sizeof(float32_t));
            XT_SSIP(z10re, pz2, sizeof(float32_t));  XT_SSXP(z10im, pz2, (N*2-3)*(int)sizeof(float32_t));
            XT_SSIP(z01re, pz2, sizeof(float32_t));  XT_SSIP(z01im, pz2, sizeof(float32_t));
            XT_SSIP(z11re, pz2, sizeof(float32_t));  XT_SSXP(z11im, pz2, (N*2-3)*(int)sizeof(float32_t));
        }
    }
}
#endif
