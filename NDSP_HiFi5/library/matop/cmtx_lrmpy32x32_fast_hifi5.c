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
    32-bit fixed-point complex data variant
    C code optimized for HiFi5 core
  IntegrIT, 2006-2019
*/

/* Portable data types. */
#include "NatureDSP_types.h"
/* Common utility and macros declarations. */
#include "common.h"
/* Signal Processing Library API. */
#include "NatureDSP_Signal_matop.h"

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
#if (0)
DISCARD_FUN(void,cmtx_lrmpy32x32_fast,( void* pScr,
                        complex_fract32* z,
                  const complex_fract32* x, 
                  const complex_fract32* y,
                  int N, int rsh, int lsh ))
#else
void cmtx_lrmpy32x32_fast ( void* pScr,
                        complex_fract32* z,
                  const complex_fract32* x, 
                  const complex_fract32* y,
                  int N, int rsh, int lsh )
{
    complex_fract32 * u;
    const ae_int32x4 * restrict px0;
    const ae_int32x4 * restrict py0;
    const ae_int32x4 * restrict py1;
    ae_int32x4 * restrict pz0;
    ae_int32x4 * restrict pz1;
    ae_int32x2 x00, x01, x02, x03, x10, x11, x12, x13;
    ae_int32x2 y00, y01, y02, y03, y10, y11, y12, y13;
    ae_int32x2 z00, z01, z02, z03, z10, z11, z12, z13;
    ae_int64 Z00re, Z00im, Z01re, Z01im;
    ae_int64 Z02re, Z02im, Z03re, Z03im;
    ae_int64 Z10re, Z10im, Z11re, Z11im;
    ae_int64 Z12re, Z12im, Z13re, Z13im;
    int m, n, p;

    NASSERT_ALIGN(pScr, HIFI_SIMD_WIDTH);
    NASSERT_ALIGN(z, HIFI_SIMD_WIDTH);
    NASSERT_ALIGN(x, HIFI_SIMD_WIDTH);
    NASSERT_ALIGN(y, HIFI_SIMD_WIDTH);
    NASSERT(0==(N%32));
    NASSERT(rsh>=0);
    NASSERT(-31<=lsh && lsh<=31);

    if (N <= 0) return;
    u = (complex_fract32 *)pScr;

    WUR_AE_SAR(rsh);
    pz1 = (ae_int32x4 *)(u);
    __Pragma("loop_count min=1");
    for (m = 0; m < N; m+=2)
    {
        pz0 = pz1;
        pz1 = (ae_int32x4 *)((complex_fract32 *)pz0 + N);
        __Pragma("loop_count min=1");
        for (p = 0; p < N; p+=4)
        {
            px0 = (const ae_int32x4 *)(x + m*N);
            py0 = (const ae_int32x4 *)(y + p);
            py1 = py0 + 1;

            Z00re = Z00im = Z01re = Z01im = AE_ZERO64();
            Z02re = Z02im = Z03re = Z03im = AE_ZERO64();
            Z10re = Z10im = Z11re = Z11im = AE_ZERO64();
            Z12re = Z12im = Z13re = Z13im = AE_ZERO64();

            __Pragma("loop_count min=2, factor=2");
            for (n = 0; n < (N >> 1); n++)
            {
                /* load data */
                AE_L32X2X2_XP(x00, x01, px0, N*sizeof(complex_fract32));
                AE_L32X2X2_XP(x10, x11, px0, (-N+2)*sizeof(complex_fract32));

                AE_L32X2X2_XP(y00, y01, py0, N*sizeof(complex_fract32));
                AE_L32X2X2_XP(y02, y03, py1, N*sizeof(complex_fract32));
                AE_L32X2X2_XP(y10, y11, py0, N*sizeof(complex_fract32));
                AE_L32X2X2_XP(y12, y13, py1, N*sizeof(complex_fract32));

                /* perform multiplications */
                AE_MULAFC32RA(Z00re, Z00im, x00, y00);    AE_MULAFC32RA(Z00re, Z00im, x01, y10);
                AE_MULAFC32RA(Z01re, Z01im, x00, y01);    AE_MULAFC32RA(Z01re, Z01im, x01, y11);
                AE_MULAFC32RA(Z02re, Z02im, x00, y02);    AE_MULAFC32RA(Z02re, Z02im, x01, y12);
                AE_MULAFC32RA(Z03re, Z03im, x00, y03);    AE_MULAFC32RA(Z03re, Z03im, x01, y13);
                AE_MULAFC32RA(Z10re, Z10im, x10, y00);    AE_MULAFC32RA(Z10re, Z10im, x11, y10);
                AE_MULAFC32RA(Z11re, Z11im, x10, y01);    AE_MULAFC32RA(Z11re, Z11im, x11, y11);
                AE_MULAFC32RA(Z12re, Z12im, x10, y02);    AE_MULAFC32RA(Z12re, Z12im, x11, y12);
                AE_MULAFC32RA(Z13re, Z13im, x10, y03);    AE_MULAFC32RA(Z13re, Z13im, x11, y13);
            }
            Z00re = AE_SRAS64(Z00re);  Z00im = AE_SRAS64(Z00im);
            Z01re = AE_SRAS64(Z01re);  Z01im = AE_SRAS64(Z01im);
            Z02re = AE_SRAS64(Z02re);  Z02im = AE_SRAS64(Z02im);
            Z03re = AE_SRAS64(Z03re);  Z03im = AE_SRAS64(Z03im);
            Z10re = AE_SRAS64(Z10re);  Z10im = AE_SRAS64(Z10im);
            Z11re = AE_SRAS64(Z11re);  Z11im = AE_SRAS64(Z11im);
            Z12re = AE_SRAS64(Z12re);  Z12im = AE_SRAS64(Z12im);
            Z13re = AE_SRAS64(Z13re);  Z13im = AE_SRAS64(Z13im);
            z00 = AE_ROUND32X2F48SASYM(Z00re, Z00im);
            z01 = AE_ROUND32X2F48SASYM(Z01re, Z01im);
            z02 = AE_ROUND32X2F48SASYM(Z02re, Z02im);
            z03 = AE_ROUND32X2F48SASYM(Z03re, Z03im);
            z10 = AE_ROUND32X2F48SASYM(Z10re, Z10im);
            z11 = AE_ROUND32X2F48SASYM(Z11re, Z11im);
            z12 = AE_ROUND32X2F48SASYM(Z12re, Z12im);
            z13 = AE_ROUND32X2F48SASYM(Z13re, Z13im);
            /* save values */
            AE_S32X2X2_IP(z00, z01, pz0, 2*sizeof(complex_fract32));
            AE_S32X2X2_IP(z02, z03, pz0, 2*sizeof(complex_fract32));
            AE_S32X2X2_IP(z10, z11, pz1, 2*sizeof(complex_fract32));
            AE_S32X2X2_IP(z12, z13, pz1, 2*sizeof(complex_fract32));
        }
    }

    WUR_AE_SAR(lsh);
    pz1 = (ae_int32x4 *)(z);
    __Pragma("loop_count min=1");
    for (m = 0; m < N; m+=2)
    {
        pz0 = pz1;
        pz1 = (ae_int32x4 *)((complex_fract32 *)pz0 + N);
        __Pragma("loop_count min=1");
        for (p = 0; p < N; p+=4)
        {
            px0 = (const ae_int32x4 *)(x + p*N);
            py0 = (const ae_int32x4 *)(u + m*N);

            Z00re = Z00im = Z01re = Z01im = AE_ZERO64();
            Z02re = Z02im = Z03re = Z03im = AE_ZERO64();
            Z10re = Z10im = Z11re = Z11im = AE_ZERO64();
            Z12re = Z12im = Z13re = Z13im = AE_ZERO64();

            __Pragma("loop_count min=2, factor=2");
            for (n = 0; n < (N >> 1); n++)
            {
                /* load data */
                AE_L32X2X2_XP(y00, y01, py0, N*sizeof(complex_fract32));
                AE_L32X2X2_XP(y10, y11, py0, (-N+2)*sizeof(complex_fract32));

                AE_L32X2X2_XP(x00, x10, px0, N*sizeof(complex_fract32));
                AE_L32X2X2_XP(x01, x11, px0, N*sizeof(complex_fract32));
                AE_L32X2X2_XP(x02, x12, px0, N*sizeof(complex_fract32));
                AE_L32X2X2_XP(x03, x13, px0, (-3*N+2)*(int)sizeof(complex_fract32));

                /* perform multiplications */
                AE_MULAAF2D32RA_HH_LL(Z00re, Z01re, x00, x01, y00, y00);  AE_MULASF2D32RA_HL_LH(Z00im, Z01im, x00, x01, y00, y00);
                AE_MULAAF2D32RA_HH_LL(Z00re, Z01re, x10, x11, y01, y01);  AE_MULASF2D32RA_HL_LH(Z00im, Z01im, x10, x11, y01, y01);
                AE_MULAAF2D32RA_HH_LL(Z10re, Z11re, x00, x01, y10, y10);  AE_MULASF2D32RA_HL_LH(Z10im, Z11im, x00, x01, y10, y10);
                AE_MULAAF2D32RA_HH_LL(Z10re, Z11re, x10, x11, y11, y11);  AE_MULASF2D32RA_HL_LH(Z10im, Z11im, x10, x11, y11, y11);
                AE_MULAAF2D32RA_HH_LL(Z02re, Z03re, x02, x03, y00, y00);  AE_MULASF2D32RA_HL_LH(Z02im, Z03im, x02, x03, y00, y00);
                AE_MULAAF2D32RA_HH_LL(Z02re, Z03re, x12, x13, y01, y01);  AE_MULASF2D32RA_HL_LH(Z02im, Z03im, x12, x13, y01, y01);
                AE_MULAAF2D32RA_HH_LL(Z12re, Z13re, x02, x03, y10, y10);  AE_MULASF2D32RA_HL_LH(Z12im, Z13im, x02, x03, y10, y10);
                AE_MULAAF2D32RA_HH_LL(Z12re, Z13re, x12, x13, y11, y11);  AE_MULASF2D32RA_HL_LH(Z12im, Z13im, x12, x13, y11, y11);
            }
            Z00re = AE_SLAS64S(Z00re);  Z00im = AE_SLAS64S(Z00im);
            Z01re = AE_SLAS64S(Z01re);  Z01im = AE_SLAS64S(Z01im);
            Z02re = AE_SLAS64S(Z02re);  Z02im = AE_SLAS64S(Z02im);
            Z03re = AE_SLAS64S(Z03re);  Z03im = AE_SLAS64S(Z03im);
            Z10re = AE_SLAS64S(Z10re);  Z10im = AE_SLAS64S(Z10im);
            Z11re = AE_SLAS64S(Z11re);  Z11im = AE_SLAS64S(Z11im);
            Z12re = AE_SLAS64S(Z12re);  Z12im = AE_SLAS64S(Z12im);
            Z13re = AE_SLAS64S(Z13re);  Z13im = AE_SLAS64S(Z13im);
            z00 = AE_ROUND32X2F48SASYM(Z00re, Z00im);
            z01 = AE_ROUND32X2F48SASYM(Z01re, Z01im);
            z02 = AE_ROUND32X2F48SASYM(Z02re, Z02im);
            z03 = AE_ROUND32X2F48SASYM(Z03re, Z03im);
            z10 = AE_ROUND32X2F48SASYM(Z10re, Z10im);
            z11 = AE_ROUND32X2F48SASYM(Z11re, Z11im);
            z12 = AE_ROUND32X2F48SASYM(Z12re, Z12im);
            z13 = AE_ROUND32X2F48SASYM(Z13re, Z13im);
            /* save values */
            AE_S32X2X2_IP(z00, z01, pz0, 2*sizeof(complex_fract32));
            AE_S32X2X2_IP(z02, z03, pz0, 2*sizeof(complex_fract32));
            AE_S32X2X2_IP(z10, z11, pz1, 2*sizeof(complex_fract32));
            AE_S32X2X2_IP(z12, z13, pz1, 2*sizeof(complex_fract32));
        }
    }
} /* cmtx_lrmpy32x32_fast() */
#endif

size_t cmtx_lrmpy32x32_fast_getScratchSize(int N)
{
    NASSERT(0==(N%32));
    return N>0 ? N*N*sizeof(complex_fract32) : 0;
}
