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
    Matrix Multiply Transpose
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
#if (0)
DISCARD_FUN(void,cmtx_mpyt32x32_fast,( void  * pScr,
                 complex_fract32 * z, 
           const complex_fract32 * x, 
           const complex_fract32 * yt,
           int M, int N, int P, int lsh ))
#else
void cmtx_mpyt32x32_fast ( void  * pScr,
                 complex_fract32 * z, 
           const complex_fract32 * x, 
           const complex_fract32 * yt,
           int M, int N, int P, int lsh )
{
    const ae_int32x4 * restrict px0;
    const ae_int32x4 * restrict py0;
    ae_int32x4 * restrict pz0;
    ae_int32x4 * restrict pz1;

    int m, n, p;

    NASSERT_ALIGN(pScr, HIFI_SIMD_WIDTH);
    NASSERT_ALIGN(z, HIFI_SIMD_WIDTH);
    NASSERT_ALIGN(x, HIFI_SIMD_WIDTH);
    NASSERT_ALIGN(yt, HIFI_SIMD_WIDTH);
    NASSERT((0==(M%32)) && (0==(N%32)) && (0==(P%32)));
    NASSERT(-31<=lsh && lsh<=31);

    if ((M <= 0) || (P <= 0)) return;
    if (N <= 0)
    {
        for (m = 0; m < ((M * P)>>1); m++) AE_S32X2X2_IP(AE_ZERO32(),AE_ZERO32(),castxcc(ae_int32x4,z),sizeof(ae_int32x4));
        return;
    }

    WUR_AE_SAR(lsh);

    pz1 = (ae_int32x4 *)(z);
    __Pragma("loop_count min=1");
    for (m = 0; m < M; m+=2)
    {
        ae_int32x2 x00, x01, x10, x11;
        ae_int32x2 y00, y01, y02, y03, y10, y11, y12, y13;
        ae_int32x2 z00, z01, z02, z03, z10, z11, z12, z13;
        ae_int64 Z00re, Z00im, Z01re, Z01im;
        ae_int64 Z02re, Z02im, Z03re, Z03im;
        ae_int64 Z10re, Z10im, Z11re, Z11im;
        ae_int64 Z12re, Z12im, Z13re, Z13im;

        pz0 = pz1;
        pz1 = (ae_int32x4 *)((complex_fract32 *)pz0 + P);
        __Pragma("loop_count min=1");
        for (p = 0; p < P; p+=4)
        {
            px0 = (const ae_int32x4 *)(x + m*N);
            py0 = (const ae_int32x4 *)(yt + p*N);

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

                AE_L32X2X2_XP(y00, y10, py0, N*sizeof(complex_fract32));
                AE_L32X2X2_XP(y01, y11, py0, N*sizeof(complex_fract32));
                AE_L32X2X2_XP(y02, y12, py0, N*sizeof(complex_fract32));
                AE_L32X2X2_XP(y03, y13, py0, (-3*N+2)*(int)sizeof(complex_fract32));

                /* perform multiplications */
                AE_MULAAF2D32RA_HH_LL(Z00re, Z01re, y00, y01, x00, x00);  AE_MULASF2D32RA_HL_LH(Z00im, Z01im, y00, y01, x00, x00);
                AE_MULAAF2D32RA_HH_LL(Z00re, Z01re, y10, y11, x01, x01);  AE_MULASF2D32RA_HL_LH(Z00im, Z01im, y10, y11, x01, x01);
                AE_MULAAF2D32RA_HH_LL(Z10re, Z11re, y00, y01, x10, x10);  AE_MULASF2D32RA_HL_LH(Z10im, Z11im, y00, y01, x10, x10);
                AE_MULAAF2D32RA_HH_LL(Z10re, Z11re, y10, y11, x11, x11);  AE_MULASF2D32RA_HL_LH(Z10im, Z11im, y10, y11, x11, x11);
                AE_MULAAF2D32RA_HH_LL(Z02re, Z03re, y02, y03, x00, x00);  AE_MULASF2D32RA_HL_LH(Z02im, Z03im, y02, y03, x00, x00);
                AE_MULAAF2D32RA_HH_LL(Z02re, Z03re, y12, y13, x01, x01);  AE_MULASF2D32RA_HL_LH(Z02im, Z03im, y12, y13, x01, x01);
                AE_MULAAF2D32RA_HH_LL(Z12re, Z13re, y02, y03, x10, x10);  AE_MULASF2D32RA_HL_LH(Z12im, Z13im, y02, y03, x10, x10);
                AE_MULAAF2D32RA_HH_LL(Z12re, Z13re, y12, y13, x11, x11);  AE_MULASF2D32RA_HL_LH(Z12im, Z13im, y12, y13, x11, x11);
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
} /* cmtx_mpyt32x32_fast() */
#endif

size_t cmtx_mpyt32x32_fast_getScratchSize (int M, int N, int P)
{
    NASSERT((0==(M%32)) && (0==(N%32)) && (0==(P%32)));
    (void)M; (void)N; (void)P;
    return 0;
}
