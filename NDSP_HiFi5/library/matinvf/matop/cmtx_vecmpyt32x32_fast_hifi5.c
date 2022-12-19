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
#if (0)
DISCARD_FUN(void,cmtx_vecmpyt32x32_fast,(  complex_fract32* y,
                        const complex_fract32* x,
                        int N, int lsh ))
#else
void cmtx_vecmpyt32x32_fast ( complex_fract32* y,
                        const complex_fract32* x,
                        int N, int lsh )
{
    const ae_int32x4 * restrict px0;
    const ae_int32x4 * restrict py0;
    ae_int32x4 * restrict pz0;
    ae_int32x4 * restrict pz1;
    ae_int32x2 conj_reg;
    int m, n;

    NASSERT_ALIGN(x, HIFI_SIMD_WIDTH);
    NASSERT_ALIGN(y, HIFI_SIMD_WIDTH);
    NASSERT(0==(N%32));
    NASSERT(-31<=lsh && lsh<=31);

    if (N <= 0) return;

    conj_reg = AE_MOVDA32X2(1, -1);
    px0 = (const ae_int32x4 *)(x);
    __Pragma("loop_count min=1");
    for (m = 0; m < N; m+=4)
    {
        ae_int32x2 x0, x1, x2, x3;
        ae_int32x2 y0, y1;
        ae_int32x2 z00, z01, z02, z03, z10, z11, z12, z13;
        ae_int32x2 z20, z21, z22, z23, z30, z31, z32, z33;
        ae_int64 Z00re, Z00im, Z01re, Z01im, Z02re, Z02im, Z03re, Z03im;
        ae_int64 Z10re, Z10im, Z11re, Z11im, Z12re, Z12im, Z13re, Z13im;
        ae_int64 Z20re, Z20im, Z21re, Z21im, Z22re, Z22im, Z23re, Z23im;
        ae_int64 Z30re, Z30im, Z31re, Z31im              , Z33re, Z33im;

        py0 = px0;
        pz0 = (ae_int32x4 *)(y + m*N + m);
        pz1 = (ae_int32x4 *)((complex_float *)pz0 + N*4);

        /* load data */
        AE_L32X2X2_IP(x0, x1, px0, 2*sizeof(complex_fract32));
        AE_L32X2X2_IP(x2, x3, px0, 2*sizeof(complex_fract32));

        {
            py0 += 2;

            /* perform multiplications */
            AE_MULZAAF2D32RA_HH_LL(Z00re, Z01re, x0, x1, x0, x0);  AE_MULZASF2D32RA_HL_LH(Z00im, Z01im, x0, x1, x0, x0);
            AE_MULZAAF2D32RA_HH_LL(Z02re, Z03re, x2, x3, x0, x0);  AE_MULZASF2D32RA_HL_LH(Z02im, Z03im, x2, x3, x0, x0);
            AE_MULZAAF2D32RA_HH_LL(Z11re, Z12re, x1, x2, x1, x1);  AE_MULZASF2D32RA_HL_LH(Z11im, Z12im, x1, x2, x1, x1);
            AE_MULZAAF2D32RA_HH_LL(Z13re, Z22re, x3, x2, x1, x2);  AE_MULZASF2D32RA_HL_LH(Z13im, Z22im, x3, x2, x1, x2);
            AE_MULZAAF2D32RA_HH_LL(Z23re, Z33re, x3, x3, x2, x3);  AE_MULZASF2D32RA_HL_LH(Z23im, Z33im, x3, x3, x2, x3);

            z00 = AE_TRUNCA32X2F64S(Z00re, Z00im, 16+lsh);
            z01 = AE_TRUNCA32X2F64S(Z01re, Z01im, 16+lsh);
            z02 = AE_TRUNCA32X2F64S(Z02re, Z02im, 16+lsh);
            z03 = AE_TRUNCA32X2F64S(Z03re, Z03im, 16+lsh);
            z11 = AE_TRUNCA32X2F64S(Z11re, Z11im, 16+lsh);
            z12 = AE_TRUNCA32X2F64S(Z12re, Z12im, 16+lsh);
            z13 = AE_TRUNCA32X2F64S(Z13re, Z13im, 16+lsh);
            z22 = AE_TRUNCA32X2F64S(Z22re, Z22im, 16+lsh);
            z23 = AE_TRUNCA32X2F64S(Z23re, Z23im, 16+lsh);
            z33 = AE_TRUNCA32X2F64S(Z33re, Z33im, 16+lsh);
            AE_MUL2P32X4S(z10, z20, z01, z02, conj_reg, conj_reg);
            AE_MUL2P32X4S(z21, z30, z12, z03, conj_reg, conj_reg);
            AE_MUL2P32X4S(z31, z32, z13, z23, conj_reg, conj_reg);

            /* save values */
            AE_S32X2X2_XP(z00, z01, pz0, N*sizeof(complex_fract32));
            AE_S32X2X2_XP(z10, z11, pz0, N*sizeof(complex_fract32));
            AE_S32X2X2_XP(z20, z21, pz0, N*sizeof(complex_fract32));
            AE_S32X2X2_XP(z30, z31, pz0, (-3*N+2)*(int)sizeof(complex_fract32));
            AE_S32X2X2_XP(z02, z03, pz0, N*sizeof(complex_fract32));
            AE_S32X2X2_XP(z12, z13, pz0, N*sizeof(complex_fract32));
            AE_S32X2X2_XP(z22, z23, pz0, N*sizeof(complex_fract32));
            AE_S32X2X2_XP(z32, z33, pz0, (-3*N+2)*(int)sizeof(complex_fract32));
        }
        __Pragma("loop_count factor=2");
        for (n = 0; n < ((N-m-4) >> 1); n++)
        {
            /* load data */
            AE_L32X2X2_IP(y0, y1, py0, 2*sizeof(complex_fract32));

            /* perform multiplications */
            AE_MULZAAF2D32RA_HH_LL(Z00re, Z01re, y0, y1, x0, x0);  AE_MULZASF2D32RA_HL_LH(Z00im, Z01im, y0, y1, x0, x0);
            AE_MULZAAF2D32RA_HH_LL(Z10re, Z11re, y0, y1, x1, x1);  AE_MULZASF2D32RA_HL_LH(Z10im, Z11im, y0, y1, x1, x1);
            AE_MULZAAF2D32RA_HH_LL(Z20re, Z21re, y0, y1, x2, x2);  AE_MULZASF2D32RA_HL_LH(Z20im, Z21im, y0, y1, x2, x2);
            AE_MULZAAF2D32RA_HH_LL(Z30re, Z31re, y0, y1, x3, x3);  AE_MULZASF2D32RA_HL_LH(Z30im, Z31im, y0, y1, x3, x3);

            z00 = AE_TRUNCA32X2F64S(Z00re, Z00im, 16+lsh);
            z01 = AE_TRUNCA32X2F64S(Z01re, Z01im, 16+lsh);
            z10 = AE_TRUNCA32X2F64S(Z10re, Z10im, 16+lsh);
            z11 = AE_TRUNCA32X2F64S(Z11re, Z11im, 16+lsh);
            z20 = AE_TRUNCA32X2F64S(Z20re, Z20im, 16+lsh);
            z21 = AE_TRUNCA32X2F64S(Z21re, Z21im, 16+lsh);
            z30 = AE_TRUNCA32X2F64S(Z30re, Z30im, 16+lsh);
            z31 = AE_TRUNCA32X2F64S(Z31re, Z31im, 16+lsh);

            /* save values */
            AE_S32X2X2_XP(z00, z01, pz0, N*sizeof(complex_fract32));
            AE_S32X2X2_XP(z10, z11, pz0, N*sizeof(complex_fract32));
            AE_S32X2X2_XP(z20, z21, pz0, N*sizeof(complex_fract32));
            AE_S32X2X2_XP(z30, z31, pz0, (-3*N+2)*(int)sizeof(complex_fract32));
            AE_MUL2P32X4S(z00, z01, z00, z01, conj_reg, conj_reg);
            AE_MUL2P32X4S(z10, z11, z10, z11, conj_reg, conj_reg);
            AE_MUL2P32X4S(z20, z21, z20, z21, conj_reg, conj_reg);
            AE_MUL2P32X4S(z30, z31, z30, z31, conj_reg, conj_reg);
            AE_S32X2X2_IP(z00, z10, pz1, 2*sizeof(complex_fract32));
            AE_S32X2X2_XP(z20, z30, pz1, (N-2)*(int)sizeof(complex_fract32));
            AE_S32X2X2_IP(z01, z11, pz1, 2*sizeof(complex_fract32));
            AE_S32X2X2_XP(z21, z31, pz1, (N-2)*(int)sizeof(complex_fract32));
        }
    }
} /* cmtx_vecmpyt32x32_fast() */
#endif
