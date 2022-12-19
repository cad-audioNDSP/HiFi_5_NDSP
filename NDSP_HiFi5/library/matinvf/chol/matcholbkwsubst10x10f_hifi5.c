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
#include "NatureDSP_types.h"
#include "common.h"
#include "NatureDSP_Signal_matinv.h"
#include <math.h>
#include "common_fpu.h"
#include "cholnf_common.h"
/*
code optimized for HiFi4 with VFPU
*/

#if (HAVE_VFPU || HAVE_FPU)
#define SZ_F32 (int)(sizeof(float32_t))
/*-------------------------------------------------------------------------
Cholesky Backward Substitution for Pseudo-inversion
These functions make backward substitution stage of pseudo-inversion. They
use Cholesky decomposition of original matrices and results of forward
substitution.

Precision:
f     single precision floating point
32x32 32-bit inputs/outputs

Input:
R[((N+1)*N)/2]
            Cholesky upper-triangle matrix R. For fixed point API, the 
            representation is Q(qR)
y[N*P]      Results of forward substitution stage. For fixed point API, 
            the representation is Q(qY)
D[N]        sequence of reciprocals of main diagonal R. NOTE: for the 
            fixed point API, these data are stored internally in special 
            format with separate mantissa and exponent for better accuracy 
            and dynamic range control. So, even for the real data, they 
            stored as pairs of 2 integers and packed to the complex_fract32 
            format
qXYR        combination of fixed point representation (matrices R, x and y) 
            qXYR=qX-qY+qR (for fixed point API only)
Output:
x[N*P]      Decision matrix x. For fixed point API, the representation is Q(qX)

N = 4, 6, 8, 10
P = 1

Restrictions:
All matrices and scratch should not overlap and be aligned on
16-bytes boundary
---------------------------------------------------------------------------*/
void  matcholbkwsubst10x10f(void * pScr,
    float32_t * x,
    const float32_t * R,
    const float32_t * D,
    const float32_t * y)
{
    NASSERT(x);
    NASSERT(R);
    NASSERT(D);
    NASSERT(y);
    NASSERT_ALIGN(x, HIFI_SIMD_WIDTH);
    NASSERT_ALIGN(R, HIFI_SIMD_WIDTH);
    NASSERT_ALIGN(D, HIFI_SIMD_WIDTH);
    NASSERT_ALIGN(y, HIFI_SIMD_WIDTH);
    //NASSERT_ALIGN(pScr, HIFI_SIMD_WIDTH);

    const xtfloat * restrict pR9 = (const xtfloat*)(R + (10 * 9) / 2 + 0); //N*(N-1)/2
    const xtfloat * restrict pR8 = (const xtfloat*)(R + (10 * 9) / 2 + 1);
    const xtfloat * restrict pR7 = (const xtfloat*)(R + (10 * 9) / 2 + 2);
    const xtfloat * restrict pR6 = (const xtfloat*)(R + (10 * 9) / 2 + 3);
    const xtfloat * restrict pR5 = (const xtfloat*)(R + (10 * 9) / 2 + 4);
    const xtfloat * restrict pR4 = (const xtfloat*)(R + (10 * 9) / 2 + 5);
    const xtfloat * restrict pR3 = (const xtfloat*)(R + (10 * 9) / 2 + 6);
    const xtfloat * restrict pR2 = (const xtfloat*)(R + (10 * 9) / 2 + 7);
    const xtfloat * restrict pR1 = (const xtfloat*)(R + (10 * 9) / 2 + 8);
    const xtfloat * restrict pY = (const xtfloat*)(y + 9);
    const xtfloat * restrict pD = (const xtfloat*)(D + 9);
    xtfloat * restrict pX = (xtfloat*)(x + 9);
    xtfloat D0;
    xtfloat R1, R2, R3, R4, R5, R6, R7, R8, R9;
    xtfloat X0, X1, X2, X3, X4, X5, X6, X7, X8, X9;

    // X0
    XT_LSXP(D0, pD, -SZ_F32);
    XT_LSXP(X0, pY, -SZ_F32);
    X0 = XT_MUL_S(X0, D0);
    XT_SSXP(X0, pX, -SZ_F32);

    // X1
    XT_LSXP(D0, pD, -SZ_F32);
    XT_LSXP(X1, pY, -SZ_F32);
    XT_LSIP(R1, pR1, 0);
    XT_MSUB_S(X1, X0, R1);
    X1 = XT_MUL_S(X1, D0);
    XT_SSXP(X1, pX, -SZ_F32);

    // X2
    XT_LSXP(D0, pD, -SZ_F32);
    XT_LSXP(X2, pY, -SZ_F32);
    XT_LSXP(R2, pR2, -SZ_F32 * 9); // N-1
    XT_MSUB_S(X2, X0, R2);
    XT_LSIP(R2, pR2, 0);
    XT_MSUB_S(X2, X1, R2);
    X2 = XT_MUL_S(X2, D0);
    XT_SSXP(X2, pX, -SZ_F32);

    // X3
    XT_LSXP(D0, pD, -SZ_F32);
    XT_LSXP(X3, pY, -SZ_F32);
    XT_LSXP(R3, pR3, -SZ_F32 * 9);
    XT_MSUB_S(X3, X0, R3);
    XT_LSXP(R3, pR3, -SZ_F32 * 8);
    XT_MSUB_S(X3, X1, R3);
    XT_LSIP(R3, pR3, 0);
    XT_MSUB_S(X3, X2, R3);
    X3 = XT_MUL_S(X3, D0);
    XT_SSXP(X3, pX, -SZ_F32);

    // X4
    XT_LSXP(D0, pD, -SZ_F32);
    XT_LSXP(X4, pY, -SZ_F32);
    XT_LSXP(R4, pR4, -SZ_F32 * 9);
    XT_MSUB_S(X4, X0, R4);
    XT_LSXP(R4, pR4, -SZ_F32 * 8);
    XT_MSUB_S(X4, X1, R4);
    XT_LSXP(R4, pR4, -SZ_F32 * 7);
    XT_MSUB_S(X4, X2, R4);
    XT_LSIP(R4, pR4, 0);
    XT_MSUB_S(X4, X3, R4);
    X4 = XT_MUL_S(X4, D0);
    XT_SSXP(X4, pX, -SZ_F32);

    // X5
    XT_LSXP(D0, pD, -SZ_F32);
    XT_LSXP(X5, pY, -SZ_F32);
    XT_LSXP(R5, pR5, -SZ_F32 * 9);
    XT_MSUB_S(X5, X0, R5);
    XT_LSXP(R5, pR5, -SZ_F32 * 8);
    XT_MSUB_S(X5, X1, R5);
    XT_LSXP(R5, pR5, -SZ_F32 * 7);
    XT_MSUB_S(X5, X2, R5);
    XT_LSXP(R5, pR5, -SZ_F32 * 6);
    XT_MSUB_S(X5, X3, R5);
    XT_LSIP(R5, pR5, 0);
    XT_MSUB_S(X5, X4, R5);
    X5 = XT_MUL_S(X5, D0);
    XT_SSXP(X5, pX, -SZ_F32);

    // X6
    XT_LSXP(D0, pD, -SZ_F32);
    XT_LSXP(X6, pY, -SZ_F32);
    XT_LSXP(R6, pR6, -SZ_F32 * 9);
    XT_MSUB_S(X6, X0, R6);
    XT_LSXP(R6, pR6, -SZ_F32 * 8);
    XT_MSUB_S(X6, X1, R6);
    XT_LSXP(R6, pR6, -SZ_F32 * 7);
    XT_MSUB_S(X6, X2, R6);
    XT_LSXP(R6, pR6, -SZ_F32 * 6);
    XT_MSUB_S(X6, X3, R6);
    XT_LSXP(R6, pR6, -SZ_F32 * 5);
    XT_MSUB_S(X6, X4, R6);
    XT_LSIP(R6, pR6, 0);
    XT_MSUB_S(X6, X5, R6);
    X6 = XT_MUL_S(X6, D0);
    XT_SSXP(X6, pX, -SZ_F32);

    // X7
    XT_LSXP(D0, pD, -SZ_F32);
    XT_LSXP(X7, pY, -SZ_F32);
    XT_LSXP(R7, pR7, -SZ_F32 * 9);
    XT_MSUB_S(X7, X0, R7);
    XT_LSXP(R7, pR7, -SZ_F32 * 8);
    XT_MSUB_S(X7, X1, R7);
    XT_LSXP(R7, pR7, -SZ_F32 * 7);
    XT_MSUB_S(X7, X2, R7);
    XT_LSXP(R7, pR7, -SZ_F32 * 6);
    XT_MSUB_S(X7, X3, R7);
    XT_LSXP(R7, pR7, -SZ_F32 * 5);
    XT_MSUB_S(X7, X4, R7);
    XT_LSXP(R7, pR7, -SZ_F32 * 4);
    XT_MSUB_S(X7, X5, R7);
    XT_LSIP(R7, pR7, 0);
    XT_MSUB_S(X7, X6, R7);
    X7 = XT_MUL_S(X7, D0);
    XT_SSXP(X7, pX, -SZ_F32);

    // X8
    XT_LSXP(D0, pD, -SZ_F32);
    XT_LSXP(X8, pY, -SZ_F32);
    XT_LSXP(R8, pR8, -SZ_F32 * 9);
    XT_MSUB_S(X8, X0, R8);
    XT_LSXP(R8, pR8, -SZ_F32 * 8);
    XT_MSUB_S(X8, X1, R8);
    XT_LSXP(R8, pR8, -SZ_F32 * 7);
    XT_MSUB_S(X8, X2, R8);
    XT_LSXP(R8, pR8, -SZ_F32 * 6);
    XT_MSUB_S(X8, X3, R8);
    XT_LSXP(R8, pR8, -SZ_F32 * 5);
    XT_MSUB_S(X8, X4, R8);
    XT_LSXP(R8, pR8, -SZ_F32 * 4);
    XT_MSUB_S(X8, X5, R8);
    XT_LSXP(R8, pR8, -SZ_F32 * 3);
    XT_MSUB_S(X8, X6, R8);
    XT_LSIP(R8, pR8, 0);
    XT_MSUB_S(X8, X7, R8);
    X8 = XT_MUL_S(X8, D0);
    XT_SSXP(X8, pX, -SZ_F32);

    // X9
    XT_LSXP(D0, pD, -SZ_F32);
    XT_LSXP(X9, pY, -SZ_F32);
    XT_LSXP(R9, pR9, -SZ_F32 * 9);
    XT_MSUB_S(X9, X0, R9);
    XT_LSXP(R9, pR9, -SZ_F32 * 8);
    XT_MSUB_S(X9, X1, R9);
    XT_LSXP(R9, pR9, -SZ_F32 * 7);
    XT_MSUB_S(X9, X2, R9);
    XT_LSXP(R9, pR9, -SZ_F32 * 6);
    XT_MSUB_S(X9, X3, R9);
    XT_LSXP(R9, pR9, -SZ_F32 * 5);
    XT_MSUB_S(X9, X4, R9);
    XT_LSXP(R9, pR9, -SZ_F32 * 4);
    XT_MSUB_S(X9, X5, R9);
    XT_LSXP(R9, pR9, -SZ_F32 * 3);
    XT_MSUB_S(X9, X6, R9);
    XT_LSXP(R9, pR9, -SZ_F32 * 2);
    XT_MSUB_S(X9, X7, R9);
    XT_LSIP(R9, pR9, 0);
    XT_MSUB_S(X9, X8, R9);
    X9 = XT_MUL_S(X9, D0);
    XT_SSXP(X9, pX, -SZ_F32);
}
/* scratch allocation functions */
size_t  matcholbkwsubst10x10f_getScratchSize()
{
    return 0;
}

#else
DISCARD_FUN(void, matcholbkwsubst10x10f, (void * pScr,
    float32_t * x,
    const float32_t * R,
    const float32_t * D,
    const float32_t * y))

    size_t  matcholbkwsubst10x10f_getScratchSize()
{
    return 0;
}
#endif
