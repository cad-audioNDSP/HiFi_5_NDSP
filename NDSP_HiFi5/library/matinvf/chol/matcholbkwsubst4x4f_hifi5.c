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
void  matcholbkwsubst4x4f(void * pScr,
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

    const xtfloat * restrict pR3 = (const xtfloat*)(R + (4 * 3) / 2 + 0); //N*(N-1)/2
    const xtfloat * restrict pR2 = (const xtfloat*)(R + (4 * 3) / 2 + 1);
    const xtfloat * restrict pR1 = (const xtfloat*)(R + (4 * 3) / 2 + 2);
    const xtfloat * restrict pY = (const xtfloat*)(y + 3);
    const xtfloat * restrict pD = (const xtfloat*)(D + 3);
    xtfloat * restrict pX = (xtfloat*)(x + 3);
    xtfloat D0;
    xtfloat R1, R2, R3;
    xtfloat X0, X1, X2, X3;

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
    XT_LSXP(R2, pR2, -SZ_F32 * 3); // N-1
    XT_MSUB_S(X2, X0, R2);
    XT_LSIP(R2, pR2, 0);
    XT_MSUB_S(X2, X1, R2);
    X2 = XT_MUL_S(X2, D0);
    XT_SSXP(X2, pX, -SZ_F32);

    // X3
    XT_LSXP(D0, pD, -SZ_F32);
    XT_LSXP(X3, pY, -SZ_F32);
    XT_LSXP(R3, pR3, -SZ_F32 * 3);
    XT_MSUB_S(X3, X0, R3);
    XT_LSXP(R3, pR3, -SZ_F32 * 2);
    XT_MSUB_S(X3, X1, R3);
    XT_LSIP(R3, pR3, 0);
    XT_MSUB_S(X3, X2, R3);
    X3 = XT_MUL_S(X3, D0);
    XT_SSXP(X3, pX, -SZ_F32);
}
/* scratch allocation functions */
size_t  matcholbkwsubst4x4f_getScratchSize()
{
    return 0;
}

#else
DISCARD_FUN(void, matcholbkwsubst4x4f, (void * pScr,
    float32_t * x,
    const float32_t * R,
    const float32_t * D,
    const float32_t * y))

    size_t  matcholbkwsubst4x4f_getScratchSize()
{
    return 0;
}
#endif
