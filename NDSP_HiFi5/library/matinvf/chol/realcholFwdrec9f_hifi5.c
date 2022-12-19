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
#include "common_fpu.h"
#include "cholnf_common.h"
/*
code optimized for HiFi4 with VFPU
*/

#if (HAVE_VFPU || HAVE_FPU)
#define SZ_F32 (sizeof(float32_t))
/*--------------------------------------------------------
Forward recursion P=1
Input:
R[((N+1)*N)/2]  upper triangular matrix R
stride          width of matrix A'*B
Z[N*1]          column in matrix A'*B
D[N]            reciprocals of main diagonal
Output:
y[N*1]		Decision matrix y
N = 1, 2, 3, 4, 5, 6, 7, 8, 9, 10
--------------------------------------------------------*/
void realcholFwdrec9f(float32_t * y,
    const float32_t * R,
    const float32_t * D,
    const float32_t * Z, int stride)
{
    xtfloat * restrict pY = (xtfloat *)y;
    const xtfloat * restrict pZ = (const xtfloat *)Z;
    const xtfloat * restrict pR1 = (const xtfloat *)(R + 1);/* + ((n*(n+1))/2) */
    const xtfloat * restrict pR2 = (const xtfloat *)(R + 3);
    const xtfloat * restrict pR3 = (const xtfloat *)(R + 6);
    const xtfloat * restrict pR4 = (const xtfloat *)(R + 10);
    const xtfloat * restrict pR5 = (const xtfloat *)(R + 15);
    const xtfloat * restrict pR6 = (const xtfloat *)(R + 21);
    const xtfloat * restrict pR7 = (const xtfloat *)(R + 28);
    const xtfloat * restrict pR8 = (const xtfloat *)(R + 36);
    const xtfloat * restrict pD = (const xtfloat *)D;
    xtfloat Z0, Z1, Z2, Z3, Z4, Z5, Z6, Z7, Z8;
    xtfloat D0, D1, D2, D3, D4, D5, D6, D7, D8;
    xtfloat Y0, Y1, Y2, Y3, Y4, Y5, Y6, Y7, Y8;
    xtfloat R0, R1, R2, R3, R4, R5, R6, R7;

    // Y0
    XT_LSXP(Z0, pZ, stride*SZ_F32);
    XT_LSIP(D0, pD, SZ_F32);
    Y0 = XT_MUL_S(Z0, D0);
    XT_SSIP(Y0, pY, SZ_F32);

    // Y1
    XT_LSXP(Z1, pZ, stride*SZ_F32);
    XT_LSIP(D1, pD, SZ_F32);
    XT_LSIP(R0, pR1, 0);
    XT_MSUB_S(Z1, R0, Y0);
    Y1 = XT_MUL_S(Z1, D1);
    XT_SSIP(Y1, pY, SZ_F32);

    // Y2
    XT_LSXP(Z2, pZ, stride*SZ_F32);
    XT_LSIP(D2, pD, SZ_F32);
    XT_LSIP(R0, pR2, SZ_F32);
    XT_LSIP(R1, pR2, 0);
    XT_MSUB_S(Z2, R0, Y0);
    XT_MSUB_S(Z2, R1, Y1);
    Y2 = XT_MUL_S(Z2, D2);
    XT_SSIP(Y2, pY, SZ_F32);

    // Y3
    XT_LSXP(Z3, pZ, stride*SZ_F32);
    XT_LSIP(D3, pD, SZ_F32);
    XT_LSIP(R0, pR3, SZ_F32);
    XT_LSIP(R1, pR3, SZ_F32);
    XT_LSIP(R2, pR3, 0);
    XT_MSUB_S(Z3, R0, Y0);
    XT_MSUB_S(Z3, R1, Y1);
    XT_MSUB_S(Z3, R2, Y2);
    Y3 = XT_MUL_S(Z3, D3);
    XT_SSIP(Y3, pY, SZ_F32);

    // Y4
    XT_LSXP(Z4, pZ, stride*SZ_F32);
    XT_LSIP(D4, pD, SZ_F32);
    XT_LSIP(R0, pR4, SZ_F32);
    XT_LSIP(R1, pR4, SZ_F32);
    XT_LSIP(R2, pR4, SZ_F32);
    XT_LSIP(R3, pR4, 0);
    XT_MSUB_S(Z4, R0, Y0);
    XT_MSUB_S(Z4, R1, Y1);
    XT_MSUB_S(Z4, R2, Y2);
    XT_MSUB_S(Z4, R3, Y3);
    Y4 = XT_MUL_S(Z4, D4);
    XT_SSIP(Y4, pY, SZ_F32);

    // Y5
    XT_LSXP(Z5, pZ, stride*SZ_F32);
    XT_LSIP(D5, pD, SZ_F32);
    XT_LSIP(R0, pR5, SZ_F32);
    XT_LSIP(R1, pR5, SZ_F32);
    XT_LSIP(R2, pR5, SZ_F32);
    XT_LSIP(R3, pR5, SZ_F32);
    XT_LSIP(R4, pR5, 0);
    XT_MSUB_S(Z5, R0, Y0);
    XT_MSUB_S(Z5, R1, Y1);
    XT_MSUB_S(Z5, R2, Y2);
    XT_MSUB_S(Z5, R3, Y3);
    XT_MSUB_S(Z5, R4, Y4);
    Y5 = XT_MUL_S(Z5, D5);
    XT_SSIP(Y5, pY, SZ_F32);

    // Y6
    XT_LSXP(Z6, pZ, stride*SZ_F32);
    XT_LSIP(D6, pD, SZ_F32);
    XT_LSIP(R0, pR6, SZ_F32);
    XT_LSIP(R1, pR6, SZ_F32);
    XT_LSIP(R2, pR6, SZ_F32);
    XT_LSIP(R3, pR6, SZ_F32);
    XT_LSIP(R4, pR6, SZ_F32);
    XT_LSIP(R5, pR6, 0);
    XT_MSUB_S(Z6, R0, Y0);
    XT_MSUB_S(Z6, R1, Y1);
    XT_MSUB_S(Z6, R2, Y2);
    XT_MSUB_S(Z6, R3, Y3);
    XT_MSUB_S(Z6, R4, Y4);
    XT_MSUB_S(Z6, R5, Y5);
    Y6 = XT_MUL_S(Z6, D6);
    XT_SSIP(Y6, pY, SZ_F32);

    // Y7
    XT_LSXP(Z7, pZ, stride*SZ_F32);
    XT_LSIP(D7, pD, SZ_F32);
    XT_LSIP(R0, pR7, SZ_F32);
    XT_LSIP(R1, pR7, SZ_F32);
    XT_LSIP(R2, pR7, SZ_F32);
    XT_LSIP(R3, pR7, SZ_F32);
    XT_LSIP(R4, pR7, SZ_F32);
    XT_LSIP(R5, pR7, SZ_F32);
    XT_LSIP(R6, pR7, 0);
    XT_MSUB_S(Z7, R0, Y0);
    XT_MSUB_S(Z7, R1, Y1);
    XT_MSUB_S(Z7, R2, Y2);
    XT_MSUB_S(Z7, R3, Y3);
    XT_MSUB_S(Z7, R4, Y4);
    XT_MSUB_S(Z7, R5, Y5);
    XT_MSUB_S(Z7, R6, Y6);
    Y7 = XT_MUL_S(Z7, D7);
    XT_SSIP(Y7, pY, SZ_F32);

    // Y8
    XT_LSXP(Z8, pZ, stride*SZ_F32);
    XT_LSIP(D8, pD, SZ_F32);
    XT_LSIP(R0, pR8, SZ_F32);
    XT_LSIP(R1, pR8, SZ_F32);
    XT_LSIP(R2, pR8, SZ_F32);
    XT_LSIP(R3, pR8, SZ_F32);
    XT_LSIP(R4, pR8, SZ_F32);
    XT_LSIP(R5, pR8, SZ_F32);
    XT_LSIP(R6, pR8, SZ_F32);
    XT_LSIP(R7, pR8, 0);
    XT_MSUB_S(Z8, R0, Y0);
    XT_MSUB_S(Z8, R1, Y1);
    XT_MSUB_S(Z8, R2, Y2);
    XT_MSUB_S(Z8, R3, Y3);
    XT_MSUB_S(Z8, R4, Y4);
    XT_MSUB_S(Z8, R5, Y5);
    XT_MSUB_S(Z8, R6, Y6);
    XT_MSUB_S(Z8, R7, Y7);
    Y8 = XT_MUL_S(Z8, D8);
    XT_SSIP(Y8, pY, SZ_F32);
}
#endif
