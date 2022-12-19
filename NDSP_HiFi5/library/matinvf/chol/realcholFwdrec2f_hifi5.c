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
void realcholFwdrec2f(float32_t * y,
    const float32_t * R,
    const float32_t * D,
    const float32_t * Z, int stride)
{
    xtfloat * restrict pY = (xtfloat *)y;
    const xtfloat * restrict pZ = (const xtfloat *)Z;
    const xtfloat * restrict pR1 = (const xtfloat *)(R + 1);/* + ((n*(n+1))/2) */
    const xtfloat * restrict pD = (const xtfloat *)D;
    xtfloat Z0, Z1;
    xtfloat D0, D1;
    xtfloat Y0, Y1;
    xtfloat R0;

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
}
#endif
