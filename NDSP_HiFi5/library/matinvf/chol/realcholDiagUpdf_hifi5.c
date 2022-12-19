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
update N-th diagonal element
Input:
stride      width of matrix Z
Z[SZ]       N-th column in convolutions matrix Z
Input/output:
y           pointer to the begining of column in matrix R[L][SR] (N+1 elements is written)
Output:
D[SD]       reciprocals of main diagonal (pointer to the N-th element
--------------------------------------------------------*/
void realcholDiagUpdf(float32_t* y, float32_t* D, const float32_t* Z, int N)
{
    xtfloat * restrict pY = (xtfloat *)y;
    xtfloat * restrict pD = (xtfloat *)D;
    const xtfloat * restrict pZ = (xtfloat *)(Z); // points to diagonal element;
    xtfloat Z0, Y0, D0;
    xtfloat Yre;
    int m;
    XT_LSIP(Z0, pZ, 0);
    for (m = 0; m<N; m++)
    {
        XT_LSIP(Yre, pY, SZ_F32);
        XT_MSUB_S(Z0, Yre, Yre);
    }

    D0 = XT_RSQRT_S(Z0);
    XT_SSI(D0, pD, 0);
    Y0 = XT_MUL_S(D0, Z0);
    XT_SSI(Y0, pY, 0);
}
#endif
