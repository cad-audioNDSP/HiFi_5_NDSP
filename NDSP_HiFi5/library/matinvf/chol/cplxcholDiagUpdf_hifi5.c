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

#include "math.h"/*DEL ME*/
/*
code optimized for HiFi4 with VFPU
*/

#if (HAVE_VFPU)
#define SZ_F32 (sizeof(float32_t))
#define SZ_CF32 (2*SZ_F32)

/* Load 32bit and replicate to xtfloatx2*/
#define _L32_SX2_IP(a,b,c) \
{ \
ae_int32x2 tmp; \
AE_L32_IP(tmp, castxcc(ae_int32, b), c); \
a = XT_AE_MOVXTFLOATX2_FROMINT32X2(tmp); }

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
void cplxcholDiagUpdf(complex_float* y, complex_float* D, const complex_float* Z, int N)
{
    xtfloatx2 * restrict pY = (xtfloatx2 *)y;
    xtfloatx2 * restrict pD = (xtfloatx2 *)D;
    const xtfloatx2 * restrict pZ = (xtfloatx2 *)(Z); // points to diagonal element;
    xtfloatx2 Z0, Zt, Y0, D0;
    xtfloatx2 Yre, Yim;
    Zt = XT_CONST_S(0);
    int m;
    _L32_SX2_IP(Z0, pZ, 0);
    for (m = 0; m<N; m++)
    {
        _L32_SX2_IP(Yre, pY, SZ_F32);
        _L32_SX2_IP(Yim, pY, SZ_F32);
        MSUB_SX2X2(Z0, Zt, Yre, Yim, Yre, Yim);
    }
    
    Z0 = XT_ADD_SX2(Z0, Zt);
    D0 = XT_RSQRT_SX2(Z0);
    XT_SSX2I(D0, pD, 0);
    Y0 = XT_MUL_SX2(D0, Z0);
    Y0 = XT_SEL32_HH_SX2(Y0, XT_CONST_S(0));
    XT_SSX2I(Y0, pY, 0);
}
#elif (HAVE_FPU)
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
void cplxcholDiagUpdf(complex_float* y, complex_float* D, const complex_float* Z, int N)
{
    int m;
    xtfloat * restrict pD = (xtfloat *)D;
    xtfloat * restrict pY = (xtfloat *)y;
    const xtfloat * restrict pZ = (const xtfloat *)Z;
    xtfloat B0 = XT_LSI(pZ,0);
    xtfloat D0, Y0;
    xtfloat Y_re, Y_im;
    for (m = 0; m<N; m++)
    {
        XT_LSIP(Y_re, pY, SZ_F32);
        XT_LSIP(Y_im, pY, SZ_F32);
        XT_MSUB_S(B0, Y_re, Y_re);
        XT_MSUB_S(B0, Y_im, Y_im);
    }
    D0 = XT_RSQRT_S(B0);
    XT_SSI(D0, pD, 0);
    XT_SSI(D0, pD, SZ_F32);
    Y0 = XT_MUL_S(D0, B0);
    XT_SSI(Y0, pY, 0);
    XT_SSI(XT_CONST_S(0), pY, SZ_F32);
}
#endif
