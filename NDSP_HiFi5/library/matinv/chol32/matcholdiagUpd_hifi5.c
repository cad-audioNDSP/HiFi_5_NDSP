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
#include "chol32x32_common.h"

/*--------------------------------------------------------
   update N-th diagonal element
   Input:
   Z[L][N+1][C]  convolutions in N-th column
   C             1 for real, 2 for complex
   Input/output:
   y             pointer to the begining of column in matrix 
                 R[L][SR] (N+1 elements is written)
   Output:
   D[L][SD]      reciprocals of main diagonal (pointer to 
                 the N-th element
--------------------------------------------------------*/
void matcholdiagUpd(int32_t* y, int32_t* D, const int64_t* Z, int N)
{
	ae_int64 a, A_re, t;
	ae_int32 d;
	int d_exp;
	int m;
	ae_ep ep;

	t = AE_L64_I((const ae_int64*)Z, 0);
	AE_SLAI72(ep, a, t, 4);
	for (m = 0; m < N; m++)
	{
		ae_int32x2 r;
		AE_L32_IP(r, castxcc(ae_int32, y), sizeof(ae_int32));
        AE_MULS32EP_HH(ep, a, r, r);
	}
	A_re = AE_SRAI72(ep, a, 4);
	invSqrt_hifi(D, (int64_t)A_re);
	d = D[0];
	d_exp = D[1];
	AE_MUL32USEP_LH(ep, a, AE_MOVINT32X2_FROMINT64(A_re), d);
	t = AE_SRAI72(ep, a, 31);
	AE_MULAF32S_HH(t, AE_MOVINT32X2_FROMINT64(A_re), d);
	A_re = t;
	d = AE_TRUNCA32F64S(A_re, d_exp + 5);
	AE_S32_L_I(d, (ae_int32*)y, 0);
}
/* diagUpd() */
