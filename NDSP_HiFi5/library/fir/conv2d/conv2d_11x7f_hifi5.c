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
  NatureDSP Signal Processing Library. FIR part
  2D Convolution  
  IntegrIT, 2006-2020
*/

/* Portable data types. */
#include "NatureDSP_types.h"
/* Signal Processing Library API. */
#include "NatureDSP_Signal_fir.h"
#include "common_fpu.h"


/*-------------------------------------------------------------------------
  2D convolution
  Functions compute the two-dimensional convolution of input matrix x[M][N]
  and y[P][Q] and store the result in matrix z[M+P-1][N+Q-1]
  Additional parameter rsh allows to contro l fixed point representation of 
  output data.
  Two versions of functions available: 
  - generic version with _gen_ suffix. 
    These functions work with arbitrary arguments.
  - fast version with no _gen_ suffix. 
    These functions expose some additional restrictions on argument


  Precision: 
  8x8      8-bit coefficients, 8-bit data, 8-bit output, Q7
  8x16     8-bit coefficients Q7, 16-bit data, 16-bit output, Q15
  16x16    16-bit coefficients, 16-bit data, 16-bit output, Q15
  f        single precision floating point data
  fp16     half precision floating point data

  Input:
  x[M][N]   input data Q15, Q7, floating point
  y[P][Q]   input data Q15, Q7, floating point
  M         number of rows in the matrix x
  N         number of columns in the matrix x
  P         number of rows in the matrix y
  Q         number of columns in the matrix y
  rsh       additional right shift (for fixed point API only)

  Output:
  z	[M+P-1][N+Q-1] output data, Q(7-rsh), Q(15-rsh)

  Temporary:
  pScr     scratch data. Should have size at least as requested by 
           corresponding scratch allocation function

  Restrictions:
  For regular routines:
  x,y,z        should not overlap
  pScr         aligned on a 16-bytes boundary
  P, Q	       >0

  For fast routines:
  x,y,z        should not overlap
  x,y,z,pScr   aligned on a 16-bytes boundary
  P, Q	       >0 and multiplies of 8
-------------------------------------------------------------------------*/

#if (HAVE_FPU==0 && HAVE_VFPU==0)
DISCARD_FUN(void,conv2d_11x7f,(void* pScr, float32_t *z, const float32_t * x, const float32_t * y, int P, int Q))
size_t conv2d_11x7f_getScratchSize (int P, int Q) 
{ 
    (void)P,(void)Q;
    return 0;
}
#elif HAVE_VFPU

void conv2d_11x7f     (void* pScr, float32_t *z, const float32_t * x, const float32_t * y, int P, int Q)
{
#   define M 11
#   define N 7
	int i, j, m, y_start, y_end, m_end, m_start;

	const xtfloatx2 * restrict pX;
	const xtfloatx4 * restrict pY;
	const xtfloatx4 * restrict pW;
	xtfloatx4 * restrict pWwr;
	xtfloatx4 * restrict pZ;
	ae_valign aX;
	ae_valignx2 alZ;
    xtfloatx2 tmp;
	xtfloatx2 Y01, Y23, Y45, Y67, Y89;
	xtfloatx2 W0z, W10, W21, W32, W43, W54, W65, Wz6;

	NASSERT(x);
	NASSERT(y);
	NASSERT(z);
	NASSERT(pScr);
    NASSERT_ALIGN(x,XCHAL_DATA_WIDTH);
    NASSERT_ALIGN(y,XCHAL_DATA_WIDTH);
    NASSERT_ALIGN(z,XCHAL_DATA_WIDTH);
    NASSERT_ALIGN(pScr,XCHAL_DATA_WIDTH);
    NASSERT(P>=0 && P%8==0);
	NASSERT(Q >= 0 && Q % 8 == 0);

	if (P <= 0 || Q <= 0) return;

	/* Store coefficients in the next order:
	*     0  x[6]
	*   x[5] x[4]
    *   x[6] x[5]
	*   x[3] x[2]
    *   x[4] x[3]
	*   x[1] x[0]
	*   x[2] x[1]
    *   X[0]   0
	* Start from the last row
	* Double up data to reduce stalls
	*/
#define WLEN 4
	pWwr = (xtfloatx4 *)pScr;
	pX = (xtfloatx2 *)(x + (M - 1) * N + 6);
	for (i = M - 1; i >= 0; i--)
	{
		ae_int32x2 v_i;							

        AE_L32_XP(v_i, castxcc(const ae_int32, pX), -(int)sizeof(float32_t));
        Wz6 = XT_AE_MOVXTFLOATX2_FROMINT32X2(v_i);//6

		Wz6 = XT_SEL32_LH_SX2((xtfloatx2)(0.0f), Wz6); // z 6 

		aX = AE_LA64_PP(pX);

		XT_LASX2RIP(W54, aX, pX);	// 5 4
		AE_SSX2X2_IP(Wz6, W54, pWwr, sizeof(xtfloatx4));

        W65 = XT_SEL32_LH_SX2(Wz6, W54);

        XT_LASX2RIP(W32, aX, pX); // 3 2
        AE_SSX2X2_IP(W65, W32, pWwr, sizeof(xtfloatx4));

        W43 = XT_SEL32_LH_SX2(W54, W32);

        XT_LASX2RIP(W10, aX, pX); // 1 0
        AE_SSX2X2_IP(W43, W10, pWwr, sizeof(xtfloatx4));

        W21 = XT_SEL32_LH_SX2(W32, W10);
        W0z = XT_SEL32_LH_SX2(W10, (xtfloatx2)(0.0f));

        AE_SSX2X2_IP(W21, W0z, pWwr, sizeof(xtfloatx4));

	}
	/*
	* Processing of convolution
	*/
	pZ = (xtfloatx4 *)(z);
	alZ = AE_ZALIGN128();
	for (i = 0; i < M + P - 1; i++)
	{
		y_start = XT_MAX(i + 1 - M, 0);
		y_end = XT_MIN(i + 1, P);
		m_end = XT_MIN(i + 1, M);
		m_start = m_end - (y_end - y_start);
		
		WAE_CBEGIN0((uintptr_t)(y + y_start*Q+4));
		WAE_CEND0((uintptr_t)(y + y_end*Q));

		WAE_CBEGIN1((uintptr_t)((float32_t *)pScr + (M - m_end) * WLEN*4 ));
		WAE_CEND1((uintptr_t)((float32_t *)pScr + (M - m_start) * WLEN*4 ));

		pW = (const xtfloatx4 *)pScr + (M - m_end) * WLEN;

		/* First N + 1 samples of the i-th row */
		{
			xtfloatx2 S0, S2, S4, S6;
			xtfloatx2 S1, S3, S5, S7;
			xtfloatx2 S2_0, S3_0;
			xtfloatx2 S4_0, S5_0, S4_1, S5_1;
			xtfloatx2 S6_0, S7_0, S6_1, S7_1, S6_2, S7_2;

			pY = (const xtfloatx4 *)(y + y_start * Q);

			S0 = S2 = S4 = S5 = (xtfloatx2)(0.0f);
			S1 = S3 = S6 = S7 = (xtfloatx2)(0.0f);
			S2_0 = S3_0 = (xtfloatx2)(0.0f);
			S4_0 = S5_0 = S4_1 = S5_1 = (xtfloatx2)(0.0f);
			S6_0 = S7_0 = S6_1 = S7_1 = S6_2 = S7_2 = (xtfloatx2)(0.0f);

			__Pragma("loop_count min=1, max=11");
			for (m = m_start; m < m_end; m++)
			{

				AE_LSX2X2_IP(Y01, Y23, pY, sizeof(xtfloatx4));
				AE_LSX2X2_XP(Y45, Y67, pY, (Q-4) * sizeof(float32_t));
                
				AE_LSX2X2_IP(Wz6, W54, pW, sizeof(xtfloatx4));
				AE_LSX2X2_IP(W65, W32, pW, sizeof(xtfloatx4));
				AE_LSX2X2_IP(W43, W10, pW, sizeof(xtfloatx4));
				AE_LSX2X2_XC1(W21, W0z, pW, sizeof(xtfloatx4));


                MADD_SX2X2(S6, S7, W65, Wz6, Y01, Y01);
                MADD_SX2X2(S4, S5, W43, W54, Y01, Y01);
                MADD_SX2X2(S2, S3, W21, W32, Y01, Y01);
                MADD_SX2X2(S0, S1, W0z, W10, Y01, Y01);

                MADD_SX2X2(S6_0, S7_0, W43, W54, Y23, Y23);
                MADD_SX2X2(S4_0, S5_0, W21, W32, Y23, Y23);
                MADD_SX2X2(S2_0, S3_0, W0z, W10, Y23, Y23);

                MADD_SX2X2(S6_1, S7_1, W21, W32, Y45, Y45);
                MADD_SX2X2(S4_1, S5_1, W0z, W10, Y45, Y45);

                MADD_SX2X2(S6_2, S7_2, W0z, W10, Y67, Y67);

			}
            
            ADD_SX2X2(S2, S3, S2, S3, S2_0, S3_0);

            ADD_SX2X2(S4, S5, S4, S5, S4_0, S5_0);
            ADD_SX2X2(S4, S5, S4, S5, S4_1, S5_1);

            ADD_SX2X2(S6, S7, S6, S7, S6_0, S7_0);
            ADD_SX2X2(S6, S7, S6, S7, S6_1, S7_1);
            ADD_SX2X2(S6, S7, S6, S7, S6_2, S7_2);


            tmp = CONST_SX2(1);

            MADDMUX_SX2X2 (S0, S4, tmp, tmp, S0,S4, 5);
            MADDMUX_SX2X2 (S1, S5, tmp, tmp, S1,S5, 5);
            MADDMUX_SX2X2 (S2, S6, tmp, tmp, S2,S6, 5);
            MADDMUX_SX2X2 (S3, S7, tmp, tmp, S3,S7, 5);

			S0 = XT_SEL32_HL_SX2(S0, S1);
			S2 = XT_SEL32_HL_SX2(S2, S3);
			S4 = XT_SEL32_HL_SX2(S4, S5);
			S6 = XT_SEL32_HL_SX2(S6, S7);


			AE_SASX2X2_IP(S0, S2, alZ, pZ);
			AE_SASX2X2_IP(S4, S6, alZ, pZ);

		}

		/* Next samples */
		pY = (const xtfloatx4 *)(y + y_start * Q + 2);

		for (j = N + 1; j < Q; j += 4)
		{
			xtfloatx2 S0_0, S1_0, S2_0, S3_0;
			xtfloatx2 S0_1, S1_1, S2_1, S3_1;
			xtfloatx2 S0_2, S1_2, S2_2, S3_2;
			xtfloatx2 S0_3, S1_3, S2_3, S3_3;


            {
                AE_LSX2X2_IP(Wz6, W54, pW, sizeof(xtfloatx4));
                AE_LSX2X2_IP(W65, W32, pW, sizeof(xtfloatx4));
                AE_LSX2X2_IP(W43, W10, pW, sizeof(xtfloatx4));
                AE_LSX2X2_XC1(W21, W0z, pW, sizeof(xtfloatx4));

                AE_LSX2IP(Y01, castxcc(xtfloatx2, pY), sizeof(xtfloatx2));
                AE_LSX2X2_IP(Y23, Y45, pY, sizeof(xtfloatx4));
                AE_LSX2X2_XC(Y67, Y89, pY, (Q - 6) * sizeof(float32_t));


                MUL_SX2X2(S0_0, S1_0, W65, Wz6, Y01, Y01);
                MUL_SX2X2(S2_0, S3_0, W65, Wz6, Y23, Y23);
                MUL_SX2X2(S0_1, S1_1, W43, W54, Y23, Y23);
                MUL_SX2X2(S2_1, S3_1, W43, W54, Y45, Y45);
                MUL_SX2X2(S0_2, S1_2, W21, W32, Y45, Y45);
                MUL_SX2X2(S2_2, S3_2, W21, W32, Y67, Y67);
                MUL_SX2X2(S0_3, S1_3, W0z, W10, Y67, Y67);
                MUL_SX2X2(S2_3, S3_3, W0z, W10, Y89, Y89);
            }
			__Pragma("loop_count max=10");
			for (m = m_start; m < m_end-1; m++)
			{
                AE_LSX2X2_IP(Wz6, W54, pW, sizeof(xtfloatx4));
                AE_LSX2X2_IP(W65, W32, pW, sizeof(xtfloatx4));
                AE_LSX2X2_IP(W43, W10, pW, sizeof(xtfloatx4));
                AE_LSX2X2_XC1(W21, W0z, pW, sizeof(xtfloatx4));

                AE_LSX2IP(Y01, castxcc(xtfloatx2,pY), sizeof(xtfloatx2));
                AE_LSX2X2_IP(Y23, Y45, pY, sizeof(xtfloatx4));
                AE_LSX2X2_XC(Y67, Y89, pY, (Q - 6) * sizeof(float32_t));

                MADD_SX2X2(S0_0, S1_0, W65, Wz6, Y01, Y01);
                MADD_SX2X2(S2_0, S3_0, W65, Wz6, Y23, Y23);
                MADD_SX2X2(S0_1, S1_1, W43, W54, Y23, Y23);
                MADD_SX2X2(S2_1, S3_1, W43, W54, Y45, Y45);
                MADD_SX2X2(S0_2, S1_2, W21, W32, Y45, Y45);
                MADD_SX2X2(S2_2, S3_2, W21, W32, Y67, Y67);
                MADD_SX2X2(S0_3, S1_3, W0z, W10, Y67, Y67);
                MADD_SX2X2(S2_3, S3_3, W0z, W10, Y89, Y89);
			}

            ADD_SX2X2(S0_0, S2_0, S0_0, S2_0, S0_1, S2_1);
            ADD_SX2X2(S1_0, S3_0, S1_0, S3_0, S1_1, S3_1);

            ADD_SX2X2(S0_2, S2_2, S0_2, S2_2, S0_3, S2_3);
            ADD_SX2X2(S1_2, S3_2, S1_2, S3_2, S1_3, S3_3);

            ADD_SX2X2(S0_0, S2_0, S0_0, S2_0, S0_2, S2_2);
            ADD_SX2X2(S1_0, S3_0, S1_0, S3_0, S1_2, S3_2);

            tmp = CONST_SX2(1);

            MADDMUX_SX2X2 (S0_0, S1_0, tmp, tmp, S0_0, S1_0, 5);
            MADDMUX_SX2X2 (S2_0, S3_0, tmp, tmp, S2_0, S3_0, 5);

			S0_0 = XT_SEL32_HL_SX2(S0_0, S1_0);
			S2_0 = XT_SEL32_HL_SX2(S2_0, S3_0);

			AE_SASX2X2_IP(S0_0, S2_0, alZ, pZ);
		}

		AE_SA128POS_FP(alZ, pZ);
		/* Last N-1 samples of the i-th row */
		{
			xtfloatx2 S0, S2, S4, S0_1;
			xtfloatx2 S1, S3, S5, S1_1;
			xtfloatx2 S0_0, S1_0, S2_0, S3_0;

			S0 = S2 = S4 = S0_1 = (xtfloatx2)(0.0f);
			S1 = S3 = S5 = S1_1 = (xtfloatx2)(0.0f);
            S0_0 = S1_0 = S2_0 = S3_0 = (xtfloatx2)(0.0f);

            pY = (const xtfloatx4*)(y + (y_start+1) * (Q) - 6);

			__Pragma("loop_count min=1, max=11");
			for (m = m_start; m < m_end; m++)
			{
                AE_LSX2IP(Y01, castxcc(xtfloatx2,pY), sizeof(xtfloatx2));
                AE_LSX2X2_XP(Y23, Y45, pY, (Q - 2) * sizeof(float32_t));

                AE_LSX2X2_IP(Wz6, W54, pW, sizeof(xtfloatx4));
                AE_LSX2X2_IP(W65, W32, pW, sizeof(xtfloatx4));
                AE_LSX2X2_IP(W43, W10, pW, sizeof(xtfloatx4));
                AE_LSX2X2_XC1(W21, W0z, pW, sizeof(xtfloatx4));


                MADD_SX2X2(S0, S1, W65, Wz6, Y01, Y01);
                MADD_SX2X2(S2, S3, W65, Wz6, Y23, Y23);
                MADD_SX2X2(S4, S5, W65, Wz6, Y45, Y45);

                MADD_SX2X2(S0_1, S1_1, W43, W54, Y23, Y23);

                MADD_SX2X2(S2_0, S3_0, W43, W54, Y45, Y45);
                MADD_SX2X2(S0_0, S1_0, W21, W32, Y45, Y45);


			}
            ADD_SX2X2(S0, S1, S0, S1, S0_1, S1_1);
            ADD_SX2X2(S0, S1, S0, S1, S0_0, S1_0);
            ADD_SX2X2(S2, S3, S2, S3, S2_0, S3_0);

            tmp = CONST_SX2(1);

            MADDMUX_SX2X2 (S0, S3, tmp, tmp, S0,S3, 5);
            MADDMUX_SX2X2 (S1, S4, tmp, tmp, S1,S4, 5);
            MADDMUX_SX2X2 (S2, S5, tmp, tmp, S2,S5, 5);

			S0 = XT_SEL32_HL_SX2(S0, S1);
			S2 = XT_SEL32_HL_SX2(S2, S3);
			S4 = XT_SEL32_HL_SX2(S4, S5);

			AE_SSX2IP(S0, castxcc(xtfloatx2, pZ), sizeof(xtfloatx2));
			AE_SSX2IP(S2, castxcc(xtfloatx2, pZ), sizeof(xtfloatx2));
			AE_SSX2IP(S4, castxcc(xtfloatx2, pZ), sizeof(xtfloatx2));
		}
	}
#   undef M
#   undef N
}

size_t conv2d_11x7f_getScratchSize(int P, int Q)
{
	(void)P, (void)Q;
	return 11 * WLEN * sizeof(xtfloatx4);
} // MxN=11x7


#else 

#define MAX(a,b)   ( (a)>(b) ? (a) : (b) )
#define MIN(a,b)   ( (a)<(b) ? (a) : (b) )

void conv2d_11x7f     (void* pScr, float32_t *z, const float32_t * x, const float32_t * y, int P, int Q)
{
    NASSERT_ALIGN(x,XCHAL_DATA_WIDTH);
    NASSERT_ALIGN(y,XCHAL_DATA_WIDTH);
    NASSERT_ALIGN(z,XCHAL_DATA_WIDTH);
    NASSERT_ALIGN(pScr,XCHAL_DATA_WIDTH);
    NASSERT(P>=0 && P%8==0);
    NASSERT(Q>=0 && Q%8==0);
    int i, j, m, n, n0, n1, m0, m1;
    float32_t * pwr = (float32_t*)pScr;
    xtfloat * restrict pWwr = (xtfloat*)pwr;
    const xtfloat * restrict pY;
    xtfloat * restrict pX;
    xtfloat * restrict pW;
    xtfloat * restrict pZ;
    xtfloat w0,w1, y0, y1, y2, y3, y4;
    xtfloat s0, s1, s2, s3;
    if (P <= 0 || Q <= 0) return;

#define N 7
#define M 11
    if (P <= 11)
    {
        for (i = 0; i < M + P - 1; i++)
        {
            for (j = 0; j < N + Q - 1; j++)
            {
                float32_t s;
                m0 = XT_MAX(i - P + 1, 0);
                m1 = XT_MIN(i + 1, M);
                n0 = XT_MAX(j - Q + 1, 0);
                n1 = XT_MIN(j + 1, N);
                s = 0;
                for (n = n0; n < n1; n++)
                for (m = m0; m < m1; m++)
                {
                    s += x[m*N + n] * y[(i - m)*Q + (j - n)];
                }
                z[i*(N + Q - 1) + j] = s;
            }
        }
    }
    else
    {
        
        pX = (xtfloat*)x;
        pW = (xtfloat*)pScr + M*N - 1;
        for (i = 0; i < M*N; i++)
        {
            XT_LSIP(w0, pX, sizeof(xtfloat));
            XT_SSXP(w0, pW, -1 * ((int)sizeof(xtfloat)));
        }

        /* first M-1 rows */
        for (i = 0; i < M - 1; i++)
        {
            for (j = 0; j < N - 1; j++)
            {
                float32_t s;
                m0 = XT_MAX(i - P + 1, 0);
                m1 = XT_MIN(i + 1, M);
                n0 = XT_MAX(j - Q + 1, 0);
                n1 = XT_MIN(j + 1, N);
                s = 0;
                for (n = n0; n < n1; n++)
                for (m = m0; m < m1; m++)
                {
                    s += x[m*N + n] * y[(i - m)*Q + (j - n)];
                }
                z[i*(N + Q - 1) + j] = s;
            }
            pZ = (xtfloat *)(z + i*(N + Q - 1) + (N - 1));
            for (j = 0; j < (Q - (N - 1)) >> 2; j++)
            {
                xtfloat w2, y5;
                s0 = XT_MOV_S(0.0f);
                s1 = XT_MOV_S(0.0f);
                s2 = XT_MOV_S(0.0f);
                s3 = XT_MOV_S(0.0f);
                pWwr = (xtfloat*)(pwr + ((M - 1) - i)*N);
                n0 = 4 * j;
                for (m = 0; m < i+1; m++)
                {
                    pY = (const xtfloat *)(y + (m)*Q + n0);
                    for (n = 0; n < 2; n++)
                    {
                        XT_LSIP(w0, pWwr, sizeof(xtfloat));
                        XT_LSIP(w1, pWwr, sizeof(xtfloat));
                        XT_LSIP(w2, pWwr, sizeof(xtfloat));
                        XT_LSIP(y0, pY, sizeof(xtfloat));
                        XT_LSIP(y1, pY, sizeof(xtfloat));
                        XT_LSIP(y2, pY, sizeof(xtfloat));
                        y3 = XT_LSI(pY, 0 * sizeof(xtfloat));
                        y4 = XT_LSI(pY, 1 * sizeof(xtfloat));
                        y5 = XT_LSI(pY, 2 * sizeof(xtfloat));
                        XT_MADD_S(s0, w0, y0);
                        XT_MADD_S(s1, w0, y1);
                        XT_MADD_S(s2, w0, y2);
                        XT_MADD_S(s3, w0, y3);

                        XT_MADD_S(s0, w1, y1);
                        XT_MADD_S(s1, w1, y2);
                        XT_MADD_S(s2, w1, y3);
                        XT_MADD_S(s3, w1, y4);

                        XT_MADD_S(s0, w2, y2);
                        XT_MADD_S(s1, w2, y3);
                        XT_MADD_S(s2, w2, y4);
                        XT_MADD_S(s3, w2, y5);
                    }
                    {
                        XT_LSIP(w0, pWwr, sizeof(xtfloat));
                        XT_LSIP(y0, pY, sizeof(xtfloat));
                        y1 = XT_LSI(pY, 0 * sizeof(xtfloat));
                        y2 = XT_LSI(pY, 1 * sizeof(xtfloat));
                        y3 = XT_LSI(pY, 2 * sizeof(xtfloat));
                        XT_MADD_S(s0, w0, y0);
                        XT_MADD_S(s1, w0, y1);
                        XT_MADD_S(s2, w0, y2);
                        XT_MADD_S(s3, w0, y3);
                    }
                }
                XT_SSIP(s0, pZ, sizeof(xtfloat));
                XT_SSIP(s1, pZ, sizeof(xtfloat));
                XT_SSIP(s2, pZ, sizeof(xtfloat));
                XT_SSIP(s3, pZ, sizeof(xtfloat));
            }
            for (j = (N - 1) + ((Q - (N - 1))&(~3)); j < Q; j++)
            {
                xtfloat s = XT_MOV_S(0.0f);
                pWwr = (xtfloat*)(pwr + ((M - 1) - i)*N);
                n0 = (j - (N - 1));
                for (m = 0; m < i + 1; m++)
                {
                    pY = (const xtfloat *)(y + (m)*Q + n0);
                    for (n = 0; n < N; n++)
                    {
                        XT_LSIP(w0, pWwr, sizeof(xtfloat));
                        XT_LSIP(y0, pY, sizeof(xtfloat));
                        XT_MADD_S(s, w0, y0);
                    }
                }
                XT_SSIP(s, pZ, sizeof(xtfloat));
            }

            for (j = Q; j < N + Q - 1; j++)
            {
                float32_t s;
                m0 = XT_MAX(i - P + 1, 0);
                m1 = XT_MIN(i + 1, M);
                n0 = XT_MAX(j - Q + 1, 0);
                n1 = XT_MIN(j + 1, N);
                s = 0;
                for (n = n0; n < n1; n++)
                for (m = m0; m < m1; m++)
                {
                    s += x[m*N + n] * y[(i - m)*Q + (j - n)];
                }
                z[i*(N + Q - 1) + j] = s;
            }

        }

        for (i = M - 1; i < P; i++)
        {
            pZ = (xtfloat *)(z + i*(N + Q - 1));
            for (j = 0; j < N - 1; j++)
            {
                m0 = i - (M - 1);
                n0 = (N - 1) - j;
                xtfloat s = XT_MOV_S(0.0f);
                for (m = 0; m < M; m++)
                {
                    pY = (const xtfloat *)(y + (m0 + m)*Q);
                    pWwr = (xtfloat *)(pwr + m*N + n0);
                    for (n = 0; n < j + 1; n++)
                    {
                        XT_LSIP(w0, pWwr, sizeof(xtfloat));
                        XT_LSIP(y0, pY, sizeof(xtfloat));
                        XT_MADD_S(s, w0, y0);
                    }
                }
                XT_SSIP(s, pZ, sizeof(xtfloat));
            }

            /* Main section */
            for (j = 0; j < (Q - (N - 1)) >> 2; j++)
            {
                xtfloat w2, y5;
                s0 = XT_MOV_S(0.0f);
                s1 = XT_MOV_S(0.0f);
                s2 = XT_MOV_S(0.0f);
                s3 = XT_MOV_S(0.0f);
                pWwr = (xtfloat*)pwr;
                m0 = i - (M - 1);
                n0 = 4*j;
                for (m = 0; m < M; m++)
                {
                    pY = (const xtfloat *)(y + (m0 + m)*Q + n0);
                    for (n = 0; n < 2; n++)
                    {
                        XT_LSIP(w0, pWwr, sizeof(xtfloat));
                        XT_LSIP(w1, pWwr, sizeof(xtfloat));
                        XT_LSIP(w2, pWwr, sizeof(xtfloat));
                        XT_LSIP(y0, pY, sizeof(xtfloat));
                        XT_LSIP(y1, pY, sizeof(xtfloat));
                        XT_LSIP(y2, pY, sizeof(xtfloat));
                        y3 = XT_LSI(pY, 0*sizeof(xtfloat));
                        y4 = XT_LSI(pY, 1 * sizeof(xtfloat));
                        y5 = XT_LSI(pY, 2 * sizeof(xtfloat));
                        XT_MADD_S(s0, w0, y0);
                        XT_MADD_S(s1, w0, y1);
                        XT_MADD_S(s2, w0, y2);
                        XT_MADD_S(s3, w0, y3);

                        XT_MADD_S(s0, w1, y1);
                        XT_MADD_S(s1, w1, y2);
                        XT_MADD_S(s2, w1, y3);
                        XT_MADD_S(s3, w1, y4);

                        XT_MADD_S(s0, w2, y2);
                        XT_MADD_S(s1, w2, y3);
                        XT_MADD_S(s2, w2, y4);
                        XT_MADD_S(s3, w2, y5);
                    }
                    {
                        XT_LSIP(w0, pWwr, sizeof(xtfloat));
                        XT_LSIP(y0, pY, sizeof(xtfloat));
                        y1 = XT_LSI(pY, 0 * sizeof(xtfloat));
                        y2 = XT_LSI(pY, 1 * sizeof(xtfloat));
                        y3 = XT_LSI(pY, 2 * sizeof(xtfloat));
                        XT_MADD_S(s0, w0, y0);
                        XT_MADD_S(s1, w0, y1);
                        XT_MADD_S(s2, w0, y2);
                        XT_MADD_S(s3, w0, y3);
                    }
                }
                XT_SSIP(s0, pZ, sizeof(xtfloat));
                XT_SSIP(s1, pZ, sizeof(xtfloat));
                XT_SSIP(s2, pZ, sizeof(xtfloat));
                XT_SSIP(s3, pZ, sizeof(xtfloat));
            }
            for (j = (N-1) + ((Q- (N-1))&(~3)); j < Q; j++)
            {
                xtfloat s = XT_MOV_S(0.0f);
                pWwr = (xtfloat*)pwr;
                m0 = i - (M - 1);
                n0 = (j - (N - 1));
                for (m = 0; m < M; m++)
                {
                    pY = (const xtfloat *)(y + (m0 + m)*Q + n0);
                    for (n = 0; n < N; n++)
                    {
                        XT_LSIP(w0, pWwr, sizeof(xtfloat));
                        XT_LSIP(y0, pY, sizeof(xtfloat));
                        XT_MADD_S(s, w0, y0);
                    }
                }
                XT_SSIP(s, pZ, sizeof(xtfloat));
            }

            for (j = Q; j < N + Q - 1; j++)
            {
                xtfloat s = XT_MOV_S(0.0f);
                n1 = N + Q - j - 1;
                m0 = i - (M - 1);
                n0 = j - (N - 1);
                for (m = 0; m < M; m++)
                {
                    pWwr = (xtfloat *)(pwr + m*N);
                    pY = (const xtfloat *)(y + (m0 + m)*Q + n0);
                    for (n = 0; n < n1; n++)
                    {
                        XT_LSIP(w0, pWwr, sizeof(xtfloat));
                        XT_LSIP(y0, pY, sizeof(xtfloat));
                        XT_MADD_S(s, w0, y0);
                    }
                }
                XT_SSIP(s, pZ, sizeof(xtfloat));
            }
        }

        /* Last M-1 rows */
        for (i = P; i < M + P - 1; i++)
        {
            for (j = 0; j < N - 1; j++)
            {
                float32_t s;
                m0 = XT_MAX(i - P + 1, 0);
                m1 = XT_MIN(i + 1, M);
                n0 = XT_MAX(j - Q + 1, 0);
                n1 = XT_MIN(j + 1, N);
                s = 0;
                for (n = n0; n < n1; n++)
                for (m = m0; m < m1; m++)
                {
                    s += x[m*N + n] * y[(i - m)*Q + (j - n)];
                }
                z[i*(N + Q - 1) + j] = s;
            }

            pZ = (xtfloat *)(z + i*(N + Q - 1) + (N - 1));
            for (j = 0; j < (Q - (N - 1)) >> 2; j++)
            {
                xtfloat w2, y5;
                s0 = XT_MOV_S(0.0f);
                s1 = XT_MOV_S(0.0f);
                s2 = XT_MOV_S(0.0f);
                s3 = XT_MOV_S(0.0f);
                pWwr = (xtfloat*)pwr;
                m1 = (M + P - i - 1);
                m0 = i - M + 1;
                n0 = 4*j;
                for (m = 0; m < m1; m++)
                {
                    pY = (const xtfloat *)(y + (m0 + m)*Q + n0);
                    for (n = 0; n < 2; n++)
                    {
                        XT_LSIP(w0, pWwr, sizeof(xtfloat));
                        XT_LSIP(w1, pWwr, sizeof(xtfloat));
                        XT_LSIP(w2, pWwr, sizeof(xtfloat));
                        XT_LSIP(y0, pY, sizeof(xtfloat));
                        XT_LSIP(y1, pY, sizeof(xtfloat));
                        XT_LSIP(y2, pY, sizeof(xtfloat));
                        y3 = XT_LSI(pY, 0*sizeof(xtfloat));
                        y4 = XT_LSI(pY, 1 * sizeof(xtfloat));
                        y5 = XT_LSI(pY, 2 * sizeof(xtfloat));
                        XT_MADD_S(s0, w0, y0);
                        XT_MADD_S(s1, w0, y1);
                        XT_MADD_S(s2, w0, y2);
                        XT_MADD_S(s3, w0, y3);

                        XT_MADD_S(s0, w1, y1);
                        XT_MADD_S(s1, w1, y2);
                        XT_MADD_S(s2, w1, y3);
                        XT_MADD_S(s3, w1, y4);

                        XT_MADD_S(s0, w2, y2);
                        XT_MADD_S(s1, w2, y3);
                        XT_MADD_S(s2, w2, y4);
                        XT_MADD_S(s3, w2, y5);
                    }
                    {
                        XT_LSIP(w0, pWwr, sizeof(xtfloat));
                        XT_LSIP(y0, pY, sizeof(xtfloat));
                        y1 = XT_LSI(pY, 0 * sizeof(xtfloat));
                        y2 = XT_LSI(pY, 1 * sizeof(xtfloat));
                        y3 = XT_LSI(pY, 2 * sizeof(xtfloat));
                        XT_MADD_S(s0, w0, y0);
                        XT_MADD_S(s1, w0, y1);
                        XT_MADD_S(s2, w0, y2);
                        XT_MADD_S(s3, w0, y3);
                    }
                }
                XT_SSIP(s0, pZ, sizeof(xtfloat));
                XT_SSIP(s1, pZ, sizeof(xtfloat));
                XT_SSIP(s2, pZ, sizeof(xtfloat));
                XT_SSIP(s3, pZ, sizeof(xtfloat));
            }
            for (j = (N - 1) + ((Q - (N - 1))&(~3)); j < Q; j++)
            {
                xtfloat s = XT_MOV_S(0.0f);
                pWwr = (xtfloat*)pwr;
                n0 = (j - (N - 1));
                m1 = (M + P - i - 1);
                m0 = i - M + 1;
                for (m = 0; m < m1; m++)
                {
                    pY = (const xtfloat *)(y + (m0 + m)*Q + n0);
                    for (n = 0; n < N; n++)
                    {
                        XT_LSIP(w0, pWwr, sizeof(xtfloat));
                        XT_LSIP(y0, pY, sizeof(xtfloat));
                        XT_MADD_S(s, w0, y0);
                    }
                }
                XT_SSIP(s, pZ, sizeof(xtfloat));
            }
            for (j = Q; j < N + Q - 1; j++)
            {
                float32_t s;
                m0 = XT_MAX(i - P + 1, 0);
                m1 = XT_MIN(i + 1, M);
                n0 = XT_MAX(j - Q + 1, 0);
                n1 = XT_MIN(j + 1, N);
                s = 0;
                for (n = n0; n < n1; n++)
                for (m = m0; m < m1; m++)
                {
                    s += x[m*N + n] * y[(i - m)*Q + (j - n)];
                }
                z[i*(N + Q - 1) + j] = s;
            }
        }
    }
#undef M
#undef N
}

size_t conv2d_11x7f_getScratchSize(int P, int Q)
{ 
    return 11 * 7 * sizeof(xtfloat);
} // MxN=11x7

#endif
