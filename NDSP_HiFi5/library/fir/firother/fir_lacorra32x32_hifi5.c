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
    Real data linear auto-correlation, 32x32-bit, no requirements on vectors
    length and alignment.
    C code optimized for HiFi5
  IntegrIT, 2006-2019
*/

/* Portable data types. */
#include "NatureDSP_types.h"
/* Signal Processing Library API. */
#include "NatureDSP_Signal_fir.h"
/* Common utility and macros declarations. */
#include "common.h"

/*-------------------------------------------------------------------------
  Linear Autocorrelation 
  Functions estimate the linear auto-correlation of vector x. Returns 
  autocorrelation of length N.

  Precision: 
  16x16   16-bit data, 16-bit outputs
  32x32   32-bit data, 32-bit outputs

  Input:
  s[]       scratch area of
            FIR_LACORRA16X16_SCRATCH_SIZE( N )
            FIR_LACORRA32X32_SCRATCH_SIZE( N ) bytes
            
  x[N]      input data Q15, Q31 
  N         length of x
  Output:
  r[N]      output data, Q15, Q31

  Restrictions:
  x,r,s   should not overlap
  N       >0
  s       aligned on an 16-bytes boundary
-------------------------------------------------------------------------*/

void fir_lacorra32x32(void    * restrict s,
                      int32_t * restrict r,
                const int32_t * restrict x, int N)
{
    void       * s_ptr;
    int32_t    * x_buf;
    const ae_int32x4 * restrict S;
          ae_int32x4 * restrict D;
    ae_valignx2 S_va;

    const ae_int32x2 *          pX;
    const ae_int32x2 *          S0;
    const ae_int32x2 *          S1;
    const ae_int32x2 *          S2;
    const ae_int32x2 *          pY;
          ae_int32x4 * restrict pR;
    ae_valignx2 ar2;
    int n, m, m_iter;

    NASSERT(s);
    NASSERT(r);
    NASSERT(x);
    NASSERT_ALIGN(s, 16);
    NASSERT(N > 0);

    //----------------------------------------------------------------------------
    // Partition the scratch memory area.
    s_ptr = s;
    x_buf = (int32_t*)s_ptr;
    s_ptr = x_buf + N + 5;
    ASSERT((int8_t *)s_ptr - (int8_t *)s <= (int)FIR_LACORRA32X32_SCRATCH_SIZE(N));

    //----------------------------------------------------------------------------
    // Copy x[N].
    S = (const ae_int32x4 *)x;
    D = (      ae_int32x4 *)x_buf;
    S_va = AE_LA128_PP(S);
    for (n = 0; n < (N >> 2); n++)
    {
        ae_int32x2 x0, x1;
        AE_LA32X2X2_IP(x0, x1, S_va, S);
        AE_S32X2X2_IP(x0, x1, D, 4 * sizeof(int32_t));
    }
    for (n = (N&~3); n < N; n++)
    {
        ae_int32x2 x0;
        AE_L32_IP(x0, castxcc(ae_int32, S), sizeof(int32_t));
        AE_S32_L_IP(x0, castxcc(ae_int32, D), sizeof(int32_t));
    }
    AE_S32_L_IP(AE_ZERO32(), castxcc(ae_int32, D), sizeof(int32_t));
    AE_S32_L_IP(AE_ZERO32(), castxcc(ae_int32, D), sizeof(int32_t));
    AE_S32_L_IP(AE_ZERO32(), castxcc(ae_int32, D), sizeof(int32_t));
    AE_S32_L_IP(AE_ZERO32(), castxcc(ae_int32, D), sizeof(int32_t));
    AE_S32_L_IP(AE_ZERO32(), castxcc(ae_int32, D), sizeof(int32_t));

    m_iter = ((N - 3) >> 2);

    pX = (const ae_int32x2 *)x_buf;
    pR = (      ae_int32x4 *)r;
    ar2 = AE_ZALIGN128();

    for (n = 0; n < (N >> 3); n++, m_iter -= 2)
    {
        ae_f64 q0, q1, q2, q3, q4, q5, q6, q7;
        ae_int32x2 x0, x1, x2, x3, x4, x5;
        ae_int32x2 y0, y1;

        pY = (const ae_int32x2 *)x_buf;
        S0 = pX;
        S1 = pX + 2;
        pX += 4;
        S2 = pX;

        q0 = q1 = q2 = q3 = q4 = q5 = q6 = q7 = AE_ZERO64();

        __Pragma("loop_count min=1");
        for (m = 0; m < m_iter; m++)
        {
            AE_L32X2_IP(y0, pY, 2 * sizeof(int32_t));
            AE_L32X2_IP(y1, pY, 2 * sizeof(int32_t));

            AE_L32X2X2_IP(x0, x1, castxcc(ae_int32x4, S0), 4 * sizeof(int32_t));
            AE_L32X2X2_IP(x2, x3, castxcc(ae_int32x4, S1), 4 * sizeof(int32_t));
            AE_L32X2X2_IP(x4, x5, castxcc(ae_int32x4, S2), 4 * sizeof(int32_t));

            AE_MULAFD32X2RA_FIR_H(q0, q1, x0, x1, y0);
            AE_MULAFD32X2RA_FIR_H(q2, q3, x1, x2, y0);
            AE_MULAFD32X2RA_FIR_H(q4, q5, x2, x3, y0);
            AE_MULAFD32X2RA_FIR_H(q6, q7, x3, x4, y0);

            AE_MULAFD32X2RA_FIR_H(q0, q1, x1, x2, y1);
            AE_MULAFD32X2RA_FIR_H(q2, q3, x2, x3, y1);
            AE_MULAFD32X2RA_FIR_H(q4, q5, x3, x4, y1);
            AE_MULAFD32X2RA_FIR_H(q6, q7, x4, x5, y1);
        }
        AE_L32X2_IP(y0, pY, 2 * sizeof(int32_t));
        AE_L32X2_IP(y1, pY, 2 * sizeof(int32_t));

        AE_L32X2X2_IP(x0, x1, castxcc(ae_int32x4, S0), 4 * sizeof(int32_t));
        AE_L32X2X2_IP(x2, x3, castxcc(ae_int32x4, S1), 4 * sizeof(int32_t));

        AE_MULAFD32X2RA_FIR_H(q0, q1, x0, x1, y0);
        AE_MULAFD32X2RA_FIR_H(q2, q3, x1, x2, y0);
        AE_MULAFD32X2RA_FIR_H(q4, q5, x2, x3, y0);

        AE_MULAFD32X2RA_FIR_H(q0, q1, x1, x2, y1);
        AE_MULAFD32X2RA_FIR_H(q2, q3, x2, x3, y1);

        AE_L32X2_IP(y0, pY, 2 * sizeof(int32_t));
        AE_MULAFD32X2RA_FIR_H(q0, q1, x2, x3, y0);

        x0 = AE_ROUND32X2F48SASYM(q0, q1);
        x1 = AE_ROUND32X2F48SASYM(q2, q3);
        x2 = AE_ROUND32X2F48SASYM(q4, q5);
        x3 = AE_ROUND32X2F48SASYM(q6, q7);
        AE_SA32X2X2_IP(x0, x1, ar2, pR);
        AE_SA32X2X2_IP(x2, x3, ar2, pR);
    }

    //pX = (const ae_int32x2 *)(x_buf + (N&~7));
    //m_iter = ((N & 7) - 1) >> 1;
    if (N & 4)
    {
        ae_f64 q0, q1, q2, q3;
        ae_int32x2 x0, x1, x2, y0;

        pY = (const ae_int32x2 *)x_buf;
        AE_L32X2_IP(y0, pY, 2 * sizeof(int32_t));

        AE_L32X2_IP(x0, pX, 2 * sizeof(int32_t));
        AE_L32X2_IP(x1, pX, 2 * sizeof(int32_t));
        S0 = pX;
        AE_L32X2_IP(x2, S0, 2 * sizeof(int32_t));

        AE_MULFD32X2RA_FIR_H(q0, q1, x0, x1, y0);
        AE_MULFD32X2RA_FIR_H(q2, q3, x1, x2, y0);

        if (N & 3)//m_iter>1
        {
            AE_L32X2_IP(y0, pY, 2 * sizeof(int32_t));
            x0 = x1; x1 = x2;
            AE_L32X2_IP(x2, S0, 2 * sizeof(int32_t));
            AE_MULAFD32X2RA_FIR_H(q0, q1, x0, x1, y0);
            AE_MULAFD32X2RA_FIR_H(q2, q3, x1, x2, y0);
        }
        if ((N & 3) == 3)//m_iter>2
        {
            AE_L32X2_IP(y0, pY, 2 * sizeof(int32_t));
            x0 = x1; x1 = x2;
            AE_L32X2_IP(x2, S0, 2 * sizeof(int32_t));
            AE_MULAFD32X2RA_FIR_H(q0, q1, x0, x1, y0);
            AE_MULAFD32X2RA_FIR_H(q2, q3, x1, x2, y0);
        }

        AE_L32X2_IP(y0, pY, 2 * sizeof(int32_t));
        AE_MULAFD32X2RA_FIR_H(q0, q1, x1, x2, y0);

        x0 = AE_ROUND32X2F48SASYM(q0, q1);
        x1 = AE_ROUND32X2F48SASYM(q2, q3);
        AE_SA32X2X2_IP(x0, x1, ar2, pR);
        //m_iter -= 2;
    }
    AE_SA128POS_FP(ar2, pR);

    if (N & 2)
    {
        ae_f64 q0, q1;
        ae_int32x2 x0, x1, x2, y0;
        pY = (const ae_int32x2 *)x_buf;
        AE_L32X2_IP(y0, pY, 2 * sizeof(int32_t));
        AE_L32X2_IP(x0, pX, 2 * sizeof(int32_t));
        x1 = AE_L32X2_I(pX, 0);
        AE_MULFD32X2RA_FIR_H(q0, q1, x0, x1, y0);
        if (N & 1)//m_iter==1
        {
            AE_L32X2_IP(y0, pY, 2 * sizeof(int32_t));
            x2 = AE_L32X2_I(pX, 2 * sizeof(int32_t));
            AE_MULAFD32X2RA_FIR_H(q0, q1, x1, x2, y0);
        }
        AE_S32RA64S_IP(q0, castxcc(ae_int32, pR), sizeof(int32_t));
        AE_S32RA64S_IP(q1, castxcc(ae_int32, pR), sizeof(int32_t));
    }

    if (N & 1)
    {
        ae_int32x2 x0, y0;
        y0 = AE_L32_I((const ae_int32 *)x_buf, 0);
        x0 = AE_L32_I((const ae_int32 *)pX, 0);
        x0 = AE_MULFP32X2RAS(x0, y0);
        AE_S32_L_I(x0, (ae_int32 *)pR, 0);
    }
} /* fir_lacorra32x32() */
