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
#include "chol32x32_common.h"

/* --------------------------------------------------------
   reversing R matrices for easier readings by rows (diagonal 
   elements are omitted):
   original R    transformed R
   0 1 3 6 a     d 8 c 4 7 b 1 3 6 a
     2 4 7 b
       5 8 c
         9 d
           e

   Input:
   R[N*(N+1)/2]    input matrix
   Output:
   Rt[N*(N-1)/2]   transposed matrix
--------------------------------------------------------*/
void cmatcholtransformR(complex_fract32 *Rt, const complex_fract32* R, int N)
#if 0
{
    int n, m;
    const complex_fract32* pR;
    for (n = 0; n < N; n++)
    {
        pR = R + (((N - n)*(N - n + 3) - 2)>>1);
        for (m = 0; m < n; m++)
        {
            Rt[0] = pR[0];
            pR += (N - n + m + 1);
            Rt += 1;
        }
    }
} 
#else
{
    int n, m,p,d;
    const ae_int32x2* pR;
    R += (N*(N + 3) - 2)>>1;
    d=N+1;
    for (n = 0; n < N; n++)
    {
        pR = (const ae_int32x2*)R;
        R-=d;
        d--;
        p=(N - n + 1)*sizeof(ae_int32x2);

        for (m = 0; m < n; m++)
        {
            ae_int32x2 t;
            AE_L32X2_XP(t,pR,p);
            p+=sizeof(ae_int32x2);
            AE_S32X2_IP(t,castxcc(ae_int32x2,Rt),sizeof(ae_int32x2));
        }
    }
} 
#endif
