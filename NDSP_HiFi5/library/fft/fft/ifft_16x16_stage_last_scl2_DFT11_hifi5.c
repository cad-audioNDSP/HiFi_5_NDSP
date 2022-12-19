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
    NatureDSP Signal Processing Library. FFT part
    Integrit, 2006-2020

    HiFi5 code
*/
#include "NatureDSP_Signal_fft.h"
#include "NatureDSP_Signal_vector.h"
#include "common.h"
/* Twiddle factor tables and FFT descriptor structure. */
#include "fft_x16_common.h"
#include "fft_16x16_stages.h"

int ifft_16x16_stage_last_scl2_DFT11(const int16_t *tw, const int16_t *x, int16_t *y, int N, int *v, int tw_step, int *bexp)
#if 0
{
    int i, res;
    NASSERT_ALIGN16(x);
    NASSERT_ALIGN16(y);
    NASSERT(N%8==0);
    res = fft_16x16_stage_last_scl2_DFT11(tw, x, y, N, v, tw_step, bexp);
    for (i = 0; i < 2 * N; i += 2)
    {
        int16_t re, im;
        re = y[i];
        im = y[i + 1];
        y[i] = im;
        y[i + 1] = re;
    }
    return res;
}
#else
{
    const ae_int16x8 *pX=(const ae_int16x8 *)y;
          ae_int32x4 *pY=(      ae_int32x4 *)y;
    int i, res;
    NASSERT_ALIGN16(x);
    NASSERT_ALIGN16(y);
    NASSERT(N%8==0);
    res = fft_16x16_stage_last_scl2_DFT11(tw, x, y, N, v, tw_step, bexp);
    // copy using implicit endianess conversion
    __Pragma("loop_count factor=2")
    for (i = 0; i < (N>>2); i++)
    {
        ae_int16x4 x0,x1;
        AE_L16X4X2_IP(x0,x1,pX,sizeof(ae_int16x8));
        AE_S32X2X2_IP(AE_MOVINT32X2_FROMINT16X4(x0),
                      AE_MOVINT32X2_FROMINT16X4(x1),
                      pY,sizeof(ae_int16x8));
    }
    return res;
}
#endif
