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
    2D-IDCT, 16-bit input & 8-bit output with no scaling
    C code optimized for HiFi4
    Integrit, 2006-2019
*/

/* Portable data types. */
#include "NatureDSP_types.h"
/* Signal Processing Library API. */
#include "NatureDSP_Signal_fft.h"
/* Common utility and macros declarations. */
#include "common.h"
/* Twiddle factor tables for DCTs. */
#include "dct2_twd.h"

/*-------------------------------------------------------------------------
  2-D Inverse Discrete Cosine Transform.
  These functions apply inverse DCT (Type II) to the series of L input 
  blocks of NxN pixels. Algorithm uses ITU-T T.81 (JPEG compression) DCT-II 
  definition with bias 128 and left-to-right, top-to-bottom orientation.

  Scaling:
      +-----------------------+--------------------------------------+
      |      Function         |           Scaling options            |
      +-----------------------+--------------------------------------+
      |       idct2d_16x8     |           0 - no scaling             |
      +-----------------------+--------------------------------------+
  Notes:
  N - IDCT size (depends on selected IDCT handle)

  Precision: 
  16x8  16-bit signed input, 8-bit unsigned output

  Input:
  x[N*N*L]    input data: L NxN blocks
  h           DCT handle
  L           number of input blocks
  scalingOpt  scaling option (see table above), should be 0

  Output:
  y[N*N*L]    output pixels: L NxN blocks

  Returned value: 0
  Restriction:
  x,y         should not overlap
  x,y         aligned on 16-bytes boundary

-------------------------------------------------------------------------*/

/*
 * 2D-DCT transform matrix, Q15
 * NOTE: coeffs are multiplied by sqrt(2)
 */
static const int16_t ALIGN(32) Cuv[] =
{
    16384,   22725,   21407,   19266,   16384,   12873,    8867,    4520,
    16384,   19266,    8867,   -4520,  -16384,  -22725,  -21407,  -12873,
    16384,   12873,   -8867,  -22725,  -16384,    4520,   21407,   19266,
    16384,    4520,  -21407,  -12873,   16384,   19266,   -8867,  -22725,
    16384,   -4520,  -21407,   12873,   16384,  -19266,   -8867,   22725,
    16384,  -12873,   -8867,   22725,  -16384,   -4520,   21407,  -19266,
    16384,  -19266,    8867,    4520,  -16384,   22725,  -21407,   12873,
    16384,  -22725,   21407,  -19266,   16384,  -12873,    8867,   -4520,
};
static const tdct2_twd twd={0,8,NULL,(void *)Cuv};
const dct_handle_t idct2d_16_8=(const dct_handle_t*)&twd;

int idct2d_16x8(uint8_t * y, int16_t * x, dct_handle_t h, int L, int scalingOpt)
{
    int16_t ALIGN(16) rows[8*8];
    const tdct2_twd *ptwd=(const tdct2_twd *)h;
          ae_int8x8  * restrict pY;
    const ae_int16x8 * restrict pX;
          ae_int16x8 * pR_wr;
    const ae_int16x8 * pR_rd;
    const ae_int16x8 * restrict pC;
    ae_int64 ACC0, ACC1, ACC2, ACC3, ACC4, ACC5, ACC6, ACC7;
    ae_int16x4 c0, c1, y0, y1;
    ae_int8x8 out;
    ae_int16x4 x00, x01, x10, x11, x20, x21, x30, x31,
               x40, x41, x50, x51, x60, x61, x70, x71;
    int l, j;

    NASSERT_ALIGN16(x);
    NASSERT_ALIGN16(y);
    NASSERT(scalingOpt == 0);
    NASSERT((ptwd->N) == 8);

    pX = (const ae_int16x8 *)(x);
    pY = (ae_int8x8 *)(y);

    for (l=0; l<L; l++)
    {
        /* Process rows */
        pR_wr = (ae_int16x8 *)(rows);
        pC=(const ae_int16x8 *)(ptwd->fft_twd);

        AE_L16X4X2_IP(x00, x01, pX, sizeof(ae_int16x8));
        AE_L16X4X2_IP(x10, x11, pX, sizeof(ae_int16x8));
        AE_L16X4X2_IP(x20, x21, pX, sizeof(ae_int16x8));
        AE_L16X4X2_IP(x30, x31, pX, sizeof(ae_int16x8));
        AE_L16X4X2_IP(x40, x41, pX, sizeof(ae_int16x8));
        AE_L16X4X2_IP(x50, x51, pX, sizeof(ae_int16x8));
        AE_L16X4X2_IP(x60, x61, pX, sizeof(ae_int16x8));
        AE_L16X4X2_IP(x70, x71, pX, sizeof(ae_int16x8));
        __Pragma("no_unroll");
        for (j=0; j<8; j++)
        {
            AE_L16X4X2_IP(c0, c1, pC, sizeof(ae_int16x8));

            /* Q15 */
            AE_MULZAAAA2Q16(ACC0, ACC1, c0, c0, x00, x10);  AE_MULAAAA2Q16(ACC0, ACC1, c1, c1, x01, x11);
            AE_MULZAAAA2Q16(ACC2, ACC3, c0, c0, x20, x30);  AE_MULAAAA2Q16(ACC2, ACC3, c1, c1, x21, x31);
            AE_MULZAAAA2Q16(ACC4, ACC5, c0, c0, x40, x50);  AE_MULAAAA2Q16(ACC4, ACC5, c1, c1, x41, x51);
            AE_MULZAAAA2Q16(ACC6, ACC7, c0, c0, x60, x70);  AE_MULAAAA2Q16(ACC6, ACC7, c1, c1, x61, x71);

            /* Q0 <- Q15 - 15 w/ rounding and saturation */
            AE_PKSR16(y0, ACC0, 1);  AE_PKSR16(y0, ACC1, 1);
            AE_PKSR16(y0, ACC2, 1);  AE_PKSR16(y0, ACC3, 1);
            AE_PKSR16(y1, ACC4, 1);  AE_PKSR16(y1, ACC5, 1);
            AE_PKSR16(y1, ACC6, 1);  AE_PKSR16(y1, ACC7, 1);

            AE_S16X4X2_IP(y0, y1, pR_wr, sizeof(ae_int16x8));
        }

        /* Process columns */
        pR_rd = (const ae_int16x8 *)(rows);
        pC=(const ae_int16x8 *)(ptwd->fft_twd);

        AE_L16X4X2_IP(x00, x01, pR_rd, sizeof(ae_int16x8));
        AE_L16X4X2_IP(x10, x11, pR_rd, sizeof(ae_int16x8));
        AE_L16X4X2_IP(x20, x21, pR_rd, sizeof(ae_int16x8));
        AE_L16X4X2_IP(x30, x31, pR_rd, sizeof(ae_int16x8));
        AE_L16X4X2_IP(x40, x41, pR_rd, sizeof(ae_int16x8));
        AE_L16X4X2_IP(x50, x51, pR_rd, sizeof(ae_int16x8));
        AE_L16X4X2_IP(x60, x61, pR_rd, sizeof(ae_int16x8));
        AE_L16X4X2_IP(x70, x71, pR_rd, sizeof(ae_int16x8));

        for (j=0; j<8; j++)
        {
            AE_L16X4X2_IP(c0, c1, pC, sizeof(ae_int16x8));

            /* Q15 <- Q15*Q0 */
            AE_MULZAAAA2Q16(ACC0, ACC1, c0, c0, x00, x10);  AE_MULAAAA2Q16(ACC0, ACC1, c1, c1, x01, x11);
            AE_MULZAAAA2Q16(ACC2, ACC3, c0, c0, x20, x30);  AE_MULAAAA2Q16(ACC2, ACC3, c1, c1, x21, x31);
            AE_MULZAAAA2Q16(ACC4, ACC5, c0, c0, x40, x50);  AE_MULAAAA2Q16(ACC4, ACC5, c1, c1, x41, x51);
            AE_MULZAAAA2Q16(ACC6, ACC7, c0, c0, x60, x70);  AE_MULAAAA2Q16(ACC6, ACC7, c1, c1, x61, x71);

            /* Q0 <- Q15-1 - 15 w/ rounding and saturation */
            AE_PKSR16(y0, ACC0, 0);  AE_PKSR16(y0, ACC1, 0);
            AE_PKSR16(y0, ACC2, 0);  AE_PKSR16(y0, ACC3, 0);
            AE_PKSR16(y1, ACC4, 0);  AE_PKSR16(y1, ACC5, 0);
            AE_PKSR16(y1, ACC6, 0);  AE_PKSR16(y1, ACC7, 0);
            out = AE_SAT8X8X16(y0, y1);
            /* Add bias */
            out = AE_ADD8(out, AE_MOVDA8(128));
            AE_S8X8_IP(out, pY, 8*sizeof(int8_t));
        }
    }

    return 0;
} /* idct2d_16x8() */
