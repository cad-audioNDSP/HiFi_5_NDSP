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
    2D-DCT, 8-bit input & 16-bit output with no scaling
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
  2-D Discrete Cosine Transform.
  These functions apply DCT (Type II) to the series of L input blocks 
  of NxN pixels. Algorithm uses ITU-T T.81 (JPEG compression) DCT-II 
  definition with bias 128 and left-to-right, top-to-bottom orientation.

  Scaling:
      +-----------------------+--------------------------------------+
      |      Function         |           Scaling options            |
      +-----------------------+--------------------------------------+
      |       dct2d_8x16      |           0 - no scaling             |
      +-----------------------+--------------------------------------+
  Notes:
  N - DCT size (depends on selected DCT handle)

  Precision: 
  8x16  8-bit unsigned input, 16-bit signed output

  Input:
  x[N*N*L]    input pixels: L NxN blocks
  h           DCT handle
  L           number of input blocks
  scalingOpt  scaling option (see table above), should be 0

  Output:
  y[N*N*L]    output of transform: L NxN blocks
  
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
    16384,   16384,   16384,   16384,   16384,   16384,   16384,   16384,
    22725,   19266,   12873,    4520,   -4520,  -12873,  -19266,  -22725,
    21407,    8867,   -8867,  -21407,  -21407,   -8867,    8867,   21407,
    19266,   -4520,  -22725,  -12873,   12873,   22725,    4520,  -19266,
    16384,  -16384,  -16384,   16384,   16384,  -16384,  -16384,   16384,
    12873,  -22725,    4520,   19266,  -19266,   -4520,   22725,  -12873,
     8867,  -21407,   21407,   -8867,   -8867,   21407,  -21407,    8867,
     4520,  -12873,   19266,  -22725,   22725,  -19266,   12873,   -4520,
};

static const tdct2_twd twd={0,8,NULL,(void *)Cuv};
const dct_handle_t dct2d_16_8=(const dct_handle_t*)&twd;

#ifndef AE_MULUS8Q8X16
int dct2d_8x16(int16_t* y, uint8_t * x, dct_handle_t h, int L, int scalingOpt)
{  
    int16_t ALIGN(16) rows[8*8];
    const tdct2_twd *ptwd=(const tdct2_twd *)h;
    const int8_t * pX;
          ae_int16x8 * pY;
          ae_int16x8 * pR_wr;
    const ae_int16x8 * pR_rd;
    const ae_int16x8 * restrict pC;
    ae_int64 ACC0, ACC1, ACC2, ACC3, ACC4, ACC5, ACC6, ACC7;
    ae_int32x2 Y0, Y1, Y2, Y3;
    ae_int16x4 c0, c1, y0, y1;
    ae_int16x4 x00, x01, x10, x11, x20, x21, x30, x31,
               x40, x41, x50, x51, x60, x61, x70, x71;
    int l, j;

    NASSERT_ALIGN16(x);
    NASSERT_ALIGN16(y);
    NASSERT(scalingOpt == 0);
    NASSERT((ptwd->N) == 8);

    pX = (const int8_t *)(x);
    pY = (ae_int16x8 *)(y);

    for (l=0; l<L; l++)
    {
        /* Process rows */
        pR_wr = (ae_int16x8 *)(rows);
        pC=(const ae_int16x8 *)(ptwd->fft_twd);

        AE_L8X4U_IP(x00, pX, 4*sizeof(uint8_t));    AE_L8X4U_IP(x01, pX, 4*sizeof(uint8_t));
        AE_L8X4U_IP(x10, pX, 4*sizeof(uint8_t));    AE_L8X4U_IP(x11, pX, 4*sizeof(uint8_t));
        AE_L8X4U_IP(x20, pX, 4*sizeof(uint8_t));    AE_L8X4U_IP(x21, pX, 4*sizeof(uint8_t));
        AE_L8X4U_IP(x30, pX, 4*sizeof(uint8_t));    AE_L8X4U_IP(x31, pX, 4*sizeof(uint8_t));
        AE_L8X4U_IP(x40, pX, 4*sizeof(uint8_t));    AE_L8X4U_IP(x41, pX, 4*sizeof(uint8_t));
        AE_L8X4U_IP(x50, pX, 4*sizeof(uint8_t));    AE_L8X4U_IP(x51, pX, 4*sizeof(uint8_t));
        AE_L8X4U_IP(x60, pX, 4*sizeof(uint8_t));    AE_L8X4U_IP(x61, pX, 4*sizeof(uint8_t));
        AE_L8X4U_IP(x70, pX, 4*sizeof(uint8_t));    AE_L8X4U_IP(x71, pX, 4*sizeof(uint8_t));
        for (j=0; j<8; j++)
        {
            AE_L16X4X2_IP(c0, c1, pC, sizeof(ae_int16x8));

            /* Q15 */
            AE_MULZAAAA2Q16(ACC0, ACC1, c0, c0, x00, x10);  AE_MULAAAA2Q16(ACC0, ACC1, c1, c1, x01, x11);
            AE_MULZAAAA2Q16(ACC2, ACC3, c0, c0, x20, x30);  AE_MULAAAA2Q16(ACC2, ACC3, c1, c1, x21, x31);
            AE_MULZAAAA2Q16(ACC4, ACC5, c0, c0, x40, x50);  AE_MULAAAA2Q16(ACC4, ACC5, c1, c1, x41, x51);
            AE_MULZAAAA2Q16(ACC6, ACC7, c0, c0, x60, x70);  AE_MULAAAA2Q16(ACC6, ACC7, c1, c1, x61, x71);
            /* Q16 <- Q15 + 1 w/ saturation */
            Y0 = AE_TRUNCA32X2F64S(ACC0, ACC1, 32+1);
            Y1 = AE_TRUNCA32X2F64S(ACC2, ACC3, 32+1);
            Y2 = AE_TRUNCA32X2F64S(ACC4, ACC5, 32+1);
            Y3 = AE_TRUNCA32X2F64S(ACC6, ACC7, 32+1);
            /* Q0 <- Q16 - 16 w/ rounding and saturation */
            y0 = AE_ROUND16X4F32SASYM(Y0, Y1);
            y1 = AE_ROUND16X4F32SASYM(Y2, Y3);
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

        {
            AE_L16X4X2_IP(c0, c1, pC, sizeof(ae_int16x8));

            /* Q15 <- Q15*Q0 */
            AE_MULZAAAA2Q16(ACC0, ACC1, c0, c0, x00, x10);  AE_MULAAAA2Q16(ACC0, ACC1, c1, c1, x01, x11);
            AE_MULZAAAA2Q16(ACC2, ACC3, c0, c0, x20, x30);  AE_MULAAAA2Q16(ACC2, ACC3, c1, c1, x21, x31);
            AE_MULZAAAA2Q16(ACC4, ACC5, c0, c0, x40, x50);  AE_MULAAAA2Q16(ACC4, ACC5, c1, c1, x41, x51);
            AE_MULZAAAA2Q16(ACC6, ACC7, c0, c0, x60, x70);  AE_MULAAAA2Q16(ACC6, ACC7, c1, c1, x61, x71);
            /* Q15 <- Q15 w/ saturation (convert 64 to 32-bit) */
            Y0 = AE_TRUNCA32X2F64S(ACC0, ACC1, 32);
            Y1 = AE_TRUNCA32X2F64S(ACC2, ACC3, 32);
            Y2 = AE_TRUNCA32X2F64S(ACC4, ACC5, 32);
            Y3 = AE_TRUNCA32X2F64S(ACC6, ACC7, 32);
            /* Remove DC bias */
            Y0 = AE_SUB32S(Y0, AE_MOVDA32X2((8*128)<<16, 0));
            /* Q0 <- Q15-1 - 15 w/ rounding and saturation */
            y0 = AE_ROUND16X4F32SASYM(Y0, Y1);
            y1 = AE_ROUND16X4F32SASYM(Y2, Y3);
            AE_S16X4X2_IP(y0, y1, pY, sizeof(ae_int16x8));
        }
        for (j=1; j<8; j++)
        {
            AE_L16X4X2_IP(c0, c1, pC, sizeof(ae_int16x8));

            /* Q15 <- Q15*Q0 */
            AE_MULZAAAA2Q16(ACC0, ACC1, c0, c0, x00, x10);  AE_MULAAAA2Q16(ACC0, ACC1, c1, c1, x01, x11);
            AE_MULZAAAA2Q16(ACC2, ACC3, c0, c0, x20, x30);  AE_MULAAAA2Q16(ACC2, ACC3, c1, c1, x21, x31);
            AE_MULZAAAA2Q16(ACC4, ACC5, c0, c0, x40, x50);  AE_MULAAAA2Q16(ACC4, ACC5, c1, c1, x41, x51);
            AE_MULZAAAA2Q16(ACC6, ACC7, c0, c0, x60, x70);  AE_MULAAAA2Q16(ACC6, ACC7, c1, c1, x61, x71);
            /* Q15 <- Q15 w/ saturation (convert 64 to 32-bit) */
            Y0 = AE_TRUNCA32X2F64S(ACC0, ACC1, 32);
            Y1 = AE_TRUNCA32X2F64S(ACC2, ACC3, 32);
            Y2 = AE_TRUNCA32X2F64S(ACC4, ACC5, 32);
            Y3 = AE_TRUNCA32X2F64S(ACC6, ACC7, 32);
            /* Q0 <- Q15-1 - 15 w/ rounding and saturation */
            y0 = AE_ROUND16X4F32SASYM(Y0, Y1);
            y1 = AE_ROUND16X4F32SASYM(Y2, Y3);
            AE_S16X4X2_IP(y0, y1, pY, sizeof(ae_int16x8));
        }
    }

    return 0;
} /* dct2d_8x16() */
#else
int dct2d_8x16(int16_t* y, uint8_t * x, dct_handle_t h, int L, int scalingOpt)
{  
    int16_t ALIGN(16) rows[8*8];
    const tdct2_twd *ptwd=(const tdct2_twd *)h;
    const ae_int8x16 * pX;
          ae_int16x8 * pY;
          ae_int16x8 * pR_wr;
    const ae_int16x8 * pR_rd;
    const ae_int16x8 * restrict pC;
    ae_int64 ACC0, ACC1, ACC2, ACC3, ACC4, ACC5, ACC6, ACC7;
    ae_int32x2 Y0, Y1, Y2, Y3;
    ae_int16x4 c0, c1, y0, y1;
    ae_int16x4 x00, x01, x10, x11, x20, x21, x30, x31,
               x40, x41, x50, x51, x60, x61, x70, x71;
    ae_int8x8 x0, x1, x2, x3, x4, x5, x6, x7;
    int l, j;

    NASSERT_ALIGN16(x);
    NASSERT_ALIGN16(y);
    NASSERT(scalingOpt == 0);
    NASSERT((ptwd->N) == 8);

    pX = (const ae_int8x16 *)(x);
    pY = (      ae_int16x8 *)(y);

    for (l=0; l<L; l++)
    {
        /* Process rows */
        pR_wr = (ae_int16x8 *)(rows);
        pC=(const ae_int16x8 *)(ptwd->fft_twd);

        AE_L8X8X2_IP(x0, x1, pX, 2*8*sizeof(uint8_t));
        AE_L8X8X2_IP(x2, x3, pX, 2*8*sizeof(uint8_t));
        AE_L8X8X2_IP(x4, x5, pX, 2*8*sizeof(uint8_t));
        AE_L8X8X2_IP(x6, x7, pX, 2*8*sizeof(uint8_t));
        for (j=0; j<8; j++)
        {
            AE_L16X4X2_IP(c0, c1, pC, sizeof(ae_int16x8));

            /* Q15 */
            AE_MULUS8Q8X16(Y0, Y1, x0, x1, x2, x3, c0, c1);
            AE_MULUS8Q8X16(Y2, Y3, x4, x5, x6, x7, c0, c1);
            /* Q16 <- Q15 + 1 w/ saturation */
            AE_MUL2P32X4S(Y0, Y1, Y0, Y1, AE_MOVDA32(2), AE_MOVDA32(2));
            AE_MUL2P32X4S(Y2, Y3, Y2, Y3, AE_MOVDA32(2), AE_MOVDA32(2));
            /* Q0 <- Q16 - 16 w/ rounding and saturation */
            y0 = AE_ROUND16X4F32SASYM(Y0, Y1);
            y1 = AE_ROUND16X4F32SASYM(Y2, Y3);
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

        {
            AE_L16X4X2_IP(c0, c1, pC, sizeof(ae_int16x8));

            /* Q15 <- Q15*Q0 */
            AE_MULZAAAA2Q16(ACC0, ACC1, c0, c0, x00, x10);  AE_MULAAAA2Q16(ACC0, ACC1, c1, c1, x01, x11);
            AE_MULZAAAA2Q16(ACC2, ACC3, c0, c0, x20, x30);  AE_MULAAAA2Q16(ACC2, ACC3, c1, c1, x21, x31);
            AE_MULZAAAA2Q16(ACC4, ACC5, c0, c0, x40, x50);  AE_MULAAAA2Q16(ACC4, ACC5, c1, c1, x41, x51);
            AE_MULZAAAA2Q16(ACC6, ACC7, c0, c0, x60, x70);  AE_MULAAAA2Q16(ACC6, ACC7, c1, c1, x61, x71);
            /* Q15 <- Q15 w/ saturation (convert 64 to 32-bit) */
            Y0 = AE_TRUNCA32X2F64S(ACC0, ACC1, 32);
            Y1 = AE_TRUNCA32X2F64S(ACC2, ACC3, 32);
            Y2 = AE_TRUNCA32X2F64S(ACC4, ACC5, 32);
            Y3 = AE_TRUNCA32X2F64S(ACC6, ACC7, 32);
            /* Remove DC bias */
            Y0 = AE_SUB32S(Y0, AE_MOVDA32X2((8*128)<<16, 0));
            /* Q0 <- Q15-1 - 15 w/ rounding and saturation */
            y0 = AE_ROUND16X4F32SASYM(Y0, Y1);
            y1 = AE_ROUND16X4F32SASYM(Y2, Y3);
            AE_S16X4X2_IP(y0, y1, pY, sizeof(ae_int16x8));
        }
        __Pragma("no_unroll");
        for (j=1; j<8; j++)
        {
            AE_L16X4X2_IP(c0, c1, pC, sizeof(ae_int16x8));

            /* Q15 <- Q15*Q0 */
            AE_MULZAAAA2Q16(ACC0, ACC1, c0, c0, x00, x10);  AE_MULAAAA2Q16(ACC0, ACC1, c1, c1, x01, x11);
            AE_MULZAAAA2Q16(ACC2, ACC3, c0, c0, x20, x30);  AE_MULAAAA2Q16(ACC2, ACC3, c1, c1, x21, x31);
            AE_MULZAAAA2Q16(ACC4, ACC5, c0, c0, x40, x50);  AE_MULAAAA2Q16(ACC4, ACC5, c1, c1, x41, x51);
            AE_MULZAAAA2Q16(ACC6, ACC7, c0, c0, x60, x70);  AE_MULAAAA2Q16(ACC6, ACC7, c1, c1, x61, x71);
            /* Q15 <- Q15 w/ saturation (convert 64 to 32-bit) */
            Y0 = AE_TRUNCA32X2F64S(ACC0, ACC1, 32);
            Y1 = AE_TRUNCA32X2F64S(ACC2, ACC3, 32);
            Y2 = AE_TRUNCA32X2F64S(ACC4, ACC5, 32);
            Y3 = AE_TRUNCA32X2F64S(ACC6, ACC7, 32);
            /* Q0 <- Q15-1 - 15 w/ rounding and saturation */
            y0 = AE_ROUND16X4F32SASYM(Y0, Y1);
            y1 = AE_ROUND16X4F32SASYM(Y2, Y3);
            AE_S16X4X2_IP(y0, y1, pY, sizeof(ae_int16x8));
        }
    }

    return 0;
} /* dct2d_8x16() */
#endif
