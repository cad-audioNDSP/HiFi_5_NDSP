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
    Complex data circular convolution, 32x16-bit
    C code optimized for HiFi5
  IntegrIT, 2006-2019
*/

/* Portable data types. */
#include "NatureDSP_types.h"
/* Signal Processing Library API. */
#include "NatureDSP_Signal_fir.h"
/* Common utility and macros declarations. */
#include "common.h"

/* Align address on a specified boundary. */
#define ALIGNED_ADDR( addr, align ) \
  (void*)(((uintptr_t)(addr)+((align)-1)) & ~((align)-1))

/*-------------------------------------------------------------------------
  Circular Convolution
  Performs circular convolution between vectors x (of length N) and y (of 
  length M) resulting in vector r of length N.
  These functions implement the circular convolution algorithm with no 
  limitations on x and y vectors length and alignment at the cost of 
  increased processing complexity. In addition, this implementation variant
  requires scratch memory area.

  Precision: 
  16x16     16x16-bit data, 16-bit outputs
  32x16     32x16-bit data, 32-bit outputs 
  32x32     32x32-bit data, 32-bit outputs
  32x32ep   the same as above but using 72-bit accumulator for intermediate 
            computations
  f         floating point

  Input:
  s[]       scratch area, 
              FIR_CONVOLA16X16_SCRATCH_SIZE(N,M) or
              FIR_CONVOLA32X16_SCRATCH_SIZE(N,M) or
              CXFIR_CONVOLA32X16_SCRATCH_SIZE(N,M) or
              FIR_CONVOLA32X32_SCRATCH_SIZE(N,M) or
              FIR_CONVOLA32X32EP_SCRATCH_SIZE(N,M) or
              FIR_CONVOLAF_SCRATCH_SIZE(N,M) bytes

  x[N]      input data Q15, Q31 or floating point
  y[M]      input data Q15, Q31 or floating point
  N         length of x
  M         length of y
  Output:
  r[N]      output data, Q15, Q31 or floating point

  Restriction:
  x,y,r,s   should not overlap
  s         must be aligned on an 16-bytes boundary
  N,M       must be >0
  N >= M-1  minimum allowed length of vector x is the length of y minus one
-------------------------------------------------------------------------*/

void cxfir_convola32x16(void            * restrict s,
                        complex_fract32 * restrict r,
                  const complex_fract32 * restrict x,
                  const complex_fract16 * restrict y,
                  int N, int M)
{
    void    * s_ptr;
    int32_t * x_buf;
    int16_t * y_buf;

    const ae_int32x2 *          S;
          ae_int32x2 * restrict D;
    const ae_int16x4 *          S16;
          ae_int32x2 * restrict D32;
          ae_int16   * restrict D16_s;
    const ae_int32x2 *          X;
    const ae_int32x2 *          SH;
          ae_int16x4 *          Y;
          ae_int32x2 * restrict R;

    ae_valign S_va, D_va, R_va;

    ae_f64     q0r, q1r, q2r, q3r, q0i, q1i, q2i, q3i;
    ae_f32x2   x0, x1, x2, x3, x4, x5;
    ae_f32x2   x0_conj, x1_conj, x2_conj, x3_conj;
    ae_f16x4   y0, y1;
    ae_int32x2 y_tmp;
    ae_int32x2 t;
    ae_int16x4 p;
    ae_int32x2 tmp;
    int _M;
    int n, m;

    NASSERT(s);
    NASSERT(r);
    NASSERT(x);
    NASSERT(y);
    NASSERT_ALIGN(s, 16);
    NASSERT(M > 0 && N > 0);
    NASSERT(N >= M - 1);

    //----------------------------------------------------------------------------
    // Partition the scratch memory area.
    s_ptr = s;
    x_buf = (int32_t*)ALIGNED_ADDR(s_ptr, 8);
    s_ptr = x_buf + 2 * (M - 1 + N + 3);
    y_buf = (int16_t*)ALIGNED_ADDR(s_ptr, 8);
    s_ptr = y_buf + 2 * (M + 3);
    ASSERT((int8_t *)s_ptr - (int8_t *)s <= (int)CXFIR_CONVOLA32X16_SCRATCH_SIZE(N, M));

    Y = (ae_int16x4 *)(x_buf + 2 * (M - 1 + N + 3) - 6);
    p = AE_ZERO16();
    AE_S16X4_IP(p, Y, +8);
    AE_S16X4_IP(p, Y, +8);
    AE_S16X4_IP(p, Y, +8);
    //----------------------------------------------------------------------------
    // Copy x[N] data into the aligned buffer in a way that simplifies the
    // convolution calculation:
    //  x_buf[M-1+N+3] = { x[N-(M-1)]..x[N-1] x[0] x[1]..x[N-1] X X X }
    // Three X locations are reserved to allow for block-4 processing.
    _M = (M + 1) >> 1;
    //
    // Copy last M-1 entries of x{N].
    //
    S = (const ae_int32x2 *)(x + (N - (M - 1)));
    D = (ae_int32x2 *)x_buf;
    S_va = AE_LA64_PP(S);
    for (m = 0; m < M - (-M & 1); m++)
    {
        AE_LA32X2_IP(t, S_va, S);
        AE_S32X2_IP(t, D, +8);
    }
    //
    // Copy x[N].
    //
    S = (const ae_int32x2 *)x;
    D = (ae_int32x2 *)((int32_t *)D - 2 * (!(M & 1)));
    S_va = AE_LA64_PP(S);
    D_va = AE_ZALIGN64();
    for (n = 0; n < N; n++)
    {
        AE_LA32X2_IP(t, S_va, S);
        AE_SA32X2_IP(t, D_va, D);
    }
    AE_SA64POS_FP(D_va, D);

    //----------------------------------------------------------------------------
    // Copy reverted y[M] data into the aligned buffer and append 3 zeros:
    //  y_buf[M+3] = { y[M-1]..y[0] 0 0 0 }
    //
    // Copy y[M] in reverted order.
    //
    S16 = (const ae_int16x4 *)(((const int16_t*)y) + 2 * M - 1);
    D32 = (ae_int32x2 *)y_buf;
    S_va = AE_LA64_PP(S16);
    for (m = 0; m < _M; m++)
    {
        AE_LA16X4_RIP(p, S_va, S16);
        tmp = AE_MOVINT32X2_FROMF16X4(p);
        AE_S32X2_IP(tmp, D32, 8);
    }
    //
    // Append three zeros to allow for block-4 processing.
    //
    D16_s = (ae_int16 *)((int16_t *)D32 - 2 * (-M & 1));
    p = AE_ZERO16();
    AE_S16_0_IP(p, D16_s, +2);
    AE_S16_0_IP(p, D16_s, +2);

    //----------------------------------------------------------------------------
    // Compute (N&~3) convolution results.
    X = (const ae_int32x2 *)x_buf;
    R = (      ae_int32x2 *)r;
    R_va = AE_ZALIGN64();
    //
    // Process vector x data in 8-entries blocks.
    //
    for (n = 0; n < (N >> 2); n++)
    {
        // Load 4 left-most unprocessed x[] entries, the first is x[8*n-(M-1)].
        // Q31
        AE_L32X2X2_IP(x0, x1, castxcc(ae_int32x4, X), 2 * sizeof(complex_fract32));
        AE_L32X2X2_IP(x2, x3, castxcc(ae_int32x4, X), 2 * sizeof(complex_fract32));
        x0_conj = AE_MUL32JS(x0);
        x1_conj = AE_MUL32JS(x1);
        x2_conj = AE_MUL32JS(x2);
        x3_conj = AE_MUL32JS(x3);

        // Use the shuttle pointer when computing the convolution. Preserve the X
        // pointer for the next iteration.
        SH = X;

        //
        // Inner loop prologue: process first 4 y[] entries for 8 accumulators.
        //
        Y = (ae_int16x4 *)y_buf;

        //
        // Inner loop kernel: process 4 y[] entries for 8 accumulators. 12-entries
        // register delay line is updated similarly to the loop prologue.
        //
#if 0
        q0r = q1r = q2r = q3r = AE_ZERO64();
        q0i = q1i = q2i = q3i = AE_ZERO64();
        __Pragma("loop_count min=1");
        for (m = 0; m < _M; m++)
        {
            // Q31
            AE_L32X2_IP(x4, SH, sizeof(complex_fract32));
            AE_L32X2_IP(x5, SH, sizeof(complex_fract32));

            // Load y[(M-1)-4*(m+1)]..y[(M-1)-4*(m+1)-3].
            // Q15
            AE_L16X4_IP(y0, Y, 2 * sizeof(complex_fract16));

            AE_MULAFC32X16W_H(q0r, q0i, x0, y0);
            AE_MULAFC32X16W_L(q0r, q0i, x1, y0);
            AE_MULAFC32X16W_H(q1r, q1i, x1, y0);
            AE_MULAFC32X16W_L(q1r, q1i, x2, y0);
            AE_MULAFC32X16W_H(q2r, q2i, x2, y0);
            AE_MULAFC32X16W_L(q2r, q2i, x3, y0);
            AE_MULAFC32X16W_H(q3r, q3i, x3, y0);
            AE_MULAFC32X16W_L(q3r, q3i, x4, y0);

            // 4 x[] entries are done, shift them out of the registers.
            x0 = x2; x1 = x3; x2 = x4; x3 = x5;
    }
#else
        AE_L32X2_IP(y_tmp, castxcc(ae_int32x2, Y), 2 * sizeof(complex_fract16));
        y0 = AE_MOVF16X4_FROMF32X2(y_tmp);

        q0r = AE_MULZAAAAQ32X16(x0_conj, x1_conj, y0);
        q0i = AE_MULZAAAAQ32X16(x0, x1, y0);
        q1r = AE_MULZAAAAQ32X16(x1_conj, x2_conj, y0);
        q1i = AE_MULZAAAAQ32X16(x1, x2, y0);
        q2r = q3r = AE_ZERO64();
        q2i = q3i = AE_ZERO64();

        for (m = 0; m < _M - 1; m++)
        {
            AE_L32X2X2_IP(x4, x5, castxcc(ae_int32x4, SH), 2 * sizeof(complex_fract32));

            AE_L32X2_IP(y_tmp, castxcc(ae_int32x2, Y), 2 * sizeof(complex_fract16));
            y1 = AE_MOVF16X4_FROMF32X2(y_tmp);

            AE_MULAAAA2Q32X16(q2r, q0r, x2_conj, x3_conj, y0, y1);
            AE_MULAAAA2Q32X16(q2i, q0i, x2, x3, y0, y1);
            AE_MULAAAA2Q32X16(q3r, q1r, x3_conj, AE_MUL32JS(x4), y0, y1);
            AE_MULAAAA2Q32X16(q3i, q1i, x3, x4, y0, y1);

            x2_conj = AE_MUL32JS(x4);
            x3_conj = AE_MUL32JS(x5);
            x2 = x4; x3 = x5;

            y0 = y1;
        }

        AE_L32X2_IP(x4, SH, sizeof(complex_fract32));
        AE_MULAAAAQ32X16(q2r, x2_conj, x3_conj, y0);
        AE_MULAAAAQ32X16(q2i, x2, x3, y0);
        AE_MULAAAAQ32X16(q3r, x3_conj, AE_MUL32JS(x4), y0);
        AE_MULAAAAQ32X16(q3i, x3, x4, y0);

        q0r = AE_SLAI64(q0r, 1); q0i = AE_SLAI64(q0i, 1);
        q1r = AE_SLAI64(q1r, 1); q1i = AE_SLAI64(q1i, 1);
        q2r = AE_SLAI64(q2r, 1); q2i = AE_SLAI64(q2i, 1);
        q3r = AE_SLAI64(q3r, 1); q3i = AE_SLAI64(q3i, 1);
#endif

        //
        // Convert and save 8 convolution results.
        //

        // Q31 <- Q16.47 - 16 w/ rounding and saturation.
        AE_SA32X2_IP(AE_ROUND32X2F48SASYM(q0r, q0i), R_va, R);
        AE_SA32X2_IP(AE_ROUND32X2F48SASYM(q1r, q1i), R_va, R);
        AE_SA32X2_IP(AE_ROUND32X2F48SASYM(q2r, q2i), R_va, R);
        AE_SA32X2_IP(AE_ROUND32X2F48SASYM(q3r, q3i), R_va, R);
    }

    N &= 3;
    if (N)
    {
        // Load 4 left-most unprocessed x[] entries, the first is x[8*n-(M-1)].
        // Q31
        AE_L32X2_IP(x0, X, +8);
        AE_L32X2_IP(x1, X, +8);
        AE_L32X2_IP(x2, X, +8);
        AE_L32X2_IP(x3, X, +8);

        // Use the shuttle pointer when computing the convolution. Preserve the X
        // pointer for the next iteration.
        SH = X;

        //
        // Inner loop prologue: process first 4 y[] entries for 8 accumulators.
        //
        Y = (ae_int16x4 *)y_buf;

        // Load y[(M-1)-4*(m+1)]..y[(M-1)-4*(m+1)-3].
        // Q15
        AE_L16X4_IP(y0, Y, 2 * sizeof(complex_fract16));

        AE_MULFC32X16W_H(q0r, q0i, x0, y0);
        AE_MULAFC32X16W_L(q0r, q0i, x1, y0);
        AE_MULFC32X16W_H(q1r, q1i, x1, y0);
        AE_MULAFC32X16W_L(q1r, q1i, x2, y0);
        AE_MULFC32X16W_H(q2r, q2i, x2, y0);
        AE_MULAFC32X16W_L(q2r, q2i, x3, y0);

        // 4 x[] entries are done, shift them out of the registers.
        x0 = x2; x1 = x3;

        //
        // Inner loop kernel: process 4 y[] entries for 8 accumulators. 12-entries
        // register delay line is updated similarly to the loop prologue.
        //
        for (m = 0; m < _M - 1; m++)
        {
            // Q31
            AE_L32X2_IP(x2, SH, sizeof(complex_fract32));
            AE_L32X2_IP(x3, SH, sizeof(complex_fract32));

            // Load y[(M-1)-4*(m+1)]..y[(M-1)-4*(m+1)-3].
            // Q15
            AE_L16X4_IP(y0, Y, 2 * sizeof(complex_fract16));

            AE_MULAFC32X16W_H(q0r, q0i, x0, y0);
            AE_MULAFC32X16W_L(q0r, q0i, x1, y0);
            AE_MULAFC32X16W_H(q1r, q1i, x1, y0);
            AE_MULAFC32X16W_L(q1r, q1i, x2, y0);
            AE_MULAFC32X16W_H(q2r, q2i, x2, y0);
            AE_MULAFC32X16W_L(q2r, q2i, x3, y0);

            // 4 x[] entries are done, shift them out of the registers.
            x0 = x2; x1 = x3;
        }

        AE_SA32X2_IP(AE_ROUND32X2F48SASYM(q0r, q0i), R_va, R);
        if (N > 1) AE_SA32X2_IP(AE_ROUND32X2F48SASYM(q1r, q1i), R_va, R);
        if (N > 2) AE_SA32X2_IP(AE_ROUND32X2F48SASYM(q2r, q2i), R_va, R);
    }
    AE_SA64POS_FP(R_va, R);
} /* cxfir_convola32x16() */
