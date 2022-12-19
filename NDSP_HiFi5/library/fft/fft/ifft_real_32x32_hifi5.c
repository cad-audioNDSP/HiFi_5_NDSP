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
    Real-valued inverse FFT: 32-bit data, 32-bit twiddle factors
    C code optimized for HiFi4
  IntegrIT, 2006-2019
*/

/* Signal Processing Library API. */
#include "NatureDSP_Signal_fft.h"
/* Common utility macros. */
#include "common.h"
/* Real-valued FFT descriptor structure. */
#include "fft_twiddles32x32.h"

#define _CONJ32(_x) {_x = AE_SEL32_HL(_x, AE_NEG32S(_x) ); }
#define SWAP_PTR(_x, _y) {int32_t *tmp = _x; _x = _y ; _y = tmp; } 


/*-------------------------------------------------------------------------
  Inverse FFT on Real Data
  These functions make inverse FFT on half spectral data forming real
  data samples.
      Scaling: 
      +-------------------+----------------------------------------+
      |      Function     |           Scaling options              |
      +-------------------+----------------------------------------+
      |  ifft_real16x16   |  2 - 16-bit dynamic scaling            | 
      |                   |  3 - fixed scaling before each stage   | 
      |  ifft_real32x32   |  2 - 32-bit dynamic scaling            | 
      |                   |  3 - fixed scaling before each stage   | 
      |  ifft_real32x16   |  2 - 32-bit dynamic scaling            | 
      |                   |  3 - fixed scaling before each stage   | 
      +-------------------+----------------------------------------+

  NOTES:
  1. Bit-reversing reordering is done here. 
  2. IFFT runs in-place algorithm so INPUT DATA WILL APPEAR DAMAGED after
     the call.
  3. Inverse FFT function for real signal transforms the input spectrum  
     and then calls ifft_cplx() with FFT size set to N/2. 
  4. 32x32,32x16,16x16  FFTs support mixed radix transforms.
  5. N - FFT size

  Precision:
  32x32  32-bit input/outputs, 32-bit twiddles
  32x16  32-bit input/outputs, 16-bit twiddles
  16x16  16-bit input/outputs, 16-bit twiddles

  Input:
  x[(N/2+1)*2]	input spectrum. Real and imaginary data are interleaved  
                and real data goes first. The imaginary part of 0th and
                N/2th input samples should be equal to zero.
  scalingOpt	scaling option (see table above)

  Output:			
  y[N]	        real output signal

  Returned value: total number of right shifts occurred during scaling 
                  procedure

  Restrictions:
  x,y           should not overlap
  x,y           aligned on a 16-bytes boundary
  x[(0)*2+1],
  x[(N/2)*2+1]  should be equal to zero
-------------------------------------------------------------------------*/
int ifft_real32x32( int32_t * y, int32_t * x, fft_handle_t h, int scalingOpt )
{
    int shift, shiftSum;
    int n, N, N4, bexp;
    fft_real32x32_descr_t *hr = (fft_real32x32_descr_t *)h;
    fft_cplx32x32_descr_t *hc = (fft_cplx32x32_descr_t *)hr->cfft_hdl;
    ae_int32x2 vA0, vA1, vB0, vB1, tw;
    ae_int32x2 * restrict p_x0,
               * restrict p_x1,
               * restrict p_y0,
               * restrict p_y1;
    ae_int32x2 * restrict ptw;
    const ae_int32x4 * restrict pX;
    ae_int32x2 scl; 
    int shiftl, shiftr; 
    ae_int16x4 tmp0, tmp1;
    ALIGN(16) const int16_t sel_tab[4] = { 0x703, 0x602, 0x105, 0x004 };
    ae_int16x4 sel = AE_L16X4_I((ae_int16x4*)sel_tab, 0);

    NASSERT(scalingOpt == 2 || scalingOpt == 3);
    NASSERT(x != y);
    NASSERT_ALIGN16(x);
    NASSERT_ALIGN16(y);
    N = 2*hc->N;
    N4 = (N+2) >> 2; /* Works for all even N */

    bexp = 0;
    if (scalingOpt == 2)
    {
        int numComplexSamples = (N>>1) + 1;
        ae_int32x2 x0, x1, x2, x3;
        ae_int16x4 nsa0, nsa1;
        int n;
        pX = (const ae_int32x4 *)x;
        nsa0 = 31; nsa1 = 31;
        
        __Pragma("loop_count min=1");
        for (n = 0; n<(numComplexSamples >> 2); n++)
        {
            /* 2 cycles unroll = 1 */
            AE_L32X2X2_IP(x0, x1, pX, sizeof(ae_int32x4));
            AE_L32X2X2_IP(x2, x3, pX, sizeof(ae_int32x4));
            nsa0 = AE_MIN16(nsa0, AE_NSA32X4(x0, x1));
            nsa1 = AE_MIN16(nsa1, AE_NSA32X4(x2, x3));
        }

        for (n = 0; n<(numComplexSamples & 3); n++)
        {
            AE_L32X2_IP(x0, castxcc(ae_int32x2, pX), sizeof(ae_int32x2));
            nsa0 = AE_MIN16(nsa0, AE_NSA32X4(x0, x0));
        }
        bexp = AE_RMIN16X4(AE_MIN16(nsa0, nsa1));
    }

    shiftSum = 0;
    AE_CALCRNG3();
    WUR_AE_SAR(0);

 /* 
    The real - to - complex spectrum conversion 
    MATLAB code:
    % N - Size of transform
    % X - input complex vector 1 x (N/2+1 )
    % x = real(N*ifft([X, conj(wrev(X(2:N/2))] ) )- output real vector 1xN
    
    twd = exp(-2*pi*1j*(0:N4-1)/N);
    a0 = X(1:N4);
    a1 = X(N/2+1:-1:N/2-N4+2); 
    b0 = a0+conj(a1);
    b1 = (a0-conj(a1))*1j.*conj(twd);
    a0 = b0+b1;
    a1 = conj(b0-b1);
    if(mod(N,4))
        x = [a0,  wrev(a1(2:N4))]; % N/2 complex samples
    else
        x = [a0,2*conj(X(N4+1)),wrev(a1(2:N4))]; % N/2 complex samples
    end

    tmp = N/2*ifft(x); 
    x = zeros(1, N);
    x(1:2:end) = real(tmp);
    x(2:2:end) = imag(tmp);
 */
    shift = (scalingOpt == 3) ? 2 : 2 - bexp;
    ASSERT(shift>-32 && shift<32);
    shiftl = XT_MAX(0, -shift);
    shiftr = XT_MAX(0, shift);
    scl = 1 << shiftl;
    WUR_AE_SAR(shiftr);

    p_x0 = (ae_int32x2 *)x;
    p_x1 = (ae_int32x2 *)x + N/2;
    p_y0 = (ae_int32x2 *)y;
    p_y1 = (ae_int32x2 *)y + N/2 - 1;
    ptw  = (ae_int32x2 *)hr->twd + 1;

    AE_L32X2_IP(vB0, p_x0,       sizeof(complex_fract32));
    AE_L32X2_XP(vB1, p_x1, -(int)sizeof(complex_fract32));

    vB0 = AE_SRAA32(vB0, shift);
    vB1 = AE_SRAA32(vB1, shift);
    AE_ADDANDSUB32S(vA0, vA1, vB0, vB1);

    vA0 = AE_SEL32_HH(vA0, vA1);
    AE_S32X2RNG_IP(vA0, p_y0, sizeof(complex_fract32)); 
    if ((N4 & 1) != 0 || (N & 3) != 0)
    {
        if (scalingOpt == 2)
        {
            /* Spectrum converter with prescaling */
            __Pragma("loop_count min=1"); 
            for (n = 1; n < N4; n++)
            {
                /*4 cycles unroll=1*/

                AE_L32X2_IP(tw, ptw, sizeof(complex_fract32));
                AE_L32X2_IP(vB0, p_x0, sizeof(complex_fract32));
                AE_L32X2_XP(vB1, p_x1, -(int)sizeof(complex_fract32));

                vB0 = AE_MULP32X2(vB0, scl);
                vB1 = AE_MULP32X2(vB1, scl);
                AE_ADDANDSUBRNG32(vA0, vA1, vB0, vB1);

                AE_DSEL16X4(tmp0, tmp1,
                            AE_MOVINT16X4_FROMINT32X2(vA0),
                            AE_MOVINT16X4_FROMINT32X2(vA1), sel);
                vB0 = AE_MOVINT32X2_FROMINT16X4(tmp0);
                vB1 = AE_MOVINT32X2_FROMINT16X4(tmp1);

                tw = AE_MUL32JS(tw);
                vB1 = AE_MULFCJ32RAS(tw, vB1);

                vA0 = AE_SUBADD32S(vB0, vB1);
                vA1 = AE_ADDSUB32S(vB1, vB0);

                AE_S32X2RNG_IP(vA0, p_y0, sizeof(complex_fract32));
                AE_S32X2RNG_XP(vA1, p_y1, -(int)sizeof(complex_fract32));
            }
        }
        else /* if (scalingOpt == 3) */
        {
            __Pragma("loop_count min=1");
            for (n = 1; n < N4; n++)
            {
                /*11 cycles unroll=4*/
                AE_L32X2_IP(tw, ptw, sizeof(complex_fract32));

                AE_L32X2_IP(vB0, p_x0, sizeof(complex_fract32));
                AE_L32X2_XP(vB1, p_x1, -(int)sizeof(complex_fract32));

                AE_ADDANDSUBRNG32(vA0, vA1, vB0, vB1);

                AE_DSEL16X4(tmp0, tmp1,
                            AE_MOVINT16X4_FROMINT32X2(vA0),
                            AE_MOVINT16X4_FROMINT32X2(vA1), sel);
                vB0 = AE_MOVINT32X2_FROMINT16X4(tmp0);
                vB1 = AE_MOVINT32X2_FROMINT16X4(tmp1);

                tw = AE_MUL32JS(tw);
                vB1 = AE_MULFCJ32RAS(tw, vB1);

                vA0 = AE_SUBADD32S(vB0, vB1);
                vA1 = AE_ADDSUB32S(vB1, vB0);

                AE_S32X2RNG_IP(vA0, p_y0, sizeof(complex_fract32));
                AE_S32X2RNG_XP(vA1, p_y1, -(int)sizeof(complex_fract32));
            }
        } /* if (scalingOpt == 3) else..  */
    } 
    else  /* if ((N4 & 1) != 0 || (N & 3) != 0) */
    {
        //n=1 for (n = 1; n < N4; n++)
        {
            AE_L32X2_IP(tw, ptw, sizeof(complex_fract32));

            AE_L32X2_IP(vB0, p_x0,       sizeof(complex_fract32));
            AE_L32X2_XP(vB1, p_x1, -(int)sizeof(complex_fract32));

            vB0 = AE_MULP32X2(vB0, scl);
            vB1 = AE_MULP32X2(vB1, scl);

            // ADD/SUBB
            AE_ADDANDSUBRNG32(vA0, vA1, vB0, vB1);

            /* vB0 = AE_SEL32_HL(vA0, vA1);
               vB1 = AE_SEL32_HL(vA1, vA0); */

            AE_DSEL16X4(tmp0, tmp1,
                AE_MOVINT16X4_FROMINT32X2(vA0),
                AE_MOVINT16X4_FROMINT32X2(vA1), sel);
            vB0 = AE_MOVINT32X2_FROMINT16X4(tmp0);
            vB1 = AE_MOVINT32X2_FROMINT16X4(tmp1);

            tw = AE_MUL32JS(tw);
            vB1 = AE_MULFCJ32RAS(tw, vB1);
        
            vA0 = AE_SUBADD32S(vB0, vB1);
            vA1 = AE_ADDSUB32S(vB1, vB0);

            AE_S32X2RNG_IP(vA0, p_y0,       sizeof(complex_fract32));
            AE_S32X2RNG_XP(vA1, p_y1, -(int)sizeof(complex_fract32));
        }
        if (scalingOpt == 2)
        {
            __Pragma("loop_count min=1");
            for (n = 2; n < N4; n+=2)
            {
                /* 15 cycles unroll = 2 */
                ae_int32x2 tw0, tw1;
                ae_int32x2 vA00, vA10, vA01, vA11, vB00, vB01, vB11, vB10;
                AE_L32X2X2_IP(vB00, vB01, castxcc(ae_int32x4, p_x0), 2 * sizeof(complex_fract32));
                AE_L32X2X2_IP(tw0,  tw1,  castxcc(ae_int32x4, ptw),  2 * sizeof(complex_fract32));

                AE_L32X2_XP(vB10, p_x1, -(int)sizeof(complex_fract32));
                AE_L32X2_XP(vB11, p_x1, -(int)sizeof(complex_fract32));

                vB00 = AE_MULP32X2(vB00, scl);
                vB10 = AE_MULP32X2(vB10, scl);
                vB01 = AE_MULP32X2(vB01, scl);
                vB11 = AE_MULP32X2(vB11, scl);

                AE_ADDANDSUBRNG32(vA00, vA10, vB00, vB10);
                AE_ADDANDSUBRNG32(vA01, vA11, vB01, vB11);

                AE_DSEL16X4(tmp0, tmp1,
                    AE_MOVINT16X4_FROMINT32X2(vA00),
                    AE_MOVINT16X4_FROMINT32X2(vA10), sel);
                vB00 = AE_MOVINT32X2_FROMINT16X4(tmp0);
                vB10 = AE_MOVINT32X2_FROMINT16X4(tmp1);

                AE_DSEL16X4(tmp0, tmp1,
                    AE_MOVINT16X4_FROMINT32X2(vA01),
                    AE_MOVINT16X4_FROMINT32X2(vA11), sel);
                vB01 = AE_MOVINT32X2_FROMINT16X4(tmp0);
                vB11 = AE_MOVINT32X2_FROMINT16X4(tmp1);

                tw0 = AE_MUL32JS(tw0);
                tw1 = AE_MUL32JS(tw1);

                vB10 = AE_MULFCJ32RAS(tw0, vB10);
                vB11 = AE_MULFCJ32RAS(tw1, vB11);

                vA00 = AE_SUBADD32S(vB00, vB10);
                vA10 = AE_ADDSUB32S(vB10, vB00);
                vA01 = AE_SUBADD32S(vB01, vB11);
                vA11 = AE_ADDSUB32S(vB11, vB01);
#if 0
                // Doesn't work! 
                AE_S32X2X2RNG_IP(vA00, vA01, castxcc(ae_int32x4, p_y0), 2 * sizeof(complex_fract32));
#else
                AE_S32X2RNG_IP(vA00, p_y0, sizeof(complex_fract32));
                AE_S32X2RNG_IP(vA01, p_y0, sizeof(complex_fract32));
#endif
                AE_S32X2RNG_XP(vA10, p_y1, -(int)sizeof(complex_fract32));
                AE_S32X2RNG_XP(vA11, p_y1, -(int)sizeof(complex_fract32));
            }
        }
        else /* if (scalingOpt == 2) */
        {
            
            __Pragma("loop_count min=1");
            for (n = 2; n < N4; n += 2)
            {
                /* 6 cycles unroll = 1 */
                ae_int32x2 tw0, tw1;
                ae_int32x2 vA00, vA10, vA01, vA11, vB00, vB01, vB11, vB10;
                AE_L32X2X2_IP(vB00, vB01, castxcc(ae_int32x4, p_x0), 2 * sizeof(complex_fract32));
                AE_L32X2X2_IP(tw0, tw1, castxcc(ae_int32x4, ptw), 2 * sizeof(complex_fract32));

                AE_L32X2_XP(vB10, p_x1, -(int)sizeof(complex_fract32));
                AE_L32X2_XP(vB11, p_x1, -(int)sizeof(complex_fract32));

                AE_ADDANDSUBRNG32(vA00, vA10, vB00, vB10);
                AE_ADDANDSUBRNG32(vA01, vA11, vB01, vB11);

                /*
                vB00 = AE_SEL32_HL(vA00, vA10);
                vB10 = AE_SEL32_HL(vA10, vA00);
                vB01 = AE_SEL32_HL(vA01, vA11);
                vB11 = AE_SEL32_HL(vA11, vA01);*/
                AE_DSEL16X4(tmp0, tmp1,
                    AE_MOVINT16X4_FROMINT32X2(vA00),
                    AE_MOVINT16X4_FROMINT32X2(vA10), sel);
                vB00 = AE_MOVINT32X2_FROMINT16X4(tmp0);
                vB10 = AE_MOVINT32X2_FROMINT16X4(tmp1);

                AE_DSEL16X4(tmp0, tmp1,
                    AE_MOVINT16X4_FROMINT32X2(vA01),
                    AE_MOVINT16X4_FROMINT32X2(vA11), sel);
                vB01 = AE_MOVINT32X2_FROMINT16X4(tmp0);
                vB11 = AE_MOVINT32X2_FROMINT16X4(tmp1);

                tw0 = AE_MUL32JS(tw0);
                tw1 = AE_MUL32JS(tw1);

                vB10 = AE_MULFCJ32RAS(tw0, vB10);
                vB11 = AE_MULFCJ32RAS(tw1, vB11);

                vA00 = AE_SUBADD32S(vB00, vB10);
                vA10 = AE_ADDSUB32S(vB10, vB00);
                vA01 = AE_SUBADD32S(vB01, vB11);
                vA11 = AE_ADDSUB32S(vB11, vB01);

                AE_S32X2X2RNG_IP(vA00, vA01, castxcc(ae_int32x4, p_y0), 2 * sizeof(complex_fract32));

                AE_S32X2RNG_XP(vA10, p_y1, -(int)sizeof(complex_fract32));
                AE_S32X2RNG_XP(vA11, p_y1, -(int)sizeof(complex_fract32));
            }
        } /* if (scalingOpt == 2) .. else ... */
    } /* if ((N4 & 1) != 0 || (N & 3) != 0) .. else ...  */

    /* When N is a multiple of 4 */
    if (0==(N&3))
    {
        /* 2*conj(x(N/4+1)) */
        vB0 = AE_L32X2_I(p_x0, 0);
        vB0 = AE_SRAA32(vB0, shift-1);
        _CONJ32(vB0); 
        AE_S32X2RNG_I(vB0, p_y0, 0);
    }
    shiftSum += shift;

    /*
     * half-sized complex-valued inverse FFT.
     */
    {
        int32_t *pdest = y; 
        const int  N = hc->N; 
        const int *tw_step = hc->tw_step;
        const cint32ptr_t *tw_tab = hc->twd;
        const fft_cplx32x32_stage_t *stg_fn = (scalingOpt == 2) ? hc->stages_s2 : hc->stages_s3;
        int v = 1;

        bexp = 0;
        if (scalingOpt==2)
        {
            AE_CALCRNG3();
            bexp = 3 - RUR_AE_SAR();
        }
        
        n = 0;
        SWAP_PTR(x, y);
        while ( stg_fn[n+1] )
        {
            shiftSum += stg_fn[n](tw_tab[n], x, y, N, &v, tw_step[n], &bexp); 
            SWAP_PTR(x, y); 
            n++; 
        }
        if (y != pdest)
        {  
            /* Execute the last stage inplace */
            y = x;
        }
        /* Last stage */
        shiftSum += stg_fn[n](tw_tab[n], x, y, N, &v, tw_step[n], &bexp);
    }

    return shiftSum;
} /* ifft_real32x32() */
