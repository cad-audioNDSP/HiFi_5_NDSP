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
  NatureDSP Signal Processing Library. IIR part
    Bi-quad Real Block IIR, 16x16-bit, Direct Form I
    Code optimized for HiFi4
  IntegrIT, 2006-2018
*/

/* Portable data types. */
#include "NatureDSP_types.h"
/* Signal Processing Library API. */
#include "NatureDSP_Signal_iir.h"
/* Common utility and macros declarations. */
#include "common.h"
/* Filter instance structure. */
#include "stereo_bqriir16x16_df1_nd_common.h"

/*-------------------------------------------------------------------------
  Bi-quad Real Block IIR
  Computes a real IIR filter (cascaded IIR direct form I or II using 5 
  coefficients per bi-quad + gain term). Real input data are stored
  in vector x. The filter output result is stored in vector r. The filter 
  calculates N output samples using SOS and G matrices.
  NOTE:
  1. Bi-quad coefficients may be derived from standard SOS and G matrices 
  generated by MATLAB. However, typically biquad stages have big peaks 
  in their step response which may cause undesirable overflows at the 
  intermediate outputs. To avoid that the additional scale factors 
  coef_g[M] may be applied. These per-section scale factors may require 
  some tuning to find a compromise between quantization noise and possible 
  overflows. Output of the last section is directed to an additional 
  multiplier, with the gain factor being a power of two, either negative 
  or non-negative. It is specified through the total gain shift amount 
  parameter gain of each filter initialization function.
  2. 16x16 filters may suffer more from accumulation of the roundoff errors,
  so filters should be properly designed to match noise requirements
  3. Due to the performance reasons, IIR biquad filters may introduce 
  additional algorithmic delay of several sampless. Amount of that delay
  might be requested by the  xxx_groupDelay API. For sensitive applications
  all the filters have delayless implementations (with  _nd  suffix in the name).
  Formally, the xxx_groupDelay APIs is also implemented for that kind of filters,
  but return zero.
  
  Precision: 
  16x16         16-bit data, 16-bit coefficients, 16-bit intermediate 
                stage outputs (DF1, DF1 stereo, DF II form)
  32x16         32-bit data, 16-bit coefficients, 32-bit intermediate 
                (DF1, DF1 stereo, DF II form) stage outputs
  32x32         32-bit data, 32-bit coefficients, 32-bit intermediate 
                (DF I, DF1 stereo,  DF II form) stage outputs 
  f             floating point (DF I, DF1 stereo, DF II and DF IIt)

  Input:
  N             length of input sample block
  M             number of bi-quad sections
  S             1 for mono, 2 for stereo API
  s[]           scratch memory area (for fixed-point functions only). 
                Minimum number of bytes depends on selected filter structure 
                and precision:
                  BQRIIR16X16_DF1_SCRATCH_SIZE( N,  M ), or
                  BQRIIR16X16_DF2_SCRATCH_SIZE( N,  M ), or
                  BQRIIR32X16_DF1_SCRATCH_SIZE( N,  M ), or
                  BQRIIR32X16_DF2_SCRATCH_SIZE( N,  M ), or
                  BQRIIR32X32_DF1_SCRATCH_SIZE( N,  M ), or
                  BQRIIR32X32_DF2_SCRATCH_SIZE( N,  M ),
                  STEREO_BQRIIR16X16_DF1_SCRATCH_SIZE( N, M ) , or
                  STEREO_BQRIIR32X32_DF1_SCRATCH_SIZE( N, M ) , or
                  STEREO_BQRIIRF_DF1_SCRATCH_SIZE    ( N, M )
                 If a particular macro returns zero, then the corresponding
                 IIR doesn't require a scratch area and parameter s may 
                 hold zero.
  coef_sos[M*5]  filter coefficients stored in blocks of 5 numbers: 
                 b0 b1 b2 a1 a2. 
                 For fixed-point funcions, fixed point format of filter 
                 coefficients is Q1.14 for 16x16 and 32x16, or Q1.30 for 32x32 
  coef_sosl[M*5] filter coefficients for the left channel (stereo filters only)
  coef_sosr[M*5] filter coefficients for the left channel (stereo filters only)
  coef_g[M]      scale factor for each section, Q15 (for fixed-point 
                 functions only)
  coef_gl[M]     scale factor for the left channel (stereo filters only)
  coef_gr[M]     scale factor for the right channel (stereo filters only)
  gain           total gain shift amount applied to output signal of the
                 last section, -48..15
  gainl          total gain shift amount  for the left channel (stereo filters 
                 only)
  gainr          total gain shift amount for the right channel (stereo filters 
                 only)

  x[N*S]         input samples, Q15, Q31 or floating point. Stereo samples 
                 go in interleaved order (left, right)
  Output:
  r[N*S]         output data, Q15, Q31 or floating point. Stereo samples go 
                 in interleaved order (left, right) 

  Restriction:
  x,r,s,coef_g,coef_sos should not overlap
  N   - must be a multiple of 2
  s[] - whenever supplied must be aligned on an 16-bytes boundary
-------------------------------------------------------------------------*/

#define sz_i16 sizeof(int16_t)
#define MAX_BUFFER_SZ ((int)(MAX_ALLOCA_SZ/(sizeof(int16_t)*2)))

void stereo_bqriir16x16_df1_nd(stereo_bqriir16x16_df1_nd_handle_t _bqriir,
                             void    * restrict       s,
                             int16_t * restrict       r,
                       const int16_t *                x, int N )
{
    stereo_bqriir16x16_df1_nd_ptr_t bqriir = (stereo_bqriir16x16_df1_nd_ptr_t)_bqriir;

    const ae_int16x4 * restrict X;
    ae_int16x4 * restrict R;
    ae_int16x8 * restrict St;
    const ae_int16x8 * restrict Cf;

    ae_int16x4 inx, y0, y1;
    ae_int16x4 ALsxsr, ALsx, ALsr;
    ae_int16x4 ARsxsr, ARsx, ARsr;
    ae_int16x4 BLsxsr, BLsx, BLsr;
    ae_int16x4 BRsxsr, BRsx, BRsr;
    ae_int16x4 ALb0b1b2_0, AL_0b0b1b2, ARb0b1b2_0, AR_0b0b1b2;
    ae_int16x4 BLb0_0b2b1, BLb1b0_0b2, BRb0_0b2b1, BRb1b0_0b2;
    ae_int16x4 AL_0_0a2a1, AR_0_0a2a1, BL_0_0a2a1, BR_0_0a2a1;
    ae_int16x4 ALgARg, BLgBRg;
    ae_int16x4 cf0, cf1, cf2, cf3;
    ae_int16x4 _0;
    ae_int32x2 ty0, ty1;
    ae_int64   accL0, accL1, accR0, accR1;
    ae_valign al_X, al_R;
    int M, gainL, gainR;
    int n, m;

    ASSERT( bqriir && bqriir->magic == STEREO_BQRIIR16X16_DF1_ND_MAGIC && r && x );
    NASSERT(N%2==0);
    if(N<=0) return;

    M = bqriir->M;
    St = (ae_int16x8 *)bqriir->state;
    Cf = (ae_int16x8 *)bqriir->coef;
    X = (const ae_int16x4 *)x;
    gainL = bqriir->gainL;
    gainR = bqriir->gainR;
    _0 = AE_ZERO16();
    al_R = AE_ZALIGN64();

    // Perform data block processing for each section, by 2 sections simulateously.
    // Use the output array r[N] for temporal storage of inter-section signal.
    for ( m=0; m<((M-1)>>1); m++ )
    {
        //
        // Load 16-bit section coefficients.
        //
        AE_L16X4X2_IP(cf0, cf1, Cf, 8*sz_i16);
        AE_L16X4X2_IP(cf2, cf3, Cf, 8*sz_i16);
        ALb0b1b2_0 = cf1;
        AL_0b0b1b2 = AE_SEL16_4321(cf1, cf1);
        AL_0_0a2a1 = AE_SEL16_5410(_0, cf0);
        ARb0b1b2_0 = cf3;
        AR_0b0b1b2 = AE_SEL16_4321(cf3, cf3);
        AR_0_0a2a1 = AE_SEL16_5410(_0, cf2);
        ALgARg = AE_SEL16_7632(cf0, cf2);

        AE_L16X4X2_IP(cf0, cf1, Cf, 8*sz_i16);
        AE_L16X4X2_IP(cf2, cf3, Cf, 8*sz_i16);
        cf1 = AE_SHORTSWAP(cf1);
        BLb0_0b2b1 = AE_SEL16_4321(cf1, cf1);
        BLb1b0_0b2 = AE_SEL16_5432(cf1, cf1);
        BL_0_0a2a1 = AE_SEL16_5410(_0, cf0);
        cf3 = AE_SHORTSWAP(cf3);
        BRb0_0b2b1 = AE_SEL16_4321(cf3, cf3);
        BRb1b0_0b2 = AE_SEL16_5432(cf3, cf3);
        BR_0_0a2a1 = AE_SEL16_5410(_0, cf2);
        BLgBRg = AE_SEL16_7632(cf0, cf2);

        //
        // Load 16-bit section's state elements.
        //
        AE_L16X4X2_I(ALsxsr, ARsxsr, St, 0);
        AE_L16X4X2_I(BLsxsr, BRsxsr, St, 8*sz_i16);
        ALsx = ALsxsr;
        ALsr = ALsxsr;
        ARsx = ARsxsr;
        ARsr = ARsxsr;
        BLsx = BLsxsr;
        BLsr = BLsxsr;
        BRsx = BRsxsr;
        BRsr = BRsxsr;

        R = (ae_int16x4 *)r;

        //
        // Pass the block of input samples through 2 sections. n-th sample at the
        // output of a biquad section:
        //   r[n] = g*x[n]*b0 + g*x[n-1]*b1 + g*x[n-2]*b2 - r[n-1]*a1 - r[n-2]*a2
        //
        al_X = AE_LA64_PP(X);
        for ( n=0; n<(N>>1); n++ )
        {
            // Load 2 input samples.
            AE_LA16X4_IP(inx, al_X, X);
            // deinterleave samples
            inx = AE_SEL16I(inx, inx, 14); // SEL pattern: 5342
                                           // inx <- g*x[n]
                                           // y0  <- g*dy0, dy0 - delayed output from previous section.
                                           // Q15 <- Q15*Q15 - 15 w/ rounding. No saturation is necessary.
            inx = AE_MULFP16X4RAS(inx, ALgARg);
            ALsx = AE_SEL16_7632(inx, ALsx);
            ARsx = AE_SEL16_5432(inx, ARsx);

            // y[n] = g*x[n]*b0 + g*x[n-1]*b1 + g*x[n-2]*b2
            // y[n] = y[n] - r[n-1]*a1 - r[n-2]*a2
            // Q29 <- Q15*Q14
            // Left Section 0 / sample 0
            // Right Section 0 / sample 0
            AE_MULZAAAA2Q16(accL0, accR0, ALsx, ARsx, AL_0b0b1b2, AR_0b0b1b2);
            AE_MULAAAA2Q16(accL0, accR0, ALsr, ARsr, AL_0_0a2a1, AR_0_0a2a1);
            AE_PKSR16(ALsr, accL0, 2);
            AE_PKSR16(ARsr, accR0, 2);

            // Left Section 0 / sample 1
            // Right Section 0 / sample 1
            AE_MULZAAAA2Q16(accL1, accR1, ALsx, ARsx, ALb0b1b2_0, ARb0b1b2_0);
            AE_MULAAAA2Q16(accL1, accR1, ALsr, ARsr, AL_0_0a2a1, AR_0_0a2a1);
            AE_PKSR16(ALsr, accL1, 2);
            AE_PKSR16(ARsr, accR1, 2);

            // pass computed samples from the 1st section
            // to the input of the 2nd section.
            y0 = ALsr;
            y1 = ARsr;
            y0 = AE_SEL16_5410(y0, y1);
            y0 = AE_MULFP16X4RAS(y0, BLgBRg);
            BLsx = AE_SEL16_7632(y0, BLsx);
            BRsx = AE_SEL16_5432(y0, BRsx);

            // Left Section 1 / sample 0
            // Right Section 1 / sample 0
            AE_MULZAAAA2Q16(accL0, accR0, BLsx, BRsx, BLb0_0b2b1, BRb0_0b2b1);
            AE_MULAAAA2Q16(accL0, accR0, BLsr, BRsr, BL_0_0a2a1, BR_0_0a2a1);
            AE_PKSR16(BLsr, accL0, 2);
            AE_PKSR16(BRsr, accR0, 2);

            // Left Section 1 / sample 1
            // Right Section 1 / sample 1
            AE_MULZAAAA2Q16(accL1, accR1, BLsx, BRsx, BLb1b0_0b2, BRb1b0_0b2);
            AE_MULAAAA2Q16(accL1, accR1, BLsr, BRsr, BL_0_0a2a1, BR_0_0a2a1);
            AE_PKSR16(BLsr, accL1, 2);
            AE_PKSR16(BRsr, accR1, 2);

            y0 = BLsr;
            y1 = BRsr;

            // Save output
            y0 = AE_SEL16_5140(y1, y0);
            y0 = AE_SEL16_2301(y0, y0);
            AE_SA16X4_IP(y0, al_R, R);
        }
        AE_SA64POS_FP(al_R, R);

        //
        // Save section's state elements.
        //
        ALsxsr = AE_SEL16_7610(ALsx, ALsr);
        ARsxsr = AE_SEL16_7610(ARsx, ARsr);
        BLsxsr = AE_SEL16_7610(BLsx, BLsr);
        BRsxsr = AE_SEL16_7610(BRsx, BRsr);
        AE_S16X4X2_IP(ALsxsr, ARsxsr, St, 8*sz_i16);
        AE_S16X4X2_IP(BLsxsr, BRsxsr, St, 8*sz_i16);

        //
        // Next sections are fed with output samples of the preceding biquad.
        //
        X = (const ae_int16x4 *)r;
    }
    //------------------------------------------------------------------------
    // Pass signal through the last biquads (1 or 2) and apply
    // the total gain.
    //------------------------------------------------------------------------
    M = (M-1)&1;

    if (M == 1) // Last 2 sections
    {
        //
        // Load 16-bit section coefficients.
        //
        AE_L16X4X2_IP(cf0, cf1, Cf, 8*sz_i16);
        AE_L16X4X2_IP(cf2, cf3, Cf, 8*sz_i16);
        ALb0b1b2_0 = cf1;
        AL_0b0b1b2 = AE_SEL16_4321(cf1, cf1);
        AL_0_0a2a1 = AE_SEL16_5410(_0, cf0);
        ARb0b1b2_0 = cf3;
        AR_0b0b1b2 = AE_SEL16_4321(cf3, cf3);
        AR_0_0a2a1 = AE_SEL16_5410(_0, cf2);
        ALgARg = AE_SEL16_7632(cf0, cf2);

        AE_L16X4X2_IP(cf0, cf1, Cf, 8*sz_i16);
        AE_L16X4X2_IP(cf2, cf3, Cf, 8*sz_i16);
        cf1 = AE_SHORTSWAP(cf1);
        BLb0_0b2b1 = AE_SEL16_4321(cf1, cf1);
        BLb1b0_0b2 = AE_SEL16_5432(cf1, cf1);
        BL_0_0a2a1 = AE_SEL16_5410(_0, cf0);
        cf3 = AE_SHORTSWAP(cf3);
        BRb0_0b2b1 = AE_SEL16_4321(cf3, cf3);
        BRb1b0_0b2 = AE_SEL16_5432(cf3, cf3);
        BR_0_0a2a1 = AE_SEL16_5410(_0, cf2);
        BLgBRg = AE_SEL16_7632(cf0, cf2);

        //
        // Load 16-bit section's state elements.
        //
        AE_L16X4X2_I(ALsxsr, ARsxsr, St, 0);
        AE_L16X4X2_I(BLsxsr, BRsxsr, St, 8*sz_i16);
        ALsx = ALsxsr;
        ALsr = ALsxsr;
        ARsx = ARsxsr;
        ARsr = ARsxsr;
        BLsx = BLsxsr;
        BLsr = BLsxsr;
        BRsx = BRsxsr;
        BRsr = BRsxsr;

        R = (ae_int16x4 *)r;

        //
        // Pass the block of input samples through 4 sections. n-th sample at the
        // output of a biquad section:
        //   r[n] = g*x[n]*b0 + g*x[n-1]*b1 + g*x[n-2]*b2 - r[n-1]*a1 - r[n-2]*a2
        //
        al_X = AE_LA64_PP(X);
        __Pragma("loop_count min=1");
        for ( n=0; n<(N>>1); n++ )
        {
            // Load 2 input samples.
            AE_LA16X4_IP(inx, al_X, X);
            // deinterleave samples
            inx = AE_SEL16I(inx, inx, 14); // SEL pattern: 5342
                                           // inx <- g*x[n]
                                           // y0  <- g*dy0, dy0 - delayed output from previous section.
                                           // Q15 <- Q15*Q15 - 15 w/ rounding. No saturation is necessary.
            inx = AE_MULFP16X4RAS(inx, ALgARg);
            ALsx = AE_SEL16_7632(inx, ALsx);
            ARsx = AE_SEL16_5432(inx, ARsx);

            // y[n] = g*x[n]*b0 + g*x[n-1]*b1 + g*x[n-2]*b2
            // y[n] = y[n] - r[n-1]*a1 - r[n-2]*a2
            // Q29 <- Q15*Q14

            // Left Section 0 / sample 0
            // Right Section 0 / sample 0
            AE_MULZAAAA2Q16(accL0, accR0, ALsx, ARsx, AL_0b0b1b2, AR_0b0b1b2);
            AE_MULAAAA2Q16(accL0, accR0, ALsr, ARsr, AL_0_0a2a1, AR_0_0a2a1);
            AE_PKSR16(ALsr, accL0, 2);
            AE_PKSR16(ARsr, accR0, 2);


            // Left Section 0 / sample 1
            // Right Section 0 / sample 1
            AE_MULZAAAA2Q16(accL1, accR1, ALsx, ARsx, ALb0b1b2_0, ARb0b1b2_0);
            AE_MULAAAA2Q16(accL1, accR1, ALsr, ARsr, AL_0_0a2a1, AR_0_0a2a1);
            AE_PKSR16(ALsr, accL1, 2);
            AE_PKSR16(ARsr, accR1, 2);

            // pass computed samples from the 1st section
            // to the input of the 2nd section.
            y0 = ALsr;
            y1 = ARsr;
            y0 = AE_SEL16_5410(y0, y1);
            y0 = AE_MULFP16X4RAS(y0, BLgBRg);
            BLsx = AE_SEL16_7632(y0, BLsx);
            BRsx = AE_SEL16_5432(y0, BRsx);

            // Left Section 1 / sample 0
            // Right Section 1 / sample 0
            AE_MULZAAAA2Q16(accL0, accR0, BLsx, BRsx, BLb0_0b2b1, BRb0_0b2b1);
            AE_MULAAAA2Q16(accL0, accR0, BLsr, BRsr, BL_0_0a2a1, BR_0_0a2a1);
            AE_PKSR16(BLsr, accL0, 2);
            AE_PKSR16(BRsr, accR0, 2);

            // Left Section 1 / sample 1
            // Right Section 1 / sample 1
            AE_MULZAAAA2Q16(accL1, accR1, BLsx, BRsx, BLb1b0_0b2, BRb1b0_0b2);
            AE_MULAAAA2Q16(accL1, accR1, BLsr, BRsr, BL_0_0a2a1, BR_0_0a2a1);
            AE_PKSR16(BLsr, accL1, 2);
            AE_PKSR16(BRsr, accR1, 2);

            // Apply the total gain shift and format outputs y[n], y[n+1]
            // Q(15+gain) <- Q(29 - 14 + gain) w/ rounding and saturation
            ty0 = AE_TRUNCA32X2F64S(accL0, accL1, 32+2+gainL);
            ty1 = AE_TRUNCA32X2F64S(accR0, accR1, 32+2+gainR);
            y0 = AE_ROUND16X4F32SASYM(ty0, ty1);
            // Save output
            y0 = AE_SEL16_7520(y0, y0);
            AE_SA16X4_IP(y0, al_R, R);
        }
        AE_SA64POS_FP(al_R, R);

        //
        // Save section's state elements.
        //
        ALsxsr = AE_SEL16_7610(ALsx, ALsr);
        ARsxsr = AE_SEL16_7610(ARsx, ARsr);
        BLsxsr = AE_SEL16_7610(BLsx, BLsr);
        BRsxsr = AE_SEL16_7610(BRsx, BRsr);
        AE_S16X4X2_IP(ALsxsr, ARsxsr, St, 8*sz_i16);
        AE_S16X4X2_IP(BLsxsr, BRsxsr, St, 8*sz_i16);
    }
    else // Last 1 section
    {
        //
        // Load 16-bit section coefficients.
        //
        AE_L16X4X2_IP(cf0, cf1, Cf, 8*sz_i16);
        AE_L16X4X2_IP(cf2, cf3, Cf, 8*sz_i16);
        ALb0b1b2_0 = cf1;
        AL_0b0b1b2 = AE_SEL16_4321(cf1, cf1);
        AL_0_0a2a1 = AE_SEL16_5410(_0, cf0);
        ARb0b1b2_0 = cf3;
        AR_0b0b1b2 = AE_SEL16_4321(cf3, cf3);
        AR_0_0a2a1 = AE_SEL16_5410(_0, cf2);
        ALgARg = AE_SEL16_7632(cf0, cf2);

        //
        // Load 16-bit section's state elements.
        //
        AE_L16X4X2_I(ALsxsr, ARsxsr, St, 0);
        ALsx = ALsxsr;
        ALsr = ALsxsr;
        ARsx = ARsxsr;
        ARsr = ARsxsr;

        R = (ae_int16x4 *)r;

        //
        // Pass the block of input samples through 4 sections. n-th sample at the
        // output of a biquad section:
        //   r[n] = g*x[n]*b0 + g*x[n-1]*b1 + g*x[n-2]*b2 - r[n-1]*a1 - r[n-2]*a2
        //
        al_X = AE_LA64_PP(X);
        __Pragma("loop_count min=1");
        for ( n=0; n<(N>>1); n++ )
        {
            // Load 2 input samples.
            AE_LA16X4_IP(inx, al_X, X);
            // deinterleave samples
            inx = AE_SEL16I(inx, inx, 14); // SEL pattern: 5342
                                           // inx <- g*x[n]
                                           // y0  <- g*dy0, dy0 - delayed output from previous section.
                                           // Q15 <- Q15*Q15 - 15 w/ rounding. No saturation is necessary.
            inx = AE_MULFP16X4RAS(inx, ALgARg);
            ALsx = AE_SEL16_7632(inx, ALsx);
            ARsx = AE_SEL16_5432(inx, ARsx);

            // y[n] = g*x[n]*b0 + g*x[n-1]*b1 + g*x[n-2]*b2
            // y[n] = y[n] - r[n-1]*a1 - r[n-2]*a2
            // Q29 <- Q15*Q14
            // Left Section 0 / sample 0
            // Right Section 0 / sample 0
            AE_MULZAAAA2Q16(accL0, accR0, ALsx, ARsx, AL_0b0b1b2, AR_0b0b1b2);
            AE_MULAAAA2Q16(accL0, accR0, ALsr, ARsr, AL_0_0a2a1, AR_0_0a2a1);
            AE_PKSR16(ALsr, accL0, 2);
            AE_PKSR16(ARsr, accR0, 2);
            // Left Section 0 / sample 1
            // Right Section 0 / sample 1
            AE_MULZAAAA2Q16(accL1, accR1, ALsx, ARsx, ALb0b1b2_0, ARb0b1b2_0);
            AE_MULAAAA2Q16(accL1, accR1, ALsr, ARsr, AL_0_0a2a1, AR_0_0a2a1);
            AE_PKSR16(ALsr, accL1, 2);
            AE_PKSR16(ARsr, accR1, 2);

            // Apply the total gain shift and format outputs y[n], y[n+1]
            // Q(15+gain) <- Q(29 - 14 + gain) w/ rounding and saturation
            ty0 = AE_TRUNCA32X2F64S(accL0, accL1, 32+2+gainL);
            ty1 = AE_TRUNCA32X2F64S(accR0, accR1, 32+2+gainR);
            y0 = AE_ROUND16X4F32SASYM(ty0, ty1);
            // Save output
            y0 = AE_SEL16_7520(y0, y0);
            AE_SA16X4_IP(y0, al_R, R);
        }
        AE_SA64POS_FP(al_R, R);

        //
        // Save section's state elements.
        //
        ALsxsr = AE_SEL16_7610(ALsx, ALsr);
        ARsxsr = AE_SEL16_7610(ARsx, ARsr);
        AE_S16X4X2_IP(ALsxsr, ARsxsr, St, 8*sz_i16);
    }
} /* stereo_bqriir16x16_df1_nd() */
