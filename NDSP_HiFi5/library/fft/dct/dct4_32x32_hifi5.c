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
    DCT-IV 32x32 with scaling option 3
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
#include "dct4_twd.h"

/*-------------------------------------------------------------------------
  Discrete Cosine Transform.
  These functions apply DCT (Type II, Type IV) to input.
  Scaling:
      +-----------------------+--------------------------------------+
      |      Function         |           Scaling options            |
      +-----------------------+--------------------------------------+
      |       dct_16x16       |  3 - fixed scaling before each stage |
      |       dct_32x16       |  3 - fixed scaling before each stage |
      |       dct_32x32       |  3 - fixed scaling before each stage |
      |       dct4_32x16      |  3 - fixed scaling before each stage |
      |       dct4_32x32      |  3 - fixed scaling before each stage |
      +-----------------------+--------------------------------------+
  NOTES:
     1. DCT runs in-place algorithm so INPUT DATA WILL APPEAR DAMAGED after 
     the call.
     2. N - DCT size (depends on selected DCT handle)

  Precision: 
  16x16  16-bit input/outputs, 16-bit twiddles
  32x16  32-bit input/outputs, 16-bit twiddles
  32x32  32-bit input/outputs, 32-bit twiddles
  f      floating point

  Input:
  x[N]        input signal
  h           DCT handle
  scalingOpt  scaling option (see table above) 
              not applicable to the floating point function
  Output:
  y[N]        transform output
  
  Returned value:
              total number of right shifts occurred during scaling 
              procedure 
  Restriction:
  x,y         should not overlap
  x,y         aligned on 16-bytes boundary
-------------------------------------------------------------------------*/

#define LOG2_SZ_CI32P 4/* log2(2*sizeof(complex_fract32)) */

#ifndef AE_DSEL32X2_HH_LL
/*
   Equal to:
   a = AE_SEL32_HH(c, d)
   b = AE_SEL32_LL(c, d)
*/
#define AE_DSEL32X2_HH_LL(a,b,c,d) \
{\
    ae_int16x4 aa,bb,cc,dd,sel; \
    sel = AE_MOVINT16X4_FROMINT32X2(AE_MOVDA32X2((7L<<24)|(5L<<16)|(6<<8)|(4), (3L<<24)|(1L<<16)|(2<<8)|(0) )); \
    cc = AE_MOVINT16X4_FROMINT32X2(c); \
    dd = AE_MOVINT16X4_FROMINT32X2(d); \
    AE_DSEL16X4(aa,bb,cc,dd,sel); \
    a = AE_MOVINT32X2_FROMINT16X4(aa); \
    b = AE_MOVINT32X2_FROMINT16X4(bb); \
}
#endif

#ifndef AE_DSEL32X2_HL_LH
/*
   Equal to:
   a = AE_SEL32_HL(c, d)
   b = AE_SEL32_LH(c, d)
*/
#define AE_DSEL32X2_HL_LH(a,b,c,d) \
{\
    ae_int16x4 aa,bb,cc,dd,sel; \
    sel = AE_MOVINT16X4_FROMINT32X2(AE_MOVDA32X2((7L<<24)|(5L<<16)|(6<<8)|(4), (1L<<24)|(3L<<16)|(0<<8)|(2) )); \
    cc = AE_MOVINT16X4_FROMINT32X2(c); \
    dd = AE_MOVINT16X4_FROMINT32X2(d); \
    AE_DSEL16X4(aa,bb,cc,dd,sel); \
    a = AE_MOVINT32X2_FROMINT16X4(aa); \
    b = AE_MOVINT32X2_FROMINT16X4(bb); \
}
#endif

static void dct3p(uint64_t *y,uint64_t *x, int N, const tdct4_twd_fr32 *pdct4_twd);
static void dct3p_N16(uint64_t *y, uint64_t *x, const tdct4_twd_fr32 *pdct4_twd);
static void fft_cplx32x32_pair (uint64_t *y, uint64_t *x, const complex_fract32 *twd, int N);

int dct4_32x32( int32_t * y, int32_t * x, dct_handle_t h, int scalingOpt)
{  
          ae_int32x4 * restrict px0;
          ae_int32x4 * restrict px1;
          ae_int32x2 * restrict py0;
          ae_int32x2 * restrict py1;
          ae_int32x2 * restrict py2;
          ae_int32x2 * restrict py3;
    const ae_int32x4 * restrict ptwd0;
    const ae_int32x4 * restrict ptwd1;
    const tdct4_twd_fr32 *ptwd=(const tdct4_twd_fr32 *)h;
    ae_int32x2 Y0, Y1, Y2, Y3;
    ae_int32x2 WV0, WV1, WV2, WV3;
    ae_int32x2 tw0, tw1, tw2, tw3;
    ae_f64 ACC0, ACC1, ACC2, ACC3, ACC4, ACC5, ACC6, ACC7;
    int N,k;

    NASSERT_ALIGN16(x);
    NASSERT_ALIGN16(y);
    NASSERT(scalingOpt == 3);

    N=ptwd->N;

    if (N==32)
    {
        dct3p_N16((uint64_t*)y,(uint64_t*)x,ptwd);
    }
    else
    {
        dct3p((uint64_t*)y,(uint64_t*)x,N>>1,ptwd);
    }

    /* split part */
    ptwd0 = (const ae_int32x4 *)(ptwd->split);
    ptwd1 = (const ae_int32x4 *)(ptwd->split+N/4);
    px0 = (ae_int32x4 *)(x);
    px1 = (ae_int32x4 *)(x+N-2*2);
    py0 = (ae_int32x2 *)(y      );
    py1 = (ae_int32x2 *)(y+N/2-2);
    py2 = (ae_int32x2 *)(y+N/2  );
    py3 = (ae_int32x2 *)(y+N  -2);
    __Pragma("loop_count min=1");
    for (k=0; k<(N>>3); k++)
    {
        AE_L32X2X2_IP(WV0, WV2, px0,  4*sizeof(int32_t));
        AE_L32X2X2_IP(WV3, WV1, px1, -4*(int)sizeof(int32_t));

        AE_L32X2X2_IP(tw0, tw1, ptwd0, 4*sizeof(int32_t));
        AE_L32X2X2_IP(tw2, tw3, ptwd1, 4*sizeof(int32_t));

        AE_MULZASF2D32S_HH_LL(ACC0, ACC2, WV0, WV1, tw0, tw2);
        AE_MULZSSF2D32S_HL_LH(ACC1, ACC3, WV1, WV0, tw2, tw0);
        AE_MULZAAF2D32S_HH_LL(ACC4, ACC6, WV2, WV3, tw1, tw3);
        AE_MULZSAF2D32S_HL_LH(ACC5, ACC7, WV3, WV2, tw3, tw1);
        Y0 = AE_ROUND32X2F64SASYM(ACC0, ACC4);
        Y1 = AE_ROUND32X2F64SASYM(ACC5, ACC1);
        Y2 = AE_ROUND32X2F64SASYM(ACC2, ACC6);
        Y3 = AE_ROUND32X2F64SASYM(ACC7, ACC3);

        AE_S32X2_IP(Y0, py0,       2*sizeof(int32_t));
        AE_S32X2_XP(Y1, py1, -2*(int)sizeof(int32_t));
        AE_S32X2_IP(Y2, py2,       2*sizeof(int32_t));
        AE_S32X2_XP(Y3, py3, -2*(int)sizeof(int32_t));
    }
    return 30-NSA(N);
} /* dct4_32x32() */


/*
    DCT-III on pairs of data
    Input/output:
    x[N]
    Temporary:
    y[N]
*/
static void dct3p(uint64_t *y, uint64_t *x, int N, const tdct4_twd_fr32* pdct4_twd)
{
          ae_int32x4 * restrict px0;
          ae_int32x4 * restrict px1;
          ae_int32x4 * restrict py0;
          ae_int32x4 * restrict py1;
    const ae_int32x2 * restrict ptwd;
    ae_int32x2 X0, X1, X2, X3, Y0, Y1, Y2, Y3;
    ae_int32x2 a0re, a0im, b0re, b0im;
    ae_int32x2 a1re, a1im, b1re, b1im;
    ae_int32x2 coef, T0;
    int k;
    ae_valign al_x0, al_x1;

    /*
      Two processing loops are combined into one:

      1. Rearrange input data to compute DCT-IV via DCT-III. Even samples are w, odd - v:
        N2 = N/2;
        w(1) = x(1);
        v(N2) = x(N);
        for k=1:(N2-1)
          w(k+1) = x(2*k+1) + x(2*k);
          v(k)   = x(2*k)   - x(2*k+1);
        end

      2. Transform samples to compute DCT-III via real IFFT:
        y(1)     = x(1) * exp(1i*pi*0/(2*N) );
        y(N/2+1) = (x(N/2+1))*sqrt(2)/2;
        for k=1:(N/2-1)
          y(k+1)   = 0.5*    ( complex(x(k+1), -x(N-k+1)) * exp(1i*pi*k/(2*N)) );
          y(N-k+1) = 0.5*conj( complex(x(k+1), -x(N-k+1)) * exp(1i*pi*k/(2*N)) );
        end
    */

    ptwd = (const ae_int32x2 *)(pdct4_twd->dct3+1);
    WUR_AE_SAR(2);
    px0 = (ae_int32x4 *)((int32_t *)x);
    px1 = (ae_int32x4 *)((int32_t *)x+N*2-1);
    py0 = (ae_int32x4 *)((int32_t *)y);
    {
        X0 = AE_L32_I((ae_int32 *)px0, 0);
        T0 = AE_L32_I((ae_int32 *)px1, 0);
        X0 = AE_SEL32_LL(X0, T0);
        X1 = AE_L32_X((ae_int32 *)px0, (N-1)*sizeof(int32_t));
        T0 = AE_L32_X((ae_int32 *)px0, (N)*sizeof(int32_t));
        X1 = AE_SEL32_LL(T0, X1);

        AE_MULF2P32X16X4RAS(X0, X1, X0, X1, AE_MOVDA16(0x4000));/* shift right by 1 bit */
        X1 = AE_ADDSUB32S_HL_LH(X1, X1);

        X1 = AE_MULFP32X2RAS(X1, AE_MOVDA32(1518500250L));
        AE_DSEL32X2_HH_LL(Y0, Y1, X0, X1);
        AE_S32X2X2_IP(Y0, Y1, py0, 2*sizeof(ae_int32x2));
    }

    px0 = (ae_int32x4 *)((int32_t *)x+1);
    px1 = (ae_int32x4 *)((int32_t *)x+N*2-2);
    al_x0 = AE_LA64_PP(px0);
    al_x1 = AE_LA64_PP(px1);
    __Pragma("loop_count min=2");
    __Pragma("no_unroll");
    for (k=1; k<(N>>1); k++)
    {
        AE_LA32X2_IP (X0, al_x0, castxcc(ae_int32x2,px0));
        AE_LA32X2_RIP(X1, al_x1, castxcc(ae_int32x2,px1));
        AE_DSEL32X2_HH_LL(Y0, Y1, X0, X1);

        AE_ADDANDSUBRNG32(X0, X1, Y0, Y1);
        X1 = AE_MUL32JS(X1);

        AE_L32X2_RIP(coef, ptwd, 2*sizeof(int32_t));
        Y0 = AE_MULFC32RAS(X0, coef);
        Y1 = AE_MULFC32RAS(X1, coef);

        AE_S32X2X2_IP(Y0, Y1, py0, 2*sizeof(ae_int32x2));
    }
    __Pragma("no_reorder");
    /* real IFFT (Nyquist sample is packed to the imaginary part of y[] */
    WUR_AE_SAR(1);
    px0 = (ae_int32x4 *)(y);
    px1 = (ae_int32x4 *)(y+N-2);
    py0 = (ae_int32x4 *)(y);
    py1 = (ae_int32x4 *)(y+N-2);
    ptwd = (const ae_int32x2 *)(pdct4_twd->rfft);
    {
        AE_L32X2X2_IP(a0re, a1re, px0, 2*sizeof(ae_int32x2));
        AE_MULF2P32X16X4RAS(a0re, a1re, a0re, a1re, AE_MOVDA16(0x4000));/* shift right by 1 bit */
        AE_DSEL32X2_HH_LL(X0, Y0, a0re, a0re);
        AE_DSEL32X2_HH_LL(X1, Y1, a1re, a1re);
        X0 = AE_SUBADD32(X0, Y0);
        X1 = AE_SUBADD32(X1, Y1);
        AE_S32X2X2_IP(X0, X1, py0, 2*sizeof(ae_int32x2));
    }
    __Pragma("loop_count min=1");
    for (k=1; k<(N>>2); k++)
    {
        AE_L32X2X2_IP(a0re, a0im, px0, 2*sizeof(ae_int32x2));
        AE_L32X2X2_IP(a1re, a1im, px1, -2*(int)sizeof(ae_int32x2));

        AE_ADDANDSUBRNG32(Y0, Y1, a0re, a1re);
        AE_ADDANDSUBRNG32(X0, X1, a0im, a1im);
        AE_DSEL32X2_HL_LH(b0re, b1re, Y1, Y0);
        AE_DSEL32X2_HL_LH(b0im, b1im, X1, X0);

        AE_L32X2_IP(coef, ptwd, 2*sizeof(int32_t));
        a0re = AE_MULFC32RAS(b1re, coef);
        a0im = AE_MULFC32RAS(b1im, coef);

        X0 = b0re;
        X1 = b0im;
        Y0 = a0re;
        Y1 = a0im;

        a0re = AE_ADDSUB32(X0, Y0);
        a0im = AE_ADDSUB32(X1, Y1);
        a1re = AE_SUBADD32(Y0, X0);
        a1im = AE_SUBADD32(Y1, X1);

        AE_S32X2X2_IP(a0re, a0im, py0, 2*sizeof(ae_int32x2));
        AE_S32X2X2_IP(a1re, a1im, py1, -2*(int)sizeof(ae_int32x2));
    }
    {
        X0 = AE_L32X2_I((ae_int32x2 *)px0, 0);
        X1 = AE_L32X2_I((ae_int32x2 *)px1, sizeof(ae_int32x2));
        T0 = AE_SEL32_HH(X0, X1);
        T0 = AE_NEG32S(T0);
        Y0 = AE_SEL32_HL(T0, X0);
        Y1 = AE_SEL32_LL(T0, X1);
        AE_S32X2_I(Y0, (ae_int32x2 *)py0, 0);
        AE_S32X2_I(Y1, (ae_int32x2 *)py1, sizeof(ae_int32x2));
    }
    /* IFFT and final permutation */

    /* Real and imaginary parts are swapped on the first and last stages to
     * inverse the FFT:
     * conj(x) == -j*swap(x) =>
     * ifft(x) == conj(fft(conj(x)) == swap(fft(swap(x)))
     */
    fft_cplx32x32_pair (x,y,pdct4_twd->fft, N/2);
    py0 = (ae_int32x4 *)(y);
    py1 = (ae_int32x4 *)(y+N-2);
    px0 = (ae_int32x4 *)(x);
    __Pragma("loop_count min=1");
    for (k=0; k<(N>>2); k++)
    {
        AE_L32X2X2_IP(X0, X1, castxcc(ae_int32x4,py0),       2*sizeof(ae_int32x2));
        AE_L32X2X2_IP(X3, X2, castxcc(ae_int32x4,py1), -2*(int)sizeof(ae_int32x2));

        AE_DSEL32X2_HH_LL(Y1, Y0, X0, X1);
        AE_DSEL32X2_HH_LL(Y2, Y3, X3, X2);

        AE_S32X2X2_IP(Y0, Y2, castxcc(ae_int32x4,px0), 2*sizeof(ae_int32x2));
        AE_S32X2X2_IP(Y1, Y3, castxcc(ae_int32x4,px0), 2*sizeof(ae_int32x2));
    }
}

/*
    paired fft: operates with pairs of data packed in 64-bit words
    input/output:
    x[N]    - data
    Input:
    twd[3/4*N]  - twiddles
    N           - FFT size
    Temporary:
    y[N]
*/
static void fft_cplx32x32_pair (uint64_t* y,uint64_t* x, const complex_fract32* twd,  int N)
{
          ae_int32x4 * restrict px0;
          ae_int32x4 * restrict px1;
          ae_int32x4 * restrict px2;
          ae_int32x4 * restrict px3;
          ae_int32x4 * restrict py0;
          ae_int32x4 * restrict py1;
          ae_int32x4 * restrict py2;
          ae_int32x4 * restrict py3;
    const ae_int32x2 * restrict ptwd;
    int logN, idx, bitrevstride, stride;
    int m, n;
    int twdstep;

    ae_int32x2 a00, a10, a20, a30, b00, b10, b20, b30;
    ae_int32x2 a01, a11, a21, a31, b01, b11, b21, b31;
    ae_int32x2 tw1, tw2, tw3;

    NASSERT( x );
    NASSERT( y );
    NASSERT( twd );
    NASSERT( x != y );
    NASSERT_ALIGN( x, 16 );
    NASSERT_ALIGN( y, 16 );
    NASSERT( N>=8 && 0 == (N&(N-1)) );

    twdstep = 1;
    logN = 30 - NSA( N );

    /*----------------------------------------------------------------------------*
     * Perform the first stage. We use DIF, all permutations are deferred until   *
     * the last stage.                                                            */
    stride = N/4;
    ptwd = (const ae_int32x2 *)(twd);
    px0 = (ae_int32x4 *)(x);
    px1 = px0 + stride;
    px2 = px1 + stride;
    px3 = px2 + stride;
    py0 = (ae_int32x4 *)(y);
    py1 = py0 + stride;
    py2 = py1 + stride;
    py3 = py2 + stride;

    for ( n=0; n<stride; n++ )
    {
        AE_L32X2_IP(tw1, ptwd, 2*sizeof(int32_t));
        AE_L32X2_IP(tw2, ptwd, 2*sizeof(int32_t));
        AE_L32X2_IP(tw3, ptwd, 2*sizeof(int32_t));

        /* Real and imaginary parts are swapped on the first and last stages to
         * inverse the FFT:
         * conj(x) == -j*swap(x) =>
         * ifft(x) == conj(fft(conj(x)) == swap(fft(swap(x)))
         * Just in case, divide data by the FFT size. */
        AE_L32X2X2_IP(a00, a01, px0, 2*sizeof(ae_int32x2));
        AE_L32X2X2_IP(a10, a11, px1, 2*sizeof(ae_int32x2));
        AE_L32X2X2_IP(a20, a21, px2, 2*sizeof(ae_int32x2));
        AE_L32X2X2_IP(a30, a31, px3, 2*sizeof(ae_int32x2));

        AE_ADDANDSUBRNG32(b00, b20, a00, a20);
        AE_ADDANDSUBRNG32(b10, b30, a10, a30);
        AE_ADDANDSUBRNG32(b01, b21, a01, a21);
        AE_ADDANDSUBRNG32(b11, b31, a11, a31);
        b30 = AE_MUL32JS(b30);
        b31 = AE_MUL32JS(b31);
        AE_ADDANDSUBRNG32(a00, a20, b00, b10);
        AE_ADDANDSUBRNG32(a01, a21, b01, b11);
        AE_ADDANDSUBRNG32(a30, a10, b20, b30);
        AE_ADDANDSUBRNG32(a31, a11, b21, b31);

        b00 = a00;
        b01 = a01;
        b10 = AE_MULFC32RAS(a10, tw1);
        b11 = AE_MULFC32RAS(a11, tw1);
        b20 = AE_MULFC32RAS(a20, tw2);
        b21 = AE_MULFC32RAS(a21, tw2);
        b30 = AE_MULFC32RAS(a30, tw3);
        b31 = AE_MULFC32RAS(a31, tw3);

        /* Two middle quartiles are swapped on all but the last stage to use the bit reversal
         * permutation instead of the digit reverse. */
        AE_S32X2X2_IP(b00, b01, py0, 2*sizeof(ae_int32x2));
        AE_S32X2X2_IP(b20, b21, py1, 2*sizeof(ae_int32x2));
        AE_S32X2X2_IP(b10, b11, py2, 2*sizeof(ae_int32x2));
        AE_S32X2X2_IP(b30, b31, py3, 2*sizeof(ae_int32x2));
    }

    /*----------------------------------------------
      Perform second through the next to last stages.*/
    for ( stride/=4; stride>1; stride/=4 )
    {
        twdstep *= 4;

        for ( m=0; m*(4*stride)<N; m++ )
        {
            ptwd = (const ae_int32x2 *)(twd);

            px0 = (ae_int32x4 *)(y+4*m*2*stride);
            px1 = px0 + stride;
            px2 = px1 + stride;
            px3 = px2 + stride;
            py0 = (ae_int32x4 *)(y+4*m*2*stride);
            py1 = py0 + stride;
            py2 = py1 + stride;
            py3 = py2 + stride;

            for ( n=0; n<stride; n++ )
            {
                AE_L32X2X2_IP(tw1, tw2, castxcc(ae_int32x4,ptwd), 4*sizeof(int32_t));
                AE_L32X2_XP(tw3, ptwd, (twdstep*3-2)*2*sizeof(int32_t));

                AE_L32X2X2_IP(a00, a01, px0, 2*sizeof(ae_int32x2));
                AE_L32X2X2_IP(a10, a11, px1, 2*sizeof(ae_int32x2));
                AE_L32X2X2_IP(a20, a21, px2, 2*sizeof(ae_int32x2));
                AE_L32X2X2_IP(a30, a31, px3, 2*sizeof(ae_int32x2));

                AE_ADDANDSUBRNG32(b00, b20, a00, a20);
                AE_ADDANDSUBRNG32(b01, b21, a01, a21);
                AE_ADDANDSUBRNG32(b10, b30, a10, a30);
                AE_ADDANDSUBRNG32(b11, b31, a11, a31);
                b30 = AE_MUL32JS(b30);
                b31 = AE_MUL32JS(b31);
                AE_ADDANDSUBRNG32(a00, a20, b00, b10);
                AE_ADDANDSUBRNG32(a01, a21, b01, b11);
                AE_ADDANDSUBRNG32(a30, a10, b20, b30);
                AE_ADDANDSUBRNG32(a31, a11, b21, b31);

                b00 = a00;
                b01 = a01;
                b10 = AE_MULFC32RAS(a10, tw1);
                b11 = AE_MULFC32RAS(a11, tw1);
                b20 = AE_MULFC32RAS(a20, tw2);
                b21 = AE_MULFC32RAS(a21, tw2);
                b30 = AE_MULFC32RAS(a30, tw3);
                b31 = AE_MULFC32RAS(a31, tw3);

                /* Two middle quartiles are swapped on all but the last stage to use the bit reversal
                 * permutation instead of the digit reverse. */
                AE_S32X2X2_IP(b00, b01, py0, 2*sizeof(ae_int32x2));
                AE_S32X2X2_IP(b20, b21, py1, 2*sizeof(ae_int32x2));
                AE_S32X2X2_IP(b10, b11, py2, 2*sizeof(ae_int32x2));
                AE_S32X2X2_IP(b30, b31, py3, 2*sizeof(ae_int32x2));
            }
        }
    }

    /*----------------------------------------------------------------------------
    Last stage (radix-4 or radix-2 for odd powers of two) with bit reversal
    permutation. */
    idx = 0;
    bitrevstride = 0x80000000U >> (logN-3+LOG2_SZ_CI32P);

    if ( stride == 1 )
    {
        py0 = (ae_int32x4 *)(y);
        px0 = (ae_int32x4 *)(x+2*0*N/4);
        px1 = (ae_int32x4 *)(x+2*1*N/4);
        px2 = (ae_int32x4 *)(x+2*2*N/4);
        px3 = (ae_int32x4 *)(x+2*3*N/4);
        for ( n=0; n<N/4; n++ )
        {
            AE_L32X2X2_IP(a00, a01, py0, 2*sizeof(ae_int32x2));
            AE_L32X2X2_IP(a10, a11, py0, 2*sizeof(ae_int32x2));
            AE_L32X2X2_IP(a20, a21, py0, 2*sizeof(ae_int32x2));
            AE_L32X2X2_IP(a30, a31, py0, 2*sizeof(ae_int32x2));

            AE_ADDANDSUBRNG32(b00, b20, a00, a20);
            AE_ADDANDSUBRNG32(b01, b21, a01, a21);
            AE_ADDANDSUBRNG32(b10, b30, a10, a30);
            AE_ADDANDSUBRNG32(b11, b31, a11, a31);
            b30 = AE_MUL32JS(b30);
            b31 = AE_MUL32JS(b31);
            AE_ADDANDSUBRNG32(a00, a20, b00, b10);
            AE_ADDANDSUBRNG32(a01, a21, b01, b11);
            AE_ADDANDSUBRNG32(a30, a10, b20, b30);
            AE_ADDANDSUBRNG32(a31, a11, b21, b31);

            /* Real and imaginary parts are swapped on the first and last stages to
             * inverse the FFT:
             * conj(x) == -j*swap(x) =>
             * ifft(x) == conj(fft(conj(x)) == swap(fft(swap(x))) */
            AE_S32X2X2_X(a00, a01, px0, idx);
            AE_S32X2X2_X(a10, a11, px1, idx);
            AE_S32X2X2_X(a20, a21, px2, idx);
            AE_S32X2X2_X(a30, a31, px3, idx);

            idx = AE_ADDBRBA32(idx, bitrevstride);
        }
    }
    else
    {
        bitrevstride >>= 1;

        py0 = (ae_int32x4 *)(y);
        px0 = (ae_int32x4 *)(x);
        px2 = (ae_int32x4 *)(x+N);
        for ( n=0; n<N/2; n++ )
        {
            AE_L32X2X2_IP(a00, a01, py0, 2*sizeof(ae_int32x2));
            AE_L32X2X2_IP(a10, a11, py0, 2*sizeof(ae_int32x2));

            AE_ADDANDSUBRNG32(b00, b10, a00, a10);
            AE_ADDANDSUBRNG32(b01, b11, a01, a11);

            /* Real and imaginary parts are swapped on the first and last stages to
             * inverse the FFT:
             * conj(x) == -j*swap(x) =>
             * ifft(x) == conj(fft(conj(x)) == swap(fft(swap(x))) */
            AE_S32X2X2_X(b00, b01, px0, idx);
            AE_S32X2X2_X(b10, b11, px2, idx);

            idx = AE_ADDBRBA32(idx, bitrevstride);
        }
    }
} /* fft_cplx32x32_pair() */

/*
    DCT-III on pairs of data, N==16
    Input/output:
    x[16]
    Temporary:
    y[16]
*/
static void dct3p_N16(uint64_t *y, uint64_t *x, const tdct4_twd_fr32* pdct4_twd)
{
    const int N = 16;
          ae_int32x4 * restrict px0;
          ae_int32x4 * restrict px1;
          ae_int32x4 * restrict py0;
    const ae_int32x4 * restrict ptwd0;
    const ae_int32x2 * restrict ptwd1;

    /*
      Below three processing loops are combined into one:

      1. Rearrange input data to compute DCT-IV via DCT-III. Even samples are w, odd - v:
        N2 = N/2;
        w(1) = x(1);
        v(N2) = x(N);
        for k=1:(N2-1)
          w(k+1) = x(2*k+1) + x(2*k);
          v(k)   = x(2*k)   - x(2*k+1);
        end

      2. Transform samples to compute DCT-III via real IFFT:
        y(1)     = x(1) * exp(1i*pi*0/(2*N) );
        y(N/2+1) = (x(N/2+1))*sqrt(2)/2;
        for k=1:(N/2-1)
          y(k+1)   = 0.5*    ( complex(x(k+1), -x(N-k+1)) * exp(1i*pi*k/(2*N)) );
          y(N-k+1) = 0.5*conj( complex(x(k+1), -x(N-k+1)) * exp(1i*pi*k/(2*N)) );
        end

      3. Transform samples to compute real IFFT via complex FFT
        twd = exp(-2*pi*1j*(0:N/4-1)/N);
        a0 = x(1:N/4);
        a1 = wrev(x(N/4+2:N/2+1));
        b0 = a0+conj(a1);
        b1 = (a0-conj(a1))*1j.*conj(twd);
        a0 = 0.5*(b0+b1);
        a1 = 0.5*conj(b0-b1);
        y = [a0, conj(x(N/4+1)),wrev(a1(2:N/4))]; % N/2 complex samples

      NOTE: loop is fully unrolled
    */
    AE_MOVSARA7X2(2, 1);
    ptwd0 = (const ae_int32x4 *)(pdct4_twd->dct3+1);
    ptwd1 = (const ae_int32x2 *)(pdct4_twd->rfft);
    px0 = (ae_int32x4 *)((int32_t *)x);
    px1 = (ae_int32x4 *)((int32_t *)x+N*2-1);
    py0 = (ae_int32x4 *)((int32_t *)y);
    {
        ae_int32x2 X0, X1, X2, X3, X4, X5, X6, X7, X8, X9, XA, XB, XC, XD, XE, XF;
        ae_int32x2 Y0, Y1, Y2, Y3, Y4, Y5, Y6, Y7, Y8, Y9, YA, YB, YC, YD, YE, YF;
        ae_int32x2 T0;
        ae_int32x2 cs0, cs1;
        ae_int64 t64, u64;
        ae_valign al_x0, al_x1;

        /* 1st sample */
        X0 = AE_L32_I((ae_int32 *)px0, 0);
        T0 = AE_L32_I((ae_int32 *)px1, 0);
        X0 = AE_SEL32_LL(X0, T0);
        X1 = AE_L32_X((ae_int32 *)px0, (N-1)*sizeof(int32_t));
        T0 = AE_L32_X((ae_int32 *)px0, (N)*sizeof(int32_t));
        X1 = AE_SEL32_LL(T0, X1);

        AE_MULF2P32X16X4RAS(X0, X1, X0, X1, AE_MOVDA16(0x4000));/* shift right by 1 bit */
        X1 = AE_ADDSUB32S_HL_LH(X1, X1);

        X1 = AE_MULFP32X2RAS(X1, AE_MOVDA32(1518500250L));
        AE_DSEL32X2_HH_LL(Y0, Y1, X0, X1);

        AE_MULF2P32X16X4RAS(X0, X1, Y0, Y1, AE_MOVDA16(0x4000));/* shift right by 1 bit */
        AE_DSEL32X2_HH_LL(X0, Y0, X0, X0);
        AE_DSEL32X2_HH_LL(X1, Y1, X1, X1);
        Y0 = AE_SUBADD32(X0, Y0);
        Y1 = AE_SUBADD32(X1, Y1);

        AE_S32X2X2_IP(Y0, Y1, py0, 2*sizeof(ae_int32x2));

        px0 = (ae_int32x4 *)((int32_t *)x+1);
        px1 = (ae_int32x4 *)((int32_t *)x+N*2-2);
        al_x0 = AE_LA64_PP(px0);
        al_x1 = AE_LA64_PP(px1);
        AE_LA32X2_IP (X2, al_x0, castxcc(ae_int32x2,px0));
        AE_LA32X2_IP (X4, al_x0, castxcc(ae_int32x2,px0));
        AE_LA32X2_IP (X6, al_x0, castxcc(ae_int32x2,px0));
        AE_LA32X2_IP (X8, al_x0, castxcc(ae_int32x2,px0));
        AE_LA32X2_IP (XA, al_x0, castxcc(ae_int32x2,px0));
        AE_LA32X2_IP (XC, al_x0, castxcc(ae_int32x2,px0));
        AE_LA32X2_IP (XE, al_x0, castxcc(ae_int32x2,px0));
        AE_LA32X2_RIP(X3, al_x1, castxcc(ae_int32x2,px1));
        AE_LA32X2_RIP(X5, al_x1, castxcc(ae_int32x2,px1));
        AE_LA32X2_RIP(X7, al_x1, castxcc(ae_int32x2,px1));
        AE_LA32X2_RIP(X9, al_x1, castxcc(ae_int32x2,px1));
        AE_LA32X2_RIP(XB, al_x1, castxcc(ae_int32x2,px1));
        AE_LA32X2_RIP(XD, al_x1, castxcc(ae_int32x2,px1));
        AE_LA32X2_RIP(XF, al_x1, castxcc(ae_int32x2,px1));
        AE_DSEL32X2_HH_LL(Y2, Y3, X2, X3);
        AE_DSEL32X2_HH_LL(Y4, Y5, X4, X5);
        AE_DSEL32X2_HH_LL(Y6, Y7, X6, X7);
        AE_DSEL32X2_HH_LL(Y8, Y9, X8, X9);
        AE_DSEL32X2_HH_LL(YA, YB, XA, XB);
        AE_DSEL32X2_HH_LL(YC, YD, XC, XD);
        AE_DSEL32X2_HH_LL(YE, YF, XE, XF);

        AE_ADDANDSUBRNG32_H(X2, X3, Y2, Y3);
        X3 = AE_MUL32JS(X3);
        AE_ADDANDSUBRNG32_H(X4, X5, Y4, Y5);
        X5 = AE_MUL32JS(X5);
        AE_ADDANDSUBRNG32_H(X6, X7, Y6, Y7);
        X7 = AE_MUL32JS(X7);
        AE_ADDANDSUBRNG32_H(X8, X9, Y8, Y9);
        X9 = AE_MUL32JS(X9);
        AE_ADDANDSUBRNG32_H(XA, XB, YA, YB);
        XB = AE_MUL32JS(XB);
        AE_ADDANDSUBRNG32_H(XC, XD, YC, YD);
        XD = AE_MUL32JS(XD);
        AE_ADDANDSUBRNG32_H(XE, XF, YE, YF);
        XF = AE_MUL32JS(XF);

        AE_L64_IP(t64, castxcc(ae_int64,ptwd0), 2*sizeof(int32_t));
        cs0 = AE_MOVINT32X2_FROMINT64(t64);
        Y2 = AE_MULFC32RAS(X2, cs0);
        Y3 = AE_MULFC32RAS(X3, cs0);

        AE_L64X2_IP(t64, u64, castxcc(ae_int64x2,ptwd0), 4*sizeof(int32_t));
        cs0 = AE_MOVINT32X2_FROMINT64(t64); cs1 = AE_MOVINT32X2_FROMINT64(u64);
        Y4 = AE_MULFC32RAS(X4, cs0);
        Y5 = AE_MULFC32RAS(X5, cs0);
        Y6 = AE_MULFC32RAS(X6, cs1);
        Y7 = AE_MULFC32RAS(X7, cs1);
        
        AE_L64X2_IP(t64, u64, castxcc(ae_int64x2,ptwd0), 4*sizeof(int32_t));
        cs0 = AE_MOVINT32X2_FROMINT64(t64); cs1 = AE_MOVINT32X2_FROMINT64(u64);
        Y8 = AE_MULFC32RAS(X8, cs0);
        Y9 = AE_MULFC32RAS(X9, cs0);
        YA = AE_MULFC32RAS(XA, cs1);
        YB = AE_MULFC32RAS(XB, cs1);

        AE_L64X2_IP(t64, u64, castxcc(ae_int64x2,ptwd0), 4*sizeof(int32_t));
        cs0 = AE_MOVINT32X2_FROMINT64(t64); cs1 = AE_MOVINT32X2_FROMINT64(u64);
        YC = AE_MULFC32RAS(XC, cs0);
        YD = AE_MULFC32RAS(XD, cs0);
        YE = AE_MULFC32RAS(XE, cs1);
        YF = AE_MULFC32RAS(XF, cs1);

        /* real IFFT (Nyquist sample is packed to the imaginary part of y[]) */
        AE_ADDANDSUBRNG32_L(X2, XE, Y2, YE);
        AE_ADDANDSUBRNG32_L(X3, XF, Y3, YF);
        AE_ADDANDSUBRNG32_L(X4, XC, Y4, YC);
        AE_ADDANDSUBRNG32_L(X5, XD, Y5, YD);
        AE_ADDANDSUBRNG32_L(X6, XA, Y6, YA);
        AE_ADDANDSUBRNG32_L(X7, XB, Y7, YB);
        AE_DSEL32X2_HL_LH(YE, Y2, XE, X2);
        AE_DSEL32X2_HL_LH(YF, Y3, XF, X3);
        AE_DSEL32X2_HL_LH(YC, Y4, XC, X4);
        AE_DSEL32X2_HL_LH(YD, Y5, XD, X5);
        AE_DSEL32X2_HL_LH(YA, Y6, XA, X6);
        AE_DSEL32X2_HL_LH(YB, Y7, XB, X7);

        AE_L32X2_IP(cs0, ptwd1, 2*sizeof(int32_t));
        Y2 = AE_MULFC32RAS(Y2, cs0);
        Y3 = AE_MULFC32RAS(Y3, cs0);
        AE_L32X2_IP(cs0, ptwd1, 2*sizeof(int32_t));
        Y4 = AE_MULFC32RAS(Y4, cs0);
        Y5 = AE_MULFC32RAS(Y5, cs0);
        AE_L32X2_IP(cs0, ptwd1, 2*sizeof(int32_t));
        Y6 = AE_MULFC32RAS(Y6, cs0);
        Y7 = AE_MULFC32RAS(Y7, cs0);

        X2 = AE_ADDSUB32(YE, Y2);
        X3 = AE_ADDSUB32(YF, Y3);
        X4 = AE_ADDSUB32(YC, Y4);
        X5 = AE_ADDSUB32(YD, Y5);
        X6 = AE_ADDSUB32(YA, Y6);
        X7 = AE_ADDSUB32(YB, Y7);
        XA = AE_SUBADD32(Y6, YA);
        XB = AE_SUBADD32(Y7, YB);
        XC = AE_SUBADD32(Y4, YC);
        XD = AE_SUBADD32(Y5, YD);
        XE = AE_SUBADD32(Y2, YE);
        XF = AE_SUBADD32(Y3, YF);

        AE_S32X2X2_IP(X2, X3, py0, 2*sizeof(ae_int32x2));
        AE_S32X2X2_IP(X4, X5, py0, 2*sizeof(ae_int32x2));
        AE_S32X2X2_IP(X6, X7, py0, 2*sizeof(ae_int32x2));

        T0 = AE_SEL32_HH(Y8, Y9);
        T0 = AE_NEG32(T0);
        X8 = AE_SEL32_HL(T0, Y8);
        X9 = AE_SEL32_LL(T0, Y9);
        AE_S32X2X2_IP(X8, X9, py0, 2*sizeof(ae_int32x2));

        AE_S32X2X2_IP(XA, XB, py0, 2*sizeof(ae_int32x2));
        AE_S32X2X2_IP(XC, XD, py0, 2*sizeof(ae_int32x2));
        AE_S32X2X2_IP(XE, XF, py0, 2*sizeof(ae_int32x2));
    }
    __Pragma("no_reorder");
    /* IFFT (on pairs of data, N=8) and final permutation */
    WUR_AE_SAR(1);
    {
        ae_int32x2 a00, a01, a10, a11, a20, a21, a30, a31;
        ae_int32x2 a40, a41, a50, a51, a60, a61, a70, a71;
        ae_int32x2 b00, b01, b10, b11, b20, b21, b30, b31;
        ae_int32x2 b40, b41, b50, b51, b60, b61, b70, b71;
        ae_int32x2 tw1, tw2, tw3, tw5, tw6, tw7;

        /*----------------------------------------------------------------------------*
         * Perform the first stage. We use DIF, all permutations are deferred until   *
         * the last stage.                                                            */
        ptwd0 = (const ae_int32x4 *)(pdct4_twd->fft);
        px0 = (ae_int32x4 *)(x);
        py0 = (ae_int32x4 *)(y);

        /* Real and imaginary parts are swapped on the first and last stages to
         * inverse the FFT:
         * conj(x) == -j*swap(x) =>
         * ifft(x) == conj(fft(conj(x)) == swap(fft(swap(x)))
         * Just in case, divide data by the FFT size. */

        AE_L32X2X2_IP(a01, a00, py0, 2*sizeof(ae_int32x2));
        AE_L32X2X2_IP(a41, a40, py0, 2*sizeof(ae_int32x2));
        AE_L32X2X2_IP(a11, a10, py0, 2*sizeof(ae_int32x2));
        AE_L32X2X2_IP(a51, a50, py0, 2*sizeof(ae_int32x2));
        AE_L32X2X2_IP(a21, a20, py0, 2*sizeof(ae_int32x2));
        AE_L32X2X2_IP(a61, a60, py0, 2*sizeof(ae_int32x2));
        AE_L32X2X2_IP(a31, a30, py0, 2*sizeof(ae_int32x2));
        AE_L32X2X2_IP(a71, a70, py0, 2*sizeof(ae_int32x2));

        AE_ADDANDSUBRNG32(b00, b20, a00, a20);
        AE_ADDANDSUBRNG32(b01, b21, a01, a21);
        AE_ADDANDSUBRNG32(b10, b30, a10, a30);
        AE_ADDANDSUBRNG32(b11, b31, a11, a31);
        b30 = AE_MUL32JS(b30);
        b31 = AE_MUL32JS(b31);
        AE_ADDANDSUBRNG32(a00, a20, b00, b10);
        AE_ADDANDSUBRNG32(a01, a21, b01, b11);
        AE_ADDANDSUBRNG32(a30, a10, b20, b30);
        AE_ADDANDSUBRNG32(a31, a11, b21, b31);

        AE_ADDANDSUBRNG32(b40, b60, a40, a60);
        AE_ADDANDSUBRNG32(b41, b61, a41, a61);
        AE_ADDANDSUBRNG32(b50, b70, a50, a70);
        AE_ADDANDSUBRNG32(b51, b71, a51, a71);
        b70 = AE_MUL32JS(b70);
        b71 = AE_MUL32JS(b71);
        AE_ADDANDSUBRNG32(a40, a60, b40, b50);
        AE_ADDANDSUBRNG32(a41, a61, b41, b51);
        AE_ADDANDSUBRNG32(a70, a50, b60, b70);
        AE_ADDANDSUBRNG32(a71, a51, b61, b71);

        AE_L32X2X2_IP(tw1, tw2, ptwd0, 2*sizeof(ae_int32x2));
        AE_L32X2X2_IP(tw3, tw5, ptwd0, 2*sizeof(ae_int32x2));
        AE_L32X2X2_IP(tw6, tw7, ptwd0, 2*sizeof(ae_int32x2));
        b00 = a00;
        b01 = a01;
        b40 = a40;
        b41 = a41;
        b10 = AE_MULFC32RAS(a10, tw1);
        b11 = AE_MULFC32RAS(a11, tw1);
        b20 = AE_MULFC32RAS(a20, tw2);
        b21 = AE_MULFC32RAS(a21, tw2);
        b30 = AE_MULFC32RAS(a30, tw3);
        b31 = AE_MULFC32RAS(a31, tw3);
        b50 = AE_MULFC32RAS(a50, tw5);
        b51 = AE_MULFC32RAS(a51, tw5);
        b60 = AE_MULFC32RAS(a60, tw6);
        b61 = AE_MULFC32RAS(a61, tw6);
        b70 = AE_MULFC32RAS(a70, tw7);
        b71 = AE_MULFC32RAS(a71, tw7);

        /*-------------------------------------------------------------------
          Last stage (radix-2) with bit reversal permutation.
          Final DCT-3 algorithm permutation is combined with FFT permutation.
        -------------------------------------------------------------------*/
        a00 = b00; a01 = b01;
        a10 = b40; a11 = b41;
        a20 = b20; a21 = b21;
        a30 = b60; a31 = b61;
        a40 = b10; a41 = b11;
        a50 = b50; a51 = b51;
        a60 = b30; a61 = b31;
        a70 = b70; a71 = b71;

        AE_ADDANDSUBRNG32(b00, b10, a00, a10);
        AE_ADDANDSUBRNG32(b01, b11, a01, a11);
        AE_ADDANDSUBRNG32(b20, b30, a20, a30);
        AE_ADDANDSUBRNG32(b21, b31, a21, a31);
        AE_ADDANDSUBRNG32(b40, b50, a40, a50);
        AE_ADDANDSUBRNG32(b41, b51, a41, a51);
        AE_ADDANDSUBRNG32(b60, b70, a60, a70);
        AE_ADDANDSUBRNG32(b61, b71, a61, a71);

        AE_DSEL32X2_HH_LL(a00, a01, b01, b00);
        AE_DSEL32X2_HH_LL(a10, a11, b11, b10);
        AE_DSEL32X2_HH_LL(a20, a21, b21, b20);
        AE_DSEL32X2_HH_LL(a30, a31, b31, b30);
        AE_DSEL32X2_HH_LL(a40, a41, b41, b40);
        AE_DSEL32X2_HH_LL(a50, a51, b51, b50);
        AE_DSEL32X2_HH_LL(a60, a61, b61, b60);
        AE_DSEL32X2_HH_LL(a70, a71, b71, b70);

        /* Real and imaginary parts are swapped on the first and last stages to
         * inverse the FFT:
         * conj(x) == -j*swap(x) =>
         * ifft(x) == conj(fft(conj(x)) == swap(fft(swap(x))) */
        AE_S32X2X2_I(a01, a70, px0, 0*2*sizeof(ae_int32x2));
        AE_S32X2X2_I(a00, a71, px0, 1*2*sizeof(ae_int32x2));
        AE_S32X2X2_I(a41, a30, px0, 2*2*sizeof(ae_int32x2));
        AE_S32X2X2_I(a40, a31, px0, 3*2*sizeof(ae_int32x2));
        AE_S32X2X2_I(a21, a50, px0, 4*2*sizeof(ae_int32x2));
        AE_S32X2X2_I(a20, a51, px0, 5*2*sizeof(ae_int32x2));
        AE_S32X2X2_I(a61, a10, px0, 6*2*sizeof(ae_int32x2));
        AE_S32X2X2_I(a60, a11, px0, 7*2*sizeof(ae_int32x2));
    } /* IFFT with permutation */
}
