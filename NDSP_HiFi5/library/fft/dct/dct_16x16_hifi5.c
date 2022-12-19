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
    DCT 16x16 with scaling option 3
    C code optimized for HiFi4
    Integrit, 2006-2019
*/

/* Portable data types. */
#include "NatureDSP_types.h"
/* Signal Processing Library API. */
#include "NatureDSP_Signal_fft.h"
/* Common utility and macros declarations. */
#include "common.h"
#include "dct2_twd.h"

static const int16_t ALIGN(32) dsel_dct_split_tbl[] = {
    (3<<8)|7, (6<<8)|2, (1<<8)|5, (4<<8)|0
};

#ifndef AE_MULFC16RAS
inline_ ae_f16x4 func_AE_MULFC16RAS(ae_f16x4 *a0, ae_f16x4 *a1)
{
    ae_f16x4 vA, vB;
    ae_f32x2 hi, lo, zhi, zlo;
    ae_f16x4 prod;

    vA = *a0; vB = *a1;
    hi=AE_CVT32X2F16_32(vA);
    lo=AE_CVT32X2F16_10(vA);
    zhi=AE_MULFC32X16RAS_H(hi, vB);
    zlo=AE_MULFC32X16RAS_L(lo, vB);

    prod = AE_ROUND16X4F32SASYM(zhi, zlo);
    return prod;
}
#define AE_MULFC16RAS(a,b) func_AE_MULFC16RAS(&(a),&(b))
#endif

/*
   scaled fft with reordering
   NOTE: y is input and output
*/
void fft16_16x16(int16_t *y, const int16_t *ptwd);/* N=16 */
void fft32_16x16(int16_t *y, const int16_t *ptwd);/* N=32 */
/* pointer to complex fft function with reordering */
typedef void(*cfftProc_func)(int16_t *y, const int16_t *ptwd);

/*
	in-place split part of DCT:
	y[2*N]	input
	z[2*N]	output
	N       size of split algorithm
*/
static void dct_split_16x16(int16_t * z, const int16_t * y, const tdct2_twd *handle);

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
int dct_16x16(int16_t* y, int16_t* x, dct_handle_t h, int scalingOpt)
{
    const cfftProc_func cfftTbl[] =
    {
        fft16_16x16,
        fft32_16x16
    };
    int N, cfftIx;
    const tdct2_twd *descr=(const tdct2_twd *)h;
    NASSERT_ALIGN16(x);
    NASSERT_ALIGN16(y);
    NASSERT(scalingOpt==3);
    NASSERT(descr->magic==MAGIC_DCT2_16);
    N = descr->N;
    NASSERT(N>=32 && 0==(N&(N-1)));

    /* fft of half-size with reordering */
    cfftIx = 25-NSA(N);
    cfftIx &= 1;
    cfftTbl[cfftIx](x, (const int16_t *)descr->fft_twd);

    dct_split_16x16(y,x,descr);

    return 30-NSA(N);
} /* dct_16x16() */

/*
	in-place split part of DCT:
	y[2*N]	input
	z[2*N]	output
	N		size of FFT
*/
void dct_split_16x16(int16_t * z, const int16_t * y, const tdct2_twd *pdescr)
{
    int N;
    ae_int16   * restrict pzk;
    ae_int32   * restrict pzN_k_1;
    ae_int16   * restrict pzN2k;
    ae_int32   * restrict pzN2_k_1;
    const ae_int16x4 * restrict py0;
    const ae_int16x4 * restrict py1;
    const ae_int16x4 * restrict ptw0;
    const ae_int16x8 * restrict ptw1;
    int k;  
    ae_int16x4 y0, y1, z0, z1;
    ae_int16x4 t0, t1, tw0, tw1, tw2, tmptw;
    ae_int16x4 dsel_tbl;
    ae_valign aly0;

    NASSERT_ALIGN16(z);
    NASSERT_ALIGN16(y);
    N = pdescr->N;
    NASSERT( (N>=32) && 0==(N&(N-1)) );

    pzk      = (ae_int16 *)&z[1];
    pzN_k_1  = (ae_int32 *)&z[N - 2];
    pzN2k    = (ae_int16 *)&z[N/2 + 1]; 
    pzN2_k_1 = (ae_int32 *)&z[N/2 - 2];
    py0  = (const ae_int16x4 *)(&y[2]);
    py1  = (const ae_int16x4 *)(&y[N - 4]);
    ptw0 = (const ae_int16x4 *)(pdescr->rfft_split_twd);
    ptw1 = (const ae_int16x8 *)(pdescr->dct_twd);

    NASSERT_ALIGN8(py1); 
    dsel_tbl = AE_L16X4_I((const ae_int16x4 *)dsel_dct_split_tbl, 0);
    WUR_AE_SAR(0x81);

    tmptw = AE_L16X4_X((const ae_int16x4 *)ptw1, N*sizeof(int16_t)-2*sizeof(ae_int16x4));
    tw0 = AE_SEL16_5410(tmptw, tmptw);
    /* z[0]   = S_add_ss(y[0], y[1]);                     *
     * z[N/2] = S_mpy_ss(S_sub_ss(y[0], y[1]), pcosi[0]); */
    y0 = AE_L16_I((ae_int16*)y, 0);
    y1 = AE_L16_I((ae_int16*)y, 2);
    z0 = AE_ADD16S(y0, y1);
    z1 = AE_SUB16S(y0, y1);
    z1 = AE_MULFP16X4RAS(z1, tw0);
    AE_S16_0_I(z0, (ae_int16 *)z, 0);
    AE_S16_0_X(z1, (ae_int16 *)z, N/2*sizeof(int16_t));

    aly0 = AE_LA64_PP(py0);
    __Pragma("loop_count min=1");
    __Pragma("no_unroll");
    for (k=1; k < N/4-1; k+=2)
    {
        AE_LA16X4_IP(y0, aly0, py0);
        AE_L16X4_XP(y1, py1, -(int)sizeof(ae_int16x4));
        y1 = AE_SEL16_5432(y1, y1); 

        AE_ADDANDSUBRNG16RAS_S2(y0, y1);
        AE_DSEL16X4(t0, t1, y0, y1, dsel_tbl);

        // COSI = cosi(k);
        // Y0 = (imag(T0)*imag(COSI) - real(T0)*real(COSI)) + ...
        //    i*(real(T0)*imag(COSI) + imag(T0)*real(COSI));
        AE_L16X4_IP(tw0, ptw0, sizeof(ae_int16x4));
        t0 = AE_MULFC16RAS(t0, tw0);

        AE_L16X4X2_IP(tw1, tw2, ptw1, 2*sizeof(ae_int16x4));
        AE_ADDANDSUBRNG16RAS_S1(t1, t0);
        z0 = AE_MULFC16RAS(t0, tw1);
        z1 = AE_MULFC16RAS(t1, tw2);

        // z[k  ] = AE_MOVAD16_3(z0);
        // z[k+1] = AE_MOVAD16_1(z0);
        t0 = AE_SHORTSWAP(z0);
        t1 = AE_SEL16_4321(z0, z0);
        AE_S16_0_IP(t0, pzk, sizeof(int16_t));
        AE_S16_0_IP(t1, pzk, sizeof(int16_t));

        // z[N-k  ] = AE_MOVAD16_2(z0);
        // z[N-k-1] = AE_MOVAD16_0(z0);
        t0 = AE_SEL16_7520(z0, z0);
        t0 = AE_NEG16S(t0);
        AE_S32_L_IP(AE_MOVINT32X2_FROMINT16X4(t0), pzN_k_1, -(int)sizeof(int32_t));

        // z[N/2-k]   = AE_MOVAD16_3(z1);
        // z[N/2-k-1] = AE_MOVAD16_1(z1);
        t1 = AE_SEL16_7531(z1, z1);
        AE_S32_L_IP(AE_MOVINT32X2_FROMINT16X4(t1), pzN2_k_1, -(int)sizeof(int32_t));

        // z[N/2+k]   = AE_MOVAD16_2(z1);
        // z[N/2+k+1] = AE_MOVAD16_0(z1);
        t0 = AE_SEL16_5432(z1, z1);
        AE_S16_0_IP(t0, pzN2k, sizeof(int16_t));
        AE_S16_0_IP(z1, pzN2k, sizeof(int16_t));

    }
    {
        AE_LA16X4_IP(y0, aly0, py0);
        AE_L16X4_XP(y1, py1, -(int)sizeof(ae_int16x4));
        y1 = AE_SEL16_5432(y1, y1);

        AE_ADDANDSUBRNG16RAS_S2(y0, y1);
        AE_DSEL16X4(t0, t1, y0, y1, dsel_tbl);

        // COSI = cosi(k);
        // Y0 = (imag(T0)*imag(COSI) - real(T0)*real(COSI)) + ...
        //    i*(real(T0)*imag(COSI) + imag(T0)*real(COSI));
        tw0 = AE_L16X4_I(ptw0, 0);
        t0 = AE_MULFC16RAS(t0, tw0);

        AE_L16X4X2_IP(tw1, tw2, ptw1, 2*sizeof(ae_int16x4));
        AE_ADDANDSUBRNG16RAS_S1(t1, t0);
        z0 = AE_MULFC16RAS(t0, tw1);
        z1 = AE_MULFC16RAS(t1, tw2);

        t0 = AE_NEG16S(z0);
        z[    k] = AE_MOVAD16_3(z0);
        z[N  -k] = AE_MOVAD16_2(t0);
        z[N/2-k] = AE_MOVAD16_3(z1);
        z[N/2+k] = AE_MOVAD16_2(z1);
    }
    {
        /*
        int16_t t0_re = y[2 * N / 4 + 0];
        int16_t t0_im = y[2 * N / 4 + 1];
        int16_t tw1r = *pw1++;
        int16_t tw1i = *pw1++;

        z[N / 4] = cmulre(t0_re, t0_im, tw1r, tw1i);
        z[3 * N / 4] = cmulim(t0_re, t0_im, tw1r, tw1i);
        */
        t0 = AE_L16X4_X((ae_int16x4*)y, 2*(N/4)*sizeof(int16_t));
        t0 = AE_SEL16_5432(t0, t0);
        t0 = AE_MULFC16RAS(t0, tw2);
        z[1*N/4] = AE_MOVAD16_1(t0);
        z[3*N/4] = AE_MOVAD16_0(t0);
    }
}
