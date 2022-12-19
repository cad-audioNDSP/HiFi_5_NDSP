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
    Complex-valued forward FFT: 16-bit data, 16-bit twiddle factors
    C code optimized for HiFi5
  IntegrIT, 2006-2019
*/

/* Library API. */
#include "NatureDSP_Signal_fft.h"
/* Common utility macros. */
#include "common.h"
/* Twiddle factor tables and FFT descriptor structure. */
#include "fft_x16_common.h"
#include "fft_16x16_stages.h"

#define SWAP_PTR(_x, _y) {int16_t *tmp = _x; _x = _y ; _y = tmp; } 

/*-------------------------------------------------------------------------
  FFT on Complex Data
  These functions make FFT on complex data.
    Scaling: 
      +-------------------+----------------------------------------+
      |      Function     |           Scaling options              |
      +-------------------+----------------------------------------+
      |  fft_cplx16x16    |  2 - 16-bit dynamic scaling            | 
      |                   |  3 - fixed scaling before each stage   | 
      |  fft_cplx32x32    |  2 - 32-bit dynamic scaling            | 
      |                   |  3 - fixed scaling before each stage   | 
      |  fft_cplx32x16    |  2 - 32-bit dynamic scaling            | 
      |                   |  3 - fixed scaling before each stage   | 
      +-------------------+----------------------------------------+
  NOTES:
  1. Bit-reversing permutation is done here. 
  2. FFT runs in-place algorithm so INPUT DATA WILL APPEAR DAMAGED after 
     the call
  3. 32x32,32x16,16x16 FFTs support mixed radix transforms 
  4. N - FFT size

  Precision: 
  32x16  32-bit input/outputs, 16-bit twiddles
  32x32  32-bit input/outputs, 32-bit twiddles
  16x16  16-bit input/outputs, 16-bit twiddles
 
  Input:
  x[2*N]     complex input signal. Real and imaginary data are interleaved 
             and real data goes first
  scalingOpt scaling option (see table above)
  Output:
  y[2*N]     output spectrum. Real and imaginary data are interleaved and 
             real data goes first

  Returned value: total number of right shifts occurred during scaling 
                  procedure

  Restrictions:
  x,y        should not overlap
  x,y        aligned on a 16-bytes boundary

-------------------------------------------------------------------------*/
int fft_cplx16x16( int16_t* y, int16_t* x, fft_handle_t h, int scalingOption )
{
  const ae_int16x4 * restrict px;
  const fft_cplx_x16_descr_t* pdescr = (const fft_cplx_x16_descr_t*)h;
  int bexp = -1;
  int v = 1;
  int s = 0;
  int shift = 0;

  int16_t *pdest = y; 
  const int  N = pdescr->N; 
  const int *tw_step = pdescr->tw_step;
  const cint16ptr_t_fft *tw_tab = pdescr->twd; 
  const fft_cplx16x16_stage_t*fnstg_Tbl = (const fft_cplx16x16_stage_t*)((scalingOption == 2) ? pdescr->fnstages_16x16_s2 : pdescr->fnstages_16x16_s3);

  NASSERT(fnstg_Tbl!=NULL);
  NASSERT_ALIGN16(x); 
  NASSERT_ALIGN16(y);
  NASSERT(x!=y);
  NASSERT(scalingOption == 2 || scalingOption == 3); 
  NASSERT(N%8==0);

  if  (scalingOption==2)
  {
    int i;
#if 0
    ae_int16x4 acc = AE_MOVINT16X4_FROMINT32X2(AE_MOVI(0)), tmp;
    px = (ae_int16x4*)x;
    __Pragma("loop_count min=4 factor=4");
    for (i = 0; i < (N >> 1); i++)   
    {
      AE_L16X4_IP(tmp, px, sizeof(*px));
      tmp = AE_ABS16S(tmp); 
      acc = AE_OR16(acc, tmp);
    }
    acc = AE_OR16(acc, AE_SEL16_5432(acc, acc));
    acc = AE_OR16(acc, AE_SHORTSWAP(acc));

    i = AE_MOVAD16_0(acc);
    bexp = NSA(i) - 16;
    XT_MOVEQZ(bexp, 0, i);
#else
    ae_int16x4 acc,a0,a1;
    a0=a1= AE_ZERO16();
    px = (ae_int16x4*)x; 
    __Pragma("loop_count min=2 factor=2");
    for (i = 0; i < (N >> 2); i++)   
    {
        ae_int16x4 t0,t1;
        AE_L16X4X2_IP(t0,t1, castxcc(ae_int16x8,px), 2*sizeof(ae_int16x4));
        a0=AE_MAX16(a0,AE_ABS16S(t0));
        a1=AE_MAX16(a1,AE_ABS16S(t1));
    }
    acc = AE_MAX16(a0,a1);
    i = AE_RMAX16X4(acc);
    bexp = NSA(i) - 16;
    XT_MOVEQZ(bexp, 0, i);
#endif
  }

  while ( tw_tab[s] != NULL )
  {
#if 1//(fnstg_Tbl!=NULL)
    {
       NASSERT(fnstg_Tbl[0]!=NULL);
       shift += fnstg_Tbl[0](tw_tab[s], x, y, N, &v, tw_step[s], &bexp);
        fnstg_Tbl++;
    }
#else
    {
        fft_cplx16x16_stage_t stg_fn;
        int stg_type;
        stg_type = (int)stg_typeTbl[s];
        NASSERT(stg_type < NUM_FFT_STAGE_TYPES);
        stg_fn = stg_Tbl[stg_type];
        NASSERT(stg_fn != NULL);
        shift += stg_fn(tw_tab[s], x, y, N, &v, tw_step[s], &bexp);
    }
#endif
    SWAP_PTR(x, y);
    s++;
  }


  if (y != pdest)
  {  
    /* Execute the last stage inplace */
    y = x;
  }

  /* Last stage */
#if 1
        NASSERT(fnstg_Tbl[0]);
        shift += fnstg_Tbl[0](tw_tab[s], x, y, N, &v, tw_step[s], &bexp);
#else    
    {
        fft_cplx16x16_stage_t stg_fn;
        int stg_type;
        stg_type = (int)stg_typeTbl[s];
        NASSERT(stg_type < NUM_FFT_STAGE_TYPES);
        stg_fn = stg_Tbl[stg_type];
        NASSERT(stg_fn != NULL);
        shift += stg_fn(tw_tab[s], x, y, N, &v, tw_step[s], &bexp);
    }
#endif
    return shift;
} /* fft_cplx16x16() */
