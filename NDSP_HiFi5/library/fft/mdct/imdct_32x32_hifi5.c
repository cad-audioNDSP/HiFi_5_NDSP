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
    IMDCT 32x32 with scaling option 3
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
  Modified Discrete Cosine Transform.
  These functions apply Modified DCT to input (convert 2N real data to N 
  spectral components) and make inverse conversion forming 2N numbers from 
  N inputs. Normally, combination of MDCT and DCT is invertible if applied 
  to subsequent data blocks with overlapping.
  Scaling:
      +-----------------------+--------------------------------------+
      |      Function         |           Scaling options            |
      +-----------------------+--------------------------------------+
      |       mdct_32x16      |  3 - fixed scaling before each stage |
      |       mdct_32x32      |  3 - fixed scaling before each stage |
      |      imdct_32x16      |  3 - fixed scaling before each stage |
      |      imdct_32x32      |  3 - fixed scaling before each stage |
      +-----------------------+--------------------------------------+
  NOTES:
     1. MDCT/IMDCT runs in-place algorithm so INPUT DATA WILL APPEAR DAMAGED 
     after the call.
     2. N - MDCT size (depends on selected MDCT handle)

  Precision: 
  32x16  32-bit input/outputs, 16-bit twiddles
  32x32  32-bit input/outputs, 32-bit twiddles

  -------------------------------------------------------------------------
  For MDCT:
  Input:
  x[2*N]      input signal
  h           MDCT handle
  scalingOpt  scaling option (see table above)
  Output:
  y[N]        output of transform 
  -------------------------------------------------------------------------
  For IMDCT:
  Input:
  x[N]        input signal
  h           IMDCT handle
  scalingOpt  scaling option (see table above)
  Output:
  y[2*N]      output of transform
  -------------------------------------------------------------------------
  Returned value:
              total number of right shifts occurred during scaling 
              procedure 
  Restriction:
  x,y         should not overlap
  x,y         aligned on 16-bytes boundary
-------------------------------------------------------------------------*/
int imdct_32x32( int32_t * y, int32_t * x, dct_handle_t h, int scalingOpt)
{  
    const tdct4_twd_fr16 *ptwd=(const tdct4_twd_fr16 *)h;
    const ae_int32x4 * restrict pyrd_pos0;
    const ae_int32x4 * restrict pyrd_pos1;
    const ae_int64x2 * restrict pyrd_neg0;
    const ae_int64x2 * restrict pyrd_neg1;
          ae_int32x4 * restrict pywr_pos0;
          ae_int32x4 * restrict pywr_pos1;
          ae_int64x2 * restrict pywr_neg0;
          ae_int64x2 * restrict pywr_neg1;
    ae_int32x2 t0, t1, t2, t3, t4, t5, t6, t7, r2, r3, r6, r7;
    ae_int64   T0, T1, T2, T3, T4, T5, T6, T7;
    int N,n,scl;
    NASSERT_ALIGN16(x);
    NASSERT_ALIGN16(y);
    NASSERT(scalingOpt == 3);
    N=ptwd->N;

    scl=dct4_32x32(y, x, h, scalingOpt);

    pyrd_pos0 = (const ae_int32x4 *)(y);
    pyrd_pos1 = (const ae_int32x4 *)(y+N/2);
    pyrd_neg0 = (const ae_int64x2 *)(y+N/2-4);
    pyrd_neg1 = (const ae_int64x2 *)(y+N-4);
    pywr_pos0 = (      ae_int32x4 *)(y+N);
    pywr_pos1 = (      ae_int32x4 *)(y);
    pywr_neg0 = (      ae_int64x2 *)(y+N/2-4);
    pywr_neg1 = (      ae_int64x2 *)(y+N-4);

    /* Transform 0...N*2 samples */
    __Pragma("loop_count min=1");
    for(n=0;n<(N>>4);n++)
    {
        AE_L32X2X2_IP(t0, t4, pyrd_pos0, 4*sizeof(int32_t));
        AE_L32X2X2_IP(t2, t6, pyrd_pos1, 4*sizeof(int32_t));
        AE_L64X2_IP  (T5, T1, pyrd_neg0, -4*(int)sizeof(int32_t));
        AE_L64X2_IP  (T7, T3, pyrd_neg1, -4*(int)sizeof(int32_t));
        t1 = AE_MOVINT32X2_FROMINT64(T1);
        t3 = AE_MOVINT32X2_FROMINT64(T3);
        t5 = AE_MOVINT32X2_FROMINT64(T5);
        t7 = AE_MOVINT32X2_FROMINT64(T7);

        AE_MULF2P32X16X4RAS(t0, t1, t0, t1, AE_MOVDA16(0xc000));/* -0.5 */
        AE_MULF2P32X16X4RAS(r2, r3, t2, t3, AE_MOVDA16(0xc000));/* -0.5 */
        AE_MULF2P32X16X4RAS(t2, t3, t2, t3, AE_MOVDA16(0x4000));/*  0.5 */
        AE_MULF2P32X16X4RAS(t4, t5, t4, t5, AE_MOVDA16(0xc000));/* -0.5 */
        AE_MULF2P32X16X4RAS(r6, r7, t6, t7, AE_MOVDA16(0xc000));/* -0.5 */
        AE_MULF2P32X16X4RAS(t6, t7, t6, t7, AE_MOVDA16(0x4000));/*  0.5 */

        T0 = AE_MOVINT64_FROMINT32X2(t0);
        T1 = AE_MOVINT64_FROMINT32X2(t1);
        T2 = AE_MOVINT64_FROMINT32X2(r2);
        T3 = AE_MOVINT64_FROMINT32X2(t3);
        T4 = AE_MOVINT64_FROMINT32X2(t4);
        T5 = AE_MOVINT64_FROMINT32X2(t5);
        T6 = AE_MOVINT64_FROMINT32X2(r6);
        T7 = AE_MOVINT64_FROMINT32X2(t7);
        __Pragma("no_reorder");
        AE_S32X2X2_X (t0, t4, pywr_pos0, N/2*sizeof(int32_t));
        AE_S32X2X2_IP(t1, t5, pywr_pos0, 4*sizeof(int32_t));
        AE_S32X2X2_X (r3, r7, pywr_pos1, N/2*sizeof(int32_t));
        AE_S32X2X2_IP(t2, t6, pywr_pos1, 4*sizeof(int32_t));
        AE_S64X2_X   (T4, T0, pywr_neg0, N*sizeof(int32_t));
        AE_S64X2_IP  (T7, T3, pywr_neg0, -4*(int)sizeof(int32_t));
        AE_S64X2_X   (T5, T1, pywr_neg1, N*sizeof(int32_t));
        AE_S64X2_IP  (T6, T2, pywr_neg1, -4*(int)sizeof(int32_t));
    }

    return scl+1;

}/* imdct_32x32() */
