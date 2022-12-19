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
    MDCT 32x32 with scaling option 3
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
int mdct_32x32( int32_t * y, int32_t * x, dct_handle_t h, int scalingOpt)
{  
    const tdct4_twd_fr16 *ptwd=(const tdct4_twd_fr16 *)h;
    const ae_int32x4 * restrict pxrd_pos0;
    const ae_int64x2 * restrict pxrd_neg0;
    const ae_int32x4 * restrict pxrd_pos1;
    const ae_int64x2 * restrict pxrd_neg1;
    ae_int32x4 * pxwr_pos;
    ae_int64x2 * pxwr_neg;
    ae_int32x2 t0, t1, t2, t3, t4, t5, t6, t7;
    ae_int32x2 t8, t9, tA, tB, tC, tD, tE, tF;
    ae_int64   t64, u64;
    int N,n,scl;
    NASSERT_ALIGN16(x);
    NASSERT_ALIGN16(y);
    NASSERT(scalingOpt == 3);
    N=ptwd->N;

    pxrd_pos0 = (const ae_int32x4 *)(x           );
    pxrd_neg0 = (const ae_int64x2 *)(x  +(N>>1)-4);
    pxrd_pos1 = (const ae_int32x4 *)(x+N         );
    pxrd_neg1 = (const ae_int64x2 *)(x+N+(N>>1)-4);
    pxwr_pos  = (      ae_int32x4 *)(x           );
    pxwr_neg  = (      ae_int64x2 *)(x  +(N>>1)-4);
    WUR_AE_SAR(1);

    /* Transform 2*N samples into N for use in DCT-IV */
    __Pragma("loop_count min=1");
    for(n=0;n<(N>>4);n++)
    {
        AE_L32X2X2_X (t7, tF, pxrd_pos0, (N>>1)*sizeof(int32_t));
        AE_L32X2X2_X (t1, t9, pxrd_pos1, (N>>1)*sizeof(int32_t));
        AE_L32X2X2_IP(t4, tC, pxrd_pos0, 4*sizeof(int32_t));
        AE_L32X2X2_IP(t2, tA, pxrd_pos1, 4*sizeof(int32_t));

        AE_L64X2_X (u64, t64, pxrd_neg1, (N>>1)*sizeof(int32_t));
        t3 = AE_MOVINT32X2_FROMINT64(t64); tB = AE_MOVINT32X2_FROMINT64(u64);
        AE_L64X2_X (u64, t64, pxrd_neg0, (N>>1)*sizeof(int32_t));
        t5 = AE_MOVINT32X2_FROMINT64(t64); tD = AE_MOVINT32X2_FROMINT64(u64);
        AE_L64X2_IP(u64, t64, pxrd_neg1, -4*(int)sizeof(int32_t));
        t0 = AE_MOVINT32X2_FROMINT64(t64); t8 = AE_MOVINT32X2_FROMINT64(u64);
        AE_L64X2_IP(u64, t64, pxrd_neg0, -4*(int)sizeof(int32_t));
        t6 = AE_MOVINT32X2_FROMINT64(t64); tE = AE_MOVINT32X2_FROMINT64(u64);

        t0 = AE_ADDRNG32(t0, t1);
        t2 = AE_ADDRNG32(t2, t3);
        t4 = AE_SUBRNG32(t4, t5);
        t6 = AE_SUBRNG32(t6, t7);
        t8 = AE_ADDRNG32(t8, t9);
        tA = AE_ADDRNG32(tA, tB);
        tC = AE_SUBRNG32(tC, tD);
        tE = AE_SUBRNG32(tE, tF);
        t0 = AE_NEG32S(t0);
        t2 = AE_NEG32S(t2);
        t8 = AE_NEG32S(t8);
        tA = AE_NEG32S(tA);

        /* make N/2...N-1 samples */
        AE_S32X2X2_X(t4, tC, pxwr_pos, (N>>1)*sizeof(int32_t));
        t64 = AE_MOVINT64_FROMINT32X2(t6);
        u64 = AE_MOVINT64_FROMINT32X2(tE);
        AE_S64X2_X(u64, t64, pxwr_neg, (N>>1)*sizeof(int32_t));
        /* make 0...N/2-1 samples */
        AE_S32X2X2_IP(t0, t8, pxwr_pos, 4*sizeof(int32_t));
        t64 = AE_MOVINT64_FROMINT32X2(t2);
        u64 = AE_MOVINT64_FROMINT32X2(tA);
        AE_S64X2_IP(u64, t64, pxwr_neg, -4*(int)sizeof(int32_t));
    }

    scl=dct4_32x32(y, x, h, scalingOpt);

    return scl+1;

}/* mdct_32x32() */
