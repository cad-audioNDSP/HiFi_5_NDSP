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

/* Library API */
#include "NatureDSP_Signal_vector.h"
#include "common.h"


/*===========================================================================
  Vector matematics:
  vec_dot              Vector Dot Product
===========================================================================*/

/*-------------------------------------------------------------------------
  Vector Dot product
  These routines take two vectors and calculates their dot product.
  Two versions of routines are available: regular versions (vec_dot64x32,
  vec_dot64x64, vec_dot64x64i, vec_dot32x16, vec_dot32x32,vec_dot16x16, 
  vec_dotf) work with arbitrary arguments, faster versions (vec_dot64x32_fast, 
  vec_dot64x64_fast, vec_dot64x64i_fast, vec_dot32x16_fast, vec_dot32x32_fast,
  vec_dot16x16_fast) apply some restrictions.  
  NOTE:
  vec_dot16x16_fast utilizes 32-bit saturating accumulator, so input data 
  should be scaled properly to avoid erroneous results.

  Precision: 
  64x32  64x32-bit data, 64-bit output (fractional multiply Q63xQ31->Q63)
  64x64  64x64-bit data, 64-bit output (fractional multiply Q63xQ63->Q63)
  64x64i 64x64-bit data, 64-bit output (low 64 bit of integer multiply)
  32x32  32x32-bit data, 64-bit output
  32x16  32x16-bit data, 64-bit output
  16x16  16x16-bit data, 64-bit output for regular version and 32-bit for 
                        fast version
  f      single precision floating point

  Input:
  x[N]  input data, Q15, Q31, Q63 or floating point
  y[N]  input data, Q15, Q31, Q63 or floating point
  N	    length of vectors
  Returns:
  dot product of all data pairs, Q31, Q63 or floating point

  Restrictions:
  Regular versions:
    none
  Faster versions:
    x,y - aligned on 16-byte boundary
    N   - multiple of 4
-------------------------------------------------------------------------*/
int64_t vec_dot64x32_fast (const int64_t * restrict x,const int32_t * restrict y,int N)
{
    int n;
    ae_int32x2 y0, x0,x1,x2,x3,y1,yy0;
    xtbool2 b=1;
    ae_int64 xw0,xw1;                                     
    ae_int64 W0,W1,W2, ACC;
    ae_ep e0,e1,e2;
    const ae_int64x2 *  restrict  px = (const ae_int64x2 *) x;
    const ae_int32x4 *  restrict  py = (const ae_int32x4 *) y;
    NASSERT_ALIGN16(x);
    NASSERT_ALIGN16(y);
    NASSERT((N&3)==0);
    if (N<=0) return 0;
    ACC=AE_ZERO64();
    W0=W1=W2=AE_ZERO64();
    e0=e1=e2=AE_MOVEA(0);
    __Pragma("loop_count min=1")
    for (n=0; n<((N>>2)); n++)
    {
        AE_L64X2_IP(xw0,xw1,px,sizeof(ae_int32x4));x0 = AE_MOVINT32X2_FROMINT64(xw0);x1 = AE_MOVINT32X2_FROMINT64(xw1);
        AE_L64X2_IP(xw0,xw1,px,sizeof(ae_int32x4));x2 = AE_MOVINT32X2_FROMINT64(xw0);x3 = AE_MOVINT32X2_FROMINT64(xw1);
        AE_L32X2X2_IP(y0,y1,py,sizeof(ae_int32x4));
        // tricky usage of HiFi5 multiplies - it uses the fact that AE_MULAAD32USEP_HL_LH is 
        // scheduled better than than AE_MULA32USEP_LL
        yy0=y0; AE_MOVT32X2(yy0,0,b);   
        AE_MULAAD32USEP_HL_LH(e0,W0,x0,yy0);
        AE_MULA32USEP_LL(e2,W2,x1,y0);
        AE_MULA32USEP_LH(e0,W0,x2,y1);
        AE_MULA32USEP_LL(e2,W0,x3,y1);

        x0=AE_SEL32_HH(x0,x1);
        x2=AE_SEL32_HH(x2,x3);
        AE_MULAAD32EP_HH_LL(e1,W1,x0,y0);
        AE_MULAAD32EP_HH_LL(e1,W1,x2,y1);
    }

    AE_ADD72(e0,W0,e2,W2);
    AE_ADD72(e1,W1,e1,W1);
    W0 = AE_SRAI72(e0,W0,31);
    AE_ADD72X64(e1,W1,W0);
    ACC = AE_SAT64S(e1,W1);
    return ACC;
}
