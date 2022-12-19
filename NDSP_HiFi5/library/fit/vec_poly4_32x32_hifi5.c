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

#include "common.h"
#include "NatureDSP_types.h"
/* Library API */
#include "NatureDSP_Signal_fit.h"
// code optimized for HiFi5

/*-------------------------------------------------------------------------
  Polynomial approximation
  Functions calculate polynomial approximation for all values from given 
  vector. Fixed point functions take polynomial coefficients in Q31 precision. 
  NOTE:
  approximation is calculated like Taylor series that is why overflow may 
  potentially occur if cumulative sum of coefficients given from the last 
  to the first coefficient is bigger that 1. To avoid this negative effect,
  all the coefficients may be scaled down and result will be shifted left 
  after all intermediate computations.

  Precision: 
  32x32  32-bit inputs, 32-bit coefficients, 32-bit output.
  f      floating point

  Input:
  x[N]    input data, Q31 or floating point
  N       length of vector
  lsh     additional left shift for result
  c[M+1]  coefficients (M=4 coefficients for vec_poly4_xxx 
          and M=8 coefficients for vec_poly8_xxx), Q31 or floating point
  Output:			
  z[N]    result, Q31 or floating point

  Restriction:
  x,c,z should not overlap
  lsh   should be in range 0...31

  PERFORMANCE NOTE:
  for optimum performance follow rules:
  x aligned on 16-byte boundary
  N   - multiple of 2
-------------------------------------------------------------------------*/
void vec_poly4_32x32 (int32_t * restrict z, const int32_t * restrict x, const int32_t * restrict c, int lsh,int N )
{
    ae_int32x2 scale=AE_SLAA32S(1,lsh);
          ae_int32x4* restrict pZ=(      ae_int32x4*)z;
    const ae_int32x4* restrict pX=(const ae_int32x4*)x;
    const ae_int32  * restrict pC=(const ae_int32  *)c;
          ae_int32x4* restrict pT;
    ae_valignx2 aX,aZ;
    ae_int32x2  x0,x1,y0,y1;
    ae_int32x2  x2,x3,y2,y3;
    ae_int32x2  t01,t00,t11,t10,t21,t20,t31,t30;
    ae_int32x2  t03,t02,t13,t12,t23,t22,t33,t32;
    ae_int32x2  c0w,c1w,c2w,c3w,c4w;
    int n;
    if (N<=0) return;
    c0w=AE_L32_I(pC,  0);
    c1w=AE_L32_I(pC,  4);
    c2w=AE_L32_I(pC,  8);
    c3w=AE_L32_I(pC, 12);
    c4w=AE_L32_I(pC, 16);
    if (N>=8)
    {
        aX=AE_LA128_PP(pX);
        aZ=AE_ZALIGN128();
        for (n=0;n<(N>>3);n++)
        {
            AE_LA32X2X2_IP(x0,x1,aX, pX); 
            AE_LA32X2X2_IP(x2,x3,aX, pX); 
            AE_L32_IP(c0w,pC,0*sizeof(ae_int32));
            c1w=AE_L32_I(pC,1*sizeof(ae_int32));
            c2w=AE_L32_I(pC,2*sizeof(ae_int32));
            c3w=AE_L32_I(pC,3*sizeof(ae_int32));
      
            AE_MOVD32X4(t31,t30,c3w,c3w);AE_MOVD32X4(t33,t22,c3w,c2w);
            AE_MOVD32X4(t21,t20,c2w,c2w);AE_MOVD32X4(t23,t32,c2w,c3w);
            AE_MOVD32X4(t11,t10,c1w,c1w);AE_MOVD32X4(t13,t02,c1w,c0w);
            AE_MOVD32X4(t01,t00,c0w,c0w);AE_MOVD32X4(t03,t12,c0w,c1w);

            AE_MULAF2P32X4RAS(t30,t31, x0,x1, c4w,c4w);
            AE_MULAF2P32X4RAS(t32,t33, x2,x3, c4w,c4w);
            AE_MULAF2P32X4RAS(t20,t21, x0,x1, t30,t31);
            AE_MULAF2P32X4RAS(t22,t23, x2,x3, t32,t33);
            AE_MULAF2P32X4RAS(t10,t11, x0,x1, t20,t21);
            AE_MULAF2P32X4RAS(t12,t13, x2,x3, t22,t23);
            AE_MULAF2P32X4RAS(t00,t01, x0,x1, t10,t11);
            AE_MULAF2P32X4RAS(t02,t03, x2,x3, t12,t13);
            AE_MUL2P32X4S(y1,y0,t01,t00,scale,scale);
            AE_MUL2P32X4S(y3,y2,t03,t02,scale,scale);

            AE_SA32X2X2_IP(y0,y1, aZ, pZ);
            AE_SA32X2X2_IP(y2,y3, aZ, pZ);
        }
        AE_SA128POS_FP(aZ,pZ);
    }
    N&=7;
    if (N)
    {
        int32_t ALIGN(16) buf[8];
        pT=(ae_int32x4*)buf; AE_S32X2X2_I(0,0,pT,sizeof(ae_int32x4));
        __Pragma("loop_count min=1,max=7")
        for (n=0; n<N; n++) 
        {
            AE_L32_IP(x0,castxcc(ae_int32,pX),sizeof(ae_int32));
            AE_S32_L_IP(x0,castxcc(ae_int32,pT),sizeof(ae_int32));
        }
        pT=(ae_int32x4*)buf;
        AE_L32X2X2_I (x0,x1,pT,0*sizeof(ae_int32x4)); 
        AE_L32X2X2_I (x2,x3,pT,1*sizeof(ae_int32x4)); 
      
        AE_MOVD32X4(t31,t30,c3w,c3w);AE_MOVD32X4(t33,t22,c3w,c2w);
        AE_MOVD32X4(t21,t20,c2w,c2w);AE_MOVD32X4(t23,t32,c2w,c3w);
        AE_MOVD32X4(t11,t10,c1w,c1w);AE_MOVD32X4(t13,t02,c1w,c0w);
        AE_MOVD32X4(t01,t00,c0w,c0w);AE_MOVD32X4(t03,t12,c0w,c1w);
        AE_MULAF2P32X4RAS(t30,t31, x0,x1, c4w,c4w);
        AE_MULAF2P32X4RAS(t32,t33, x2,x3, c4w,c4w);
        AE_MULAF2P32X4RAS(t20,t21, x0,x1, t30,t31);
        AE_MULAF2P32X4RAS(t22,t23, x2,x3, t32,t33);
        AE_MULAF2P32X4RAS(t10,t11, x0,x1, t20,t21);
        AE_MULAF2P32X4RAS(t12,t13, x2,x3, t22,t23);
        AE_MULAF2P32X4RAS(t00,t01, x0,x1, t10,t11);
        AE_MULAF2P32X4RAS(t02,t03, x2,x3, t12,t13);
        AE_MUL2P32X4S(y1,y0,t01,t00,scale,scale);
        AE_MUL2P32X4S(y3,y2,t03,t02,scale,scale);

        AE_S32X2X2_I (y0,y1, pT,0*sizeof(ae_int32x4)); 
        AE_S32X2X2_I (y2,y3, pT,1*sizeof(ae_int32x4)); 
        __Pragma("loop_count min=1,max=7")
        for (n=0; n<N; n++) 
        {
            AE_L32_IP(x0,castxcc(ae_int32,pT),sizeof(ae_int32));
            AE_S32_L_IP(x0,castxcc(ae_int32,pZ),sizeof(ae_int32));
        }
    }
}
