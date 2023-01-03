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
/* Code optimized for HiFi5 core */
#include "NatureDSP_Signal_matop.h"
#include "NatureDSP_types.h"
#include "common.h"
/*-------------------------------------------------------------------------
  Matrix by Vector Multiply
  These functions compute the expression z = 2^lsh * x * y for the matrices 
  x and vector y. 
  NOTE: lsh factor is not relevant for floating point routines.

  Two versions of functions available: regular version (mtx_vecmpy32x32,  
  mtx_vecmpy16x16, mtx_vecmpy8x8, mtx_vecmpy8x16, mtx_vecmpyf) with arbitrary 
  arguments and faster version (mtx_vecmpy32x32_fast, mtx_vecmpy16x16_fast, 
  mtx_vecmpy8x8_fast, mtx_vecmpy8x16_fast,  mtx_vecmpyf_fast) that apply 
  some restrictions

  Precision: 
  32x32 32-bit input, 32-bit output
  16x16 16-bit input, 16-bit output
  8x8   8-bit inputs, 8-bit output
  8x16  8/16-bit inputs, 16-bit output
  f     floating point

  Input:
  x[M*N] input matrix,Q31,Q15 or floating point
  y[N]   input vector,Q31,Q15 or floating point
  M      number of rows in matrix x
  N      number of columns in matrix x
  lsh    additional left shift(applied to the fixed-
         point functions only) 
  Output:
  z[M]   output vector,Q31,Q15 or floating point

  Restriction:
  For regular routines (mtx_vecmpy32x32, mtx_vecmpy16x16, mtx_vecmpy8x8,
  mtx_vecmpy8x16,  mtx_vecmpyf)
  x,y,z should not overlap

  For faster routines  (mtx_vecmpy32x32_fast, mtx_vecmpy16x16_fast, 
  mtx_vecmpy8x8_fast, mtx_vecmpy8x16_fast, mtx_vecmpyf_fast)
  x,y,z   should not overlap
  x,y     aligned on 16-byte boundary
  N, M    multiples of 4
  lsh     should be in range:
          -31...31 for mtx_vecmpy32x32, mtx_vecmpy32x32_fast
          -15...15 for mtx_vecmpy16x16, mtx_vecmpy16x16_fast, 
                   mtx_vecmpy8x8_fast, mtx_vecmpy8x16_fast   
-------------------------------------------------------------------------*/
void mtx_vecmpy16x16_fast(int16_t* restrict z,
    const int16_t* restrict x,
    const int16_t* restrict y,
    int M, int N, int lsh)
#if 1
{

    ae_valign az;
    const ae_int16x4 * restrict px0;
    const ae_int16x4 * restrict px1;
    const ae_int16x4 * restrict px2;
    const ae_int16x4 * restrict py;
    const ae_int16x8 * restrict px0_;
    const ae_int16x8 * restrict px1_;
    const ae_int16x8 * restrict px3_;
    const ae_int16x8 * restrict py_;
    ae_int16x4 * restrict pz;
    ae_int16x4 X0, X1, X2, X3, Y0, Y1;
    ae_int16x4 X4, X5, X6, X7;
    ae_int32x2 Z0, Z1;
    NASSERT_ALIGN(x, HIFI_SIMD_WIDTH);
    NASSERT_ALIGN(y, HIFI_SIMD_WIDTH);
    NASSERT((N & 3) == 0);
    NASSERT((M & 3) == 0);
    int m, n;
    if (N <= 0 || M <= 0)    /* exceptional situation */
    {
        for (m = 0; m < (M>>2); m++) AE_S16X4_IP(0,castxcc(ae_int16x4,z),sizeof(ae_int16x4));
        return;
    }
    py = (const ae_int16x4 *)y;
    pz = (ae_int16x4 *)z;
    az = AE_ZALIGN64();

    px0 = (const ae_int16x4 *)(x);
    px1 = px0 + (N >> 2);
    px2 = px0 + (N >> 1);
    if (N & 4)
    {
        /* optimized for low-bit ordering: N is a multiple of 4, but not a multiple of 8
        Here, 2 successive raws begin in different banks
        */
        ae_valignx2 v0, v1;

        px2 = px1;
        __Pragma("ymemory(px2)")
            __Pragma("loop_count min=1")
            for (m = 0; m < (M >> 2); m++)
            {
                ae_int64 B0, B1, B2, B3;
                py_ = (const ae_int16x8 *)y;

                px0_ = (const ae_int16x8 *)(x+4*m*N);
                px1_ = (const ae_int16x8 *)(x + 4*m*N+N);
                px3_ = (const ae_int16x8 *)(x + 4 * m*N +3*N);
                v0 = AE_LA128_PP(px1_);
                v1 = AE_LA128_PP(px3_);
                
                B0 = B1 = B2 = B3 = AE_ZERO64();
                for (n = 0; n < (N >> 3); n++)
                {
                    AE_L16X4X2_IP(Y0, Y1, py_, 16);

                    AE_L16X4X2_X (X4, X5, px0_, 4 * N);
                    AE_L16X4X2_IP(X0, X1, px0_, 16);
                    AE_LA16X4X2_IP(X2, X3, v0, px1_);
                    AE_LA16X4X2_IP(X6, X7, v1, px3_);

                    AE_MULAAAA2Q16(B0, B1, Y0, Y0, X0, X2);
                    AE_MULAAAA2Q16(B0, B1, Y1, Y1, X1, X3);
                    AE_MULAAAA2Q16(B2, B3, Y0, Y0, X4, X6);
                    AE_MULAAAA2Q16(B2, B3, Y1, Y1, X5, X7);
                }
                py = (const ae_int16x4 *)py_;
                px0 = (const ae_int16x4 *)px0_;
                px1 = (const ae_int16x4 *)px1_;
                X2 = AE_L16X4_X(px0, 4 * N);
                X3 = AE_L16X4_X(px1, 4 * N);
                AE_L16X4_IP(Y0, py, 8);
                AE_L16X4_IP(X0, px0, 8);
                AE_L16X4_IP(X1, px1, 8);
                AE_MULAAAA2Q16(B0, B1, Y0, Y0, X0, X1);
                AE_MULAAAA2Q16(B2, B3, Y0, Y0, X2, X3);

                Z0 = AE_TRUNCA32X2F64S(B0, B1, lsh + 33);
                Z1 = AE_TRUNCA32X2F64S(B2, B3, lsh + 33);
                Y0 = AE_ROUND16X4F32SASYM(Z0, Z1);
                AE_SA16X4_IP(Y0, az, pz);
                px0 += 3 * (N >> 2);
                px2 += 3 * (N >> 2);
            }
    }
    else
    {
        /*        optimized for low-bit ordering: N is a multiple of 8    */
        px0_ = (const ae_int16x8 *)(x);
        px1_ = px0_ + (N >> 3);

        __Pragma("loop_count min=1")
            for (m = 0; m < (M >> 2); m++)
            {
                ae_int64 B0, B1, B2, B3;
                py_ = (const ae_int16x8 *)y;
                B0 = B1 = B2 = B3 = AE_ZERO64();
                __Pragma("loop_count min=1")
                    for (n = 0; n < (N >> 3); n++)
                    {
                        AE_L16X4X2_IP(Y0, Y1, py_, 16);
                        AE_L16X4X2_X(X4, X5, px0_, 4 * N);
                        AE_L16X4X2_X(X6, X7, px1_, 4 * N);
                        AE_L16X4X2_IP(X0, X1, px0_, 16);
                        AE_L16X4X2_IP(X2, X3, px1_, 16);

                        AE_MULAAAA2Q16(B0, B1, Y0, Y0, X0, X2);
                        AE_MULAAAA2Q16(B2, B3, Y0, Y0, X4, X6);
                        AE_MULAAAA2Q16(B0, B1, Y1, Y1, X1, X3);
                        AE_MULAAAA2Q16(B2, B3, Y1, Y1, X5, X7);
                    }
                Z0 = AE_TRUNCA32X2F64S(B0, B1, lsh + 33);
                Z1 = AE_TRUNCA32X2F64S(B2, B3, lsh + 33);
                Y0 = AE_ROUND16X4F32SASYM(Z0, Z1);
                AE_SA16X4_IP(Y0, az, pz);
                px0_ += 3 * (N >> 3);
                px1_ += 3 * (N >> 3);
            }
    }
    AE_SA64POS_FP(az, pz);
}
#else
// this variant is bit better for very big N but for moderate N it is slower due to the loop overheads
{
    const ae_int16x8* restrict pX0;
    const ae_int16x8* restrict pX1;
    const ae_int16x8* restrict pW0;
    const ae_int16x8* restrict pW1;
    const ae_int16x8* restrict pW3;
    const ae_int16x8* restrict pW5;
    const ae_int16x8* restrict pW7;
    const ae_int16x8* restrict pY;
    ae_valign aZ;
    ae_int16x4 *pZ;

    int m,n;
    NASSERT_ALIGN(x, HIFI_SIMD_WIDTH);
    NASSERT_ALIGN(y, HIFI_SIMD_WIDTH);
    NASSERT((N & 3) == 0);
    NASSERT((M & 3) == 0);
    if (N <= 0 || M <= 0)    /* exceptional situation */
    {
        for (m = 0; m < (M>>2); m++) AE_S16X4_IP(0,castxcc(ae_int16x4,z),sizeof(ae_int16x4));
        return;
    }

    NASSERT(lsh >= -15 && lsh <= 15);
    pZ=(ae_int16x4*)z;
    switch (N&12)
    {
    default:
        aZ=AE_ZALIGN64();
        for(m=0; m<(M>>3); m++)
        {
            ae_int64 a0,a1,a2,a3,a4,a5,a6,a7;
            ae_int16x4 y0,y1;
            ae_int16x4 x0,x1,x2,x3,x4,x5,x6,x7;
            ae_int16x4 x8,x9,xa,xb,xc,xd,xe,xf;
            ae_valignx2 aW1,aW3,aW5,aW7;
            pY =(const ae_int16x8*)(y);
            pW0=(const ae_int16x8*)(x);
            pW1=(const ae_int16x8*)XT_ADDX2(N,(uintptr_t)pW0);
            pW3=(const ae_int16x8*)XT_ADDX4(N,(uintptr_t)pW1);
            pW5=(const ae_int16x8*)XT_ADDX8(N,(uintptr_t)pW1);
            pW7=(const ae_int16x8*)XT_ADDX8(N,(uintptr_t)pW3);
            AE_MOVDX2(a0,a1,0,0); AE_MOVDX2(a2,a3,0,0); 
            AE_MOVDX2(a4,a5,0,0); AE_MOVDX2(a6,a7,0,0);
            aW1=AE_LA128_PP(pW1);
            aW3=AE_LA128_PP(pW3);
            aW5=AE_LA128_PP(pW5);
            aW7=AE_LA128_PP(pW7);
            for(n=0; n<(N>>3); n++)
            {
                AE_L16X4X2_IP (y0,y1,pY ,sizeof(ae_int16x8));
                AE_LA16X4X2_IP(xe,xf,aW7,pW7);
                AE_L16X4X2_X  (xc,xd,pW0,6*N*sizeof(int16_t));
                AE_LA16X4X2_IP(xa,xb,aW5,pW5);
                AE_L16X4X2_X  (x8,x9,pW0,4*N*sizeof(int16_t));
                AE_LA16X4X2_IP(x6,x7,aW3,pW3);
                AE_L16X4X2_X  (x4,x5,pW0,2*N*sizeof(int16_t));
                AE_LA16X4X2_IP(x2,x3,aW1,pW1);
                AE_L16X4X2_IP (x0,x1,pW0,sizeof(ae_int16x8));

                AE_MULAAAA2Q16(a0,a1,x0,x2,y0,y0);
                AE_MULAAAA2Q16(a0,a1,x1,x3,y1,y1);
                AE_MULAAAA2Q16(a2,a3,x4,x6,y0,y0);
                AE_MULAAAA2Q16(a2,a3,x5,x7,y1,y1);
                AE_MULAAAA2Q16(a4,a5,x8,xa,y0,y0);
                AE_MULAAAA2Q16(a4,a5,x9,xb,y1,y1);
                AE_MULAAAA2Q16(a6,a7,xc,xe,y0,y0);
                AE_MULAAAA2Q16(a6,a7,xd,xf,y1,y1);
            }
            AE_L16X4_IP (y0,castxcc(ae_int16x4,pY ),sizeof(ae_int16x4));
            AE_L16X4_IP (xe,castxcc(ae_int16x4,pW7),sizeof(ae_int16x4));
            xc=AE_L16X4_X  ((const ae_int16x4*)pW0 ,6*N*sizeof(int16_t));
            AE_L16X4_IP (xa,castxcc(ae_int16x4,pW5),sizeof(ae_int16x4));
            x8=AE_L16X4_X  ((const ae_int16x4*)pW0 ,4*N*sizeof(int16_t));
            AE_L16X4_IP (x6,castxcc(ae_int16x4,pW3),sizeof(ae_int16x4));
            x4=AE_L16X4_X  ((const ae_int16x4*)pW0 ,2*N*sizeof(int16_t));
            AE_L16X4_IP (x2,castxcc(ae_int16x4,pW1),sizeof(ae_int16x4));
            AE_L16X4_IP (x0,castxcc(ae_int16x4,pW0),sizeof(ae_int16x4));

            AE_MULAAAA2Q16(a0,a1,x0,x2,y0,y0);
            AE_MULAAAA2Q16(a2,a3,x4,x6,y0,y0);
            AE_MULAAAA2Q16(a4,a5,x8,xa,y0,y0);
            AE_MULAAAA2Q16(a6,a7,xc,xe,y0,y0);
            AE_SA16X4_IP(AE_ROUND16X4F32SASYM(AE_TRUNCA32X2F64S(a0, a1, lsh + 33), AE_TRUNCA32X2F64S(a2, a3, lsh + 33)),aZ,pZ);
            AE_SA16X4_IP(AE_ROUND16X4F32SASYM(AE_TRUNCA32X2F64S(a4, a5, lsh + 33), AE_TRUNCA32X2F64S(a6, a7, lsh + 33)),aZ,pZ);
            x+=8*N;
        }
        if(M&4)
        {
            ae_int64 a0,a1,a2,a3;
            ae_int16x4 y0,y1;
            ae_int16x4 x0,x1,x2,x3,x4,x5,x6,x7;
            ae_valignx2 aW1,aW3;
            pY =(const ae_int16x8*)(y);
            pW0=(const ae_int16x8*)(x);
            pW1=(const ae_int16x8*)XT_ADDX2(N,(uintptr_t)pW0);
            pW3=(const ae_int16x8*)XT_ADDX4(N,(uintptr_t)pW1);
            AE_MOVDX2(a0,a1,0,0); AE_MOVDX2(a2,a3,0,0); 
            aW1=AE_LA128_PP(pW1);
            aW3=AE_LA128_PP(pW3);
            for(n=0; n<(N>>3); n++)
            {
                AE_L16X4X2_IP (y0,y1,pY ,sizeof(ae_int16x8));
                AE_LA16X4X2_IP(x6,x7,aW3,pW3);
                AE_L16X4X2_X  (x4,x5,pW0,2*N*sizeof(int16_t));
                AE_LA16X4X2_IP(x2,x3,aW1,pW1);
                AE_L16X4X2_IP (x0,x1,pW0,sizeof(ae_int16x8));

                AE_MULAAAA2Q16(a0,a1,x0,x2,y0,y0);
                AE_MULAAAA2Q16(a0,a1,x1,x3,y1,y1);
                AE_MULAAAA2Q16(a2,a3,x4,x6,y0,y0);
                AE_MULAAAA2Q16(a2,a3,x5,x7,y1,y1);
            }
            AE_L16X4_IP (y0,castxcc(ae_int16x4,pY) ,sizeof(ae_int16x4));
            AE_L16X4_IP (x6,castxcc(ae_int16x4,pW3),sizeof(ae_int16x4));
            x4=AE_L16X4_X  ((const ae_int16x4*)pW0 ,2*N*sizeof(int16_t));
            AE_L16X4_IP (x2,castxcc(ae_int16x4,pW1),sizeof(ae_int16x4));
            AE_L16X4_IP (x0,castxcc(ae_int16x4,pW0),sizeof(ae_int16x4));
            AE_MULAAAA2Q16(a0,a1,x0,x2,y0,y0);
            AE_MULAAAA2Q16(a2,a3,x4,x6,y0,y0);
            AE_SA16X4_IP(AE_ROUND16X4F32SASYM(AE_TRUNCA32X2F64S(a0, a1, lsh + 33), AE_TRUNCA32X2F64S(a2, a3, lsh + 33)),aZ,pZ);
            x+=4*N;

        }
        AE_SA64POS_FP(aZ,pZ);
        break;
    case 8:
        // N is a multiple of 8
        NASSERT(N%8==0);
        aZ=AE_ZALIGN64();
        __Pragma("loop_count min=1")
        for(m=0; m<(M>>2); m++)
        {
            ae_int16x4 y0,y1,y2,y3;
            ae_int16x4 x0,x1,x2,x3,x4,x5,x6,x7;
            ae_int16x4 z0,z1,z2,z3,z4,z5,z6,z7;
            ae_int64 a0,a1,a2,a3,b0,b1,b2,b3;
            pX0=(const ae_int16x8*)(x);
            pY =(const ae_int16x8*)(y);
            AE_L16X4X2_IP(y0,y1,pY ,sizeof(ae_int16x8));
            AE_L16X4X2_X (x6,x7,pX0,3*N*sizeof(int16_t));
            AE_L16X4X2_X (x4,x5,pX0,2*N*sizeof(int16_t));
            AE_L16X4X2_X (x2,x3,pX0,1*N*sizeof(int16_t));
            AE_L16X4X2_IP(x0,x1,pX0,sizeof(ae_int16x8));
            AE_MOVDX2(b0,b1,0,0);AE_MOVDX2(b2,b3,0,0);
            AE_MULZAAAA2Q16(a0,a1,x0,x2,y0,y0);
            AE_MULZAAAA2Q16(a2,a3,x4,x6,y0,y0);
            AE_MULAAAA2Q16(a0,a1,x1,x3,y1,y1);
            AE_MULAAAA2Q16(a2,a3,x5,x7,y1,y1);
            pX1=pX0+1;
            __Pragma("ymemory(pX1)");
            for(n=0; n<(N>>4); n++)
            {
                AE_L16X4X2_I (y2,y3,pY ,  sizeof(ae_int16x8));
                AE_L16X4X2_IP(y0,y1,pY ,2*sizeof(ae_int16x8));
                AE_L16X4X2_X (x6,x7,pX0,3*N*sizeof(int16_t));
                AE_L16X4X2_X (z6,z7,pX1,3*N*sizeof(int16_t));
                AE_L16X4X2_X (x4,x5,pX0,2*N*sizeof(int16_t));
                AE_L16X4X2_X (z4,z5,pX1,2*N*sizeof(int16_t));
                AE_L16X4X2_X (x2,x3,pX0,1*N*sizeof(int16_t));
                AE_L16X4X2_X (z2,z3,pX1,1*N*sizeof(int16_t));
                AE_L16X4X2_IP(x0,x1,pX0,2*sizeof(ae_int16x8));
                AE_L16X4X2_IP(z0,z1,pX1,2*sizeof(ae_int16x8));
                AE_MULAAAA2Q16(a2,a3,x5,x7,y1,y1);
                AE_MULAAAA2Q16(b2,b3,z5,z7,y3,y3);
                AE_MULAAAA2Q16(a2,a3,x4,x6,y0,y0);
                AE_MULAAAA2Q16(b2,b3,z4,z6,y2,y2);
                AE_MULAAAA2Q16(a0,a1,x1,x3,y1,y1);
                AE_MULAAAA2Q16(b0,b1,z1,z3,y3,y3);
                AE_MULAAAA2Q16(a0,a1,x0,x2,y0,y0);
                AE_MULAAAA2Q16(b0,b1,z0,z2,y2,y2);
            }
            AE_SA16X4_IP(AE_ROUND16X4F32SASYM(AE_TRUNCA32X2F64S(a0+b0, a1+b1, lsh + 33), AE_TRUNCA32X2F64S(a2+b2, a3+b3, lsh + 33)),aZ,pZ);
            x+=4*N;
        }
        AE_SA64POS_FP(aZ,pZ);
        break;
    case 0:
        // N is a multiple of 16
        NASSERT(N%16==0);
        aZ=AE_ZALIGN64();
        __Pragma("loop_count min=1")
        for(m=0; m<(M>>2); m++)
        {
            int64_t A0,A1,A2,A3;
            ae_int64 a0,a1,a2,a3;
            A0=A1=A2=A3=0;
            pX0=(const ae_int16x8*)(x);
            pX1=pX0+1;
            pY =(const ae_int16x8*)(y);
            AE_MOVDX2(a0,a1,0,0); AE_MOVDX2(a2,a3,0,0);
            __Pragma("ymemory(pX1)");
            __Pragma("loop_count min=1")
            for(n=0; n<(N>>4); n++)
            {
                ae_int16x4 y0,y1,y2,y3;
                ae_int16x4 x0,x1,x2,x3,x4,x5,x6,x7;
                ae_int16x4 z0,z1,z2,z3,z4,z5,z6,z7;
                AE_L16X4X2_I (y2,y3,pY ,  sizeof(ae_int16x8));
                AE_L16X4X2_IP(y0,y1,pY ,2*sizeof(ae_int16x8));
                AE_L16X4X2_X (x6,x7,pX0,3*N*sizeof(int16_t));
                AE_L16X4X2_X (x4,x5,pX0,2*N*sizeof(int16_t));
                AE_L16X4X2_X (x2,x3,pX0,1*N*sizeof(int16_t));
                AE_L16X4X2_IP(x0,x1,pX0,2*sizeof(ae_int16x8));
                AE_L16X4X2_X (z6,z7,pX1,3*N*sizeof(int16_t));
                AE_L16X4X2_X (z4,z5,pX1,2*N*sizeof(int16_t));
                AE_L16X4X2_X (z2,z3,pX1,1*N*sizeof(int16_t));
                AE_L16X4X2_IP(z0,z1,pX1,2*sizeof(ae_int16x8));
                AE_MULAAAA2Q16(a2,a3,x5,x7,y1,y1);
                AE_MULAAAA2Q16(a2,a3,x4,x6,y0,y0);
                AE_MULAAAA2Q16(a2,a3,z5,z7,y3,y3);
                AE_MULAAAA2Q16(a2,a3,z4,z6,y2,y2);
                AE_MULAAAA2Q16(a0,a1,x1,x3,y1,y1);
                AE_MULAAAA2Q16(a0,a1,x0,x2,y0,y0);
                AE_MULAAAA2Q16(a0,a1,z1,z3,y3,y3);
                AE_MULAAAA2Q16(a0,a1,z0,z2,y2,y2);
            }
            AE_SA16X4_IP(AE_ROUND16X4F32SASYM(AE_TRUNCA32X2F64S(a0, a1, lsh + 33), AE_TRUNCA32X2F64S(a2, a3, lsh + 33)),aZ,pZ);
            x+=4*N;
        }
        AE_SA64POS_FP(aZ,pZ);
        break;
    }

}
#endif
