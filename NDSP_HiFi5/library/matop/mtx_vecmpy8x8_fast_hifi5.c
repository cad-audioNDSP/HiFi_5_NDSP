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
#include "NatureDSP_Signal_matop.h"
#include "NatureDSP_types.h"
#include "common.h"
/* Code optimized for HiFi5 core */

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
#if !defined(AE_MULA8Q8X8)
void mtx_vecmpy8x8_fast( int8_t* restrict z,
                 const int8_t* restrict x,
                 const int8_t* restrict y,
                 int M, int N, int lsh)
{
    const int8_t* restrict pY;
    const int8_t* restrict pX0;
    const int8_t* restrict pX1;
    const int8_t* restrict pX2;
    const int8_t* restrict pX3;
    int sa;
    int m,n;
    NASSERT_ALIGN(z,HIFI_SIMD_WIDTH);
    NASSERT_ALIGN(x,HIFI_SIMD_WIDTH);
    NASSERT_ALIGN(y,HIFI_SIMD_WIDTH);
    NASSERT(M%4==0);
    NASSERT(N%4==0);
    NASSERT(lsh >= -15 && lsh <= 15);
    if (M<=0) return;
    if (N<=0)    /* exceptional situation */
    {
        for (m=0; m<M; m++) z[m]=0;
        return ;
    }
    sa=49+lsh;
    if ((N&15)==0)
    {   // N is a multiple of 16
        pX0=x;
        for (m=0; m<M;m+=4)
        {
            ae_int8x8 x0,x1,x2,x3,x4,x5,x6,x7,y0,y1;
            ae_int64 A0,A1,A2,A3;
            ae_int64 B0,B1,B2,B3;
            ae_int8x8 r;
            ae_int32x2 res;
            pY=y;
            AE_MOVDX2(A0,B0,0,0);
            AE_MOVDX2(A1,B1,0,0);
            AE_MOVDX2(A2,B2,0,0);
            AE_MOVDX2(A3,B3,0,0);
            __Pragma("loop_count min=1")
            for (n=0; n<(N>>4); n++)
            {
                AE_L8X8X2_X(x2,x3,(const ae_int8x16*)pX0,1*N);
                AE_L8X8X2_X(x4,x5,(const ae_int8x16*)pX0,2*N);
                AE_L8X8X2_X(x6,x7,(const ae_int8x16*)pX0,3*N);
                AE_L8X8X2_IP(x0,x1,castxcc(ae_int8x16,pX0),sizeof(ae_int8x16));
                AE_L8X8X2_IP(y0,y1,castxcc(ae_int8x16,pY),sizeof(ae_int8x16));
                AE_MULAAAA2Q8(A0,B0,x0,y0);
                AE_MULAAAA2Q8(A1,B1,x2,y0);
                AE_MULAAAA2Q8(A2,B2,x4,y0);
                AE_MULAAAA2Q8(A3,B3,x6,y0);
                AE_MULAAAA2Q8(A0,B0,x1,y1);
                AE_MULAAAA2Q8(A1,B1,x3,y1);
                AE_MULAAAA2Q8(A2,B2,x5,y1);
                AE_MULAAAA2Q8(A3,B3,x7,y1);
            }
            pX0+=3*N;
            r=AE_ROUND8X4F32SASYM_L(AE_TRUNCA32X2F64S(A3+B3, A2+B2, sa),AE_TRUNCA32X2F64S(A1+B1, A0+B0, sa));
            res=AE_MOVINT32X2_FROMINT8X8(r);
            AE_S32_L_IP(res,castxcc(ae_int32,z),sizeof(int32_t));
        }
        return;
    }
    if ((N&7)==0)
    {   // N is a multiple of 8
        pX0=x+0*N;
        pX1=x+1*N;
        pX2=x+2*N;
        pX3=x+3*N;
        for (m=0; m<M;m+=4)
        {
            ae_int8x8 x0,x1,x2,x3,y0;
            ae_int64 A0,A1,A2,A3;
            ae_int64 B0,B1,B2,B3;
            ae_int8x8 r;
            ae_int32x2 res;
            pY=y;
            AE_MOVDX2(A0,B0,0,0);
            AE_MOVDX2(A1,B1,0,0);
            AE_MOVDX2(A2,B2,0,0);
            AE_MOVDX2(A3,B3,0,0);
            __Pragma("loop_count min=1")
            for (n=0; n<(N>>3); n++)
            {
                AE_L8X8_IP(y0,castxcc(ae_int8x8,pY),sizeof(ae_int8x8));
                AE_L8X8_IP(x0,castxcc(ae_int8x8,pX0),sizeof(ae_int8x8));
                AE_L8X8_IP(x1,castxcc(ae_int8x8,pX1),sizeof(ae_int8x8));
                AE_L8X8_IP(x2,castxcc(ae_int8x8,pX2),sizeof(ae_int8x8));
                AE_L8X8_IP(x3,castxcc(ae_int8x8,pX3),sizeof(ae_int8x8));
                AE_MULAAAA2Q8(A0,B0,x0,y0);
                AE_MULAAAA2Q8(A1,B1,x1,y0);
                AE_MULAAAA2Q8(A2,B2,x2,y0);
                AE_MULAAAA2Q8(A3,B3,x3,y0);
            }
            pX0+=3*N;
            pX1+=3*N;
            pX2+=3*N;
            pX3+=3*N;
            r=AE_ROUND8X4F32SASYM_L(AE_TRUNCA32X2F64S(A3+B3, A2+B2, sa),AE_TRUNCA32X2F64S(A1+B1, A0+B0, sa));
            res=AE_MOVINT32X2_FROMINT8X8(r);
            AE_S32_L_IP(res,castxcc(ae_int32,z),sizeof(int32_t));
        }
        return;
    }
#if 0
    {   // N is a multiple of 4
        pX0=x+0*N;
        pX1=x+1*N;
        pX2=x+2*N;
        pX3=x+3*N;
        for (m=0; m<M;m+=4)
        {
            ae_int8x8 x0,x1,x2,x3,y0;
            ae_int64 A0,A1,A2,A3;
            ae_int64 B0,B1,B2,B3;
            ae_int8x8 r;
            ae_int32x2 res;
            ae_valign aX1,aX3;
            pY=y;
            AE_MOVDX2(A0,B0,0,0);
            AE_MOVDX2(A1,B1,0,0);
            AE_MOVDX2(A2,B2,0,0);
            AE_MOVDX2(A3,B3,0,0);
            NASSERT_ALIGN8(pX0);
            NASSERT_ALIGN8(pX2);
            aX1=AE_LA64_PP(pX1);
            aX3=AE_LA64_PP(pX3);
            for (n=0; n<(N>>3); n++)
            {
                AE_L8X8_IP(y0,castxcc(ae_int8x8,pY),sizeof(ae_int8x8));
                AE_L8X8_IP(x0,castxcc(ae_int8x8,pX0),sizeof(ae_int8x8));
                AE_LA8X8_IP(x1,aX1,castxcc(ae_int8x8,pX1));
                AE_L8X8_IP(x2,castxcc(ae_int8x8,pX2),sizeof(ae_int8x8));
                AE_LA8X8_IP(x3,aX3,castxcc(ae_int8x8,pX3));
                AE_MULAAAA2Q8(A0,B0,x0,y0);
                AE_MULAAAA2Q8(A1,B1,x1,y0);
                AE_MULAAAA2Q8(A2,B2,x2,y0);
                AE_MULAAAA2Q8(A3,B3,x3,y0);
            }
            if (N&4)
            {
                ae_int32x2 t;
                t=AE_L32_I((const ae_int32*)pY,0); AE_MOVT32X2(t,0,AE_MOVAB2(2)); y0=AE_MOVINT8X8_FROMINT32X2(t);  // mask half of vector
                AE_L32_IP(t,castxcc(ae_int32,pX0),sizeof(ae_int32)); x0=AE_MOVINT8X8_FROMINT32X2(t);
                AE_L32_IP(t,castxcc(ae_int32,pX1),sizeof(ae_int32)); x1=AE_MOVINT8X8_FROMINT32X2(t);
                AE_L32_IP(t,castxcc(ae_int32,pX2),sizeof(ae_int32)); x2=AE_MOVINT8X8_FROMINT32X2(t);
                AE_L32_IP(t,castxcc(ae_int32,pX3),sizeof(ae_int32)); x3=AE_MOVINT8X8_FROMINT32X2(t);
                AE_MULAAAA2Q8(A0,B0,x0,y0);
                AE_MULAAAA2Q8(A1,B1,x1,y0);
                AE_MULAAAA2Q8(A2,B2,x2,y0);
                AE_MULAAAA2Q8(A3,B3,x3,y0);
            }
            pX0+=3*N;
            pX1+=3*N;
            pX2+=3*N;
            pX3+=3*N;
            r=AE_ROUND8X4F32SASYM_L(AE_TRUNCA32X2F64S(A3+B3, A2+B2, sa),AE_TRUNCA32X2F64S(A1+B1, A0+B0, sa));
            res=AE_MOVINT32X2_FROMINT8X8(r);
            AE_S32_L_IP(res,castxcc(ae_int32,z),sizeof(int32_t));
        }
    }
#else
    {   // N is a multiple of 4
        pX0=x+0*N;
        pX1=x+1*N;
        pX2=x+2*N;
        pX3=x+3*N;
        for (m=0; m<M;m+=4)
        {
            ae_int8x8 x0,x1,x2,x3,y0;
            ae_int64 A0,A1,A2,A3;
            ae_int64 B0,B1,B2,B3;
            ae_int8x8 r;
            ae_int32x2 res;
            pY=y;
            AE_MOVDX2(A0,B0,0,0);
            AE_MOVDX2(A1,B1,0,0);
            AE_MOVDX2(A2,B2,0,0);
            AE_MOVDX2(A3,B3,0,0);
            NASSERT_ALIGN8(pX0);
            NASSERT_ALIGN8(pX2);
            for (n=0; n<(N>>3); n++)
            {
                ae_int32x2 t,t0,t1;
                AE_L32X2_IP(t,castxcc(ae_int32x2,pY),sizeof(ae_int8x8));  y0 =AE_MOVINT8X8_FROMINT32X2(t);
                AE_L32X2_IP(t,castxcc(ae_int32x2,pX0),sizeof(ae_int8x8)); x0 =AE_MOVINT8X8_FROMINT32X2(t);
                t1=AE_L32_I (  (const ae_int32*)pX1 ,sizeof(ae_int32));
                AE_L32_IP  (t0,castxcc(ae_int32,pX1),sizeof(ae_int8x8));
                t0=AE_SEL32_HH(t0,t1);x1 =AE_MOVINT8X8_FROMINT32X2(t0);
                AE_L32X2_IP(t,castxcc(ae_int32x2,pX2),sizeof(ae_int8x8)); x2 =AE_MOVINT8X8_FROMINT32X2(t);
                t1=AE_L32_I (  (const ae_int32*)pX3 ,sizeof(ae_int32));
                AE_L32_IP  (t0,castxcc(ae_int32,pX3),sizeof(ae_int8x8));
                t0=AE_SEL32_HH(t0,t1);x3 =AE_MOVINT8X8_FROMINT32X2(t0);
                AE_MULAAAA2Q8(A0,B0,x0,y0);
                AE_MULAAAA2Q8(A1,B1,x1,y0);
                AE_MULAAAA2Q8(A2,B2,x2,y0);
                AE_MULAAAA2Q8(A3,B3,x3,y0);
            }
            if (N&4)
            {
                ae_int32x2 t;
                t=AE_L32_I((const ae_int32*)pY,0); AE_MOVT32X2(t,0,AE_MOVAB2(2)); y0=AE_MOVINT8X8_FROMINT32X2(t);  // mask half of vector
                AE_L32_IP(t,castxcc(ae_int32,pX0),sizeof(ae_int32)); x0=AE_MOVINT8X8_FROMINT32X2(t);
                AE_L32_IP(t,castxcc(ae_int32,pX1),sizeof(ae_int32)); x1=AE_MOVINT8X8_FROMINT32X2(t);
                AE_L32_IP(t,castxcc(ae_int32,pX2),sizeof(ae_int32)); x2=AE_MOVINT8X8_FROMINT32X2(t);
                AE_L32_IP(t,castxcc(ae_int32,pX3),sizeof(ae_int32)); x3=AE_MOVINT8X8_FROMINT32X2(t);
                AE_MULAAAA2Q8(A0,B0,x0,y0);
                AE_MULAAAA2Q8(A1,B1,x1,y0);
                AE_MULAAAA2Q8(A2,B2,x2,y0);
                AE_MULAAAA2Q8(A3,B3,x3,y0);
            }
            pX0+=3*N;
            pX1+=3*N;
            pX2+=3*N;
            pX3+=3*N;
            r=AE_ROUND8X4F32SASYM_L(AE_TRUNCA32X2F64S(A3+B3, A2+B2, sa),AE_TRUNCA32X2F64S(A1+B1, A0+B0, sa));
            res=AE_MOVINT32X2_FROMINT8X8(r);
            AE_S32_L_IP(res,castxcc(ae_int32,z),sizeof(int32_t));
        }
    }
#endif
}
#else
void mtx_vecmpy8x8_fast( int8_t* restrict z,
                 const int8_t* restrict x,
                 const int8_t* restrict y,
                 int M, int N, int lsh)
{
    ae_valignx2 alX1, alX2, alX3;
    int8_t * restrict pZ0;
    const int8_t * restrict pX0;
    const int8_t * restrict pX1;
    const int8_t * restrict pX2;
    const int8_t * restrict pX3;
    const int8_t * restrict pY0;
    ae_int32x2 Z0, Z1;
    ae_int16x4 t0;
    ae_int8x8 x00, x10, x20, x30, x01, x11, x21, x31;
    ae_int8x8 y0, y1, z0;
    int m,n;

    NASSERT(lsh >= -15 && lsh <= 15);
    NASSERT_ALIGN(z,HIFI_SIMD_WIDTH);
    NASSERT_ALIGN(x,HIFI_SIMD_WIDTH);
    NASSERT_ALIGN(y,HIFI_SIMD_WIDTH);
    NASSERT(M%4==0);
    NASSERT(N%4==0);
    if (M<=0) return;
    if (N<=0)    /* exceptional situation */
    {
        for (m=0; m<M; m++) *z++ = 0;
        return;
    }

    pZ0 = z;
    if (N & 15) /* N is not a multiple of 16 */
    {
        __Pragma("loop_count min=1")
        for (m=0; m<(M>>2);m++)
        {
            pX0 = x + m*N*4;
            pX1 = pX0 + N;
            pX2 = pX1 + N;
            pX3 = pX2 + N;

            pY0 = y;

            alX1 = AE_LA128_PP(pX1);
            alX2 = AE_LA128_PP(pX2);
            alX3 = AE_LA128_PP(pX3);
            Z0 = Z1 = AE_ZERO32();

            for (n=0; n<(N>>4); n++)
            {
                /* load x matrix, 4x16 values, 8-bit */
                AE_L8X8X2_IP (x00, x01,       castxcc(ae_int8x16,pX0), 16);
                AE_LA8X8X2_IP(x10, x11, alX1, castxcc(ae_int8x16,pX1));
                AE_LA8X8X2_IP(x20, x21, alX2, castxcc(ae_int8x16,pX2));
                AE_LA8X8X2_IP(x30, x31, alX3, castxcc(ae_int8x16,pX3));

                /* load y vector, 16 values, 8-bit */
                AE_L8X8X2_IP(y0, y1, castxcc(ae_int8x16,pY0), 16);

                /* make multiply, using 4-way 32-bit output octal mac, 8X8 bit */
                AE_MULA8Q8X8(Z0, Z1, x10, x00, x30, x20, y0);
                AE_MULA8Q8X8(Z0, Z1, x11, x01, x31, x21, y1);
            }
            {
                ae_valignx2 alX0, alY0;
                alX0 = AE_LA128_PP(pX0);
                alY0 = AE_LA128_PP(pY0);
                /* load x matrix, 4x16 values, 8-bit */
                AE_LAV8X8X2_XP(x00, x01, alX0, castxcc(ae_int8x16,pX0), N&15);
                AE_LAV8X8X2_XP(x10, x11, alX1, castxcc(ae_int8x16,pX1), N&15);
                AE_LAV8X8X2_XP(x20, x21, alX2, castxcc(ae_int8x16,pX2), N&15);
                AE_LAV8X8X2_XP(x30, x31, alX3, castxcc(ae_int8x16,pX3), N&15);

                /* load y vector, 4...12 values, 8-bit */
                AE_LAV8X8X2_XP(y0, y1, alY0, castxcc(ae_int8x16,pY0), N&15);

                /* make multiply, using 4-way 32-bit output octal mac, 8X8 bit */
                AE_MULA8Q8X8(Z0, Z1, x10, x00, x30, x20, y0);
                AE_MULA8Q8X8(Z0, Z1, x11, x01, x31, x21, y1);
            }
            /* Q15 + lsh <- Q14 + 1 + lsh w/ truncation and saturation */
            t0 = AE_TRUNCA16X4F32S(Z1, Z0, 16+1+lsh);
            /* Q7 + lsh <- Q15 + lsh -8 w/ rounding and saturation */
            z0 = AE_ROUND8X8F16SASYM(t0, t0);
            AE_S32_H_IP(AE_MOVINT32X2_FROMINT8X8(z0), castxcc(ae_int32,pZ0), 4);
        }
    }
    else
    {
        __Pragma("loop_count min=1")
        for (m=0; m<(M>>2);m++)
        {
            pX0 = x + m*N*4;
            pY0 = y;

            Z0 = Z1 = AE_ZERO32();

            __Pragma("loop_count min=1")
            for (n=0; n<(N>>4); n++)
            {
                /* load x matrix, 4x16 values, 8-bit */
                AE_L8X8X2_XP(x00, x01, castxcc(ae_int8x16,pX0), N);
                AE_L8X8X2_XP(x10, x11, castxcc(ae_int8x16,pX0), N);
                AE_L8X8X2_XP(x20, x21, castxcc(ae_int8x16,pX0), N);
                AE_L8X8X2_XP(x30, x31, castxcc(ae_int8x16,pX0), -3*N+16);

                /* load y vector, 16 values, 8-bit */
                AE_L8X8X2_IP(y0, y1, castxcc(ae_int8x16,pY0), 16);

                /* make multiply, using 4-way 32-bit output octal mac, 8X8 bit */
                AE_MULA8Q8X8(Z0, Z1, x10, x00, x30, x20, y0);
                AE_MULA8Q8X8(Z0, Z1, x11, x01, x31, x21, y1);
            }
            /* Q15 + lsh <- Q14 + 1 + lsh w/ truncation and saturation */
            t0 = AE_TRUNCA16X4F32S(Z1, Z0, 16+1+lsh);
            /* Q7 + lsh <- Q15 + lsh -8 w/ rounding and saturation */
            z0 = AE_ROUND8X8F16SASYM(t0, t0);
            AE_S32_H_IP(AE_MOVINT32X2_FROMINT8X8(z0), castxcc(ae_int32,pZ0), 4);
        }
    }
}
#endif
