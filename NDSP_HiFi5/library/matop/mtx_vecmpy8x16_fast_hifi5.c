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
#if !defined(AE_MULA8Q8X16)
void mtx_vecmpy8x16_fast( int16_t* restrict z,
                 const int8_t* restrict x,
                 const int16_t* restrict y,
                 int M, int N, int lsh)
{
    static const uint64_t dsel_rephhll_tbl=0x40516273c8d9eafbULL;
    const ae_int8x8 * restrict pX0;
    const ae_int8x8 * restrict pX1;
    const ae_int8x8 * restrict pX2;
    const ae_int8x8 * restrict pX3;
    const ae_int16x4 * restrict pY;
    ae_valign aX1,aX3;
    int m,n;
    int sa=41+lsh;
    ae_int8x8 dsel_rephhll;
    NASSERT_ALIGN(z,HIFI_SIMD_WIDTH);
    NASSERT_ALIGN(x,HIFI_SIMD_WIDTH);
    NASSERT_ALIGN(y,HIFI_SIMD_WIDTH);
    NASSERT(M%4==0);
    NASSERT(N%4==0);
    NASSERT(lsh >= -15 && lsh <= 15);
    if (N<=0)    /* exceptional situation */
    {
        if (M<=0) return;
        for (m=0; m<(M>>3); m++) AE_S16X4X2_IP(0,0,castxcc(ae_int16x8,z),sizeof(ae_int16x8));
        if (M&4) AE_S16X4_IP(0,castxcc(ae_int16x4,z),sizeof(ae_int16x4));
        return ;
    }
    dsel_rephhll=AE_L8X8_I((const ae_int8x8*)&dsel_rephhll_tbl,0);
    pX0=(const ae_int8x8*)( x );
    pX1=(const ae_int8x8*)( (uintptr_t)x + N);
    pX2=(const ae_int8x8*)XT_ADDX2(N, (uintptr_t)pX0);
    pX3=(const ae_int8x8*)XT_ADDX2(N, (uintptr_t)pX1);
    if (N&4)
    {   
        // separate variant for N not a multiple of 8
        for (m=0; m<(M>>2); m++)
        {
            ae_int16x4 r;
            ae_int64 A0,A1,A2,A3;
            AE_MOVDX2(A0,A1,0,0); AE_MOVDX2(A2,A3,0,0);
            pY=(const ae_int16x4*)y;
            aX1=AE_LA64_PP(pX1);
            aX3=AE_LA64_PP(pX3);
            for (n=0; n<(N>>3); n++)
            {
                ae_int16x4 y0,y1;
                ae_int8x8 x0,x1,x2,x3,x00,x11,x22,x33;
                AE_L16X4X2_IP(y0,y1,castxcc(ae_int16x8,pY),sizeof(ae_int16x8));
                AE_L8X8_IP (x0,pX0,sizeof(ae_int8x8));
                AE_LA8X8_IP(x1,aX1,pX1);
                AE_L8X8_IP (x2,pX2,sizeof(ae_int8x8));
                AE_LA8X8_IP(x3,aX3,pX3);
                AE_DSEL8X8(x00,x11,x0,x1,dsel_rephhll);
                AE_DSEL8X8(x22,x33,x2,x3,dsel_rephhll);
                AE_MULAAAA2Q16X8(A0, A1, y0, y0, x00);
                AE_MULAAAA2Q16X8(A0, A1, y1, y1, x11);
                AE_MULAAAA2Q16X8(A2, A3, y0, y0, x22);
                AE_MULAAAA2Q16X8(A2, A3, y1, y1, x33);
            }
            {
                ae_int64 q;
                ae_int32x2 t;
                ae_int16x4 y0;
                ae_int8x8 x0,x1,x2,x3,x00,x11,x22,x33;
                // convert endianness for 16 bit data because 8-bit data also read by 32-bit 
                // instructions
                AE_L64_IP(q,castxcc(ae_int64,pY),sizeof(ae_int64)); y0=AE_MOVINT16X4_FROMINT64(q);
                AE_L32_IP (t,castxcc(ae_int32,pX0),4); x0=AE_MOVINT8X8_FROMINT32X2(t);
                AE_L32_IP (t,castxcc(ae_int32,pX1),4); x1=AE_MOVINT8X8_FROMINT32X2(t);
                AE_L32_IP (t,castxcc(ae_int32,pX2),4); x2=AE_MOVINT8X8_FROMINT32X2(t);
                AE_L32_IP (t,castxcc(ae_int32,pX3),4); x3=AE_MOVINT8X8_FROMINT32X2(t);
                AE_DSEL8X8(x00,x11,x0,x1,dsel_rephhll);
                AE_DSEL8X8(x22,x33,x2,x3,dsel_rephhll);
                AE_MULAAAA2Q16X8(A0, A1, y0, y0, x00);
                AE_MULAAAA2Q16X8(A2, A3, y0, y0, x22);
            }
            r=AE_ROUND16X4F32SASYM(AE_TRUNCA32X2F64S(A0,A1, sa),AE_TRUNCA32X2F64S(A2,A3, sa));
            AE_S16X4_IP(r,castxcc(ae_int16x4,z),sizeof(ae_int16x4));
            pX0=(const ae_int8x8*)( (uintptr_t)pX0 + 3*N);
            pX1=(const ae_int8x8*)( (uintptr_t)pX1 + 3*N);
            pX2=(const ae_int8x8*)( (uintptr_t)pX2 + 3*N);
            pX3=(const ae_int8x8*)( (uintptr_t)pX3 + 3*N);
        }
        return;
    }
    if (N&8)
    {
        // for N a multiple of 8 but not a multiple of 16
        // in that case (pX0,pX1), (pX2,pX3) are belonging to the different banks
        for (m=0; m<M;m+=4)
        {
            ae_int16x4 r;
            ae_int64 A0,A1,A2,A3;
            AE_MOVDX2(A0,A1,0,0); AE_MOVDX2(A2,A3,0,0);
            pY=(const ae_int16x4*)y;
            __Pragma("ymemory(pX1)");
            __Pragma("ymemory(pX3)");
            __Pragma("loop_count min=1")
            for (n=0; n<(N>>3); n++)
            {
                ae_int16x4 y0,y1;
                ae_int8x8 x0,x1,x2,x3,x00,x11,x22,x33;
                AE_L16X4X2_IP(y0,y1,castxcc(ae_int16x8,pY),sizeof(ae_int16x8));
                AE_L8X8_IP (x0,pX0,sizeof(ae_int8x8));
                AE_L8X8_IP (x1,pX1,sizeof(ae_int8x8));
                AE_L8X8_IP (x2,pX2,sizeof(ae_int8x8));
                AE_L8X8_IP (x3,pX3,sizeof(ae_int8x8));
                AE_DSEL8X8(x00,x11,x0,x1,dsel_rephhll);
                AE_DSEL8X8(x22,x33,x2,x3,dsel_rephhll);
                AE_MULAAAA2Q16X8(A0, A1, y0, y0, x00);
                AE_MULAAAA2Q16X8(A0, A1, y1, y1, x11);
                AE_MULAAAA2Q16X8(A2, A3, y0, y0, x22);
                AE_MULAAAA2Q16X8(A2, A3, y1, y1, x33);
            }
            r=AE_ROUND16X4F32SASYM(AE_TRUNCA32X2F64S(A0,A1, sa),AE_TRUNCA32X2F64S(A2,A3, sa));
            AE_S16X4_IP(r,castxcc(ae_int16x4,z),sizeof(ae_int16x4));
            pX0=(const ae_int8x8*)( (uintptr_t)pX0 + 3*N);
            pX1=(const ae_int8x8*)( (uintptr_t)pX1 + 3*N);
            pX2=(const ae_int8x8*)( (uintptr_t)pX2 + 3*N);
            pX3=(const ae_int8x8*)( (uintptr_t)pX3 + 3*N);
        }
    }
    else
    {
        // for N a multiple of 16
        NASSERT(N%16==0);
        for (m=0; m<(M>>2);m++)
        {
            ae_int16x4 r;
            ae_int64 A0,A1,A2,A3;
            AE_MOVDX2(A0,A1,0,0); AE_MOVDX2(A2,A3,0,0);
            pY=(const ae_int16x4*)y;
            __Pragma("loop_count min=1")
            for (n=0; n<(N>>4); n++)
            {
                ae_int16x4 y0,y1,y2,y3;
                ae_int8x8 x0,x1,x2,x3,x00,x11,x22,x33;
                ae_int8x8 x4,x5,x6,x7,x44,x55,x66,x77;
                AE_L16X4X2_IP(y0,y1,castxcc(ae_int16x8,pY),sizeof(ae_int16x8));
                AE_L16X4X2_IP(y2,y3,castxcc(ae_int16x8,pY),sizeof(ae_int16x8));
                AE_L8X8X2_IP (x0,x4,castxcc(ae_int8x16,pX0),sizeof(ae_int8x16));
                AE_L8X8X2_IP (x1,x5,castxcc(ae_int8x16,pX1),sizeof(ae_int8x16));
                AE_L8X8X2_IP (x2,x6,castxcc(ae_int8x16,pX2),sizeof(ae_int8x16));
                AE_L8X8X2_IP (x3,x7,castxcc(ae_int8x16,pX3),sizeof(ae_int8x16));
                AE_DSEL8X8(x00,x11,x0,x1,dsel_rephhll);
                AE_DSEL8X8(x22,x33,x2,x3,dsel_rephhll);
                AE_DSEL8X8(x44,x55,x4,x5,dsel_rephhll);
                AE_DSEL8X8(x66,x77,x6,x7,dsel_rephhll);
                AE_MULAAAA2Q16X8(A0, A1, y0, y0, x00);
                AE_MULAAAA2Q16X8(A0, A1, y1, y1, x11);
                AE_MULAAAA2Q16X8(A2, A3, y0, y0, x22);
                AE_MULAAAA2Q16X8(A2, A3, y1, y1, x33);
                AE_MULAAAA2Q16X8(A0, A1, y2, y2, x44);
                AE_MULAAAA2Q16X8(A0, A1, y3, y3, x55);
                AE_MULAAAA2Q16X8(A2, A3, y2, y2, x66);
                AE_MULAAAA2Q16X8(A2, A3, y3, y3, x77);
            }
            r=AE_ROUND16X4F32SASYM(AE_TRUNCA32X2F64S(A0,A1, sa),AE_TRUNCA32X2F64S(A2,A3, sa));
            AE_S16X4_IP(r,castxcc(ae_int16x4,z),sizeof(ae_int16x4));
            pX0=(const ae_int8x8*)( (uintptr_t)pX0 + 3*N);
            pX1=(const ae_int8x8*)( (uintptr_t)pX1 + 3*N);
            pX2=(const ae_int8x8*)( (uintptr_t)pX2 + 3*N);
            pX3=(const ae_int8x8*)( (uintptr_t)pX3 + 3*N);
        }
    }
}
#else
void mtx_vecmpy8x16_fast( int16_t* restrict z,
                 const int8_t* restrict x,
                 const int16_t* restrict y,
                 int M, int N, int lsh)
{
    ae_valignx2 alX1, alX2, alX3;
    int16_t * restrict pZ0;
    const int8_t * restrict pX0;
    const int8_t * restrict pX1;
    const int8_t * restrict pX2;
    const int8_t * restrict pX3;
    const int16_t * restrict pY0;
    ae_int32x2 Z0, Z1;
    ae_int8x8 x00, x10, x20, x30, x01, x11, x21, x31;
    ae_int16x4 y0, y1, y2, y3, z0;
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

    WUR_AE_SAR(9+lsh);
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

                /* load y vector, 16 values, 16-bit */
                AE_L16X4X2_IP(y0, y1, castxcc(ae_int16x8,pY0), 8*sizeof(int16_t));
                AE_L16X4X2_IP(y2, y3, castxcc(ae_int16x8,pY0), 8*sizeof(int16_t));

                /* make multiply, using 4-way 32-bit output octal mac, 8X8 bit */
                AE_MULA8Q8X16(Z0, Z1, x00, x10, x20, x30, y0, y1);
                AE_MULA8Q8X16(Z0, Z1, x01, x11, x21, x31, y2, y3);
            }
            {
                int offs1, offs2;
                ae_valignx2 alX0;
                alX0 = AE_LA128_PP(pX0);
                /* load x matrix, 4x16 values, 8-bit */
                AE_LAV8X8X2_XP(x00, x01, alX0, castxcc(ae_int8x16,pX0), N&15);
                AE_LAV8X8X2_XP(x10, x11, alX1, castxcc(ae_int8x16,pX1), N&15);
                AE_LAV8X8X2_XP(x20, x21, alX2, castxcc(ae_int8x16,pX2), N&15);
                AE_LAV8X8X2_XP(x30, x31, alX3, castxcc(ae_int8x16,pX3), N&15);

                /* Load last 4...12 values in a special manner to avoid outbound memory accesses */
                offs1 = (N&15)>4 ? sizeof(ae_int16x4) : 0;
                offs2 = (N&15)>8 ? sizeof(ae_int16x4)*2 : 0;
                y0 = AE_L16X4_I((ae_int16x4 *)pY0, 0);
                y1 = AE_L16X4_X((ae_int16x4 *)pY0, offs1);
                y2 = AE_L16X4_X((ae_int16x4 *)pY0, offs2);
                y3 = AE_XOR16(y3, y3);
                /* make multiply, using 4-way 32-bit output octal mac, 8X8 bit */
                AE_MULA8Q8X16(Z0, Z1, x00, x10, x20, x30, y0, y1);
                AE_MULA8Q8X16(Z0, Z1, x01, x11, x21, x31, y2, y3);
            }
            /* Q15 + lsh <- Q22 - 7 + lsh w/ rounding and saturation */
            Z0 = AE_SLAS32S(Z0);
            Z1 = AE_SLAS32S(Z1);
            z0 = AE_ROUND16X4F32SASYM(Z0, Z1);
            AE_S16X4_IP(z0, castxcc(ae_int16x4,pZ0), sizeof(ae_int16x4));
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
                AE_L16X4X2_IP(y0, y1, castxcc(ae_int16x8,pY0), 8*sizeof(int16_t));
                AE_L16X4X2_IP(y2, y3, castxcc(ae_int16x8,pY0), 8*sizeof(int16_t));

                /* make multiply, using 4-way 32-bit output octal mac, 8X8 bit */
                AE_MULA8Q8X16(Z0, Z1, x00, x10, x20, x30, y0, y1);
                AE_MULA8Q8X16(Z0, Z1, x01, x11, x21, x31, y2, y3);
            }
            /* Q15 + lsh <- Q22 - 7 + lsh w/ rounding and saturation */
            Z0 = AE_SLAS32S(Z0);
            Z1 = AE_SLAS32S(Z1);
            z0 = AE_ROUND16X4F32SASYM(Z0, Z1);
            AE_S16X4_IP(z0, castxcc(ae_int16x4,pZ0), sizeof(ae_int16x4));
        }
    }
}
#endif
