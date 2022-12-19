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
  Matrix Multiply
  These functions compute the expression z = 2^lsh * x * y for the matrices 
  x and y. The columnar dimension of x must match the row dimension of y. 
  The resulting matrix has the same number of rows as x and the same number 
  of columns as y.
  Transposing API allows to interpret input yt as transposed matrix y.

  NOTE: lsh factor is not relevant for floating point routines.

  Functions require scratch memory for storing intermediate data. This 
  scratch memory area should be aligned on 16 byte boundary and its size is 
  calculated by dedicated scratch allocation functions.

  Two versions of functions available: regular version (mtx_mpy[t]32x32, 
  mtx_mpy[t]16x16, mtx_mpy[t]8x16, mtx_mpy[t]8x8, mtx[t]_mpyf) with 
  arbitrary arguments and faster version (mtx_mpy[t]32x32_fast, 
  mtx_mpy[t]16x16_fast, mtx_mpy[t]8x16_fast, mtx_mpy[t]8x8_fast, 
  mtx_mpy[t]f_fast, cntx_mpyt32x32_fast, cntx_mpytf_fast) that apply 
  some restrictions

  Precision:
  32x32 32-bit inputs, 32-bit output
  16x16 16-bit inputs, 16-bit output
  8x8   8-bit inputs, 8-bit output
  8x16  8/16-bit inputs, 16-bit output
  f     floating point

  Input:
  x[M*N]      input matrix x, Q7, Q15, Q31 or floating point
  y[N*P]      input matrix y, Q7, Q15, Q31 or floating point
  yt[P*N]     transposed input matrix y. Q31,Q15, Q7 floating point. (for 
              transposing API only)
  M           number of rows in matrix x and z
  N           number of columns in matrix x and number of rows in matrix y
  P           number of columns in matrices y and z
  lsh         left shift applied to the result (applied to the fixed-
              point functions only) 
  Output:
  z[M*P]      output matrix z, Q7, Q15, Q31 or floating point 
  Scratch:
  pScr        size in bytes defined by corresponding scratch allocation 
              functions

  Restrictions:
  For regular routines mpy[t]32x32, mtx_mpy[t]16x16, mtx_mpy[t]8x8, 
  mtx_mpy[t]8x16, mtx_mpy[t]f):
  pScr    aligned on 16-byte boundary
  x,y,z   should not overlap

  For faster routines (mtx_mpy[t]32x32_fast, mtx_mpy[t]16x16_fast, 
  mtx_mpy[t]8x8_fast, mtx_mpy[t]8x16_fast, 
  mtx_mpy[t]f_fast):
  x,y,z       should not overlap
  x,y,z,pScr  aligned on 16-byte boundary
  M,N,P       multiplies of 4 for mtx_mpy[t]32x32_fast, mtx_mpy[t]16x16_fast, 
              mtx_mpy[t]8x8_fast, mtx_mpy[t]8x16_fast, mtx_mpy[t]f_fast
              multiplies of 32 for cntx_mpyt32x32_fast, cntx_mpytf_fast
  lsh         should be in range:
              -31...31 for mtx_mpy32x32, mtx_mpy32x32_fast, cntx_mpyt32x32_fast, 
                       cntx_mpytf_fast
              -15...15 for mtx_mpy16x16, mtx_mpy16x16_fast, mtx_mpy[t]8x8, 
                       mtx_mpy[t]8x8_fast, mtx_mpy[t]8x16, 
                       mtx_mpy[t]8x16_fast 

-------------------------------------------------------------------------*/
#if !defined(AE_MULA8Q8X8) || !defined(AE_MULA4O8X16)
void mtx_mpyt8x8_fast (  void* pScr,
                     int8_t* restrict z,
               const int8_t* restrict x,
               const int8_t* restrict y,
               int M, int N, int P, int lsh )
#if 0
{
    int8_t* restrict pZ;
    const int8_t* restrict pX0;
    const int8_t* restrict pX1;
    const int8_t* restrict pY0;
    const int8_t* restrict pY1;
    int sa=49+lsh;
    int m,n,p;
    NASSERT(lsh >= -15 && lsh <= 15);
    NASSERT_ALIGN(pScr,HIFI_SIMD_WIDTH);
    (void)pScr;
    NASSERT_ALIGN(pScr,HIFI_SIMD_WIDTH);
    NASSERT_ALIGN(z,HIFI_SIMD_WIDTH);
    NASSERT_ALIGN(x,HIFI_SIMD_WIDTH);
    NASSERT_ALIGN(y,HIFI_SIMD_WIDTH);
    NASSERT(M%4==0);
    NASSERT(N%4==0);
    NASSERT(P%4==0);
    (void)pScr;
    if (P<=0 || M<=0) return;
    if (N<=0)    /* exceptional situation */
    {
        for (m=0; m<((M*P)>>4); m++) AE_S16X4X2_IP(AE_ZERO16(),AE_ZERO16(),castxcc(ae_int16x8,z),sizeof(ae_int16x8));
        return ;
    }
    pZ=z;
    pX0=x;
    pX1=x+2*N;
    __Pragma("loop_count min=1")
    for (m=0; m<(M>>2);m++)
    {
        pY0=y;
        pY1=(const int8_t*)XT_ADDX2(N,(uintptr_t)pY0);
        __Pragma("loop_count min=1")
        for (p=0; p<(P>>2); p++)
        {
            ae_int32x2 res0,res1,res2,res3;
            ae_int8x8  r0,r1,r2,r3;
            ae_int64 a00,a01,a02,a03,a10,a11,a12,a13;
            ae_int64 a20,a21,a22,a23,a30,a31,a32,a33;
            AE_MOVDX2(a00,a01,0,0); AE_MOVDX2(a02,a03,0,0); AE_MOVDX2(a10,a11,0,0); AE_MOVDX2(a12,a13,0,0); 
            AE_MOVDX2(a20,a21,0,0); AE_MOVDX2(a22,a23,0,0); AE_MOVDX2(a30,a31,0,0); AE_MOVDX2(a32,a33,0,0); 
            __Pragma("loop_count min=1")
            for (n=0; n<(N>>2); n++)
            {
                 ae_int32x2 x0,x1,x2,x3,y0,y1,y2,y3;
                 ae_int8x8 y01,y23;
                 x1=AE_L32_X((const ae_int32*)pX0,N);
                 AE_L32_IP(x0,castxcc(ae_int32,pX0),4);
                 x3=AE_L32_X((const ae_int32*)pX1,N);
                 AE_L32_IP(x2,castxcc(ae_int32,pX1),4);
                 y1=AE_L32_X((const ae_int32*)pY0,N);
                 AE_L32_IP(y0,castxcc(ae_int32,pY0),4);
                 y3=AE_L32_X((const ae_int32*)pY1,N);
                 AE_L32_IP(y2,castxcc(ae_int32,pY1),4);

                 y01=AE_MOVINT8X8_FROMINT32X2(AE_SEL32_HH(y0,y1));
                 y23=AE_MOVINT8X8_FROMINT32X2(AE_SEL32_HH(y2,y3));
                 AE_MULAAAA2Q8(a00,a01,AE_MOVINT8X8_FROMINT32X2(x0),y01);
                 AE_MULAAAA2Q8(a02,a03,AE_MOVINT8X8_FROMINT32X2(x0),y23);
                 AE_MULAAAA2Q8(a10,a11,AE_MOVINT8X8_FROMINT32X2(x1),y01);
                 AE_MULAAAA2Q8(a12,a13,AE_MOVINT8X8_FROMINT32X2(x1),y23);
                 AE_MULAAAA2Q8(a20,a21,AE_MOVINT8X8_FROMINT32X2(x2),y01);
                 AE_MULAAAA2Q8(a22,a23,AE_MOVINT8X8_FROMINT32X2(x2),y23);
                 AE_MULAAAA2Q8(a30,a31,AE_MOVINT8X8_FROMINT32X2(x3),y01);
                 AE_MULAAAA2Q8(a32,a33,AE_MOVINT8X8_FROMINT32X2(x3),y23);
            }
            r0=AE_ROUND8X4F32SASYM_L(AE_TRUNCA32X2F64S(a03, a02, sa),AE_TRUNCA32X2F64S(a01, a00, sa));
            r1=AE_ROUND8X4F32SASYM_L(AE_TRUNCA32X2F64S(a13, a12, sa),AE_TRUNCA32X2F64S(a11, a10, sa));
            r2=AE_ROUND8X4F32SASYM_L(AE_TRUNCA32X2F64S(a23, a22, sa),AE_TRUNCA32X2F64S(a21, a20, sa));
            r3=AE_ROUND8X4F32SASYM_L(AE_TRUNCA32X2F64S(a33, a32, sa),AE_TRUNCA32X2F64S(a31, a30, sa));
            res0=AE_MOVINT32X2_FROMINT8X8(r0);
            res1=AE_MOVINT32X2_FROMINT8X8(r1);
            res2=AE_MOVINT32X2_FROMINT8X8(r2);
            res3=AE_MOVINT32X2_FROMINT8X8(r3);
            AE_S32_L_X (res3,(ae_int32*)pZ       ,3*P);
            AE_S32_L_X (res2,(ae_int32*)pZ       ,2*P);
            AE_S32_L_X (res1,(ae_int32*)pZ       ,1*P);
            AE_S32_L_IP(res0,castxcc(ae_int32,pZ),4);
            pX0-=N;
            pX1-=N;
            pY0+=3*N;
            pY1+=3*N;
        }
        pZ+=3*P;
        pX0=(const int8_t*)XT_ADDX4(N,(uintptr_t)pX0);
        pX1=(const int8_t*)XT_ADDX4(N,(uintptr_t)pX1);
    }
}
#else
{
    static const uint64_t dsel_rephhll_tbl=0x40516273c8d9eafbULL;
    ae_int8x8 dsel_rephhll;
    ae_valign aY1,aY3;
    ae_valign aX1,aX3;
    int8_t* restrict pZ;
    const int8_t* restrict pX0;
    const int8_t* restrict pX1;
    const int8_t* restrict pX2;
    const int8_t* restrict pX3;
    const int8_t* restrict pY0;
    const int8_t* restrict pY1;
    const int8_t* restrict pY2;
    const int8_t* restrict pY3;
    int sa=49+lsh;
    int m,n,p;
    NASSERT(lsh >= -15 && lsh <= 15);
    NASSERT_ALIGN(pScr,HIFI_SIMD_WIDTH);
    (void)pScr;
    NASSERT_ALIGN(pScr,HIFI_SIMD_WIDTH);
    NASSERT_ALIGN(z,HIFI_SIMD_WIDTH);
    NASSERT_ALIGN(x,HIFI_SIMD_WIDTH);
    NASSERT_ALIGN(y,HIFI_SIMD_WIDTH);
    NASSERT(M%4==0);
    NASSERT(N%4==0);
    NASSERT(P%4==0);
    (void)pScr;
    if (P<=0 || M<=0) return;
    if (N<=0)    /* exceptional situation */
    {
        for (m=0; m<((M*P)>>4); m++) AE_S16X4X2_IP(AE_ZERO16(),AE_ZERO16(),castxcc(ae_int16x8,z),sizeof(ae_int16x8));
        return ;
    }
    dsel_rephhll=AE_L8X8_I((const ae_int8x8*)&dsel_rephhll_tbl,0);
    pZ=z;
    pX0=x;
    pX1=x+1*N;
    pX2=x+2*N;
    pX3=x+3*N;
    if (N&4)
    {
        /*
            this variant for N not a multiple of 8 but multiple of 4. 
            In that case each even row will be unaligned
        */
        __Pragma("loop_count min=1")
        for (m=0; m<(M>>2);m++)
        {
            pY0=y;
            pY1=y+N;
            pY2=(const int8_t*)XT_ADDX2(N,(uintptr_t)pY0);
            pY3=(const int8_t*)XT_ADDX2(N,(uintptr_t)pY1);
            __Pragma("loop_count min=1")
            for (p=0; p<(P>>2); p++)
            {
                ae_int8x8 x0,x1,x2,x3,x4,x5,x6,x7;
                ae_int32x2 y0,y1,y2,y3;
                ae_int32x2 x04,x15,x26,x37;
                ae_int8x8 y01,y23,y45,y67;
                ae_int32x2 res0,res1,res2,res3;
                ae_int8x8  r0,r1,r2,r3;
                ae_int64 a00,a01,a02,a03,a10,a11,a12,a13;
                ae_int64 a20,a21,a22,a23,a30,a31,a32,a33;
                AE_MOVDX2(a00,a01,0,0); AE_MOVDX2(a02,a03,0,0); AE_MOVDX2(a10,a11,0,0); AE_MOVDX2(a12,a13,0,0); 
                AE_MOVDX2(a20,a21,0,0); AE_MOVDX2(a22,a23,0,0); AE_MOVDX2(a30,a31,0,0); AE_MOVDX2(a32,a33,0,0); 
                NASSERT_ALIGN(pX0,8);
                NASSERT_ALIGN(pX2,8);
                NASSERT_ALIGN(pY0,8);
                NASSERT_ALIGN(pY2,8);
                aY1=AE_LA64_PP(pY1);
                aY3=AE_LA64_PP(pY3);
                aX1=AE_LA64_PP(pX1);
                aX3=AE_LA64_PP(pX3);
                for (n=0; n<(N>>3); n++)
                {
                    AE_L32X2_IP (x04,    castxcc(ae_int32x2,pX0),8);
                    AE_L32X2_IP (x26,    castxcc(ae_int32x2,pX2),8);
                    AE_LA32X2_IP(x15,aX1,castxcc(ae_int32x2,pX1)  );
                    AE_LA32X2_IP(x37,aX3,castxcc(ae_int32x2,pX3)  );
                    AE_L32X2_IP (y0,     castxcc(ae_int32x2,pY0),8);
                    AE_LA32X2_IP(y1,aY1, castxcc(ae_int32x2,pY1)  );
                    AE_L32X2_IP (y2,     castxcc(ae_int32x2,pY2),8);
                    AE_LA32X2_IP(y3,aY3, castxcc(ae_int32x2,pY3)  );

                    AE_DSEL8X8(x0,x4,AE_MOVINT8X8_FROMINT32X2(x04),AE_MOVINT8X8_FROMINT32X2(x04),dsel_rephhll);
                    AE_DSEL8X8(x2,x6,AE_MOVINT8X8_FROMINT32X2(x26),AE_MOVINT8X8_FROMINT32X2(x26),dsel_rephhll);
                    AE_DSEL8X8(x1,x5,AE_MOVINT8X8_FROMINT32X2(x15),AE_MOVINT8X8_FROMINT32X2(x15),dsel_rephhll);
                    AE_DSEL8X8(x3,x7,AE_MOVINT8X8_FROMINT32X2(x37),AE_MOVINT8X8_FROMINT32X2(x37),dsel_rephhll);
                    AE_DSEL8X8(y01,y45,AE_MOVINT8X8_FROMINT32X2(y0),AE_MOVINT8X8_FROMINT32X2(y1),dsel_rephhll);
                    AE_DSEL8X8(y23,y67,AE_MOVINT8X8_FROMINT32X2(y2),AE_MOVINT8X8_FROMINT32X2(y3),dsel_rephhll);

                    AE_MULAAAA2Q8(a00,a01,x0,y01);
                    AE_MULAAAA2Q8(a02,a03,x0,y23);
                    AE_MULAAAA2Q8(a10,a11,x1,y01);
                    AE_MULAAAA2Q8(a12,a13,x1,y23);
                    AE_MULAAAA2Q8(a20,a21,x2,y01);
                    AE_MULAAAA2Q8(a22,a23,x2,y23);
                    AE_MULAAAA2Q8(a30,a31,x3,y01);
                    AE_MULAAAA2Q8(a32,a33,x3,y23);
                    AE_MULAAAA2Q8(a00,a01,x4,y45);
                    AE_MULAAAA2Q8(a02,a03,x4,y67);
                    AE_MULAAAA2Q8(a10,a11,x5,y45);
                    AE_MULAAAA2Q8(a12,a13,x5,y67);
                    AE_MULAAAA2Q8(a20,a21,x6,y45);
                    AE_MULAAAA2Q8(a22,a23,x6,y67);
                    AE_MULAAAA2Q8(a30,a31,x7,y45);
                    AE_MULAAAA2Q8(a32,a33,x7,y67);
                }
                {
                    ae_int32x2 t;
                    AE_L32_IP(t,castxcc(ae_int32,pX0),4); x0=AE_MOVINT8X8_FROMINT32X2( t);
                    AE_L32_IP(t,castxcc(ae_int32,pX1),4); x1=AE_MOVINT8X8_FROMINT32X2( t);
                    AE_L32_IP(t,castxcc(ae_int32,pX2),4); x2=AE_MOVINT8X8_FROMINT32X2( t);
                    AE_L32_IP(t,castxcc(ae_int32,pX3),4); x3=AE_MOVINT8X8_FROMINT32X2( t);
                    AE_L32_IP(y0,castxcc(ae_int32,pY0),4);
                    AE_L32_IP(y1,castxcc(ae_int32,pY1),4);
                    AE_L32_IP(y2,castxcc(ae_int32,pY2),4);
                    AE_L32_IP(y3,castxcc(ae_int32,pY3),4);

                    y01=AE_MOVINT8X8_FROMINT32X2(AE_SEL32_HH(y0,y1));
                    y23=AE_MOVINT8X8_FROMINT32X2(AE_SEL32_HH(y2,y3));
                    AE_MULAAAA2Q8(a00,a01,x0,y01);
                    AE_MULAAAA2Q8(a02,a03,x0,y23);
                    AE_MULAAAA2Q8(a10,a11,x1,y01);
                    AE_MULAAAA2Q8(a12,a13,x1,y23);
                    AE_MULAAAA2Q8(a20,a21,x2,y01);
                    AE_MULAAAA2Q8(a22,a23,x2,y23);
                    AE_MULAAAA2Q8(a30,a31,x3,y01);
                    AE_MULAAAA2Q8(a32,a33,x3,y23);
                }
                r0=AE_ROUND8X4F32SASYM_L(AE_TRUNCA32X2F64S(a03, a02, sa),AE_TRUNCA32X2F64S(a01, a00, sa));
                r1=AE_ROUND8X4F32SASYM_L(AE_TRUNCA32X2F64S(a13, a12, sa),AE_TRUNCA32X2F64S(a11, a10, sa));
                r2=AE_ROUND8X4F32SASYM_L(AE_TRUNCA32X2F64S(a23, a22, sa),AE_TRUNCA32X2F64S(a21, a20, sa));
                r3=AE_ROUND8X4F32SASYM_L(AE_TRUNCA32X2F64S(a33, a32, sa),AE_TRUNCA32X2F64S(a31, a30, sa));
                res0=AE_MOVINT32X2_FROMINT8X8(r0);
                res1=AE_MOVINT32X2_FROMINT8X8(r1);
                res2=AE_MOVINT32X2_FROMINT8X8(r2);
                res3=AE_MOVINT32X2_FROMINT8X8(r3);
                AE_S32_L_X (res3,(ae_int32*)pZ       ,3*P);
                AE_S32_L_X (res2,(ae_int32*)pZ       ,2*P);
                AE_S32_L_X (res1,(ae_int32*)pZ       ,1*P);
                AE_S32_L_IP(res0,castxcc(ae_int32,pZ),4);
                pX0-=N;
                pX1-=N;
                pX2-=N;
                pX3-=N;
                pY0+=3*N;
                pY1+=3*N;
                pY2+=3*N;
                pY3+=3*N;
            }
            pZ+=3*P;
            pX0=(const int8_t*)XT_ADDX4(N,(uintptr_t)pX0);
            pX1=(const int8_t*)XT_ADDX4(N,(uintptr_t)pX1);
            pX2=(const int8_t*)XT_ADDX4(N,(uintptr_t)pX2);
            pX3=(const int8_t*)XT_ADDX4(N,(uintptr_t)pX3);
        }
    }
    else
    {
        /*
            separate code for N a multiple of 8 
        */
        __Pragma("loop_count min=1")
        for (m=0; m<(M>>2);m++)
        {
            pY0=y;
            pY1=y+N;
            pY2=(const int8_t*)XT_ADDX2(N,(uintptr_t)pY0);
            pY3=(const int8_t*)XT_ADDX2(N,(uintptr_t)pY1);
            __Pragma("loop_count min=1")
            for (p=0; p<(P>>2); p++)
            {
                ae_int8x8 x0,x1,x2,x3,x4,x5,x6,x7;
                ae_int32x2 y0,y1,y2,y3;
                ae_int32x2 x04,x15,x26,x37;
                ae_int8x8 y01,y23,y45,y67;
                ae_int32x2 res0,res1,res2,res3;
                ae_int8x8  r0,r1,r2,r3;
                ae_int64 a00,a01,a02,a03,a10,a11,a12,a13;
                ae_int64 a20,a21,a22,a23,a30,a31,a32,a33;
                AE_MOVDX2(a00,a01,0,0); AE_MOVDX2(a02,a03,0,0); AE_MOVDX2(a10,a11,0,0); AE_MOVDX2(a12,a13,0,0); 
                AE_MOVDX2(a20,a21,0,0); AE_MOVDX2(a22,a23,0,0); AE_MOVDX2(a30,a31,0,0); AE_MOVDX2(a32,a33,0,0); 
                NASSERT_ALIGN(pX0,8);
                NASSERT_ALIGN(pX2,8);
                NASSERT_ALIGN(pY0,8);
                NASSERT_ALIGN(pY2,8);
                __Pragma("loop_count min=1")
                for (n=0; n<(N>>3); n++)
                {
                    x15=AE_L32X2_X ((const ae_int32x2*)pX0,N);
                    x37=AE_L32X2_X ((const ae_int32x2*)pX2,N);
                    AE_L32X2_IP (x04,castxcc(ae_int32x2,pX0),8);
                    AE_L32X2_IP (x26,castxcc(ae_int32x2,pX2),8);

                    y1=AE_L32X2_X  ((const ae_int32x2*)pY0,N);
                    y3=AE_L32X2_X  ((const ae_int32x2*)pY2,N);
                    AE_L32X2_IP    (y0, castxcc(ae_int32x2,pY0),8);
                    AE_L32X2_IP    (y2, castxcc(ae_int32x2,pY2),8);

                    AE_DSEL8X8(x0,x4,AE_MOVINT8X8_FROMINT32X2(x04),AE_MOVINT8X8_FROMINT32X2(x04),dsel_rephhll);
                    AE_DSEL8X8(x2,x6,AE_MOVINT8X8_FROMINT32X2(x26),AE_MOVINT8X8_FROMINT32X2(x26),dsel_rephhll);
                    AE_DSEL8X8(x1,x5,AE_MOVINT8X8_FROMINT32X2(x15),AE_MOVINT8X8_FROMINT32X2(x15),dsel_rephhll);
                    AE_DSEL8X8(x3,x7,AE_MOVINT8X8_FROMINT32X2(x37),AE_MOVINT8X8_FROMINT32X2(x37),dsel_rephhll);
                    AE_DSEL8X8(y01,y45,AE_MOVINT8X8_FROMINT32X2(y0),AE_MOVINT8X8_FROMINT32X2(y1),dsel_rephhll);
                    AE_DSEL8X8(y23,y67,AE_MOVINT8X8_FROMINT32X2(y2),AE_MOVINT8X8_FROMINT32X2(y3),dsel_rephhll);

                    AE_MULAAAA2Q8(a00,a01,x0,y01);
                    AE_MULAAAA2Q8(a02,a03,x0,y23);
                    AE_MULAAAA2Q8(a10,a11,x1,y01);
                    AE_MULAAAA2Q8(a12,a13,x1,y23);
                    AE_MULAAAA2Q8(a20,a21,x2,y01);
                    AE_MULAAAA2Q8(a22,a23,x2,y23);
                    AE_MULAAAA2Q8(a30,a31,x3,y01);
                    AE_MULAAAA2Q8(a32,a33,x3,y23);
                    AE_MULAAAA2Q8(a00,a01,x4,y45);
                    AE_MULAAAA2Q8(a02,a03,x4,y67);
                    AE_MULAAAA2Q8(a10,a11,x5,y45);
                    AE_MULAAAA2Q8(a12,a13,x5,y67);
                    AE_MULAAAA2Q8(a20,a21,x6,y45);
                    AE_MULAAAA2Q8(a22,a23,x6,y67);
                    AE_MULAAAA2Q8(a30,a31,x7,y45);
                    AE_MULAAAA2Q8(a32,a33,x7,y67);
                }
                r0=AE_ROUND8X4F32SASYM_L(AE_TRUNCA32X2F64S(a03, a02, sa),AE_TRUNCA32X2F64S(a01, a00, sa));
                r1=AE_ROUND8X4F32SASYM_L(AE_TRUNCA32X2F64S(a13, a12, sa),AE_TRUNCA32X2F64S(a11, a10, sa));
                r2=AE_ROUND8X4F32SASYM_L(AE_TRUNCA32X2F64S(a23, a22, sa),AE_TRUNCA32X2F64S(a21, a20, sa));
                r3=AE_ROUND8X4F32SASYM_L(AE_TRUNCA32X2F64S(a33, a32, sa),AE_TRUNCA32X2F64S(a31, a30, sa));
                res0=AE_MOVINT32X2_FROMINT8X8(r0);
                res1=AE_MOVINT32X2_FROMINT8X8(r1);
                res2=AE_MOVINT32X2_FROMINT8X8(r2);
                res3=AE_MOVINT32X2_FROMINT8X8(r3);
                AE_S32_L_X (res3,(ae_int32*)pZ       ,3*P);
                AE_S32_L_X (res2,(ae_int32*)pZ       ,2*P);
                AE_S32_L_X (res1,(ae_int32*)pZ       ,1*P);
                AE_S32_L_IP(res0,castxcc(ae_int32,pZ),4);
                pX0-=N;
                pX2-=N;
                pY0+=3*N;
                pY2+=3*N;
            }
            pZ+=3*P;
            pX0=(const int8_t*)XT_ADDX4(N,(uintptr_t)pX0);
            pX2=(const int8_t*)XT_ADDX4(N,(uintptr_t)pX2);
        }
    }
}
#endif

#else
/* Variant with using the 'Neural Network' instruction extension */
void mtx_mpyt8x8_fast (  void* pScr,
                     int8_t* restrict z,
               const int8_t* restrict x,
               const int8_t* restrict y,
               int M, int N, int P, int lsh )
{
    ae_valign alX1, alX3, alY1, alY3;
    int8_t * restrict pZ0;
    const int8_t * restrict pX0;
    const int8_t * restrict pX1;
    const int8_t * restrict pX2;
    const int8_t * restrict pX3;
    const int8_t * restrict pY0;
    const int8_t * restrict pY1;
    const int8_t * restrict pY2;
    const int8_t * restrict pY3;
    ae_int32x2 Z00, Z01, Z10, Z11, Z20, Z21, Z30, Z31;
    ae_int16x4 t0, t1, t2, t3;
    ae_int16x4 X0, X1, X2, X3, Y0, Y1, Y2, Y3;
    ae_int8x8 x0, x1, x2, x3;
    ae_int8x8 y0, y1, y2, y3;
    ae_int8x8 z0, z1;
    int m,n,p;

    NASSERT(lsh >= -15 && lsh <= 15);
    NASSERT_ALIGN(pScr,HIFI_SIMD_WIDTH);
    NASSERT_ALIGN(z,HIFI_SIMD_WIDTH);
    NASSERT_ALIGN(x,HIFI_SIMD_WIDTH);
    NASSERT_ALIGN(y,HIFI_SIMD_WIDTH);
    NASSERT(M%4==0);
    NASSERT(N%4==0);
    NASSERT(P%4==0);
    (void)pScr;
    if (P<=0 || M<=0) return;
    if (N<=0)    /* exceptional situation */
    {
        for (m=0; m<((M*P)>>4); m++) AE_S16X4X2_IP(AE_ZERO16(),AE_ZERO16(),castxcc(ae_int16x8,z),sizeof(ae_int16x8));
        return ;
    }

    if (N&4)
    {
        /*
            this variant for N not a multiple of 8 but multiple of 4. 
            In that case each even row will be unaligned
        */
        __Pragma("loop_count min=1")
        for (m=0; m<(M>>2);m++)
        {
            pZ0 = z + m*P*4;
            __Pragma("loop_count min=1")
            for (p=0; p<(P>>2); p++)
            {
                pX0 = x + m*N*4;
                pX1 = pX0 + N;
                pX2 = pX1 + N;
                pX3 = pX2 + N;

                pY0 = y + p*N*4;
                pY1 = pY0 + N;
                pY2 = pY1 + N;
                pY3 = pY2 + N;

                alX1 = AE_LA64_PP(pX1);
                alX3 = AE_LA64_PP(pX3);
                alY1 = AE_LA64_PP(pY1);
                alY3 = AE_LA64_PP(pY3);
                Z00 = Z01 = Z10 = Z11 = Z20 = Z21 = Z30 = Z31 = AE_ZERO32();

                for (n=0; n<(N>>3); n++)
                {
                    /* load x matrix, 4x8 values, 8-bit */
                    AE_L8X8_IP (x0,       castxcc(ae_int8x8,pX0), 8);
                    AE_LA8X8_IP(x1, alX1, castxcc(ae_int8x8,pX1));
                    AE_L8X8_IP (x2,       castxcc(ae_int8x8,pX2), 8);
                    AE_LA8X8_IP(x3, alX3, castxcc(ae_int8x8,pX3));

                    /* load y matrix, 4x8 values, 8-bit */
                    AE_L8X8_IP (y0,       castxcc(ae_int8x8,pY0), 8);
                    AE_LA8X8_IP(y1, alY1, castxcc(ae_int8x8,pY1));
                    AE_L8X8_IP (y2,       castxcc(ae_int8x8,pY2), 8);
                    AE_LA8X8_IP(y3, alY3, castxcc(ae_int8x8,pY3));

                    /* make multiply, using 4-way 32-bit output octal mac, 8X8 bit */
                    AE_MULA8Q8X8(Z00, Z01, y1, y0, y3, y2, x0);
                    AE_MULA8Q8X8(Z10, Z11, y1, y0, y3, y2, x1);
                    AE_MULA8Q8X8(Z20, Z21, y1, y0, y3, y2, x2);
                    AE_MULA8Q8X8(Z30, Z31, y1, y0, y3, y2, x3);
                }
                {
                    /* load x matrix, 4x4 values, cvt 8->16-bit */
                    X0 = AE_L8X4S_I(pX0, 0);
                    X1 = AE_L8X4S_I(pX1, 0);
                    X2 = AE_L8X4S_I(pX2, 0);
                    X3 = AE_L8X4S_I(pX3, 0);

                    /* load y matrix, 4x4 values, 8-bit */
                    Y0 = AE_L8X4S_I(pY0, 0);
                    Y1 = AE_L8X4S_I(pY1, 0);
                    Y2 = AE_L8X4S_I(pY2, 0);
                    Y3 = AE_L8X4S_I(pY3, 0);
                    y0 = AE_SAT8X8X16(Y0, Y0);
                    y1 = AE_SAT8X8X16(Y1, Y1);
                    y2 = AE_SAT8X8X16(Y2, Y2);
                    y3 = AE_SAT8X8X16(Y3, Y3);

                    /* make multiply, using 8-way 32-bit output quad mac, 8X16 bit */
                    AE_MULA4O8X16(Z00, Z01, Z10, Z11, y1, y0, y3, y2, X0, X1);
                    AE_MULA4O8X16(Z20, Z21, Z30, Z31, y1, y0, y3, y2, X2, X3);
                }
                /* Q15 + lsh <- Q14 + 1 + lsh w/ truncation and saturation */
                t0 = AE_TRUNCA16X4F32S(Z01, Z00, 16+1+lsh);
                t1 = AE_TRUNCA16X4F32S(Z11, Z10, 16+1+lsh);
                t2 = AE_TRUNCA16X4F32S(Z21, Z20, 16+1+lsh);
                t3 = AE_TRUNCA16X4F32S(Z31, Z30, 16+1+lsh);
                /* Q7 + lsh <- Q15 + lsh -8 w/ rounding and saturation */
                z0 = AE_ROUND8X8F16SASYM(t0, t1);
                z1 = AE_ROUND8X8F16SASYM(t2, t3);
                AE_S32_H_XP(AE_MOVINT32X2_FROMINT8X8(z0), castxcc(ae_int32,pZ0), P);
                AE_S32_L_XP(AE_MOVINT32X2_FROMINT8X8(z0), castxcc(ae_int32,pZ0), P);
                AE_S32_H_XP(AE_MOVINT32X2_FROMINT8X8(z1), castxcc(ae_int32,pZ0), P);
                AE_S32_L_XP(AE_MOVINT32X2_FROMINT8X8(z1), castxcc(ae_int32,pZ0), -3*P+4);
            }
        }
    }
    else
    {
        __Pragma("loop_count min=1")
        for (m=0; m<(M>>2);m++)
        {
            pZ0 = z + m*P*4;
            __Pragma("loop_count min=1")
            for (p=0; p<(P>>2); p++)
            {
                pX0 = x + m*N*4;
                pY0 = y + p*N*4;

                Z00 = Z01 = Z10 = Z11 = Z20 = Z21 = Z30 = Z31 = AE_ZERO32();

                __Pragma("loop_count min=1");
                for (n=0; n<(N>>3); n++)
                {
                    /* load x matrix, 4x8 values, 8-bit */
                    AE_L8X8_XP(x0, castxcc(ae_int8x8,pX0), N);
                    AE_L8X8_XP(x1, castxcc(ae_int8x8,pX0), N);
                    AE_L8X8_XP(x2, castxcc(ae_int8x8,pX0), N);
                    AE_L8X8_XP(x3, castxcc(ae_int8x8,pX0), -3*N+8);

                    /* load y matrix, 4x8 values, 8-bit */
                    AE_L8X8_XP(y0, castxcc(ae_int8x8,pY0), N);
                    AE_L8X8_XP(y1, castxcc(ae_int8x8,pY0), N);
                    AE_L8X8_XP(y2, castxcc(ae_int8x8,pY0), N);
                    AE_L8X8_XP(y3, castxcc(ae_int8x8,pY0), -3*N+8);

                    /* make multiply, using 4-way 32-bit output octal mac, 8X8 bit */
                    AE_MULA8Q8X8(Z00, Z01, y1, y0, y3, y2, x0);
                    AE_MULA8Q8X8(Z10, Z11, y1, y0, y3, y2, x1);
                    AE_MULA8Q8X8(Z20, Z21, y1, y0, y3, y2, x2);
                    AE_MULA8Q8X8(Z30, Z31, y1, y0, y3, y2, x3);
                }
                /* Q15 + lsh <- Q14 + 1 + lsh w/ truncation and saturation */
                t0 = AE_TRUNCA16X4F32S(Z01, Z00, 16+1+lsh);
                t1 = AE_TRUNCA16X4F32S(Z11, Z10, 16+1+lsh);
                t2 = AE_TRUNCA16X4F32S(Z21, Z20, 16+1+lsh);
                t3 = AE_TRUNCA16X4F32S(Z31, Z30, 16+1+lsh);
                /* Q7 + lsh <- Q15 + lsh -8 w/ rounding and saturation */
                z0 = AE_ROUND8X8F16SASYM(t0, t1);
                z1 = AE_ROUND8X8F16SASYM(t2, t3);
                AE_S32_H_XP(AE_MOVINT32X2_FROMINT8X8(z0), castxcc(ae_int32,pZ0), P);
                AE_S32_L_XP(AE_MOVINT32X2_FROMINT8X8(z0), castxcc(ae_int32,pZ0), P);
                AE_S32_H_XP(AE_MOVINT32X2_FROMINT8X8(z1), castxcc(ae_int32,pZ0), P);
                AE_S32_L_XP(AE_MOVINT32X2_FROMINT8X8(z1), castxcc(ae_int32,pZ0), -3*P+4);
            }
        }
    }
}
#endif

size_t mtx_mpyt8x8_fast_getScratchSize ( int M, int N, int P)
{
    (void)M;(void)N;(void)P;
    return 0;
}
