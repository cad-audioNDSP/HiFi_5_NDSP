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
/* code optimimized for HiFi5 core */

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
#if !( defined(AE_MULA4O8X16) && defined(AE_MULA8Q8X8) )
void mtx_mpy8x8_fast (  void* pScr,
                     int8_t* restrict z,
               const int8_t* restrict x,
               const int8_t* restrict y,
               int M, int N, int P, int lsh )
{
    int m,n,p;
          int8_t * restrict pZ;
    const int8_t * restrict pX0;
    const int8_t * restrict pX1;
    const int8_t * restrict pY0;
    const int8_t * restrict pY1;
    int sa=49+lsh;
    ae_int8x8 dsel0,dsel1;
    static const uint64_t ALIGN(16) dseltbl[]=
    {
        0x40c851d962ea73fbULL,  // interleave1 {0,8,1,9,2,10,3,11,4,12,5,13,6,14,7,15}
        0x4051c8d96273eafbULL   // interleave2 {0,1,8,9,2,3,10,11,4,5,12,13,6,7,14,15}
    };
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
    dsel0=AE_L8X8_I((const ae_int8x8*)dseltbl,0*sizeof(uint64_t));
    dsel1=AE_L8X8_I((const ae_int8x8*)dseltbl,1*sizeof(uint64_t));
#if 0
    if (N==4)
    {
        pX0=x;
        pX1=x+2*4;
        pZ=z;
        __Pragma("loop_count min=1");
        for (m=0; m<(M>>2);m++)
        {
            ae_int8x8 x0,x1,x2,x3;
            ae_int32x2 t;
            pY0=y;
            pY1=(const int8_t*)XT_ADDX2(P,(uintptr_t)pY0);
            t=AE_L32_I((const ae_int32*)pX0,4)    ;x1=AE_MOVINT8X8_FROMINT32X2(t);
            AE_L32_IP(t,castxcc(ae_int32,pX0),4*4);x0=AE_MOVINT8X8_FROMINT32X2(t);
            t=AE_L32_I((const ae_int32*)pX1,4)    ;x3=AE_MOVINT8X8_FROMINT32X2(t);
            AE_L32_IP(t,castxcc(ae_int32,pX1),4*4);x2=AE_MOVINT8X8_FROMINT32X2(t);
            __Pragma("loop_count min=1");
            for (p=0; p<P; p+=4)
            {
                ae_int8x8 y0,y1,y2,y3,z0,z1,z2,z3;
                ae_int8x8 r0,r1,r2,r3;
                ae_int32x2  res0,res1,res2,res3;
                ae_int64 a00,a01,a02,a03,a10,a11,a12,a13;
                ae_int64 a20,a21,a22,a23,a30,a31,a32,a33;
                t=AE_L32_X ((const ae_int32*)pY1 ,1*P);y0=AE_MOVINT8X8_FROMINT32X2(t);
                AE_L32_IP(t,castxcc(ae_int32,pY1),4  );y1=AE_MOVINT8X8_FROMINT32X2(t);
                t=AE_L32_X ((const ae_int32*)pY0 ,1*P);y2=AE_MOVINT8X8_FROMINT32X2(t);
                AE_L32_IP(t,castxcc(ae_int32,pY0),4  );y3=AE_MOVINT8X8_FROMINT32X2(t);
                // transpose y
                AE_DSEL8X8(z0,z1,y0,y1,dsel0);
                AE_DSEL8X8(z2,z3,y2,y3,dsel0);
                AE_DSEL8X8(y0,y1,z0,z2,dsel1);
                AE_DSEL8X8(y2,y3,z1,z3,dsel1);

                AE_MULZAAAA2Q8(a03,a02,x0,y0);
                AE_MULZAAAA2Q8(a13,a12,x1,y0);
                AE_MULZAAAA2Q8(a01,a00,x0,y1);
                AE_MULZAAAA2Q8(a11,a10,x1,y1);

                AE_MULZAAAA2Q8(a23,a22,x2,y0);
                AE_MULZAAAA2Q8(a33,a32,x3,y0);
                AE_MULZAAAA2Q8(a21,a20,x2,y1);
                AE_MULZAAAA2Q8(a31,a30,x3,y1);

                r0=AE_ROUND8X4F32SASYM_L(AE_TRUNCA32X2F64S(a03, a02, sa),AE_TRUNCA32X2F64S(a01, a00, sa));
                r1=AE_ROUND8X4F32SASYM_L(AE_TRUNCA32X2F64S(a13, a12, sa),AE_TRUNCA32X2F64S(a11, a10, sa));
                r2=AE_ROUND8X4F32SASYM_L(AE_TRUNCA32X2F64S(a23, a22, sa),AE_TRUNCA32X2F64S(a21, a20, sa));
                r3=AE_ROUND8X4F32SASYM_L(AE_TRUNCA32X2F64S(a33, a32, sa),AE_TRUNCA32X2F64S(a31, a30, sa));
                res0=AE_MOVINT32X2_FROMINT8X8(r0);
                res1=AE_MOVINT32X2_FROMINT8X8(r1);
                res2=AE_MOVINT32X2_FROMINT8X8(r2);
                res3=AE_MOVINT32X2_FROMINT8X8(r3);
                AE_S32_L_X(res3,(ae_int32*)pZ,3*P);
                AE_S32_L_X(res2,(ae_int32*)pZ,2*P);
                AE_S32_L_X(res1,(ae_int32*)pZ,1*P);
                AE_S32_L_IP(res0,castxcc(ae_int32,pZ),4);
            }
            pZ+=3*P;
        }
    }
    else
    {
        pX0=x;
        pX1=x+2*N;
        pZ=z;
        __Pragma("loop_count min=1");
        for (m=0; m<(M>>2);m++)
        {
            __Pragma("loop_count min=1");
            for (p=0; p<P; p+=4)
            {
                ae_int8x8 r0,r1,r2,r3;
                ae_int32x2  res0,res1,res2,res3;
                ae_int64 a00,a01,a02,a03,a10,a11,a12,a13;
                ae_int64 a20,a21,a22,a23,a30,a31,a32,a33;
                AE_MOVDX2(a00,a01,0,0); AE_MOVDX2(a02,a03,0,0); AE_MOVDX2(a10,a11,0,0); AE_MOVDX2(a12,a13,0,0); 
                AE_MOVDX2(a20,a21,0,0); AE_MOVDX2(a22,a23,0,0); AE_MOVDX2(a30,a31,0,0); AE_MOVDX2(a32,a33,0,0); 
                pY0=&y[p];
                pY1=(const int8_t*)XT_ADDX2(P,(uintptr_t)pY0);
                __Pragma("loop_count min=2");
                for (n=0; n<(N>>2); n++)
                {
                    ae_int8x8 y0,y1,y2,y3,z0,z1,z2,z3;
                    ae_int8x8 x0,x1,x2,x3;
                    ae_int32x2 t;
                    t=AE_L32_X ((const ae_int32*)pY1 ,1*P);y0=AE_MOVINT8X8_FROMINT32X2(t);
                    AE_L32_XP(t,castxcc(ae_int32,pY1),4*P);y1=AE_MOVINT8X8_FROMINT32X2(t);
                    t=AE_L32_X ((const ae_int32*)pY0 ,1*P);y2=AE_MOVINT8X8_FROMINT32X2(t);
                    AE_L32_XP(t,castxcc(ae_int32,pY0),4*P);y3=AE_MOVINT8X8_FROMINT32X2(t);
                    t=AE_L32_X((const ae_int32*)pX0,N)  ;x1=AE_MOVINT8X8_FROMINT32X2(t);
                    AE_L32_IP(t,castxcc(ae_int32,pX0),4);x0=AE_MOVINT8X8_FROMINT32X2(t);
                    t=AE_L32_X((const ae_int32*)pX1,N)  ;x3=AE_MOVINT8X8_FROMINT32X2(t);
                    AE_L32_IP(t,castxcc(ae_int32,pX1),4);x2=AE_MOVINT8X8_FROMINT32X2(t);
                    // transpose y
                    AE_DSEL8X8(z0,z1,y0,y1,dsel0);
                    AE_DSEL8X8(z2,z3,y2,y3,dsel0);
                    AE_DSEL8X8(y0,y1,z0,z2,dsel1);
                    AE_DSEL8X8(y2,y3,z1,z3,dsel1);

                    AE_MULAAAA2Q8(a03,a02,x0,y0);
                    AE_MULAAAA2Q8(a13,a12,x1,y0);
                    AE_MULAAAA2Q8(a01,a00,x0,y1);
                    AE_MULAAAA2Q8(a11,a10,x1,y1);

                    AE_MULAAAA2Q8(a23,a22,x2,y0);
                    AE_MULAAAA2Q8(a33,a32,x3,y0);
                    AE_MULAAAA2Q8(a21,a20,x2,y1);
                    AE_MULAAAA2Q8(a31,a30,x3,y1);
                }
                r0=AE_ROUND8X4F32SASYM_L(AE_TRUNCA32X2F64S(a03, a02, sa),AE_TRUNCA32X2F64S(a01, a00, sa));
                r1=AE_ROUND8X4F32SASYM_L(AE_TRUNCA32X2F64S(a13, a12, sa),AE_TRUNCA32X2F64S(a11, a10, sa));
                r2=AE_ROUND8X4F32SASYM_L(AE_TRUNCA32X2F64S(a23, a22, sa),AE_TRUNCA32X2F64S(a21, a20, sa));
                r3=AE_ROUND8X4F32SASYM_L(AE_TRUNCA32X2F64S(a33, a32, sa),AE_TRUNCA32X2F64S(a31, a30, sa));
                res0=AE_MOVINT32X2_FROMINT8X8(r0);
                res1=AE_MOVINT32X2_FROMINT8X8(r1);
                res2=AE_MOVINT32X2_FROMINT8X8(r2);
                res3=AE_MOVINT32X2_FROMINT8X8(r3);
                AE_S32_L_X(res3,(ae_int32*)pZ,3*P);
                AE_S32_L_X(res2,(ae_int32*)pZ,2*P);
                AE_S32_L_X(res1,(ae_int32*)pZ,1*P);
                AE_S32_L_IP(res0,castxcc(ae_int32,pZ),4);
                pX0-=N;
                pX1-=N;
            }
            pX0=(const int8_t*)XT_ADDX4(N,(uintptr_t)pX0);
            pX1=(const int8_t*)XT_ADDX4(N,(uintptr_t)pX1);
            pZ+=3*P;
        }
    }
#else
    {
        pX0=x;
        pX1=x+2*N;
        pZ=z;
        __Pragma("loop_count min=1");
        for (m=0; m<(M>>2);m++)
        {
            __Pragma("loop_count min=1");
            for (p=0; p<P; p+=4)
            {
                ae_int8x8 r0,r1,r2,r3;
                ae_int32x2  res0,res1,res2,res3;
                ae_int64 a00,a01,a02,a03,a10,a11,a12,a13;
                ae_int64 a20,a21,a22,a23,a30,a31,a32,a33;
                AE_MOVDX2(a00,a01,0,0); AE_MOVDX2(a02,a03,0,0); AE_MOVDX2(a10,a11,0,0); AE_MOVDX2(a12,a13,0,0); 
                AE_MOVDX2(a20,a21,0,0); AE_MOVDX2(a22,a23,0,0); AE_MOVDX2(a30,a31,0,0); AE_MOVDX2(a32,a33,0,0); 
                pY0=&y[p];
                pY1=(const int8_t*)XT_ADDX2(P,(uintptr_t)pY0);
                __Pragma("loop_count min=1");
                for (n=0; n<(N>>2); n++)
                {
                    ae_int8x8 y0,y1,y2,y3,z0,z1,z2,z3;
                    ae_int8x8 x0,x1,x2,x3;
                    ae_int32x2 t;
                    t=AE_L32_X ((const ae_int32*)pY1 ,1*P);y0=AE_MOVINT8X8_FROMINT32X2(t);
                    AE_L32_XP(t,castxcc(ae_int32,pY1),4*P);y1=AE_MOVINT8X8_FROMINT32X2(t);
                    t=AE_L32_X ((const ae_int32*)pY0 ,1*P);y2=AE_MOVINT8X8_FROMINT32X2(t);
                    AE_L32_XP(t,castxcc(ae_int32,pY0),4*P);y3=AE_MOVINT8X8_FROMINT32X2(t);
                    t=AE_L32_X((const ae_int32*)pX0,N)  ;x1=AE_MOVINT8X8_FROMINT32X2(t);
                    AE_L32_IP(t,castxcc(ae_int32,pX0),4);x0=AE_MOVINT8X8_FROMINT32X2(t);
                    t=AE_L32_X((const ae_int32*)pX1,N)  ;x3=AE_MOVINT8X8_FROMINT32X2(t);
                    AE_L32_IP(t,castxcc(ae_int32,pX1),4);x2=AE_MOVINT8X8_FROMINT32X2(t);
                    // transpose y
                    AE_DSEL8X8(z0,z1,y0,y1,dsel0);
                    AE_DSEL8X8(z2,z3,y2,y3,dsel0);
                    AE_DSEL8X8(y0,y1,z0,z2,dsel1);

                    AE_MULAAAA2Q8(a03,a02,x0,y0);
                    AE_MULAAAA2Q8(a13,a12,x1,y0);
                    AE_MULAAAA2Q8(a01,a00,x0,y1);
                    AE_MULAAAA2Q8(a11,a10,x1,y1);

                    AE_MULAAAA2Q8(a23,a22,x2,y0);
                    AE_MULAAAA2Q8(a33,a32,x3,y0);
                    AE_MULAAAA2Q8(a21,a20,x2,y1);
                    AE_MULAAAA2Q8(a31,a30,x3,y1);
                }
                r0=AE_ROUND8X4F32SASYM_L(AE_TRUNCA32X2F64S(a03, a02, sa),AE_TRUNCA32X2F64S(a01, a00, sa));
                r1=AE_ROUND8X4F32SASYM_L(AE_TRUNCA32X2F64S(a13, a12, sa),AE_TRUNCA32X2F64S(a11, a10, sa));
                r2=AE_ROUND8X4F32SASYM_L(AE_TRUNCA32X2F64S(a23, a22, sa),AE_TRUNCA32X2F64S(a21, a20, sa));
                r3=AE_ROUND8X4F32SASYM_L(AE_TRUNCA32X2F64S(a33, a32, sa),AE_TRUNCA32X2F64S(a31, a30, sa));
                res0=AE_MOVINT32X2_FROMINT8X8(r0);
                res1=AE_MOVINT32X2_FROMINT8X8(r1);
                res2=AE_MOVINT32X2_FROMINT8X8(r2);
                res3=AE_MOVINT32X2_FROMINT8X8(r3);
                AE_S32_L_X(res3,(ae_int32*)pZ,3*P);
                AE_S32_L_X(res2,(ae_int32*)pZ,2*P);
                AE_S32_L_X(res1,(ae_int32*)pZ,1*P);
                AE_S32_L_IP(res0,castxcc(ae_int32,pZ),4);
                pX0-=N;
                pX1-=N;
            }
            pX0=(const int8_t*)XT_ADDX4(N,(uintptr_t)pX0);
            pX1=(const int8_t*)XT_ADDX4(N,(uintptr_t)pX1);
            pZ+=3*P;
        }
    }
#endif
}
#else
/* Variant with using the 'Neural Network' instruction extension */
void mtx_mpy8x8_fast (  void* pScr,
                     int8_t* restrict z,
               const int8_t* restrict x,
               const int8_t* restrict y,
               int M, int N, int P, int lsh )
{
    int m,n,p;
          int8_t * restrict pZ0;
          int8_t * restrict pZ1;
          int8_t * restrict pZ2;
          int8_t * restrict pZ3;
    const int8_t * restrict pX0;
    const int8_t * restrict pX1;
    const int8_t * restrict pX2;
    const int8_t * restrict pX3;
    const int8_t * restrict pY0;
    ae_int32x2 Z00, Z01, Z02, Z03, Z10, Z11, Z12, Z13;
    ae_int32x2 Z20, Z21, Z22, Z23, Z30, Z31, Z32, Z33;
    ae_int16x4 t0, t1, t2, t3, t4, t5, t6, t7;
    ae_int16x4 X0, X1, X2, X3;
    ae_int8x8 x0, x1, x2, x3;
    ae_int8x8 y0, y1, y2, y3;
    ae_int8x8 z0, z1, z2, z3;
    ae_int8x8 dsel_idx, dsel_cvt;
    static const uint8_t ALIGN(16) dsel_tbl[] =
    {
        (15<<4)|13, (7<<4)|5, (14<<4)|12, (6<<4)|4, (11<<4)|9, (3<<4)|1, (10<<4)|8, (2<<4)|0,
        (15<<4)|14, (13<<4)|12, (11<<4)|10, (9<<4)|8, (7<<4)|6, (5<<4)|4, (3<<4)|2, (1<<4)|0,
        (14<<4)|6, (12<<4)|4, (10<<4)|2, (8<<4)|0, (14<<4)|6, (12<<4)|4, (10<<4)|2, (8<<4)|0
    };

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
        return;
    }

    if (P == 4)
    {
        ae_valign alX1, alX3;

        dsel_idx = AE_L8X8_I((const ae_int8x8 *)dsel_tbl, 8);
        /* Process multiplication by 4x4 blocks */
        for (m = 0; m < (M>>2); m++)
        {
            pZ0 = z + m*4*P;

            pX0 = x + m*4*N;
            pX1 = pX0 + N;
            pX2 = pX1 + N;
            pX3 = pX2 + N;
            pY0 = y;

            alX1 = AE_LA64_PP(pX1);
            alX3 = AE_LA64_PP(pX3);
            Z00 = Z01 = Z10 = Z11 = Z20 = Z21 = Z30 = Z31 = AE_ZERO32();

            for (n = 0; n < (N>>3); n++)
            {
                /* load x matrix, 4x8 values, 8-bit */
                AE_L8X8_IP (x0,       castxcc(ae_int8x8,pX0), 8);
                AE_LA8X8_IP(x1, alX1, castxcc(ae_int8x8,pX1));
                AE_L8X8_IP (x2,       castxcc(ae_int8x8,pX2), 8);
                AE_LA8X8_IP(x3, alX3, castxcc(ae_int8x8,pX3));

                /* load y matrix, 8x4 values, 8-bit */
                AE_L8X8_IP(y0, castxcc(ae_int8x8,pY0), 8);
                AE_L8X8_IP(y1, castxcc(ae_int8x8,pY0), 8);
                AE_L8X8_IP(y2, castxcc(ae_int8x8,pY0), 8);
                AE_L8X8_IP(y3, castxcc(ae_int8x8,pY0), 8);
                /* transpose y matrix */
                AE_DSEL8X8(z0, z1, y0, y1, dsel_idx);
                AE_DSEL8X8(z2, z3, y2, y3, dsel_idx);
                AE_DSEL8X8(y0, y2, z0, z2, dsel_idx);
                AE_DSEL8X8(y1, y3, z1, z3, dsel_idx);

                /* make multiply, using 4-way 32-bit output octal mac, 8X8 bit */
                AE_MULA8Q8X8(Z00, Z01, y0, y1, y2, y3, x0);
                AE_MULA8Q8X8(Z10, Z11, y0, y1, y2, y3, x1);
                AE_MULA8Q8X8(Z20, Z21, y0, y1, y2, y3, x2);
                AE_MULA8Q8X8(Z30, Z31, y0, y1, y2, y3, x3);
            }
            if (N & 4)
            {
                ae_int32x2 T0, T1, T2, T3, T4, T5, T6, T7;
                /* load x matrix, 4x4 values, cvt 8->16-bit */
                X0 = AE_L8X4S_I(pX0, 0);
                X1 = AE_L8X4S_I(pX1, 0);
                X2 = AE_L8X4S_I(pX2, 0);
                X3 = AE_L8X4S_I(pX3, 0);

                /* load y matrix, 4x4 values, 8-bit */
                AE_L8X8_IP(y0, castxcc(ae_int8x8,pY0), 8);
                AE_L8X8_IP(y1, castxcc(ae_int8x8,pY0), 8);
                /* transpose y matrix */
                z0 = AE_SEL8X8I(y0, y1, 26);/* AE_SELI_8B_EXTRACT_ODD  */
                z1 = AE_SEL8X8I(y0, y1, 25);/* AE_SELI_8B_EXTRACT_EVEN */
                y0 = AE_SEL8X8I(z0, z0, 28);/* AE_SELI_8B_EXTRACT_1_OF_2_ODD_EVEN  */
                y1 = AE_SEL8X8I(z1, z1, 28);/* AE_SELI_8B_EXTRACT_1_OF_2_ODD_EVEN */

                /* make multiply, using 8-way 32-bit output quad mac, 8X16 bit */
                AE_MULA4O8X16(Z00, T0, Z01, T1, y0, y1, y0, y1, X0, X0);
                AE_MULA4O8X16(Z10, T2, Z11, T3, y0, y1, y0, y1, X1, X1);
                AE_MULA4O8X16(Z20, T4, Z21, T5, y0, y1, y0, y1, X2, X2);
                AE_MULA4O8X16(Z30, T6, Z31, T7, y0, y1, y0, y1, X3, X3);
            }
            /* Q15 + lsh <- Q14 + 1 + lsh w/ truncation and saturation */
            t0 = AE_TRUNCA16X4F32S(Z00, Z01, 16+1+lsh);
            t1 = AE_TRUNCA16X4F32S(Z10, Z11, 16+1+lsh);
            t2 = AE_TRUNCA16X4F32S(Z20, Z21, 16+1+lsh);
            t3 = AE_TRUNCA16X4F32S(Z30, Z31, 16+1+lsh);
            /* Q7 + lsh <- Q15 + lsh - 8 w/ rounding and saturation */
            z0 = AE_ROUND8X8F16SASYM(t0, t1);
            z1 = AE_ROUND8X8F16SASYM(t2, t3);
            AE_S8X8_IP(z0, castxcc(ae_int8x8,pZ0), 8);
            AE_S8X8_IP(z1, castxcc(ae_int8x8,pZ0), 8);
        }
    }
    else if (P & 4)/* P is more than 4 and is a multiple of 4 */
    {
        ae_valign alZ1, alZ3;
        alZ1 = alZ3 = AE_ZALIGN64();
        dsel_idx = AE_L8X8_I((const ae_int8x8 *)dsel_tbl, 0);
        dsel_cvt = AE_L8X8_I((const ae_int8x8 *)dsel_tbl, 16);
        /* Process multiplication by 4x8 blocks */
        for (m = 0; m < (M>>2); m++)
        {
            pZ0 = z + m*4*P;
            pZ1 = pZ0 + P;
            pZ2 = pZ1 + P;
            pZ3 = pZ2 + P;

            for (p = 0; p < (P>>3); p++)
            {
                pX0 = x + m*4*N;
                pY0 = y + p*8;
                Z00 = Z01 = Z02 = Z03 = Z10 = Z11 = Z12 = Z13 = AE_ZERO32();
                Z20 = Z21 = Z22 = Z23 = Z30 = Z31 = Z32 = Z33 = AE_ZERO32();
                __Pragma("loop_count min=1");
                for (n = 0; n < (N>>2); n++)
                {
                    ae_int16x4 Y10, Y11, Y30, Y31;
                    /* load x matrix, 4x4 values, cvt 8->16-bit */
                    AE_L8X4S_XP(X0, pX0, N);
                    AE_L8X4S_XP(X1, pX0, N);
                    AE_L8X4S_XP(X2, pX0, N);
                    AE_L8X4S_XP(X3, pX0, -N*3 + 4);

                    /* load y matrix, 4x8 values, 8-bit */
                    AE_L8X8_XP(y0, castxcc(ae_int8x8,pY0), P);
                    AE_L8X4S_IP(Y10, pY0, 4);  AE_L8X4S_XP(Y11, pY0, P-4);  y1 = AE_SAT8X8X16(Y10, Y11);
                    AE_L8X8_XP(y2, castxcc(ae_int8x8,pY0), P);
                    AE_L8X4S_IP(Y30, pY0, 4);  AE_L8X4S_XP(Y31, pY0, P-4);  y3 = AE_SAT8X8X16(Y30, Y31);
                    /* transpose y matrix */
                    AE_DSEL8X8(z0, z1, y0, y2, dsel_idx);
                    AE_DSEL8X8(z2, z3, y1, y3, dsel_idx);
                    AE_DSEL8X8(y0, y1, z0, z2, dsel_idx);
                    AE_DSEL8X8(y2, y3, z1, z3, dsel_idx);

                    /* make multiply, using 8-way 32-bit output quad mac, 8X16 bit */
                    AE_MULA4O8X16(Z00, Z01, Z02, Z03, y0, y1, y2, y3, X0, X0);
                    AE_MULA4O8X16(Z10, Z11, Z12, Z13, y0, y1, y2, y3, X1, X1);
                    AE_MULA4O8X16(Z20, Z21, Z22, Z23, y0, y1, y2, y3, X2, X2);
                    AE_MULA4O8X16(Z30, Z31, Z32, Z33, y0, y1, y2, y3, X3, X3);
                }
                /* Q15 + lsh <- Q14 + 1 + lsh w/ truncation and saturation */
                t0 = AE_TRUNCA16X4F32S(Z00, Z01, 16+1+lsh);
                t1 = AE_TRUNCA16X4F32S(Z02, Z03, 16+1+lsh);
                t2 = AE_TRUNCA16X4F32S(Z10, Z11, 16+1+lsh);
                t3 = AE_TRUNCA16X4F32S(Z12, Z13, 16+1+lsh);
                t4 = AE_TRUNCA16X4F32S(Z20, Z21, 16+1+lsh);
                t5 = AE_TRUNCA16X4F32S(Z22, Z23, 16+1+lsh);
                t6 = AE_TRUNCA16X4F32S(Z30, Z31, 16+1+lsh);
                t7 = AE_TRUNCA16X4F32S(Z32, Z33, 16+1+lsh);
                /* Q7 + lsh <- Q15 + lsh -8 w/ rounding and saturation */
                z0 = AE_ROUND8X8F16SASYM(t0, t1);
                z1 = AE_ROUND8X8F16SASYM(t2, t3);
                z2 = AE_ROUND8X8F16SASYM(t4, t5);
                z3 = AE_ROUND8X8F16SASYM(t6, t7);
                AE_S8X8_IP (z0,       castxcc(ae_int8x8,pZ0), 8);
                AE_SA8X8_IP(z1, alZ1, castxcc(ae_int8x8,pZ1));
                AE_S8X8_IP (z2,       castxcc(ae_int8x8,pZ2), 8);
                AE_SA8X8_IP(z3, alZ3, castxcc(ae_int8x8,pZ3));
            }
            AE_SA64POS_FP(alZ1, pZ1);
            AE_SA64POS_FP(alZ3, pZ3);
            /* Process last 4 columns */
            {
                pX0 = x + m*4*N;
                pY0 = y + p*8;
                Z00 = Z01 = Z10 = Z11 = Z20 = Z21 = Z30 = Z31 = AE_ZERO32();
                __Pragma("loop_count min=1");
                for (n = 0; n < (N>>2); n++)
                {
                    ae_int16x4 Y0, Y1, Y2, Y3;
                    /* load x matrix, 4x4 values, cvt 8->16-bit */
                    AE_L8X4S_XP(X0, pX0, N);
                    AE_L8X4S_XP(X1, pX0, N);
                    AE_L8X4S_XP(X2, pX0, N);
                    AE_L8X4S_XP(X3, pX0, -N*3 + 4);

                    /* load y matrix, 4x4 values, 8-bit */
                    AE_L8X4S_XP(Y0, pY0, P);
                    AE_L8X4S_XP(Y1, pY0, P);
                    AE_L8X4S_XP(Y2, pY0, P);
                    AE_L8X4S_XP(Y3, pY0, P);
                    AE_DSEL8X8(y0, y1, AE_MOVINT8X8_FROMINT16X4(Y0), AE_MOVINT8X8_FROMINT16X4(Y1), dsel_cvt);
                    AE_DSEL8X8(y2, y3, AE_MOVINT8X8_FROMINT16X4(Y2), AE_MOVINT8X8_FROMINT16X4(Y3), dsel_cvt);
                    /* transpose y matrix */
                    AE_DSEL8X8(z0, z1, y0, y2, dsel_idx);
                    AE_DSEL8X8(z2, z3, y1, y3, dsel_idx);
                    AE_DSEL8X8(y0, y1, z0, z2, dsel_idx);
                    AE_DSEL8X8(y2, y3, z1, z3, dsel_idx);

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
                AE_S32_H_I(AE_MOVINT32X2_FROMINT8X8(z0), (ae_int32 *)pZ0, 0);
                AE_S32_L_I(AE_MOVINT32X2_FROMINT8X8(z0), (ae_int32 *)pZ1, 0);
                AE_S32_H_I(AE_MOVINT32X2_FROMINT8X8(z1), (ae_int32 *)pZ2, 0);
                AE_S32_L_I(AE_MOVINT32X2_FROMINT8X8(z1), (ae_int32 *)pZ3, 0);
            }
        }
    }
    else/* P is a multiple of 8 */
    {
        dsel_idx = AE_L8X8_I((const ae_int8x8 *)dsel_tbl, 0);
        /* Process multiplication by 4x8 blocks */
        for (m = 0; m < (M>>2); m++)
        {
            pZ0 = z + m*4*P;

            for (p = 0; p < (P>>3); p++)
            {
                pX0 = x + m*4*N;
                pY0 = y + p*8;
                Z00 = Z01 = Z02 = Z03 = Z10 = Z11 = Z12 = Z13 = AE_ZERO32();
                Z20 = Z21 = Z22 = Z23 = Z30 = Z31 = Z32 = Z33 = AE_ZERO32();
                __Pragma("loop_count min=1");
                for (n = 0; n < (N>>2); n++)
                {
                    /* load x matrix, 4x4 values, cvt 8->16-bit */
                    AE_L8X4S_XP(X0, pX0, N);
                    AE_L8X4S_XP(X1, pX0, N);
                    AE_L8X4S_XP(X2, pX0, N);
                    AE_L8X4S_XP(X3, pX0, -N*3 + 4);

                    /* load y matrix, 4x8 values, 8-bit */
                    AE_L8X8_XP(y0, castxcc(ae_int8x8,pY0), P);
                    AE_L8X8_XP(y1, castxcc(ae_int8x8,pY0), P);
                    AE_L8X8_XP(y2, castxcc(ae_int8x8,pY0), P);
                    AE_L8X8_XP(y3, castxcc(ae_int8x8,pY0), P);
                    /* transpose y matrix */
                    AE_DSEL8X8(z0, z1, y0, y2, dsel_idx);
                    AE_DSEL8X8(z2, z3, y1, y3, dsel_idx);
                    AE_DSEL8X8(y0, y1, z0, z2, dsel_idx);
                    AE_DSEL8X8(y2, y3, z1, z3, dsel_idx);

                    /* make multiply, using 8-way 32-bit output quad mac, 8X16 bit */
                    AE_MULA4O8X16(Z00, Z01, Z02, Z03, y0, y1, y2, y3, X0, X0);
                    AE_MULA4O8X16(Z10, Z11, Z12, Z13, y0, y1, y2, y3, X1, X1);
                    AE_MULA4O8X16(Z20, Z21, Z22, Z23, y0, y1, y2, y3, X2, X2);
                    AE_MULA4O8X16(Z30, Z31, Z32, Z33, y0, y1, y2, y3, X3, X3);
                }
                /* Q15 + lsh <- Q14 + 1 + lsh w/ truncation and saturation */
                t0 = AE_TRUNCA16X4F32S(Z00, Z01, 16+1+lsh);
                t1 = AE_TRUNCA16X4F32S(Z02, Z03, 16+1+lsh);
                t2 = AE_TRUNCA16X4F32S(Z10, Z11, 16+1+lsh);
                t3 = AE_TRUNCA16X4F32S(Z12, Z13, 16+1+lsh);
                t4 = AE_TRUNCA16X4F32S(Z20, Z21, 16+1+lsh);
                t5 = AE_TRUNCA16X4F32S(Z22, Z23, 16+1+lsh);
                t6 = AE_TRUNCA16X4F32S(Z30, Z31, 16+1+lsh);
                t7 = AE_TRUNCA16X4F32S(Z32, Z33, 16+1+lsh);
                /* Q7 + lsh <- Q15 + lsh -8 w/ rounding and saturation */
                z0 = AE_ROUND8X8F16SASYM(t0, t1);
                z1 = AE_ROUND8X8F16SASYM(t2, t3);
                z2 = AE_ROUND8X8F16SASYM(t4, t5);
                z3 = AE_ROUND8X8F16SASYM(t6, t7);
                AE_S8X8_XP(z0, castxcc(ae_int8x8,pZ0), P);
                AE_S8X8_XP(z1, castxcc(ae_int8x8,pZ0), P);
                AE_S8X8_XP(z2, castxcc(ae_int8x8,pZ0), P);
                AE_S8X8_XP(z3, castxcc(ae_int8x8,pZ0), -P*3 + 8);
            }
        }
    }
}
#endif

size_t mtx_mpy8x8_fast_getScratchSize ( int M, int N, int P)
{
    (void)M;(void)N;(void)P;
    return 0;
}
