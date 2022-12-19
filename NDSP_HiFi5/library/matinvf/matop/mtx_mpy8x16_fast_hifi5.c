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
#if !defined(AE_MULA4O8X16) || !defined(AE_MULA8Q8X16)
void mtx_mpy8x16_fast(  void* pScr,
                     int16_t* restrict z,
               const  int8_t* restrict x,
               const int16_t* restrict y,
               int M, int N, int P, int lsh )
{
    static const int16_t ALIGN(16) dsel_interleave1_tbl[]={6|(7<<8), 4|(5<<8), 2|(3<<8), 0|(1<<8)};
    ae_int16x4 dsel_interleave1;
    int m, n, p,sa=lsh+33;
    const ae_int16x4 * restrict pX;
    const ae_int16x4 * restrict pX0;
    const ae_int16x4 * restrict pX1;
    const ae_int16x4 * restrict pX2;
    const ae_int16x4 * restrict pX3;
    const ae_int32   * restrict pY0;
    const ae_int32   * restrict pY1;
          ae_int32   * restrict pZ;
    ae_int16x4 x0, x1, x2, x3;

    (void)pScr;
    NASSERT(N % 4 == 0);
    NASSERT(M % 4 == 0);
    NASSERT(P % 4 == 0);
    NASSERT_ALIGN(x,HIFI_SIMD_WIDTH);
    NASSERT_ALIGN(y,HIFI_SIMD_WIDTH);
    NASSERT_ALIGN(z,HIFI_SIMD_WIDTH);
    if ((M <= 0) || (P <= 0)) return;
    if (N <= 0)
    {
        for (m = 0; m < ((M * P)>>3); m++) AE_S16X4X2_IP(0,0,castxcc(ae_int16x8,z),sizeof (ae_int16x8));
        return;
    }
    dsel_interleave1=AE_L16X4_I((const ae_int16x4*)dsel_interleave1_tbl,0);

    __Pragma("loop_count min=1");
    for (p = 0; p < (P>>1); p+=2)
    {
        pX = (const ae_int16x4 *)x;
        pZ = (      ae_int32   *)z;
        __Pragma("ymemory(pY0)")
        __Pragma("ymemory(pY1)")
        __Pragma("loop_count min=1");
        for (m = 0; m < (M>>2); m++)
        {
            ae_int16x4 t04,t15,t26,t37,s0,s1,s2,s3,v0,v1,v2,v3;
            ae_int64 B0, B1, B2, B3, B4, B5, B6, B7;
            ae_int64 C0, C1, C2, C3, C4, C5, C6, C7;
            pY0 = (const ae_int32   *)(y);
            pY1 = (const ae_int32   *)(y+2*P);
            pX0 = (const ae_int16x4 *)pX;
            pX1 = (const ae_int16x4 *)XT_ADD(N, (uintptr_t)pX0);
            pX2 = (const ae_int16x4 *)XT_ADD(N, (uintptr_t)pX1);
            pX3 = (const ae_int16x4 *)XT_ADD(N, (uintptr_t)pX2);

            t15=AE_L16X4_X ((const ae_int16x4*)pY0 ,  P*sizeof(int16_t));
            AE_L16X4_XP(t04,castxcc(ae_int16x4,pY0),4*P*sizeof(int16_t));
            t37=AE_L16X4_X ((const ae_int16x4*)pY1 ,  P*sizeof(int16_t));
            AE_L16X4_XP(t26,castxcc(ae_int16x4,pY1),4*P*sizeof(int16_t));

            AE_DSEL16X4(s0,s1,t04,t15,dsel_interleave1);
            AE_DSEL16X4(s2,s3,t26,t37,dsel_interleave1);
            AE_DSEL16X4(v0,v2,s0,s2,dsel_interleave1);
            AE_DSEL16X4(v1,v3,s1,s3,dsel_interleave1);

            AE_L8X4F_IP(x0, castxcc(int8_t,pX0), 4*sizeof(int8_t));
            AE_L8X4F_IP(x1, castxcc(int8_t,pX1), 4*sizeof(int8_t));
            AE_L8X4F_IP(x2, castxcc(int8_t,pX2), 4*sizeof(int8_t));
            AE_L8X4F_IP(x3, castxcc(int8_t,pX3), 4*sizeof(int8_t));

            AE_MULZAAAA2Q16(B0,C0, x0,x0, v0,v2);
            AE_MULZAAAA2Q16(B1,C1, x0,x0, v1,v3);
            AE_MULZAAAA2Q16(B2,C2, x1,x1, v0,v2);
            AE_MULZAAAA2Q16(B3,C3, x1,x1, v1,v3);
            AE_MULZAAAA2Q16(B4,C4, x2,x2, v0,v2);
            AE_MULZAAAA2Q16(B5,C5, x2,x2, v1,v3);
            AE_MULZAAAA2Q16(B6,C6, x3,x3, v0,v2);
            AE_MULZAAAA2Q16(B7,C7, x3,x3, v1,v3);
            for (n = 0; n < (N>>2)-1; n++)
            {
                t15=AE_L16X4_X ((const ae_int16x4*)pY0 ,  P*sizeof(int16_t));
                AE_L16X4_XP(t04,castxcc(ae_int16x4,pY0),4*P*sizeof(int16_t));
                t37=AE_L16X4_X ((const ae_int16x4*)pY1 ,  P*sizeof(int16_t));
                AE_L16X4_XP(t26,castxcc(ae_int16x4,pY1),4*P*sizeof(int16_t));
                AE_L8X4F_IP(x0, castxcc(int8_t,pX0), 4*sizeof(int8_t));
                AE_L8X4F_IP(x1, castxcc(int8_t,pX1), 4*sizeof(int8_t));
                AE_L8X4F_IP(x2, castxcc(int8_t,pX2), 4*sizeof(int8_t));
                AE_L8X4F_IP(x3, castxcc(int8_t,pX3), 4*sizeof(int8_t));

                AE_DSEL16X4(s0,s1,t04,t15,dsel_interleave1);
                AE_DSEL16X4(s2,s3,t26,t37,dsel_interleave1);
                AE_DSEL16X4(v0,v2,s0,s2,dsel_interleave1);
                AE_DSEL16X4(v1,v3,s1,s3,dsel_interleave1);

                AE_MULAAAA2Q16(B0,C0, x0,x0, v0,v2);
                AE_MULAAAA2Q16(B1,C1, x0,x0, v1,v3);
                AE_MULAAAA2Q16(B2,C2, x1,x1, v0,v2);
                AE_MULAAAA2Q16(B3,C3, x1,x1, v1,v3);
                AE_MULAAAA2Q16(B4,C4, x2,x2, v0,v2);
                AE_MULAAAA2Q16(B5,C5, x2,x2, v1,v3);
                AE_MULAAAA2Q16(B6,C6, x3,x3, v0,v2);
                AE_MULAAAA2Q16(B7,C7, x3,x3, v1,v3);
            }
            AE_S16X4_XP(AE_ROUND16X4F32SASYM(AE_TRUNCA32X2F64S(B0, B1, sa), AE_TRUNCA32X2F64S(C0, C1, sa)), castxcc(ae_int16x4,pZ), P*sizeof(int16_t));
            AE_S16X4_XP(AE_ROUND16X4F32SASYM(AE_TRUNCA32X2F64S(B2, B3, sa), AE_TRUNCA32X2F64S(C2, C3, sa)), castxcc(ae_int16x4,pZ), P*sizeof(int16_t));
            AE_S16X4_XP(AE_ROUND16X4F32SASYM(AE_TRUNCA32X2F64S(B4, B5, sa), AE_TRUNCA32X2F64S(C4, C5, sa)), castxcc(ae_int16x4,pZ), P*sizeof(int16_t));
            AE_S16X4_XP(AE_ROUND16X4F32SASYM(AE_TRUNCA32X2F64S(B6, B7, sa), AE_TRUNCA32X2F64S(C6, C7, sa)), castxcc(ae_int16x4,pZ), P*sizeof(int16_t));
            pX = (const ae_int16x4 *)XT_ADDX4(N, (uintptr_t)pX);
        }
        z += 4;
        y += 4;
    }
}
#else
void mtx_mpy8x16_fast(  void* pScr,
                     int16_t* restrict z,
               const  int8_t* restrict x,
               const int16_t* restrict y,
               int M, int N, int P, int lsh )
{
    static const int16_t ALIGN(16) dsel_interleave1_tbl[]={6|(7<<8), 4|(5<<8), 2|(3<<8), 0|(1<<8)};
    int m, n, p;
    const int8_t  * restrict pX0;
    const int8_t  * restrict pX1;
    const int8_t  * restrict pX2;
    const int8_t  * restrict pX3;
    const int16_t * restrict pY0;
          int16_t * restrict pZ0;
    ae_valign alX1, alX3;
    ae_int8x8 x0, x1, x2, x3;
    ae_int16x4 y0, y1, y2, y3, y4, y5, y6, y7;
    ae_int32x2 Z00, Z01, Z10, Z11, Z20, Z21, Z30, Z31;
    ae_int16x4 z0, z1, z2, z3;
    ae_int16x4 dsel_interleave1;

    (void)pScr;
    NASSERT(N % 4 == 0);
    NASSERT(M % 4 == 0);
    NASSERT(P % 4 == 0);
    NASSERT_ALIGN(x,HIFI_SIMD_WIDTH);
    NASSERT_ALIGN(y,HIFI_SIMD_WIDTH);
    NASSERT_ALIGN(z,HIFI_SIMD_WIDTH);
    if ((M <= 0) || (P <= 0)) return;
    if (N <= 0)
    {
        for (m = 0; m < ((M * P)>>3); m++) AE_S16X4X2_IP(0,0,castxcc(ae_int16x8,z),sizeof (ae_int16x8));
        return;
    }

    WUR_AE_SAR(9+lsh);
    dsel_interleave1=AE_L16X4_I((const ae_int16x4*)dsel_interleave1_tbl,0);

    if (N & 4)
    {
        __Pragma("loop_count min=1");
        for (m = 0; m < (M>>2); m++)
        {
            pZ0 = z + m*P*4;
            __Pragma("loop_count min=1");
            for (p = 0; p < (P>>2); p++)
            {
                pX0 = x + m*N*4;
                pX1 = pX0 + N;
                pX2 = pX1 + N;
                pX3 = pX2 + N;
                pY0 = y + p*4;
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

                    /* load y matrix, 8x4 values, 16-bit */
                    AE_L16X4_XP(y0, castxcc(ae_int16x4,pY0), P*sizeof(int16_t));
                    AE_L16X4_XP(y1, castxcc(ae_int16x4,pY0), P*sizeof(int16_t));
                    AE_L16X4_XP(y2, castxcc(ae_int16x4,pY0), P*sizeof(int16_t));
                    AE_L16X4_XP(y3, castxcc(ae_int16x4,pY0), P*sizeof(int16_t));
                    AE_L16X4_XP(y4, castxcc(ae_int16x4,pY0), P*sizeof(int16_t));
                    AE_L16X4_XP(y5, castxcc(ae_int16x4,pY0), P*sizeof(int16_t));
                    AE_L16X4_XP(y6, castxcc(ae_int16x4,pY0), P*sizeof(int16_t));
                    AE_L16X4_XP(y7, castxcc(ae_int16x4,pY0), P*sizeof(int16_t));

                    AE_DSEL16X4(y0, y1, y0, y1, dsel_interleave1);  
                    AE_DSEL16X4(y2, y3, y2, y3, dsel_interleave1);  
                    AE_DSEL16X4(y0, y2, y0, y2, dsel_interleave1);  
                    AE_DSEL16X4(y1, y3, y1, y3, dsel_interleave1);  
                    AE_DSEL16X4(y4, y5, y4, y5, dsel_interleave1);  
                    AE_DSEL16X4(y6, y7, y6, y7, dsel_interleave1);  
                    AE_DSEL16X4(y4, y6, y4, y6, dsel_interleave1);  
                    AE_DSEL16X4(y5, y7, y5, y7, dsel_interleave1);  

                    /* make multiply, using 4-way 32-bit output octal mac, 8X16->32 bit */
                    AE_MULA8Q8X16(Z00, Z01, x0, x1, x2, x3, y0, y4);
                    AE_MULA8Q8X16(Z10, Z11, x0, x1, x2, x3, y1, y5);
                    AE_MULA8Q8X16(Z20, Z21, x0, x1, x2, x3, y2, y6);
                    AE_MULA8Q8X16(Z30, Z31, x0, x1, x2, x3, y3, y7);
                }
                /* Last 4 multiplications */
                {
                    ae_int16x4 X0, X1, X2, X3;
                    /* load x matrix, 4x4 values, 8-bit */
                    X0 = AE_L8X4S_I(pX0, 0);
                    X1 = AE_L8X4S_I(pX1, 0);
                    X2 = AE_L8X4S_I(pX2, 0);
                    X3 = AE_L8X4S_I(pX3, 0);
                    x0 = AE_SAT8X8X16(X0, X0);
                    x1 = AE_SAT8X8X16(X1, X1);
                    x2 = AE_SAT8X8X16(X2, X2);
                    x3 = AE_SAT8X8X16(X3, X3);

                    /* load y matrix, 4x4 values, 16-bit */
                    AE_L16X4_XP(y0, castxcc(ae_int16x4,pY0), P*sizeof(int16_t));
                    AE_L16X4_XP(y1, castxcc(ae_int16x4,pY0), P*sizeof(int16_t));
                    AE_L16X4_XP(y2, castxcc(ae_int16x4,pY0), P*sizeof(int16_t));
                    AE_L16X4_XP(y3, castxcc(ae_int16x4,pY0), P*sizeof(int16_t));

                    AE_DSEL16X4(y0, y1, y0, y1, dsel_interleave1);  
                    AE_DSEL16X4(y2, y3, y2, y3, dsel_interleave1);  
                    AE_DSEL16X4(y0, y2, y0, y2, dsel_interleave1);  
                    AE_DSEL16X4(y1, y3, y1, y3, dsel_interleave1);  

                    /* make multiply, using 4-way 32-bit output octal mac, 8X16->32 bit */
                    AE_MULA4O8X16(Z00, Z01, Z10, Z11, x0, x1, x2, x3, y0, y1);
                    AE_MULA4O8X16(Z20, Z21, Z30, Z31, x0, x1, x2, x3, y2, y3);
                }
                Z00 = AE_SLAS32S(Z00);
                Z01 = AE_SLAS32S(Z01);
                Z10 = AE_SLAS32S(Z10);
                Z11 = AE_SLAS32S(Z11);
                Z20 = AE_SLAS32S(Z20);
                Z21 = AE_SLAS32S(Z21);
                Z30 = AE_SLAS32S(Z30);
                Z31 = AE_SLAS32S(Z31);
                z0 = AE_ROUND16X4F32SASYM(Z00, Z10);
                z1 = AE_ROUND16X4F32SASYM(Z20, Z30);
                z2 = AE_ROUND16X4F32SASYM(Z01, Z11);
                z3 = AE_ROUND16X4F32SASYM(Z21, Z31);
                AE_DSEL16X4(z0, z1, z0, z1, dsel_interleave1);  
                AE_DSEL16X4(z2, z3, z2, z3, dsel_interleave1);  
                AE_S16X4_XP(z0, castxcc(ae_int16x4,pZ0), P*sizeof(int16_t));
                AE_S16X4_XP(z1, castxcc(ae_int16x4,pZ0), P*sizeof(int16_t));
                AE_S16X4_XP(z2, castxcc(ae_int16x4,pZ0), P*sizeof(int16_t));
                AE_S16X4_XP(z3, castxcc(ae_int16x4,pZ0), (-3*P+4)*(int)sizeof(int16_t));
            }
        }
    }
    else
    {
        __Pragma("loop_count min=1");
        for (m = 0; m < (M>>2); m++)
        {
            pZ0 = z + m*P*4;
            __Pragma("loop_count min=1");
            for (p = 0; p < (P>>2); p++)
            {
                pX0 = x + m*N*4;
                pY0 = y + p*4;

                Z00 = Z01 = Z10 = Z11 = Z20 = Z21 = Z30 = Z31 = AE_ZERO32();
                __Pragma("loop_count min=1");
                for (n = 0; n < (N>>3); n++)
                {
                    /* load x matrix, 4x8 values, 8-bit */
                    AE_L8X8_XP(x0, castxcc(ae_int8x8,pX0), N);
                    AE_L8X8_XP(x1, castxcc(ae_int8x8,pX0), N);
                    AE_L8X8_XP(x2, castxcc(ae_int8x8,pX0), N);
                    AE_L8X8_XP(x3, castxcc(ae_int8x8,pX0), -3*N+8);

                    /* load y matrix, 8x4 values, 16-bit */
                    AE_L16X4_XP(y0, castxcc(ae_int16x4,pY0), P*sizeof(int16_t));
                    AE_L16X4_XP(y1, castxcc(ae_int16x4,pY0), P*sizeof(int16_t));
                    AE_L16X4_XP(y2, castxcc(ae_int16x4,pY0), P*sizeof(int16_t));
                    AE_L16X4_XP(y3, castxcc(ae_int16x4,pY0), P*sizeof(int16_t));
                    AE_L16X4_XP(y4, castxcc(ae_int16x4,pY0), P*sizeof(int16_t));
                    AE_L16X4_XP(y5, castxcc(ae_int16x4,pY0), P*sizeof(int16_t));
                    AE_L16X4_XP(y6, castxcc(ae_int16x4,pY0), P*sizeof(int16_t));
                    AE_L16X4_XP(y7, castxcc(ae_int16x4,pY0), P*sizeof(int16_t));

                    AE_DSEL16X4(y0, y1, y0, y1, dsel_interleave1);  
                    AE_DSEL16X4(y2, y3, y2, y3, dsel_interleave1);  
                    AE_DSEL16X4(y0, y2, y0, y2, dsel_interleave1);  
                    AE_DSEL16X4(y1, y3, y1, y3, dsel_interleave1);  
                    AE_DSEL16X4(y4, y5, y4, y5, dsel_interleave1);  
                    AE_DSEL16X4(y6, y7, y6, y7, dsel_interleave1);  
                    AE_DSEL16X4(y4, y6, y4, y6, dsel_interleave1);  
                    AE_DSEL16X4(y5, y7, y5, y7, dsel_interleave1);  

                    /* make multiply, using 4-way 32-bit output octal mac, 8X16->32 bit */
                    AE_MULA8Q8X16(Z00, Z01, x0, x1, x2, x3, y0, y4);
                    AE_MULA8Q8X16(Z10, Z11, x0, x1, x2, x3, y1, y5);
                    AE_MULA8Q8X16(Z20, Z21, x0, x1, x2, x3, y2, y6);
                    AE_MULA8Q8X16(Z30, Z31, x0, x1, x2, x3, y3, y7);
                }
                Z00 = AE_SLAS32S(Z00);
                Z01 = AE_SLAS32S(Z01);
                Z10 = AE_SLAS32S(Z10);
                Z11 = AE_SLAS32S(Z11);
                Z20 = AE_SLAS32S(Z20);
                Z21 = AE_SLAS32S(Z21);
                Z30 = AE_SLAS32S(Z30);
                Z31 = AE_SLAS32S(Z31);
                z0 = AE_ROUND16X4F32SASYM(Z00, Z10);
                z1 = AE_ROUND16X4F32SASYM(Z20, Z30);
                z2 = AE_ROUND16X4F32SASYM(Z01, Z11);
                z3 = AE_ROUND16X4F32SASYM(Z21, Z31);
                AE_DSEL16X4(z0, z1, z0, z1, dsel_interleave1);  
                AE_DSEL16X4(z2, z3, z2, z3, dsel_interleave1);  
                AE_S16X4_XP(z0, castxcc(ae_int16x4,pZ0), P*sizeof(int16_t));
                AE_S16X4_XP(z1, castxcc(ae_int16x4,pZ0), P*sizeof(int16_t));
                AE_S16X4_XP(z2, castxcc(ae_int16x4,pZ0), P*sizeof(int16_t));
                AE_S16X4_XP(z3, castxcc(ae_int16x4,pZ0), (-3*P+4)*(int)sizeof(int16_t));
            }
        }
    }
}
#endif

size_t mtx_mpy8x16_fast_getScratchSize(int M, int N, int P)
{
    (void)M;(void)N;(void)P;
    return 0;
}
