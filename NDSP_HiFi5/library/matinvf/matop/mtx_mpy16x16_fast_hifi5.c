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
void mtx_mpy16x16_fast(void* pScr,
    int16_t* restrict z,
    const int16_t* restrict x,
    const int16_t* restrict y,
    int M, int N, int P, int lsh)
{
    const ae_int16x4* restrict pY0;
    const ae_int16x8* restrict pX0;
    const ae_int16x8* restrict pX1;
    const ae_int16x8* restrict pX2;
    const ae_int16x8* restrict pX3;
    ae_int16x4* restrict pZ;
    int m, n, p;
    static const int16_t ALIGN(16) dsel_interleave1_tbl[]={6|(7<<8), 4|(5<<8), 2|(3<<8), 0|(1<<8)};
    ae_int16x4 dsel_interleave1=AE_L16X4_I((const ae_int16x4*)dsel_interleave1_tbl,0);

    NASSERT(lsh >= -15 && lsh <= 15);
    NASSERT_ALIGN(pScr,HIFI_SIMD_WIDTH);
    NASSERT_ALIGN(z,HIFI_SIMD_WIDTH);
    NASSERT_ALIGN(x,HIFI_SIMD_WIDTH);
    NASSERT_ALIGN(y,HIFI_SIMD_WIDTH);
    ASSERT((N & 3) == 0);
    ASSERT((M & 3) == 0);
    ASSERT((P & 3) == 0);
    if ((M <= 0) || (P <= 0)) return;
    if (N <= 0)
    {
        for (m = 0; m < ((M * P)>>3); m++) AE_S16X4X2_IP(0,0,castxcc(ae_int16x8,z),sizeof(ae_int16x8));
        return;
    }

    if (N&4)
    {
        // N is not a multiple of 8
        NASSERT((N+4)%8==0);
        __Pragma("loop_count min=1")
        for (p = 0; p<(P>>2); p++)
        {
            pZ=(ae_int16x4*)z;
            pX3=(const ae_int16x8*)x;
            __Pragma("loop_count min=1")
            for (m = 0; m<(M>>2); m++)
            {
                ae_valignx2 aX0,aX2;
                ae_int64 a0,a1,a2,a3,b0,b1,b2,b3,c0,c1,c2,c3,d0,d1,d2,d3;
                ae_int16x4 x0,x1,x2,x3,x4,x5,x6,x7,y0,y1,y2,y3;
                pY0=(const ae_int16x4*)(y);
                pX0=(const ae_int16x8*)pX3;
                pX1=(const ae_int16x8*)XT_ADDX2(N,(uintptr_t)pX0);
                pX2=(const ae_int16x8*)XT_ADDX4(N,(uintptr_t)pX0);
                pX3=(const ae_int16x8*)XT_ADDX4(N,(uintptr_t)pX1);

                AE_L16X4_IP(x0,castxcc(ae_int16x4,pX0),sizeof(ae_int16x4));
                AE_L16X4_IP(x1,castxcc(ae_int16x4,pX1),sizeof(ae_int16x4));
                AE_L16X4_IP(x2,castxcc(ae_int16x4,pX2),sizeof(ae_int16x4));
                AE_L16X4_IP(x3,castxcc(ae_int16x4,pX3),sizeof(ae_int16x4));

                AE_L16X4_XP(y0,pY0,P*sizeof(int16_t));
                AE_L16X4_XP(y1,pY0,P*sizeof(int16_t));
                AE_L16X4_XP(y2,pY0,P*sizeof(int16_t));
                AE_L16X4_XP(y3,pY0,P*sizeof(int16_t));

                AE_DSEL16X4(y0,y1,y0,y1,dsel_interleave1);  
                AE_DSEL16X4(y2,y3,y2,y3,dsel_interleave1);  
                AE_DSEL16X4(y0,y2,y0,y2,dsel_interleave1);  
                AE_DSEL16X4(y1,y3,y1,y3,dsel_interleave1);  

                AE_MULZAAAA2Q16(a0,a1,x0,x0,y0,y1);
                AE_MULZAAAA2Q16(a2,a3,x0,x0,y2,y3);
                AE_MULZAAAA2Q16(b0,b1,x1,x1,y0,y1);
                AE_MULZAAAA2Q16(b2,b3,x1,x1,y2,y3);
                AE_MULZAAAA2Q16(c0,c1,x2,x2,y0,y1);
                AE_MULZAAAA2Q16(c2,c3,x2,x2,y2,y3);
                AE_MULZAAAA2Q16(d0,d1,x3,x3,y0,y1);
                AE_MULZAAAA2Q16(d2,d3,x3,x3,y2,y3);
                aX0=AE_LA128_PP(pX0);
                aX2=AE_LA128_PP(pX2);
                for (n = 0; n<(N>>3); n++)
                {
                    AE_LA16X4X2_IP(x0,x4,aX0,pX0);
                    AE_L16X4X2_IP (x1,x5,pX1,sizeof(ae_int16x8));
                    AE_LA16X4X2_IP (x2,x6,aX2,pX2);
                    AE_L16X4X2_IP(x3,x7,pX3,sizeof(ae_int16x8));

                    AE_L16X4_XP(y0,pY0,P*sizeof(int16_t));
                    AE_L16X4_XP(y1,pY0,P*sizeof(int16_t));
                    AE_L16X4_XP(y2,pY0,P*sizeof(int16_t));
                    AE_L16X4_XP(y3,pY0,P*sizeof(int16_t));

                    AE_DSEL16X4(y0,y1,y0,y1,dsel_interleave1);  
                    AE_DSEL16X4(y2,y3,y2,y3,dsel_interleave1);  
                    AE_DSEL16X4(y0,y2,y0,y2,dsel_interleave1);  
                    AE_DSEL16X4(y1,y3,y1,y3,dsel_interleave1);  

                    AE_MULAAAA2Q16(a0,a1,x0,x0,y0,y1);
                    AE_MULAAAA2Q16(a2,a3,x0,x0,y2,y3);
                    AE_MULAAAA2Q16(b0,b1,x1,x1,y0,y1);
                    AE_MULAAAA2Q16(b2,b3,x1,x1,y2,y3);
                    AE_MULAAAA2Q16(c0,c1,x2,x2,y0,y1);
                    AE_MULAAAA2Q16(c2,c3,x2,x2,y2,y3);
                    AE_MULAAAA2Q16(d0,d1,x3,x3,y0,y1);
                    AE_MULAAAA2Q16(d2,d3,x3,x3,y2,y3);

                    AE_L16X4_XP(y0,pY0,P*sizeof(int16_t));
                    AE_L16X4_XP(y1,pY0,P*sizeof(int16_t));
                    AE_L16X4_XP(y2,pY0,P*sizeof(int16_t));
                    AE_L16X4_XP(y3,pY0,P*sizeof(int16_t));

                    AE_DSEL16X4(y0,y1,y0,y1,dsel_interleave1);  
                    AE_DSEL16X4(y2,y3,y2,y3,dsel_interleave1);  
                    AE_DSEL16X4(y0,y2,y0,y2,dsel_interleave1);  
                    AE_DSEL16X4(y1,y3,y1,y3,dsel_interleave1);  

                    AE_MULAAAA2Q16(a0,a1,x4,x4,y0,y1);
                    AE_MULAAAA2Q16(a2,a3,x4,x4,y2,y3);
                    AE_MULAAAA2Q16(b0,b1,x5,x5,y0,y1);
                    AE_MULAAAA2Q16(b2,b3,x5,x5,y2,y3);
                    AE_MULAAAA2Q16(c0,c1,x6,x6,y0,y1);
                    AE_MULAAAA2Q16(c2,c3,x6,x6,y2,y3);
                    AE_MULAAAA2Q16(d0,d1,x7,x7,y0,y1);
                    AE_MULAAAA2Q16(d2,d3,x7,x7,y2,y3);
                }
                AE_S16X4_XP(AE_ROUND16X4F32SASYM(AE_TRUNCA32X2F64S(a0, a1, lsh + 33), AE_TRUNCA32X2F64S(a2, a3, lsh + 33)),pZ,P*sizeof(int16_t));
                AE_S16X4_XP(AE_ROUND16X4F32SASYM(AE_TRUNCA32X2F64S(b0, b1, lsh + 33), AE_TRUNCA32X2F64S(b2, b3, lsh + 33)),pZ,P*sizeof(int16_t));
                AE_S16X4_XP(AE_ROUND16X4F32SASYM(AE_TRUNCA32X2F64S(c0, c1, lsh + 33), AE_TRUNCA32X2F64S(c2, c3, lsh + 33)),pZ,P*sizeof(int16_t));
                AE_S16X4_XP(AE_ROUND16X4F32SASYM(AE_TRUNCA32X2F64S(d0, d1, lsh + 33), AE_TRUNCA32X2F64S(d2, d3, lsh + 33)),pZ,P*sizeof(int16_t));
            }
            z+=4;
            y+=4;
        }
    }
    else
    {
        // N is a multiple of 8
        NASSERT(N%8==0);
        __Pragma("loop_count min=1")
        for (p = 0; p<(P>>2); p++)
        {
            pZ=(ae_int16x4*)z;
            pX3=(const ae_int16x8*)x;
            __Pragma("loop_count min=1")
            for (m = 0; m<(M>>2); m++)
            {
                ae_int64 a0,a1,a2,a3,b0,b1,b2,b3,c0,c1,c2,c3,d0,d1,d2,d3;
                pY0=(const ae_int16x4*)(y);
                pX0=(const ae_int16x8*)pX3;
                pX1=(const ae_int16x8*)XT_ADDX2(N,(uintptr_t)pX0);
                pX2=(const ae_int16x8*)XT_ADDX4(N,(uintptr_t)pX0);
                pX3=(const ae_int16x8*)XT_ADDX4(N,(uintptr_t)pX1);

                AE_MOVDX2(a0,a1,0,0); AE_MOVDX2(a2,a3,0,0); 
                AE_MOVDX2(b0,b1,0,0); AE_MOVDX2(b2,b3,0,0); 
                AE_MOVDX2(c0,c1,0,0); AE_MOVDX2(c2,c3,0,0); 
                AE_MOVDX2(d0,d1,0,0); AE_MOVDX2(d2,d3,0,0);
                __Pragma("loop_count min=1")
                for (n = 0; n<(N>>3); n++)
                {
                    ae_int16x4 x0,x1,x2,x3,x4,x5,x6,x7,y0,y1,y2,y3;

                    AE_L16X4X2_IP (x0,x4,pX0,sizeof(ae_int16x8));
                    AE_L16X4X2_IP (x1,x5,pX1,sizeof(ae_int16x8));
                    AE_L16X4X2_IP (x2,x6,pX2,sizeof(ae_int16x8));
                    AE_L16X4X2_IP (x3,x7,pX3,sizeof(ae_int16x8));

                    AE_L16X4_XP(y0,pY0,P*sizeof(int16_t));
                    AE_L16X4_XP(y1,pY0,P*sizeof(int16_t));
                    AE_L16X4_XP(y2,pY0,P*sizeof(int16_t));
                    AE_L16X4_XP(y3,pY0,P*sizeof(int16_t));

                    AE_DSEL16X4(y0,y1,y0,y1,dsel_interleave1);  
                    AE_DSEL16X4(y2,y3,y2,y3,dsel_interleave1);  
                    AE_DSEL16X4(y0,y2,y0,y2,dsel_interleave1);  
                    AE_DSEL16X4(y1,y3,y1,y3,dsel_interleave1);  

                    AE_MULAAAA2Q16(a0,a1,x0,x0,y0,y1);
                    AE_MULAAAA2Q16(a2,a3,x0,x0,y2,y3);
                    AE_MULAAAA2Q16(b0,b1,x1,x1,y0,y1);
                    AE_MULAAAA2Q16(b2,b3,x1,x1,y2,y3);
                    AE_MULAAAA2Q16(c0,c1,x2,x2,y0,y1);
                    AE_MULAAAA2Q16(c2,c3,x2,x2,y2,y3);
                    AE_MULAAAA2Q16(d0,d1,x3,x3,y0,y1);
                    AE_MULAAAA2Q16(d2,d3,x3,x3,y2,y3);

                    AE_L16X4_XP(y0,pY0,P*sizeof(int16_t));
                    AE_L16X4_XP(y1,pY0,P*sizeof(int16_t));
                    AE_L16X4_XP(y2,pY0,P*sizeof(int16_t));
                    AE_L16X4_XP(y3,pY0,P*sizeof(int16_t));

                    AE_DSEL16X4(y0,y1,y0,y1,dsel_interleave1);  
                    AE_DSEL16X4(y2,y3,y2,y3,dsel_interleave1);  
                    AE_DSEL16X4(y0,y2,y0,y2,dsel_interleave1);  
                    AE_DSEL16X4(y1,y3,y1,y3,dsel_interleave1);  

                    AE_MULAAAA2Q16(a0,a1,x4,x4,y0,y1);
                    AE_MULAAAA2Q16(a2,a3,x4,x4,y2,y3);
                    AE_MULAAAA2Q16(b0,b1,x5,x5,y0,y1);
                    AE_MULAAAA2Q16(b2,b3,x5,x5,y2,y3);
                    AE_MULAAAA2Q16(c0,c1,x6,x6,y0,y1);
                    AE_MULAAAA2Q16(c2,c3,x6,x6,y2,y3);
                    AE_MULAAAA2Q16(d0,d1,x7,x7,y0,y1);
                    AE_MULAAAA2Q16(d2,d3,x7,x7,y2,y3);
                }
                AE_S16X4_XP(AE_ROUND16X4F32SASYM(AE_TRUNCA32X2F64S(a0, a1, lsh + 33), AE_TRUNCA32X2F64S(a2, a3, lsh + 33)),pZ,P*sizeof(int16_t));
                AE_S16X4_XP(AE_ROUND16X4F32SASYM(AE_TRUNCA32X2F64S(b0, b1, lsh + 33), AE_TRUNCA32X2F64S(b2, b3, lsh + 33)),pZ,P*sizeof(int16_t));
                AE_S16X4_XP(AE_ROUND16X4F32SASYM(AE_TRUNCA32X2F64S(c0, c1, lsh + 33), AE_TRUNCA32X2F64S(c2, c3, lsh + 33)),pZ,P*sizeof(int16_t));
                AE_S16X4_XP(AE_ROUND16X4F32SASYM(AE_TRUNCA32X2F64S(d0, d1, lsh + 33), AE_TRUNCA32X2F64S(d2, d3, lsh + 33)),pZ,P*sizeof(int16_t));
            }
            z+=4;
            y+=4;
        }
    }
} 

size_t mtx_mpy16x16_fast_getScratchSize(int M, int N, int P)
{
    (void)M; (void)N; (void)P;
    return 0;
}
