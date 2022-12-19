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
#include "mtx_mpy8x8_helper.h"

#if USE_NN_EXTENSION_8X8==0
static const uint64_t dsel_rephhll_tbl=0x40516273c8d9eafbULL;

/* special code for multiplication of small matrices without scratch */
static void mtx_mpyt8x8_small(
                     int8_t* restrict z,
               const int8_t* restrict x,
               const int8_t* restrict y,
               int M, int N, int P, int lsh )
{
    int sa=49+lsh;
    ae_int8x8 dsel_rephhll;
    dsel_rephhll=AE_L8X8_I((const ae_int8x8*)&dsel_rephhll_tbl,0);
    int m,n,p;
    const ae_int8 * restrict pX0;
    const ae_int8 * restrict pX1;
    const ae_int8 * restrict pY0;
    ae_int64 A00,A10,A01,A11,t0,t1,t2,t3;

    AE_MOVDX2(t0,t1,0,0);
    AE_MOVDX2(t2,t3,0,0);
    pX0=(const ae_int8*)&x[0*N];
    pX1=(const ae_int8*)&x[1*N];
    for (m=0; m<(M>>1);m++)
    {
        pY0=(const ae_int8*)&y[0];
        for (p=0; p<(P>>1); p++)
        {
            ae_int8x8 z00,z01,z10,z11;
            AE_MOVDX2(A00,A01,0,0);
            AE_MOVDX2(A10,A11,0,0);
            __Pragma("loop_count min=1")
            for (n=0; n<N; n++)
            {
                ae_int8x8 x01,t01,x0,x1,y0,y1;
                AE_L8_IP(x0,castxcc(ae_int8,pX0),1);
                AE_L8_IP(x1,castxcc(ae_int8,pX1),1);
                y1=AE_L8_X ((const ae_int8*)pY0 ,N);
                AE_L8_IP(y0,castxcc(ae_int8,pY0),1);
                AE_DSEL8X8(x01,t01,x0,x1,dsel_rephhll);
                AE_MULAAAA2Q8(A00,A10,x01,y0);
                AE_MULAAAA2Q8(A01,A11,x01,y1);
            }
            z00=AE_ROUND8X4F32SASYM_L(AE_TRUNCA32X2F64S(A00,A00, sa-2),AE_TRUNCA32X2F64S(A00,A00, sa-2));
            z01=AE_ROUND8X4F32SASYM_L(AE_TRUNCA32X2F64S(A01,A01, sa-2),AE_TRUNCA32X2F64S(A01,A01, sa-2));
            z10=AE_ROUND8X4F32SASYM_L(AE_TRUNCA32X2F64S(A10,A10, sa-2),AE_TRUNCA32X2F64S(A10,A10, sa-2));
            z11=AE_ROUND8X4F32SASYM_L(AE_TRUNCA32X2F64S(A11,A11, sa-2),AE_TRUNCA32X2F64S(A11,A11, sa-2));
            AE_S8_0_X (z10,(ae_int8*)z,P  );
            AE_S8_0_X (z11,(ae_int8*)z,P+1);
            AE_S8_0_I (z01,(ae_int8*) z, 1);
            AE_S8_0_IP(z00,castxcc(ae_int8,z),2);
            pY0+=N;
            pX0-=N;
            pX1-=N;
        }
        if (P&1)
        {
            ae_int8x8 z00,z10;
            AE_MOVDX2(A00,A10,0,0);
            __Pragma("loop_count min=1")
            for (n=0; n<N; n++)
            {
                ae_int8x8 x0,x1,y0;
                ae_int8x8 x01,t01;
                AE_L8_IP(x0,castxcc(ae_int8,pX0),1);
                AE_L8_IP(x1,castxcc(ae_int8,pX1),1);
                AE_L8_IP(y0,castxcc(ae_int8,pY0),1);
                AE_DSEL8X8(x01,t01,x0,x1,dsel_rephhll);
                AE_MULAAAA2Q8(A00,A10,x01,y0);
            }
            z00=AE_ROUND8X4F32SASYM_L(AE_TRUNCA32X2F64S(A00,A00, sa-2),AE_TRUNCA32X2F64S(A00,A00, sa-2));
            z10=AE_ROUND8X4F32SASYM_L(AE_TRUNCA32X2F64S(A10,A10, sa-2),AE_TRUNCA32X2F64S(A10,A10, sa-2));
            AE_S8_0_X (z10,(ae_int8*)z,P  );
            AE_S8_0_IP(z00,castxcc(ae_int8,z),1);
            pX0-=N;
            pX1-=N;
        }
        z+=P;
        pX0=(const ae_int8*)XT_ADDX2(N,(uintptr_t)pX0);
        pX1=(const ae_int8*)XT_ADDX2(N,(uintptr_t)pX1);
    }
    if (M&1)
    {
        pY0=(const ae_int8*)&y[0];
        for (p=0; p<P; p++)
        {
            ae_int8x8 z00;
            A00=0;
            __Pragma("loop_count min=1")
            for (n=0; n<N; n++)
            {
                ae_int8x8 x0,y0;
                AE_L8_IP(x0,castxcc(ae_int8,pX0),1);
                AE_L8_IP(y0,castxcc(ae_int8,pY0),1);
                AE_MULAAAA2Q8(A00,t0,x0,y0);
            }
            z00=AE_ROUND8X4F32SASYM_L(AE_TRUNCA32X2F64S(A00,A00, sa-2),AE_TRUNCA32X2F64S(A00,A00, sa-2));
            AE_S8_0_IP(z00,castxcc(ae_int8,z),1);
            pX0-=N;
        }
    }
}
#endif

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

void mtx_mpyt8x8 (  void* pScr,
                     int8_t* restrict z,
               const int8_t* restrict x,
               const int8_t* restrict y,
               int M, int N, int P, int lsh )
{
#if USE_NN_EXTENSION_8X8==0
    ae_int8 * restrict pZ;
    int m,n,p;
    NASSERT(lsh >= -15 && lsh <= 15);
    NASSERT_ALIGN(pScr,HIFI_SIMD_WIDTH);
    (void)pScr;
    int8_t *w=(int8_t *)pScr;
    int8_t *v;
    int N0,P0;

    if (P<=0 || M<=0) return;
    if (N<=0)    /* exceptional situation */
    {
        ae_valign aZ;
        for (m=0; m<((M*P)&7); m++) *z++=0;
        aZ=AE_ZALIGN64();
        for (m=0; m<((M*P)>>3); m++) AE_SA8X8_IP(AE_MOVINT8X8_FROMINT32X2(AE_ZERO32()),aZ,castxcc(ae_int8x8,z));
        AE_SA64POS_FP(aZ,z);
        return ;
    }
    if (M*N*P<128 || M<2)  
    {
        mtx_mpyt8x8_small(z,x,y,M,N,P,lsh);
        return;
    }
    N0=((N+7)&~7);
    P0=((P+3)&~3);
    v=w+N0*P0;
#if 0
    for (p=0; p<P; p++)
    {
        for (n=0; n<N ; n++) w[p*N0+n]=y[p*N +n];
        for (   ; n<N0; n++) w[p*N0+n]=0;

    }
    for (   ; p<P0; p++)
    for (n=0; n<N0; n++) w[p*N0+n]=0;
#else
    pZ=(ae_int8*)w;
    for (p=0; p<((P0*N0)>>4); p++) AE_S32X2X2_IP(0,0,castxcc(ae_int32x4,pZ),sizeof(ae_int32x4));
    pZ=(ae_int8*)w;
    if (N<8)
    {
        for (p=0; p<P; p++)
        {
            ae_int8x8 x;
            __Pragma("loop_count min=1,max=7")
            for (n=0; n<N ; n++)
            {
                AE_L8_XP(x,castxcc(ae_int8,y),1);
                AE_S8_0_IP(x,pZ,1);
            }
            pZ+=N0-N;
        }
    }
    else
    {
        for (p=0; p<P; p++)
        {
            ae_int8x8 x;
            ae_valign aY;
            aY=AE_LA64_PP(y);
            for (n=0; n<(N>>3) ; n++)
            {
                AE_LA8X8_IP(x,aY,castxcc(ae_int8x8,y));
                AE_S8X8_IP(x,castxcc(ae_int8x8,pZ),sizeof(ae_int8x8));
            }
            __Pragma("loop_count max=7")
            for (n=0; n<(N&7) ; n++)
            {
                AE_L8_XP(x,castxcc(ae_int8,y),1);
                AE_S8_0_IP(x,pZ,1);
            }
            pZ+=N0-N;
        }
    }
#endif
    y=w;
    /* for very small N we simply copy input matrix x to Mx8 */
    if (N<8)
    {
        int8_t *newx=v+((M+3)&~3)*P0;
        mtx_mpy8x8_copysmallx(newx,x,M, N);
        x=newx;
        N=8; 
    }

    mtx_mpyt8x8_partiallyaligned (v,x,y,M, N,P, lsh );
    mtx_mpyt8x8_copyz( z,v,M,P);
#else
    int N0;
    int m,n,p;
    int8_t *yp;
    ae_valignx2 alX0, alX1, alX2, alX3;
    ae_valignx2 alY0;
    ae_valignx2 alZ0, alZ1, alZ2, alZ3;
    int8_t * restrict pZ0;
    int8_t * restrict pZ1;
    int8_t * restrict pZ2;
    int8_t * restrict pZ3;
    const int8_t * restrict pX0;
    const int8_t * restrict pX1;
    const int8_t * restrict pX2;
    const int8_t * restrict pX3;
    const int8_t * restrict pY0;

    ae_int32x2 Z00, Z01, Z10, Z11, Z20, Z21, Z30, Z31;
    ae_int16x4 t0, t1, t2, t3;
    ae_int8x8 x00, x01, x10, x11, x20, x21, x30, x31;
    ae_int8x8 y00, y01, y10, y11, y20, y21, y30, y31;
    ae_int8x8 z0, z1, z2, z3;
    ae_int8x8 zero;

    NASSERT(lsh >= -15 && lsh <= 15);
    NASSERT_ALIGN(pScr,HIFI_SIMD_WIDTH);

    if (M<=0 || P<=0) return;
    zero = AE_MOVINT8X8_FROMINT64(AE_ZERO());
    if (N<=0)    /* exceptional situation */
    {
        pZ0 = z;
        alZ0 = AE_ZALIGN128();
        for (m=0; m<((M*P)>>4); m++)  AE_SA8X8X2_IP(zero, zero, alZ0, castxcc(ae_int8x16,pZ0));
        AE_SAV8X8X2_XP(zero, zero, alZ0, castxcc(ae_int8x16,pZ0), (M*P)&15);
        AE_SA128POS_FP(alZ0, pZ0);
        return;
    }

    /* Copy matrix y to buffer with zero-padded columns and rows; *
     * number of columns is rounded to a next multiple of 16,     *
     * number of rows is rounded to a next multiple of 4          */
    yp = (int8_t *)pScr;
    N0 = (N+15)&~15;
    {
        int P0 = (P+3)&~3;
        pY0 = y;
        pZ0 = yp;
        alY0 = AE_LA128_PP(pY0);
        for (p = 0; p < P; p++)
        {
            for (n = 0; n < ((N-1)>>4); n++)
            {
                AE_LA8X8X2_IP(y00, y01, alY0, castxcc(ae_int8x16,pY0));
                AE_S8X8X2_IP(y00, y01, castxcc(ae_int8x16,pZ0), sizeof(ae_int8x16));
            }
            {
                AE_LAV8X8X2_XP(y00, y01, alY0, castxcc(ae_int8x16,pY0), ((N-1)&15)+1);
                AE_S8X8X2_IP(y00, y01, castxcc(ae_int8x16,pZ0), sizeof(ae_int8x16));
            }
        }
        for (p = 0; p < ((P0-P)*N0)>>4; p++)
        {
            AE_S8X8X2_IP(zero, zero, castxcc(ae_int8x16,pZ0), sizeof(ae_int8x16));
        }
    }

    pZ3 = z;
    for (m=0; m<(M>>2);m++)
    {
        pZ0 = pZ3;
        pZ1 = pZ0 + P;
        pZ2 = pZ1 + P;
        pZ3 = pZ2 + P;
        for (p=0; p < P; p+=4)
        {
            pX0 = x + m*N*4;
            pX1 = pX0 + N;
            pX2 = pX1 + N;
            pX3 = pX2 + N;
            pY0 = yp + p*N0;
            alX0 = AE_LA128_PP(pX0);
            alX1 = AE_LA128_PP(pX1);
            alX2 = AE_LA128_PP(pX2);
            alX3 = AE_LA128_PP(pX3);

            Z00 = Z01 = Z10 = Z11 = Z20 = Z21 = Z30 = Z31 = AE_ZERO32();

            for (n=0; n<((N-1)>>4); n++)
            {
                /* load x matrix, 4x16 values, 8-bit */
                AE_LA8X8X2_IP(x00, x01, alX0, castxcc(ae_int8x16,pX0));
                AE_LA8X8X2_IP(x10, x11, alX1, castxcc(ae_int8x16,pX1));
                AE_LA8X8X2_IP(x20, x21, alX2, castxcc(ae_int8x16,pX2));
                AE_LA8X8X2_IP(x30, x31, alX3, castxcc(ae_int8x16,pX3));

                /* load y matrix, 4x16 values, 8-bit */
                AE_L8X8X2_XP(y00, y01, castxcc(ae_int8x16,pY0), N0);
                AE_L8X8X2_XP(y10, y11, castxcc(ae_int8x16,pY0), N0);
                AE_L8X8X2_XP(y20, y21, castxcc(ae_int8x16,pY0), N0);
                AE_L8X8X2_XP(y30, y31, castxcc(ae_int8x16,pY0), -3*N0+16);

                /* make multiply, using 4-way 32-bit output octal mac, 8X8 bit */
                AE_MULA8Q8X8(Z00, Z01, y00, y10, y20, y30, x00);
                AE_MULA8Q8X8(Z10, Z11, y00, y10, y20, y30, x10);
                AE_MULA8Q8X8(Z20, Z21, y00, y10, y20, y30, x20);
                AE_MULA8Q8X8(Z30, Z31, y00, y10, y20, y30, x30);
                AE_MULA8Q8X8(Z00, Z01, y01, y11, y21, y31, x01);
                AE_MULA8Q8X8(Z10, Z11, y01, y11, y21, y31, x11);
                AE_MULA8Q8X8(Z20, Z21, y01, y11, y21, y31, x21);
                AE_MULA8Q8X8(Z30, Z31, y01, y11, y21, y31, x31);
            }
            /* Process last 1...16 multiplications of each output value */
            {
                /* load x matrix, 4x16 values, 8-bit */
                AE_LAV8X8X2_XP(x00, x01, alX0, castxcc(ae_int8x16,pX0), ((N-1)&15)+1);
                AE_LAV8X8X2_XP(x10, x11, alX1, castxcc(ae_int8x16,pX1), ((N-1)&15)+1);
                AE_LAV8X8X2_XP(x20, x21, alX2, castxcc(ae_int8x16,pX2), ((N-1)&15)+1);
                AE_LAV8X8X2_XP(x30, x31, alX3, castxcc(ae_int8x16,pX3), ((N-1)&15)+1);

                /* load y matrix, 4x16 values, 8-bit */
                AE_L8X8X2_XP(y00, y01, castxcc(ae_int8x16,pY0), N0);
                AE_L8X8X2_XP(y10, y11, castxcc(ae_int8x16,pY0), N0);
                AE_L8X8X2_XP(y20, y21, castxcc(ae_int8x16,pY0), N0);
                AE_L8X8X2_XP(y30, y31, castxcc(ae_int8x16,pY0), -3*N0+16);

                /* make multiply, using 4-way 32-bit output octal mac, 8X8 bit */
                AE_MULA8Q8X8(Z00, Z01, y00, y10, y20, y30, x00);
                AE_MULA8Q8X8(Z10, Z11, y00, y10, y20, y30, x10);
                AE_MULA8Q8X8(Z20, Z21, y00, y10, y20, y30, x20);
                AE_MULA8Q8X8(Z30, Z31, y00, y10, y20, y30, x30);
                AE_MULA8Q8X8(Z00, Z01, y01, y11, y21, y31, x01);
                AE_MULA8Q8X8(Z10, Z11, y01, y11, y21, y31, x11);
                AE_MULA8Q8X8(Z20, Z21, y01, y11, y21, y31, x21);
                AE_MULA8Q8X8(Z30, Z31, y01, y11, y21, y31, x31);
            }
            /* Q15 + lsh <- Q14 + 1 + lsh w/ truncation and saturation */
            t0 = AE_TRUNCA16X4F32S(Z00, Z01, 16+1+lsh);
            t1 = AE_TRUNCA16X4F32S(Z10, Z11, 16+1+lsh);
            t2 = AE_TRUNCA16X4F32S(Z20, Z21, 16+1+lsh);
            t3 = AE_TRUNCA16X4F32S(Z30, Z31, 16+1+lsh);
            /* Q7 + lsh <- Q15 + lsh -8 w/ rounding and saturation */
            z0 = AE_ROUND8X8F16SASYM(t0, t0);
            z1 = AE_ROUND8X8F16SASYM(t1, t1);
            z2 = AE_ROUND8X8F16SASYM(t2, t2);
            z3 = AE_ROUND8X8F16SASYM(t3, t3);
            alZ0 = alZ1 = alZ2 = alZ3 = AE_ZALIGN128();
            AE_SAV8X8X2_XP(z0, z0, alZ0, castxcc(ae_int8x16,pZ0), XT_MIN(4, P-p));
            AE_SAV8X8X2_XP(z1, z1, alZ1, castxcc(ae_int8x16,pZ1), XT_MIN(4, P-p));
            AE_SAV8X8X2_XP(z2, z2, alZ2, castxcc(ae_int8x16,pZ2), XT_MIN(4, P-p));
            AE_SAV8X8X2_XP(z3, z3, alZ3, castxcc(ae_int8x16,pZ3), XT_MIN(4, P-p));
            AE_SA128POS_FP(alZ0, pZ0);
            AE_SA128POS_FP(alZ1, pZ1);
            AE_SA128POS_FP(alZ2, pZ2);
            AE_SA128POS_FP(alZ3, pZ3);
        }
    }
    /* Process last M%4 rows */
    pZ0 = pZ3;
    for (m = (M&~3); m < M; m++)
    {
        for (p=0; p < P; p+=4)
        {
            pX0 = x + m*N;
            pY0 = yp + p*N0;
            alX0 = AE_LA128_PP(pX0);

            Z00 = Z01 = AE_ZERO32();

            for (n=0; n<((N-1)>>4); n++)
            {
                /* load x matrix, 1x16 values, 8-bit */
                AE_LA8X8X2_IP(x00, x01, alX0, castxcc(ae_int8x16,pX0));

                /* load y matrix, 4x16 values, 8-bit */
                AE_L8X8X2_XP(y00, y01, castxcc(ae_int8x16,pY0), N0);
                AE_L8X8X2_XP(y10, y11, castxcc(ae_int8x16,pY0), N0);
                AE_L8X8X2_XP(y20, y21, castxcc(ae_int8x16,pY0), N0);
                AE_L8X8X2_XP(y30, y31, castxcc(ae_int8x16,pY0), -3*N0+16);

                /* make multiply, using 4-way 32-bit output octal mac, 8X8 bit */
                AE_MULA8Q8X8(Z00, Z01, y00, y10, y20, y30, x00);
                AE_MULA8Q8X8(Z00, Z01, y01, y11, y21, y31, x01);
            }
            /* Process last 1...16 multiplications of each output value */
            {
                /* load x matrix, 1x16 values, 8-bit */
                AE_LAV8X8X2_XP(x00, x01, alX0, castxcc(ae_int8x16,pX0), ((N-1)&15)+1);

                /* load y matrix, 4x16 values, 8-bit */
                AE_L8X8X2_XP(y00, y01, castxcc(ae_int8x16,pY0), N0);
                AE_L8X8X2_XP(y10, y11, castxcc(ae_int8x16,pY0), N0);
                AE_L8X8X2_XP(y20, y21, castxcc(ae_int8x16,pY0), N0);
                AE_L8X8X2_XP(y30, y31, castxcc(ae_int8x16,pY0), -3*N0+16);

                /* make multiply, using 4-way 32-bit output octal mac, 8X8 bit */
                AE_MULA8Q8X8(Z00, Z01, y00, y10, y20, y30, x00);
                AE_MULA8Q8X8(Z00, Z01, y01, y11, y21, y31, x01);
            }
            /* Q15 + lsh <- Q14 + 1 + lsh w/ truncation and saturation */
            t0 = AE_TRUNCA16X4F32S(Z00, Z01, 16+1+lsh);
            /* Q7 + lsh <- Q15 + lsh -8 w/ rounding and saturation */
            z0 = AE_ROUND8X8F16SASYM(t0, t0);
            alZ0 = AE_ZALIGN128();
            AE_SAV8X8X2_XP(z0, z0, alZ0, castxcc(ae_int8x16,pZ0), XT_MIN(4, P-p));
            AE_SA128POS_FP(alZ0, pZ0);
        }
    }
#endif
}
size_t mtx_mpyt8x8_getScratchSize ( int M, int N, int P)
{
#if USE_NN_EXTENSION_8X8==0
    /* allocate matrices of size (N+7)&~7 x (P+3)&~3 and ((M+3)&~3) x (P+3)&~3 */
    size_t ysz,zsz,xsz;
    if (M<=0 || N<=0 || P<=0) return 0;
    int P0=((P+3)&~3);
    int N0=((N+7)&~7);
    int M0=((M+3)&~3);
    ysz=N0 * P0;
    zsz=M0 * P0;
    xsz=(N<8) ? M*8 : 0;
    return ysz+zsz+xsz;
#else
    if (M<=0 || N<=0 || P<=0) return 0;
    P = (P+3)&~3;
    N = (N+15)&~15;
    return P*N;
#endif
}
