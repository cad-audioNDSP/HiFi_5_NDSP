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
#include "mtx_mpy8x8_helper.h"

/* code optimized for HiFi5 core */

#if USE_NN_EXTENSION_8X8==0
static const uint64_t dsel_rephhll_tbl=0x40516273c8d9eafbULL;

/* special code for multiplication of small matrices without scratch */
static void mtx_mpy8x8_small(
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
    for (m=0; m<(M&~1);m+=2)
    {
        pY0=(const ae_int8*)&y[0];
        for (p=0; p<(P&~1); p+=2)
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
                y1=AE_L8_I ((const ae_int8*)pY0 ,1);
                AE_L8_XP(y0,castxcc(ae_int8,pY0),P);
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
            pY0-=(P*N-2);
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
                AE_L8_XP(y0,castxcc(ae_int8,pY0),P);
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
                AE_L8_XP(y0,castxcc(ae_int8,pY0),P);
                AE_MULAAAA2Q8(A00,t0,x0,y0);
            }
            z00=AE_ROUND8X4F32SASYM_L(AE_TRUNCA32X2F64S(A00,A00, sa-2),AE_TRUNCA32X2F64S(A00,A00, sa-2));
            AE_S8_0_IP(z00,castxcc(ae_int8,z),1);
            pX0-=N;
            pY0-=(P*N-1);
        }
    }
}
#else
static const uint8_t ALIGN(16) dsel_tbl[] =
{
    (15<<4)|13, (7<<4)|5, (14<<4)|12, (6<<4)|4, (11<<4)|9, (3<<4)|1, (10<<4)|8, (2<<4)|0
};

/* special code for multiplication of matrices with small P */
static void mtx_mpy8x8_smallP(
                int8_t* restrict z,
               const int8_t* restrict x,
               const int8_t* restrict y,
               int M, int N, int P, int lsh )
{
    int m, n;
    ae_valignx2 alX0, alX1, alX2, alX3, alY0, alY1, alY2, alY3;
    ae_valignx2 alZ0;
    int8_t * restrict pZ0;
    const int8_t * restrict pX0;
    const int8_t * restrict pX1;
    const int8_t * restrict pX2;
    const int8_t * restrict pX3;
    const int8_t * restrict pY0;
    const int8_t * restrict pY1;
    const int8_t * restrict pY2;
    const int8_t * restrict pY3;

    ae_int32x2 Z00, Z01, Z02, Z03,
               Z10, Z11, Z12, Z13,
               Z20, Z21, Z22, Z23,
               Z30, Z31, Z32, Z33;
    ae_int16x4 t0, t1, t2, t3;
    ae_int8x8 x00, x10, x20, x30;
    ae_int8x8 y00, y10, y20, y30;
    ae_int8x8 t00, t10, t20, t30;
    ae_int8x8 z0, z1;
    ae_int8x8 zero, dsel_idx, tmp;

    NASSERT(P < 9);

    WUR_AE_CBEGIN0((uintptr_t)y);
    WUR_AE_CEND0((uintptr_t)(y + N*P));
    zero = AE_MOVINT8X8_FROMINT64(AE_ZERO());
    dsel_idx = AE_L8X8_I((const ae_int8x8 *)dsel_tbl, 0);

    for (m=0; m<(M>>2); m++)
    {
        pZ0 = z + m*P*4;
        pX0 = x + m*N*4;
        pX1 = pX0 + N;
        pX2 = pX1 + N;
        pX3 = pX2 + N;
        pY0 = y;
        alX0 = AE_LA128_PP(pX0);
        alX1 = AE_LA128_PP(pX1);
        alX2 = AE_LA128_PP(pX2);
        alX3 = AE_LA128_PP(pX3);
        alY0 = AE_LA128_PP(pY0);

        Z00 = Z01 = Z02 = Z03 = Z10 = Z11 = Z12 = Z13 = AE_ZERO32();
        Z20 = Z21 = Z22 = Z23 = Z30 = Z31 = Z32 = Z33 = AE_ZERO32();

        for (n=0; n<((N-1)>>2); n++)
        {
            /* load x matrix, 4x4 values, 8-bit */
            alX0 = AE_LA128_PP(pX0);
            alX1 = AE_LA128_PP(pX1);
            alX2 = AE_LA128_PP(pX2);
            alX3 = AE_LA128_PP(pX3);
            AE_LAV8X8X2_XP(x00, tmp, alX0, castxcc(ae_int8x16,pX0), 4);
            AE_LAV8X8X2_XP(x10, tmp, alX1, castxcc(ae_int8x16,pX1), 4);
            AE_LAV8X8X2_XP(x20, tmp, alX2, castxcc(ae_int8x16,pX2), 4);
            AE_LAV8X8X2_XP(x30, tmp, alX3, castxcc(ae_int8x16,pX3), 4);
            x00 = AE_SEL8X8I(x00, x00, 1);/* AE_SELI_8B_INTERLEAVE_4_HH */
            x10 = AE_SEL8X8I(x10, x10, 1);/* AE_SELI_8B_INTERLEAVE_4_HH */
            x20 = AE_SEL8X8I(x20, x20, 1);/* AE_SELI_8B_INTERLEAVE_4_HH */
            x30 = AE_SEL8X8I(x30, x30, 1);/* AE_SELI_8B_INTERLEAVE_4_HH */

            /* load y matrix, 4x8 values, 8-bit */
            AE_LAV8X8X2_XP(y00, tmp, alY0, castxcc(ae_int8x16,pY0), P);
            AE_LAV8X8X2_XP(y10, tmp, alY0, castxcc(ae_int8x16,pY0), P);
            AE_LAV8X8X2_XP(y20, tmp, alY0, castxcc(ae_int8x16,pY0), P);
            AE_LAV8X8X2_XP(y30, tmp, alY0, castxcc(ae_int8x16,pY0), P);
            AE_DSEL8X8(t00, t10, y00, y20, dsel_idx);
            AE_DSEL8X8(t20, t30, y10, y30, dsel_idx);
            AE_DSEL8X8(y00, y10, t00, t20, dsel_idx);
            AE_DSEL8X8(y20, y30, t10, t30, dsel_idx);

            /* make multiply, using 8-way 32-bit output quad mac, 8X8 bit */
            AE_MULA4O8X8(Z00, Z01, Z02, Z03, y00, y10, y20, y30, x00);
            AE_MULA4O8X8(Z10, Z11, Z12, Z13, y00, y10, y20, y30, x10);
            AE_MULA4O8X8(Z20, Z21, Z22, Z23, y00, y10, y20, y30, x20);
            AE_MULA4O8X8(Z30, Z31, Z32, Z33, y00, y10, y20, y30, x30);
        }
        /* Process last 1...4 multiplications of each output value */
        {
            /* load x matrix, 4x4 values, 8-bit */
            alX0 = AE_LA128_PP(pX0);
            alX1 = AE_LA128_PP(pX1);
            alX2 = AE_LA128_PP(pX2);
            alX3 = AE_LA128_PP(pX3);
            AE_LAV8X8X2_XP(x00, tmp, alX0, castxcc(ae_int8x16,pX0), ((N-1)&3) + 1);
            AE_LAV8X8X2_XP(x10, tmp, alX1, castxcc(ae_int8x16,pX1), ((N-1)&3) + 1);
            AE_LAV8X8X2_XP(x20, tmp, alX2, castxcc(ae_int8x16,pX2), ((N-1)&3) + 1);
            AE_LAV8X8X2_XP(x30, tmp, alX3, castxcc(ae_int8x16,pX3), ((N-1)&3) + 1);
            x00 = AE_SEL8X8I(x00, x00, 1);/* AE_SELI_8B_INTERLEAVE_4_HH */
            x10 = AE_SEL8X8I(x10, x10, 1);/* AE_SELI_8B_INTERLEAVE_4_HH */
            x20 = AE_SEL8X8I(x20, x20, 1);/* AE_SELI_8B_INTERLEAVE_4_HH */
            x30 = AE_SEL8X8I(x30, x30, 1);/* AE_SELI_8B_INTERLEAVE_4_HH */

            /* load y matrix, 4x8 values, 8-bit */
            pY1 = pY0;  AE_ADDCIRC_XC(castxcc(ae_int64,pY1), P); alY1 = AE_LA128_PP(pY1);
            pY2 = pY1;  AE_ADDCIRC_XC(castxcc(ae_int64,pY2), P); alY2 = AE_LA128_PP(pY2);
            pY3 = pY2;  AE_ADDCIRC_XC(castxcc(ae_int64,pY3), P); alY3 = AE_LA128_PP(pY3);
            AE_LAV8X8X2_XP(y00, tmp, alY0, castxcc(ae_int8x16,pY0), P);
            AE_LAV8X8X2_XP(y10, tmp, alY1, castxcc(ae_int8x16,pY1), P);
            AE_LAV8X8X2_XP(y20, tmp, alY2, castxcc(ae_int8x16,pY2), P);
            AE_LAV8X8X2_XP(y30, tmp, alY3, castxcc(ae_int8x16,pY3), P);
            AE_DSEL8X8(t00, t10, y00, y20, dsel_idx);
            AE_DSEL8X8(t20, t30, y10, y30, dsel_idx);
            AE_DSEL8X8(y00, y10, t00, t20, dsel_idx);
            AE_DSEL8X8(y20, y30, t10, t30, dsel_idx);

            /* make multiply, using 8-way 32-bit output quad mac, 8X8 bit */
            AE_MULA4O8X8(Z00, Z01, Z02, Z03, y00, y10, y20, y30, x00);
            AE_MULA4O8X8(Z10, Z11, Z12, Z13, y00, y10, y20, y30, x10);
            AE_MULA4O8X8(Z20, Z21, Z22, Z23, y00, y10, y20, y30, x20);
            AE_MULA4O8X8(Z30, Z31, Z32, Z33, y00, y10, y20, y30, x30);
        }
        /* Q15 + lsh <- Q14 + 1 + lsh w/ truncation and saturation */
        t0 = AE_TRUNCA16X4F32S(Z00, Z01, 16+1+lsh);
        t1 = AE_TRUNCA16X4F32S(Z02, Z03, 16+1+lsh);
        t2 = AE_TRUNCA16X4F32S(Z10, Z11, 16+1+lsh);
        t3 = AE_TRUNCA16X4F32S(Z12, Z13, 16+1+lsh);
        /* Q7 + lsh <- Q15 + lsh -8 w/ rounding and saturation */
        z0 = AE_ROUND8X8F16SASYM(t0, t1);
        z1 = AE_ROUND8X8F16SASYM(t2, t3);
        alZ0 = AE_ZALIGN128();
        AE_SAV8X8X2_XP(z0, z0, alZ0, castxcc(ae_int8x16,pZ0), P);
        AE_SAV8X8X2_XP(z1, z1, alZ0, castxcc(ae_int8x16,pZ0), P);
        AE_SA128POS_FP(alZ0, pZ0);

        /* Save next 2 rows */
        t0 = AE_TRUNCA16X4F32S(Z20, Z21, 16+1+lsh);
        t1 = AE_TRUNCA16X4F32S(Z22, Z23, 16+1+lsh);
        t2 = AE_TRUNCA16X4F32S(Z30, Z31, 16+1+lsh);
        t3 = AE_TRUNCA16X4F32S(Z32, Z33, 16+1+lsh);
        z0 = AE_ROUND8X8F16SASYM(t0, t1);
        z1 = AE_ROUND8X8F16SASYM(t2, t3);
        alZ0 = AE_ZALIGN128();
        AE_SAV8X8X2_XP(z0, z0, alZ0, castxcc(ae_int8x16,pZ0), P);
        AE_SAV8X8X2_XP(z1, z1, alZ0, castxcc(ae_int8x16,pZ0), P);
        AE_SA128POS_FP(alZ0, pZ0);
    }
    /* Process last M%4 rows */
    for (m = (M&~3); m < M; m++)
    {
        pZ0 = z + m*P;
        pX0 = x + m*N;
        pY0 = y;
        alX0 = AE_LA128_PP(pX0);
        alY0 = AE_LA128_PP(pY0);

        Z00 = Z01 = Z02 = Z03 = AE_ZERO32();

        for (n=0; n<((N-1)>>2); n++)
        {
            /* load x matrix, 1x4 values, 8-bit */
            AE_LAV8X8X2_XP(x00, tmp, alX0, castxcc(ae_int8x16,pX0), 4);
            x00 = AE_SEL8X8I(x00, x00, 1);/* AE_SELI_8B_INTERLEAVE_4_HH */

            /* load y matrix, 4x8 values, 8-bit */
            AE_LAV8X8X2_XP(y00, tmp, alY0, castxcc(ae_int8x16,pY0), P);
            AE_LAV8X8X2_XP(y10, tmp, alY0, castxcc(ae_int8x16,pY0), P);
            AE_LAV8X8X2_XP(y20, tmp, alY0, castxcc(ae_int8x16,pY0), P);
            AE_LAV8X8X2_XP(y30, tmp, alY0, castxcc(ae_int8x16,pY0), P);
            AE_DSEL8X8(t00, t10, y00, y20, dsel_idx);
            AE_DSEL8X8(t20, t30, y10, y30, dsel_idx);
            AE_DSEL8X8(y00, y10, t00, t20, dsel_idx);
            AE_DSEL8X8(y20, y30, t10, t30, dsel_idx);

            /* make multiply, using 8-way 32-bit output quad mac, 8X8 bit */
            AE_MULA4O8X8(Z00, Z01, Z02, Z03, y00, y10, y20, y30, x00);
        }
        /* Process last 1...4 multiplications of each output value */
        {
            /* load x matrix, 1x4 values, 8-bit */
            AE_LAV8X8X2_XP(x00, tmp, alX0, castxcc(ae_int8x16,pX0), ((N-1)&3) + 1);
            x00 = AE_SEL8X8I(x00, x00, 1);/* AE_SELI_8B_INTERLEAVE_4_HH */

            /* load y matrix, 4x8 values, 8-bit */
            /* apply circular addressing to prevent out-of-bound memory accesses */
            pY1 = pY0;  AE_ADDCIRC_XC(castxcc(ae_int64,pY1), P); alY1 = AE_LA128_PP(pY1);
            pY2 = pY1;  AE_ADDCIRC_XC(castxcc(ae_int64,pY2), P); alY2 = AE_LA128_PP(pY2);
            pY3 = pY2;  AE_ADDCIRC_XC(castxcc(ae_int64,pY3), P); alY3 = AE_LA128_PP(pY3);
            AE_LAV8X8X2_XP(y00, tmp, alY0, castxcc(ae_int8x16,pY0), P);
            AE_LAV8X8X2_XP(y10, tmp, alY1, castxcc(ae_int8x16,pY1), P);
            AE_LAV8X8X2_XP(y20, tmp, alY2, castxcc(ae_int8x16,pY2), P);
            AE_LAV8X8X2_XP(y30, tmp, alY3, castxcc(ae_int8x16,pY3), P);
            AE_DSEL8X8(t00, t10, y00, y20, dsel_idx);
            AE_DSEL8X8(t20, t30, y10, y30, dsel_idx);
            AE_DSEL8X8(y00, y10, t00, t20, dsel_idx);
            AE_DSEL8X8(y20, y30, t10, t30, dsel_idx);

            /* make multiply, using 8-way 32-bit output quad mac, 8X8 bit */
            AE_MULA4O8X8(Z00, Z01, Z02, Z03, y00, y10, y20, y30, x00);
        }
        /* Q15 + lsh <- Q14 + 1 + lsh w/ truncation and saturation */
        t0 = AE_TRUNCA16X4F32S(Z00, Z01, 16+1+lsh);
        t1 = AE_TRUNCA16X4F32S(Z02, Z03, 16+1+lsh);
        /* Q7 + lsh <- Q15 + lsh -8 w/ rounding and saturation */
        z0 = AE_ROUND8X8F16SASYM(t0, t1);
        alZ0 = AE_ZALIGN128();
        AE_SAV8X8X2_XP(z0, z0, alZ0, castxcc(ae_int8x16,pZ0), P);
        AE_SA128POS_FP(alZ0, pZ0);
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
void mtx_mpy8x8 (  void* pScr,
                     int8_t* restrict z,
               const int8_t* restrict x,
               const int8_t* restrict y,
               int M, int N, int P, int lsh )
{
#if USE_NN_EXTENSION_8X8==0
    ae_int8   * restrict pZ;

    int N0,P0;
    int8_t *w=(int8_t *)pScr;
    int8_t *v;
    int m,n,p;
    ae_int8x8 dsel_rephhll;

    NASSERT(lsh >= -15 && lsh <= 15);
    NASSERT_ALIGN(pScr,HIFI_SIMD_WIDTH);
    (void)pScr;
    if (M<=0 || P<=0) return ;
    if (N<=0)    /* exceptional situation */
    {
        ae_valign aZ;
        if (M<=0 || P<=0) return;
        for (m=0; m<((M*P)&7); m++) *z++=0;
        aZ=AE_ZALIGN64();
        for (m=0; m<((M*P)>>3); m++) AE_SA8X8_IP(AE_MOVINT8X8_FROMINT32X2(AE_ZERO32()),aZ,castxcc(ae_int8x8,z));
        AE_SA64POS_FP(aZ,z);
        return ;
    }

    /* if matrix sizes are very small it is better to multiply them element-wise */
    if (M*N*P<128 || M<4)  
    {
        return mtx_mpy8x8_small(z,x,y,M,N,P,lsh);
    }

    /* copy y to the scratch with padding and transposing */
    dsel_rephhll=AE_L8X8_I((const ae_int8x8*)&dsel_rephhll_tbl,0);
    N0=((N+7)&~7);
    P0=((P+3)&~3);
    v=w+N0*P0;
#if 0
    for (p=0; p<P; p++)
    {
        for (n=0; n<N ; n++) w[p*N0+n]=y[n*P+p];
        for (   ; n<N0; n++) w[p*N0+n]=0;

    }
    for (   ; p<P0; p++)
    for (n=0; n<N0; n++) w[p*N0+n]=0;
#else
    pZ=(ae_int8*)w;
    for (p=0; p<((P0*N0)>>4); p++) AE_S32X2X2_IP(0,0,castxcc(ae_int32x4,pZ),sizeof(ae_int32x4));
    pZ=(ae_int8*)w;
    for (p=0; p<P; p++)
    {
        ae_int8x8 x;
        for (n=0; n<N ; n++)
        {
            AE_L8_XP(x,castxcc(ae_int8,y),P);
            AE_S8_0_IP(x,pZ,1);
        }
        y-=P*N-1;
        pZ+=N0-N;
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
    int N0, P0;
    int8_t *ytp;
    int m,n,p;
    ae_valignx2 alX0, alX1, alX2, alX3;
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
    const int8_t * restrict pY1;
    const int8_t * restrict pY2;
    const int8_t * restrict pY3;

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
    pZ0 = z;
    if (N<=0)    /* exceptional situation */
    {
        alZ0 = AE_ZALIGN128();
        for (m=0; m<((M*P)>>4); m++)  AE_SA8X8X2_IP(zero, zero, alZ0, castxcc(ae_int8x16,pZ0));
        AE_SAV8X8X2_XP(zero, zero, alZ0, castxcc(ae_int8x16,pZ0), (M*P)&15);
        AE_SA128POS_FP(alZ0, pZ0);
        return;
    }

    /* Call specialized function when number of columns P is small */
    if (P <= 8)
    {
        return mtx_mpy8x8_smallP(z, x, y, M, N, P, lsh);
    }

    P0 = (P+15) & ~15;
    N0 = (N+15) & ~15;
    /* Copy matrix y to buffer with zero-padded columns and rows and transposed; *
     * number of columns is rounded to a next multiple of 16,                    *
     * number of rows is rounded to a next multiple of 16                        */
    WUR_AE_CBEGIN0((uintptr_t)(y));
    WUR_AE_CEND0((uintptr_t)(y + N*P));
    ytp = (int8_t *)pScr;
    {
        ae_valignx2 alY0, alY1, alY2, alY3;
        ae_int8x8 dsel_idx;
        dsel_idx = AE_L8X8_I((const ae_int8x8 *)dsel_tbl, 0);
        /* Clear last columns of padded matrix */
        {
            pZ0 = ytp + N0-16;
            __Pragma("loop_count min=1");
            for (p = 0; p < P0; p++)
            {
                AE_S8X8X2_XP(zero, zero, castxcc(ae_int8x16,pZ0), N0);
            }
        }
        for (n = 0; n < N; n+=4)
        {
            pY0 = y + n*P;
            pY1 = pY0;  AE_ADDCIRC_XC(castxcc(ae_int64,pY1), P);
            pY2 = pY1;  AE_ADDCIRC_XC(castxcc(ae_int64,pY2), P);
            pY3 = pY2;  AE_ADDCIRC_XC(castxcc(ae_int64,pY3), P);
            pZ0 = ytp + n;
            alY0 = AE_LA128_PP(pY0);
            alY1 = AE_LA128_PP(pY1);
            alY2 = AE_LA128_PP(pY2);
            alY3 = AE_LA128_PP(pY3);
            for (p = 0; p < P; p+=16)
            {
                AE_LAV8X8X2_XP(y00, y01, alY0, castxcc(ae_int8x16,pY0), P-p);
                AE_LAV8X8X2_XP(y10, y11, alY1, castxcc(ae_int8x16,pY1), P-p);
                AE_LAV8X8X2_XP(y20, y21, alY2, castxcc(ae_int8x16,pY2), P-p);
                AE_LAV8X8X2_XP(y30, y31, alY3, castxcc(ae_int8x16,pY3), P-p);
                /* transpose 4x16 to 16x4 */
                AE_DSEL8X8(x00, x10, y30, y10, dsel_idx);
                AE_DSEL8X8(x20, x30, y20, y00, dsel_idx);
                AE_DSEL8X8(y00, y10, x00, x20, dsel_idx);
                AE_DSEL8X8(y20, y30, x10, x30, dsel_idx);
                AE_DSEL8X8(x01, x11, y31, y11, dsel_idx);
                AE_DSEL8X8(x21, x31, y21, y01, dsel_idx);
                AE_DSEL8X8(y01, y11, x01, x21, dsel_idx);
                AE_DSEL8X8(y21, y31, x11, x31, dsel_idx);
                /* save result */
                AE_S32_H_XP(AE_MOVINT32X2_FROMINT8X8(y00), castxcc(ae_int32,pZ0), N0);
                AE_S32_H_XP(AE_MOVINT32X2_FROMINT8X8(y10), castxcc(ae_int32,pZ0), N0);
                AE_S32_H_XP(AE_MOVINT32X2_FROMINT8X8(y20), castxcc(ae_int32,pZ0), N0);
                AE_S32_H_XP(AE_MOVINT32X2_FROMINT8X8(y30), castxcc(ae_int32,pZ0), N0);
                AE_S32_L_XP(AE_MOVINT32X2_FROMINT8X8(y00), castxcc(ae_int32,pZ0), N0);
                AE_S32_L_XP(AE_MOVINT32X2_FROMINT8X8(y10), castxcc(ae_int32,pZ0), N0);
                AE_S32_L_XP(AE_MOVINT32X2_FROMINT8X8(y20), castxcc(ae_int32,pZ0), N0);
                AE_S32_L_XP(AE_MOVINT32X2_FROMINT8X8(y30), castxcc(ae_int32,pZ0), N0);
                AE_S32_H_XP(AE_MOVINT32X2_FROMINT8X8(y01), castxcc(ae_int32,pZ0), N0);
                AE_S32_H_XP(AE_MOVINT32X2_FROMINT8X8(y11), castxcc(ae_int32,pZ0), N0);
                AE_S32_H_XP(AE_MOVINT32X2_FROMINT8X8(y21), castxcc(ae_int32,pZ0), N0);
                AE_S32_H_XP(AE_MOVINT32X2_FROMINT8X8(y31), castxcc(ae_int32,pZ0), N0);
                AE_S32_L_XP(AE_MOVINT32X2_FROMINT8X8(y01), castxcc(ae_int32,pZ0), N0);
                AE_S32_L_XP(AE_MOVINT32X2_FROMINT8X8(y11), castxcc(ae_int32,pZ0), N0);
                AE_S32_L_XP(AE_MOVINT32X2_FROMINT8X8(y21), castxcc(ae_int32,pZ0), N0);
                AE_S32_L_XP(AE_MOVINT32X2_FROMINT8X8(y31), castxcc(ae_int32,pZ0), N0);
            }
        }
    }

    pZ3 = z;
    for (m=0; m<(M>>2);m++)
    {
        pZ0 = pZ3;
        pZ1 = pZ0 + P;
        pZ2 = pZ1 + P;
        pZ3 = pZ2 + P;
        pY3 = ytp;
        for (p=0; p < P; p+=4)
        {
            pX0 = x + m*N*4;
            pX1 = pX0 + N;
            pX2 = pX1 + N;
            pX3 = pX2 + N;
            pY0 = pY3;
            pY1 = pY0 + N0;
            pY2 = pY1 + N0;
            pY3 = pY2 + N0;
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
                AE_L8X8X2_IP(y00, y01, castxcc(ae_int8x16,pY0), sizeof(ae_int8x16));
                AE_L8X8X2_IP(y10, y11, castxcc(ae_int8x16,pY1), sizeof(ae_int8x16));
                AE_L8X8X2_IP(y20, y21, castxcc(ae_int8x16,pY2), sizeof(ae_int8x16));
                AE_L8X8X2_IP(y30, y31, castxcc(ae_int8x16,pY3), sizeof(ae_int8x16));

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
                AE_LAV8X8X2_XP(x00, x01, alX0, castxcc(ae_int8x16,pX0), ((N-1)&15) + 1);
                AE_LAV8X8X2_XP(x10, x11, alX1, castxcc(ae_int8x16,pX1), ((N-1)&15) + 1);
                AE_LAV8X8X2_XP(x20, x21, alX2, castxcc(ae_int8x16,pX2), ((N-1)&15) + 1);
                AE_LAV8X8X2_XP(x30, x31, alX3, castxcc(ae_int8x16,pX3), ((N-1)&15) + 1);

                /* load y matrix, 4x16 values, 8-bit */
                AE_L8X8X2_IP(y00, y01, castxcc(ae_int8x16,pY0), sizeof(ae_int8x16));
                AE_L8X8X2_IP(y10, y11, castxcc(ae_int8x16,pY1), sizeof(ae_int8x16));
                AE_L8X8X2_IP(y20, y21, castxcc(ae_int8x16,pY2), sizeof(ae_int8x16));
                AE_L8X8X2_IP(y30, y31, castxcc(ae_int8x16,pY3), sizeof(ae_int8x16));

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
        alZ0 = AE_ZALIGN128();
        pY3 = ytp;
        for (p=0; p < P; p+=4)
        {
            pX0 = x + m*N;
            pY0 = pY3;
            pY1 = pY0 + N0;
            pY2 = pY1 + N0;
            pY3 = pY2 + N0;
            alX0 = AE_LA128_PP(pX0);

            Z00 = Z01 = AE_ZERO32();

            for (n=0; n<((N-1)>>4); n++)
            {
                /* load x matrix, 1x16 values, 8-bit */
                AE_LA8X8X2_IP(x00, x01, alX0, castxcc(ae_int8x16,pX0));

                /* load y matrix, 4x16 values, 8-bit */
                AE_L8X8X2_IP(y00, y01, castxcc(ae_int8x16,pY0), sizeof(ae_int8x16));
                AE_L8X8X2_IP(y10, y11, castxcc(ae_int8x16,pY1), sizeof(ae_int8x16));
                AE_L8X8X2_IP(y20, y21, castxcc(ae_int8x16,pY2), sizeof(ae_int8x16));
                AE_L8X8X2_IP(y30, y31, castxcc(ae_int8x16,pY3), sizeof(ae_int8x16));

                /* make multiply, using 4-way 32-bit output octal mac, 8X8 bit */
                AE_MULA8Q8X8(Z00, Z01, y00, y10, y20, y30, x00);
                AE_MULA8Q8X8(Z00, Z01, y01, y11, y21, y31, x01);
            }
            /* Process last 1...16 multiplications of each output value */
            {
                /* load x matrix, 1x16 values, 8-bit */
                AE_LAV8X8X2_XP(x00, x01, alX0, castxcc(ae_int8x16,pX0), ((N-1)&15) + 1);

                /* load y matrix, 4x16 values, 8-bit */
                AE_L8X8X2_IP(y00, y01, castxcc(ae_int8x16,pY0), sizeof(ae_int8x16));
                AE_L8X8X2_IP(y10, y11, castxcc(ae_int8x16,pY1), sizeof(ae_int8x16));
                AE_L8X8X2_IP(y20, y21, castxcc(ae_int8x16,pY2), sizeof(ae_int8x16));
                AE_L8X8X2_IP(y30, y31, castxcc(ae_int8x16,pY3), sizeof(ae_int8x16));

                /* make multiply, using 4-way 32-bit output octal mac, 8X8 bit */
                AE_MULA8Q8X8(Z00, Z01, y00, y10, y20, y30, x00);
                AE_MULA8Q8X8(Z00, Z01, y01, y11, y21, y31, x01);
            }
            /* Q15 + lsh <- Q14 + 1 + lsh w/ truncation and saturation */
            t0 = AE_TRUNCA16X4F32S(Z00, Z01, 16+1+lsh);
            /* Q7 + lsh <- Q15 + lsh -8 w/ rounding and saturation */
            z0 = AE_ROUND8X8F16SASYM(t0, t0);
            AE_SAV8X8X2_XP(z0, z0, alZ0, castxcc(ae_int8x16,pZ0), XT_MIN(4, P-p));
        }
        AE_SA128POS_FP(alZ0, pZ0);
    }
#endif
}

size_t mtx_mpy8x8_getScratchSize ( int M, int N, int P)
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
    if (M<=0 || N<=0 || P<=8) return 0;
    P = (P+15)&~15;
    N = (N+15)&~15;
    return P*N;
#endif
}
