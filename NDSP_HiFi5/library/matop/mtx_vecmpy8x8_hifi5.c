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
void mtx_vecmpy8x8   ( int8_t* restrict z,
                 const int8_t* restrict x,
                 const int8_t* restrict y,
                 int M, int N, int lsh)
#if USE_NN_EXTENSION_8X8==0
{
    static const int64_t seltbl[]={
        0x2345670123456701ULL,
        0x3456701234567012ULL,
        0x4567012345670123ULL};
    ae_int8x8 sel_rotate1,sel_rotate2,sel_rotate3;
    const int8_t* restrict pX0;
    const int8_t* restrict pX1;
    const int8_t* restrict pX2;
    const int8_t* restrict pX3;
    const int8_t* restrict pY;
          ae_int8* restrict pZ=(ae_int8* )z;
    ae_valign aX0,aX1,aX2,aX3,aY;
    ae_int8x8 r;
    ae_int64 A0,A1,A2,A3;
    ae_int64 B0,B1,B2,B3;
    ae_int8x8 x0,x1,x2,x3,y0;
    int m,n;
    NASSERT(lsh >= -15 && lsh <= 15);
    if (M<=0) return;
    if (N<=0)    /* exceptional situation */
    {
        for (m=0; m<M; m++) z[m]=0;
        return ;
    }
    sel_rotate1=AE_MOVINT8X8_FROMINT64(AE_L64_I((const ae_int64*)seltbl,0*8));
    sel_rotate2=AE_MOVINT8X8_FROMINT64(AE_L64_I((const ae_int64*)seltbl,1*8));
    sel_rotate3=AE_MOVINT8X8_FROMINT64(AE_L64_I((const ae_int64*)seltbl,2*8));
    int sa=49+lsh;
    pX0=x+0*N;
    pX1=x+1*N;
    pX2=x+2*N;
    pX3=x+3*N;
    for (m=0; m<(M&~3);m+=4)
    {
        pY=y;
        AE_MOVDX2(A0,B0,0,0);
        AE_MOVDX2(A1,B1,0,0);
        AE_MOVDX2(A2,B2,0,0);
        AE_MOVDX2(A3,B3,0,0);
        aX0=AE_LA64_PP(pX0); 
        aX1=AE_LA64_PP(pX1); 
        for (n=0; n<(N>>3); n++)
        {
            AE_LA8X8_IP(x0,aX0,castxcc(ae_int8x8,pX0));
            AE_LA8X8_IP(x1,aX1,castxcc(ae_int8x8,pX1));
            aX2=AE_LA64_PP(pX2); AE_LA8X8_IP(x2,aX2,castxcc(ae_int8x8,pX2));
            aX3=AE_LA64_PP(pX3); AE_LA8X8_IP(x3,aX3,castxcc(ae_int8x8,pX3));
            aY =AE_LA64_PP(pY ); AE_LA8X8_IP(y0,aY ,castxcc(ae_int8x8,pY ));
            AE_MULAAAA2Q8(A0,B0,x0,y0);
            AE_MULAAAA2Q8(A1,B1,x1,y0);
            AE_MULAAAA2Q8(A2,B2,x2,y0);
            AE_MULAAAA2Q8(A3,B3,x3,y0);
        }
        for(n=0;n<(N&7);n++)
        {
            ae_int8x8 t,z=AE_MOVINT8X8_FROMINT32X2(AE_ZERO32());
            AE_L8_IP(x0,castxcc(ae_int8,pX0),1);
            AE_L8_IP(x1,castxcc(ae_int8,pX1),1);
            AE_L8_IP(x2,castxcc(ae_int8,pX2),1);
            AE_L8_IP(x3,castxcc(ae_int8,pX3),1);
            AE_L8_IP(y0,castxcc(ae_int8,pY ),1);
            AE_MOVT8X16_L(t, y0, z, y0, 1);
            AE_MULAAAA2Q8(A0,B0,x0,y0);
            AE_MULAAAA2Q8(A1,B1,x1,y0);
            AE_MULAAAA2Q8(A2,B2,x2,y0);
            AE_MULAAAA2Q8(A3,B3,x3,y0);
        }
        r=AE_ROUND8X4F32SASYM_L(AE_TRUNCA32X2F64S(A3+B3, A2+B2, sa),AE_TRUNCA32X2F64S(A1+B1, A0+B0, sa));
        AE_S8_0_I(AE_SEL8X8(r,r,sel_rotate1),pZ,1);
        AE_S8_0_I(AE_SEL8X8(r,r,sel_rotate2),pZ,2);
        AE_S8_0_I(AE_SEL8X8(r,r,sel_rotate3),pZ,3);
        AE_S8_0_IP(r,pZ,4);
        pX0+=3*N;
        pX1+=3*N;
        pX2+=3*N;
        pX3+=3*N;
    }
    if (M&2)
    {
        pY=y;
        AE_MOVDX2(A0,B0,0,0);
        AE_MOVDX2(A1,B1,0,0);
        aX0=AE_LA64_PP(pX0); 
        aX1=AE_LA64_PP(pX1); 
        aY =AE_LA64_PP(pY ); 
        for (n=0; n<(N>>3); n++)
        {
            AE_LA8X8_IP(x0,aX0,castxcc(ae_int8x8,pX0));
            AE_LA8X8_IP(x1,aX1,castxcc(ae_int8x8,pX1));
            AE_LA8X8_IP(y0,aY ,castxcc(ae_int8x8,pY ));
            AE_MULAAAA2Q8(A0,B0,x0,y0);
            AE_MULAAAA2Q8(A1,B1,x1,y0);
        }
        for(n=0;n<(N&7);n++)
        {
            ae_int8x8 t,z=AE_MOVINT8X8_FROMINT32X2(AE_ZERO32());
            AE_L8_IP(x0,castxcc(ae_int8,castxcc(ae_int8x8,pX0)),1);
            AE_L8_IP(x1,castxcc(ae_int8,castxcc(ae_int8x8,pX1)),1);
            AE_L8_IP(y0,castxcc(ae_int8,castxcc(ae_int8x8,pY )),1);
            AE_MOVT8X16_L(t, y0, z, y0, 1);
            AE_MULAAAA2Q8(A0,B0,x0,y0);
            AE_MULAAAA2Q8(A1,B1,x1,y0);
        }
        r=AE_ROUND8X4F32SASYM_L(0,AE_TRUNCA32X2F64S(A1+B1, A0+B0, sa));
        AE_S8_0_I(AE_SEL8X8(r,r,sel_rotate1),pZ,1);
        AE_S8_0_IP(r,pZ,2);
        pX0+=1*N;
    }
    if (M&1)
    {
        pY=y;
        AE_MOVDX2(A0,B0,0,0);
        aX0=AE_LA64_PP(pX0); 
        aY =AE_LA64_PP(pY ); 
        for (n=0; n<(N>>3); n++)
        {
            AE_LA8X8_IP(x0,aX0,castxcc(ae_int8x8,pX0));
            AE_LA8X8_IP(y0,aY ,castxcc(ae_int8x8,pY ));
            AE_MULAAAA2Q8(A0,B0,x0,y0);
        }
        for(n=0;n<(N&7);n++)
        {
            ae_int8x8 t,z=AE_MOVINT8X8_FROMINT32X2(AE_ZERO32());
            AE_L8_IP(x0,castxcc(ae_int8,pX0),1);
            AE_L8_IP(y0,castxcc(ae_int8,pY ),1);
            AE_MOVT8X16_L(t, y0, z, y0, 1);
            AE_MULAAAA2Q8(A0,B0,x0,y0);
        }
        r=AE_ROUND8X4F32SASYM_L(0,AE_TRUNCA32X2F64S(0, A0+B0, sa));
        AE_S8_0_IP(r,pZ,1);
    }
}
#elif 0
{
    static const int64_t seltbl[]={
        0x2345670123456701ULL};
    ae_int8x8 sel_rotate1;
    const int8_t* restrict pX0;
    const int8_t* restrict pX1;
    const int8_t* restrict pY;
          ae_int8* restrict pZ=(ae_int8* )z;
    ae_valign aX0,aX1,aY;
    ae_int8x8 r;
    ae_int64 A0,A1;
    ae_int64 B0,B1;
    ae_int8x8 x0,x1,y0;
    int m,n;
    NASSERT(lsh >= -15 && lsh <= 15);
    if (M<=0) return;
    if (N<=0)    /* exceptional situation */
    {
        for (m=0; m<M; m++) z[m]=0;
        return ;
    }
    sel_rotate1=AE_MOVINT8X8_FROMINT64(AE_L64_I((const ae_int64*)seltbl,0*8));
    int sa=49+lsh;
    pX0=x+0*N;
    pX1=x+1*N;
    for (m=0; m<(M&~1);m+=2)
    {
        pY=y;
        AE_MOVDX2(A0,B0,0,0);
        AE_MOVDX2(A1,B1,0,0);
        aX0=AE_LA64_PP(pX0); 
        aX1=AE_LA64_PP(pX1); 
        for (n=0; n<(N>>3); n++)
        {
            AE_LA8X8_IP(x0,aX0,castxcc(ae_int8x8,pX0));
            AE_LA8X8_IP(x1,aX1,castxcc(ae_int8x8,pX1));
            aY =AE_LA64_PP(pY ); AE_LA8X8_IP(y0,aY ,castxcc(ae_int8x8,pY ));
            AE_MULAAAA2Q8(A0,B0,x0,y0);
            AE_MULAAAA2Q8(A1,B1,x1,y0);
        }
        for(n=0;n<(N&7);n++)
        {
            ae_int8x8 t,z=AE_MOVINT8X8_FROMINT32X2(AE_ZERO32());
            AE_L8_IP(x0,castxcc(ae_int8,pX0),1);
            AE_L8_IP(x1,castxcc(ae_int8,pX1),1);
            AE_L8_IP(y0,castxcc(ae_int8,pY ),1);
            AE_MOVT8X16_L(t, y0, z, y0, 1);
            AE_MULAAAA2Q8(A0,B0,x0,y0);
            AE_MULAAAA2Q8(A1,B1,x1,y0);
        }
        r=AE_ROUND8X4F32SASYM_L(0,AE_TRUNCA32X2F64S(A1+B1, A0+B0, sa));
        AE_S8_0_I(AE_SEL8X8(r,r,sel_rotate1),pZ,1);
        AE_S8_0_IP(r,pZ,2);
        pX0+=1*N;
        pX1+=1*N;
    }
    if (M&1)
    {
        pY=y;
        AE_MOVDX2(A0,B0,0,0);
        aX0=AE_LA64_PP(pX0); 
        aY =AE_LA64_PP(pY ); 
        for (n=0; n<(N>>3); n++)
        {
            AE_LA8X8_IP(x0,aX0,castxcc(ae_int8x8,pX0));
            AE_LA8X8_IP(y0,aY ,castxcc(ae_int8x8,pY ));
            AE_MULAAAA2Q8(A0,B0,x0,y0);
        }
        for(n=0;n<(N&7);n++)
        {
            ae_int8x8 t,z=AE_MOVINT8X8_FROMINT32X2(AE_ZERO32());
            AE_L8_IP(x0,castxcc(ae_int8,pX0),1);
            AE_L8_IP(y0,castxcc(ae_int8,pY ),1);
            AE_MOVT8X16_L(t, y0, z, y0, 1);
            AE_MULAAAA2Q8(A0,B0,x0,y0);
        }
        r=AE_ROUND8X4F32SASYM_L(0,AE_TRUNCA32X2F64S(0, A0+B0, sa));
        AE_S8_0_IP(r,pZ,1);
    }
}
#else
{
    int m,n;
    ae_valignx2 alX0, alX1, alX2, alX3;
    ae_valignx2 alY0;
    ae_valignx2 alZ0;
    int8_t * restrict pZ0;
    const int8_t * restrict pX0;
    const int8_t * restrict pX1;
    const int8_t * restrict pX2;
    const int8_t * restrict pX3;
    const int8_t * restrict pY0;

    ae_int32x2 Z00, Z01;
    ae_int16x4 t0;
    ae_int8x8 x00, x01, x10, x11, x20, x21, x30, x31;
    ae_int8x8 y00, y01;
    ae_int8x8 z0;
    ae_int8x8 zero;

    NASSERT(lsh >= -15 && lsh <= 15);

    if (M<=0) return;
    zero = AE_MOVINT8X8_FROMINT64(AE_ZERO());
    if (N<=0)    /* exceptional situation */
    {
        pZ0 = z;
        alZ0 = AE_ZALIGN128();
        for (m=0; m<(M>>4); m++)  AE_SA8X8X2_IP(zero, zero, alZ0, castxcc(ae_int8x16,pZ0));
        AE_SAV8X8X2_XP(zero, zero, alZ0, castxcc(ae_int8x16,pZ0), M&15);
        AE_SA128POS_FP(alZ0, pZ0);
        return;
    }

    WUR_AE_CBEGIN0((uintptr_t)x);
    WUR_AE_CEND0((uintptr_t)(x + M*N));
    pZ0 = z;
    for (m=0; m < M; m+=4)
    {
        pX0 = x + m*N;
        pX1 = pX0;  AE_ADDCIRC_XC(castxcc(ae_int64,pX1), N);
        pX2 = pX1;  AE_ADDCIRC_XC(castxcc(ae_int64,pX2), N);
        pX3 = pX2;  AE_ADDCIRC_XC(castxcc(ae_int64,pX3), N);
        pY0 = y;
        alX2 = AE_LA128_PP(pX2);
        alX3 = AE_LA128_PP(pX3);

        Z00 = Z01 = AE_ZERO32();

        for (n=0; n<((N-1)>>4); n++)
        {
            /* load x matrix, 4x16 values, 8-bit */
            alX0 = AE_LA128_PP(pX0);  AE_LA8X8X2_IP(x00, x01, alX0, castxcc(ae_int8x16,pX0));
            alX1 = AE_LA128_PP(pX1);  AE_LA8X8X2_IP(x10, x11, alX1, castxcc(ae_int8x16,pX1));
            AE_LA8X8X2_IP(x20, x21, alX2, castxcc(ae_int8x16,pX2));
            AE_LA8X8X2_IP(x30, x31, alX3, castxcc(ae_int8x16,pX3));

            /* load y vector, 16 values, 8-bit */
            alY0 = AE_LA128_PP(pY0);  AE_LA8X8X2_IP(y00, y01, alY0, castxcc(ae_int8x16,pY0));

            /* make multiply, using 4-way 32-bit output octal mac, 8X8 bit */
            AE_MULA8Q8X8(Z00, Z01, x00, x10, x20, x30, y00);
            AE_MULA8Q8X8(Z00, Z01, x01, x11, x21, x31, y01);
        }
        /* Process last 1...16 multiplications of each output value */
        {
            /* load x matrix, 4x16 values, 8-bit */
            alX0 = AE_LA128_PP(pX0);  AE_LAV8X8X2_XP(x00, x01, alX0, castxcc(ae_int8x16,pX0), ((N-1)&15)+1);
            alX1 = AE_LA128_PP(pX1);  AE_LAV8X8X2_XP(x10, x11, alX1, castxcc(ae_int8x16,pX1), ((N-1)&15)+1);
            AE_LAV8X8X2_XP(x20, x21, alX2, castxcc(ae_int8x16,pX2), ((N-1)&15)+1);
            AE_LAV8X8X2_XP(x30, x31, alX3, castxcc(ae_int8x16,pX3), ((N-1)&15)+1);

            /* load y vector, 16 values, 8-bit */
            alY0 = AE_LA128_PP(pY0);  AE_LAV8X8X2_XP(y00, y01, alY0, castxcc(ae_int8x16,pY0), ((N-1)&15)+1);

            /* make multiply, using 4-way 32-bit output octal mac, 8X8 bit */
            AE_MULA8Q8X8(Z00, Z01, x00, x10, x20, x30, y00);
            AE_MULA8Q8X8(Z00, Z01, x01, x11, x21, x31, y01);
        }
        /* Q15 + lsh <- Q14 + 1 + lsh w/ truncation and saturation */
        t0 = AE_TRUNCA16X4F32S(Z00, Z01, 16+1+lsh);
        /* Q7 + lsh <- Q15 + lsh -8 w/ rounding and saturation */
        z0 = AE_ROUND8X8F16SASYM(t0, t0);
        alZ0 = AE_ZALIGN128();
        AE_SAV8X8X2_XP(z0, z0, alZ0, castxcc(ae_int8x16,pZ0), XT_MIN(4, M-m));
        AE_SA128POS_FP(alZ0, pZ0);
    }
}
#endif
