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
/*
   Matrix transpose
   Optimized for HiFi5
 */

/* Cross-platform data type definitions. */
#include "NatureDSP_types.h"
/* Library API */
#include "NatureDSP_Signal_matop.h"
/* Common helper macros. */
#include "common.h"

/*-------------------------------------------------------------------------
  Matrix Transpose
  These functions transpose matrices.

  Precision: 
  32x32 32-bit input, 32-bit output
  16x16 16-bit input, 16-bit output
  8x8   8-bit inputs, 8-bit output
  f     floating point

  Input:
  x[M][N] input matrix,Q31,Q15,Q7 or floating point
  M       number of rows in matrix x
  N       number of columns in matrix x
  Output:
  y[N][M] output vector,Q31,Q15,Q7 or floating point

  Restriction:
  For regular routines (mtx_transpose_32x32, mtx_transpose_16x16, 
  mtx_transpose_8x8, mtx_transposef):
  x,y should not overlap

  For faster routines (mtx_transpose 32x32_fast, mtx_transpose 16x16_fast, 
  mtx_transpose_8x8_fast, mtx_transposef_fast)
  x,y   should not overlap
  x,y   aligned on 16-byte boundary
  N and M are multiples of 4
-------------------------------------------------------------------------*/
void mtx_transpose8x8   (int8_t *   restrict y, const int8_t * restrict x, int M, int N)
#if !defined(AE_LAV8X8X2_XP) || !defined(AE_SAV8X8X2_XP)
{
          int8_t* restrict pY0;
          int8_t* restrict pY1;
          int8_t* restrict pY2;
          int8_t* restrict pY3;
          int8_t* restrict pY4;
          int8_t* restrict pY5;
          int8_t* restrict pY6;
          int8_t* restrict pY7;
    const int8_t* restrict pX0;
    const int8_t* restrict pX1;
    const int8_t* restrict pX2;
    const int8_t* restrict pX3;
    const int8_t* restrict pX4;
    const int8_t* restrict pX5;
    const int8_t* restrict pX6;
    const int8_t* restrict pX7;
    static const uint64_t ALIGN(16) dseltbl[]=
    {
        0x40c851d962ea73fbULL,  // interleave1 {0,8,1,9,2,10,3,11,4,12,5,13,6,14,7,15}
        0x4051c8d96273eafbULL,  // interleave2 {0,1,8,9,2,3,10,11,4,5,12,13,6,7,14,15}
        0x40516273c8d9eafbULL   // interleave4 {0,1,2,3,8,9,10,11,4,5,6,7,12,13,14,15}
    };
    ae_int8x8 dsel0,dsel1,dsel2;
    int m, n;
    NASSERT( x );
    NASSERT( y );
    NASSERT(x != y);
    if (M<=0 || N<=0) return;
    dsel0=AE_L8X8_I((const ae_int8x8*)dseltbl,0*sizeof(ae_int16x4));
    dsel1=AE_L8X8_I((const ae_int8x8*)dseltbl,1*sizeof(ae_int16x4));
    dsel2=AE_L8X8_I((const ae_int8x8*)dseltbl,2*sizeof(ae_int16x4));
    /* bottom right corner */
    if ((M&7) && (N&7))
    {
        for ( m=(M&~7); m<M; m++ )  
        for ( n=(N&~7); n<N; n++ )    
        {    
            y[n*M+m] = x[m*N+n];
        }
    }
    /* copy bottom border */
    if (N>=8 && (M&7))
    {
        if (M<8)
        {
            for ( n=0     ; n<(N&~7); n++ )    
            for ( m=(M&~7); m<M; m++ )  
            {    
                y[n*M+m] = x[m*N+n];
            }
        }
        else
        {
            /* if matrix height >=8 we may just copy last 8 rows */
            ae_valign aX,aY;
            m=M-8;
            pY0=(int8_t*)(y+m);
            pX0=x+(m+0)*N;
            pX1=pX0+N;
            pX2=(const int8_t*)XT_ADDX2(N,(uintptr_t)pX0);
            pX3=(const int8_t*)XT_ADDX2(N,(uintptr_t)pX1);
            pX4=(const int8_t*)XT_ADDX4(N,(uintptr_t)pX0);
            pX5=(const int8_t*)XT_ADDX4(N,(uintptr_t)pX1);
            pX6=(const int8_t*)XT_ADDX4(N,(uintptr_t)pX2);
            pX7=(const int8_t*)XT_ADDX4(N,(uintptr_t)pX3);
            aY=AE_ZALIGN64();
            __Pragma("loop_count min=1")
            for ( n=0; n<(N&~7); n+=8 )    
            {
                ae_int8x8 x0,x1,x2,x3,x4,x5,x6,x7;
                ae_int8x8 y0,y1,y2,y3,y4,y5,y6,y7;
                aX=AE_LA64_PP(pX0); AE_LA8X8_IP(x0,aX,castxcc(ae_int8x8,pX0));
                aX=AE_LA64_PP(pX1); AE_LA8X8_IP(x1,aX,castxcc(ae_int8x8,pX1));
                aX=AE_LA64_PP(pX2); AE_LA8X8_IP(x2,aX,castxcc(ae_int8x8,pX2));
                aX=AE_LA64_PP(pX3); AE_LA8X8_IP(x3,aX,castxcc(ae_int8x8,pX3));
                aX=AE_LA64_PP(pX4); AE_LA8X8_IP(x4,aX,castxcc(ae_int8x8,pX4));
                aX=AE_LA64_PP(pX5); AE_LA8X8_IP(x5,aX,castxcc(ae_int8x8,pX5));
                aX=AE_LA64_PP(pX6); AE_LA8X8_IP(x6,aX,castxcc(ae_int8x8,pX6));
                aX=AE_LA64_PP(pX7); AE_LA8X8_IP(x7,aX,castxcc(ae_int8x8,pX7));
                AE_DSEL8X8(y0,y1,x0,x1,dsel0);
                AE_DSEL8X8(y2,y3,x2,x3,dsel0);
                AE_DSEL8X8(y4,y5,x4,x5,dsel0);
                AE_DSEL8X8(y6,y7,x6,x7,dsel0);
                AE_DSEL8X8(x0,x1,y0,y2,dsel1);
                AE_DSEL8X8(x2,x3,y1,y3,dsel1);
                AE_DSEL8X8(x4,x5,y4,y6,dsel1);
                AE_DSEL8X8(x6,x7,y5,y7,dsel1);
                AE_DSEL8X8(y0,y1,x0,x4,dsel2);
                AE_DSEL8X8(y2,y3,x1,x5,dsel2);
                AE_DSEL8X8(y4,y5,x2,x6,dsel2);
                AE_DSEL8X8(y6,y7,x3,x7,dsel2);
                AE_SA8X8_IP(y0,aY,castxcc(ae_int8x8,pY0)); AE_SA64POS_FP(aY,pY0); pY0+=M-8;
                AE_SA8X8_IP(y1,aY,castxcc(ae_int8x8,pY0)); AE_SA64POS_FP(aY,pY0); pY0+=M-8;
                AE_SA8X8_IP(y2,aY,castxcc(ae_int8x8,pY0)); AE_SA64POS_FP(aY,pY0); pY0+=M-8;
                AE_SA8X8_IP(y3,aY,castxcc(ae_int8x8,pY0)); AE_SA64POS_FP(aY,pY0); pY0+=M-8;
                AE_SA8X8_IP(y4,aY,castxcc(ae_int8x8,pY0)); AE_SA64POS_FP(aY,pY0); pY0+=M-8;
                AE_SA8X8_IP(y5,aY,castxcc(ae_int8x8,pY0)); AE_SA64POS_FP(aY,pY0); pY0+=M-8;
                AE_SA8X8_IP(y6,aY,castxcc(ae_int8x8,pY0)); AE_SA64POS_FP(aY,pY0); pY0+=M-8;
                AE_SA8X8_IP(y7,aY,castxcc(ae_int8x8,pY0)); AE_SA64POS_FP(aY,pY0); pY0+=M-8;
            }
        }
    }
    /* copy right border */
    if (M>=8 && (N&7))
    {
        if (N<8)
        {   // no way to optimize on very small N
            for ( n=(N&~7); n<N; n++ )    
            for ( m=0; m<(M&~7); m++ )  
            {    
                y[n*M+m] = x[m*N+n];
            }
        }
        else
        {   // here we are able to transpose 8 last columns
            ae_valign aX,aY,aY0;
            aY0=aY=AE_ZALIGN64();
            n=N-8;
            pX0=(int8_t*)(x+n);
            pX1=(int8_t*)(pX0+N);
            pY0=y+n*M;
            pY1=pY0+M;
            pY2=(int8_t*)XT_ADDX2(M,(uintptr_t)pY0);
            pY3=(int8_t*)XT_ADDX2(M,(uintptr_t)pY1);
            pY4=(int8_t*)XT_ADDX4(M,(uintptr_t)pY0);
            pY5=(int8_t*)XT_ADDX4(M,(uintptr_t)pY1);
            pY6=(int8_t*)XT_ADDX4(M,(uintptr_t)pY2);
            pY7=(int8_t*)XT_ADDX4(M,(uintptr_t)pY3);
            __Pragma("loop_count min=1")
            for ( m=0; m<(M>>3); m++ )  
            {
                ae_int8x8 x0,x1,x2,x3,x4,x5,x6,x7;
                ae_int8x8 y0,y1,y2,y3,y4,y5,y6,y7;

                aX=AE_LA64_PP(pX0);AE_LA8X8_IP(x0,aX,castxcc(ae_int8x8,pX0)); pX0+=2*N-8;
                aX=AE_LA64_PP(pX1);AE_LA8X8_IP(x1,aX,castxcc(ae_int8x8,pX1)); pX1+=2*N-8;
                aX=AE_LA64_PP(pX0);AE_LA8X8_IP(x2,aX,castxcc(ae_int8x8,pX0)); pX0+=2*N-8;
                aX=AE_LA64_PP(pX1);AE_LA8X8_IP(x3,aX,castxcc(ae_int8x8,pX1)); pX1+=2*N-8;
                aX=AE_LA64_PP(pX0);AE_LA8X8_IP(x4,aX,castxcc(ae_int8x8,pX0)); pX0+=2*N-8;
                aX=AE_LA64_PP(pX1);AE_LA8X8_IP(x5,aX,castxcc(ae_int8x8,pX1)); pX1+=2*N-8;
                aX=AE_LA64_PP(pX0);AE_LA8X8_IP(x6,aX,castxcc(ae_int8x8,pX0)); pX0+=2*N-8;
                aX=AE_LA64_PP(pX1);AE_LA8X8_IP(x7,aX,castxcc(ae_int8x8,pX1)); pX1+=2*N-8;

                AE_DSEL8X8(y0,y1,x0,x1,dsel0);
                AE_DSEL8X8(y2,y3,x2,x3,dsel0);
                AE_DSEL8X8(y4,y5,x4,x5,dsel0);
                AE_DSEL8X8(y6,y7,x6,x7,dsel0);
                AE_DSEL8X8(x0,x1,y0,y2,dsel1);
                AE_DSEL8X8(x2,x3,y1,y3,dsel1);
                AE_DSEL8X8(x4,x5,y4,y6,dsel1);
                AE_DSEL8X8(x6,x7,y5,y7,dsel1);
                AE_DSEL8X8(y0,y1,x0,x4,dsel2);
                AE_DSEL8X8(y2,y3,x1,x5,dsel2);
                AE_DSEL8X8(y4,y5,x2,x6,dsel2);
                AE_DSEL8X8(y6,y7,x3,x7,dsel2);
                AE_SA8X8_IP(y0,aY0,castxcc(ae_int8x8,pY0)); //AE_SA64POS_FP(aY,pY0);
                AE_SA8X8_IP(y1,aY,castxcc(ae_int8x8,pY1)); AE_SA64POS_FP(aY,pY1);
                AE_SA8X8_IP(y2,aY,castxcc(ae_int8x8,pY2)); AE_SA64POS_FP(aY,pY2);
                AE_SA8X8_IP(y3,aY,castxcc(ae_int8x8,pY3)); AE_SA64POS_FP(aY,pY3);
                AE_SA8X8_IP(y4,aY,castxcc(ae_int8x8,pY4)); AE_SA64POS_FP(aY,pY4);
                AE_SA8X8_IP(y5,aY,castxcc(ae_int8x8,pY5)); AE_SA64POS_FP(aY,pY5);
                AE_SA8X8_IP(y6,aY,castxcc(ae_int8x8,pY6)); AE_SA64POS_FP(aY,pY6);
                AE_SA8X8_IP(y7,aY,castxcc(ae_int8x8,pY7)); AE_SA64POS_FP(aY,pY7);
            }
            AE_SA64POS_FP(aY0,pY0);
        }
    }
    /* transpose by 8x8 cells */
    if ((N&~7)>0 && (M&~7)>0)
    {
        for ( n=0; n<(N>>3); n++ )    
        {    
            ae_valign aX,aY,aY0;
            aY0=aY=AE_ZALIGN64();
            pX0=(int8_t*)(x); x+=8;
            pX1=(int8_t*)(pX0+N);
            pY0=y; y=(int8_t*)XT_ADDX8(M,(uintptr_t)y);
            pY1=pY0+M;
            pY2=(int8_t*)XT_ADDX2(M,(uintptr_t)pY0);
            pY3=(int8_t*)XT_ADDX2(M,(uintptr_t)pY1);
            pY4=(int8_t*)XT_ADDX4(M,(uintptr_t)pY0);
            pY5=(int8_t*)XT_ADDX4(M,(uintptr_t)pY1);
            pY6=(int8_t*)XT_ADDX4(M,(uintptr_t)pY2);
            pY7=(int8_t*)XT_ADDX4(M,(uintptr_t)pY3);
            __Pragma("loop_count min=1")
            for ( m=0; m<(M>>3); m++ )  
            {
                ae_int8x8 x0,x1,x2,x3,x4,x5,x6,x7;
                ae_int8x8 y0,y1,y2,y3,y4,y5,y6,y7;

                aX=AE_LA64_PP(pX0);AE_LA8X8_IP(x0,aX,castxcc(ae_int8x8,pX0)); pX0+=2*N-8;
                aX=AE_LA64_PP(pX1);AE_LA8X8_IP(x1,aX,castxcc(ae_int8x8,pX1)); pX1+=2*N-8;
                aX=AE_LA64_PP(pX0);AE_LA8X8_IP(x2,aX,castxcc(ae_int8x8,pX0)); pX0+=2*N-8;
                aX=AE_LA64_PP(pX1);AE_LA8X8_IP(x3,aX,castxcc(ae_int8x8,pX1)); pX1+=2*N-8;
                aX=AE_LA64_PP(pX0);AE_LA8X8_IP(x4,aX,castxcc(ae_int8x8,pX0)); pX0+=2*N-8;
                aX=AE_LA64_PP(pX1);AE_LA8X8_IP(x5,aX,castxcc(ae_int8x8,pX1)); pX1+=2*N-8;
                aX=AE_LA64_PP(pX0);AE_LA8X8_IP(x6,aX,castxcc(ae_int8x8,pX0)); pX0+=2*N-8;
                aX=AE_LA64_PP(pX1);AE_LA8X8_IP(x7,aX,castxcc(ae_int8x8,pX1)); pX1+=2*N-8;

                AE_DSEL8X8(y0,y1,x0,x1,dsel0);
                AE_DSEL8X8(y2,y3,x2,x3,dsel0);
                AE_DSEL8X8(y4,y5,x4,x5,dsel0);
                AE_DSEL8X8(y6,y7,x6,x7,dsel0);
                AE_DSEL8X8(x0,x1,y0,y2,dsel1);
                AE_DSEL8X8(x2,x3,y1,y3,dsel1);
                AE_DSEL8X8(x4,x5,y4,y6,dsel1);
                AE_DSEL8X8(x6,x7,y5,y7,dsel1);
                AE_DSEL8X8(y0,y1,x0,x4,dsel2);
                AE_DSEL8X8(y2,y3,x1,x5,dsel2);
                AE_DSEL8X8(y4,y5,x2,x6,dsel2);
                AE_DSEL8X8(y6,y7,x3,x7,dsel2);
                AE_SA8X8_IP(y0,aY0,castxcc(ae_int8x8,pY0)); //AE_SA64POS_FP(aY,pY0);
                AE_SA8X8_IP(y1,aY,castxcc(ae_int8x8,pY1)); AE_SA64POS_FP(aY,pY1);
                AE_SA8X8_IP(y2,aY,castxcc(ae_int8x8,pY2)); AE_SA64POS_FP(aY,pY2);
                AE_SA8X8_IP(y3,aY,castxcc(ae_int8x8,pY3)); AE_SA64POS_FP(aY,pY3);
                AE_SA8X8_IP(y4,aY,castxcc(ae_int8x8,pY4)); AE_SA64POS_FP(aY,pY4);
                AE_SA8X8_IP(y5,aY,castxcc(ae_int8x8,pY5)); AE_SA64POS_FP(aY,pY5);
                AE_SA8X8_IP(y6,aY,castxcc(ae_int8x8,pY6)); AE_SA64POS_FP(aY,pY6);
                AE_SA8X8_IP(y7,aY,castxcc(ae_int8x8,pY7)); AE_SA64POS_FP(aY,pY7);
            }
            AE_SA64POS_FP(aY0,pY0);
        }
    }
}
#else
{
          int8_t* restrict pY0;
          int8_t* restrict pY1;
          int8_t* restrict pY2;
          int8_t* restrict pY3;
          int8_t* restrict pY4;
          int8_t* restrict pY5;
          int8_t* restrict pY6;
          int8_t* restrict pY7;
    const int8_t* restrict pX0;
    const int8_t* restrict pX1;
    const int8_t* restrict pX2;
    const int8_t* restrict pX3;
    const int8_t* restrict pX4;
    const int8_t* restrict pX5;
    const int8_t* restrict pX6;
    static const uint64_t ALIGN(16) dseltbl[]=
    {
        0x40c851d962ea73fbULL,  // interleave1 {0,8,1,9,2,10,3,11,4,12,5,13,6,14,7,15}
    };
    ae_int8x8 dsel0;
    int m, n;
    NASSERT( x );
    NASSERT( y );
    NASSERT(x != y);
    if (M<=0 || N<=0) return;
    dsel0=AE_L8X8_I((const ae_int8x8*)dseltbl,0*sizeof(ae_int16x4));

    /* bottom right corner */
    if ((M&7) && (N&7))
    {
        ae_int8x8 x0;
        for ( n=(N&~7); n<N; n++ )
        {
            pX0 = x + n + (M&~7)*N;
            pY0 = y + n*M + (M&~7);
            __Pragma("loop_count min=1");
            for ( m=(M&~7); m<M; m++ )
            {
                AE_L8_XP(x0, castxcc(ae_int8, pX0), N);
                AE_S8_0_IP(x0, castxcc(ae_int8, pY0), 1);
            }
        }
    }
    /* copy bottom border */
    if (N>7 && (M&7))
    {
        ae_int8x8 x0,x1,x2,x3,x4,x5,x6,x7;
        ae_int8x8 y0,y1,y2,y3,y4,y5,y6,y7;
        ae_valign aX;
        ae_valignx2 aY;

        m = (M&~7);
        WUR_AE_CBEGIN0((uintptr_t)(x + m*N));
        WUR_AE_CEND0((uintptr_t)(x + M*N));
        pX0 = x + m*N;
        pX1 = pX0;  AE_ADDCIRC_XC(castxcc(ae_int64,pX1), N);
        pX2 = pX1;  AE_ADDCIRC_XC(castxcc(ae_int64,pX2), N);
        pX3 = pX2;  AE_ADDCIRC_XC(castxcc(ae_int64,pX3), N);
        pX4 = pX3;  AE_ADDCIRC_XC(castxcc(ae_int64,pX4), N);
        pX5 = pX4;  AE_ADDCIRC_XC(castxcc(ae_int64,pX5), N);
        pX6 = pX5;  AE_ADDCIRC_XC(castxcc(ae_int64,pX6), N);
        pY0 = y + m;
        aY = AE_ZALIGN128();
        __Pragma("loop_count min=1")
        for ( n=0; n<(N>>3); n++ )
        {
            aX=AE_LA64_PP(pX0); AE_LA8X8_IP(x0,aX,castxcc(ae_int8x8,pX0));
            aX=AE_LA64_PP(pX1); AE_LA8X8_IP(x1,aX,castxcc(ae_int8x8,pX1));
            aX=AE_LA64_PP(pX2); AE_LA8X8_IP(x2,aX,castxcc(ae_int8x8,pX2));
            aX=AE_LA64_PP(pX3); AE_LA8X8_IP(x3,aX,castxcc(ae_int8x8,pX3));
            aX=AE_LA64_PP(pX4); AE_LA8X8_IP(x4,aX,castxcc(ae_int8x8,pX4));
            aX=AE_LA64_PP(pX5); AE_LA8X8_IP(x5,aX,castxcc(ae_int8x8,pX5));
            aX=AE_LA64_PP(pX6); AE_LA8X8_IP(x6,aX,castxcc(ae_int8x8,pX6));
            x7 = AE_INT8X8_XOR_INT8X8(x7, x7);
            AE_DSEL8X8(y0,y4,x0,x4,dsel0);
            AE_DSEL8X8(y1,y5,x1,x5,dsel0);
            AE_DSEL8X8(y2,y6,x2,x6,dsel0);
            AE_DSEL8X8(y3,y7,x3,x7,dsel0);
            AE_DSEL8X8(x0,x2,y0,y2,dsel0);
            AE_DSEL8X8(x1,x3,y1,y3,dsel0);
            AE_DSEL8X8(x4,x6,y4,y6,dsel0);
            AE_DSEL8X8(x5,x7,y5,y7,dsel0);
            AE_DSEL8X8(y0,y1,x0,x1,dsel0);
            AE_DSEL8X8(y2,y3,x2,x3,dsel0);
            AE_DSEL8X8(y4,y5,x4,x5,dsel0);
            AE_DSEL8X8(y6,y7,x6,x7,dsel0);
            AE_SAV8X8X2_XP(y0,y0,aY,castxcc(ae_int8x16,pY0), M&7); AE_SA128POS_FP(aY,pY0); pY0+=(M&~7);
            AE_SAV8X8X2_XP(y1,y1,aY,castxcc(ae_int8x16,pY0), M&7); AE_SA128POS_FP(aY,pY0); pY0+=(M&~7);
            AE_SAV8X8X2_XP(y2,y2,aY,castxcc(ae_int8x16,pY0), M&7); AE_SA128POS_FP(aY,pY0); pY0+=(M&~7);
            AE_SAV8X8X2_XP(y3,y3,aY,castxcc(ae_int8x16,pY0), M&7); AE_SA128POS_FP(aY,pY0); pY0+=(M&~7);
            AE_SAV8X8X2_XP(y4,y4,aY,castxcc(ae_int8x16,pY0), M&7); AE_SA128POS_FP(aY,pY0); pY0+=(M&~7);
            AE_SAV8X8X2_XP(y5,y5,aY,castxcc(ae_int8x16,pY0), M&7); AE_SA128POS_FP(aY,pY0); pY0+=(M&~7);
            AE_SAV8X8X2_XP(y6,y6,aY,castxcc(ae_int8x16,pY0), M&7); AE_SA128POS_FP(aY,pY0); pY0+=(M&~7);
            AE_SAV8X8X2_XP(y7,y7,aY,castxcc(ae_int8x16,pY0), M&7); AE_SA128POS_FP(aY,pY0); pY0+=(M&~7);
        }
    }
    /* copy right border */
    if (M>7 && (N&7))
    {
        ae_int8x8 x0,x1,x2,x3,x4,x5,x6,x7;
        ae_int8x8 y0,y1,y2,y3,y4,y5,y6,y7;
        ae_int8x8 tmp;
        ae_valignx2 aX;
        ae_valign aY;

        n = (N&~7);
        WUR_AE_CBEGIN0((uintptr_t)(y + n*M));
        WUR_AE_CEND0((uintptr_t)(y + N*M));
        pX0 = x + n;
        pY0 = y + n*M;
        pY1 = pY0;  AE_ADDCIRC_XC(castxcc(ae_int64,pY1), M);
        pY2 = pY1;  AE_ADDCIRC_XC(castxcc(ae_int64,pY2), M);
        pY3 = pY2;  AE_ADDCIRC_XC(castxcc(ae_int64,pY3), M);
        pY4 = pY3;  AE_ADDCIRC_XC(castxcc(ae_int64,pY4), M);
        pY5 = pY4;  AE_ADDCIRC_XC(castxcc(ae_int64,pY5), M);
        pY6 = pY5;  AE_ADDCIRC_XC(castxcc(ae_int64,pY6), M);
        aY = AE_ZALIGN64();
        __Pragma("loop_count min=1")
        for ( m=0; m<(M>>3); m++ )
        {
            aX=AE_LA128_PP(pX0); AE_LAV8X8X2_XP(x0,tmp,aX,castxcc(ae_int8x16,pX0), N&7); pX0 += (N&~7);
            aX=AE_LA128_PP(pX0); AE_LAV8X8X2_XP(x1,tmp,aX,castxcc(ae_int8x16,pX0), N&7); pX0 += (N&~7);
            aX=AE_LA128_PP(pX0); AE_LAV8X8X2_XP(x2,tmp,aX,castxcc(ae_int8x16,pX0), N&7); pX0 += (N&~7);
            aX=AE_LA128_PP(pX0); AE_LAV8X8X2_XP(x3,tmp,aX,castxcc(ae_int8x16,pX0), N&7); pX0 += (N&~7);
            aX=AE_LA128_PP(pX0); AE_LAV8X8X2_XP(x4,tmp,aX,castxcc(ae_int8x16,pX0), N&7); pX0 += (N&~7);
            aX=AE_LA128_PP(pX0); AE_LAV8X8X2_XP(x5,tmp,aX,castxcc(ae_int8x16,pX0), N&7); pX0 += (N&~7);
            aX=AE_LA128_PP(pX0); AE_LAV8X8X2_XP(x6,tmp,aX,castxcc(ae_int8x16,pX0), N&7); pX0 += (N&~7);
            aX=AE_LA128_PP(pX0); AE_LAV8X8X2_XP(x7,tmp,aX,castxcc(ae_int8x16,pX0), N&7); pX0 += (N&~7);
            AE_DSEL8X8(y0,y4,x0,x4,dsel0);
            AE_DSEL8X8(y1,y5,x1,x5,dsel0);
            AE_DSEL8X8(y2,y6,x2,x6,dsel0);
            AE_DSEL8X8(y3,y7,x3,x7,dsel0);
            AE_DSEL8X8(x0,x2,y0,y2,dsel0);
            AE_DSEL8X8(x1,x3,y1,y3,dsel0);
            AE_DSEL8X8(x4,x6,y4,y6,dsel0);
            AE_DSEL8X8(x5,x7,y5,y7,dsel0);
            AE_DSEL8X8(y0,y1,x0,x1,dsel0);
            AE_DSEL8X8(y2,y3,x2,x3,dsel0);
            AE_DSEL8X8(y4,y5,x4,x5,dsel0);
            AE_DSEL8X8(y6,y7,x6,x7,dsel0);
            AE_SA8X8_IP(y6,aY,castxcc(ae_int8x8,pY6)); AE_SA64POS_FP(aY,pY6);
            AE_SA8X8_IP(y5,aY,castxcc(ae_int8x8,pY5)); AE_SA64POS_FP(aY,pY5);
            AE_SA8X8_IP(y4,aY,castxcc(ae_int8x8,pY4)); AE_SA64POS_FP(aY,pY4);
            AE_SA8X8_IP(y3,aY,castxcc(ae_int8x8,pY3)); AE_SA64POS_FP(aY,pY3);
            AE_SA8X8_IP(y2,aY,castxcc(ae_int8x8,pY2)); AE_SA64POS_FP(aY,pY2);
            AE_SA8X8_IP(y1,aY,castxcc(ae_int8x8,pY1)); AE_SA64POS_FP(aY,pY1);
            AE_SA8X8_IP(y0,aY,castxcc(ae_int8x8,pY0)); AE_SA64POS_FP(aY,pY0);
        }
    }
    /* transpose by 8x8 cells */
    if (N>7 && M>7)
    {
        for ( n=0; n<(N>>3); n++ )    
        {    
            ae_valign aX,aY,aY0;
            aY0=aY=AE_ZALIGN64();
            pX0=(int8_t*)(x); x+=8;
            pX1=(int8_t*)(pX0+N);
            pY0=y; y=(int8_t*)XT_ADDX8(M,(uintptr_t)y);
            pY1=pY0+M;
            pY2=(int8_t*)XT_ADDX2(M,(uintptr_t)pY0);
            pY3=(int8_t*)XT_ADDX2(M,(uintptr_t)pY1);
            pY4=(int8_t*)XT_ADDX4(M,(uintptr_t)pY0);
            pY5=(int8_t*)XT_ADDX4(M,(uintptr_t)pY1);
            pY6=(int8_t*)XT_ADDX4(M,(uintptr_t)pY2);
            pY7=(int8_t*)XT_ADDX4(M,(uintptr_t)pY3);
            __Pragma("loop_count min=1")
            for ( m=0; m<(M>>3); m++ )  
            {
                ae_int8x8 x0,x1,x2,x3,x4,x5,x6,x7;
                ae_int8x8 y0,y1,y2,y3,y4,y5,y6,y7;

                aX=AE_LA64_PP(pX0);AE_LA8X8_IP(x0,aX,castxcc(ae_int8x8,pX0)); pX0+=2*N-8;
                aX=AE_LA64_PP(pX1);AE_LA8X8_IP(x1,aX,castxcc(ae_int8x8,pX1)); pX1+=2*N-8;
                aX=AE_LA64_PP(pX0);AE_LA8X8_IP(x2,aX,castxcc(ae_int8x8,pX0)); pX0+=2*N-8;
                aX=AE_LA64_PP(pX1);AE_LA8X8_IP(x3,aX,castxcc(ae_int8x8,pX1)); pX1+=2*N-8;
                aX=AE_LA64_PP(pX0);AE_LA8X8_IP(x4,aX,castxcc(ae_int8x8,pX0)); pX0+=2*N-8;
                aX=AE_LA64_PP(pX1);AE_LA8X8_IP(x5,aX,castxcc(ae_int8x8,pX1)); pX1+=2*N-8;
                aX=AE_LA64_PP(pX0);AE_LA8X8_IP(x6,aX,castxcc(ae_int8x8,pX0)); pX0+=2*N-8;
                aX=AE_LA64_PP(pX1);AE_LA8X8_IP(x7,aX,castxcc(ae_int8x8,pX1)); pX1+=2*N-8;
                AE_DSEL8X8(y0,y4,x0,x4,dsel0);
                AE_DSEL8X8(y1,y5,x1,x5,dsel0);
                AE_DSEL8X8(y2,y6,x2,x6,dsel0);
                AE_DSEL8X8(y3,y7,x3,x7,dsel0);
                AE_DSEL8X8(x0,x2,y0,y2,dsel0);
                AE_DSEL8X8(x1,x3,y1,y3,dsel0);
                AE_DSEL8X8(x4,x6,y4,y6,dsel0);
                AE_DSEL8X8(x5,x7,y5,y7,dsel0);
                AE_DSEL8X8(y0,y1,x0,x1,dsel0);
                AE_DSEL8X8(y2,y3,x2,x3,dsel0);
                AE_DSEL8X8(y4,y5,x4,x5,dsel0);
                AE_DSEL8X8(y6,y7,x6,x7,dsel0);
                AE_SA8X8_IP(y0,aY0,castxcc(ae_int8x8,pY0)); //AE_SA64POS_FP(aY,pY0);
                AE_SA8X8_IP(y1,aY,castxcc(ae_int8x8,pY1)); AE_SA64POS_FP(aY,pY1);
                AE_SA8X8_IP(y2,aY,castxcc(ae_int8x8,pY2)); AE_SA64POS_FP(aY,pY2);
                AE_SA8X8_IP(y3,aY,castxcc(ae_int8x8,pY3)); AE_SA64POS_FP(aY,pY3);
                AE_SA8X8_IP(y4,aY,castxcc(ae_int8x8,pY4)); AE_SA64POS_FP(aY,pY4);
                AE_SA8X8_IP(y5,aY,castxcc(ae_int8x8,pY5)); AE_SA64POS_FP(aY,pY5);
                AE_SA8X8_IP(y6,aY,castxcc(ae_int8x8,pY6)); AE_SA64POS_FP(aY,pY6);
                AE_SA8X8_IP(y7,aY,castxcc(ae_int8x8,pY7)); AE_SA64POS_FP(aY,pY7);
            }
            AE_SA64POS_FP(aY0,pY0);
        }
    }
}
#endif
