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
   Code optimized for HiFi5 core
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
void mtx_transpose16x16(int16_t* restrict   y, const int16_t* restrict  x, int M, int N)
{
#if !defined(AE_LAV16X4X2_XP) || !defined(AE_SAV16X4X2_XP)
    ae_int16* restrict pY0;
    ae_int16* restrict pY1;
    ae_int16* restrict pY2;
    ae_int16* restrict pY3;
    ae_int16* restrict pY4;
    ae_int16* restrict pY5;
    ae_int16* restrict pY6;
    ae_int16* restrict pY7;
    const ae_int16* restrict pX0;
    const ae_int16* restrict pX1;
    const ae_int16* restrict pX2;
    const ae_int16* restrict pX3;

    int m, n;
    static const int16_t ALIGN(16) dsel_interleave1_tbl[]={6|(7<<8), 4|(5<<8), 2|(3<<8), 0|(1<<8)};
    ae_int16x4 dsel_interleave1=AE_L16X4_I((const ae_int16x4*)dsel_interleave1_tbl,0);

    NASSERT( x );
    NASSERT( y );
    NASSERT(x != y);
    if (M <= 0 || N <= 0) return;

    if (M>7)
    {
        ae_valignx2 aX,aY,aY0;
        aY0=aY=AE_ZALIGN128();
        // by 8x8 blocks
        for ( n=0; n<(N&~7); n+=8 )    
        {    
            pX0=(const ae_int16*)(x);
            pX1=(const ae_int16*)XT_ADDX2(N,(uintptr_t)pX0);
            pY0=(ae_int16*)(y);
            pY1=(ae_int16*)XT_ADDX2(M,(uintptr_t)pY0);
            pY2=(ae_int16*)XT_ADDX4(M,(uintptr_t)pY0);
            pY4=(ae_int16*)XT_ADDX8(M,(uintptr_t)pY0);
            pY3=(ae_int16*)XT_ADDX4(M,(uintptr_t)pY1);
            pY5=(ae_int16*)XT_ADDX8(M,(uintptr_t)pY1);
            pY6=(ae_int16*)XT_ADDX8(M,(uintptr_t)pY2);
            pY7=(ae_int16*)XT_ADDX8(M,(uintptr_t)pY3);
            __Pragma("loop_count min=1")
            for ( m=0; m<(M>>3); m++)  
            {      
                ae_int16x4 x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,xa,xb,xc,xd,xe,xf;
                aX=AE_LA128_PP(pX0);AE_LA16X4X2_IP(x0,x1,aX,castxcc(ae_int16x8,pX0)); pX0+=2*N-8;
                aX=AE_LA128_PP(pX1);AE_LA16X4X2_IP(x2,x3,aX,castxcc(ae_int16x8,pX1)); pX1+=2*N-8;
                aX=AE_LA128_PP(pX0);AE_LA16X4X2_IP(x4,x5,aX,castxcc(ae_int16x8,pX0)); pX0+=2*N-8;
                aX=AE_LA128_PP(pX1);AE_LA16X4X2_IP(x6,x7,aX,castxcc(ae_int16x8,pX1)); pX1+=2*N-8;
                aX=AE_LA128_PP(pX0);AE_LA16X4X2_IP(x8,x9,aX,castxcc(ae_int16x8,pX0)); pX0+=2*N-8;
                aX=AE_LA128_PP(pX1);AE_LA16X4X2_IP(xa,xb,aX,castxcc(ae_int16x8,pX1)); pX1+=2*N-8;
                aX=AE_LA128_PP(pX0);AE_LA16X4X2_IP(xc,xd,aX,castxcc(ae_int16x8,pX0)); pX0+=2*N-8;
                aX=AE_LA128_PP(pX1);AE_LA16X4X2_IP(xe,xf,aX,castxcc(ae_int16x8,pX1)); pX1+=2*N-8;

                // transpose 4 4x4 blocks
                AE_DSEL16X4(x0,x2,x0,x2,dsel_interleave1);  
                AE_DSEL16X4(x4,x6,x4,x6,dsel_interleave1);  
                AE_DSEL16X4(x0,x4,x0,x4,dsel_interleave1);  
                AE_DSEL16X4(x2,x6,x2,x6,dsel_interleave1);  

                AE_DSEL16X4(x1,x3,x1,x3,dsel_interleave1);  
                AE_DSEL16X4(x5,x7,x5,x7,dsel_interleave1);  
                AE_DSEL16X4(x1,x5,x1,x5,dsel_interleave1);  
                AE_DSEL16X4(x3,x7,x3,x7,dsel_interleave1);  

                AE_DSEL16X4(x8,xa,x8,xa,dsel_interleave1);  
                AE_DSEL16X4(xc,xe,xc,xe,dsel_interleave1);  
                AE_DSEL16X4(x8,xc,x8,xc,dsel_interleave1);  
                AE_DSEL16X4(xa,xe,xa,xe,dsel_interleave1);  

                AE_DSEL16X4(x9,xb,x9,xb,dsel_interleave1);  
                AE_DSEL16X4(xd,xf,xd,xf,dsel_interleave1);  
                AE_DSEL16X4(x9,xd,x9,xd,dsel_interleave1);  
                AE_DSEL16X4(xb,xf,xb,xf,dsel_interleave1);  

                AE_SA16X4X2_IP(x0,x8,aY0,castxcc(ae_int16x8,pY0)); //AE_SA128POS_FP(aY0,pY0);
                AE_SA16X4X2_IP(x2,xa,aY ,castxcc(ae_int16x8,pY1)); AE_SA128POS_FP(aY,pY1);
                AE_SA16X4X2_IP(x4,xc,aY ,castxcc(ae_int16x8,pY2)); AE_SA128POS_FP(aY,pY2);
                AE_SA16X4X2_IP(x6,xe,aY ,castxcc(ae_int16x8,pY3)); AE_SA128POS_FP(aY,pY3);
                AE_SA16X4X2_IP(x1,x9,aY ,castxcc(ae_int16x8,pY4)); AE_SA128POS_FP(aY,pY4);
                AE_SA16X4X2_IP(x3,xb,aY ,castxcc(ae_int16x8,pY5)); AE_SA128POS_FP(aY,pY5);
                AE_SA16X4X2_IP(x5,xd,aY ,castxcc(ae_int16x8,pY6)); AE_SA128POS_FP(aY,pY6);
                AE_SA16X4X2_IP(x7,xf,aY ,castxcc(ae_int16x8,pY7)); AE_SA128POS_FP(aY,pY7);
            }  
            AE_SA128POS_FP(aY0,pY0);
            x+=8;
            y+=8*M;
        }
        if (N&4)
        {    // by 4 columns
            ae_valign aX;
            pX0=(ae_int16*)(x);
            pX1=(const ae_int16*)XT_ADDX2(N,(uintptr_t)pX0);
            pY0=(ae_int16*)(y);
            pY1=(ae_int16*)XT_ADDX2(M,(uintptr_t)pY0);
            pY2=(ae_int16*)XT_ADDX4(M,(uintptr_t)pY0);
            pY3=(ae_int16*)XT_ADDX4(M,(uintptr_t)pY1);
            __Pragma("loop_count min=1")
            for ( m=0; m<(M>>3); m++ )  
            {      
                ae_int16x4 x0,x2,x4,x6,x8,xa,xc,xe;
                aX=AE_LA64_PP(pX0);AE_LA16X4_IP(x0,aX,castxcc(ae_int16x4,pX0)); pX0+=2*N-4;
                aX=AE_LA64_PP(pX1);AE_LA16X4_IP(x2,aX,castxcc(ae_int16x4,pX1)); pX1+=2*N-4;
                aX=AE_LA64_PP(pX0);AE_LA16X4_IP(x4,aX,castxcc(ae_int16x4,pX0)); pX0+=2*N-4;
                aX=AE_LA64_PP(pX1);AE_LA16X4_IP(x6,aX,castxcc(ae_int16x4,pX1)); pX1+=2*N-4;
                aX=AE_LA64_PP(pX0);AE_LA16X4_IP(x8,aX,castxcc(ae_int16x4,pX0)); pX0+=2*N-4;
                aX=AE_LA64_PP(pX1);AE_LA16X4_IP(xa,aX,castxcc(ae_int16x4,pX1)); pX1+=2*N-4;
                aX=AE_LA64_PP(pX0);AE_LA16X4_IP(xc,aX,castxcc(ae_int16x4,pX0)); pX0+=2*N-4;
                aX=AE_LA64_PP(pX1);AE_LA16X4_IP(xe,aX,castxcc(ae_int16x4,pX1)); pX1+=2*N-4;

                // transpose 4 4x4 blocks
                AE_DSEL16X4(x0,x2,x0,x2,dsel_interleave1);  
                AE_DSEL16X4(x4,x6,x4,x6,dsel_interleave1);  
                AE_DSEL16X4(x0,x4,x0,x4,dsel_interleave1);  
                AE_DSEL16X4(x2,x6,x2,x6,dsel_interleave1);  

                AE_DSEL16X4(x8,xa,x8,xa,dsel_interleave1);  
                AE_DSEL16X4(xc,xe,xc,xe,dsel_interleave1);  
                AE_DSEL16X4(x8,xc,x8,xc,dsel_interleave1);  
                AE_DSEL16X4(xa,xe,xa,xe,dsel_interleave1);  

                AE_SA16X4X2_IP(x0,x8,aY0,castxcc(ae_int16x8,pY0)); //AE_SA128POS_FP(aY0,pY0);
                AE_SA16X4X2_IP(x2,xa,aY ,castxcc(ae_int16x8,pY1)); AE_SA128POS_FP(aY,pY1);
                AE_SA16X4X2_IP(x4,xc,aY ,castxcc(ae_int16x8,pY2)); AE_SA128POS_FP(aY,pY2);
                AE_SA16X4X2_IP(x6,xe,aY ,castxcc(ae_int16x8,pY3)); AE_SA128POS_FP(aY,pY3);
            }  
            AE_SA128POS_FP(aY0,pY0);
            n+=4;
            x+=4;
            y+=4*M;
        }
        if (N&3)
        {   // column by column
            for ( n=(N&~3); n<N; n++ )    
            {    
                pY0=(ae_int16*)(y);
                pX0=(ae_int16*)(x);
                pX1=(const ae_int16*)XT_ADDX2(N,(uintptr_t)pX0);
                __Pragma("loop_count min=1")
                for ( m=0; m<(M>>3); m++ )  
                {      
                    ae_int16x4 x0,x1,x2,x3,x4,x5,x6,x7;
                    x2=AE_L16_X (pX0,2*N*sizeof(int16_t)); AE_L16_XP(x0,pX0,4*N*sizeof(int16_t));
                    x3=AE_L16_X (pX1,2*N*sizeof(int16_t)); AE_L16_XP(x1,pX1,4*N*sizeof(int16_t));
                    x6=AE_L16_X (pX0,2*N*sizeof(int16_t)); AE_L16_XP(x4,pX0,4*N*sizeof(int16_t));
                    x7=AE_L16_X (pX1,2*N*sizeof(int16_t)); AE_L16_XP(x5,pX1,4*N*sizeof(int16_t));

                    AE_DSEL16X4(x0,x1,x0,x1,dsel_interleave1);  
                    AE_DSEL16X4(x2,x3,x2,x3,dsel_interleave1);  
                    AE_DSEL16X4(x0,x2,x0,x2,dsel_interleave1);  
                    AE_DSEL16X4(x1,x3,x1,x3,dsel_interleave1);  

                    AE_DSEL16X4(x4,x5,x4,x5,dsel_interleave1);  
                    AE_DSEL16X4(x6,x7,x6,x7,dsel_interleave1);  
                    AE_DSEL16X4(x4,x6,x4,x6,dsel_interleave1);  
                    AE_DSEL16X4(x5,x7,x5,x7,dsel_interleave1);  
                    AE_SA16X4X2_IP(x0,x4,aY0,castxcc(ae_int16x8,pY0)); //AE_SA128POS_FP(aY0,pY0);
                }  
                AE_SA128POS_FP(aY0,pY0);
                x++;
                y+=M;
            }
        }
        x-=N;
        y-=M*N;
    }

    if (N>7)
    {
        m=(M&~7);
        // last 4 rows if any
        if (M&4)
        {
            ae_valign aY,aY0;
            aY=aY0=AE_ZALIGN64();
            pX0=(const ae_int16*)(x+m*N);
            pX1=(const ae_int16*)XT_ADDX2(N,(uintptr_t)pX0);
            pX2=(const ae_int16*)XT_ADDX4(N,(uintptr_t)pX0);
            pX3=(const ae_int16*)XT_ADDX4(N,(uintptr_t)pX1);

            pY0=(ae_int16*)(y+m);
            pY1=(ae_int16*)XT_ADDX2(M,(uintptr_t)pY0);
            __Pragma("loop_count min=1")
            for ( n=0; n<(N>>3); n++ )    
            {      
                ae_valignx2 aX;
                ae_int16x4 x0,x1,x2,x3,x4,x5,x6,x7;

                aX=AE_LA128_PP(pX0);AE_LA16X4X2_IP(x0,x1,aX,castxcc(ae_int16x8,pX0)); 
                aX=AE_LA128_PP(pX1);AE_LA16X4X2_IP(x2,x3,aX,castxcc(ae_int16x8,pX1)); 
                aX=AE_LA128_PP(pX2);AE_LA16X4X2_IP(x4,x5,aX,castxcc(ae_int16x8,pX2)); 
                aX=AE_LA128_PP(pX3);AE_LA16X4X2_IP(x6,x7,aX,castxcc(ae_int16x8,pX3)); 

                // transpose 2 4x4 blocks
                AE_DSEL16X4(x0,x2,x0,x2,dsel_interleave1);  
                AE_DSEL16X4(x4,x6,x4,x6,dsel_interleave1);  
                AE_DSEL16X4(x0,x4,x0,x4,dsel_interleave1);  
                AE_DSEL16X4(x2,x6,x2,x6,dsel_interleave1);  
 
                AE_DSEL16X4(x1,x3,x1,x3,dsel_interleave1);  
                AE_DSEL16X4(x5,x7,x5,x7,dsel_interleave1);  
                AE_DSEL16X4(x1,x5,x1,x5,dsel_interleave1);  
                AE_DSEL16X4(x3,x7,x3,x7,dsel_interleave1);  

                AE_SA16X4_IP(x0,aY ,castxcc(ae_int16x4,pY0)); AE_SA64POS_FP(aY,pY0); pY0+=2*M-4;
                AE_SA16X4_IP(x2,aY ,castxcc(ae_int16x4,pY1)); AE_SA64POS_FP(aY,pY1); pY1+=2*M-4;
                AE_SA16X4_IP(x4,aY ,castxcc(ae_int16x4,pY0)); AE_SA64POS_FP(aY,pY0); pY0+=2*M-4;
                AE_SA16X4_IP(x6,aY ,castxcc(ae_int16x4,pY1)); AE_SA64POS_FP(aY,pY1); pY1+=2*M-4;
                AE_SA16X4_IP(x1,aY ,castxcc(ae_int16x4,pY0)); AE_SA64POS_FP(aY,pY0); pY0+=2*M-4;
                AE_SA16X4_IP(x3,aY ,castxcc(ae_int16x4,pY1)); AE_SA64POS_FP(aY,pY1); pY1+=2*M-4;
                AE_SA16X4_IP(x5,aY ,castxcc(ae_int16x4,pY0)); AE_SA64POS_FP(aY,pY0); pY0+=2*M-4;
                AE_SA16X4_IP(x7,aY ,castxcc(ae_int16x4,pY1)); AE_SA64POS_FP(aY,pY1); pY1+=2*M-4;
            }  
            m+=4;
        }
        // last 1..3 rows
        if (M&3)
        {
            for (; m<M; m++ )  
            {    
                ae_int16x4 x0,x1;
                ae_valignx2 aX;
                pX0=(const ae_int16*)(x+(m+0)*N);
                pY0=(ae_int16*)(y+m);
                pY1=(ae_int16*)XT_ADDX2(M,(uintptr_t)pY0);
                aX=AE_LA128_PP(pX0);
                __Pragma("loop_count min=1")
                for ( n=0; n<(N>>3); n++ )    
                {      
                    AE_LA16X4X2_IP(x0,x1,aX,castxcc(ae_int16x8,pX0)); 
                    AE_S16_0_XP(AE_SEL16_6543(x0,x0),pY0,2*M*sizeof(int16_t)); 
                    AE_S16_0_XP(AE_SEL16_5432(x0,x0),pY1,2*M*sizeof(int16_t)); 
                    AE_S16_0_XP(AE_SEL16_4321(x0,x0),pY0,2*M*sizeof(int16_t)); 
                    AE_S16_0_XP(x0,pY1,2*M*sizeof(int16_t)); 
                    AE_S16_0_XP(AE_SEL16_6543(x1,x1),pY0,2*M*sizeof(int16_t)); 
                    AE_S16_0_XP(AE_SEL16_5432(x1,x1),pY1,2*M*sizeof(int16_t)); 
                    AE_S16_0_XP(AE_SEL16_4321(x1,x1),pY0,2*M*sizeof(int16_t)); 
                    AE_S16_0_XP(x1,pY1,2*M*sizeof(int16_t)); 
                }  
            }
        }
    }

    // transpose small rectangle at the right bottom if needed
    if ((M&7) && (N&7))
    {
        for ( n=(N&~7); n<N; n++ )    
        {    
            __Pragma("no_unroll")
            for ( m=(M&~7); m<M; m++ )  
            {      
                y[n*M+m] = x[m*N+n];    
            }  
        }
    }
#else
    ae_int16* restrict pY0;
    ae_int16* restrict pY1;
    ae_int16* restrict pY2;
    ae_int16* restrict pY3;
    ae_int16* restrict pY4;
    ae_int16* restrict pY5;
    ae_int16* restrict pY6;
    ae_int16* restrict pY7;
    const ae_int16* restrict pX0;
    const ae_int16* restrict pX1;
    const ae_int16* restrict pX2;
    const ae_int16* restrict pX3;
    const ae_int16* restrict pX4;
    const ae_int16* restrict pX5;
    const ae_int16* restrict pX6;

    int m, n;
    static const int16_t ALIGN(16) dsel_interleave1_tbl[]={6|(7<<8), 4|(5<<8), 2|(3<<8), 0|(1<<8)};
    ae_int16x4 dsel_interleave1=AE_L16X4_I((const ae_int16x4*)dsel_interleave1_tbl,0);

    NASSERT( x );
    NASSERT( y );
    NASSERT(x != y);
    if (M <= 0 || N <= 0) return;

    /* copy right border */
    if (M>7 && (N&7))
    {
        ae_int16x4 x00,x10,x20,x30,x40,x50,x60,x70;
        ae_int16x4 x01,x11,x21,x31,x41,x51,x61,x71;
        ae_valignx2 aX, aY;

        n = (N&~7);
        WUR_AE_CBEGIN0((uintptr_t)(y + n*M));
        WUR_AE_CEND0((uintptr_t)(y + N*M));
        pX0 = (const ae_int16 *)(x + n);
        pY0 = (      ae_int16 *)(y + n*M);
        pY1 = pY0;  AE_ADDCIRC_XC(castxcc(ae_int64,pY1), M*sizeof(int16_t));
        pY2 = pY1;  AE_ADDCIRC_XC(castxcc(ae_int64,pY2), M*sizeof(int16_t));
        pY3 = pY2;  AE_ADDCIRC_XC(castxcc(ae_int64,pY3), M*sizeof(int16_t));
        pY4 = pY3;  AE_ADDCIRC_XC(castxcc(ae_int64,pY4), M*sizeof(int16_t));
        pY5 = pY4;  AE_ADDCIRC_XC(castxcc(ae_int64,pY5), M*sizeof(int16_t));
        pY6 = pY5;  AE_ADDCIRC_XC(castxcc(ae_int64,pY6), M*sizeof(int16_t));
        aY = AE_ZALIGN128();
        __Pragma("loop_count min=1")
        for ( m=0; m<(M>>3); m++ )
        {
            aX=AE_LA128_PP(pX0); AE_LAV16X4X2_XP(x00, x01, aX,castxcc(ae_int16x8,pX0), (N&7)*sizeof(int16_t)); pX0 += (N&~7);
            aX=AE_LA128_PP(pX0); AE_LAV16X4X2_XP(x10, x11, aX,castxcc(ae_int16x8,pX0), (N&7)*sizeof(int16_t)); pX0 += (N&~7);
            aX=AE_LA128_PP(pX0); AE_LAV16X4X2_XP(x20, x21, aX,castxcc(ae_int16x8,pX0), (N&7)*sizeof(int16_t)); pX0 += (N&~7);
            aX=AE_LA128_PP(pX0); AE_LAV16X4X2_XP(x30, x31, aX,castxcc(ae_int16x8,pX0), (N&7)*sizeof(int16_t)); pX0 += (N&~7);
            aX=AE_LA128_PP(pX0); AE_LAV16X4X2_XP(x40, x41, aX,castxcc(ae_int16x8,pX0), (N&7)*sizeof(int16_t)); pX0 += (N&~7);
            aX=AE_LA128_PP(pX0); AE_LAV16X4X2_XP(x50, x51, aX,castxcc(ae_int16x8,pX0), (N&7)*sizeof(int16_t)); pX0 += (N&~7);
            aX=AE_LA128_PP(pX0); AE_LAV16X4X2_XP(x60, x61, aX,castxcc(ae_int16x8,pX0), (N&7)*sizeof(int16_t)); pX0 += (N&~7);
            aX=AE_LA128_PP(pX0); AE_LAV16X4X2_XP(x70, x71, aX,castxcc(ae_int16x8,pX0), (N&7)*sizeof(int16_t)); pX0 += (N&~7);

            AE_DSEL16X4(x00, x10, x00, x10, dsel_interleave1);
            AE_DSEL16X4(x20, x30, x20, x30, dsel_interleave1);
            AE_DSEL16X4(x00, x20, x00, x20, dsel_interleave1);
            AE_DSEL16X4(x10, x30, x10, x30, dsel_interleave1);
            AE_DSEL16X4(x40, x50, x40, x50, dsel_interleave1);
            AE_DSEL16X4(x60, x70, x60, x70, dsel_interleave1);
            AE_DSEL16X4(x40, x60, x40, x60, dsel_interleave1);
            AE_DSEL16X4(x50, x70, x50, x70, dsel_interleave1);

            AE_DSEL16X4(x01, x11, x01, x11, dsel_interleave1);
            AE_DSEL16X4(x21, x31, x21, x31, dsel_interleave1);
            AE_DSEL16X4(x01, x21, x01, x21, dsel_interleave1);
            AE_DSEL16X4(x11, x31, x11, x31, dsel_interleave1);
            AE_DSEL16X4(x41, x51, x41, x51, dsel_interleave1);
            AE_DSEL16X4(x61, x71, x61, x71, dsel_interleave1);
            AE_DSEL16X4(x41, x61, x41, x61, dsel_interleave1);
            AE_DSEL16X4(x51, x71, x51, x71, dsel_interleave1);

            AE_SA16X4X2_IP(x21, x61, aY, castxcc(ae_int16x8,pY6)); AE_SA128POS_FP(aY, pY6);
            AE_SA16X4X2_IP(x11, x51, aY, castxcc(ae_int16x8,pY5)); AE_SA128POS_FP(aY, pY5);
            AE_SA16X4X2_IP(x01, x41, aY, castxcc(ae_int16x8,pY4)); AE_SA128POS_FP(aY, pY4);
            AE_SA16X4X2_IP(x30, x70, aY, castxcc(ae_int16x8,pY3)); AE_SA128POS_FP(aY, pY3);
            AE_SA16X4X2_IP(x20, x60, aY, castxcc(ae_int16x8,pY2)); AE_SA128POS_FP(aY, pY2);
            AE_SA16X4X2_IP(x10, x50, aY, castxcc(ae_int16x8,pY1)); AE_SA128POS_FP(aY, pY1);
            AE_SA16X4X2_IP(x00, x40, aY, castxcc(ae_int16x8,pY0)); AE_SA128POS_FP(aY, pY0);
        }
    }

    /* copy bottom border */
    if (N>7 && (M&7))
    {
        ae_int16x4 x00,x10,x20,x30,x40,x50,x60,x70;
        ae_int16x4 x01,x11,x21,x31,x41,x51,x61,x71;
        ae_valignx2 aX, aY;

        m = (M&~7);
        WUR_AE_CBEGIN0((uintptr_t)(x + m*N));
        WUR_AE_CEND0((uintptr_t)(x + M*N));
        pX0 = (const ae_int16 *)(x + m*N);
        pX1 = pX0;  AE_ADDCIRC_XC(castxcc(ae_int64,pX1), N*sizeof(int16_t));
        pX2 = pX1;  AE_ADDCIRC_XC(castxcc(ae_int64,pX2), N*sizeof(int16_t));
        pX3 = pX2;  AE_ADDCIRC_XC(castxcc(ae_int64,pX3), N*sizeof(int16_t));
        pX4 = pX3;  AE_ADDCIRC_XC(castxcc(ae_int64,pX4), N*sizeof(int16_t));
        pX5 = pX4;  AE_ADDCIRC_XC(castxcc(ae_int64,pX5), N*sizeof(int16_t));
        pX6 = pX5;  AE_ADDCIRC_XC(castxcc(ae_int64,pX6), N*sizeof(int16_t));
        pY0 = (ae_int16 *)(y + m);
        aY = AE_ZALIGN128();
        __Pragma("loop_count min=1")
        for ( n=0; n<(N>>3); n++ )
        {
            aX=AE_LA128_PP(pX0); AE_LA16X4X2_IP(x00, x01, aX,castxcc(ae_int16x8,pX0));
            aX=AE_LA128_PP(pX1); AE_LA16X4X2_IP(x10, x11, aX,castxcc(ae_int16x8,pX1));
            aX=AE_LA128_PP(pX2); AE_LA16X4X2_IP(x20, x21, aX,castxcc(ae_int16x8,pX2));
            aX=AE_LA128_PP(pX3); AE_LA16X4X2_IP(x30, x31, aX,castxcc(ae_int16x8,pX3));
            aX=AE_LA128_PP(pX4); AE_LA16X4X2_IP(x40, x41, aX,castxcc(ae_int16x8,pX4));
            aX=AE_LA128_PP(pX5); AE_LA16X4X2_IP(x50, x51, aX,castxcc(ae_int16x8,pX5));
            aX=AE_LA128_PP(pX6); AE_LA16X4X2_IP(x60, x61, aX,castxcc(ae_int16x8,pX6));
            x70 = AE_XOR16(x70, x70);
            x71 = AE_XOR16(x71, x71);

            AE_DSEL16X4(x00, x10, x00, x10, dsel_interleave1);
            AE_DSEL16X4(x20, x30, x20, x30, dsel_interleave1);
            AE_DSEL16X4(x00, x20, x00, x20, dsel_interleave1);
            AE_DSEL16X4(x10, x30, x10, x30, dsel_interleave1);
            AE_DSEL16X4(x40, x50, x40, x50, dsel_interleave1);
            AE_DSEL16X4(x60, x70, x60, x70, dsel_interleave1);
            AE_DSEL16X4(x40, x60, x40, x60, dsel_interleave1);
            AE_DSEL16X4(x50, x70, x50, x70, dsel_interleave1);

            AE_DSEL16X4(x01, x11, x01, x11, dsel_interleave1);
            AE_DSEL16X4(x21, x31, x21, x31, dsel_interleave1);
            AE_DSEL16X4(x01, x21, x01, x21, dsel_interleave1);
            AE_DSEL16X4(x11, x31, x11, x31, dsel_interleave1);
            AE_DSEL16X4(x41, x51, x41, x51, dsel_interleave1);
            AE_DSEL16X4(x61, x71, x61, x71, dsel_interleave1);
            AE_DSEL16X4(x41, x61, x41, x61, dsel_interleave1);
            AE_DSEL16X4(x51, x71, x51, x71, dsel_interleave1);

            AE_SAV16X4X2_XP(x00, x40, aY,castxcc(ae_int16x8,pY0), (M&7)*sizeof(int16_t)); AE_SA128POS_FP(aY,pY0); pY0+=(M&~7);
            AE_SAV16X4X2_XP(x10, x50, aY,castxcc(ae_int16x8,pY0), (M&7)*sizeof(int16_t)); AE_SA128POS_FP(aY,pY0); pY0+=(M&~7);
            AE_SAV16X4X2_XP(x20, x60, aY,castxcc(ae_int16x8,pY0), (M&7)*sizeof(int16_t)); AE_SA128POS_FP(aY,pY0); pY0+=(M&~7);
            AE_SAV16X4X2_XP(x30, x70, aY,castxcc(ae_int16x8,pY0), (M&7)*sizeof(int16_t)); AE_SA128POS_FP(aY,pY0); pY0+=(M&~7);
            AE_SAV16X4X2_XP(x01, x41, aY,castxcc(ae_int16x8,pY0), (M&7)*sizeof(int16_t)); AE_SA128POS_FP(aY,pY0); pY0+=(M&~7);
            AE_SAV16X4X2_XP(x11, x51, aY,castxcc(ae_int16x8,pY0), (M&7)*sizeof(int16_t)); AE_SA128POS_FP(aY,pY0); pY0+=(M&~7);
            AE_SAV16X4X2_XP(x21, x61, aY,castxcc(ae_int16x8,pY0), (M&7)*sizeof(int16_t)); AE_SA128POS_FP(aY,pY0); pY0+=(M&~7);
            AE_SAV16X4X2_XP(x31, x71, aY,castxcc(ae_int16x8,pY0), (M&7)*sizeof(int16_t)); AE_SA128POS_FP(aY,pY0); pY0+=(M&~7);
        }
    }

    /* transpose small rectangle at the right bottom if needed */
    if ((M&7) && (N&7))
    {
        ae_int16x4 x0;
        for ( n=(N&~7); n<N; n++ )
        {
            pX0 = (const ae_int16 *)(x + n + (M&~7)*N);
            pY0 = (      ae_int16 *)(y + n*M + (M&~7));
            __Pragma("loop_count min=1");
            for ( m=(M&~7); m<M; m++ )
            {
                AE_L16_XP(x0, pX0, N*sizeof(int16_t));
                AE_S16_0_IP(x0, pY0, sizeof(int16_t));
            }
        }
    }

    /* transpose by 8x8 cells */
    if (M>7)
    {
        ae_valignx2 aX,aY,aY0;
        aY0=aY=AE_ZALIGN128();

        for ( n=0; n<(N>>3); n++ )    
        {    
            pX0=(const ae_int16*)(x);
            pX1=(const ae_int16*)XT_ADDX2(N,(uintptr_t)pX0);
            pY0=(ae_int16*)(y);
            pY1=(ae_int16*)XT_ADDX2(M,(uintptr_t)pY0);
            pY2=(ae_int16*)XT_ADDX4(M,(uintptr_t)pY0);
            pY4=(ae_int16*)XT_ADDX8(M,(uintptr_t)pY0);
            pY3=(ae_int16*)XT_ADDX4(M,(uintptr_t)pY1);
            pY5=(ae_int16*)XT_ADDX8(M,(uintptr_t)pY1);
            pY6=(ae_int16*)XT_ADDX8(M,(uintptr_t)pY2);
            pY7=(ae_int16*)XT_ADDX8(M,(uintptr_t)pY3);
            __Pragma("loop_count min=1")
            for ( m=0; m<(M>>3); m++)  
            {
                ae_int16x4 x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,xa,xb,xc,xd,xe,xf;
                aX=AE_LA128_PP(pX0);AE_LA16X4X2_IP(x0,x1,aX,castxcc(ae_int16x8,pX0)); pX0+=2*N-8;
                aX=AE_LA128_PP(pX1);AE_LA16X4X2_IP(x2,x3,aX,castxcc(ae_int16x8,pX1)); pX1+=2*N-8;
                aX=AE_LA128_PP(pX0);AE_LA16X4X2_IP(x4,x5,aX,castxcc(ae_int16x8,pX0)); pX0+=2*N-8;
                aX=AE_LA128_PP(pX1);AE_LA16X4X2_IP(x6,x7,aX,castxcc(ae_int16x8,pX1)); pX1+=2*N-8;
                aX=AE_LA128_PP(pX0);AE_LA16X4X2_IP(x8,x9,aX,castxcc(ae_int16x8,pX0)); pX0+=2*N-8;
                aX=AE_LA128_PP(pX1);AE_LA16X4X2_IP(xa,xb,aX,castxcc(ae_int16x8,pX1)); pX1+=2*N-8;
                aX=AE_LA128_PP(pX0);AE_LA16X4X2_IP(xc,xd,aX,castxcc(ae_int16x8,pX0)); pX0+=2*N-8;
                aX=AE_LA128_PP(pX1);AE_LA16X4X2_IP(xe,xf,aX,castxcc(ae_int16x8,pX1)); pX1+=2*N-8;

                // transpose 4 4x4 blocks
                AE_DSEL16X4(x0,x2,x0,x2,dsel_interleave1);
                AE_DSEL16X4(x4,x6,x4,x6,dsel_interleave1);
                AE_DSEL16X4(x0,x4,x0,x4,dsel_interleave1);
                AE_DSEL16X4(x2,x6,x2,x6,dsel_interleave1);

                AE_DSEL16X4(x1,x3,x1,x3,dsel_interleave1);
                AE_DSEL16X4(x5,x7,x5,x7,dsel_interleave1);
                AE_DSEL16X4(x1,x5,x1,x5,dsel_interleave1);
                AE_DSEL16X4(x3,x7,x3,x7,dsel_interleave1);

                AE_DSEL16X4(x8,xa,x8,xa,dsel_interleave1);
                AE_DSEL16X4(xc,xe,xc,xe,dsel_interleave1);
                AE_DSEL16X4(x8,xc,x8,xc,dsel_interleave1);
                AE_DSEL16X4(xa,xe,xa,xe,dsel_interleave1);

                AE_DSEL16X4(x9,xb,x9,xb,dsel_interleave1);
                AE_DSEL16X4(xd,xf,xd,xf,dsel_interleave1);
                AE_DSEL16X4(x9,xd,x9,xd,dsel_interleave1);
                AE_DSEL16X4(xb,xf,xb,xf,dsel_interleave1);

                AE_SA16X4X2_IP(x0,x8,aY0,castxcc(ae_int16x8,pY0)); //AE_SA128POS_FP(aY0,pY0);
                AE_SA16X4X2_IP(x2,xa,aY ,castxcc(ae_int16x8,pY1)); AE_SA128POS_FP(aY,pY1);
                AE_SA16X4X2_IP(x4,xc,aY ,castxcc(ae_int16x8,pY2)); AE_SA128POS_FP(aY,pY2);
                AE_SA16X4X2_IP(x6,xe,aY ,castxcc(ae_int16x8,pY3)); AE_SA128POS_FP(aY,pY3);
                AE_SA16X4X2_IP(x1,x9,aY ,castxcc(ae_int16x8,pY4)); AE_SA128POS_FP(aY,pY4);
                AE_SA16X4X2_IP(x3,xb,aY ,castxcc(ae_int16x8,pY5)); AE_SA128POS_FP(aY,pY5);
                AE_SA16X4X2_IP(x5,xd,aY ,castxcc(ae_int16x8,pY6)); AE_SA128POS_FP(aY,pY6);
                AE_SA16X4X2_IP(x7,xf,aY ,castxcc(ae_int16x8,pY7)); AE_SA128POS_FP(aY,pY7);
            }
            AE_SA128POS_FP(aY0,pY0);
            x+=8;
            y+=8*M;
        }
    }
#endif
}
