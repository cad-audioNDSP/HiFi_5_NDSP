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
void mtx_transpose16x16_fast(int16_t  *  y, const int16_t*     x, int M, int N)
{
    const ae_int16x8* restrict pX0;
    const ae_int16x8* restrict pX1;
    const ae_int16x8* restrict pX3;
    ae_int16x8* restrict pY ;
    ae_int16x8* restrict pY1;
    ae_int16x8* restrict pY3;
    ae_int16x8* restrict pY5;
    ae_int16x8* restrict pY7;
    int m, n;
    static const int16_t ALIGN(16) dsel_interleave1_tbl[]={6|(7<<8), 4|(5<<8), 2|(3<<8), 0|(1<<8)};
    ae_int16x4 dsel_interleave1=AE_L16X4_I((const ae_int16x4*)dsel_interleave1_tbl,0);
    NASSERT(x);
    NASSERT(y);
    NASSERT(x != y);
    NASSERT(N % 4 == 0 && M % 4 == 0);
    if (M <= 0 || N <= 0) return;
    NASSERT_ALIGN(x, HIFI_SIMD_WIDTH);
    NASSERT_ALIGN(y, HIFI_SIMD_WIDTH);
    if (M>7)
    {
        ae_valignx2 aY1,aY3,aY5,aY7;
        aY1=aY3=aY5=aY7=AE_ZALIGN128();
        for ( n=0; n<(N>>3); n++ )    
        {    
            pY =(ae_int16x8*)(y);
            pY1=(ae_int16x8*)XT_ADDX2(M,(uintptr_t)pY );
            pY3=(ae_int16x8*)XT_ADDX4(M,(uintptr_t)pY1);
            pY5=(ae_int16x8*)XT_ADDX8(M,(uintptr_t)pY1);
            pY7=(ae_int16x8*)XT_ADDX8(M,(uintptr_t)pY3);

            pX0=(const ae_int16x8*)(x);
            pX1=(ae_int16x8*)XT_ADDX2(N,(uintptr_t)pX0);

            __Pragma("loop_count min=1")
            for ( m=0; m<(M>>3); m++)  
            {   
                ae_int16x4 x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,xa,xb,xc,xd,xe,xf;
                AE_L16X4X2_XP(x0,x1,pX0,2*N*sizeof(int16_t));
                AE_L16X4X2_XP(x4,x5,pX0,2*N*sizeof(int16_t));
                AE_L16X4X2_XP(x8,x9,pX0,2*N*sizeof(int16_t));
                AE_L16X4X2_XP(xc,xd,pX0,2*N*sizeof(int16_t));
                x3=AE_L16X4_I((const ae_int16x4*)pX1,sizeof(ae_int16x4));  AE_L16X4_XP(x2,castxcc(ae_int16x4,pX1),2*N*sizeof(int16_t));
                x7=AE_L16X4_I((const ae_int16x4*)pX1,sizeof(ae_int16x4));  AE_L16X4_XP(x6,castxcc(ae_int16x4,pX1),2*N*sizeof(int16_t));
                xb=AE_L16X4_I((const ae_int16x4*)pX1,sizeof(ae_int16x4));  AE_L16X4_XP(xa,castxcc(ae_int16x4,pX1),2*N*sizeof(int16_t));
                xf=AE_L16X4_I((const ae_int16x4*)pX1,sizeof(ae_int16x4));  AE_L16X4_XP(xe,castxcc(ae_int16x4,pX1),2*N*sizeof(int16_t));
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

                AE_S16X4X2_XP (x0,x8,pY,2*M*sizeof(int16_t));
                AE_SA16X4X2_IP(x2,xa,aY1,pY1);
                AE_S16X4X2_XP (x4,xc,pY,2*M*sizeof(int16_t));
                AE_SA16X4X2_IP(x6,xe,aY3,pY3);
                AE_S16X4X2_XP (x1,x9,pY,2*M*sizeof(int16_t));
                AE_SA16X4X2_IP(x3,xb,aY5,pY5);
                AE_S16X4X2_XP (x5,xd,pY,sizeof(ae_int16x8)-6*M*sizeof(int16_t));
                AE_SA16X4X2_IP(x7,xf,aY7,pY7);
            }  
            AE_SA128POS_FP(aY1,pY1);
            AE_SA128POS_FP(aY3,pY3);
            AE_SA128POS_FP(aY5,pY5);
            AE_SA128POS_FP(aY7,pY7);
            x+=8;
            y+=8*M;
        }
        if (N&4)
        {   // last 4 columns if any    
            ae_valignx2 aY1,aY3;
            pY =(ae_int16x8*)(y);
            pY1=(ae_int16x8*)XT_ADDX2(M,(uintptr_t)pY );
            pY3=(ae_int16x8*)XT_ADDX4(M,(uintptr_t)pY1);
            pX0=(const ae_int16x8*)(x);
            pX1=(ae_int16x8*)XT_ADDX2(N,(uintptr_t)pX0);
            aY1=aY3=AE_ZALIGN128();
            __Pragma("loop_count min=1")
            for ( m=0; m<(M>>3); m++)  
            {   
                ae_int16x4 x0,x1,x2,x3,x4,x5,x6,x7;
                AE_L16X4_XP(x0,castxcc(ae_int16x4,pX0),2*N*sizeof(int16_t));
                AE_L16X4_XP(x1,castxcc(ae_int16x4,pX1),2*N*sizeof(int16_t));
                AE_L16X4_XP(x2,castxcc(ae_int16x4,pX0),2*N*sizeof(int16_t));
                AE_L16X4_XP(x3,castxcc(ae_int16x4,pX1),2*N*sizeof(int16_t));
                AE_L16X4_XP(x4,castxcc(ae_int16x4,pX0),2*N*sizeof(int16_t));
                AE_L16X4_XP(x5,castxcc(ae_int16x4,pX1),2*N*sizeof(int16_t));
                AE_L16X4_XP(x6,castxcc(ae_int16x4,pX0),2*N*sizeof(int16_t));
                AE_L16X4_XP(x7,castxcc(ae_int16x4,pX1),2*N*sizeof(int16_t));
                AE_DSEL16X4(x0,x1,x0,x1,dsel_interleave1);  
                AE_DSEL16X4(x2,x3,x2,x3,dsel_interleave1);  
                AE_DSEL16X4(x0,x2,x0,x2,dsel_interleave1);  
                AE_DSEL16X4(x1,x3,x1,x3,dsel_interleave1);  
                AE_DSEL16X4(x4,x5,x4,x5,dsel_interleave1);  
                AE_DSEL16X4(x6,x7,x6,x7,dsel_interleave1);  
                AE_DSEL16X4(x4,x6,x4,x6,dsel_interleave1);  
                AE_DSEL16X4(x5,x7,x5,x7,dsel_interleave1);  
                AE_S16X4X2_XP (x0,x4,pY,2*M*sizeof(int16_t));
                AE_SA16X4X2_IP(x1,x5,aY1,pY1);
                AE_S16X4X2_XP (x2,x6,pY,sizeof(ae_int16x8)-2*M*sizeof(int16_t));
                AE_SA16X4X2_IP(x3,x7,aY3,pY3);
            }  
           AE_SA128POS_FP(aY1,pY1);
           AE_SA128POS_FP(aY3,pY3);
        }
        x-=(N&~7);
        y-=(N&~7)*M;
    }
    if (M&4)
    {    // last 4 rows if any
        ae_valignx2 aX1,aX3;
        pX0=(const ae_int16x8*)(x+(M&~7)*N);
        pX1=(const ae_int16x8*)XT_ADDX2(N,(uintptr_t)pX0);
        pX3=(const ae_int16x8*)XT_ADDX4(N,(uintptr_t)pX1);
        pY = (ae_int16x8*)(y+(M&~7));
        pY1=(ae_int16x8*)XT_ADDX2(M,(uintptr_t)pY);
        aX1=AE_LA128_PP(pX1);
        aX3=AE_LA128_PP(pX3);
        for ( n=0; n<(N>>3); n++ )    
        {   
            ae_int16x4 x0,x1,x2,x3,x4,x5,x6,x7;
            AE_L16X4X2_X  (x4,x5,pX0,2*N*sizeof(int16_t));
            AE_L16X4X2_IP (x0,x1,pX0,sizeof(ae_int16x8));
            AE_LA16X4X2_IP(x2,x3,aX1,pX1);
            AE_LA16X4X2_IP(x6,x7,aX3,pX3);
            AE_DSEL16X4(x0,x2,x0,x2,dsel_interleave1);  
            AE_DSEL16X4(x4,x6,x4,x6,dsel_interleave1);  
            AE_DSEL16X4(x0,x4,x0,x4,dsel_interleave1);  
            AE_DSEL16X4(x2,x6,x2,x6,dsel_interleave1);  
            AE_DSEL16X4(x1,x3,x1,x3,dsel_interleave1);  
            AE_DSEL16X4(x5,x7,x5,x7,dsel_interleave1);  
            AE_DSEL16X4(x1,x5,x1,x5,dsel_interleave1);  
            AE_DSEL16X4(x3,x7,x3,x7,dsel_interleave1);  
            AE_S16X4_XP(x0,castxcc(ae_int16x4,pY ),2*M*sizeof(int16_t));
            AE_S16X4_XP(x2,castxcc(ae_int16x4,pY1),2*M*sizeof(int16_t));
            AE_S16X4_XP(x4,castxcc(ae_int16x4,pY ),2*M*sizeof(int16_t));
            AE_S16X4_XP(x6,castxcc(ae_int16x4,pY1),2*M*sizeof(int16_t));
            AE_S16X4_XP(x1,castxcc(ae_int16x4,pY ),2*M*sizeof(int16_t));
            AE_S16X4_XP(x3,castxcc(ae_int16x4,pY1),2*M*sizeof(int16_t));
            AE_S16X4_XP(x5,castxcc(ae_int16x4,pY ),2*M*sizeof(int16_t));
            AE_S16X4_XP(x7,castxcc(ae_int16x4,pY1),2*M*sizeof(int16_t));
        }  
        if (N&4)
        {      
            ae_int16x4 x0,x2,x4,x6;
            x4=AE_L16X4_X  ((const ae_int16x4*)pX0 ,2*N*sizeof(int16_t));
            AE_L16X4_IP (x0,castxcc(ae_int16x4,pX0),sizeof(ae_int16x4));
            AE_L16X4_IP (x2,castxcc(ae_int16x4,pX1),sizeof(ae_int16x4));
            AE_L16X4_IP (x6,castxcc(ae_int16x4,pX3),sizeof(ae_int16x4));
            AE_DSEL16X4(x0,x2,x0,x2,dsel_interleave1);  
            AE_DSEL16X4(x4,x6,x4,x6,dsel_interleave1);  
            AE_DSEL16X4(x0,x4,x0,x4,dsel_interleave1);  
            AE_DSEL16X4(x2,x6,x2,x6,dsel_interleave1);  
            AE_S16X4_XP(x0,castxcc(ae_int16x4,pY ),2*M*sizeof(int16_t));
            AE_S16X4_XP(x2,castxcc(ae_int16x4,pY1),2*M*sizeof(int16_t));
            AE_S16X4_XP(x4,castxcc(ae_int16x4,pY ),2*M*sizeof(int16_t));
            AE_S16X4_XP(x6,castxcc(ae_int16x4,pY1),2*M*sizeof(int16_t));
        }  
    }
}
