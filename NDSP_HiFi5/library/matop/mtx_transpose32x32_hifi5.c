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
void mtx_transpose32x32(int32_t*    y, const int32_t*     x, int M, int N)
{
    const ae_int32x4 * restrict pX0;
    const ae_int32x4 * restrict pX1;
    const ae_int32x4 * restrict pX2;
    const ae_int32x4 * restrict pX3;
          ae_int32x4 * restrict pY0;
          ae_int32x4 * restrict pY1;
          ae_int32x4 * restrict pY2;
          ae_int32x4 * restrict pY3;
    int m, n;
    static const int16_t ALIGN(32) dsel_ind[4] = { 1797, 1540, 769, 512 };
    ae_int16x4 ind;
    ind = AE_L16X4_I((const ae_int16x4*)(dsel_ind), 0);
    NASSERT(x);
    NASSERT(y);
    NASSERT(x != y);
    if (M <= 0 || N <= 0) return;
    if (N&~3)
    {
        for ( m=0; m<(M>>2); m++ )  
        {
            ae_valignx2 ax0,ax1,ax2,ax3;
            pX0=(const ae_int32x4*)(x+4*m*N);
            pX1=(      ae_int32x4*)XT_ADDX4(N,(uintptr_t)pX0);
            pX2=(      ae_int32x4*)XT_ADDX8(N,(uintptr_t)pX0);
            pX3=(      ae_int32x4*)XT_ADDX8(N,(uintptr_t)pX1);
            pY0=(      ae_int32x4*)(y);
            pY1=(      ae_int32x4*)XT_ADDX4(M,(uintptr_t)y);
            pY2=(      ae_int32x4*)XT_ADDX8(M,(uintptr_t)y);
            pY3=(      ae_int32x4*)XT_ADDX8(M,(uintptr_t)pY1);

            ax2=AE_LA128_PP(pX2);
            ax3=AE_LA128_PP(pX3);
            __Pragma("loop_count min=1")
            for ( n=0; n<(N>>2); n++ )    
            {   
                ae_int32x2 x0,x1,x2,x3,x4,x5,x6,x7;
                ae_int32x2 y0,y1,y2,y3,y4,y5,y6,y7;
                ae_int16x4 C0,C1,C2,C3,C4,C5,C6,C7;
                ae_valignx2 ay;

                ax0=AE_LA128_PP(pX0); AE_LA32X2X2_IP(x0,x1,ax0,pX0);
                ax1=AE_LA128_PP(pX1); AE_LA32X2X2_IP(x2,x3,ax1,pX1);
                                      AE_LA32X2X2_IP(x4,x5,ax2,pX2);
                                      AE_LA32X2X2_IP(x6,x7,ax3,pX3);

                AE_DSEL16X4(C0, C1, AE_MOVINT16X4_FROMINT32X2(x0), AE_MOVINT16X4_FROMINT32X2(x2), ind);
                AE_DSEL16X4(C2, C3, AE_MOVINT16X4_FROMINT32X2(x1), AE_MOVINT16X4_FROMINT32X2(x3), ind);
                AE_DSEL16X4(C4, C5, AE_MOVINT16X4_FROMINT32X2(x4), AE_MOVINT16X4_FROMINT32X2(x6), ind);
                AE_DSEL16X4(C6, C7, AE_MOVINT16X4_FROMINT32X2(x5), AE_MOVINT16X4_FROMINT32X2(x7), ind);
                y0 = AE_MOVINT32X2_FROMINT16X4(C0);
                y1 = AE_MOVINT32X2_FROMINT16X4(C1);
                y2 = AE_MOVINT32X2_FROMINT16X4(C2);
                y3 = AE_MOVINT32X2_FROMINT16X4(C3);
                y4 = AE_MOVINT32X2_FROMINT16X4(C4);
                y5 = AE_MOVINT32X2_FROMINT16X4(C5);
                y6 = AE_MOVINT32X2_FROMINT16X4(C6);
                y7 = AE_MOVINT32X2_FROMINT16X4(C7);

                ay=AE_ZALIGN128(); 
                AE_SA32X2X2_IP(y0, y4,ay,pY0); AE_SA128POS_FP(ay,pY0); pY0 = (ae_int32x4*)( (uintptr_t)pY0+M*4*sizeof(int32_t)-sizeof(ae_int32x4));
                AE_SA32X2X2_IP(y1, y5,ay,pY1); AE_SA128POS_FP(ay,pY1); pY1 = (ae_int32x4*)( (uintptr_t)pY1+M*4*sizeof(int32_t)-sizeof(ae_int32x4));
                AE_SA32X2X2_IP(y2, y6,ay,pY2); AE_SA128POS_FP(ay,pY2); pY2 = (ae_int32x4*)( (uintptr_t)pY2+M*4*sizeof(int32_t)-sizeof(ae_int32x4));
                AE_SA32X2X2_IP(y3, y7,ay,pY3); AE_SA128POS_FP(ay,pY3); pY3 = (ae_int32x4*)( (uintptr_t)pY3+M*4*sizeof(int32_t)-sizeof(ae_int32x4));
            }  
            y+=4;
        }
        y-=(M&~3);
    }

    if (M>3 && (N&3))
    {
        pX0=(const ae_int32x4*)(x+0*N+(N&~3));
        pX1=(      ae_int32x4*)XT_ADDX4(N,(uintptr_t)pX0);
        pX2=(      ae_int32x4*)XT_ADDX8(N,(uintptr_t)pX0);
        pX3=(      ae_int32x4*)XT_ADDX8(N,(uintptr_t)pX1);
        for ( n=(N&~3); n<N; n++ )    
        {    
            ae_valignx2 ay;
            ae_int32x2 x0,x1,x2,x3;
            pY0=(ae_int32x4*)(y+n*M);
            ay=AE_ZALIGN128();
            for ( m=0; m<(M>>2); m++)  
            {      
                AE_L32_XP(x0,castxcc(ae_int32,pX0),N*4*sizeof(int32_t));
                AE_L32_XP(x1,castxcc(ae_int32,pX1),N*4*sizeof(int32_t));
                AE_L32_XP(x2,castxcc(ae_int32,pX2),N*4*sizeof(int32_t));
                AE_L32_XP(x3,castxcc(ae_int32,pX3),N*4*sizeof(int32_t));
                AE_SA32X2X2_IP(AE_SEL32_HH(x0,x1),AE_SEL32_HH(x2,x3),ay,pY0);
            }  
            AE_SA128POS_FP(ay,pY0);
            pX0=(const ae_int32x4*)(((uintptr_t)pX0)+sizeof(int32_t)-((M&~3)*N)*sizeof(int32_t));
            pX1=(const ae_int32x4*)(((uintptr_t)pX1)+sizeof(int32_t)-((M&~3)*N)*sizeof(int32_t));
            pX2=(const ae_int32x4*)(((uintptr_t)pX2)+sizeof(int32_t)-((M&~3)*N)*sizeof(int32_t));
            pX3=(const ae_int32x4*)(((uintptr_t)pX3)+sizeof(int32_t)-((M&~3)*N)*sizeof(int32_t));
        }
    }
    y+=(M&~3);
    pX0=(const ae_int32x4*)(x+(M&~3)*N);
    for (m=(M&~3); m<M; m++ )  
    {    
        ae_valignx2 ax;
        pY0=(ae_int32x4*)(y);
        pY1=(ae_int32x4*)XT_ADDX4(M,(uintptr_t)pY0);
        pY2=(ae_int32x4*)XT_ADDX8(M,(uintptr_t)pY0);
        pY3=(ae_int32x4*)XT_ADDX8(M,(uintptr_t)pY1);
        ax=AE_LA128_PP(pX0);
        for ( n=0; n<(N>>2); n++)    
        {      
            ae_int32x2 x0,x1;
            AE_LA32X2X2_IP(x0,x1,ax,pX0);
            AE_S32_H_XP(x0,castxcc(ae_int32,pY0),4*M*sizeof(int32_t));
            AE_S32_L_XP(x0,castxcc(ae_int32,pY1),4*M*sizeof(int32_t));
            AE_S32_H_XP(x1,castxcc(ae_int32,pY2),4*M*sizeof(int32_t));
            AE_S32_L_XP(x1,castxcc(ae_int32,pY3),4*M*sizeof(int32_t));
        }  
        __Pragma("no_unroll");
        __Pragma("loop_count max=3");
        for (n=0; n<(N&3); n++ )    
        {      
            ae_int32x2 x0;
            AE_L32_IP  (x0,castxcc(ae_int32,pX0),sizeof(int32_t));
            AE_S32_H_XP(x0,castxcc(ae_int32,pY0),M*sizeof(int32_t));
        }  
        y++;
    }
}
