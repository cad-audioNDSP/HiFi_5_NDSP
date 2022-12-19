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
void mtx_transpose8x8_fast   (int8_t   *  y, const int8_t *     x, int M, int N)
{
    ae_int8x8 dsel0,dsel1;
    static const uint64_t ALIGN(16) dseltbl[]=
    {
        0x40c851d962ea73fbULL,  // interleave1 {0,8,1,9,2,10,3,11,4,12,5,13,6,14,7,15}
        0x4051c8d96273eafbULL   // interleave2 {0,1,8,9,2,3,10,11,4,5,12,13,6,7,14,15}
    };
    const int8_t* restrict pX0;
    const int8_t* restrict pX1;
    const int8_t* restrict pX2;
    const int8_t* restrict pX3;
          int32_t* restrict pY;
    int m, n;
    ae_int8x8 tx0,tx1,tx2,tx3;
    ae_int8x8 ty0,ty1,ty2,ty3;
    NASSERT( x );
    NASSERT( y );
    NASSERT(x != y);
    NASSERT(N%4==0 && M%4==0);
    NASSERT_ALIGN(x,HIFI_SIMD_WIDTH);
    NASSERT_ALIGN(y,HIFI_SIMD_WIDTH);
    if (N<=0 || M<=0) return;
    dsel0=AE_L8X8_I((const ae_int8x8*)dseltbl,0*sizeof(ae_int16x4));
    dsel1=AE_L8X8_I((const ae_int8x8*)dseltbl,1*sizeof(ae_int16x4));
    pX0=x+0*N;
    pX1=x+1*N;
    pX2=x+2*N;
    pX3=x+3*N;
    if (N&4)
    {
        /* case with N a multiple of 8 plus 4 */
        for ( m=0; m<M; m+=4 )  
        {    
            ae_int32x2 t;
            ae_valign aX0,aX1,aX2,aX3;
            pY=(int32_t*)(y+m);
            aX0=AE_LA64_PP(pX0);
            aX1=AE_LA64_PP(pX1);
            aX2=AE_LA64_PP(pX2);
            aX3=AE_LA64_PP(pX3);
            for ( n=0; n<(N>>3); n++ )    
            {      
                AE_LA8X8_IP(tx0,aX0,castxcc(ae_int8x8,pX0));
                AE_LA8X8_IP(tx1,aX1,castxcc(ae_int8x8,pX1));
                AE_LA8X8_IP(tx2,aX2,castxcc(ae_int8x8,pX2));
                AE_LA8X8_IP(tx3,aX3,castxcc(ae_int8x8,pX3));
                AE_DSEL8X8(ty0,ty1,tx1,tx0,dsel0);
                AE_DSEL8X8(ty2,ty3,tx3,tx2,dsel0);
                AE_DSEL8X8(tx0,tx1,ty2,ty0,dsel1);
                AE_DSEL8X8(tx2,tx3,ty3,ty1,dsel1);
                AE_S32_H_XP( AE_MOVINT32X2_FROMINT8X8(tx0),castxcc(ae_int32,pY),1*M);
                AE_S32_L_XP( AE_MOVINT32X2_FROMINT8X8(tx0),castxcc(ae_int32,pY),1*M);
                AE_S32_H_XP( AE_MOVINT32X2_FROMINT8X8(tx1),castxcc(ae_int32,pY),1*M);
                AE_S32_L_XP( AE_MOVINT32X2_FROMINT8X8(tx1),castxcc(ae_int32,pY),1*M);
                AE_S32_H_XP( AE_MOVINT32X2_FROMINT8X8(tx2),castxcc(ae_int32,pY),1*M);
                AE_S32_L_XP( AE_MOVINT32X2_FROMINT8X8(tx2),castxcc(ae_int32,pY),1*M);
                AE_S32_H_XP( AE_MOVINT32X2_FROMINT8X8(tx3),castxcc(ae_int32,pY),1*M);
                AE_S32_L_XP( AE_MOVINT32X2_FROMINT8X8(tx3),castxcc(ae_int32,pY),1*M);
            }  
            AE_L32_XP(t,castxcc(ae_int32,pX0),1*sizeof(ae_int32)+3*N); tx0=AE_MOVINT8X8_FROMINT32X2(t);
            AE_L32_XP(t,castxcc(ae_int32,pX1),1*sizeof(ae_int32)+3*N); tx1=AE_MOVINT8X8_FROMINT32X2(t);
            AE_L32_XP(t,castxcc(ae_int32,pX2),1*sizeof(ae_int32)+3*N); tx2=AE_MOVINT8X8_FROMINT32X2(t);
            AE_L32_XP(t,castxcc(ae_int32,pX3),1*sizeof(ae_int32)+3*N); tx3=AE_MOVINT8X8_FROMINT32X2(t);
            AE_DSEL8X8(ty0,ty1,tx1,tx0,dsel0);
            AE_DSEL8X8(ty2,ty3,tx3,tx2,dsel0);
            AE_DSEL8X8(tx0,tx1,ty2,ty0,dsel1);
            AE_DSEL8X8(tx2,tx3,ty3,ty1,dsel1);
            AE_S32_L_XP( AE_MOVINT32X2_FROMINT8X8(tx1),castxcc(ae_int32,pY),1*M);
            AE_S32_H_XP( AE_MOVINT32X2_FROMINT8X8(tx1),castxcc(ae_int32,pY),1*M);
            AE_S32_L_XP( AE_MOVINT32X2_FROMINT8X8(tx0),castxcc(ae_int32,pY),1*M);
            AE_S32_H_XP( AE_MOVINT32X2_FROMINT8X8(tx0),castxcc(ae_int32,pY),1*M);
        }
    }
    else
    {
        /* case with N a multiple of 8 */
        for ( m=0; m<M; m+=4 )  
        {    
            pY=(int32_t*)(y+m);
            NASSERT_ALIGN(pX0,8);
            NASSERT_ALIGN(pX1,8);
            NASSERT_ALIGN(pX2,8);
            NASSERT_ALIGN(pX3,8);
            __Pragma("loop_count min=1")
            for ( n=0; n<(N>>3); n++ )    
            {      
                AE_L8X8_IP(tx0,castxcc(ae_int8x8,pX0),sizeof(ae_int8x8));
                AE_L8X8_IP(tx1,castxcc(ae_int8x8,pX1),sizeof(ae_int8x8));
                AE_L8X8_IP(tx2,castxcc(ae_int8x8,pX2),sizeof(ae_int8x8));
                AE_L8X8_IP(tx3,castxcc(ae_int8x8,pX3),sizeof(ae_int8x8));
                AE_DSEL8X8(ty0,ty1,tx1,tx0,dsel0);
                AE_DSEL8X8(ty2,ty3,tx3,tx2,dsel0);
                AE_DSEL8X8(tx0,tx1,ty2,ty0,dsel1);
                AE_DSEL8X8(tx2,tx3,ty3,ty1,dsel1);
                AE_S32_H_XP( AE_MOVINT32X2_FROMINT8X8(tx0),castxcc(ae_int32,pY),1*M);
                AE_S32_L_XP( AE_MOVINT32X2_FROMINT8X8(tx0),castxcc(ae_int32,pY),1*M);
                AE_S32_H_XP( AE_MOVINT32X2_FROMINT8X8(tx1),castxcc(ae_int32,pY),1*M);
                AE_S32_L_XP( AE_MOVINT32X2_FROMINT8X8(tx1),castxcc(ae_int32,pY),1*M);
                AE_S32_H_XP( AE_MOVINT32X2_FROMINT8X8(tx2),castxcc(ae_int32,pY),1*M);
                AE_S32_L_XP( AE_MOVINT32X2_FROMINT8X8(tx2),castxcc(ae_int32,pY),1*M);
                AE_S32_H_XP( AE_MOVINT32X2_FROMINT8X8(tx3),castxcc(ae_int32,pY),1*M);
                AE_S32_L_XP( AE_MOVINT32X2_FROMINT8X8(tx3),castxcc(ae_int32,pY),1*M);
            }  
            pX0+=3*N;
            pX1+=3*N;
            pX2+=3*N;
            pX3+=3*N;
        }
    }
}
