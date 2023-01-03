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
void mtx_mpy32x32(void* pScr,
    int32_t* restrict z,
    const int32_t* restrict x,
    const int32_t* restrict y,
    int M, int N, int P, int lsh)
{
    int m, n, p;
    int K,K1;
    int32_t * restrict Y0_;
    int32_t * restrict Y1_;
    int32_t * restrict Y2_;
    int32_t * restrict Y3_;
    const int32_t *b_y;
    const ae_int32x4 * restrict pY;
    const ae_int32x2 * restrict pY_;
    ae_int32x4 * restrict pY0;
    ae_int32x4 * restrict pY1;
    ae_int32x4 * restrict pY2;
    ae_int32x4 * restrict pY3;
    const ae_int32x4 * restrict px;
    const ae_int32x4 * restrict pX;
    const ae_int32x4 * restrict pX1;
    const ae_int32x4 * restrict pX2;
    const ae_int32x4 * restrict pX3;
    ae_int32   * restrict pZ;
    ae_valignx2 aX, aX1, aX2, aX3;

    ae_int16x4 ind;
    xtbool2 bmask0, bmask1;
    static const int16_t ALIGN(16) dsel_ind[4] = { 1797, 1540, 769, 512 };
    // set up circular buffers for y to avoid reading after the end of input
    WUR_AE_CBEGIN0(((uintptr_t)y)&~15);
    WUR_AE_CEND0  ((((uintptr_t)(y+N*P))+15)&~15);
    WUR_AE_CBEGIN1(((uintptr_t)y)&~7);
    WUR_AE_CEND1  ((((uintptr_t)(y+N*P))+7)&~7);

    NASSERT(lsh >= -31 && lsh <= 31);
    /* exceptional situations */
    if ((M <= 0) || (P <= 0)) return;
    if (N <= 0)
    {
        int m;
        for (m = 0; m < M * P; m++) z[m] = 0;
        return;
    }
    ind = AE_L16X4_I((ae_int16x4*)(dsel_ind), 0);
    bmask0 = AE_MOVBA2((N & 2) + ((N & 2) >> 1) + 2 * (int)((N & 3) == 1));  // 2*(N%4) if N%4<2, else 3 
    bmask1 = AE_MOVBA2(((int)((N & 3) == 3)) << 1);  // 2 if (N%4)=3, else 0

    K = ((N + 3)&(~3));
    K1= ((N + 7)&(~7))+4;
    Y0_ = (int32_t *)pScr;
#if 0
    Y1_ = Y0_ + K;
    Y2_ = Y1_ + K;
    Y3_ = Y2_ + K;
#else
    {   // place Y0,Y2 and Y1,Y3 to the different banks
        Y1_ = Y0_ + K1;
        Y2_ = Y1_ + K1;
        Y3_ = Y2_ + K1;
    }
#endif
    /*
    main loop unrolled 4x over P
    */
    for (p = 0; p < (P >> 2); p++, z += 4, y += 4)
    {
        NASSERT_ALIGN(Y0_, HIFI_SIMD_WIDTH);
        NASSERT_ALIGN(Y1_, HIFI_SIMD_WIDTH);
        pZ = (ae_int32*)(z);
        /* copy 4 input columns to the srcatch memory with zero padding */
        pY0 = (ae_int32x4 *)Y0_;
        pY1 = (ae_int32x4 *)Y1_;
        pY2 = (ae_int32x4 *)Y2_;
        pY3 = (ae_int32x4 *)Y3_;
        b_y = y ;
        for (n = 0; n < (N >> 2); n++)
        {
            ae_int32x2 Y0, Y1, Y2, Y3, Y4, Y5, Y6, Y7;
            ae_int16x4 C0, C1, C2, C3, C4, C5, C6, C7;
            ae_valignx2 aY0,aY1;

            pY = (const ae_int32x4 *)(b_y);
            pY_= (const ae_int32x2 *)XT_ADDX2(2*P,(uintptr_t)b_y);
            aY0 = AE_LA128_PP(pY ); AE_LA32X2X2_IP(Y0, Y1, aY0, pY );
            aY1 = AE_LA128_PP(pY_); AE_LA32X2X2_IP(Y2, Y3, aY1, castxcc(ae_int32x4,pY_));
            pY = (const ae_int32x4 *)XT_ADDX4(2*P,(uintptr_t)b_y);
            pY_= (const ae_int32x2 *)(b_y + 3 * P);
            aY0 = AE_LA128_PP(pY ); AE_LA32X2X2_IP(Y4, Y5, aY0, pY );
            aY1 = AE_LA128_PP(pY_); AE_LA32X2X2_IP(Y6, Y7, aY1, castxcc(ae_int32x4,pY_));
            AE_DSEL16X4(C0, C1, AE_MOVINT16X4_FROMINT32X2(Y0), AE_MOVINT16X4_FROMINT32X2(Y2), ind);
            AE_DSEL16X4(C2, C3, AE_MOVINT16X4_FROMINT32X2(Y4), AE_MOVINT16X4_FROMINT32X2(Y6), ind);
            AE_DSEL16X4(C4, C5, AE_MOVINT16X4_FROMINT32X2(Y1), AE_MOVINT16X4_FROMINT32X2(Y3), ind);
            AE_DSEL16X4(C6, C7, AE_MOVINT16X4_FROMINT32X2(Y5), AE_MOVINT16X4_FROMINT32X2(Y7), ind);

            AE_S32X2X2_IP(AE_MOVINT32X2_FROMINT16X4(C0), AE_MOVINT32X2_FROMINT16X4(C2), pY0, 2 * sizeof(ae_int32x2));
            AE_S32X2X2_IP(AE_MOVINT32X2_FROMINT16X4(C1), AE_MOVINT32X2_FROMINT16X4(C3), pY1, 2 * sizeof(ae_int32x2));
            AE_S32X2X2_IP(AE_MOVINT32X2_FROMINT16X4(C4), AE_MOVINT32X2_FROMINT16X4(C6), pY2, 2 * sizeof(ae_int32x2));
            AE_S32X2X2_IP(AE_MOVINT32X2_FROMINT16X4(C5), AE_MOVINT32X2_FROMINT16X4(C7), pY3, 2 * sizeof(ae_int32x2));
            b_y=(int32_t*)XT_ADDX8(2*P,(uintptr_t)b_y);
        }

        //handle tail
        if (N & 3)
        {
            ae_int32x2 Y0, Y1, Y2, Y3, Y4, Y5, Y6, Y7;
            ae_int16x4 C0, C1, C2, C3, C4, C5, C6, C7;
            ae_valignx2 aY;
            int P1,P2;

            pY = (const ae_int32x4 *)(b_y);
            AE_LA32X2X2POS_PC(aY,pY);
            AE_LA32X2X2_IC(Y0, Y1, aY, pY);
            P1=P; XT_MOVGEZ(P1,0,1-(N&3));
            P2=P; XT_MOVGEZ(P2,0,2-(N&3));
            pY = (const ae_int32x4 *)XT_ADDX4(P1,(uintptr_t)(b_y));
            AE_LA32X2X2POS_PC(aY,pY);
            AE_LA32X2X2_IC(Y2, Y3, aY, pY);
            pY =(const ae_int32x4 *)XT_ADDX8(P2,(uintptr_t)(b_y));
            AE_LA32X2X2POS_PC(aY,pY);
            AE_LA32X2X2_IC(Y4, Y5, aY, pY);
            Y6=Y7=0;

            AE_DSEL16X4(C0, C1, AE_MOVINT16X4_FROMINT32X2(Y0), AE_MOVINT16X4_FROMINT32X2(Y2), ind);
            AE_DSEL16X4(C2, C3, AE_MOVINT16X4_FROMINT32X2(Y4), AE_MOVINT16X4_FROMINT32X2(Y6), ind);
            AE_DSEL16X4(C4, C5, AE_MOVINT16X4_FROMINT32X2(Y1), AE_MOVINT16X4_FROMINT32X2(Y3), ind);
            AE_DSEL16X4(C6, C7, AE_MOVINT16X4_FROMINT32X2(Y5), AE_MOVINT16X4_FROMINT32X2(Y7), ind);
            Y0 = AE_MOVINT32X2_FROMINT16X4(C0);
            Y1 = AE_MOVINT32X2_FROMINT16X4(C1);
            Y2 = AE_MOVINT32X2_FROMINT16X4(C2);
            Y3 = AE_MOVINT32X2_FROMINT16X4(C3);
            Y4 = AE_MOVINT32X2_FROMINT16X4(C4);
            Y5 = AE_MOVINT32X2_FROMINT16X4(C5);
            Y6 = AE_MOVINT32X2_FROMINT16X4(C6);
            Y7 = AE_MOVINT32X2_FROMINT16X4(C7);

            AE_MOVF32X2(Y0, AE_ZERO32(), bmask0);
            AE_MOVF32X2(Y1, AE_ZERO32(), bmask0);
            AE_MOVF32X2(Y4, AE_ZERO32(), bmask0);
            AE_MOVF32X2(Y5, AE_ZERO32(), bmask0);

            AE_MOVF32X2(Y2, AE_ZERO32(), bmask1);
            AE_MOVF32X2(Y3, AE_ZERO32(), bmask1);
            AE_MOVF32X2(Y6, AE_ZERO32(), bmask1);
            AE_MOVF32X2(Y7, AE_ZERO32(), bmask1);

            AE_S32X2X2_IP(Y0, Y2, pY0, 2 * sizeof(ae_int32x2));
            AE_S32X2X2_IP(Y1, Y3, pY1, 2 * sizeof(ae_int32x2));
            AE_S32X2X2_IP(Y4, Y6, pY2, 2 * sizeof(ae_int32x2));
            AE_S32X2X2_IP(Y5, Y7, pY3, 2 * sizeof(ae_int32x2));
        }

        /* compute 4*M outputs  */
        px = (ae_int32x4 *)x;

        for (m = 0; m < (M >> 2); m++)
        {
            ae_valignx2 aZ;
            ae_int64 q0, q1, q2, q3;
            ae_int64 q4, q5, q6, q7;
            ae_int64 w0, w1, w2, w3;
            ae_int64 w4, w5, w6, w7;
            pY0 = (ae_int32x4 *)Y0_;
            pY1 = (ae_int32x4 *)Y1_;
            pY2 = (ae_int32x4 *)Y2_;
            pY3 = (ae_int32x4 *)Y3_;
            pX = px;
            pX1 = (const ae_int32x4 *)XT_ADDX4(N, (uintptr_t)pX );
            pX2 = (const ae_int32x4 *)XT_ADDX8(N, (uintptr_t)pX );
            pX3 = (const ae_int32x4 *)XT_ADDX4(N, (uintptr_t)pX2);
            px  = (const ae_int32x4 *)XT_ADDX8(N, (uintptr_t)pX2);
            aX  = AE_LA128_PP(pX );
            aX1 = AE_LA128_PP(pX1);
            aX2 = AE_LA128_PP(pX2);
            aX3 = AE_LA128_PP(pX3);
            AE_MOVDX2(q0,q1,0,0);AE_MOVDX2(q2,q3,0,0);
            AE_MOVDX2(q4,q5,0,0);AE_MOVDX2(q6,q7,0,0);
            AE_MOVDX2(w0,w1,0,0);AE_MOVDX2(w2,w3,0,0);
            AE_MOVDX2(w4,w5,0,0);AE_MOVDX2(w6,w7,0,0);
            __Pragma("loop_count min=1")
            for (n = 0; n < (K >> 2); n++)
            {
                ae_int32x2 X0, X1, X2, X3, X4,X5,X6,X7, Y0, Y1, Y2, Y3, Y4, Y5, Y6, Y7;
                AE_LA32X2X2_IP(X0, X1, aX, pX);
                AE_LA32X2X2_IP(X2, X3, aX1, pX1);
                AE_LA32X2X2_IP(X4, X5, aX2, pX2);
                AE_LA32X2X2_IP(X6, X7, aX3, pX3);
                AE_L32X2X2_X (Y2, Y3, pY0, (K1)*sizeof(int32_t));
                AE_L32X2X2_IP(Y0, Y1, pY0, sizeof(ae_int32x4));
                AE_L32X2X2_X (Y6, Y7, pY2, (K1)*sizeof(int32_t));
                AE_L32X2X2_IP(Y4, Y5, pY2, sizeof(ae_int32x4));

                AE_MULAAF2D32RA_HH_LL(q0, q1, X0, X0, Y0, Y2); AE_MULAAF2D32RA_HH_LL(w0, w1, X4, X4, Y0, Y2);
                AE_MULAAF2D32RA_HH_LL(q0, q1, X1, X1, Y1, Y3); AE_MULAAF2D32RA_HH_LL(w0, w1, X5, X5, Y1, Y3);
                AE_MULAAF2D32RA_HH_LL(q2, q3, X0, X0, Y4, Y6); AE_MULAAF2D32RA_HH_LL(w2, w3, X4, X4, Y4, Y6);
                AE_MULAAF2D32RA_HH_LL(q2, q3, X1, X1, Y5, Y7); AE_MULAAF2D32RA_HH_LL(w2, w3, X5, X5, Y5, Y7);
                AE_MULAAF2D32RA_HH_LL(q4, q5, X2, X2, Y0, Y2); AE_MULAAF2D32RA_HH_LL(w4, w5, X6, X6, Y0, Y2);
                AE_MULAAF2D32RA_HH_LL(q4, q5, X3, X3, Y1, Y3); AE_MULAAF2D32RA_HH_LL(w4, w5, X7, X7, Y1, Y3);
                AE_MULAAF2D32RA_HH_LL(q6, q7, X2, X2, Y4, Y6); AE_MULAAF2D32RA_HH_LL(w6, w7, X6, X6, Y4, Y6);
                AE_MULAAF2D32RA_HH_LL(q6, q7, X3, X3, Y5, Y7); AE_MULAAF2D32RA_HH_LL(w6, w7, X7, X7, Y5, Y7);
            }

            aZ=AE_ZALIGN128();
            AE_SA32X2X2_IP(AE_TRUNCA32X2F64S(q0, q1, 16 + lsh),AE_TRUNCA32X2F64S(q2, q3, 16 + lsh),aZ,castxcc(ae_int32x4,pZ));  AE_SA128POS_FP(aZ,pZ);
            pZ=(ae_int32*)(P*sizeof(int32_t)-sizeof(ae_int32x4)+(uintptr_t)pZ);
            AE_SA32X2X2_IP(AE_TRUNCA32X2F64S(q4, q5, 16 + lsh),AE_TRUNCA32X2F64S(q6, q7, 16 + lsh),aZ,castxcc(ae_int32x4,pZ));  AE_SA128POS_FP(aZ,pZ);
            pZ=(ae_int32*)(P*sizeof(int32_t)-sizeof(ae_int32x4)+(uintptr_t)pZ);
            AE_SA32X2X2_IP(AE_TRUNCA32X2F64S(w0, w1, 16 + lsh),AE_TRUNCA32X2F64S(w2, w3, 16 + lsh),aZ,castxcc(ae_int32x4,pZ));  AE_SA128POS_FP(aZ,pZ);
            pZ=(ae_int32*)(P*sizeof(int32_t)-sizeof(ae_int32x4)+(uintptr_t)pZ);
            AE_SA32X2X2_IP(AE_TRUNCA32X2F64S(w4, w5, 16 + lsh),AE_TRUNCA32X2F64S(w6, w7, 16 + lsh),aZ,castxcc(ae_int32x4,pZ));  AE_SA128POS_FP(aZ,pZ);
            pZ=(ae_int32*)(P*sizeof(int32_t)-sizeof(ae_int32x4)+(uintptr_t)pZ);
        }
        if (M&2)
        {
            ae_valignx2 aZ=AE_ZALIGN128();
            ae_int64 q0, q1, q2, q3;
            ae_int64 q4, q5, q6, q7;
            pY0 = (ae_int32x4 *)Y0_;
            pY1 = (ae_int32x4 *)Y1_;
            pY2 = (ae_int32x4 *)Y2_;
            pY3 = (ae_int32x4 *)Y3_;
            pX = px;
            pX1 = (const ae_int32x4 *)XT_ADDX4(N, (uintptr_t)pX);
            px = (const ae_int32x4 *)XT_ADDX8(N, (uintptr_t)px);
            aX = AE_LA128_PP(pX);
            aX1 = AE_LA128_PP(pX1);
            AE_MOVDX2(q0,q1,0,0);AE_MOVDX2(q2,q3,0,0);
            AE_MOVDX2(q4,q5,0,0);AE_MOVDX2(q6,q7,0,0);
            __Pragma("loop_count min=1")
            for (n = 0; n < (K >> 2); n++)
            {
                ae_int32x2 X0, X1, X2, X3, Y0, Y1, Y2, Y3, Y4, Y5, Y6, Y7;
                AE_LA32X2X2_IP(X0, X1, aX, pX);
                AE_LA32X2X2_IP(X2, X3, aX1, pX1);
                AE_L32X2X2_X (Y2, Y3, pY0, (K1)*sizeof(int32_t));
                AE_L32X2X2_IP(Y0, Y1, pY0, sizeof(ae_int32x4));
                AE_L32X2X2_X (Y6, Y7, pY2, (K1)*sizeof(int32_t));
                AE_L32X2X2_IP(Y4, Y5, pY2, sizeof(ae_int32x4));

                AE_MULAAF2D32RA_HH_LL(q0, q1, X0, X0, Y0, Y2);
                AE_MULAAF2D32RA_HH_LL(q0, q1, X1, X1, Y1, Y3);
                AE_MULAAF2D32RA_HH_LL(q2, q3, X0, X0, Y4, Y6);
                AE_MULAAF2D32RA_HH_LL(q4, q5, X2, X2, Y0, Y2);
                AE_MULAAF2D32RA_HH_LL(q6, q7, X2, X2, Y4, Y6);
                AE_MULAAF2D32RA_HH_LL(q2, q3, X1, X1, Y5, Y7);
                AE_MULAAF2D32RA_HH_LL(q4, q5, X3, X3, Y1, Y3);
                AE_MULAAF2D32RA_HH_LL(q6, q7, X3, X3, Y5, Y7);
            }

            AE_SA32X2X2_IP(AE_TRUNCA32X2F64S(q0, q1, 16 + lsh),AE_TRUNCA32X2F64S(q2, q3, 16 + lsh),aZ,castxcc(ae_int32x4,pZ)); 
            AE_SA128POS_FP(aZ,pZ);
            pZ=(ae_int32*)(P*sizeof(int32_t)-sizeof(ae_int32x4)+(uintptr_t)pZ);
            AE_SA32X2X2_IP(AE_TRUNCA32X2F64S(q4, q5, 16 + lsh),AE_TRUNCA32X2F64S(q6, q7, 16 + lsh),aZ,castxcc(ae_int32x4,pZ)); 
            AE_SA128POS_FP(aZ,pZ);
            pZ=(ae_int32*)(P*sizeof(int32_t)-sizeof(ae_int32x4)+(uintptr_t)pZ);
        }
        if (M & 1)
        {
            ae_int32x2 z0, z1;
            ae_int64 q0, q1, q2, q3;
            pY0 = (ae_int32x4 *)Y0_;
            pY1 = (ae_int32x4 *)Y1_;
            pY2 = (ae_int32x4 *)Y2_;
            pY3 = (ae_int32x4 *)Y3_;
            pX = px;
            aX = AE_LA128_PP(pX);
            AE_MOVDX2(q0,q1,0,0); AE_MOVDX2(q2,q3,0,0);
            __Pragma("loop_count min=1")
            for (n = 0; n < (K >> 2); n++)
            {
                ae_int32x2 Y0, Y1, Y2, Y3, Y4, Y5, Y6, Y7, X0, X1;

                AE_LA32X2X2_IP(X0, X1, aX, pX);
                AE_L32X2X2_X (Y2, Y3, pY0, (K1)*sizeof(int32_t));
                AE_L32X2X2_IP(Y0, Y1, pY0, 2 * sizeof(ae_int32x2));
                AE_L32X2X2_X (Y6, Y7, pY2, (K1)*sizeof(int32_t));
                AE_L32X2X2_IP(Y4, Y5, pY2, 2 * sizeof(ae_int32x2));

                AE_MULAAF2D32RA_HH_LL(q0, q1, X0, X0, Y0, Y2);
                AE_MULAAF2D32RA_HH_LL(q0, q1, X1, X1, Y1, Y3);

                AE_MULAAF2D32RA_HH_LL(q2, q3, X0, X0, Y4, Y6);
                AE_MULAAF2D32RA_HH_LL(q2, q3, X1, X1, Y5, Y7);
            }
            z0 = AE_TRUNCA32X2F64S(q0, q1, 16 + lsh);
            z1 = AE_TRUNCA32X2F64S(q2, q3, 16 + lsh);
            AE_S32_L_I (z0, pZ, 1 * sizeof(int32_t));
            AE_S32_H_I (z1, pZ, 2 * sizeof(int32_t));
            AE_S32_L_I (z1, pZ, 3 * sizeof(int32_t));
            AE_S32_H_XP(z0, pZ, P*sizeof(int32_t));
        }
    }
    /* process last 2 columns if any */
    if (P & 2)
    {
        NASSERT_ALIGN(Y0_, HIFI_SIMD_WIDTH);
        NASSERT_ALIGN(Y1_, HIFI_SIMD_WIDTH);
        pZ = (ae_int32*)(z);
        /* copy 2 input columns to the srcatch memory with zero padding */
        pY0 = (ae_int32x4 *)Y0_;
        pY1 = (ae_int32x4 *)Y1_;
        b_y = y;
        for (n = 0; n < (N >> 2); n++)
        {
            ae_int32x2 y0, y1, y2, y3;
            ae_valign aY0,aY1;

            pY_ = (const ae_int32x2 *)(b_y);
            pY  = (const ae_int32x4 *)XT_ADDX2(2*P,(uintptr_t)b_y);
            aY0 = AE_LA64_PP(pY_); AE_LA32X2_IP(y0, aY0, pY_);
            aY1 = AE_LA64_PP(pY ); AE_LA32X2_IP(y1, aY1, castxcc(ae_int32x2,pY ));
            pY_ = (const ae_int32x2 *)XT_ADDX4(2*P,(uintptr_t)b_y);
            pY  = (const ae_int32x4 *)(b_y + 3 * P);
            aY0 = AE_LA64_PP(pY_); AE_LA32X2_IP(y2, aY0, pY_);
            aY1 = AE_LA64_PP(pY ); AE_LA32X2_IP(y3, aY1, castxcc(ae_int32x2,pY ));
            b_y=(int32_t*)XT_ADDX8(2*P,(uintptr_t)b_y);
            AE_S32X2X2_IP(AE_SEL32_HH(y0, y1), AE_SEL32_HH(y2, y3), pY0, 2 * sizeof(ae_int32x2));
            AE_S32X2X2_IP(AE_SEL32_LL(y0, y1), AE_SEL32_LL(y2, y3), pY1, 2 * sizeof(ae_int32x2));
        }
        if (N & 3)
        {
            ae_int32x2 y0, y1, y2, y0_, y1_, y2_, y3_;
            ae_valign aY;
            int P1,P2;

            P1=P; XT_MOVGEZ(P1,0,1-(N&3));
            P2=P; XT_MOVGEZ(P2,0,2-(N&3));
            pY_ = (const ae_int32x2 *)(b_y);
            AE_LA32X2POS_PC(aY,pY_);
            AE_LA32X2_IC(y0, aY, pY_);
            pY_ = (const ae_int32x2 *)XT_ADDX4(P1,(uintptr_t)(b_y));
            AE_LA32X2POS_PC(aY,pY_);
            AE_LA32X2_IC(y1, aY, pY_);
            pY_ = (const ae_int32x2 *)XT_ADDX8(P2,(uintptr_t)(b_y));
            AE_LA32X2POS_PC(aY,pY_);
            AE_LA32X2_IC(y2, aY, pY_);

            y0_ = AE_SEL32_HH(y0, y1);
            y1_ = AE_SEL32_HH(y2,  0);
            y2_ = AE_SEL32_LL(y0, y1);
            y3_ = AE_SEL32_LL(y2,  0);

            AE_MOVF32X2(y0_, AE_ZERO32(), bmask0);
            AE_MOVF32X2(y1_, AE_ZERO32(), bmask1);
            AE_MOVF32X2(y2_, AE_ZERO32(), bmask0);
            AE_MOVF32X2(y3_, AE_ZERO32(), bmask1);

            AE_S32X2X2_IP(y0_, y1_, pY0, 2 * sizeof(ae_int32x2));
            AE_S32X2X2_IP(y2_, y3_, pY1, 2 * sizeof(ae_int32x2));
        }
        /* compute 2*M outputs  */
        px = (ae_int32x4 *)x;
        for (m = 0; m < (M >> 2); m++)
        {
            ae_int64 q0, q1, q2, q3;
            ae_int64 q4, q5, q6, q7;
            ae_int32x2 z0, z1, z2, z3;
            pY0 = (ae_int32x4 *)Y0_;
            pY1 = (ae_int32x4 *)Y1_;
            pX = px;
            pX1 = (const ae_int32x4 *)XT_ADDX4(N, (uintptr_t)pX);
            pX2 = (const ae_int32x4 *)XT_ADDX4(N, (uintptr_t)pX1);
            pX3 = (const ae_int32x4 *)XT_ADDX4(N, (uintptr_t)pX2);
            px  = (const ae_int32x4 *)XT_ADDX4(N, (uintptr_t)pX3);
            aX  = AE_LA128_PP(pX);
            aX1 = AE_LA128_PP(pX1);
            aX2 = AE_LA128_PP(pX2);
            aX3 = AE_LA128_PP(pX3);
            AE_MOVDX2(q0,q1,0,0); AE_MOVDX2(q2,q3,0,0);
            AE_MOVDX2(q4,q5,0,0); AE_MOVDX2(q6,q7,0,0);
            __Pragma("loop_count min=1")
            for (n = 0; n < (K >> 2); n++)
            {
                    ae_int32x2 X0, X1, X2, X3, X4, X5, X6, X7, Y0, Y1, Y2, Y3;
                    AE_LA32X2X2_IP(X0, X1, aX, pX);
                    AE_LA32X2X2_IP(X2, X3, aX1, pX1);
                    AE_LA32X2X2_IP(X4, X5, aX2, pX2);
                    AE_LA32X2X2_IP(X6, X7, aX3, pX3);
                    AE_L32X2X2_X (Y2, Y3, pY0, (K1)*sizeof(int32_t));
                    AE_L32X2X2_IP(Y0, Y1, pY0, 2 * sizeof(ae_int32x2));

                    AE_MULAAF2D32RA_HH_LL(q0, q1, X0, X0, Y0, Y2);
                    AE_MULAAF2D32RA_HH_LL(q0, q1, X1, X1, Y1, Y3);
                    AE_MULAAF2D32RA_HH_LL(q2, q3, X2, X2, Y0, Y2);
                    AE_MULAAF2D32RA_HH_LL(q2, q3, X3, X3, Y1, Y3);
                    AE_MULAAF2D32RA_HH_LL(q4, q5, X4, X4, Y0, Y2);
                    AE_MULAAF2D32RA_HH_LL(q4, q5, X5, X5, Y1, Y3);
                    AE_MULAAF2D32RA_HH_LL(q6, q7, X6, X6, Y0, Y2);
                    AE_MULAAF2D32RA_HH_LL(q6, q7, X7, X7, Y1, Y3);
            }
            z0 = AE_TRUNCA32X2F64S(q0, q1, 16 + lsh);
            z1 = AE_TRUNCA32X2F64S(q2, q3, 16 + lsh);
            z2 = AE_TRUNCA32X2F64S(q4, q5, 16 + lsh);
            z3 = AE_TRUNCA32X2F64S(q6, q7, 16 + lsh);
            AE_S32_L_I (z0, pZ, 1 * sizeof(int32_t));
            AE_S32_H_XP(z0, pZ, P * sizeof(int32_t));
            AE_S32_L_I (z1, pZ, 1 * sizeof(int32_t));
            AE_S32_H_XP(z1, pZ, P * sizeof(int32_t));
            AE_S32_L_I (z2, pZ, 1 * sizeof(int32_t));
            AE_S32_H_XP(z2, pZ, P * sizeof(int32_t));
            AE_S32_L_I (z3, pZ, 1 * sizeof(int32_t));
            AE_S32_H_XP(z3, pZ, P * sizeof(int32_t));
        }
        m <<= 2;
        if (M & 2)
        {
            ae_int64 q0, q1, q2, q3;
            ae_int32x2 z0, z1;
            pY0 = (ae_int32x4 *)Y0_;
            pY1 = (ae_int32x4 *)Y1_;
            pX = px;
            pX1 = (const ae_int32x4 *)XT_ADDX4(N, (uintptr_t)pX);
            px  = (const ae_int32x4 *)XT_ADDX4(N, (uintptr_t)pX1);
            aX  = AE_LA128_PP(pX);
            aX1 = AE_LA128_PP(pX1);
            AE_MOVDX2(q0,q1,0,0);
            AE_MOVDX2(q2,q3,0,0);
            __Pragma("loop_count min=1")
            for (n = 0; n < (K >> 2); n++)
            {
                ae_int32x2 y0, y1, y2, y3, x0, x1, x2, x3;
                AE_LA32X2X2_IP(x0, x1, aX, pX);
                AE_LA32X2X2_IP(x2, x3, aX1, pX1);
                AE_L32X2X2_IP(y0, y1, pY0, 2 * sizeof(ae_int32x2));
                AE_L32X2X2_IP(y2, y3, pY1, 2 * sizeof(ae_int32x2));
                AE_MULAAF2D32RA_HH_LL(q0, q1, x0, x0, y0, y2);
                AE_MULAAF2D32RA_HH_LL(q0, q1, x1, x1, y1, y3);
                AE_MULAAF2D32RA_HH_LL(q2, q3, x2, x2, y0, y2);
                AE_MULAAF2D32RA_HH_LL(q2, q3, x3, x3, y1, y3);
            }
            z0 = AE_TRUNCA32X2F64S(q0, q1, 16 + lsh);
            z1 = AE_TRUNCA32X2F64S(q2, q3, 16 + lsh);
            AE_S32_L_I (z0, pZ, 1 * sizeof(int32_t));
            AE_S32_H_XP(z0, pZ, P * sizeof(int32_t));
            AE_S32_L_I (z1, pZ, 1 * sizeof(int32_t));
            AE_S32_H_XP(z1, pZ, P * sizeof(int32_t));
            m += 2;
        }
        if (M & 1)
        {
            ae_int64 q0, q1;
            ae_int32x2 z0;
            pY0 = (ae_int32x4 *)Y0_;
            pY1 = (ae_int32x4 *)Y1_;
            pX = px;
            aX = AE_LA128_PP(pX);
            AE_MOVDX2(q0,q1,0,0);
            __Pragma("loop_count min=1")
            for (n = 0; n < (K >> 2); n++)
            {
                ae_int32x2 y0, y1, y2, y3, x1, x2;
                AE_LA32X2X2_IP(x1, x2, aX, pX);
                AE_L32X2X2_IP(y0, y1, pY0, 2 * sizeof(ae_int32x2));
                AE_L32X2X2_IP(y2, y3, pY1, 2 * sizeof(ae_int32x2));
                AE_MULAAF2D32RA_HH_LL(q0, q1, x1, x1, y0, y2);
                AE_MULAAF2D32RA_HH_LL(q0, q1, x2, x2, y1, y3);
            }
            z0 = AE_TRUNCA32X2F64S(q0, q1, 16 + lsh);
            AE_S32_L_I (z0, pZ, 1 * sizeof(int32_t));
            AE_S32_H_XP(z0, pZ, P * sizeof(int32_t));
        }
        z += 2, y += 2;  /* go to the next 2 columns */
    }
    /* process last column if any */
    if (P & 1)
    {
        NASSERT_ALIGN(Y0_, HIFI_SIMD_WIDTH);
        NASSERT_ALIGN(Y1_, HIFI_SIMD_WIDTH);
        pZ = (ae_int32*)(z);
        /* copy input column to the srcatch memory with zero padding */
        pY0 = (ae_int32x4 *)Y0_;
        pY_ = (const ae_int32x2 *)(y);
        for (n = 0; n < (N >> 2); n++)
        {
            ae_int32x2 y0, y1, y2, y3;
            AE_L32_XP(y0, castxcc(ae_int32, pY_), P*sizeof(int32_t));
            AE_L32_XP(y1, castxcc(ae_int32, pY_), P*sizeof(int32_t));
            AE_L32_XP(y2, castxcc(ae_int32, pY_), P*sizeof(int32_t));
            AE_L32_XP(y3, castxcc(ae_int32, pY_), P*sizeof(int32_t));

            y0 = AE_SEL32_LL(y0, y1);
            y1 = AE_SEL32_LL(y2, y3);

            AE_S32X2X2_IP(y0, y1, pY0, 2 * sizeof(ae_int32x2));
        }
        if (N & 3)
        {
            ae_int32x2 y0, y1, y2;
            int P1,P2;
            P1=P*sizeof(int32_t);   XT_MOVGEZ(P1,0,1-(N&3));
            P2=2*P*sizeof(int32_t); XT_MOVGEZ(P2,0,2-(N&3));
            y0=AE_L32_I ((const ae_int32*) pY_,  0);
            y1=AE_L32_X ((const ae_int32*) pY_, P1);
            y2=AE_L32_X ((const ae_int32*) pY_, P2);

            y0 = AE_SEL32_HH(y0, y1);
            y1 = AE_SEL32_HH(y2, 0);

            AE_MOVF32X2(y0, AE_ZERO32(), bmask0);
            AE_MOVF32X2(y1, AE_ZERO32(), bmask1);

            AE_S32X2X2_IP(y0, y1, pY0, 2 * sizeof(ae_int32x2));
        }

        /* compute M outputs  */
        px = (ae_int32x4 *)x;
        for (m = 0; m < M; m++)
        {
            ae_int64 q0, q1;
            ae_int32x2 z0;
            pY0 = (ae_int32x4 *)Y0_;
            pX = px;
            px = (const ae_int32x4 *)XT_ADDX4(N, (uintptr_t)px);
            aX = AE_LA128_PP(pX);
            AE_MOVDX2(q0,q1,0,0);
            __Pragma("loop_count min=1")
            for (n = 0; n < (K >> 2); n++)
            {
                ae_int32x2 y0, y1, x0, x1;
                AE_LA32X2X2_IP(x0, x1, aX, pX);
                AE_L32X2X2_IP(y0, y1, pY0, 2 * sizeof(ae_int32x2));
                AE_MULAAFD32RA_HH_LL(q0, x0, y0);
                AE_MULAAFD32RA_HH_LL(q1, x1, y1);
            }
            z0 = AE_TRUNCA32X2F64S(q0 + q1, q0 + q1, 16 + lsh);
            AE_S32_L_XP(z0, pZ, P*sizeof(int32_t));
        }
    }
}

size_t mtx_mpy32x32_getScratchSize(int M, int N, int P)
{
    int K;
    (void)M; (void)P;
    K = ((N + 7)&(~7))+4;

    return N <= 0 ? 0 : (K * 4 * sizeof(int32_t));
}
