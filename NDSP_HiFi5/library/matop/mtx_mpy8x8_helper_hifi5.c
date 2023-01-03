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

/* helper functions for matrix 8x8 multiplies  */
#include "NatureDSP_types.h"
#include "common.h"
#include "mtx_mpy8x8_helper.h"
/* code optimized for HiFi5 core */
#if USE_NN_EXTENSION_8X8==0
static const uint64_t dsel_rephhll_tbl=0x40516273c8d9eafbULL;

/*-------------------------------------------------------------------------
  Matrix Multiply

  multiplies x by y' where y is written with paddings

  input:
  x[M][N]
  y[P][N0]

  output:
  z[M][N0] -- zero padded output

  N0 - closest bigger multiple of 8 of N

  z,y - aligned
-------------------------------------------------------------------------*/
void mtx_mpyt8x8_partiallyaligned (  
                     int8_t* restrict z,
               const int8_t* restrict x,
               const int8_t* restrict y,
               int M, int N, int P, int lsh )
{
    const ae_int8x8 * restrict pX0;
    const ae_int8x8 * restrict pX1;
    const ae_int8x8 * restrict pX2;
    const ae_int8x8 * restrict pX3;
    const ae_int8x8 * restrict pY0;
    const ae_int8x8 * restrict pY2;
          ae_int8   * restrict pZ;
    int M0,N0,P0;
    int sa,m,n,p;
    ae_int8x8 dsel_rephhll;

    NASSERT(lsh >= -15 && lsh <= 15);
    dsel_rephhll=AE_L8X8_I((const ae_int8x8*)&dsel_rephhll_tbl,0);
    sa=49+lsh;
    M0=M&~3;
    N0=((N+7)&~7);
    P0=((P+3)&~3);
    __Pragma("no_reorder");
    pZ=(ae_int8 *)z;
    if (N0==8)
    {
        // fully unrolled loop for small N0 and P0
        pX0=(const ae_int8x8 *)&x[0*N];
        pX1=(const ae_int8x8 *)&x[1*N];
        pX2=(const ae_int8x8 *)XT_ADDX2(N,(uintptr_t)pX0);
        pX3=(const ae_int8x8 *)XT_ADDX2(N,(uintptr_t)pX1);
        if (P0==4)
        {
            ae_int8x8 x04,x26,x15,x37,y0,y1,y2,y3,x0,x1,x2,x3,x4,x5,x6,x7,y01,y23,y45,y67;
            ae_int8x8 r0,r1,r2,r3;
            ae_valign aX0,aX1,aX2,aX3;
            ae_int64 a00,a01,a02,a03,a10,a11,a12,a13;
            ae_int64 a20,a21,a22,a23,a30,a31,a32,a33;
            pY0=(const ae_int8x8 *)y;
            pY2=(const ae_int8x8 *)XT_ADDX2(N0,(uintptr_t)y);
            y0= AE_L8X8_I  (pY0,0);
            y1= AE_L8X8_X  (pY0,N0);
            y2= AE_L8X8_I  (pY2,0);
            y3= AE_L8X8_X  (pY2,N0);

            AE_DSEL8X8(y01,y45,y0,y1,dsel_rephhll);
            AE_DSEL8X8(y23,y67,y2,y3,dsel_rephhll);
            for (m=0; m<(M0>>2);m++)
            {
                aX0=AE_LA64_PP(pX0);
                aX1=AE_LA64_PP(pX1);
                aX2=AE_LA64_PP(pX2);
                aX3=AE_LA64_PP(pX3);
                AE_LA8X8_IP(x04,aX0,pX0 );
                AE_LA8X8_IP(x15,aX1,pX1 );
                AE_LA8X8_IP(x26,aX2,pX2 );
                AE_LA8X8_IP(x37,aX3,pX3 );
                pX0=(const ae_int8x8*)((uintptr_t)pX0+(4*N-8));
                pX1=(const ae_int8x8*)((uintptr_t)pX1+(4*N-8));
                pX2=(const ae_int8x8*)((uintptr_t)pX2+(4*N-8));
                pX3=(const ae_int8x8*)((uintptr_t)pX3+(4*N-8));
                AE_DSEL8X8(x0,x4,x04,x04,dsel_rephhll);
                AE_DSEL8X8(x2,x6,x26,x26,dsel_rephhll);
                AE_DSEL8X8(x1,x5,x15,x15,dsel_rephhll);
                AE_DSEL8X8(x3,x7,x37,x37,dsel_rephhll);

                AE_MULZAAAA2Q8(a00,a01,x0,y01);
                AE_MULZAAAA2Q8(a02,a03,x0,y23);
                AE_MULZAAAA2Q8(a10,a11,x1,y01);
                AE_MULZAAAA2Q8(a12,a13,x1,y23);
                AE_MULZAAAA2Q8(a20,a21,x2,y01);
                AE_MULZAAAA2Q8(a22,a23,x2,y23);
                AE_MULZAAAA2Q8(a30,a31,x3,y01);
                AE_MULZAAAA2Q8(a32,a33,x3,y23);
                AE_MULAAAA2Q8(a00,a01,x4,y45);
                AE_MULAAAA2Q8(a02,a03,x4,y67);
                AE_MULAAAA2Q8(a10,a11,x5,y45);
                AE_MULAAAA2Q8(a12,a13,x5,y67);
                AE_MULAAAA2Q8(a20,a21,x6,y45);
                AE_MULAAAA2Q8(a22,a23,x6,y67);
                AE_MULAAAA2Q8(a30,a31,x7,y45);
                AE_MULAAAA2Q8(a32,a33,x7,y67);
                r0=AE_ROUND8X4F32SASYM_L(AE_TRUNCA32X2F64S(a03,a02, sa),AE_TRUNCA32X2F64S(a01,a00, sa));
                r1=AE_ROUND8X4F32SASYM_L(AE_TRUNCA32X2F64S(a13,a12, sa),AE_TRUNCA32X2F64S(a11,a10, sa));
                r2=AE_ROUND8X4F32SASYM_L(AE_TRUNCA32X2F64S(a23,a22, sa),AE_TRUNCA32X2F64S(a21,a20, sa));
                r3=AE_ROUND8X4F32SASYM_L(AE_TRUNCA32X2F64S(a33,a32, sa),AE_TRUNCA32X2F64S(a31,a30, sa));
                AE_S32_L_X (AE_MOVINT32X2_FROMINT8X8(r3),(ae_int32*)pZ       ,3*P0);
                AE_S32_L_X (AE_MOVINT32X2_FROMINT8X8(r2),(ae_int32*)pZ       ,2*P0);
                AE_S32_L_X (AE_MOVINT32X2_FROMINT8X8(r1),(ae_int32*)pZ       ,1*P0);
                AE_S32_L_XP(AE_MOVINT32X2_FROMINT8X8(r0),castxcc(ae_int32,pZ),4+3*P0);
            }
        }
        else
        {
            for (m=0; m<(M0>>2);m++)
            {
                ae_int8x8 x04,x26,x15,x37,y0,y1,y2,y3,x0,x1,x2,x3,x4,x5,x6,x7,y01,y23,y45,y67;
                ae_int8x8 r0,r1,r2,r3;
                ae_valign aX0,aX1,aX2,aX3;
                ae_int64 a00,a01,a02,a03,a10,a11,a12,a13;
                ae_int64 a20,a21,a22,a23,a30,a31,a32,a33;
                pY0=(const ae_int8x8 *)y;
                pY2=(const ae_int8x8 *)XT_ADDX2(N0,(uintptr_t)y);
                aX0=AE_LA64_PP(pX0);
                aX1=AE_LA64_PP(pX1);
                aX2=AE_LA64_PP(pX2);
                aX3=AE_LA64_PP(pX3);
                AE_LA8X8_IP(x04,aX0,pX0 );
                AE_LA8X8_IP(x15,aX1,pX1 );
                AE_LA8X8_IP(x26,aX2,pX2 );
                AE_LA8X8_IP(x37,aX3,pX3 );
                pX0=(const ae_int8x8*)((uintptr_t)pX0+(4*N-8));
                pX1=(const ae_int8x8*)((uintptr_t)pX1+(4*N-8));
                pX2=(const ae_int8x8*)((uintptr_t)pX2+(4*N-8));
                pX3=(const ae_int8x8*)((uintptr_t)pX3+(4*N-8));
                AE_DSEL8X8(x0,x4,x04,x04,dsel_rephhll);
                AE_DSEL8X8(x2,x6,x26,x26,dsel_rephhll);
                AE_DSEL8X8(x1,x5,x15,x15,dsel_rephhll);
                AE_DSEL8X8(x3,x7,x37,x37,dsel_rephhll);
                __Pragma("loop_count min=1")
                for (p=0; p<(P0>>2); p++)
                {
                    y1= AE_L8X8_X  (pY0,N0);
                    AE_L8X8_XP (y0, pY0,sizeof(ae_int8x8)+3*N0);
                    y3=AE_L8X8_X   (pY2,N0);
                    AE_L8X8_XP (y2, pY2,sizeof(ae_int8x8)+3*N0);

                    AE_DSEL8X8(y01,y45,y0,y1,dsel_rephhll);
                    AE_DSEL8X8(y23,y67,y2,y3,dsel_rephhll);

                    AE_MULZAAAA2Q8(a00,a01,x0,y01);
                    AE_MULZAAAA2Q8(a02,a03,x0,y23);
                    AE_MULZAAAA2Q8(a10,a11,x1,y01);
                    AE_MULZAAAA2Q8(a12,a13,x1,y23);
                    AE_MULZAAAA2Q8(a20,a21,x2,y01);
                    AE_MULZAAAA2Q8(a22,a23,x2,y23);
                    AE_MULZAAAA2Q8(a30,a31,x3,y01);
                    AE_MULZAAAA2Q8(a32,a33,x3,y23);
                    AE_MULAAAA2Q8(a00,a01,x4,y45);
                    AE_MULAAAA2Q8(a02,a03,x4,y67);
                    AE_MULAAAA2Q8(a10,a11,x5,y45);
                    AE_MULAAAA2Q8(a12,a13,x5,y67);
                    AE_MULAAAA2Q8(a20,a21,x6,y45);
                    AE_MULAAAA2Q8(a22,a23,x6,y67);
                    AE_MULAAAA2Q8(a30,a31,x7,y45);
                    AE_MULAAAA2Q8(a32,a33,x7,y67);
                    r0=AE_ROUND8X4F32SASYM_L(AE_TRUNCA32X2F64S(a03,a02, sa),AE_TRUNCA32X2F64S(a01,a00, sa));
                    r1=AE_ROUND8X4F32SASYM_L(AE_TRUNCA32X2F64S(a13,a12, sa),AE_TRUNCA32X2F64S(a11,a10, sa));
                    r2=AE_ROUND8X4F32SASYM_L(AE_TRUNCA32X2F64S(a23,a22, sa),AE_TRUNCA32X2F64S(a21,a20, sa));
                    r3=AE_ROUND8X4F32SASYM_L(AE_TRUNCA32X2F64S(a33,a32, sa),AE_TRUNCA32X2F64S(a31,a30, sa));
                    AE_S32_L_X (AE_MOVINT32X2_FROMINT8X8(r3),(ae_int32*)pZ ,3*P0);
                    AE_S32_L_X (AE_MOVINT32X2_FROMINT8X8(r2),(ae_int32*)pZ ,2*P0);
                    AE_S32_L_X (AE_MOVINT32X2_FROMINT8X8(r1),(ae_int32*)pZ ,1*P0);
                    AE_S32_L_IP(AE_MOVINT32X2_FROMINT8X8(r0),castxcc(ae_int32,pZ),4);
                }
                pZ+=3*P0;
            }
        }
    }
    else
    {
        for (m=0; m<M0;m+=4)
        {
            pY0=(const ae_int8x8 *)y;
            pY2=(const ae_int8x8 *)XT_ADDX2(N0,(uintptr_t)y);
            for (p=0; p<P0; p+=4)
            {
                ae_int8x8 r0,r1,r2,r3;
                ae_valign aX0,aX1,aX2,aX3;
                ae_int64 a00,a01,a02,a03,a10,a11,a12,a13;
                ae_int64 a20,a21,a22,a23,a30,a31,a32,a33;
                AE_MOVDX2(a00,a01,0,0); AE_MOVDX2(a02,a03,0,0); AE_MOVDX2(a10,a11,0,0); AE_MOVDX2(a12,a13,0,0); 
                AE_MOVDX2(a20,a21,0,0); AE_MOVDX2(a22,a23,0,0); AE_MOVDX2(a30,a31,0,0); AE_MOVDX2(a32,a33,0,0); 

                pX0=(const ae_int8x8 *)&x[(m+0)*N]; aX0=AE_LA64_PP(pX0);
                pX1=(const ae_int8x8 *)&x[(m+1)*N]; aX1=AE_LA64_PP(pX1);
                pX2=(const ae_int8x8 *)&x[(m+2)*N]; aX2=AE_LA64_PP(pX2);
                pX3=(const ae_int8x8 *)&x[(m+3)*N]; aX3=AE_LA64_PP(pX3);
                __Pragma("loop_count min=1")
                for (n=0; n<(N0>>3); n++)
                {
                    ae_int8x8 x04,x26,x15,x37,y0,y1,y2,y3,x0,x1,x2,x3,x4,x5,x6,x7,y01,y23,y45,y67;
                    AE_LA8X8_IP(x04,aX0,pX0 );
                    AE_LA8X8_IP(x15,aX1,pX1 );
                    AE_LA8X8_IP(x26,aX2,pX2 );
                    AE_LA8X8_IP(x37,aX3,pX3 );
                    y1= AE_L8X8_X  (pY0,N0);
                    AE_L8X8_IP (y0, pY0,sizeof(ae_int8x8));
                    y3=AE_L8X8_X   (pY2,N0);
                    AE_L8X8_IP (y2, pY2,sizeof(ae_int8x8));

                    AE_DSEL8X8(x0,x4,x04,x04,dsel_rephhll);
                    AE_DSEL8X8(x2,x6,x26,x26,dsel_rephhll);
                    AE_DSEL8X8(x1,x5,x15,x15,dsel_rephhll);
                    AE_DSEL8X8(x3,x7,x37,x37,dsel_rephhll);
                    AE_DSEL8X8(y01,y45,y0,y1,dsel_rephhll);
                    AE_DSEL8X8(y23,y67,y2,y3,dsel_rephhll);

                    AE_MULAAAA2Q8(a00,a01,x0,y01);
                    AE_MULAAAA2Q8(a02,a03,x0,y23);
                    AE_MULAAAA2Q8(a10,a11,x1,y01);
                    AE_MULAAAA2Q8(a12,a13,x1,y23);
                    AE_MULAAAA2Q8(a20,a21,x2,y01);
                    AE_MULAAAA2Q8(a22,a23,x2,y23);
                    AE_MULAAAA2Q8(a30,a31,x3,y01);
                    AE_MULAAAA2Q8(a32,a33,x3,y23);
                    AE_MULAAAA2Q8(a00,a01,x4,y45);
                    AE_MULAAAA2Q8(a02,a03,x4,y67);
                    AE_MULAAAA2Q8(a10,a11,x5,y45);
                    AE_MULAAAA2Q8(a12,a13,x5,y67);
                    AE_MULAAAA2Q8(a20,a21,x6,y45);
                    AE_MULAAAA2Q8(a22,a23,x6,y67);
                    AE_MULAAAA2Q8(a30,a31,x7,y45);
                    AE_MULAAAA2Q8(a32,a33,x7,y67);
                }
                r0=AE_ROUND8X4F32SASYM_L(AE_TRUNCA32X2F64S(a03,a02, sa),AE_TRUNCA32X2F64S(a01,a00, sa));
                r1=AE_ROUND8X4F32SASYM_L(AE_TRUNCA32X2F64S(a13,a12, sa),AE_TRUNCA32X2F64S(a11,a10, sa));
                r2=AE_ROUND8X4F32SASYM_L(AE_TRUNCA32X2F64S(a23,a22, sa),AE_TRUNCA32X2F64S(a21,a20, sa));
                r3=AE_ROUND8X4F32SASYM_L(AE_TRUNCA32X2F64S(a33,a32, sa),AE_TRUNCA32X2F64S(a31,a30, sa));
                AE_S32_L_X (AE_MOVINT32X2_FROMINT8X8(r3),(ae_int32*)pZ       ,3*P0);
                AE_S32_L_X (AE_MOVINT32X2_FROMINT8X8(r2),(ae_int32*)pZ       ,2*P0);
                AE_S32_L_X (AE_MOVINT32X2_FROMINT8X8(r1),(ae_int32*)pZ       ,1*P0);
                AE_S32_L_IP(AE_MOVINT32X2_FROMINT8X8(r0),castxcc(ae_int32,pZ),4);
                pY0+=3*(N0>>3);
                pY2+=3*(N0>>3);
            }
            pZ+=3*P0;
        }
    }
    __Pragma("no_reorder")
    /* process last rows if number of rows is not a multiple of 4 
       in this loop, the last elements are taken by vector loads with 
       backward offset
    */
    if (M&3)
    switch(M&3)
    {
    case 2:
        {
            m=M-4;
            pZ=(ae_int8*)(z+m*P0);
            pY0=(const ae_int8x8 *)y;
            pY2=(const ae_int8x8 *)XT_ADDX2(N0,(uintptr_t)y);
            for (p=0; p<P; p+=4)
            {
                ae_int8x8 r2,r3;
                ae_valign aX2,aX3;
                ae_int64 a20,a21,a22,a23,a30,a31,a32,a33;
                AE_MOVDX2(a20,a21,0,0); AE_MOVDX2(a22,a23,0,0); AE_MOVDX2(a30,a31,0,0); AE_MOVDX2(a32,a33,0,0); 

                pX2=(const ae_int8x8 *)&x[(m+2)*N]; aX2=AE_LA64_PP(pX2);
                pX3=(const ae_int8x8 *)&x[(m+3)*N]; aX3=AE_LA64_PP(pX3);
                for (n=0; n<(N>>3); n++)
                {
                    ae_int8x8 x26,x37,y0,y1,y2,y3,x2,x3,x6,x7,y01,y23,y45,y67;
                    AE_LA8X8_IP(x26,aX2,pX2 );
                    AE_LA8X8_IP(x37,aX3,pX3 );
                    y1= AE_L8X8_X  (pY0,N0);
                    AE_L8X8_IP (y0, pY0,sizeof(ae_int8x8));
                    y3=AE_L8X8_X   (pY2,N0);
                    AE_L8X8_IP (y2, pY2,sizeof(ae_int8x8));

                    AE_DSEL8X8(x2,x6,x26,x26,dsel_rephhll);
                    AE_DSEL8X8(x3,x7,x37,x37,dsel_rephhll);
                    AE_DSEL8X8(y01,y45,y0,y1,dsel_rephhll);
                    AE_DSEL8X8(y23,y67,y2,y3,dsel_rephhll);

                    AE_MULAAAA2Q8(a20,a21,x2,y01);
                    AE_MULAAAA2Q8(a22,a23,x2,y23);
                    AE_MULAAAA2Q8(a30,a31,x3,y01);
                    AE_MULAAAA2Q8(a32,a33,x3,y23);
                    AE_MULAAAA2Q8(a20,a21,x6,y45);
                    AE_MULAAAA2Q8(a22,a23,x6,y67);
                    AE_MULAAAA2Q8(a30,a31,x7,y45);
                    AE_MULAAAA2Q8(a32,a33,x7,y67);
                }           
                if (N&7)
                {
                    int off=8-(N&7);
                    ae_int8x8 x26,x37,y0,y1,y2,y3,x2,x3,x6,x7,y01,y23,y45,y67;
                    /* take the last elements in two steps:
                       - move back pointers
                       - execute vector load
                       - rotate vectors with zero padding 
                    */
                    pX2=(const ae_int8x8*)(((uintptr_t)pX2)-off); aX2=AE_LA64_PP(pX2);
                    pX3=(const ae_int8x8*)(((uintptr_t)pX3)-off); aX3=AE_LA64_PP(pX3);
                    AE_LA8X8_IP(x26,aX2,pX2 );
                    AE_LA8X8_IP(x37,aX3,pX3 );

                    x26=AE_MOVINT8X8_FROMINT64(AE_SLAA64(AE_MOVINT64_FROMINT8X8(x26), off<<3));
                    x37=AE_MOVINT8X8_FROMINT64(AE_SLAA64(AE_MOVINT64_FROMINT8X8(x37), off<<3));

                    y1= AE_L8X8_X  (pY0,N0);
                    AE_L8X8_IP (y0, pY0,sizeof(ae_int8x8));
                    y3=AE_L8X8_X   (pY2,N0);
                    AE_L8X8_IP (y2, pY2,sizeof(ae_int8x8));

                    AE_DSEL8X8(x2,x6,x26,x26,dsel_rephhll);
                    AE_DSEL8X8(x3,x7,x37,x37,dsel_rephhll);
                    AE_DSEL8X8(y01,y45,y0,y1,dsel_rephhll);
                    AE_DSEL8X8(y23,y67,y2,y3,dsel_rephhll);

                    AE_MULAAAA2Q8(a20,a21,x2,y01);
                    AE_MULAAAA2Q8(a22,a23,x2,y23);
                    AE_MULAAAA2Q8(a30,a31,x3,y01);
                    AE_MULAAAA2Q8(a32,a33,x3,y23);
                    AE_MULAAAA2Q8(a20,a21,x6,y45);
                    AE_MULAAAA2Q8(a22,a23,x6,y67);
                    AE_MULAAAA2Q8(a30,a31,x7,y45);
                    AE_MULAAAA2Q8(a32,a33,x7,y67);
                }
                r2=AE_ROUND8X4F32SASYM_L(AE_TRUNCA32X2F64S(a23,a22, sa),AE_TRUNCA32X2F64S(a21,a20, sa));
                r3=AE_ROUND8X4F32SASYM_L(AE_TRUNCA32X2F64S(a33,a32, sa),AE_TRUNCA32X2F64S(a31,a30, sa));
                AE_S32_L_X (AE_MOVINT32X2_FROMINT8X8(r2),(ae_int32*)pZ  ,2*P0);
                AE_S32_L_X (AE_MOVINT32X2_FROMINT8X8(r3),(ae_int32*)pZ  ,3*P0);
                pZ+=4;
                pY0+=3*(N0>>3);
                pY2+=3*(N0>>3);
            }
        }
        break;
    case 1:
        {
            m=M-4;
            pZ=(ae_int8*)(z+m*P0);
            pY0=(const ae_int8x8 *)y;
            pY2=(const ae_int8x8 *)XT_ADDX2(N0,(uintptr_t)y);
            for (p=0; p<P; p+=4)
            {
                ae_int8x8 r3;
                ae_valign aX3;
                ae_int64 a30,a31,a32,a33;
                AE_MOVDX2(a30,a31,0,0); AE_MOVDX2(a32,a33,0,0); 
                pX3=(const ae_int8x8 *)&x[(m+3)*N]; aX3=AE_LA64_PP(pX3);
                for (n=0; n<(N>>3); n++)
                {
                    ae_int8x8 x37,y0,y1,y2,y3,x3,x7,y01,y23,y45,y67;
                    AE_LA8X8_IP(x37,aX3,pX3 );
                    y1= AE_L8X8_X  (pY0,N0);
                    AE_L8X8_IP (y0, pY0,sizeof(ae_int8x8));
                    y3=AE_L8X8_X   (pY2,N0);
                    AE_L8X8_IP (y2, pY2,sizeof(ae_int8x8));

                    AE_DSEL8X8(x3,x7,x37,x37,dsel_rephhll);
                    AE_DSEL8X8(y01,y45,y0,y1,dsel_rephhll);
                    AE_DSEL8X8(y23,y67,y2,y3,dsel_rephhll);

                    AE_MULAAAA2Q8(a30,a31,x3,y01);
                    AE_MULAAAA2Q8(a32,a33,x3,y23);
                    AE_MULAAAA2Q8(a30,a31,x7,y45);
                    AE_MULAAAA2Q8(a32,a33,x7,y67);
                }           
                if (N&7)
                {
                    int off=8-(N&7);
                    ae_int8x8 x37,y0,y1,y2,y3,x3,x7,y01,y23,y45,y67;
                    /* take the last elements in two steps:
                       - move back pointers
                       - execute vector load
                       - rotate vectors with zero padding 
                    */
                    pX3=(const ae_int8x8*)(((uintptr_t)pX3)-off); aX3=AE_LA64_PP(pX3);
                    AE_LA8X8_IP(x37,aX3,pX3 );

                    x37=AE_MOVINT8X8_FROMINT64(AE_SLAA64(AE_MOVINT64_FROMINT8X8(x37), off<<3));

                    y1= AE_L8X8_X  (pY0,N0);
                    AE_L8X8_IP (y0, pY0,sizeof(ae_int8x8));
                    y3=AE_L8X8_X   (pY2,N0);
                    AE_L8X8_IP (y2, pY2,sizeof(ae_int8x8));

                    AE_DSEL8X8(x3,x7,x37,x37,dsel_rephhll);
                    AE_DSEL8X8(y01,y45,y0,y1,dsel_rephhll);
                    AE_DSEL8X8(y23,y67,y2,y3,dsel_rephhll);

                    AE_MULAAAA2Q8(a30,a31,x3,y01);
                    AE_MULAAAA2Q8(a32,a33,x3,y23);
                    AE_MULAAAA2Q8(a30,a31,x7,y45);
                    AE_MULAAAA2Q8(a32,a33,x7,y67);
                }
                r3=AE_ROUND8X4F32SASYM_L(AE_TRUNCA32X2F64S(a33,a32, sa),AE_TRUNCA32X2F64S(a31,a30, sa));
                AE_S32_L_X (AE_MOVINT32X2_FROMINT8X8(r3),(ae_int32*)pZ  ,3*P0);
                pZ+=4;
                pY0+=3*(N0>>3);
                pY2+=3*(N0>>3);
            }
        }
        break;
    case 3:
        {
            m=M-4;
            pZ=(ae_int8*)(z+m*P0);
            pY0=(const ae_int8x8 *)y;
            pY2=(const ae_int8x8 *)XT_ADDX2(N0,(uintptr_t)y);
            for (p=0; p<P; p+=4)
            {
                ae_int8x8 r1,r2,r3;
                ae_valign aX1,aX2,aX3;
                ae_int64                 a10,a11,a12,a13;
                ae_int64 a20,a21,a22,a23,a30,a31,a32,a33;
                                                                AE_MOVDX2(a10,a11,0,0); AE_MOVDX2(a12,a13,0,0); 
                AE_MOVDX2(a20,a21,0,0); AE_MOVDX2(a22,a23,0,0); AE_MOVDX2(a30,a31,0,0); AE_MOVDX2(a32,a33,0,0); 

                pX1=(const ae_int8x8 *)&x[(m+1)*N]; aX1=AE_LA64_PP(pX1);
                pX2=(const ae_int8x8 *)&x[(m+2)*N]; aX2=AE_LA64_PP(pX2);
                pX3=(const ae_int8x8 *)&x[(m+3)*N]; aX3=AE_LA64_PP(pX3);
                for (n=0; n<(N>>3); n++)
                {
                    ae_int8x8 x26,x15,x37,y0,y1,y2,y3,x1,x2,x3,x5,x6,x7,y01,y23,y45,y67;
                    AE_LA8X8_IP(x15,aX1,pX1 );
                    AE_LA8X8_IP(x26,aX2,pX2 );
                    AE_LA8X8_IP(x37,aX3,pX3 );
                    y1= AE_L8X8_X  (pY0,N0);
                    AE_L8X8_IP (y0, pY0,sizeof(ae_int8x8));
                    y3=AE_L8X8_X   (pY2,N0);
                    AE_L8X8_IP (y2, pY2,sizeof(ae_int8x8));

                    AE_DSEL8X8(x2,x6,x26,x26,dsel_rephhll);
                    AE_DSEL8X8(x1,x5,x15,x15,dsel_rephhll);
                    AE_DSEL8X8(x3,x7,x37,x37,dsel_rephhll);
                    AE_DSEL8X8(y01,y45,y0,y1,dsel_rephhll);
                    AE_DSEL8X8(y23,y67,y2,y3,dsel_rephhll);

                    AE_MULAAAA2Q8(a10,a11,x1,y01);
                    AE_MULAAAA2Q8(a12,a13,x1,y23);
                    AE_MULAAAA2Q8(a20,a21,x2,y01);
                    AE_MULAAAA2Q8(a22,a23,x2,y23);
                    AE_MULAAAA2Q8(a30,a31,x3,y01);
                    AE_MULAAAA2Q8(a32,a33,x3,y23);
                    AE_MULAAAA2Q8(a10,a11,x5,y45);
                    AE_MULAAAA2Q8(a12,a13,x5,y67);
                    AE_MULAAAA2Q8(a20,a21,x6,y45);
                    AE_MULAAAA2Q8(a22,a23,x6,y67);
                    AE_MULAAAA2Q8(a30,a31,x7,y45);
                    AE_MULAAAA2Q8(a32,a33,x7,y67);
                }           
                if (N&7)
                {
                    int off=8-(N&7);
                    ae_int8x8 x26,x15,x37,y0,y1,y2,y3,x1,x2,x3,x5,x6,x7,y01,y23,y45,y67;
                    /* take the last elements in two steps:
                       - move back pointers
                       - execute vector load
                       - rotate vectors with zero padding 
                    */
                    pX1=(const ae_int8x8*)(((uintptr_t)pX1)-off); aX1=AE_LA64_PP(pX1);
                    pX2=(const ae_int8x8*)(((uintptr_t)pX2)-off); aX2=AE_LA64_PP(pX2);
                    pX3=(const ae_int8x8*)(((uintptr_t)pX3)-off); aX3=AE_LA64_PP(pX3);
                    AE_LA8X8_IP(x15,aX1,pX1 );
                    AE_LA8X8_IP(x26,aX2,pX2 );
                    AE_LA8X8_IP(x37,aX3,pX3 );

                    x15=AE_MOVINT8X8_FROMINT64(AE_SLAA64(AE_MOVINT64_FROMINT8X8(x15), off<<3));
                    x26=AE_MOVINT8X8_FROMINT64(AE_SLAA64(AE_MOVINT64_FROMINT8X8(x26), off<<3));
                    x37=AE_MOVINT8X8_FROMINT64(AE_SLAA64(AE_MOVINT64_FROMINT8X8(x37), off<<3));

                    y1= AE_L8X8_X  (pY0,N0);
                    AE_L8X8_IP (y0, pY0,sizeof(ae_int8x8));
                    y3=AE_L8X8_X   (pY2,N0);
                    AE_L8X8_IP (y2, pY2,sizeof(ae_int8x8));

                    AE_DSEL8X8(x2,x6,x26,x26,dsel_rephhll);
                    AE_DSEL8X8(x1,x5,x15,x15,dsel_rephhll);
                    AE_DSEL8X8(x3,x7,x37,x37,dsel_rephhll);
                    AE_DSEL8X8(y01,y45,y0,y1,dsel_rephhll);
                    AE_DSEL8X8(y23,y67,y2,y3,dsel_rephhll);

                    AE_MULAAAA2Q8(a10,a11,x1,y01);
                    AE_MULAAAA2Q8(a12,a13,x1,y23);
                    AE_MULAAAA2Q8(a20,a21,x2,y01);
                    AE_MULAAAA2Q8(a22,a23,x2,y23);
                    AE_MULAAAA2Q8(a30,a31,x3,y01);
                    AE_MULAAAA2Q8(a32,a33,x3,y23);
                    AE_MULAAAA2Q8(a10,a11,x5,y45);
                    AE_MULAAAA2Q8(a12,a13,x5,y67);
                    AE_MULAAAA2Q8(a20,a21,x6,y45);
                    AE_MULAAAA2Q8(a22,a23,x6,y67);
                    AE_MULAAAA2Q8(a30,a31,x7,y45);
                    AE_MULAAAA2Q8(a32,a33,x7,y67);
                }
                r1=AE_ROUND8X4F32SASYM_L(AE_TRUNCA32X2F64S(a13,a12, sa),AE_TRUNCA32X2F64S(a11,a10, sa));
                r2=AE_ROUND8X4F32SASYM_L(AE_TRUNCA32X2F64S(a23,a22, sa),AE_TRUNCA32X2F64S(a21,a20, sa));
                r3=AE_ROUND8X4F32SASYM_L(AE_TRUNCA32X2F64S(a33,a32, sa),AE_TRUNCA32X2F64S(a31,a30, sa));
                AE_S32_L_X (AE_MOVINT32X2_FROMINT8X8(r1),(ae_int32*)pZ  ,1*P0);
                AE_S32_L_X (AE_MOVINT32X2_FROMINT8X8(r2),(ae_int32*)pZ  ,2*P0);
                AE_S32_L_X (AE_MOVINT32X2_FROMINT8X8(r3),(ae_int32*)pZ  ,3*P0);
                pZ+=4;
                pY0+=3*(N0>>3);
                pY2+=3*(N0>>3);
            }
        }
        break;
    }
}

/*-------------------------------------------------------------------------
copy small x to the r with padding

input:
x[M][N]
Output:
r[M][8]

r - aligned
N<8
-------------------------------------------------------------------------*/
void mtx_mpy8x8_copysmallx(int8_t * restrict r,const int8_t * restrict x,int M,int N)
{
    int m,n;
    ae_int8* pZ;
    pZ=(ae_int8*)r;
    NASSERT(N>=1 && N<8);
    NASSERT_ALIGN8(r);
    for(m=0; m<M; m++) 
    {
        AE_S32X2_I(0,(ae_int32x2*)pZ,0);
        __Pragma("loop_count min=1,max=7")
        for(n=0; n<N; n++) 
        {
            ae_int8x8 t;
            AE_L8_IP(t,castxcc(ae_int8,x),1);
            AE_S8_0_IP(t,pZ,1);
        }
        pZ+=8-N;
    }
}


/*-------------------------------------------------------------------------
copy from v to z removing the padding

input:
v[M0][P0]
Output:
z[M][P]

v - aligned
-------------------------------------------------------------------------*/
void mtx_mpyt8x8_copyz(  
                     int8_t* restrict z,
               const int8_t* restrict v,
               int M, int P )
#if 0
{
    int m,p,P0;
    P0=P&~3;
    NASSERT(P0>0);
    for (m=0; m<M;m++)
    {
        for (p=0; p<P; p++)
        {
            z[m*P+p]=v[m*P0+p];
        }
    }
}
#else
{
    ae_valign aZ,aV;
    ae_int8 * restrict pZ;
    int m,p,P0;
    P0=(P+3)&~3;
    NASSERT(P0>0);
    pZ=(ae_int8*)v;
    aZ=AE_ZALIGN64();
    for (m=0; m<M; m++) 
    {
        ae_int8x8 x;
        pZ=(ae_int8*)(v+m*P0);
        aV=AE_LA64_PP(pZ);
        for (p=0; p<(P>>3); p++) 
        {
            AE_LA8X8_IP(x,aV,castxcc(ae_int8x8,pZ));
            AE_SA8X8_IP(x,aZ,castxcc(ae_int8x8,z ));
        }
        AE_SA64POS_FP(aZ,z);
        for (p=0; p<(P&7); p++) 
        {
            AE_L8_IP(x,castxcc(ae_int8,pZ),1);
            AE_S8_0_IP(x,castxcc(ae_int8,z ),1);
        }
    }
}
#endif
#endif /* USE_NN_EXTENSION */
