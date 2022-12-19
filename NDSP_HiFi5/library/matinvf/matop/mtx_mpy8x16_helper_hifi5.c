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

/* helper functions for matrix 8x16 multiplies  */
#include "NatureDSP_types.h"
#include "common.h"
#include "mtx_mpy8x16_helper.h"
/* code optimized for HiFi5 core */

#if USE_NN_EXTENSION_8X16==0
void mtx_mpyt8x16_copyz(  
                     int16_t* restrict z,
               const int16_t* restrict v,
               int M, int P )
#if 0
{
    int m,p,P0;
    P0=(P+3)&~3;
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
    ae_int16 * restrict pZ;
    int m,p,P0;
    P0=(P+3)&~3;
    NASSERT(P0>0);
    pZ=(ae_int16*)v;
    aZ=AE_ZALIGN64();
    for (m=0; m<M; m++) 
    {
        ae_int16x4 x;
        pZ=(ae_int16*)(v+m*P0);
        aV=AE_LA64_PP(pZ);
        for (p=0; p<(P>>2); p++) 
        {
            AE_LA16X4_IP(x,aV,castxcc(ae_int16x4,pZ));
            AE_SA16X4_IP(x,aZ,castxcc(ae_int16x4,z ));
        }
        AE_SA64POS_FP(aZ,z);
        for (p=0; p<(P&3); p++) 
        {
            AE_L16_IP(x,castxcc(ae_int16,pZ),2);
            AE_S16_0_IP(x,castxcc(ae_int16,z ),2);
        }
    }
}
#endif

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
void mtx_mpyt8x16_partiallyaligned (  
                     int16_t* restrict z,
               const int8_t * restrict x,
               const int16_t* restrict y,
               int M, int N, int P, int lsh )
{
    const ae_int8x8* restrict pX0;
    const ae_int8x8* restrict pX1;
    const ae_int8x8* restrict pX2;
    const ae_int8x8* restrict pX3;
    const ae_int16x4* restrict pY0;
    const ae_int16x4* restrict pY2;

    ae_int16x4* restrict pZ ;

    int m,n,p,P0,N0;
    int sa=lsh+41;

    /* circular buffer definition is needed to avoid reading 
       after the end of buffer during the last iterations 
       in that cases we are reading the data using circular 
       addressing instead of linear one
     */
//    WAE_CBEGIN0((uintptr_t)(x)&~(HIFI_SIMD_WIDTH-1) );
//    WAE_CEND0  ((uintptr_t)(x+M*N+(HIFI_SIMD_WIDTH-1))&~(HIFI_SIMD_WIDTH-1));
    WAE_CBEGIN0((uintptr_t)(x)&~7 );
    WAE_CEND0  ((uintptr_t)(x+M*N+7)&~7);

    N0=(N+7)&~7;
    P0=(P+3)&~3;
    pZ =(ae_int16x4*)z;
    for (m=0; m<(M&~3);m+=4)
    {
        pX0=(const ae_int8x8*)&x[(m+0)*N];
        pX1=(const ae_int8x8*)        (N+(uintptr_t)pX0);
        pX2=(const ae_int8x8*)XT_ADDX2(N,(uintptr_t)pX0);
        pX3=(const ae_int8x8*)XT_ADDX2(N,(uintptr_t)pX1);
        pY0=(const ae_int16x4*)y;
        pY2=(const ae_int16x4*)XT_ADDX4(N0,(uintptr_t)y);
        for (p=0; p<(P0>>2); p++)
        {
            ae_int64 a00,a01,a02,a03,
                     a10,a11,a12,a13,
                     a20,a21,a22,a23,
                     a30,a31,a32,a33;
            ae_valign aX0,aX1,aX2,aX3;
            ae_int16x4 r0,r1,r2,r3;
            aX0=AE_LA64_PP(pX0);
            aX1=AE_LA64_PP(pX1);
            aX2=AE_LA64_PP(pX2);
            aX3=AE_LA64_PP(pX3);
            NASSERT_ALIGN(pY0,16);
            NASSERT_ALIGN(pY2,16);
            AE_MOVDX2(a00,a01,0,0);AE_MOVDX2(a02,a03,0,0);
            AE_MOVDX2(a10,a11,0,0);AE_MOVDX2(a12,a13,0,0);
            AE_MOVDX2(a20,a21,0,0);AE_MOVDX2(a22,a23,0,0);
            AE_MOVDX2(a30,a31,0,0);AE_MOVDX2(a32,a33,0,0);
            __Pragma("loop_count min=1")
            for (n=0; n<(N0>>3); n++)
            {
                ae_int8x8 x0,x1,x2,x3;
                ae_int16x4 y00,y01,y10,y11,y20,y21,y30,y31;
                AE_L16X4X2_X (y10,y11,(const ae_int16x8*)pY0 ,N0*sizeof(int16_t));
                AE_L16X4X2_IP(y00,y01,castxcc(ae_int16x8,pY0),sizeof(ae_int16x8));
                AE_L16X4X2_X (y30,y31,(const ae_int16x8*)pY2 ,N0*sizeof(int16_t));
                AE_L16X4X2_IP(y20,y21,castxcc(ae_int16x8,pY2),sizeof(ae_int16x8));
                AE_LA8X8_IP(x0,aX0,pX0);
                AE_LA8X8_IP(x1,aX1,pX1);
                AE_LA8X8_IP(x2,aX2,pX2);
                AE_LA8X8_IP(x3,aX3,pX3);
                AE_MULAAAA2Q16X8(a00,a10,y00,y11,x0);
                AE_MULAAAA2Q16X8(a10,a00,y10,y01,x0);
                AE_MULAAAA2Q16X8(a20,a30,y20,y31,x0);
                AE_MULAAAA2Q16X8(a30,a20,y30,y21,x0);

                AE_MULAAAA2Q16X8(a01,a11,y00,y11,x1);
                AE_MULAAAA2Q16X8(a11,a01,y10,y01,x1);
                AE_MULAAAA2Q16X8(a21,a31,y20,y31,x1);
                AE_MULAAAA2Q16X8(a31,a21,y30,y21,x1);

                AE_MULAAAA2Q16X8(a02,a12,y00,y11,x2);
                AE_MULAAAA2Q16X8(a12,a02,y10,y01,x2);
                AE_MULAAAA2Q16X8(a22,a32,y20,y31,x2);
                AE_MULAAAA2Q16X8(a32,a22,y30,y21,x2);

                AE_MULAAAA2Q16X8(a03,a13,y00,y11,x3);
                AE_MULAAAA2Q16X8(a13,a03,y10,y01,x3);
                AE_MULAAAA2Q16X8(a23,a33,y20,y31,x3);
                AE_MULAAAA2Q16X8(a33,a23,y30,y21,x3);
            }
            r0=AE_ROUND16X4F32SASYM(AE_TRUNCA32X2F64S(a00,a10, sa), AE_TRUNCA32X2F64S(a20,a30, sa));
            r1=AE_ROUND16X4F32SASYM(AE_TRUNCA32X2F64S(a01,a11, sa), AE_TRUNCA32X2F64S(a21,a31, sa));
            r2=AE_ROUND16X4F32SASYM(AE_TRUNCA32X2F64S(a02,a12, sa), AE_TRUNCA32X2F64S(a22,a32, sa));
            r3=AE_ROUND16X4F32SASYM(AE_TRUNCA32X2F64S(a03,a13, sa), AE_TRUNCA32X2F64S(a23,a33, sa));
            AE_S16X4_X (r1,pZ,1*P0*sizeof(int16_t));
            AE_S16X4_X (r2,pZ,2*P0*sizeof(int16_t));
            AE_S16X4_X (r3,pZ,3*P0*sizeof(int16_t));
            AE_S16X4_IP(r0,pZ,sizeof(ae_int16x4));
            pX0=(const ae_int8x8*)((uintptr_t)pX0-N0);
            pX1=(const ae_int8x8*)((uintptr_t)pX1-N0);
            pX2=(const ae_int8x8*)((uintptr_t)pX2-N0);
            pX3=(const ae_int8x8*)((uintptr_t)pX3-N0);
            pY0=(const ae_int16x4*)((uintptr_t)pY0+3*N0*sizeof(int16_t));
            pY2=(const ae_int16x4*)((uintptr_t)pY2+3*N0*sizeof(int16_t));
        }
        pZ=(ae_int16x4*)(3*P0*sizeof(int16_t)+(uintptr_t)pZ);
    }
    if (M&2)
    {
        pX0=(const ae_int8x8*)&x[(m+0)*N];
        pX1=(const ae_int8x8*)        (N+(uintptr_t)pX0);
        pY0=(const ae_int16x4*)y;
        pY2=(const ae_int16x4*)XT_ADDX4(N0,(uintptr_t)y);
        for (p=0; p<(P0>>2); p++)
        {
            ae_int64 a00,a01,
                     a10,a11,
                     a20,a21,
                     a30,a31;
            ae_valign aX0,aX1;
            ae_int16x4 r0,r1;
            pX1=(const ae_int8x8*)        (N+(uintptr_t)pX0);
            aX0=AE_LA64_PP(pX0);
            //aX1=AE_LA64_PP(pX1);
            AE_LA8X8POS_PC(aX1,pX1);
            NASSERT_ALIGN(pY0,16);
            NASSERT_ALIGN(pY2,16);
            AE_MOVDX2(a00,a01,0,0);
            AE_MOVDX2(a10,a11,0,0);
            AE_MOVDX2(a20,a21,0,0);
            AE_MOVDX2(a30,a31,0,0);
            __Pragma("loop_count min=1")
            for (n=0; n<(N0>>3); n++)
            {
                ae_int8x8 x0,x1;
                ae_int16x4 y00,y01,y10,y11,y20,y21,y30,y31;
                AE_L16X4X2_X (y10,y11,(const ae_int16x8*)pY0 ,N0*sizeof(int16_t));
                AE_L16X4X2_IP(y00,y01,castxcc(ae_int16x8,pY0),sizeof(ae_int16x8));
                AE_L16X4X2_X (y30,y31,(const ae_int16x8*)pY2 ,N0*sizeof(int16_t));
                AE_L16X4X2_IP(y20,y21,castxcc(ae_int16x8,pY2),sizeof(ae_int16x8));
                AE_LA8X8_IP(x0,aX0,pX0);
                //AE_LA8X8_IP(x1,aX1,pX1);
                AE_LA8X8_IC(x1,aX1,pX1);
                AE_MULAAAA2Q16X8(a00,a10,y00,y11,x0);
                AE_MULAAAA2Q16X8(a10,a00,y10,y01,x0);
                AE_MULAAAA2Q16X8(a20,a30,y20,y31,x0);
                AE_MULAAAA2Q16X8(a30,a20,y30,y21,x0);

                AE_MULAAAA2Q16X8(a01,a11,y00,y11,x1);
                AE_MULAAAA2Q16X8(a11,a01,y10,y01,x1);
                AE_MULAAAA2Q16X8(a21,a31,y20,y31,x1);
                AE_MULAAAA2Q16X8(a31,a21,y30,y21,x1);

            }
            r0=AE_ROUND16X4F32SASYM(AE_TRUNCA32X2F64S(a00,a10, sa), AE_TRUNCA32X2F64S(a20,a30, sa));
            r1=AE_ROUND16X4F32SASYM(AE_TRUNCA32X2F64S(a01,a11, sa), AE_TRUNCA32X2F64S(a21,a31, sa));
            AE_S16X4_X (r1,pZ,1*P0*sizeof(int16_t));
            AE_S16X4_IP(r0,pZ,sizeof(ae_int16x4));
            pX0=(const ae_int8x8*)((uintptr_t)pX0-N0);
            pY0=(const ae_int16x4*)((uintptr_t)pY0+3*N0*sizeof(int16_t));
            pY2=(const ae_int16x4*)((uintptr_t)pY2+3*N0*sizeof(int16_t));
        }
        pZ=(ae_int16x4*)(1*P0*sizeof(int16_t)+(uintptr_t)pZ);
        m+=2;
    }
    if (M&1)
    {
        pY0=(const ae_int16x4*)y;
        pY2=(const ae_int16x4*)XT_ADDX4(N0,(uintptr_t)y);
        for (p=0; p<(P0>>2); p++)
        {
            ae_int64 a00,a10,a20,a30;
            ae_valign aX0;
            ae_int16x4 r0;
            pX0=(const ae_int8x8*)&x[(m+0)*N];
            AE_LA8X8POS_PC(aX0,pX0);
            NASSERT_ALIGN(pY0,16);
            NASSERT_ALIGN(pY2,16);
            AE_MOVDX2(a00,a10,0,0);
            AE_MOVDX2(a20,a30,0,0);
            __Pragma("loop_count min=1")
            for (n=0; n<(N0>>3); n++)
            {
                ae_int8x8 x0;
                ae_int16x4 y00,y01,y10,y11,y20,y21,y30,y31;
                AE_L16X4X2_X (y10,y11,(const ae_int16x8*)pY0 ,N0*sizeof(int16_t));
                AE_L16X4X2_IP(y00,y01,castxcc(ae_int16x8,pY0),sizeof(ae_int16x8));
                AE_L16X4X2_X (y30,y31,(const ae_int16x8*)pY2 ,N0*sizeof(int16_t));
                AE_L16X4X2_IP(y20,y21,castxcc(ae_int16x8,pY2),sizeof(ae_int16x8));
                AE_LA8X8_IC(x0,aX0,pX0);
                AE_MULAAAA2Q16X8(a00,a10,y00,y11,x0);
                AE_MULAAAA2Q16X8(a10,a00,y10,y01,x0);
                AE_MULAAAA2Q16X8(a20,a30,y20,y31,x0);
                AE_MULAAAA2Q16X8(a30,a20,y30,y21,x0);
            }
            r0=AE_ROUND16X4F32SASYM(AE_TRUNCA32X2F64S(a00,a10, sa), AE_TRUNCA32X2F64S(a20,a30, sa));
            AE_S16X4_IP(r0,pZ,sizeof(ae_int16x4));
            pY0=(const ae_int16x4*)((uintptr_t)pY0+3*N0*sizeof(int16_t));
            pY2=(const ae_int16x4*)((uintptr_t)pY2+3*N0*sizeof(int16_t));
        }
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
void mtx_mpy8x16_copysmallx(int8_t * restrict r,const int8_t * restrict x,int M,int N)
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
#endif /* USE_NN_EXTENSION_8X16 */
