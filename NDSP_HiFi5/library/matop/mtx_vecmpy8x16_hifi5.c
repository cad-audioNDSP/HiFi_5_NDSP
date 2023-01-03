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
#include "mtx_mpy8x16_helper.h"
/* Code optimized for HiFi4 core */

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
void mtx_vecmpy8x16   ( int16_t* restrict z,
                 const int8_t* restrict x,
                 const int16_t* restrict y,
                 int M, int N, int lsh)
#if 0
{
    int m,n;
    NASSERT(lsh >= -15 && lsh <= 15);
    if (M<=0) return;
    if (N<=0)    /* exceptional situation */
    {
        for (m=0; m<M; m++) z[m]=0;
        return ;
    }
    ae_int64 rnd;
    int rsh=7-lsh;
	rnd=AE_SLAA64S(1,rsh+15);

    for (m=0; m<(M&~1); m+=2)
    {
		const ae_int16x4 * pY;
		ae_valign yv;
		pY=(const ae_int16x4 *)(y);
		yv=AE_LA64_PP(pY);
        ae_int64 a0,a1;
        a0=a1=AE_ZERO64();
        for (n=0; n<(N&~3); n+=4)
        {
            ae_int32x2 t0,t1,t2,t3;
			ae_int16x4 y0;
			ae_int32x2 q0,q1;
			ae_int32x2 s0,s1;
			AE_LA16X4_IP(y0,yv,pY);
            t0=AE_MOVDA32X2((uint8_t)x[m*N+n+0*N  ],(uint8_t)x[m*N+n+0*N+1]);
            t1=AE_MOVDA32X2((uint8_t)x[m*N+n+1*N  ],(uint8_t)x[m*N+n+1*N+1]);
            t2=AE_MOVDA32X2((uint8_t)x[m*N+n+0*N+2],(uint8_t)x[m*N+n+0*N+3]);
            t3=AE_MOVDA32X2((uint8_t)x[m*N+n+1*N+2],(uint8_t)x[m*N+n+1*N+3]);
            t0=AE_SEXT32(t0,7);
            t1=AE_SEXT32(t1,7);
            t2=AE_SEXT32(t2,7);
            t3=AE_SEXT32(t3,7);

			q0=AE_SEXT32(          AE_MOVINT32X2_FROMINT16X4(y0)    ,15);
			q1=AE_SEXT32(AE_SRLI32(AE_MOVINT32X2_FROMINT16X4(y0),16),15);
			s0=AE_SEL32_HH(q1,q0);
			s1=AE_SEL32_LL(q1,q0);

            AE_MULAAD32_HH_LL(a0,t0,s0);
            AE_MULAAD32_HH_LL(a0,t2,s1);
            AE_MULAAD32_HH_LL(a1,t1,s0);
            AE_MULAAD32_HH_LL(a1,t3,s1);
        }
		if (N&2)
		{
            ae_int32x2 s0,t0,t1;
            s0=AE_MOVDA32X2((uint16_t)y[n  ],(uint16_t)y[n+1]);
            t0=AE_MOVDA32X2((uint8_t)x[m*N+n+0*N],(uint8_t)x[m*N+n+0*N+1]);
            t1=AE_MOVDA32X2((uint8_t)x[m*N+n+1*N],(uint8_t)x[m*N+n+1*N+1]);
            t0=AE_SEXT32(t0,7);
            t1=AE_SEXT32(t1,7);
            s0=AE_SEXT32(s0,15);
            AE_MULAAD32_HH_LL(a0,t0,s0);
            AE_MULAAD32_HH_LL(a1,t1,s0);
			n+=2;
		}
        if (N&1)
        {
            ae_int32x2 s0,t0;
            s0=AE_MOVDA32X2((uint16_t)y[n  ],(uint16_t)y[n  ]);
            t0=AE_MOVDA32X2((uint8_t)x[m*N+n+0*N],(uint8_t)x[m*N+n+1*N]);
            t0=AE_SEXT32(t0,7);
            s0=AE_SEXT32(s0,15);
            AE_MULA32_HH(a0,t0,s0);
            AE_MULA32_LL(a1,t0,s0);
        }
		// Pack, round, store
        {
            ae_f16x4 t0;
            t0=AE_ROUND16X4F32SASYM(AE_TRUNCA32X2F64S(a0,a1,48-rsh),AE_TRUNCA32X2F64S(AE_ZERO64(),AE_ZERO64(),48-rsh));
			z[m]=AE_MOVAD16_3(t0);
			z[m+1]=AE_MOVAD16_2(t0);
        }
    }
	if (M&1)
	{
        ae_int64 a0;
        a0=AE_ZERO64();
        for (n=0; n<(N&~1); n+=2)
        {
            ae_int32x2 s0,t0;
            s0=AE_MOVDA32X2((uint16_t)y[n+0],(uint16_t)y[n+1]);
            t0=AE_MOVDA32X2((uint8_t)x[m*N+n+0*N],(uint8_t)x[m*N+n+0*N+1]);
            t0=AE_SEXT32(t0,7);
            s0=AE_SEXT32(s0,15);
            AE_MULAAD32_HH_LL(a0,t0,s0);
        }
        if (N&1)
        {
            ae_int32x2 s0,t0;
            s0=AE_MOVDA32X2((uint16_t)y[n],(uint16_t)y[n]);
            t0=AE_MOVDA32X2((uint8_t)x[m*N+n+0*N],(uint8_t)x[m*N+n+0*N]);
            t0=AE_SEXT32(t0,7);
            s0=AE_SEXT32(s0,15);
            AE_MULA32_HH(a0,t0,s0);
        }
		// Pack, round, store
        {
            ae_f16x4 t0;
            t0=AE_ROUND16X4F32SASYM(AE_TRUNCA32X2F64S(a0,a0,48-rsh),AE_TRUNCA32X2F64S(a0,a0,48-rsh));
            z[m]=(int16_t)AE_MOVAD16_3(t0);
        }
	}
}
#elif USE_NN_EXTENSION_8X16==0
{
    ae_int16 * restrict pZ=(ae_int16 * )z;
    static const uint64_t dsel_rephhll_tbl=0x40516273c8d9eafbULL;
    const ae_int8x8 * restrict pX0;
    const ae_int8x8 * restrict pX1;
    const ae_int8x8 * restrict pX2;
    const ae_int8x8 * restrict pX3;
    const ae_int16x4 * restrict pY;
    ae_valign aX0,aX1,aX2,aX3;
    ae_valignx2 aY;
    int m,n;
    int sa=41+lsh;
    ae_int8x8 dsel_rephhll;
    NASSERT(lsh >= -15 && lsh <= 15);
    if (M<=0) return;
    if (N<=0)    /* exceptional situation */
    {
        for (m=0; m<M; m++) AE_S16_0_IP(0,pZ,sizeof(int16_t));
        return ;
    }
    dsel_rephhll=AE_L8X8_I((const ae_int8x8*)&dsel_rephhll_tbl,0);
    pX0=(const ae_int8x8*)( x );
    pX1=(const ae_int8x8*)( (uintptr_t)x + N);
    pX2=(const ae_int8x8*)XT_ADDX2(N, (uintptr_t)pX0);
    pX3=(const ae_int8x8*)XT_ADDX2(N, (uintptr_t)pX1);
    /* separate code for N<8 to avoid reads beyond the array boundary */
    if (N<8)
    {
        for (m=0; m<(M&~3);m+=4)
        {
            ae_int64 A0,A1,A2,A3;
            AE_MOVDX2(A0,A1,0,0); AE_MOVDX2(A2,A3,0,0);
            pY =(const ae_int16x4*)y;
            __Pragma("loop_count min=1,max=7")
            for(n=0; n<N; n++)
            {
                ae_int16x4 y0;
                ae_int8x8 x0,x1,x2,x3,x00,x11,x22,x33;
                AE_L16_IP(y0,castxcc(ae_int16,pY),sizeof(int16_t));
                AE_L8_IP(x0,castxcc(ae_int8,pX0),1);
                AE_L8_IP(x1,castxcc(ae_int8,pX1),1);
                AE_L8_IP(x2,castxcc(ae_int8,pX2),1);
                AE_L8_IP(x3,castxcc(ae_int8,pX3),1);
                AE_DSEL8X8(x00,x11,x0,x1,dsel_rephhll);
                AE_DSEL8X8(x22,x33,x2,x3,dsel_rephhll);
                AE_MULAAAA2Q16X8(A0, A1, y0, y0, x00);
                AE_MULAAAA2Q16X8(A2, A3, y0, y0, x22);
            }
            AE_S16_0_IP(AE_ROUND16X4F32SASYM(AE_TRUNCA32X2F64S(A0,A0,sa-2),AE_TRUNCA32X2F64S(A0,A0,sa-2)),pZ,sizeof(int16_t));
            AE_S16_0_IP(AE_ROUND16X4F32SASYM(AE_TRUNCA32X2F64S(A1,A1,sa-2),AE_TRUNCA32X2F64S(A1,A1,sa-2)),pZ,sizeof(int16_t));
            AE_S16_0_IP(AE_ROUND16X4F32SASYM(AE_TRUNCA32X2F64S(A2,A2,sa-2),AE_TRUNCA32X2F64S(A2,A2,sa-2)),pZ,sizeof(int16_t));
            AE_S16_0_IP(AE_ROUND16X4F32SASYM(AE_TRUNCA32X2F64S(A3,A3,sa-2),AE_TRUNCA32X2F64S(A3,A3,sa-2)),pZ,sizeof(int16_t));
            pX0=(const ae_int8x8*)( (uintptr_t)pX0 + 3*N);
            pX1=(const ae_int8x8*)( (uintptr_t)pX1 + 3*N);
            pX2=(const ae_int8x8*)( (uintptr_t)pX2 + 3*N);
            pX3=(const ae_int8x8*)( (uintptr_t)pX3 + 3*N);
        }
        for (;m<M; m++)
        {
            ae_int64 A0,A1;
            AE_MOVDX2(A0,A1,0,0);
            pY =(const ae_int16x4*)y;
            __Pragma("loop_count min=1,max=7")
            for(n=0; n<N; n++)
            {
                ae_int16x4 y0;
                ae_int8x8 x0;
                AE_L16_IP(y0,castxcc(ae_int16,pY),sizeof(int16_t));
                AE_L8_IP(x0,castxcc(ae_int8,pX0),1);
                AE_MULAAAA2Q16X8(A0, A1, y0, 0, x0);
            }
            A0=A0+A1;
            AE_S16_0_IP(AE_ROUND16X4F32SASYM(AE_TRUNCA32X2F64S(A0,A0,sa-2),AE_TRUNCA32X2F64S(A0,A0,sa-2)),pZ,sizeof(int16_t));
        }
        return;
    }
    for (m=0; m<(M&~3);m+=4)
    {
        ae_int64 A0,A1,A2,A3;
        AE_MOVDX2(A0,A1,0,0); AE_MOVDX2(A2,A3,0,0);
        pY =(const ae_int16x4*)y;
        aX0=AE_LA64_PP(pX0);
        aX1=AE_LA64_PP(pX1);
        aY =AE_LA128_PP(pY);
        __Pragma("loop_count min=1")
        for (n=0; n<(N>>3); n++)
        {
            ae_int16x4 y0,y1;
            ae_int8x8 x0,x1,x2,x3,x00,x11,x22,x33;
            AE_LA16X4X2_IP(y0,y1,aY,castxcc(ae_int16x8,pY));
                                 AE_LA8X8_IP(x0,aX0,pX0);
                                 AE_LA8X8_IP(x1,aX1,pX1);
            aX2=AE_LA64_PP(pX2); AE_LA8X8_IP(x2,aX2,pX2);
            aX3=AE_LA64_PP(pX3); AE_LA8X8_IP(x3,aX3,pX3);
            AE_DSEL8X8(x00,x11,x0,x1,dsel_rephhll);
            AE_DSEL8X8(x22,x33,x2,x3,dsel_rephhll);
            AE_MULAAAA2Q16X8(A0, A1, y0, y0, x00);
            AE_MULAAAA2Q16X8(A0, A1, y1, y1, x11);
            AE_MULAAAA2Q16X8(A2, A3, y0, y0, x22);
            AE_MULAAAA2Q16X8(A2, A3, y1, y1, x33);
        }
        for(n=0; n<(N&7); n++)
        {
            ae_int16x4 y0;
            ae_int8x8 x0,x1,x2,x3,x00,x11,x22,x33;
            AE_L16_IP(y0,castxcc(ae_int16,pY),sizeof(int16_t));
            AE_L8_IP(x0,castxcc(ae_int8,pX0),1);
            AE_L8_IP(x1,castxcc(ae_int8,pX1),1);
            AE_L8_IP(x2,castxcc(ae_int8,pX2),1);
            AE_L8_IP(x3,castxcc(ae_int8,pX3),1);
            AE_MOVF16X4(y0,0,AE_MOVBA4(1));
            AE_DSEL8X8(x00,x11,x0,x1,dsel_rephhll);
            AE_DSEL8X8(x22,x33,x2,x3,dsel_rephhll);
            AE_MULAAAA2Q16X8(A0, A1, y0, y0, x00);
            AE_MULAAAA2Q16X8(A2, A3, y0, y0, x22);
        }
        AE_S16_0_IP(AE_ROUND16X4F32SASYM(AE_TRUNCA32X2F64S(A0,A0,sa),AE_TRUNCA32X2F64S(A0,A0,sa)),pZ,sizeof(int16_t));
        AE_S16_0_IP(AE_ROUND16X4F32SASYM(AE_TRUNCA32X2F64S(A1,A1,sa),AE_TRUNCA32X2F64S(A1,A1,sa)),pZ,sizeof(int16_t));
        AE_S16_0_IP(AE_ROUND16X4F32SASYM(AE_TRUNCA32X2F64S(A2,A2,sa),AE_TRUNCA32X2F64S(A2,A2,sa)),pZ,sizeof(int16_t));
        AE_S16_0_IP(AE_ROUND16X4F32SASYM(AE_TRUNCA32X2F64S(A3,A3,sa),AE_TRUNCA32X2F64S(A3,A3,sa)),pZ,sizeof(int16_t));
        pX0=(const ae_int8x8*)( (uintptr_t)pX0 + 3*N);
        pX1=(const ae_int8x8*)( (uintptr_t)pX1 + 3*N);
        pX2=(const ae_int8x8*)( (uintptr_t)pX2 + 3*N);
        pX3=(const ae_int8x8*)( (uintptr_t)pX3 + 3*N);
    }
    if (M&2)
    {
        ae_int64 A0,A1;
        AE_MOVDX2(A0,A1,0,0);
        pY =(const ae_int16x4*)y;
        aX0=AE_LA64_PP(pX0);
        aX1=AE_LA64_PP(pX1);
        aY =AE_LA128_PP(pY);
        __Pragma("loop_count min=1")
        for (n=0; n<(N>>3); n++)
        {
            ae_int16x4 y0,y1;
            ae_int8x8 x0,x1,x00,x11;
            AE_LA16X4X2_IP(y0,y1,aY,castxcc(ae_int16x8,pY));
            AE_LA8X8_IP(x0,aX0,pX0);
            AE_LA8X8_IP(x1,aX1,pX1);
            AE_DSEL8X8(x00,x11,x0,x1,dsel_rephhll);
            AE_MULAAAA2Q16X8(A0, A1, y0, y0, x00);
            AE_MULAAAA2Q16X8(A0, A1, y1, y1, x11);
        }
        for(n=0; n<(N&7); n++)
        {
            ae_int16x4 y0;
            ae_int8x8 x0,x1,x00,x11;
            AE_L16_IP(y0,castxcc(ae_int16,pY),sizeof(int16_t));
            AE_L8_IP(x0,castxcc(ae_int8,pX0),1);
            AE_L8_IP(x1,castxcc(ae_int8,pX1),1);
            AE_MOVF16X4(y0,0,AE_MOVBA4(1));
            AE_DSEL8X8(x00,x11,x0,x1,dsel_rephhll);
            AE_MULAAAA2Q16X8(A0, A1, y0, y0, x00);
        }
        AE_S16_0_IP(AE_ROUND16X4F32SASYM(AE_TRUNCA32X2F64S(A0,A0,sa),AE_TRUNCA32X2F64S(A0,A0,sa)),pZ,sizeof(int16_t));
        AE_S16_0_IP(AE_ROUND16X4F32SASYM(AE_TRUNCA32X2F64S(A1,A1,sa),AE_TRUNCA32X2F64S(A1,A1,sa)),pZ,sizeof(int16_t));
        pX0=(const ae_int8x8*)( (uintptr_t)pX0 + 1*N);
        m+=2;
    }
    if (M&1)
    {
        ae_int64 A0,A1;
        AE_MOVDX2(A0,A1,0,0);
        pY =(const ae_int16x4*)y;
        aX0=AE_LA64_PP(pX0);
        aY =AE_LA128_PP(pY);
        __Pragma("loop_count min=1")
        for (n=0; n<(N>>3); n++)
        {
            ae_int16x4 y0,y1;
            ae_int8x8 x0;
            AE_LA16X4X2_IP(y0,y1,aY,castxcc(ae_int16x8,pY));
            AE_LA8X8_IP(x0,aX0,pX0);
            AE_MULAAAA2Q16X8(A0, A1, y0, y1, x0);
        }
        for(n=0; n<(N&7); n++)
        {
            ae_int16x4 y0;
            ae_int8x8 x0;
            AE_L16_IP(y0,castxcc(ae_int16,pY),sizeof(int16_t));
            AE_L8_IP(x0,castxcc(ae_int8,pX0),1);
            AE_MOVF16X4(y0,0,AE_MOVBA4(1));
            AE_MULAAAA2Q16X8(A0, A1, y0, 0, x0);
        }
        A0=A0+A1;
        AE_S16_0_IP(AE_ROUND16X4F32SASYM(AE_TRUNCA32X2F64S(A0,A0,sa),AE_TRUNCA32X2F64S(A0,A0,sa)),pZ,sizeof(int16_t));
    }
}
#else
{
    int m,n;
    ae_valignx2 alX0, alX1, alX2, alX3;
    ae_valignx2 alY0;
    ae_valignx2 alZ0;
    int16_t * restrict pZ0;
    const int8_t * restrict pX0;
    const int8_t * restrict pX1;
    const int8_t * restrict pX2;
    const int8_t * restrict pX3;
    const int16_t * restrict pY0;

    ae_int32x2 Z00, Z01;
    ae_int8x8 x00, x01, x10, x11, x20, x21, x30, x31;
    ae_int16x4 y00, y01, y02, y03;
    ae_int16x4 z0;
    ae_int16x4 zero;

    NASSERT(lsh >= -15 && lsh <= 15);

    if (M<=0) return;
    zero = AE_ZERO16();
    if (N<=0)    /* exceptional situation */
    {
        pZ0 = z;
        alZ0 = AE_ZALIGN128();
        for (m=0; m<(M>>3); m++)  AE_SA16X4X2_IP(zero, zero, alZ0, castxcc(ae_int16x8,pZ0));
        AE_SAV16X4X2_XP(zero, zero, alZ0, castxcc(ae_int16x8,pZ0), (M&7)*sizeof(int16_t));
        AE_SA128POS_FP(alZ0, pZ0);
        return;
    }

    WUR_AE_SAR(9+lsh);
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
            alY0 = AE_LA128_PP(pY0);
            AE_LA16X4X2_IP(y00, y01, alY0, castxcc(ae_int16x8,pY0));
            AE_LA16X4X2_IP(y02, y03, alY0, castxcc(ae_int16x8,pY0));

            /* make multiply, using 4-way 32-bit output octal mac, 8X8 bit */
            AE_MULA8Q8X16(Z00, Z01, x00, x10, x20, x30, y00, y01);
            AE_MULA8Q8X16(Z00, Z01, x01, x11, x21, x31, y02, y03);
        }
        /* Process last 1...16 multiplications of each output value */
        {
            int tail0, tail1;
            ae_int8x8 tmp;
            /* read last values in a different manner to prevent out-of-bound memory accesses */
            tail0 = ((N-1)&15)+1 - 8;
            tail1 = XT_MIN(((N-1)&15)+1, 8);
            /* load x matrix, 4x16 values, 8-bit */
            alX0 = AE_LA128_PP(pX0);
            alX1 = AE_LA128_PP(pX1);
            AE_LAV8X8X2_XP(x00, tmp, alX0, castxcc(ae_int8x16,pX0), tail0);
            AE_LAV8X8X2_XP(x10, tmp, alX1, castxcc(ae_int8x16,pX1), tail0);
            AE_LAV8X8X2_XP(x20, tmp, alX2, castxcc(ae_int8x16,pX2), tail0);
            AE_LAV8X8X2_XP(x30, tmp, alX3, castxcc(ae_int8x16,pX3), tail0);
            AE_LAV8X8X2_XP(x01, tmp, alX0, castxcc(ae_int8x16,pX0), tail1);
            AE_LAV8X8X2_XP(x11, tmp, alX1, castxcc(ae_int8x16,pX1), tail1);
            AE_LAV8X8X2_XP(x21, tmp, alX2, castxcc(ae_int8x16,pX2), tail1);
            AE_LAV8X8X2_XP(x31, tmp, alX3, castxcc(ae_int8x16,pX3), tail1);

            /* load y vector, 16 values, 8-bit */
            alY0 = AE_LA128_PP(pY0);
            AE_LAV16X4X2_XP(y00, y01, alY0, castxcc(ae_int16x8,pY0), tail0*sizeof(int16_t));
            AE_LAV16X4X2_XP(y02, y03, alY0, castxcc(ae_int16x8,pY0), tail1*sizeof(int16_t));

            /* make multiply, using 4-way 32-bit output octal mac, 8X8 bit */
            AE_MULA8Q8X16(Z00, Z01, x00, x10, x20, x30, y00, y01);
            AE_MULA8Q8X16(Z00, Z01, x01, x11, x21, x31, y02, y03);
        }
        /* Q15 + lsh <- Q22 - 7 + lsh w/ rounding and saturation */
        Z00 = AE_SLAS32S(Z00);
        Z01 = AE_SLAS32S(Z01);
        z0 = AE_ROUND16X4F32SASYM(Z00, Z01);
        alZ0 = AE_ZALIGN128();
        AE_SAV16X4X2_XP(z0, z0, alZ0, castxcc(ae_int16x8,pZ0), XT_MIN(4, M-m)*sizeof(int16_t));
        AE_SA128POS_FP(alZ0, pZ0);
    }
}
#endif
