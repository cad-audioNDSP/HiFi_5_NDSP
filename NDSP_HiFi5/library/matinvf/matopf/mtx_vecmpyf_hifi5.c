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
  NatureDSP Signal Processing Library. Matrix Operations
  Matrix multiply
  * Optimized code for HiFi4
  IntegrIT, 2006-2019
  */
/* Code optimized for HiFi5 core */

/* Cross-platform data type definitions. */
#include "NatureDSP_types.h"
#include "NatureDSP_Signal_matop.h"
/* Common helper macros. */
#include "common.h"
#include "common_fpu.h"

#if (HAVE_FPU==0 && HAVE_VFPU==0)
DISCARD_FUN(void ,mtx_vecmpyf,( float32_t * z, const float32_t * x,  const float32_t * y, int M, int N ))
#elif (HAVE_VFPU)
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
void mtx_vecmpyf(float32_t * z, const float32_t * x, const float32_t * y, int M, int N)
{
    const xtfloatx4 *restrict px0;
    const xtfloatx4 *restrict px1;
    const xtfloatx4 *restrict px2;
    const xtfloatx4 *restrict px3;
    const xtfloatx4 *restrict py;
    xtfloatx2 *restrict pz_2;
    xtfloatx4 *restrict pz_4;
    ae_int32x2 zero_v;
    xtfloatx2 y0, y1;
    xtfloatx2 z0, z1;
    xtfloatx2 x00, x01, x10, x11;
    xtfloatx2 x20, x21, x30, x31;

    xtfloatx2 acc00, acc01, acc10, acc11;
    xtfloatx2 acc20, acc21, acc30, acc31;
    xtfloatx2 acc02, acc03, acc12, acc13;
    xtfloatx2 acc22, acc23, acc32, acc33;
    xtfloat z0_;
    ae_valignx2 vx0, vx1, vy;
    ae_valignx2 vx2, vx3;
    ae_valign vz_2;
    ae_valignx2 vz_4;
    xtbool2 bmask0, bmask1;
    int m, n;

    NASSERT(x);
    NASSERT(y);
    NASSERT(z);
    NASSERT((z != x) && (z != y));

    if (M <= 0)    return;
    /* If N<=0 then clear output vector and return */
    if (N <= 0)
    {
        for (m = 0; m<M; m++) z[m] = 0.f;
        return;
    }
    bmask0 = AE_MOVBA2((N & 2) + ((N & 2) >> 1) + 2 * (int)((N & 3) == 1));  // 2*(N%4) if N%4<2, else 3 
    bmask1 = AE_MOVBA2(((int)((N & 3) == 3)) << 1);  // 2 if (N%4)=3, else 0

    x00 = (xtfloatx2)0.f;
    zero_v = AE_MOVINT32X2_FROMXTFLOATX2(x00);
    
    if (N&4 && M>=4 ){
        pz_4 = (xtfloatx4 *)z;

        __Pragma("loop_count min=1")
        for (m = 0; m < (M >> 2); m++)
        {
            px0 = (const xtfloatx4 *) (x);
            px1 = (const xtfloatx4 *) XT_ADDX4(N,(uintptr_t) px0);
            px2 = (const xtfloatx4 *) XT_ADDX4(N,(uintptr_t) px1);
            px3 = (const xtfloatx4 *) XT_ADDX4(N,(uintptr_t) px2);
            x   = (const float32_t *) XT_ADDX4(N,(uintptr_t) px3);
            py  = (const xtfloatx4 *)(y);

            //vx0 = AE_LA128_PP(px0);
            vx1 = AE_LA128_PP(px1);
            vx2 = AE_LA128_PP(px2);
            vx3 = AE_LA128_PP(px3);

            vy = AE_LA128_PP(py);
            acc00 = acc01 = acc10 = acc11 = 0.0f;
            acc20 = acc21 = acc30 = acc31 = 0.0f;
            acc02 = acc03 = acc12 = acc13 = 0.0f;
            acc22 = acc23 = acc32 = acc33 = 0.0f;

            for (n = 0; n < (N >> 3); n++)
            {
                vx0 = AE_LA128_PP(px0); AE_LASX2X2_IP(x00, x01, vx0, px0);
                /*vx1 = AE_LA128_PP(px1);*/ AE_LASX2X2_IP(x10, x11, vx1, px1);
                AE_LASX2X2_IP(x20, x21, vx2, px2);
                AE_LASX2X2_IP(x30, x31, vx3, px3);

                AE_LASX2X2_IP(y0, y1, vy, py);

                MADDQ_S(acc00, acc10, x00, x10, y0);
                MADDQ_S(acc01, acc11, x01, x11, y1);

                MADDQ_S(acc20, acc30, x20, x30, y0);
                MADDQ_S(acc21, acc31, x21, x31, y1);

                AE_LASX2X2_IP(x00, x01, vx0, px0);
                AE_LASX2X2_IP(x10, x11, vx1, px1);
                AE_LASX2X2_IP(x20, x21, vx2, px2);
                AE_LASX2X2_IP(x30, x31, vx3, px3);

                vy = AE_LA128_PP(py); AE_LASX2X2_IP(y0, y1, vy, py);

                MADDQ_S(acc02, acc12, x00, x10, y0);
                MADDQ_S(acc03, acc13, x01, x11, y1);

                MADDQ_S(acc22, acc32, x20, x30, y0);
                MADDQ_S(acc23, acc33, x21, x31, y1);
            }

            acc00 += acc02;  acc01 += acc03;
            acc10 += acc12;  acc11 += acc13;
            acc20 += acc22;  acc21 += acc23;
            acc30 += acc32;  acc31 += acc33;

            vx0 = AE_LA128_PP(px0); AE_LASX2X2_IP(x00, x01, vx0, px0);
            vx1 = AE_LA128_PP(px1); AE_LASX2X2_IP(x10, x11, vx1, px1);
            AE_LASX2X2_IP(x20, x21, vx2, px2);
            AE_LASX2X2_IP(x30, x31, vx3, px3);

            AE_LASX2X2_IP(y0, y1, vy, py);

            MADDQ_S(acc00, acc10, x00, x10, y0);
            MADDQ_S(acc01, acc11, x01, x11, y1);
            MADDQ_S(acc20, acc30, x20, x30, y0);
            MADDQ_S(acc21, acc31, x21, x31, y1);

            vx0 = AE_LA128_PP(px0); AE_LASX2X2_IP(x00, x01, vx0, px0);
            vx1 = AE_LA128_PP(px1); AE_LASX2X2_IP(x10, x11, vx1, px1);
            AE_LASX2X2_IP(x20, x21, vx2, px2);
            AE_LASX2X2_IP(x30, x31, vx3, px3);

            AE_LASX2X2_IP(y0, y1, vy, py);

            MOVF_SX2(y0 , zero_v, bmask0);
            MOVF_SX2(y1 , zero_v, bmask1);

            MOVF_SX2(x00, zero_v, bmask0);
            MOVF_SX2(x10, zero_v, bmask0);
            MOVF_SX2(x01, zero_v, bmask1);
            MOVF_SX2(x11, zero_v, bmask1);  

            MOVF_SX2(x20, zero_v, bmask0);
            MOVF_SX2(x30, zero_v, bmask0);
            MOVF_SX2(x21, zero_v, bmask1);
            MOVF_SX2(x31, zero_v, bmask1);

            MADDQ_S(acc00, acc10, x00, x10, y0);
            MADDQ_S(acc01, acc11, x01, x11, y1);
            MADDQ_S(acc20, acc30, x20, x30, y0);
            MADDQ_S(acc21, acc31, x21, x31, y1);
            
            acc00 = acc00 + acc01;
            acc10 = acc10 + acc11;
            acc20 = acc20 + acc21;
            acc30 = acc30 + acc31;

            z0 = XT_SEL32_HL_SX2(acc00, acc10) + XT_SEL32_LH_SX2(acc00, acc10);
            z1 = XT_SEL32_HL_SX2(acc20, acc30) + XT_SEL32_LH_SX2(acc20, acc30);

            vz_4 = AE_ZALIGN128(); AE_SASX2X2_IP(z0, z1, vz_4, pz_4);  AE_SA128POS_FP(vz_4, pz_4);
        }
        z = (float32_t *)pz_4;
    } else if ( M>=4 ){
        pz_4 = (xtfloatx4 *)z;

        __Pragma("loop_count min=1")
        for (m = 0; m < (M >> 2); m++)
        {
            px0 = (const xtfloatx4 *) (x);
            px1 = (const xtfloatx4 *) XT_ADDX4(N,(uintptr_t) px0);
            px2 = (const xtfloatx4 *) XT_ADDX4(N,(uintptr_t) px1);
            px3 = (const xtfloatx4 *) XT_ADDX4(N,(uintptr_t) px2);
            x   = (const float32_t *) XT_ADDX4(N,(uintptr_t) px3);
            py = (const xtfloatx4 *)(y);

            //vx0 = AE_LA128_PP(px0);
            vx1 = AE_LA128_PP(px1);
            vx2 = AE_LA128_PP(px2);
            vx3 = AE_LA128_PP(px3);

            vy = AE_LA128_PP(py);
            acc00 = acc01 = acc10 = acc11 = 0.0f;
            acc20 = acc21 = acc30 = acc31 = 0.0f;
            acc02 = acc03 = acc12 = acc13 = 0.0f;
            acc22 = acc23 = acc32 = acc33 = 0.0f;

            for (n = 0; n < (N >> 3); n++)
            {
                vx0 = AE_LA128_PP(px0); AE_LASX2X2_IP(x00, x01, vx0, px0);
                /*vx1 = AE_LA128_PP(px1);*/ AE_LASX2X2_IP(x10, x11, vx1, px1);
                AE_LASX2X2_IP(x20, x21, vx2, px2);
                AE_LASX2X2_IP(x30, x31, vx3, px3);

                AE_LASX2X2_IP(y0, y1, vy, py);

                MADDQ_S(acc00, acc10, x00, x10, y0);
                MADDQ_S(acc01, acc11, x01, x11, y1);

                MADDQ_S(acc20, acc30, x20, x30, y0);
                MADDQ_S(acc21, acc31, x21, x31, y1);

                AE_LASX2X2_IP(x00, x01, vx0, px0);
                AE_LASX2X2_IP(x10, x11, vx1, px1);
                AE_LASX2X2_IP(x20, x21, vx2, px2);
                AE_LASX2X2_IP(x30, x31, vx3, px3);

                vy = AE_LA128_PP(py); AE_LASX2X2_IP(y0, y1, vy, py);

                MADDQ_S(acc02, acc12, x00, x10, y0);
                MADDQ_S(acc03, acc13, x01, x11, y1);

                MADDQ_S(acc22, acc32, x20, x30, y0);
                MADDQ_S(acc23, acc33, x21, x31, y1);
            }

            acc00 += acc02;  acc01 += acc03;
            acc10 += acc12;  acc11 += acc13;
            acc20 += acc22;  acc21 += acc23;
            acc30 += acc32;  acc31 += acc33;

            vx0 = AE_LA128_PP(px0); AE_LASX2X2_IP(x00, x01, vx0, px0);
            vx1 = AE_LA128_PP(px1); AE_LASX2X2_IP(x10, x11, vx1, px1);
            AE_LASX2X2_IP(x20, x21, vx2, px2);
            AE_LASX2X2_IP(x30, x31, vx3, px3);

            AE_LASX2X2_IP(y0, y1, vy, py);

            MOVF_SX2(y0 , zero_v, bmask0);
            MOVF_SX2(y1 , zero_v, bmask1);

            MOVF_SX2(x00, zero_v, bmask0);
            MOVF_SX2(x10, zero_v, bmask0);
            MOVF_SX2(x01, zero_v, bmask1);
            MOVF_SX2(x11, zero_v, bmask1);
            MOVF_SX2(x20, zero_v, bmask0);
            MOVF_SX2(x30, zero_v, bmask0);
            MOVF_SX2(x21, zero_v, bmask1);
            MOVF_SX2(x31, zero_v, bmask1);

            MADDQ_S(acc00, acc10, x00, x10, y0);
            MADDQ_S(acc01, acc11, x01, x11, y1);
            MADDQ_S(acc20, acc30, x20, x30, y0);
            MADDQ_S(acc21, acc31, x21, x31, y1);

            acc00 = acc00 + acc01;
            acc10 = acc10 + acc11;
            acc20 = acc20 + acc21;
            acc30 = acc30 + acc31;

            z0 = XT_SEL32_HL_SX2(acc00, acc10) + XT_SEL32_LH_SX2(acc00, acc10);
            z1 = XT_SEL32_HL_SX2(acc20, acc30) + XT_SEL32_LH_SX2(acc20, acc30);

            vz_4 = AE_ZALIGN128();  AE_SASX2X2_IP(z0, z1, vz_4, pz_4); AE_SA128POS_FP(vz_4, pz_4);
        }
        z = (float32_t *)pz_4;
    }

    if (M&2)
    {
        px0 = (const xtfloatx4 *)(x);
        px1 = (const xtfloatx4 *)XT_ADDX4( N, (uintptr_t)px0);
        x   = (const float32_t *)XT_ADDX4( N, (uintptr_t)px1);
        py  = (const xtfloatx4 *)(y);

        pz_2 = (xtfloatx2 *)(z);

        vx1 = AE_LA128_PP(px1);
        vx0 = AE_LA128_PP(px0);
        vy = AE_LA128_PP(py);
        acc00 = acc01 = acc10 = acc11 = (xtfloatx2)0.0f;

        for (n = 0; n < (N >> 2); n++)
        {
            AE_LASX2X2_IP(x00, x01, vx0, px0);
            AE_LASX2X2_IP(x10, x11, vx1, px1);

            AE_LASX2X2_IP(y0, y1, vy, py);

            /* perform multiplications */
            MADDQ_S(acc00, acc10, x00, x10, y0);
            MADDQ_S(acc01, acc11, x01, x11, y1);
        }
        AE_LASX2X2_IP(x00, x01, vx0, px0);
        AE_LASX2X2_IP(x10, x11, vx1, px1);
        AE_LASX2X2_IP(y0, y1, vy, py);

        MOVF_SX2(x00, zero_v, bmask0);
        MOVF_SX2(x10, zero_v, bmask0);
        MOVF_SX2(y0 , zero_v, bmask0);

        MOVF_SX2(x01, zero_v, bmask1);
        MOVF_SX2(x11, zero_v, bmask1);
        MOVF_SX2( y1, zero_v, bmask1);

        /* perform multiplications */
        MADDQ_S(acc00, acc10, x00, x10, y0);
        MADDQ_S(acc01, acc11, x01, x11, y1);

        acc00 = acc00 + acc01;
        acc10 = acc10 + acc11;

        z0 = XT_SEL32_HL_SX2(acc00, acc10) + XT_SEL32_LH_SX2(acc00, acc10);

        vz_2 = AE_ZALIGN64(); AE_SASX2IP(z0, vz_2, pz_2); AE_SASX2POSFP(vz_2, pz_2);
        z = (float32_t *)pz_2;
    }
    

    /* Compute last (M%2) output element */
    if(M&1)
    {
        px1 = (const xtfloatx4 *)(x);
        py = (const xtfloatx4 *)(y);

        vy = AE_LA128_PP(py);
        vx1 = AE_LA128_PP(px1);
        acc00 = acc01 = (xtfloatx2)0.0f;

        for (n = 0; n < (N >> 2); n++)
        {
            AE_LASX2X2_IP(x00, x01, vx1, px1);
            AE_LASX2X2_IP(y0, y1, vy, py);

            MADD_SX2X2(acc00, acc01, x00, x01, y0, y1);
        }

        AE_LASX2X2_IP(x00, x01, vx1, px1);
        AE_LASX2X2_IP(y0, y1, vy, py);

        MOVF_SX2(y0 , zero_v, bmask0);
        MOVF_SX2(y1 , zero_v, bmask1);
        
        MOVF_SX2(x00, zero_v, bmask0);
        MOVF_SX2(x01, zero_v, bmask1);

        MADD_SX2X2(acc00, acc01, x00, x01, y0, y1);

        acc00 = acc00 + acc01;

        z0_ = XT_RADD_SX2(acc00);

        XT_SSIP(z0_, castxcc(xtfloat,z), sizeof(float32_t));
    }
} /* mtx_vecmpyf() */

#elif (HAVE_FPU)

void mtx_vecmpyf( float32_t * z, const float32_t * x,  const float32_t * y, int M, int N )
{
    int m, n;
    xtfloat * restrict pZ;
    const xtfloat * restrict pX;
    const xtfloat * restrict pY;

    NASSERT(x);
    NASSERT(y);
    NASSERT(z);
    NASSERT((z != x) && (z != y));

    if (M <= 0)    return;
    // If N<=0 then clear output vector and return 
    pZ=(xtfloat*)z;
    if (N<0)
    {
        for ( m=0; m<M; m++) XT_SSIP(XT_CONST_S(0),pZ,sizeof(xtfloat));
        return;
    }
    for ( m=0; m<(M&~3); m+=4)
    {
        xtfloat a0,a1,a2,a3,x0,x1,x2,x3,y0;
        a0=a1=a2=a3=XT_CONST_S(0);
        pX=(const xtfloat*)(x+m*N);
        pY=(const xtfloat*)(y);
        for ( n=0; n<N; n++ )
        {
            XT_LSIP(y0,pY,sizeof(xtfloat));
            x1=XT_LSX(pX,1*N*sizeof(xtfloat));
            x2=XT_LSX(pX,2*N*sizeof(xtfloat));
            x3=XT_LSX(pX,3*N*sizeof(xtfloat));
            XT_LSIP(x0,pX,sizeof(xtfloat));
            XT_MADD_S(a0,x0,y0);
            XT_MADD_S(a1,x1,y0);
            XT_MADD_S(a2,x2,y0);
            XT_MADD_S(a3,x3,y0);
        }
        XT_SSIP(a0,pZ,sizeof(xtfloat));
        XT_SSIP(a1,pZ,sizeof(xtfloat));
        XT_SSIP(a2,pZ,sizeof(xtfloat));
        XT_SSIP(a3,pZ,sizeof(xtfloat));
    }
    if (M&2)
    {
        xtfloat a0,a1,x0,x1,y0;
        a0=a1=XT_CONST_S(0);
        pX=(const xtfloat*)(x+m*N);
        pY=(const xtfloat*)(y);
        for ( n=0; n<N; n++ )
        {
            XT_LSIP(y0,pY,sizeof(xtfloat));
            x1=XT_LSX(pX,1*N*sizeof(xtfloat));
            XT_LSIP(x0,pX,sizeof(xtfloat));
            XT_MADD_S(a0,x0,y0);
            XT_MADD_S(a1,x1,y0);
        }
        XT_SSIP(a0,pZ,sizeof(xtfloat));
        XT_SSIP(a1,pZ,sizeof(xtfloat));
        m+=2;
    }
    if (M&1)
    {
        xtfloat a0,x0,y0;
        a0=XT_CONST_S(0);
        pX=(const xtfloat*)(x+m*N);
        pY=(const xtfloat*)(y);
        for ( n=0; n<N; n++ )
        {
            XT_LSIP(y0,pY,sizeof(xtfloat));
            XT_LSIP(x0,pX,sizeof(xtfloat));
            XT_MADD_S(a0,x0,y0);
        }
        XT_SSIP(a0,pZ,sizeof(xtfloat));
    }
}

#endif

