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
void mtx_mpyt32x32 (  void* pScr,
                     int32_t* restrict z,
               const int32_t* restrict x,
               const int32_t* restrict yt,
               int M, int N, int P, int lsh )
{
    const ae_int32x4 * restrict px;
    const ae_int32x4 * restrict px0;
    const ae_int32x4 * restrict px1;
    const ae_int32x4 * restrict px2;
    const ae_int32x4 * restrict px3;
    const ae_int32x4 * restrict py;
    const ae_int32x4 * restrict py0;
    const ae_int32x4 * restrict py1;
    ae_int32x4 * restrict pz_;
    ae_int32x2 * restrict pz;
    ae_int32x2 * restrict pz0;
    ae_int32x2 * restrict pz1;
    ae_valignx2 ay0, ay1, ax0, ax1;
    ae_valign az0, az1;
    ae_int32x2 vx0,vx1,vx2,vx3,  vx4,vx5,vx6,vx7, vy0, vy1,vy2,vy3, vt0, vt1,vt2,vt3, vz0,vz1;
    ae_f64   ACC00, ACC01, ACC10, ACC11;
    xtbool2 bmask0, bmask1;
    int n, m, p;

    if ((M <= 0) || (P <= 0)) return;
    if (N <= 0)
    {
        for (m = 0; m < M * P; m++) z[m] = 0;
        return;
    }

    bmask0 = AE_MOVBA2((N & 2) + ((N & 2) >> 1) + 2 * (int)((N & 3) == 1));  // 2*(N%4) if N%4<2, else 3 
    bmask1 = AE_MOVBA2(((int)((N & 3) == 3)) << 1);  // 2 if (N%4)=3, else 0

    for (p = 0; p < (P >> 1); p++)
    {
        py0 = (const ae_int32x4 *)yt;
        py1 = (const ae_int32x4 *)((uintptr_t)py0 + N*sizeof(int32_t));
        pz_ = (ae_int32x4 *)pScr;
        ay0 = AE_LA128_PP(py0);
        ay1 = AE_LA128_PP(py1);
        for (n = 0; n < (N >> 2); n++)
        {
            AE_LA32X2X2_IP(vt0, vt1, ay0, py0);
            AE_LA32X2X2_IP(vt2, vt3, ay1, py1);
            AE_S32X2X2_IP(vt0, vt1, pz_, 2*sizeof(ae_int32x2));
            AE_S32X2X2_IP(vt2, vt3, pz_, 2*sizeof(ae_int32x2));
        }
        AE_LA32X2X2_IP(vt0, vt1, ay0, py0);
        AE_LA32X2X2_IP(vt2, vt3, ay1, py1);
        AE_MOVF32X2(vt0, AE_ZERO32(), bmask0);
        AE_MOVF32X2(vt1, AE_ZERO32(), bmask1);
        AE_MOVF32X2(vt2, AE_ZERO32(), bmask0);
        AE_MOVF32X2(vt3, AE_ZERO32(), bmask1);
        AE_S32X2X2_IP(vt0, vt1, pz_, 2 * sizeof(ae_int32x2));
        AE_S32X2X2_IP(vt2, vt3, pz_, 2 * sizeof(ae_int32x2));

        pz = (ae_int32x2 *)z;
        pz0 = (ae_int32x2 *)pz;
        pz1 = (ae_int32x2 *)((uintptr_t)pz0+P*sizeof(int32_t));
        px = (const ae_int32x4 *)(x);

        for (m = 0; m < (M >> 2); m++)
        {
            ae_f64 ACC20,ACC21,ACC30,ACC31;
            ae_valignx2 ax2,ax3;
            px0 = (const ae_int32x4 *)px;
            px1 = (const ae_int32x4 *)XT_ADDX2(N*2,(uintptr_t)px);
            px2 = (const ae_int32x4 *)XT_ADDX4(N*2,(uintptr_t)px);
            px3 = (const ae_int32x4 *)XT_ADDX2(N*2,(uintptr_t)px2);
            px  = (const ae_int32x4 *)XT_ADDX8(N*2,(uintptr_t)px);
            ax0 = AE_LA128_PP(px0);
            ax1 = AE_LA128_PP(px1);
            ax2 = AE_LA128_PP(px2);
            ax3 = AE_LA128_PP(px3);
            py = (const ae_int32x4 *)pScr;
            py0 = (const ae_int32x4 *)py;

            AE_MOVDX2(ACC00,ACC01,0,0); AE_MOVDX2(ACC10,ACC11,0,0);
            AE_MOVDX2(ACC20,ACC21,0,0); AE_MOVDX2(ACC30,ACC31,0,0);
            __Pragma("loop_count min=1")
            for (n = 0; n < N; n += 4)
            {
                /* load data from 'x' */
                AE_LA32X2X2_IP(vx0, vx1, ax0, px0);
                AE_LA32X2X2_IP(vx2, vx3, ax1, px1);
                AE_LA32X2X2_IP(vx4, vx5, ax2, px2);
                AE_LA32X2X2_IP(vx6, vx7, ax3, px3);
                /* load data from 'y' */
                AE_L32X2X2_I (vy2, vy3, py0, 1 * sizeof(ae_int32x4));
                AE_L32X2X2_XP(vy0, vy1, py0, 2 * sizeof(ae_int32x4));
                /* perform multiplications */
                AE_MULAAF2D32RA_HH_LL(ACC00, ACC10, vx0, vx2, vy0, vy0);
                AE_MULAAF2D32RA_HH_LL(ACC01, ACC11, vx0, vx2, vy2, vy2);
                AE_MULAAF2D32RA_HH_LL(ACC00, ACC10, vx1, vx3, vy1, vy1);
                AE_MULAAF2D32RA_HH_LL(ACC01, ACC11, vx1, vx3, vy3, vy3);
                AE_MULAAF2D32RA_HH_LL(ACC20, ACC30, vx4, vx6, vy0, vy0);
                AE_MULAAF2D32RA_HH_LL(ACC21, ACC31, vx4, vx6, vy2, vy2);
                AE_MULAAF2D32RA_HH_LL(ACC20, ACC30, vx5, vx7, vy1, vy1);
                AE_MULAAF2D32RA_HH_LL(ACC21, ACC31, vx5, vx7, vy3, vy3);
            }
            /* save values */
            az0 = AE_ZALIGN64();
            az1 = AE_ZALIGN64();
            AE_SA32X2_IP(AE_TRUNCA32X2F64S(ACC00, ACC01, 16 + lsh), az0, pz0);
            AE_SA32X2_IP(AE_TRUNCA32X2F64S(ACC10, ACC11, 16 + lsh), az1, pz1);
            AE_SA64POS_FP(az0, pz0);
            AE_SA64POS_FP(az1, pz1);
            pz0 = (ae_int32x2 *)XT_ADDX8(P - 1, (uintptr_t)pz0);
            pz1 = (ae_int32x2 *)XT_ADDX8(P - 1, (uintptr_t)pz1);
            AE_SA32X2_IP(AE_TRUNCA32X2F64S(ACC20, ACC21, 16 + lsh), az0, pz0);
            AE_SA32X2_IP(AE_TRUNCA32X2F64S(ACC30, ACC31, 16 + lsh), az1, pz1);
            AE_SA64POS_FP(az0, pz0);
            AE_SA64POS_FP(az1, pz1);
            pz0 = (ae_int32x2 *)XT_ADDX8(P - 1, (uintptr_t)pz0);
            pz1 = (ae_int32x2 *)XT_ADDX8(P - 1, (uintptr_t)pz1);

        }
        if (M&2)
        {
            px0 = (const ae_int32x4 *)px;
            px1 = (const ae_int32x4 *)((uintptr_t)px+N*sizeof(int32_t));
            ax0 = AE_LA128_PP(px0);
            ax1 = AE_LA128_PP(px1);
            py = (const ae_int32x4 *)pScr;
            py0 = (const ae_int32x4 *)py;
            py1 = (const ae_int32x4 *)((uintptr_t)py+4*sizeof(int32_t));

            AE_MOVDX2(ACC00,ACC01,0,0); AE_MOVDX2(ACC10,ACC11,0,0);
            __Pragma("loop_count min=1")
            for (n = 0; n < N; n += 4)
            {
                /* load data from 'x' */
                AE_LA32X2X2_IP(vx0, vx1, ax0, px0);
                AE_LA32X2X2_IP(vx2, vx3, ax1, px1);
                /* load data from 'y' */
                AE_L32X2X2_IP(vy0, vy1, py0, 4 * sizeof(ae_int32x2));
                AE_L32X2X2_IP(vy2, vy3, py1, 4 * sizeof(ae_int32x2));
                /* perform multiplications */
                AE_MULAAF2D32RA_HH_LL(ACC00, ACC10, vx0, vx2, vy0, vy0);
                AE_MULAAF2D32RA_HH_LL(ACC01, ACC11, vx0, vx2, vy2, vy2);
                AE_MULAAF2D32RA_HH_LL(ACC00, ACC10, vx1, vx3, vy1, vy1);
                AE_MULAAF2D32RA_HH_LL(ACC01, ACC11, vx1, vx3, vy3, vy3);
            }
            /* format values */
            vz0 = AE_TRUNCA32X2F64S(ACC00, ACC01, 16 + lsh);
            vz1 = AE_TRUNCA32X2F64S(ACC10, ACC11, 16 + lsh);
            /* save values */
            az0 = AE_ZALIGN64();
            az1 = AE_ZALIGN64();
            AE_SA32X2_IP(vz0, az0, pz0);
            AE_SA32X2_IP(vz1, az1, pz1);
            AE_SA64POS_FP(az0, pz0);
            AE_SA64POS_FP(az1, pz1);
            pz0 = (ae_int32x2 *)XT_ADDX8(P - 1, (uintptr_t)pz0);
            pz1 = (ae_int32x2 *)XT_ADDX8(P - 1, (uintptr_t)pz1);

            px = (const ae_int32x4 *)XT_ADDX8(N, (uintptr_t)px);
        }
        if (M & 1)
        {
            px0 = (const ae_int32x4 *)px;
            ax0 = AE_LA128_PP(px0);
            py = (const ae_int32x4 *)pScr;

            ACC00 = ACC01 = AE_ZERO64();
            __Pragma("loop_count min=1")
            for (n = 0; n < N; n += 4)
            {
                /* load data from 'x' */
                AE_LA32X2X2_IP(vx0, vx1, ax0, px0);
                /* load data from 'y' */
                AE_L32X2X2_IP(vy0, vy1, py, 2*sizeof(ae_int32x2));
                AE_L32X2X2_IP(vy2, vy3, py, 2*sizeof(ae_int32x2));
                /* perform multiplications */
                AE_MULAAF2D32RA_HH_LL(ACC00, ACC01, vx0, vx0, vy0, vy2);
                AE_MULAAF2D32RA_HH_LL(ACC00, ACC01, vx1, vx1, vy1, vy3);
            }
            /* format values */
            vz0 = AE_TRUNCA32X2F64S(ACC00, ACC01, 16 + lsh);
            /* save values */
            az0 = AE_ZALIGN64();
            AE_SA32X2_IP(vz0, az0, pz0);
            AE_SA64POS_FP(az0, pz0);
        }

        yt += 2 * N;
        z += 2;
    }

    if (P & 1)
    {
        py0 = (const ae_int32x4 *)yt;
        pz_ = (ae_int32x4 *)pScr;
        ay0 = AE_LA128_PP(py0);
        for (n = 0; n < (N >> 2); n++)
        {
            AE_LA32X2X2_IP(vt0, vt1, ay0, py0);
            AE_S32X2X2_IP(vt0, vt1, pz_, 2 * sizeof(ae_int32x2));
        }
        AE_LA32X2X2_IP(vt0, vt1, ay0, py0);
        AE_MOVF32X2(vt0, AE_ZERO32(), bmask0);
        AE_MOVF32X2(vt1, AE_ZERO32(), bmask1);
        AE_S32X2X2_IP(vt0, vt1, pz_, 2 * sizeof(ae_int32x2));

        pz = (ae_int32x2 *)z;
        pz0 = (ae_int32x2 *)pz;
        pz1 = (ae_int32x2 *)XT_ADDX4(P, (uintptr_t)pz0);
        px = (const ae_int32x4 *)(x);

        for (m = 0; m < (M >> 1); m++)
        {
            px0 = (const ae_int32x4 *)px;
            px1 = (const ae_int32x4 *)XT_ADDX4(N, (uintptr_t)px);
            ax0 = AE_LA128_PP(px0);
            ax1 = AE_LA128_PP(px1);
            py = (const ae_int32x4 *)pScr;

            ACC00 = ACC10 = AE_ZERO64();
            __Pragma("loop_count min=1")
            for (n = 0; n < N; n += 4)
            {
                /* load data from 'x' */
                AE_LA32X2X2_IP(vx0, vx1, ax0, px0);
                AE_LA32X2X2_IP(vx2, vx3, ax1, px1);
                /* load data from 'y' */
                AE_L32X2X2_IP(vy0, vy1, py, 2*sizeof(ae_int32x2));
                /* perform multiplications */
                AE_MULAAF2D32RA_HH_LL(ACC00, ACC10, vx0, vx2, vy0, vy0);
                AE_MULAAF2D32RA_HH_LL(ACC00, ACC10, vx1, vx3, vy1, vy1);
            }
            /* format values */
            vz0 = AE_TRUNCA32X2F64S(ACC00, ACC00, 16 + lsh);
            vz1 = AE_TRUNCA32X2F64S(ACC10, ACC10, 16 + lsh);
            /* save values */
            AE_S32_L_I(vz0, (ae_int32 *)pz0, 0);
            AE_S32_L_I(vz1, (ae_int32 *)pz1, 0);
            pz0 = (ae_int32x2 *)XT_ADDX8(P, (uintptr_t)pz0);
            pz1 = (ae_int32x2 *)XT_ADDX8(P, (uintptr_t)pz1);

            px = (const ae_int32x4 *)XT_ADDX8(N, (uintptr_t)px);
        }
        if (M & 1)
        {
            px0 = (const ae_int32x4 *)px;
            ax0 = AE_LA128_PP(px0);
            py = (const ae_int32x4 *)pScr;

            ACC00 = ACC01 = AE_ZERO64();

            /* pre-load data from 'x' */
            AE_LA32X2X2_IP(vx0, vx1, ax0, px0);
            /* pre-load data from 'y' */
            AE_L32X2X2_IP(vy0, vy1, py, 2*sizeof(ae_int32x2));
            for (n = 0; n < N - 4; n += 4)
            {
                /* perform multiplications */
                AE_MULAAF2D32RA_HH_LL(ACC00, ACC01, vx0, vx1, vy0, vy1);
                /* load data from 'x' */
                AE_LA32X2X2_IP(vx0, vx1, ax0, px0);
                /* load data from 'y' */
                AE_L32X2X2_IP(vy0, vy1, py, 2*sizeof(ae_int32x2));
            }
            /* perform multiplications */
            AE_MULAAF2D32RA_HH_LL(ACC00, ACC01, vx0, vx1, vy0, vy1);
            /* format values */
            vz0 = AE_TRUNCA32X2F64S(ACC00+ACC01, ACC00+ACC01, 16 + lsh);
            /* save values */
            AE_S32_L_I(vz0, (ae_int32 *)pz0, 0);
        }
    }
}

size_t mtx_mpyt32x32_getScratchSize ( int M, int N, int P )
{
	return N<=0? 0: ((((N)+4))*2*sizeof(int32_t));
} /* mtx_mpyt32x32_getScratchSize() */
