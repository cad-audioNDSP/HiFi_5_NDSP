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
  NatureDSP Signal Processing Library. IIR filters
    Kalman filter update for order 1
    Single precision floating point
    C code optimized for HiFi5 with VFPU/SFPU
  IntegrIT, 2006-2019
*/

/* Portable data types. */
#include "NatureDSP_types.h"
/* Common utility and macros declarations. */
#include "common.h"
#include "common_fpu.h"
/* Signal Processing Library API. */
#include "NatureDSP_Signal_iir.h"
#include "inff_tbl.h"

#define ALIGN_SIZE      (HIFI_SIMD_WIDTH)

/*-------------------------------------------------------------------------
  Kalman Filter Update
  For input matrices U (N-by-M), H (M-by-N) and R (M-by-M), evaluate 
  U*(H*U+R)^-1, and store the resulting M-by-N matrix to the output 
  argument K.

  Precision: 
  32x32   32-bit data, 32-bit coefficients
  f       single precision floating point

  Parameter:
  M       order, equals to 1
  Input:
  N       number of states
  U[N*M]  state matrix
  H[M*N]  measurement matrix
  R[M*M]  noise estimate (measurement covariance matrix)
  qK      fixed point position for the output matrix K
  qU      fixed point position for the input matrix U
  qH      fixed point position for the input matrix H
  qR      fixed point position for the input matrix R
  Output:
  K[N*M]  Kalman gain matrix
  Temporary:
  pScr    scratch memory area of size specified by the 
          corresponding scratch allocation function (in bytes) 

  Restrictions:
  U,R,H,K  must not overlap and must be aligned by 16-bytes
  N        multiple of 32
-------------------------------------------------------------------------*/
#if (!HAVE_VFPU && !HAVE_FPU) 
DISCARD_FUN(void, kalmanupd1f, ( void * pScr, 
                                 float32_t * K, 
                           const float32_t * U, 
                           const float32_t * H,
                           const float32_t * R,
                           int N ));
#elif (HAVE_VFPU)
void kalmanupd1f( void * pScr, 
                  float32_t * K, 
            const float32_t * U, 
            const float32_t * H,
            const float32_t * R,
            int N )
{
    const xtfloatx4  * restrict pX = (const xtfloatx4  *)U;
    const xtfloatx4  * restrict pY = (const xtfloatx4  *)H;
          xtfloatx4  * restrict pZ = (      xtfloatx4  *)K;

    xtfloatx2 f,s,r;
    xtbool2 bzero,bsmall,binf;
    static const union {int32_t i; float32_t f;} invrealminf= {0x7e800000} ;// 1.f/realminf.f
    int n;
    NASSERT_ALIGN(K, ALIGN_SIZE);
    NASSERT_ALIGN(U, ALIGN_SIZE);
    NASSERT_ALIGN(H, ALIGN_SIZE);
    NASSERT_ALIGN(R, ALIGN_SIZE);
    NASSERT(0==(N%32));
    if (N<=0) return;

    xtfloatx2 a0,a1,a2,a3,a4,a5,a6,a7;
    a0=a1=a2=a3=CONST_SX2(0);
    a4=a5=a6=a7=CONST_SX2(0);
    for (n=0; n<(N>>4); n++) 
    {
        xtfloatx2 x0,x1,x2,x3,y0,y1,y2,y3;
        AE_LSX2X2_IP(x0,x1,pX,sizeof(xtfloatx4));
        AE_LSX2X2_IP(x2,x3,pX,sizeof(xtfloatx4));
        AE_LSX2X2_IP(y0,y1,pY,sizeof(xtfloatx4));
        AE_LSX2X2_IP(y2,y3,pY,sizeof(xtfloatx4));
        MADD_SX2X2(a0,a1,x0,x1,y0,y1);
        MADD_SX2X2(a2,a3,x2,x3,y2,y3);
        AE_LSX2X2_IP(x0,x1,pX,sizeof(xtfloatx4));
        AE_LSX2X2_IP(x2,x3,pX,sizeof(xtfloatx4));
        AE_LSX2X2_IP(y0,y1,pY,sizeof(xtfloatx4));
        AE_LSX2X2_IP(y2,y3,pY,sizeof(xtfloatx4));
        MADD_SX2X2(a4,a5,x0,x1,y0,y1);
        MADD_SX2X2(a6,a7,x2,x3,y2,y3);
    }
    a0=XT_ADD_SX2(a0,a4);
    a1=XT_ADD_SX2(a1,a5);
    a2=XT_ADD_SX2(a2,a6);
    a3=XT_ADD_SX2(a3,a7);
    a0=XT_ADD_SX2(a0,a1);
    a2=XT_ADD_SX2(a2,a3);
    a0=XT_ADD_SX2(a0,a2);
    a0=ADD_HL_LH_S(a0,a0);
    f=XT_LOW_S(a0);
    f=XT_ADD_S(f,R[0]);
    bzero=XT_OEQ_S(f,0);
    bsmall=XT_OLE_S(XT_ABS_S(f),realminf.f);
    binf  =XT_OEQ_S(XT_ABS_S(f),plusInff.f);
    s=XT_CONST_S(1); 
    MOVT_SX2 (s,invrealminf.f,bsmall);
    f=XT_MUL_S(f,s);
    r=XT_RECIP_S(f);
    MOVT_SX2 (r,plusInff.f   ,bzero);
    MOVT_SX2 (r,XT_CONST_S(0),binf);
    pX = (const xtfloatx4  *)U;
    for (n=0; n<(N>>4); n++) 
    {
        xtfloatx2 x0,x1,x2,x3;
        AE_LSX2X2_IP(x0,x1,pX,sizeof(xtfloatx4));
        AE_LSX2X2_IP(x2,x3,pX,sizeof(xtfloatx4));
        MULQ_S(x0,x1,x0,x1,s);
        MULQ_S(x2,x3,x2,x3,s);
        MULQ_S(x0,x1,x0,x1,r);
        MULQ_S(x2,x3,x2,x3,r);
        AE_SSX2X2_IP(x0,x1,pZ,sizeof(xtfloatx4));
        AE_SSX2X2_IP(x2,x3,pZ,sizeof(xtfloatx4));
        AE_LSX2X2_IP(x0,x1,pX,sizeof(xtfloatx4));
        AE_LSX2X2_IP(x2,x3,pX,sizeof(xtfloatx4));
        MULQ_S(x0,x1,x0,x1,s);
        MULQ_S(x2,x3,x2,x3,s);
        MULQ_S(x0,x1,x0,x1,r);
        MULQ_S(x2,x3,x2,x3,r);
        AE_SSX2X2_IP(x0,x1,pZ,sizeof(xtfloatx4));
        AE_SSX2X2_IP(x2,x3,pZ,sizeof(xtfloatx4));
    }
} /* kalmanupd1f() */
#elif (HAVE_FPU)
void kalmanupd1f( void * pScr, 
                  float32_t * K, 
            const float32_t * U, 
            const float32_t * H,
            const float32_t * R,
            int N )
{
    static const union {int32_t i; float32_t f;} invrealminf= {0x7e800000} ;// 1.f/realminf.f
    const xtfloat  * restrict pX = (const xtfloat  *)U;
    const xtfloat  * restrict pY = (const xtfloat  *)H;
          xtfloat  * restrict pZ = (      xtfloat  *)K;
    xtfloat f,r,s;
    int n;
    NASSERT_ALIGN(K, ALIGN_SIZE);
    NASSERT_ALIGN(U, ALIGN_SIZE);
    NASSERT_ALIGN(H, ALIGN_SIZE);
    NASSERT_ALIGN(R, ALIGN_SIZE);
    NASSERT(0==(N%32));
    if (N<=0) return ;

    {
      xtfloat acc0, acc1;
      xtbool bzero,bsmall,binf;
      int n;
      acc0 = R[0]; acc1 = XT_CONST_S(0);
      for (n = 0; n<(N&~1); n+=2)
      {
        xtfloat x0,y0;
        XT_LSIP(x0, pX, sizeof(xtfloat));
        XT_LSIP(y0, pY, sizeof(xtfloat));
        XT_MADD_S(acc0,x0,y0);
        XT_LSIP(x0, pX, sizeof(xtfloat));
        XT_LSIP(y0, pY, sizeof(xtfloat));
        XT_MADD_S(acc1,x0,y0);
      }

      f=XT_ADD_S(acc0,acc1);
      bzero=XT_OEQ_S(f,0);
      bsmall=XT_OLE_S(XT_ABS_S(f),realminf.f);
      binf  =XT_OEQ_S(XT_ABS_S(f),plusInff.f);
      s=XT_CONST_S(1); 
      XT_MOVT_S(s,invrealminf.f,bsmall);
      f=XT_MUL_S(f,s);
      r=XT_RECIP_S(f);
      XT_MOVT_S(r,plusInff.f   ,bzero);
      XT_MOVT_S(r,XT_CONST_S(0),binf);
    }
    pX = (const xtfloat  *)U;
    for ( n=0; n<(N>>2); n++) 
    {
        xtfloat x0,x1,x2,x3;
        XT_LSIP(x0, pX, sizeof(xtfloat));
        XT_LSIP(x1, pX, sizeof(xtfloat));
        XT_LSIP(x2, pX, sizeof(xtfloat));
        XT_LSIP(x3, pX, sizeof(xtfloat));
        XT_SSIP(XT_MUL_S(XT_MUL_S(x0,s),r),pZ, sizeof(xtfloat));
        XT_SSIP(XT_MUL_S(XT_MUL_S(x1,s),r),pZ, sizeof(xtfloat));
        XT_SSIP(XT_MUL_S(XT_MUL_S(x2,s),r),pZ, sizeof(xtfloat));
        XT_SSIP(XT_MUL_S(XT_MUL_S(x3,s),r),pZ, sizeof(xtfloat));
    }
} /* kalmanupd1f() */
#endif

/* Returns: size of scratch memory area, in bytes. */
size_t kalmanupd1f_getScratchSize(int N)
{
    NASSERT(0==(N%32));
    (void)N;
    return 0;
}
