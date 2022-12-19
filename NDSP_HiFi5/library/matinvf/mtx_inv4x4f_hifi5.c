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
 * Real Matrix Inversion
 * Optimized code for HiFi5
  IntegrIT, 2006-2019
 */

/* Cross-platform data type definitions. */
#include "NatureDSP_types.h"
/* Common helper macros. */
#include "common.h"
/* Matrix functions */
#include "NatureDSP_Signal_matinv.h"
#include "common_fpu.h"


#if (HAVE_VFPU)
#define SZ_F32 (sizeof(float32_t))
#define SZ_2F32 (2*SZ_F32)
/*-------------------------------------------------------------------------
  These functions implement in-place matrix inversion by Gauss elimination 
  with full pivoting
  NOTE: user may detect "invalid" or "divide-by-zero" exception in the CPU 
  flags which MAY indicate that inversion results are not accurate. Also 
  it's responsibility of the user to provide valid input matrix for 
  inversion.
  Fixed point version takes representation of input matrix and forms 
  representation of output matrix with proper scaling.

  Precision: 
  f     floating point
  32x32 32-bit input, 32-bit output

  Input:
  x[N*N]      input matrix
  qX          input matrix representation (for fixed point API only)
  Output:
  x[N*N]      result
  Temporary:
  pScr        scratch memory. Size in bytes is defined by corresponding 
              scratch allocation function 
  N is 2,3,4,6,8,10

  Returned value: floating functions return none, fixed point functions 
                  return fixed-point representation of inverted matrix
  Restrictions:
  none
-------------------------------------------------------------------------*/
void mtx_inv4x4f(void* pScr,float32_t* x)
{
  xtfloatx2 *restrict pk;
  xtfloatx2 *restrict pmax;
  xtfloatx2 *restrict pA;
  xtfloatx2 *restrict pA0;
  xtfloatx2 *restrict pA1;
  xtfloatx2 *restrict pA2;
  xtfloatx2 R, C0, C1, C2;
  xtfloatx2 A00, A01, A02, A03;
  xtfloatx2 A10, A11, A12, A13;
  xtfloatx2 A20, A21, A22, A23;
  xtfloatx2 Amax0, Amax1, Amax2, Amax3,
              Ak0,   Ak1,   Ak2,   Ak3;
  ae_valignx2 vX;
  int n,k;
  float32_t *A;//float32_t ALIGN(16) A[32];
    NASSERT_ALIGN16(pScr);
    A = (float32_t *)pScr;
    /* Copy the matrix to buffer and
   * initialize identity matrix:
   *
   * x00 x01 x02 x03 1.0 0.0 0.0 0.0
   * x10 x11 x12 x13 0.0 1.0 0.0 0.0
   * x20 x21 x22 x23 0.0 0.0 1.0 0.0
   * x30 x31 x32 x33 0.0 0.0 0.0 1.0
   */
  C0 = (xtfloatx2)0.0f;
  C1 = (xtfloatx2)1.0f;
  pA0 = (xtfloatx2 *)(x);
  pA1 = (xtfloatx2 *)(A);
  {
    xtfloatx2 C01, C10;
    C01 = XT_SEL32_HH_SX2(C0, C1);
    C10 = XT_SEL32_HH_SX2(C1, C0);
    /* Load input matrix */
    vX = AE_LA128_PP(pA0);
    AE_LASX2X2_IP(A00,A01, vX, castxcc(xtfloatx4,pA0));
    AE_LASX2X2_IP(A02,A03, vX, castxcc(xtfloatx4,pA0));
    AE_LASX2X2_IP(A10,A11, vX, castxcc(xtfloatx4,pA0));
    AE_LASX2X2_IP(A12,A13, vX, castxcc(xtfloatx4,pA0));

    /* Save input and identity matrix */
    /* 1-st row */
    AE_SSX2X2_IP(A00,A01, castxcc(xtfloatx4,pA1), sizeof(xtfloatx4));
    AE_SSX2X2_IP(C10,C0 , castxcc(xtfloatx4,pA1), sizeof(xtfloatx4));
    /* 2-nd row */
    AE_SSX2X2_IP(A02,A03, castxcc(xtfloatx4,pA1), sizeof(xtfloatx4));
    AE_SSX2X2_IP(C01,C0 , castxcc(xtfloatx4,pA1), sizeof(xtfloatx4));
    /* 3-rd row */
    AE_SSX2X2_IP(A10,A11, castxcc(xtfloatx4,pA1), sizeof(xtfloatx4));
    AE_SSX2X2_IP(C0 ,C10, castxcc(xtfloatx4,pA1), sizeof(xtfloatx4));
    /* 4-th row */
    AE_SSX2X2_IP(A12,A13, castxcc(xtfloatx4,pA1), sizeof(xtfloatx4));
    AE_SSX2X2_IP(C0 ,C01, castxcc(xtfloatx4,pA1), sizeof(xtfloatx4));
  }

  /* Set bounds of the buffer */
  WUR_AE_CBEGIN0((uintptr_t)(A));
  WUR_AE_CEND0  ((uintptr_t)(A+32));

  pk  = (xtfloatx2 *)(A);
  pA0 = (xtfloatx2 *)(A+8);
  pA1 = (xtfloatx2 *)(A+16);
  pA2 = (xtfloatx2 *)(A+24);
  /* Gauss elimination */
  for(k=0; k<4; k++)
  {
    xtfloatx2 amax;
    unsigned int imax;
    /* pivoting */
    imax=k;
    amax=XT_CONST_S(0);
    /* find absolute max value in the k-th column */
    pA=(xtfloatx2*)&A[2*k*4+k];
    for(n=k; n<4; n++)
    {
        xtbool2 cond;
        xtfloatx2 t;
        //XT_LSXP(t,castxcc(xtfloat,pA),2*N*sizeof(xtfloat));
        {
            ae_int32x2 tt;
            AE_L32_XP(tt,castxcc(ae_int32,pA),2*4*sizeof(xtfloat));
            t=XT_AE_MOVXTFLOATX2_FROMINT32X2(tt);
        }
        BMAXNUMABS_SX2(cond,amax,amax,t);
        XT_MOVF_SX2(R  ,  t,cond);
        XT_MOVEQZ(imax, n, AE_MOVAB2(cond));
    }
    R = XT_RECIP_SX2(R);
    
    /* swap k-th and imax-th rows */
    pmax = (xtfloatx2 *)(A+imax*8);
    AE_LSX2X2_I (Ak0,Ak1,(const xtfloatx4*)pmax, 0*sizeof(xtfloatx4));
    AE_LSX2X2_I (Ak2,Ak3,(const xtfloatx4*)pmax, 1*sizeof(xtfloatx4));
    AE_LSX2X2_I (Amax0,Amax1,(const xtfloatx4*)pk, 0*sizeof(xtfloatx4));
    AE_LSX2X2_I (Amax2,Amax3,(const xtfloatx4*)pk, 1*sizeof(xtfloatx4));
    AE_SSX2X2_IP(Amax0,Amax1,castxcc(xtfloatx4,pmax), sizeof(xtfloatx4));
    AE_SSX2X2_IP(Amax2,Amax3,castxcc(xtfloatx4,pmax), sizeof(xtfloatx4));
    __Pragma("no_reorder")
    /* multiply k-th row by the reciprocal *
     * pivot element during swapping rows  */
    MUL_SX2X2(Ak0,Ak1,Ak0,Ak1,R,R);
    MUL_SX2X2(Ak2,Ak3,Ak2,Ak3,R,R);
    AE_SSX2X2_IP(Ak0,Ak1,castxcc(xtfloatx4,pk), sizeof(xtfloatx4));
    AE_SSX2X2_IP(Ak2,Ak3,castxcc(xtfloatx4,pk), sizeof(xtfloatx4));

    /* elimination */
    /* join forward and back substitution */

    C0 = XT_LSX((xtfloat *)pA0, k*SZ_F32);
    C1 = XT_LSX((xtfloat *)pA1, k*SZ_F32);
    C2 = XT_LSX((xtfloat *)pA2, k*SZ_F32);

    AE_LSX2X2_I (A00,A01,(const xtfloatx4*)pA0, 0*sizeof(xtfloatx4));
    AE_LSX2X2_I (A02,A03,(const xtfloatx4*)pA0, 1*sizeof(xtfloatx4));
    AE_LSX2X2_I (A10,A11,(const xtfloatx4*)pA1, 0*sizeof(xtfloatx4));
    AE_LSX2X2_I (A12,A13,(const xtfloatx4*)pA1, 1*sizeof(xtfloatx4));
    AE_LSX2X2_I (A20,A21,(const xtfloatx4*)pA2, 0*sizeof(xtfloatx4));
    AE_LSX2X2_I (A22,A23,(const xtfloatx4*)pA2, 1*sizeof(xtfloatx4));

    MSUB_SX2X2(A00,A01, Ak0, Ak1, C0,C0);
    MSUB_SX2X2(A02,A03, Ak2, Ak3, C0,C0);
    MSUB_SX2X2(A10,A11, Ak0, Ak1, C1,C1);
    MSUB_SX2X2(A12,A13, Ak2, Ak3, C1,C1);
    MSUB_SX2X2(A20,A21, Ak0, Ak1, C2,C2);
    MSUB_SX2X2(A22,A23, Ak2, Ak3, C2,C2);

    AE_SSX2X2_IP (A00,A01,castxcc(xtfloatx4,pA0), sizeof(xtfloatx4));
    AE_SSX2X2_XC (A02,A03,castxcc(xtfloatx4,pA0), sizeof(xtfloatx4));
    AE_SSX2X2_IP (A10,A11,castxcc(xtfloatx4,pA1), sizeof(xtfloatx4));
    AE_SSX2X2_XC (A12,A13,castxcc(xtfloatx4,pA1), sizeof(xtfloatx4));
    AE_SSX2X2_IP (A20,A21,castxcc(xtfloatx4,pA2), sizeof(xtfloatx4));
    AE_SSX2X2_XC (A22,A23,castxcc(xtfloatx4,pA2), sizeof(xtfloatx4));
  }
  /* copy inverse matrix to x */
  pA0 = (xtfloatx2 *)(x);
  pA1 = (xtfloatx2 *)(A+4);
  vX = AE_ZALIGN128();
  {
    AE_LSX2X2_IP(A00,A01,castxcc(xtfloatx4,pA1),4*sizeof(xtfloatx2));
    AE_LSX2X2_IP(A02,A03,castxcc(xtfloatx4,pA1),4*sizeof(xtfloatx2));
    AE_LSX2X2_IP(A10,A11,castxcc(xtfloatx4,pA1),4*sizeof(xtfloatx2));
    AE_LSX2X2_IP(A12,A13,castxcc(xtfloatx4,pA1),4*sizeof(xtfloatx2));

    AE_SASX2X2_IP(A00,A01, vX, castxcc(xtfloatx4,pA0));
    AE_SASX2X2_IP(A02,A03, vX, castxcc(xtfloatx4,pA0));
    AE_SASX2X2_IP(A10,A11, vX, castxcc(xtfloatx4,pA0));
    AE_SASX2X2_IP(A12,A13, vX, castxcc(xtfloatx4,pA0));

  }
 AE_SA128POS_FP(vX, pA0);
}/* mtx_inv4x4f() */
size_t mtx_inv4x4f_getScratchSize        () 
{
    return 32*sizeof(float32_t);
}
#elif (HAVE_FPU)
// for scalar FPU
void mtx_inv4x4f(void* pScr,float32_t* x)
{
    xtfloat * restrict pAn1;
    xtfloat * restrict pAn2;
    xtfloat * restrict pAn3;
    ae_int32x2 * restrict pA0;
    ae_int32x2 * restrict pA1;
    ae_valign vX;
    int n,m,k;
    float32_t *A = (float32_t *)pScr; //float32_t ALIGN(8) A[32];
    /* Set bounds of the buffer */
    WUR_AE_CBEGIN0((uintptr_t)(A));
    WUR_AE_CEND0  ((uintptr_t)(A+32));
      /* Copy the matrix to buffer and
       * initialize identity matrix:
       *
       * x00 x01 x02 x03 1.0 0.0 0.0 0.0
       * x10 x11 x12 x13 0.0 1.0 0.0 0.0
       * x20 x21 x22 x23 0.0 0.0 1.0 0.0
       * x30 x31 x32 x33 0.0 0.0 0.0 1.0
       */
    {
        ae_int32x2 C0,C1;
        ae_int32x2 C01, C10;
        ae_int32x2 A00,A01,A02,A03,A10,A11,A12,A13;
        C0 = AE_ZERO32();
        C1 = 0x3f800000;
        pA0 = (ae_int32x2 *)(x);
        pA1 = (ae_int32x2 *)(A);
        C01 = AE_SEL32_HH(C0, C1);
        C10 = AE_SEL32_HH(C1, C0);
        /* Load input matrix */
        vX = AE_LA64_PP(pA0);
        AE_LA32X2_IP(A00, vX, pA0);
        AE_LA32X2_IP(A01, vX, pA0);
        AE_LA32X2_IP(A02, vX, pA0);
        AE_LA32X2_IP(A03, vX, pA0);
        AE_LA32X2_IP(A10, vX, pA0);
        AE_LA32X2_IP(A11, vX, pA0);
        AE_LA32X2_IP(A12, vX, pA0);
        AE_LA32X2_IP(A13, vX, pA0);
        /* Save input and identity matrix */
        /* 1-st row */
        AE_S32X2_IP(A00, pA1, 2*sizeof(float32_t));
        AE_S32X2_IP(A01, pA1, 2*sizeof(float32_t));
        AE_S32X2_IP(C10, pA1, 2*sizeof(float32_t));
        AE_S32X2_IP(C0 , pA1, 2*sizeof(float32_t));
        /* 2-nd row */
        AE_S32X2_IP(A02, pA1, 2*sizeof(float32_t));
        AE_S32X2_IP(A03, pA1, 2*sizeof(float32_t));
        AE_S32X2_IP(C01, pA1, 2*sizeof(float32_t));
        AE_S32X2_IP(C0 , pA1, 2*sizeof(float32_t));
        /* 3-rd row */
        AE_S32X2_IP(A10, pA1, 2*sizeof(float32_t));
        AE_S32X2_IP(A11, pA1, 2*sizeof(float32_t));
        AE_S32X2_IP(C0 , pA1, 2*sizeof(float32_t));
        AE_S32X2_IP(C10, pA1, 2*sizeof(float32_t));
        /* 4-th row */
        AE_S32X2_IP(A12, pA1, 2*sizeof(float32_t));
        AE_S32X2_IP(A13, pA1, 2*sizeof(float32_t));
        AE_S32X2_IP(C0 , pA1, 2*sizeof(float32_t));
        AE_S32X2_IP(C01, pA1, 2*sizeof(float32_t));
    }

    for(k=0; k<4; k++)
    {
        xtfloat * pmax;
        xtfloat * pAk;
        xtfloat amax, max,R;
        unsigned int imax;
        /* pivoting */
        imax=k;
        amax=0.0f;
        /* find absolute max value in the k-th column */
        pAk  = (xtfloat *)(A+k*9);
        for(n=k; n<4; n++)
        {
            xtbool cond;
            xtfloat a, t;
            XT_LSIP(t,pAk,8*sizeof(xtfloat));
            a = XT_ABS_S(t);

            cond = (XT_OLT_S(amax, a));
            XT_MOVT_S(amax, a, cond);
            XT_MOVT_S(max,  t, cond);
            XT_MOVT  (imax, n, cond);
        }

        R = max;
        R = XT_RECIP_S(R);

        /* swap k-th and imax-th rows */
        pmax = (xtfloat *)(A+imax*8);
        pAk  = (xtfloat *)(A+k*8);
        for(m=0; m<8; m++) 
        { 
            xtfloat t,am; 
            t=XT_LSI(pmax,0);
            am=XT_LSI(pAk,0);
            XT_SSIP(am,pmax,sizeof(xtfloat));
            XT_SSIP(XT_MUL_S(R,t),pAk,sizeof(xtfloat));
        }

        /* elimination */
        __Pragma("no_reorder")
        {
            pAk  = (xtfloat *)(A+k*8);
            xtfloat c1,c2,c3;
            pAn1=(xtfloat*)pAk    ; AE_ADDCIRC_XC(castxcc(ae_int64,pAn1),8*sizeof(xtfloat));
            pAn2=(xtfloat*)pAn1   ; AE_ADDCIRC_XC(castxcc(ae_int64,pAn2),8*sizeof(xtfloat));
            pAn3=(xtfloat*)pAn2   ; AE_ADDCIRC_XC(castxcc(ae_int64,pAn3),8*sizeof(xtfloat));
            c1=XT_LSX(pAn1,k*sizeof(xtfloat));
            c2=XT_LSX(pAn2,k*sizeof(xtfloat));
            c3=XT_LSX(pAn3,k*sizeof(xtfloat));
            for (m=0; m<8; m++)
            {
                xtfloat akm,an1,an2,an3;
                XT_LSIP(akm,pAk,sizeof(xtfloat));
                an1=XT_LSI(pAn1,0);
                an2=XT_LSI(pAn2,0);
                an3=XT_LSI(pAn3,0);
                XT_MSUB_S(an1,akm,c1);
                XT_MSUB_S(an2,akm,c2);
                XT_MSUB_S(an3,akm,c3);
                XT_SSIP(an1,pAn1,sizeof(xtfloat));
                XT_SSIP(an2,pAn2,sizeof(xtfloat));
                XT_SSIP(an3,pAn3,sizeof(xtfloat));
            }
        }
    }
    /* copy to x */
  /* copy inverse matrix to x */
    {
        __Pragma("no_reorder")
        pA0 = (ae_int32x2 *)(x);
        pA1 = (ae_int32x2 *)(A+4);
        ae_int32x2 A00,A01,A02,A03,A10,A11,A12,A13;
        vX = AE_ZALIGN64();
        AE_L32X2_IP(A00, pA1,   2*sizeof(float32_t));
        AE_L32X2_IP(A01, pA1, 3*2*sizeof(float32_t));
        AE_L32X2_IP(A02, pA1,   2*sizeof(float32_t));
        AE_L32X2_IP(A03, pA1, 3*2*sizeof(float32_t));
        AE_L32X2_IP(A10, pA1,   2*sizeof(float32_t));
        AE_L32X2_IP(A11, pA1, 3*2*sizeof(float32_t));
        AE_L32X2_IP(A12, pA1,   2*sizeof(float32_t));
        AE_L32X2_IP(A13, pA1,   0);

        AE_SA32X2_IP(A00, vX, pA0);
        AE_SA32X2_IP(A01, vX, pA0);
        AE_SA32X2_IP(A02, vX, pA0);
        AE_SA32X2_IP(A03, vX, pA0);
        AE_SA32X2_IP(A10, vX, pA0);
        AE_SA32X2_IP(A11, vX, pA0);
        AE_SA32X2_IP(A12, vX, pA0);
        AE_SA32X2_IP(A13, vX, pA0);
        AE_SA64POS_FP(vX, pA0);
    }
}
size_t mtx_inv4x4f_getScratchSize        () 
{
    return 32*sizeof(float32_t);
}
#else
DISCARD_FUN(void,mtx_inv4x4f,(void* pScr,float32_t* x))
size_t mtx_inv4x4f_getScratchSize        () 
{
    return 0;
}
#endif
