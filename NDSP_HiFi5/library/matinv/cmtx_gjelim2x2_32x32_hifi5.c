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
 * complex Matrix Gauss-Jordan Elimination for linear equation problem, 32-bit 
   fixed point API
 * Optimized code for HiFi5
 */

/* Cross-platform data type definitions. */
#include "NatureDSP_types.h"
/* Common helper macros. */
#include "common.h"
/* Matrix functions */
#include "NatureDSP_Signal_matinv.h"

/*-------------------------------------------------------------------------
  These functions implement Gauss elimination elimination process with full 
  pivoting to find solution of linear equations A*y=x
  
  Fixed point version takes representation of input matrix/vector and forms 
  representation of output vector with proper scaling.

  Precision: 
  f     floating point
  32x32 32-bit input, 32-bit output

  Input:
  A[N*N]      input matrix, representation is defined by parameter qA
  x[N]        input rigth side of equation. For fixed point API, 
              representation is defined by parameter qX
  qA          input matrix A representation (for fixed point API only)
  qX          input vector x representation (for fixed point API only)
  Output:
  y[N]        output vector
  Temporary:
  pScr        scratch memory. Size in bytes is defined by corresponding 
              scratch allocation function 
  N is 2,3,4,6,8,10

  Returned value: fixed point functions return fixed-point representation 
                  of resulted vector
  Restrictions:
  none
-------------------------------------------------------------------------*/
int  cmtx_gjelim2x2_32x32  (void* pScr, complex_fract32 *y, const complex_fract32* A,const complex_fract32 * x, int qA, int qX)
{
    const int N=2;
    xtbool2 bnegim=AE_LT32(AE_MOVDA32X2(0,-1),AE_ZERO32());
    int k;
    ae_int32x2 *B; // [N][N]
    ae_int32x2 *C; // [N][N]
    ae_int32x2 *T; // [N]
    ae_int32x2 * restrict pB; // [N][N]
    ae_int32x2 * restrict pC; // [N][N]
    ae_int32x2 * restrict pT; // [N][N]
    ae_valignx2 aX;
    const ae_int32x2 * restrict pX;
          ae_int32x2 * restrict pY;
    ae_int32x2 * restrict pBwr; // [N][N]
    ae_int32x2 * restrict pCwr; // [N]

    int qB,qC; // fixed point representations
    NASSERT_ALIGN16(pScr);
    // allocate on scratch
    B=(ae_int32x2 *)pScr;
    C=B+N*N;
    T=C+N;
    // setup circular pointers for B and T
    WUR_AE_CBEGIN0((uintptr_t)B);    WUR_AE_CEND0((uintptr_t)(B+N*N));
    WUR_AE_CBEGIN1((uintptr_t)T);    WUR_AE_CEND1((uintptr_t)(T+N));
    // copy input
    pB=B;
    pX=(const ae_int32x2*)A;
    aX=AE_LA128_PP(pX);
    {
        ae_int32x2 b0,b1,b2,b3;
        AE_LA32X2X2_IP(b0,b1,aX,castxcc(ae_int32x4,pX));
        AE_LA32X2X2_IP(b2,b3,aX,castxcc(ae_int32x4,pX));
        AE_S32X2X2_IP(b0,b1,castxcc(ae_int32x4,pB),sizeof(ae_int32x4));
        AE_S32X2X2_IP(b2,b3,castxcc(ae_int32x4,pB),sizeof(ae_int32x4));
    }
    pC=C;
    pX=(const ae_int32x2*)x;
    aX=AE_LA128_PP(pX);
    {
        ae_int32x2 c0,c1;
        AE_LA32X2X2_IP(c0,c1,aX,castxcc(ae_int32x4,pX));
        AE_S32X2X2_IP(c0,c1,castxcc(ae_int32x4,pC),sizeof(ae_int32x4));
    }

    qB=31;
    qC=qX+(31-qA); // representation of inverted matrix

    for (k=0; k<N; k++)
    {
        int imax;
        int e,expB,expC;
        // find matrix normalization scale
        {
            ae_int32x2 b0,b1,b2,b3,c0,c1;
            AE_L32X2X2_I (b0,b1,(const ae_int32x4*)B,0*sizeof(ae_int32x4));
            AE_L32X2X2_I (b2,b3,(const ae_int32x4*)B,1*sizeof(ae_int32x4));
            AE_L32X2X2_I (c0,c1,(const ae_int32x4*)C,0);
            expB=AE_RMIN16X4(AE_MIN16(AE_NSA32X4(b0,b1),AE_NSA32X4(b2,b3)));
            expC=AE_RMIN16X4(AE_NSA32X4(c0,c1));

            // pivoting
            if (k==0)
            {
                int off0,off1;
                imax=0; 
                XT_MOVT(imax, 1, AE_LE64(AE_MULZAAD32_HH_LL(b0,b0),AE_MULZAAD32_HH_LL(b2,b2)));
                off0=(int)(imax*sizeof(ae_int32x2)*N);
                off1=(int)((N-1-imax)*sizeof(ae_int32x2)*N);
                AE_S32X2X2_X(b0,b1,(ae_int32x4*)B,off0);
                AE_S32X2X2_X(b2,b3,(ae_int32x4*)B,off1);
                pC=&C[k];        
                off0=(int)(imax*sizeof(ae_int32x2));
                off1=(int)((N-1-imax)*sizeof(ae_int32x2));
                AE_S32X2_X(c0,(ae_int32x2*)C,off0);
                AE_S32X2_X(c1,(ae_int32x2*)C,off1);
            }
        }
        // find normalization factor
        {
            ae_int32x2 Tk;
            int e_den,e_num,e_rden;
            ae_int64 d2,n2;
            ae_int32x2 rden;
            pT=T;
            pB=(B+k);
            {
                ae_int32x2 b0,b1;
                AE_L32X2_XP(b0,pB,sizeof(ae_int32x2)*N);
                AE_L32X2_XP(b1,pB,sizeof(ae_int32x2)*N);
                AE_S32X2X2_IP(b0,b1,castxcc(ae_int32x4,pT),2*sizeof(ae_int32x2));
            }
            pT=T;
            Tk=AE_L32X2_X(pT,k*sizeof(ae_int32x2));
            AE_S32X2_X(AE_SLAA32S(AE_MOVDA32X2(1,0),qB),pT,k*sizeof(ae_int32x2));
            d2=AE_MULZAAD32_HH_LL(Tk,Tk);
            n2=AE_ZERO64();
            {
                ae_int32x2 Tn0,Tn1;
                ae_int64 t;
                AE_L32X2X2_IP(Tn0,Tn1,castxcc(ae_int32x4,pT),2*sizeof(ae_int32x2));
                t=AE_MULZAAD32_HH_LL(Tn0,Tn0); n2=t;
                t=AE_MULZAAD32_HH_LL(Tn1,Tn1); n2=AE_MAX64(n2,t);
            }
            n2=AE_ADD64S(AE_SRAI64(n2,1),AE_SRAI64(d2,1));
            e_num=AE_NSA64(n2);
            e_num=(e_num>>1)-1;
            e_den=AE_NSA64(d2);
            e_den=(e_den>>1);
            d2=AE_SLAA64(d2,(e_den<<1));
            rden=AE_TRUNCI32X2F64S(d2,d2,0);
            e_rden=AE_NSAZ32_L(rden);
            rden=AE_SLAA32S(rden,e_rden);
            {   // reciprocal: assume rden is positive!
                ae_int32x2 X,Y,E;
                X=rden;
                Y=AE_SUB32(0xBAEC0000,X);
                E=0x40000000; AE_MULSFP32X2RAS(E,X,Y); AE_MULAFP32X2RAS(Y,Y,AE_SLLI32(E,1));
                E=0x40000000; AE_MULSFP32X2RAS(E,X,Y); AE_MULAFP32X2RAS(Y,Y,AE_SLLI32(E,1));
                E=0x40000000; AE_MULSFP32X2RAS(E,X,Y); AE_MULAFP32X2RAS(Y,Y,AE_SLLI32(E,1));
                E=0x40000000; AE_MULSFP32X2RAS(E,X,Y); AE_MULAFP32X2RAS(Y,Y,AE_SLLI32(E,1));
                rden=Y;
            }
            Tk=AE_MULFP32X2RAS(Tk,rden);
            AE_MOVT32X2(Tk,AE_NEG32S(Tk),bnegim);
            Tk=AE_SLAA32S(Tk,e_den+e_rden-1);
            pT=T;
            {
                ae_int32x2 Tn0,Tn1;
                AE_L32X2X2_I(Tn0,Tn1,(const ae_int32x4*)pT,0);
                Tn0=AE_SLAA32S(Tn0,e_num);
                Tn1=AE_SLAA32S(Tn1,e_num);
                Tn0=AE_MULFC32RAS(Tn0,Tk);
                Tn1=AE_MULFC32RAS(Tn1,Tk);
                AE_S32X2X2_IP(Tn0,Tn1,castxcc(ae_int32x4,pT),2*sizeof(ae_int32x2));
            }
           e=e_den-e_num+1;
        }
        // scale B and C (2-way shifts possible!)
        qB=qB+expB-e;
        qC=qC+expC-e;
        // Gauss-Jordan elimination: might be made using curcular addressing on B/C
        {
            ae_int32x2 Cin;
            ae_int32x2 Ti,Ckn;
            ae_int32x2 Bk0,Bk1;
            pB=(ae_int32x2 *)B;
            pC=(ae_int32x2 *)C;
            pBwr=(ae_int32x2 *)B;
            pCwr=(ae_int32x2 *)C;
            {
                ae_int32x2 expB_coef,expC_coef;
                ae_int32x2 b0,b1,c0,c1,b2,b3;
                expB_coef=AE_SLAA32S(1,expB);
                expC_coef=AE_SLAA32S(1,expC);
                AE_L32X2X2_I(b0,b1,(const ae_int32x4*)pB,0*sizeof(ae_int32x4));
                AE_L32X2X2_I(b2,b3,(const ae_int32x4*)pB,1*sizeof(ae_int32x4));
                AE_L32X2X2_I(c0,c1,(const ae_int32x4*)pC,0*sizeof(ae_int32x4));
                AE_MUL2P32X4(b0,b1,b0,b1,expB_coef,expB_coef);
                AE_MUL2P32X4(b2,b3,b2,b3,expB_coef,expB_coef);
                AE_MUL2P32X4(c0,c1,c0,c1,expC_coef,expC_coef);
                AE_S32X2X2_I(b0,b1, (ae_int32x4*)pBwr,0*sizeof(ae_int32x4));
                AE_S32X2X2_I(b2,b3, (ae_int32x4*)pBwr,1*sizeof(ae_int32x4));
                AE_S32X2X2_I(c0,c1, (ae_int32x4*)pCwr,0*sizeof(ae_int32x4));
                __Pragma("no_reorder")
            }
            pC=(ae_int32x2 *)C;
            Ckn=pC[k];
            pB=(B+k*N);
            pT=(T+k);
            AE_L32X2X2_XC(Bk0,Bk1,castxcc(ae_int32x4,pB),sizeof(ae_int32x2)*N);
            pBwr=pB;
            AE_ADDCIRC32X2_XC1(pT,sizeof(ae_int32x2));
            {
                ae_int32x2 Bin0,Bin1;
                Ti=AE_L32X2_I(pT,0);
                Ti=AE_NEG32S(Ti);
                Cin=AE_L32X2_X (pT,(int)((uintptr_t)C-(uintptr_t)T));
                Cin=AE_SLAA32S(Cin,-e);
                AE_MULAFC32RAS(Cin,Ckn,Ti);
                AE_S32X2_X (Cin,pT,(int)((uintptr_t)C-(uintptr_t)T));
                AE_ADDCIRC32X2_XC1(pT,sizeof(ae_int32x2));
                AE_L32X2X2_I(Bin0,Bin1,(const ae_int32x4*)pB,0);
                Bin0=AE_SLAA32S(Bin0,-e);
                Bin1=AE_SLAA32S(Bin1,-e);
                AE_MULAFC32RAS(Bin0,Bk0,Ti);
                AE_MULAFC32RAS(Bin1,Bk1,Ti);
                AE_S32X2X2_XC (Bin0,Bin1,castxcc(ae_int32x4,pBwr),sizeof(ae_int32x2)*N);
            }
            Ti=AE_L32X2_I(pT,0);
            Cin=AE_MULFC32RAS(Ckn,Ti);
            AE_S32X2_X (Cin,pT,(int)((uintptr_t)C-(uintptr_t)T));
            AE_ADDCIRC32X2_XC1(pT,sizeof(ae_int32x2));
            Bk0=AE_MULFC32RAS(Bk0,Ti);
            Bk1=AE_MULFC32RAS(Bk1,Ti);
            AE_S32X2X2_XC (Bk0,Bk1,castxcc(ae_int32x4,pBwr),sizeof(ae_int32x2)*N);
            pBwr=(B+k);
            AE_S32X2_XP(AE_ZERO32(),pBwr,N*sizeof(ae_int32x2));
            AE_S32X2_XP(AE_ZERO32(),pBwr,N*sizeof(ae_int32x2));
        }
    }
    // copy back to the output
    pC=C;
    pY=(ae_int32x2 *)y;
    {
        ae_int32x2 c0,c1;
        aX=AE_ZALIGN128();
        AE_L32X2X2_IP(c0,c1,castxcc(ae_int32x4,pC),2*sizeof(ae_int32x2));
        AE_SA32X2X2_IP(c0,c1,aX,castxcc(ae_int32x4,pY));
        AE_SA128POS_FP(aX,pY);
    }
    return qC;
}
// scratch allocation
size_t cmtx_gjelim2x2_32x32_getScratchSize   () 
{
    return  (2*2+2*2)*sizeof(complex_fract32);  
}
