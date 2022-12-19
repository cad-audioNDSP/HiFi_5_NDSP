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
 * complex Matrix Gauss-Jordan Elimination for linear equation problem 8x8, 32-bit 
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
int  cmtx_gjelim8x8_32x32  (void* pScr, complex_fract32 *y, const complex_fract32* A,const complex_fract32 * x, int qA, int qX)
{
    const int N=8;
    xtbool2 bnegim=AE_LT32(AE_MOVDA32X2(0,-1),AE_ZERO32());
    int k,n;
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
    for (k=0; k<N*N/2; k++) 
    {
        ae_int32x2 b0,b1;
        AE_LA32X2X2_IP(b0,b1,aX,castxcc(ae_int32x4,pX));
        AE_S32X2X2_IP(b0,b1,castxcc(ae_int32x4,pB),sizeof(ae_int32x4));
    }
    pC=C;
    pX=(const ae_int32x2*)x;
    aX=AE_LA128_PP(pX);
    for (k=0; k<N/2; k++) 
    {
        ae_int32x2 b0,b1;
        AE_LA32X2X2_IP(b0,b1,aX,castxcc(ae_int32x4,pX));
        AE_S32X2X2_IP (b0,b1,castxcc(ae_int32x4,pC),sizeof(ae_int32x4));
    }

    qB=31;
    qC=qX+(31-qA); // representation of inverted matrix

    for (k=0; k<N; k++)
    {
        int i,imax;
        int e,expB,expC;
        // find matrix normalization scale
        pB=(ae_int32x2 *)B;
        pC=(ae_int32x2 *)C;
        {
            ae_int16x4 minnsaB0=31,minnsaB1=31,minnsaC=31;
            for(n=0; n<(N); n++) 
            {
                ae_int32x2 b0,b1,b2,b3,b4,b5,b6,b7,c0;
                AE_L32X2X2_IP(b0,b1,castxcc(ae_int32x4,pB),2*sizeof(ae_int32x2));
                AE_L32X2X2_IP(b2,b3,castxcc(ae_int32x4,pB),2*sizeof(ae_int32x2));
                AE_L32X2X2_IP(b4,b5,castxcc(ae_int32x4,pB),2*sizeof(ae_int32x2));
                AE_L32X2X2_IP(b6,b7,castxcc(ae_int32x4,pB),2*sizeof(ae_int32x2));
                AE_L32X2_IP  (c0,   castxcc(ae_int32x2,pC),1*sizeof(ae_int32x2));
                minnsaB0=AE_MIN16(minnsaB0,AE_NSA32X4(b0,b1));
                minnsaB1=AE_MIN16(minnsaB1,AE_NSA32X4(b2,b3));
                minnsaB0=AE_MIN16(minnsaB0,AE_NSA32X4(b4,b5));
                minnsaB1=AE_MIN16(minnsaB1,AE_NSA32X4(b6,b7));
                minnsaC=AE_MIN16(minnsaC,AE_NSA32X4(c0,c0));
            }
            expB=AE_RMIN16X4(AE_MIN16(minnsaB0,minnsaB1));
            expC=AE_RMIN16X4(minnsaC);
        }
        // pivoting
        {
            ae_int64 bmax64;
            imax=k; bmax64=AE_ZERO64();
            pB=&B[k*N+k];
            for (n=k; n<N; n++)
            {
                ae_int32x2 bb;
                ae_int64 b;
                xtbool bbig;
                AE_L32X2_XP(bb,pB,sizeof(ae_int32x2)*N);
                b=AE_MULZAAD32_HH_LL(bb,bb);
                bbig=AE_LE64(bmax64,b);
                AE_MOVT64(bmax64,b,bbig);
                XT_MOVT(imax, n, bbig);
            }
        }
        int off=(int)((imax-k)*sizeof(ae_int32x2)*N);
        pB=&B[k*N];
        pC=&C[k];        
        for (n=0; n<N; n++) 
        {
            ae_int32x2 bk,bi;
            bk=AE_L32X2_I(pB,0);
            bi=AE_L32X2_X(pB,off);
            AE_S32X2_X (bk,pB,off);
            AE_S32X2_IP(bi,pB,sizeof(ae_int32x2));
        }
        off=(int)((imax-k)*sizeof(ae_int32x2));
        {
            ae_int32x2 ck,ci;
            ck=AE_L32X2_I(pC,0);
            ci=AE_L32X2_X(pC,off);
            AE_S32X2_X (ck,pC,off);
            AE_S32X2_IP(ci,pC,sizeof(ae_int32x2));
        }
        // find normalization factor
        {
            ae_int32x2 Tk;
            int e_den,e_num,e_rden;
            ae_int64 d2,n2;
            ae_int32x2 rden;
            pT=T;
            pB=(B+k);
            for (n=0; n<N/2; n++)
            {
                ae_int32x2 b0,b1;
                AE_L32X2_XP(b0,pB,sizeof(ae_int32x2)*N);
                AE_L32X2_XP(b1,pB,sizeof(ae_int32x2)*N);
                AE_S32X2X2_IP(b0,b1,castxcc(ae_int32x4,pT),sizeof(ae_int32x4));
            }
            pT=T;
            Tk=AE_L32X2_X(pT,k*sizeof(ae_int32x2));
            AE_S32X2_X(AE_SLAA32S(AE_MOVDA32X2(1,0),qB),pT,k*sizeof(ae_int32x2));
            d2=AE_MULZAAD32_HH_LL(Tk,Tk);
            n2=AE_ZERO64();
            for (n=0; n<N/2; n++) 
            {
                ae_int32x2 Tn0,Tn1;
                AE_L32X2X2_IP(Tn0,Tn1,castxcc(ae_int32x4,pT),sizeof(ae_int32x4));
                n2=AE_MAX64(n2,AE_MULZAAD32_HH_LL(Tn0,Tn0));
                n2=AE_MAX64(n2,AE_MULZAAD32_HH_LL(Tn1,Tn1));
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
            for (n=0; n<N/2; n++)
            {
                ae_int32x2 Tn0,Tn1;
                AE_L32X2X2_I(Tn0,Tn1,(const ae_int32x4*)pT,0*sizeof(ae_int32x4));
                Tn0=AE_MULFC32RAS(AE_SLAA32S(Tn0,e_num),Tk);
                Tn1=AE_MULFC32RAS(AE_SLAA32S(Tn1,e_num),Tk);
                AE_S32X2X2_IP(Tn0,Tn1,castxcc(ae_int32x4,pT),sizeof(ae_int32x4));
            }
           e=e_den-e_num+1;
        }
        // scale B and C (2-way shifts possible!)
        qB=qB+expB-e;
        qC=qC+expC-e;
        // Gauss-Jordan elimination:
        if (k<4)
        {
            ae_int32x2 esh_coef,expB_coef,expC_coef;
            ae_int32x2 Cin;
            ae_int32x2 Ti,Ckn;
            ae_int32x2 Bk0,Bk1,Bk2,Bk3,Bk4,Bk5,Bk6,Bk7;
            pB=(ae_int32x2 *)B;
            pC=(ae_int32x2 *)C;
            pBwr=(ae_int32x2 *)B;
            pCwr=(ae_int32x2 *)C;
            expB_coef=AE_SLAA32S(1,expB);
            expC_coef=AE_SLAA32S(1,expC);
            for (n=0; n<N; n++)
            {
                ae_int32x2 b0,b1,b2,b3,b4,b5,b6,b7,c0;
                AE_L32X2X2_IP(b0,b1,castxcc(ae_int32x4,pB),sizeof(ae_int32x4));
                AE_L32X2X2_IP(b2,b3,castxcc(ae_int32x4,pB),sizeof(ae_int32x4));
                AE_L32X2X2_IP(b4,b5,castxcc(ae_int32x4,pB),sizeof(ae_int32x4));
                AE_L32X2X2_IP(b6,b7,castxcc(ae_int32x4,pB),sizeof(ae_int32x4));
                AE_L32X2_IP(c0,castxcc(ae_int32x2,pC),sizeof(ae_int32x2));
                AE_MUL2P32X4(b0,b1,b0,b1,expB_coef,expB_coef);
                AE_MUL2P32X4(b2,b3,b2,b3,expB_coef,expB_coef);
                AE_MUL2P32X4(b4,b5,b4,b5,expB_coef,expB_coef);
                AE_MUL2P32X4(b6,b7,b6,b7,expB_coef,expB_coef);
                c0=AE_MULP32X2(c0,expC_coef);
                AE_S32X2X2_IP(b0,b1,castxcc(ae_int32x4,pBwr),sizeof(ae_int32x4));
                AE_S32X2X2_IP(b2,b3,castxcc(ae_int32x4,pBwr),sizeof(ae_int32x4));
                AE_S32X2X2_IP(b4,b5,castxcc(ae_int32x4,pBwr),sizeof(ae_int32x4));
                AE_S32X2X2_IP(b6,b7,castxcc(ae_int32x4,pBwr),sizeof(ae_int32x4));
                AE_S32X2_IP(c0,castxcc(ae_int32x2,pCwr),sizeof(ae_int32x2));
            }
            __Pragma("no_reorder")
            pC=(ae_int32x2 *)C;
            Ckn=pC[k];
            pB=(B+k*N);
            pT=(T+k);
            AE_L32X2X2_I (Bk2,Bk3,(const ae_int32x4*)pB,1*sizeof(ae_int32x4));
            AE_L32X2X2_I (Bk4,Bk5,(const ae_int32x4*)pB,2*sizeof(ae_int32x4));
            AE_L32X2X2_I (Bk6,Bk7,(const ae_int32x4*)pB,3*sizeof(ae_int32x4));
            AE_L32X2X2_XC(Bk0,Bk1,castxcc(ae_int32x4,pB),N*sizeof(ae_int32x2));
            AE_ADDCIRC32X2_XC1(pT,sizeof(ae_int32x2));
            pBwr=pB;
            esh_coef=AE_SLAA32S(0x40000000,-e+1);
            for (i=0; i<N-1; i++)
            {
                ae_int32x2 Bin0,Bin1,Bin2,Bin3,Bin4,Bin5,Bin6,Bin7;
                Ti=AE_L32X2_I(pT,0);
                Ti=AE_NEG32S(Ti);
                Cin=AE_L32X2_X (pT,(int)((uintptr_t)C-(uintptr_t)T));
                Cin=AE_SLAA32S(Cin,-e);
                AE_MULAFC32RAS(Cin,Ckn,Ti);
                AE_S32X2_X (Cin,pT,(int)((uintptr_t)C-(uintptr_t)T));
                AE_ADDCIRC32X2_XC1(pT,sizeof(ae_int32x2));
                AE_L32X2X2_I (Bin2,Bin3,(const ae_int32x4*)pB,1*sizeof(ae_int32x4));
                AE_L32X2X2_I (Bin4,Bin5,(const ae_int32x4*)pB,2*sizeof(ae_int32x4));
                AE_L32X2X2_I (Bin6,Bin7,(const ae_int32x4*)pB,3*sizeof(ae_int32x4));
                AE_L32X2X2_XC(Bin0,Bin1,castxcc(ae_int32x4,pB),N*sizeof(ae_int32x2));
                AE_MULF2P32X4RAS(Bin0,Bin1,Bin0,Bin1,esh_coef,esh_coef);
                AE_MULF2P32X4RAS(Bin2,Bin3,Bin2,Bin3,esh_coef,esh_coef);
                AE_MULF2P32X4RAS(Bin4,Bin5,Bin4,Bin5,esh_coef,esh_coef);
                AE_MULF2P32X4RAS(Bin6,Bin7,Bin6,Bin7,esh_coef,esh_coef);
                AE_MULAFC32RAS(Bin0,Bk0,Ti);
                AE_MULAFC32RAS(Bin1,Bk1,Ti);
                AE_MULAFC32RAS(Bin2,Bk2,Ti);
                AE_MULAFC32RAS(Bin3,Bk3,Ti);
                AE_MULAFC32RAS(Bin4,Bk4,Ti);
                AE_MULAFC32RAS(Bin5,Bk5,Ti);
                AE_MULAFC32RAS(Bin6,Bk6,Ti);
                AE_MULAFC32RAS(Bin7,Bk7,Ti);
                AE_S32X2X2_I (Bin2,Bin3,(ae_int32x4*)pBwr,1*sizeof(ae_int32x4));
                AE_S32X2X2_I (Bin4,Bin5,(ae_int32x4*)pBwr,2*sizeof(ae_int32x4));
                AE_S32X2X2_I (Bin6,Bin7,(ae_int32x4*)pBwr,3*sizeof(ae_int32x4));
                AE_S32X2X2_XC(Bin0,Bin1,castxcc(ae_int32x4,pBwr),N*sizeof(ae_int32x2));
            }
            Ti=AE_L32X2_I(pT,0);
            Ckn=AE_MULFC32RAS(Ckn,Ti);
            AE_S32X2_X (Ckn,pT,(int)((uintptr_t)C-(uintptr_t)T));
            Bk0=AE_MULFC32RAS(Bk0,Ti);
            Bk1=AE_MULFC32RAS(Bk1,Ti);
            Bk2=AE_MULFC32RAS(Bk2,Ti);
            Bk3=AE_MULFC32RAS(Bk3,Ti);
            Bk4=AE_MULFC32RAS(Bk4,Ti);
            Bk5=AE_MULFC32RAS(Bk5,Ti);
            Bk6=AE_MULFC32RAS(Bk6,Ti);
            Bk7=AE_MULFC32RAS(Bk7,Ti);
            AE_S32X2X2_I (Bk2,Bk3,(ae_int32x4*)pBwr,1*sizeof(ae_int32x4));
            AE_S32X2X2_I (Bk4,Bk5,(ae_int32x4*)pBwr,2*sizeof(ae_int32x4));
            AE_S32X2X2_I (Bk6,Bk7,(ae_int32x4*)pBwr,3*sizeof(ae_int32x4));
            AE_S32X2X2_XC(Bk0,Bk1,castxcc(ae_int32x4,pBwr),N*sizeof(ae_int32x2));
            pBwr=(B+k);
            for (i=0; i<N; i++)  AE_S32X2_XP(AE_ZERO32(),pBwr,N*sizeof(ae_int32x2));
            __Pragma("no_reorder");
        }
        else    // k=4...7
        {
            ae_int32x2 esh_coef,expB_coef,expC_coef;
            ae_int32x2 Cin;
            ae_int32x2 Ti,Ckn;
            ae_int32x2 Bk4,Bk5,Bk6,Bk7;
            pB=((ae_int32x2 *)B)+4;
            pC=(ae_int32x2 *)C;
            pBwr=((ae_int32x2 *)B)+4;
            pCwr=(ae_int32x2 *)C;
            expB_coef=AE_SLAA32S(1,expB);
            expC_coef=AE_SLAA32S(1,expC);
            for (n=0; n<N; n++)
            {
                ae_int32x2 b4,b5,b6,b7,c0;
                AE_L32X2X2_IP(b4,b5,castxcc(ae_int32x4,pB),  sizeof(ae_int32x4));
                AE_L32X2X2_IP(b6,b7,castxcc(ae_int32x4,pB),3*sizeof(ae_int32x4));
                AE_L32X2_IP(c0,castxcc(ae_int32x2,pC),sizeof(ae_int32x2));
                AE_MUL2P32X4(b4,b5,b4,b5,expB_coef,expB_coef);
                AE_MUL2P32X4(b6,b7,b6,b7,expB_coef,expB_coef);
                c0=AE_MULP32X2(c0,expC_coef);
                AE_S32X2X2_IP(b4,b5,castxcc(ae_int32x4,pBwr),  sizeof(ae_int32x4));
                AE_S32X2X2_IP(b6,b7,castxcc(ae_int32x4,pBwr),3*sizeof(ae_int32x4));
                AE_S32X2_IP(c0,castxcc(ae_int32x2,pCwr),sizeof(ae_int32x2));
            }
            __Pragma("no_reorder");
            pC=(ae_int32x2 *)C;
            Ckn=pC[k];
            pB=(B+k*N+4);
            pT=(T+k);
            AE_L32X2X2_I (Bk6,Bk7,(const ae_int32x4*)pB,1*sizeof(ae_int32x4));
            AE_L32X2X2_XC(Bk4,Bk5,castxcc(ae_int32x4,pB),N*sizeof(ae_int32x2));
            AE_ADDCIRC32X2_XC1(pT,sizeof(ae_int32x2));
            pBwr=pB;
            esh_coef=AE_SLAA32S(0x40000000,-e+1);
            for (i=0; i<N-1; i++)
            {
                ae_int32x2 Bin4,Bin5,Bin6,Bin7;
                Ti=AE_L32X2_I(pT,0);
                Ti=AE_NEG32S(Ti);
                Cin=AE_L32X2_X (pT,(int)((uintptr_t)C-(uintptr_t)T));
                Cin=AE_SLAA32S(Cin,-e);
                AE_MULAFC32RAS(Cin,Ckn,Ti);
                AE_S32X2_X (Cin,pT,(int)((uintptr_t)C-(uintptr_t)T));
                AE_ADDCIRC32X2_XC1(pT,sizeof(ae_int32x2));
                AE_L32X2X2_I (Bin6,Bin7,(const ae_int32x4*)pB,1*sizeof(ae_int32x4));
                AE_L32X2X2_XC(Bin4,Bin5,castxcc(ae_int32x4,pB),N*sizeof(ae_int32x2));
                AE_MULF2P32X4RAS(Bin4,Bin5,Bin4,Bin5,esh_coef,esh_coef);
                AE_MULF2P32X4RAS(Bin6,Bin7,Bin6,Bin7,esh_coef,esh_coef);
                AE_MULAFC32RAS(Bin4,Bk4,Ti);
                AE_MULAFC32RAS(Bin5,Bk5,Ti);
                AE_MULAFC32RAS(Bin6,Bk6,Ti);
                AE_MULAFC32RAS(Bin7,Bk7,Ti);
                AE_S32X2X2_I (Bin6,Bin7,(ae_int32x4*)pBwr,1*sizeof(ae_int32x4));
                AE_S32X2X2_XC(Bin4,Bin5,castxcc(ae_int32x4,pBwr),N*sizeof(ae_int32x2));
            }
            Ti=AE_L32X2_I(pT,0);
            Ckn=AE_MULFC32RAS(Ckn,Ti);
            AE_S32X2_X (Ckn,pT,(int)((uintptr_t)C-(uintptr_t)T));
            Bk4=AE_MULFC32RAS(Bk4,Ti);
            Bk5=AE_MULFC32RAS(Bk5,Ti);
            Bk6=AE_MULFC32RAS(Bk6,Ti);
            Bk7=AE_MULFC32RAS(Bk7,Ti);
            AE_S32X2X2_I (Bk6,Bk7,(ae_int32x4*)pBwr,1*sizeof(ae_int32x4));
            AE_S32X2X2_XC(Bk4,Bk5,castxcc(ae_int32x4,pBwr),N*sizeof(ae_int32x2));
            pBwr=(B+k);
            for (i=0; i<N; i++)  AE_S32X2_XP(AE_ZERO32(),pBwr,N*sizeof(ae_int32x2));
            __Pragma("no_reorder");
        }
    }
    // copy back to the output
    pC=C;
    pY=(ae_int32x2 *)y;
    aX=AE_ZALIGN128();
    {
        ae_int32x2 c0,c1,c2,c3;
        AE_L32X2X2_IP(c0,c1,castxcc(ae_int32x4,pC),sizeof(ae_int32x4));
        AE_L32X2X2_IP(c2,c3,castxcc(ae_int32x4,pC),sizeof(ae_int32x4));
        AE_SA32X2X2_IP(c0,c1,aX,castxcc(ae_int32x4,pY));
        AE_SA32X2X2_IP(c2,c3,aX,castxcc(ae_int32x4,pY));
        AE_L32X2X2_IP(c0,c1,castxcc(ae_int32x4,pC),sizeof(ae_int32x4));
        AE_L32X2X2_IP(c2,c3,castxcc(ae_int32x4,pC),sizeof(ae_int32x4));
        AE_SA32X2X2_IP(c0,c1,aX,castxcc(ae_int32x4,pY));
        AE_SA32X2X2_IP(c2,c3,aX,castxcc(ae_int32x4,pY));
        AE_SA128POS_FP(aX,pY);
    }
    return qC;
}
// scratch allocation
size_t cmtx_gjelim8x8_32x32_getScratchSize   () 
{ 
    return  (8*8+2*8)*sizeof(complex_fract32);  
}
