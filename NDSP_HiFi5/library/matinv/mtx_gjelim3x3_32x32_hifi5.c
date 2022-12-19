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
 * Real Matrix Gauss-Jordan Elimination for linear equation problem 3x3, 32-bit 
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
int  mtx_gjelim3x3_32x32  (void* pScr, int32_t *y, const int32_t* A,const int32_t * x, int qA, int qX) 
{
    const int N=3;
    int N0=(N+1)&~1;
    int k,n;
    int32_t *B; // [N][N]
    int32_t *C; // [N0]
    int32_t *T; // [N0]
    ae_int32x2 * restrict pB; // [N][N]
    ae_int32x2 * restrict pC; // [N]
    ae_int32x2 * restrict pT; // [N]
    const ae_int32x2 * restrict pX; // [N]
          ae_int32x2 * restrict pY; // [N]
    ae_int32x2 * restrict pBwr; // [N][N]
    ae_int32x2 * restrict pCwr; // [N]

    int qB,qC; // fixed point representations
    NASSERT_ALIGN16(pScr);
    // allocate on scratch
    B=(int32_t *)pScr;
    C=B+N*N0;
    T=C+N0;
    NASSERT_ALIGN16(B);
    NASSERT_ALIGN16(C);
    NASSERT_ALIGN16(T);
    // setup circular pointers for B and T
    WUR_AE_CBEGIN2((uintptr_t)B);    WUR_AE_CEND2((uintptr_t)(B+N*N0));
    WUR_AE_CBEGIN1((uintptr_t)T);    WUR_AE_CEND1((uintptr_t)(T+N));
    WUR_AE_CBEGIN0((uintptr_t)C);    WUR_AE_CEND0((uintptr_t)(C+N));

    // copy input
    {
        ae_valign aX;
        ae_valignx2 ax;
        ae_int32x2  c0,c1;
        pB=(ae_int32x2*)B;
        pX=(ae_int32x2*)A;
        ae_int32x2 b0,b1,b2,b3,b4;
        ax=AE_LA128_PP(pX);
        AE_LA32X2X2_IP(b0,b1,ax,castxcc(ae_int32x4,pX));
        AE_LA32X2X2_IP(b2,b3,ax,castxcc(ae_int32x4,pX));
        b4=AE_L32_I((const ae_int32*)pX,0);
        AE_S32X2X2_I (b0,AE_SEL32_HH(b1,0),(ae_int32x4*)pB,0*sizeof(ae_int32x4));
        AE_S32X2X2_I (AE_SEL32_LH(b1,b2),AE_SEL32_LL(b2,0),(ae_int32x4*)pB,1*sizeof(ae_int32x4));
        AE_S32X2X2_I (b3,AE_SEL32_HH(b4,0),(ae_int32x4*)pB,2*sizeof(ae_int32x4));
        pX=(ae_int32x2*)x;
        pC=(ae_int32x2*)C;
        c1=AE_L32_I((ae_int32*)pX,2*sizeof(int32_t));
        aX=AE_LA64_PP(pX);
        AE_LA32X2_IP(c0,aX,pX);
        AE_S32X2_I(c0,pC,0);
        AE_S32X2_I(AE_SEL32_HL(c1,AE_ZERO32()),pC,1*sizeof(ae_int32x2));
    }
    qB=31;
    qC=qX+(31-qA); // representation of inverted matrix

    for (k=0; k<N; k++)
    {
        int imax;
        int e,expB,expC;
        // find matrix normalization scale
        ae_int16x4 minnsaB,minnsaC;
        pB=(ae_int32x2*)B;
        pC=(ae_int32x2*)C;
        {
            ae_int32x2 b0,b1,b2,b3,b4,b5,c0,c1;
            AE_L32X2X2_IP(b0,b1,castxcc(ae_int32x4,pB),sizeof(ae_int32x4));
            AE_L32X2X2_IP(b2,b3,castxcc(ae_int32x4,pB),sizeof(ae_int32x4));
            AE_L32X2X2_IP(b4,b5,castxcc(ae_int32x4,pB),sizeof(ae_int32x4));
            AE_L32X2X2_IP(c0,c1,castxcc(ae_int32x4,pC),sizeof(ae_int32x4));
            minnsaB=AE_NSA32X4(b0,b1);
            minnsaB=AE_MIN16(minnsaB,AE_NSA32X4(b2,b3));
            minnsaB=AE_MIN16(minnsaB,AE_NSA32X4(b4,b5));
            minnsaC=AE_NSA32X4(c0,c1);
        }
        expB=AE_RMIN16X4(minnsaB);
        expC=AE_RMIN16X4(minnsaC);
        // pivoting
        {
            int off;
            ae_int32x2 bk,bi,ck,ci;
            ae_int32x2 imax2,n2;
            ae_int32x2 bmax=AE_ZERO32();
            pB=(ae_int32x2*)(&B[k*N0+k]);
            imax2=n2=k;
            for (n=k; n<N; n++)
            {
                ae_int32x2 b;
                xtbool2 bbig;
                AE_L32_XP(b,castxcc(ae_int32,pB),N0*sizeof(int32_t));
                b=AE_ABS32S(b);
                bbig=AE_LE32(bmax,b);
                AE_MOVT32X2(imax2,n2,bbig);
                bmax=AE_MAX32(bmax,b);
                n2=AE_ADD32(n2,1);
            }
            imax=AE_MOVAD32_L(imax2);

            off=(int)((imax-k)*sizeof(int32_t)*N0);
            pB=(ae_int32x2*)&B[k*N0];
            bk=AE_L32X2_I(pB,0);
            bi=AE_L32X2_X(pB,off);
            AE_S32X2_X (bk,pB,off);
            AE_S32X2_IP(bi,pB,sizeof(ae_int32x2));
            bk=AE_L32X2_I(pB,0);
            bi=AE_L32X2_X(pB,off);
            AE_S32X2_X (bk,pB,off);
            AE_S32X2_IP(bi,pB,sizeof(ae_int32x2));
            pC=(ae_int32x2*)&C[k];        
            off=(int)((imax-k)*sizeof(int32_t));
            ck=AE_L32_I((const ae_int32*)pC,0);
            ci=AE_L32_X((const ae_int32*)pC,off);
            AE_S32_L_X (ck,(ae_int32*)pC,off);
            AE_S32_L_I (ci,(ae_int32*)pC,0);
        }

        // find normalization factor
        {
            ae_int32x2 rden,tk,maxden,maxnum;
            ae_int32x2 b;
            int e_den,e_num;
            pT=(ae_int32x2*)T;
            pB=(ae_int32x2*)&B[k];
            AE_L32_XP(b,castxcc(ae_int32,pB),N0*sizeof(int32_t));
            AE_S32_L_I(b,(ae_int32*)pT,0*sizeof(int32_t));
            AE_L32_XP(b,castxcc(ae_int32,pB),N0*sizeof(int32_t));
            AE_S32_L_I(b,(ae_int32*)pT,1*sizeof(int32_t));
            AE_L32_XP(b,castxcc(ae_int32,pB),N0*sizeof(int32_t));
            AE_S32_L_I(b,(ae_int32*)pT,2*sizeof(int32_t));
            AE_S32_L_I(AE_ZERO32(),(ae_int32*)pT,3*sizeof(int32_t));
            tk=AE_L32_X((const ae_int32*)pT,k*sizeof(int32_t));
            AE_S32_L_X(AE_SLAA32S(1,qB),(ae_int32*)pT,k*sizeof(int32_t));
            maxnum=AE_MAXABS32S(pT[0],pT[1]);
            maxnum=AE_MAX32(maxnum,AE_SEL32_LH(maxnum,maxnum));
            maxden=AE_ABS32S(tk);
            maxnum=AE_ADD32S(AE_SRAI32(maxnum,1),AE_SRAI32(maxden,1));
            e_den=AE_NSAZ32_L(maxden); 
            e_num=AE_NSAZ32_L(maxnum)-1; 
            e=e_den-e_num+1;
            tk=AE_SLAA32S(tk,e_den);
            {   // reciprocal
                ae_int32x2 x,y,e;
                xtbool2 sx;
                x=tk;
                sx=AE_LT32(x,AE_ZERO32());
                x=AE_ABS32S(x);
                y=AE_SUB32(0xBAEC0000,x);
                e=0x40000000; AE_MULSFP32X2RAS(e,x,y);  AE_MULAFP32X2RAS(y,y,AE_SLAI32(e,1));
                e=0x40000000; AE_MULSFP32X2RAS(e,x,y);  AE_MULAFP32X2RAS(y,y,AE_SLAI32(e,1));
                e=0x40000000; AE_MULSFP32X2RAS(e,x,y);  AE_MULAFP32X2RAS(y,y,AE_SLAI32(e,1));
                e=0x40000000; AE_MULSFP32X2RAS(e,x,y);  AE_MULAFP32X2RAS(y,y,AE_SLAI32(e,1));
                AE_MOVT32X2(y,AE_NEG32(y),sx);
                rden=y;
            }
            pT[0]=AE_MULFP32X2RAS( AE_SLAA32S(pT[0],e_num),rden);
            pT[1]=AE_MULFP32X2RAS( AE_SLAA32S(pT[1],e_num),rden);
        }
        // scale B and C (2-way shifts possible!)
        qB=qB+expB-e;
        qC=qC+expC-e;
        // Gauss-Jordan elimination

        {
            ae_int32x2 esh_coef,expB_coef,expC_coef;
            ae_int32x2 Ti,Ckn;
            ae_int32x2 Bk0,Bk1;
            pB=(ae_int32x2 *)B;
            pC=(ae_int32x2 *)C;
            pBwr=(ae_int32x2 *)B;
            pCwr=(ae_int32x2 *)C;
            expB_coef=AE_SLAA32S(1,expB);
            expC_coef=AE_SLAA32S(1,expC);
            {
                ae_int32x2 b0,b1,b2,b3,b4,b5,c0,c1;
                AE_L32X2X2_IP(b0,b1,castxcc(ae_int32x4,pB),sizeof(ae_int32x4));
                AE_L32X2X2_IP(b2,b3,castxcc(ae_int32x4,pB),sizeof(ae_int32x4));
                AE_L32X2X2_IP(b4,b5,castxcc(ae_int32x4,pB),sizeof(ae_int32x4));
                AE_L32X2X2_IP(c0,c1,castxcc(ae_int32x4,pC),sizeof(ae_int32x4));
                AE_MUL2P32X4(b0,b1,b0,b1,expB_coef,expB_coef);
                AE_MUL2P32X4(b2,b3,b2,b3,expB_coef,expB_coef);
                AE_MUL2P32X4(b4,b5,b4,b5,expB_coef,expB_coef);
                AE_MUL2P32X4(c0,c1,c0,c1,expC_coef,expC_coef);
                AE_S32X2X2_IP(b0,b1,castxcc(ae_int32x4,pBwr),sizeof(ae_int32x4));
                AE_S32X2X2_IP(b2,b3,castxcc(ae_int32x4,pBwr),sizeof(ae_int32x4));
                AE_S32X2X2_IP(b4,b5,castxcc(ae_int32x4,pBwr),sizeof(ae_int32x4));
                AE_S32X2X2_IP(c0,c1,castxcc(ae_int32x4,pCwr),sizeof(ae_int32x4));
            }
            __Pragma("no_reorder")
            pB=(ae_int32x2 *)(B+k*N0);
            pT=(ae_int32x2 *)(T+k);
            pC=(ae_int32x2 *)(C+k);
            AE_L32_XC     (Ckn,castxcc(ae_int32,pC),sizeof(int32_t));
            AE_L32X2X2_XC2(Bk0,Bk1,castxcc(ae_int32x4,pB),N0/2*sizeof(ae_int32x2));
            AE_L32_XC1(Ti, castxcc(ae_int32,pT),sizeof(int32_t));
            pBwr=pB;
            pCwr=pC;
            esh_coef=AE_SLAA32S(0x40000000,-e+1);
            {
                ae_int32x2 Bin0,Bin1,Bin2,Bin3,Cin0,Cin1,Ti0,Ti1;
                AE_L32_XC1     (Ti0 ,castxcc(ae_int32,pT),sizeof(int32_t));
                AE_L32_XC1     (Ti1 ,castxcc(ae_int32,pT),sizeof(int32_t));
                AE_L32_XC      (Cin0,castxcc(ae_int32,pC),sizeof(int32_t));
                AE_L32_XC      (Cin1,castxcc(ae_int32,pC),sizeof(int32_t));
                AE_L32X2X2_XC2 (Bin0,Bin1,castxcc(ae_int32x4,pB),N0/2*sizeof(ae_int32x2));
                AE_L32X2X2_XC2 (Bin2,Bin3,castxcc(ae_int32x4,pB),N0/2*sizeof(ae_int32x2));
                AE_MULF2P32X4RAS (Cin0,Cin1,Cin0,Cin1,esh_coef,esh_coef);
                AE_MULSF2P32X4RAS (Cin0,Cin1,Ckn,Ckn,Ti0,Ti1);
                AE_MULF2P32X4RAS (Bin0,Bin1,Bin0,Bin1,esh_coef,esh_coef);
                AE_MULF2P32X4RAS (Bin2,Bin3,Bin2,Bin3,esh_coef,esh_coef);
                AE_MULSF2P32X4RAS(Bin0,Bin1,Bk0,Bk1,Ti0,Ti0);
                AE_MULSF2P32X4RAS(Bin2,Bin3,Bk0,Bk1,Ti1,Ti1);
                AE_S32_L_XC    (Cin0,castxcc(ae_int32,pCwr),sizeof(int32_t));
                AE_S32_L_XC    (Cin1,castxcc(ae_int32,pCwr),sizeof(int32_t));
                AE_S32X2X2_XC2 (Bin0,Bin1,castxcc(ae_int32x4,pBwr),N0/2*sizeof(ae_int32x2));
                AE_S32X2X2_XC2 (Bin2,Bin3,castxcc(ae_int32x4,pBwr),N0/2*sizeof(ae_int32x2));
            }
            Ckn=AE_MULFP32X2RAS(Ckn,Ti);
            AE_S32_L_XC (Ckn,castxcc(ae_int32,pCwr),sizeof(int32_t));
            AE_MULF2P32X4RAS(Bk0,Bk1,Bk0,Bk1,Ti,Ti);
            AE_S32X2X2_XC2 (Bk0,Bk1,castxcc(ae_int32x4,pBwr),N0/2*sizeof(ae_int32x2));
            pBwr=(ae_int32x2*)(B+k);
            AE_S32_L_XP(AE_ZERO32(),castxcc(ae_int32,pBwr),N0*sizeof(int32_t));
            AE_S32_L_XP(AE_ZERO32(),castxcc(ae_int32,pBwr),N0*sizeof(int32_t));
            AE_S32_L_XP(AE_ZERO32(),castxcc(ae_int32,pBwr),N0*sizeof(int32_t));
            __Pragma("no_reorder");
        }
    }
    // copy back to the output
    {
        ae_valign aY=AE_ZALIGN64();
        ae_int32x2 c0,c1;
        pY=(ae_int32x2*)y;
        pC=(ae_int32x2*)C;
        AE_L32X2X2_I(c0,c1,(const ae_int32x4*)pC,0*sizeof(ae_int32x2));
        AE_SA32X2_IP(c0,aY,pY);
        AE_SA64POS_FP(aY,pY);
        AE_S32_L_I(AE_SEL32_HH(c1,c1),(ae_int32*)pY,0);
    }
    return qC;
}

// scratch allocation
size_t mtx_gjelim3x3_32x32_getScratchSize  ()
{
    const int N=3;
    int N0= (N+1)&~1;
    return (N*N0+2*N0)*sizeof(int32_t);
}
