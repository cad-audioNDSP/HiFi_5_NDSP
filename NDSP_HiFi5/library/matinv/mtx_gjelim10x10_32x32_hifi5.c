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
 * Real Matrix Gauss-Jordan Elimination for linear equation problem 10x10, 32-bit 
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
int  mtx_gjelim10x10_32x32(void* pScr, int32_t *y, const int32_t* A,const int32_t * x, int qA, int qX) 
{
    const int N=10;
    const int N0=12;
    int k,n;
    int32_t *B; // [N][N]
    int32_t *C; // [N]
    int32_t *T; // [N]
    ae_int32x2 * restrict pB; // [N][N]
    ae_int32x2 * restrict pC; // [N]
    ae_int32x2 * restrict pT; // [N]
    const ae_int32x2 * restrict pX; // [N]
          ae_int32x2 * restrict pY; // [N]
    ae_int32x2 * restrict pBwr; // [N][N]
    ae_int32x2 * restrict pCwr; // [N][N]

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
    WUR_AE_CBEGIN0((uintptr_t)B);    WUR_AE_CEND0((uintptr_t)(B+N*N0));
    WUR_AE_CBEGIN1((uintptr_t)T);    WUR_AE_CEND1((uintptr_t)(T+N));
    // copy input
    {
        ae_valign aX;
        pB=(ae_int32x2*)B;
        pX=(const ae_int32x2*)A;
        aX=AE_LA64_PP(pX);
        for (k=0; k<(N); k++) 
        {
            ae_int32x2 b0,b1,b2,b3,b4;
            AE_LA32X2_IP(b0,aX,pX);
            AE_LA32X2_IP(b1,aX,pX);
            AE_LA32X2_IP(b2,aX,pX);
            AE_LA32X2_IP(b3,aX,pX);
            AE_LA32X2_IP(b4,aX,pX);
            AE_S32X2X2_IP(b0,b1,castxcc(ae_int32x4,pB),sizeof(ae_int32x4));
            AE_S32X2X2_IP(b2,b3,castxcc(ae_int32x4,pB),sizeof(ae_int32x4));
            AE_S32X2X2_IP(b4, 0,castxcc(ae_int32x4,pB),sizeof(ae_int32x4));
        }
        pC=(ae_int32x2*)C;
        pX=(const ae_int32x2*)x;
        aX=AE_LA64_PP(pX);
        {
            ae_int32x2 c0,c1,c2,c3,c4;
            AE_LA32X2_IP(c0,aX,pX);
            AE_LA32X2_IP(c1,aX,pX);
            AE_LA32X2_IP(c2,aX,pX);
            AE_LA32X2_IP(c3,aX,pX);
            AE_LA32X2_IP(c4,aX,pX);
            AE_S32X2X2_I(c0,c1,(ae_int32x4*)pC,0*sizeof(ae_int32x4));
            AE_S32X2X2_I(c2,c3,(ae_int32x4*)pC,1*sizeof(ae_int32x4));
            AE_S32X2X2_I(c4, 0,(ae_int32x4*)pC,2*sizeof(ae_int32x4));
        }
    }
    qB=31;
    qC=qX+(31-qA); // representation of inverted matrix

    for (k=0; k<N; k++)
    {
        int i,imax;
        int e,expB,expC;
        // find matrix normalization scale
        // find matrix normalization scale
        ae_int16x4 minnsaB0=31,minnsaB1=31,minnsaC=31;
        pB=(ae_int32x2*)B;
        pC=(ae_int32x2*)C;
        for(n=0; n<((N*N0)>>3); n++) 
        {
            ae_int32x2 b0,b1,b2,b3;
            AE_L32X2X2_IP(b0,b1,castxcc(ae_int32x4,pB),sizeof(ae_int32x4));
            AE_L32X2X2_IP(b2,b3,castxcc(ae_int32x4,pB),sizeof(ae_int32x4));
            minnsaB0=AE_MIN16(minnsaB0,AE_NSA32X4(b0,b1));
            minnsaB1=AE_MIN16(minnsaB1,AE_NSA32X4(b2,b3));
        }
        for(n=0; n<(N0>>2); n++) 
        {
            ae_int32x2 c0,c1;
            AE_L32X2X2_IP(c0,c1,castxcc(ae_int32x4,pC),sizeof(ae_int32x4));
            minnsaC=AE_MIN16(minnsaC,AE_NSA32X4(c0,c1));
        }
        expB=AE_RMIN16X4(AE_MIN16(minnsaB0,minnsaB1));
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
            for (n=0; n<(N/2); n++) 
            {
                bk=AE_L32X2_I(pB,0);
                bi=AE_L32X2_X(pB,off);
                AE_S32X2_X (bk,pB,off);
                AE_S32X2_IP(bi,pB,sizeof(ae_int32x2));
            }
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
            int e_den,e_num;
            pT=(ae_int32x2*)T;
            pB=(ae_int32x2*)&B[k];
            for (n=0; n<N/2; n++) 
            {
                ae_int32x2 b0,b1;
                AE_L32_XP(b0,castxcc(ae_int32,pB),N0*sizeof(int32_t));
                AE_L32_XP(b1,castxcc(ae_int32,pB),N0*sizeof(int32_t));
                AE_S32X2_IP(AE_SEL32_HH(b0,b1),pT,sizeof(ae_int32x2));
            }
            AE_S32X2_I(0,pT,0);
            pT=(ae_int32x2*)T;
            tk=AE_L32_X((const ae_int32*)pT,k*sizeof(int32_t));
            AE_S32_L_X(AE_SLAA32S(1,qB),(ae_int32*)pT,k*sizeof(int32_t));
            {
                ae_int32x2 x,y,ee;
                xtbool2 sx;
                ae_int32x2 t0,t1,t2,t3,t4,t5;
                AE_L32X2X2_I(t0,t1,(const ae_int32x4*)pT,0*sizeof(ae_int32x4));
                AE_L32X2X2_I(t2,t3,(const ae_int32x4*)pT,1*sizeof(ae_int32x4));
                AE_L32X2X2_I(t4,t5,(const ae_int32x4*)pT,2*sizeof(ae_int32x4));
                maxnum=AE_MAX32(AE_MAX32(AE_MAXABS32S(t0,t1),AE_MAXABS32S(t2,t3)),AE_MAXABS32S(t4,t5));
                maxnum=AE_MAX32(maxnum,AE_SEL32_LH(maxnum,maxnum));
                maxden=AE_ABS32S(tk);
                maxnum=AE_ADD32S(AE_SRAI32(maxnum,1),AE_SRAI32(maxden,1));
                e_den=AE_NSAZ32_L(maxden); 
                e_num=AE_NSAZ32_L(maxnum)-1; 
                e=e_den-e_num+1;
                tk=AE_SLAA32S(tk,e_den);
               // reciprocal
                x=tk;
                sx=AE_LT32(x,AE_ZERO32());
                x=AE_ABS32S(x);
                y=AE_SUB32(0xBAEC0000,x);
                ee=0x40000000; AE_MULSFP32X2RAS(ee,x,y);  AE_MULAFP32X2RAS(y,y,AE_SLAI32(ee,1));
                ee=0x40000000; AE_MULSFP32X2RAS(ee,x,y);  AE_MULAFP32X2RAS(y,y,AE_SLAI32(ee,1));
                ee=0x40000000; AE_MULSFP32X2RAS(ee,x,y);  AE_MULAFP32X2RAS(y,y,AE_SLAI32(ee,1));
                ee=0x40000000; AE_MULSFP32X2RAS(ee,x,y);  AE_MULAFP32X2RAS(y,y,AE_SLAI32(ee,1));
                AE_MOVT32X2(y,AE_NEG32(y),sx);
                rden=y;

                t0=AE_SLAA32S(t0,e_num);
                t1=AE_SLAA32S(t1,e_num);
                t2=AE_SLAA32S(t2,e_num);
                t3=AE_SLAA32S(t3,e_num);
                t4=AE_SLAA32S(t4,e_num);
                t5=AE_SLAA32S(t5,e_num);
                AE_MULF2P32X4RAS(t0,t1,t0,t1,rden,rden);
                AE_MULF2P32X4RAS(t2,t3,t2,t3,rden,rden);
                AE_MULF2P32X4RAS(t4,t5,t4,t5,rden,rden);
                AE_S32X2X2_I(t0,t1,(ae_int32x4*)pT,0*sizeof(ae_int32x4));
                AE_S32X2X2_I(t2,t3,(ae_int32x4*)pT,1*sizeof(ae_int32x4));
                AE_S32X2X2_I(t4,t5,(ae_int32x4*)pT,2*sizeof(ae_int32x4));
            }
        }
        // scale B and C (2-way shifts possible!)
        qB=qB+expB-e;
        qC=qC+expC-e;
        // Gauss-Jordan elimination: might be made using curcular addressing on B/C
        {
            ae_int32x2 esh_coef,expB_coef,expC_coef;
            ae_int32x2 Cin;
            ae_int32x2 Ti,Ckn;
            ae_int32x2 Bk0,Bk1,Bk2,Bk3,Bk4,Bk5;
            pB=(ae_int32x2 *)B;
            pC=(ae_int32x2 *)C;
            pBwr=(ae_int32x2 *)B;
            pCwr=(ae_int32x2 *)C;
            expB_coef=AE_SLAA32S(1,expB);
            expC_coef=AE_SLAA32S(1,expC);
            for (n=0; n<((N*N0)>>3); n++)
            {
                ae_int32x2 b0,b1,b2,b3;
                AE_L32X2X2_IP (b0,b1,castxcc(ae_int32x4,pB),sizeof(ae_int32x4));
                AE_L32X2X2_IP (b2,b3,castxcc(ae_int32x4,pB),sizeof(ae_int32x4));
                AE_MUL2P32X4(b0,b1,b0,b1,expB_coef,expB_coef);
                AE_MUL2P32X4(b2,b3,b2,b3,expB_coef,expB_coef);
                AE_S32X2X2_IP (b0,b1,castxcc(ae_int32x4,pBwr),sizeof(ae_int32x4));
                AE_S32X2X2_IP (b2,b3,castxcc(ae_int32x4,pBwr),sizeof(ae_int32x4));
            }
            {
                ae_int32x2 c0,c1,c2,c3,c4,c5;
                AE_L32X2X2_IP (c0,c1,castxcc(ae_int32x4,pC),sizeof(ae_int32x4));
                AE_L32X2X2_IP (c2,c3,castxcc(ae_int32x4,pC),sizeof(ae_int32x4));
                AE_L32X2X2_IP (c4,c5,castxcc(ae_int32x4,pC),sizeof(ae_int32x4));
                AE_MUL2P32X4(c0,c1,c0,c1,expC_coef,expC_coef);
                AE_MUL2P32X4(c2,c3,c2,c3,expC_coef,expC_coef);
                AE_MUL2P32X4(c4,c5,c4,c5,expC_coef,expC_coef);
                AE_S32X2X2_IP (c0,c1,castxcc(ae_int32x4,pCwr),sizeof(ae_int32x4));
                AE_S32X2X2_IP (c2,c3,castxcc(ae_int32x4,pCwr),sizeof(ae_int32x4));
                AE_S32X2X2_IP (c4,c5,castxcc(ae_int32x4,pCwr),sizeof(ae_int32x4));
            }
             __Pragma("no_reorder");    // important! - C[] is addressed below circularly via pT
            pC=(ae_int32x2 *)C;
            Ckn=AE_L32_X((const ae_int32*)pC,k*sizeof(int32_t));
            pB=(ae_int32x2 *)(B+k*N0);
            pT=(ae_int32x2 *)(T+k);
            AE_L32X2X2_I  (Bk2,Bk3,(const ae_int32x4*)pB,2*sizeof(ae_int32x2));
            AE_L32X2X2_I  (Bk4,Bk5,(const ae_int32x4*)pB,4*sizeof(ae_int32x2));
            AE_L32X2X2_XC (Bk0,Bk1,castxcc(ae_int32x4,pB),N0/2*sizeof(ae_int32x2));
            AE_ADDCIRC32X2_XC1(pT,sizeof(int32_t));
            pBwr=pB;
            esh_coef=AE_SLAA32S(0x40000000,-e+1);
            for (i=0; i<N-1; i++)
            {
                ae_int32x2 Bin0,Bin1,Bin2,Bin3,Bin4,Bin5;
                Ti=AE_L32_I((const ae_int32*)pT,0);
                Ti=AE_NEG32S(Ti);
                Cin=AE_L32_X ((const ae_int32*)pT,(int)((uintptr_t)C-(uintptr_t)T));
                Cin=AE_SLAA32S(Cin,-e);
                AE_MULAFP32X2RAS(Cin,Ckn,Ti);
                AE_S32_L_X (Cin,(ae_int32*)pT,(int)((uintptr_t)C-(uintptr_t)T));
                AE_ADDCIRC32X2_XC1(pT,sizeof(int32_t));
                AE_L32X2X2_I  (Bin2,Bin3,(const ae_int32x4*)pB,2*sizeof(ae_int32x2));
                AE_L32X2X2_I  (Bin4,Bin5,(const ae_int32x4*)pB,4*sizeof(ae_int32x2));
                AE_L32X2X2_XC (Bin0,Bin1,castxcc(ae_int32x4,pB),N0/2*sizeof(ae_int32x2));
                AE_MULF2P32X4RAS(Bin0,Bin1,Bin0,Bin1,esh_coef,esh_coef);
                AE_MULF2P32X4RAS(Bin2,Bin3,Bin2,Bin3,esh_coef,esh_coef);
                Bin4=AE_MULFP32X2RAS(Bin4,esh_coef);
                AE_MULAF2P32X4RAS(Bin0,Bin1,Bk0,Bk1,Ti,Ti);
                AE_MULAF2P32X4RAS(Bin2,Bin3,Bk2,Bk3,Ti,Ti);
                AE_MULAFP32X2RAS(Bin4,Bk4,Ti);
                AE_S32X2X2_I  (Bin2,Bin3,(ae_int32x4*)pBwr,2*sizeof(ae_int32x2));
                AE_S32X2X2_I  (Bin4,Bin5,(ae_int32x4*)pBwr,4*sizeof(ae_int32x2));
                AE_S32X2X2_XC (Bin0,Bin1,castxcc(ae_int32x4,pBwr),N0/2*sizeof(ae_int32x2));
            }
            Ti=AE_L32_I((const ae_int32*)pT,0);
            Ckn=AE_MULFP32X2RAS(Ckn,Ti);
            AE_S32_L_X (Ckn,(ae_int32*)pT,(int)((uintptr_t)C-(uintptr_t)T));
            AE_MULF2P32X4RAS(Bk0,Bk1,Bk0,Bk1,Ti,Ti);
            AE_MULF2P32X4RAS(Bk2,Bk3,Bk2,Bk3,Ti,Ti);
            Bk4=AE_MULFP32X2RAS(Bk4,Ti);
            AE_S32X2X2_I  (Bk2,Bk3,(ae_int32x4*)pBwr,2*sizeof(ae_int32x2));
            AE_S32X2X2_I  (Bk4,Bk5,(ae_int32x4*)pBwr,4*sizeof(ae_int32x2));
            AE_S32X2X2_XC (Bk0,Bk1,castxcc(ae_int32x4,pBwr),N0/2*sizeof(ae_int32x2));
            pBwr=(ae_int32x2*)(B+k);
            for (i=0; i<N; i++)  AE_S32_L_XP(AE_ZERO32(),castxcc(ae_int32,pBwr),N0*sizeof(int32_t));
            __Pragma("no_reorder");
        }
    }
    // copy back to the output
    {
        ae_int32x2 c0,c1,c2,c3,c4,c5;
        ae_valign aY;
        pY=(ae_int32x2*)y;
        aY=AE_ZALIGN64();
        pC=(ae_int32x2*)C;
        AE_L32X2X2_I (c0,c1,(const ae_int32x4*)pC,0*sizeof(ae_int32x4));
        AE_L32X2X2_I (c2,c3,(const ae_int32x4*)pC,1*sizeof(ae_int32x4));
        AE_L32X2X2_I (c4,c5,(const ae_int32x4*)pC,2*sizeof(ae_int32x4));
        AE_SA32X2_IP(c0,aY,pY);
        AE_SA32X2_IP(c1,aY,pY);
        AE_SA32X2_IP(c2,aY,pY);
        AE_SA32X2_IP(c3,aY,pY);
        AE_SA32X2_IP(c4,aY,pY);
        AE_SA64POS_FP(aY,pY);
    }
    return qC;
}
size_t mtx_gjelim10x10_32x32_getScratchSize () 
{
    const int N=10;
    const int N0=12;
    return (N*N0+2*N0)*sizeof(int32_t);
}
