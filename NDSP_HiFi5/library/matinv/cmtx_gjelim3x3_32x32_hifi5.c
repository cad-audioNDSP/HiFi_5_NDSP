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
 * complex Matrix Gauss-Jordan Elimination for linear equation problem 3x3, 32-bit 
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
int  cmtx_gjelim3x3_32x32  (void* pScr, complex_fract32 *y, const complex_fract32* A,const complex_fract32 * x, int qA, int qX)
{
    const int N=3;
    const int N0=4;
    xtbool2 bnegim=AE_LT32(AE_MOVDA32X2(0,-1),AE_ZERO32());
    int k,n;
    ae_int32x2 *B; // [N][N0]
    ae_int32x2 *C; // [N0]
    ae_int32x2 *T; // [N]
    ae_int32x2 * restrict pB; // [N][N0]
    ae_int32x2 * restrict pC; // [N0]
    ae_int32x2 * restrict pT; // [N0]
    const ae_int32x2 * restrict pX;
          ae_int32x2 * restrict pY;
    ae_int32x2 * restrict pBwr; // [N][N0]
    ae_int32x2 * restrict pCwr; //  [N0]

    int qB,qC; // fixed point representations
    NASSERT_ALIGN16(pScr);

    // allocate on scratch
    B=(ae_int32x2 *)pScr;
    C=B+N*N0;
    T=C+N0;
    // setup circular pointers for B and T
    WUR_AE_CBEGIN2((uintptr_t)B);    WUR_AE_CEND2((uintptr_t)(B+N*N0));
    WUR_AE_CBEGIN1((uintptr_t)T);    WUR_AE_CEND1((uintptr_t)(T+N));
    WUR_AE_CBEGIN0((uintptr_t)C);    WUR_AE_CEND0((uintptr_t)(C+N));
    // copy input
    pB=B;
    pX=(const ae_int32x2*)A;
    for (k=0; k<N; k++) 
    {
        ae_int32x2 b0,b1,b2;
        AE_L32X2_IP(b0,pX,sizeof(ae_int32x2));
        AE_L32X2_IP(b1,pX,sizeof(ae_int32x2));
        AE_L32X2_IP(b2,pX,sizeof(ae_int32x2));
        AE_S32X2X2_IP(b0,b1,castxcc(ae_int32x4,pB),sizeof(ae_int32x4));
        AE_S32X2X2_IP(b2,0 ,castxcc(ae_int32x4,pB),sizeof(ae_int32x4));
    }
    pX=(const ae_int32x2*)x;
    {
        ae_int32x2 b0,b1,b2;
        b0=AE_L32X2_I(pX,0*sizeof(ae_int32x2));
        b1=AE_L32X2_I(pX,1*sizeof(ae_int32x2));
        b2=AE_L32X2_I(pX,2*sizeof(ae_int32x2));
        AE_S32X2X2_I(b0,b1,(ae_int32x4*)C,0*sizeof(ae_int32x4));
        AE_S32X2X2_I(b2, 0,(ae_int32x4*)C,1*sizeof(ae_int32x4));
    }
    qB=31;
    qC=qX+(31-qA); // representation of inverted matrix

    for (k=0; k<N; k++)
    {
        int imax;
        int e,expB,expC;
        // find matrix normalization scale
        pB=(ae_int32x2 *)B;
        pC=(ae_int32x2 *)C;
        {
            ae_int16x4 minnsaB=31,minnsaC=31;
            for(n=0; n<(N); n++) 
            {
                ae_int32x2 b0,b1,b2,b3,c0;
                AE_L32X2X2_IP(b0,b1,castxcc(ae_int32x4,pB),sizeof(ae_int32x4));
                AE_L32X2X2_IP(b2,b3,castxcc(ae_int32x4,pB),sizeof(ae_int32x4));
                AE_L32X2_IP  (c0,   castxcc(ae_int32x2,pC),sizeof(ae_int32x2));
                minnsaB=AE_MIN16(minnsaB,AE_MIN16(AE_NSA32X4(b0,b1),AE_NSA32X4(b2,b3)));
                minnsaC=AE_MIN16(minnsaC,AE_NSA32X4(c0,c0));
            }
            expB=AE_RMIN16X4(minnsaB);
            expC=AE_RMIN16X4(minnsaC);
        }
        // pivoting
        {
            ae_int64 bmax64;
            imax=k; bmax64=AE_ZERO64();
            pB=&B[k*N0+k];
            for (n=k; n<N; n++)
            {
                ae_int32x2 bb;
                ae_int64 b;
                xtbool bbig;
                AE_L32X2_XP(bb,pB,sizeof(ae_int32x2)*N0);
                b=AE_MULZAAD32_HH_LL(bb,bb);
                bbig=AE_LE64(bmax64,b);
                AE_MOVT64(bmax64,b,bbig);
                XT_MOVT(imax, n, bbig);
            }
        }
        int off=(int)((imax-k)*sizeof(ae_int32x2)*N0);
        pB=&B[k*N0];
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
            pB=(B+k);
            {
                ae_int32x2 b0,b1,b2;
                AE_L32X2_XP(b0,pB,sizeof(ae_int32x2)*N0);
                AE_L32X2_XP(b1,pB,sizeof(ae_int32x2)*N0);
                AE_L32X2_XP(b2,pB,sizeof(ae_int32x2)*N0);
                AE_S32X2X2_I(b0,b1,(ae_int32x4*)T,0*sizeof(ae_int32x2));
                AE_S32X2_I  (b2,   (ae_int32x2*)T,2*sizeof(ae_int32x2));
            }
            Tk=AE_L32X2_X(T,k*sizeof(ae_int32x2));
            AE_S32X2_X(AE_SLAA32S(AE_MOVDA32X2(1,0),qB),T,k*sizeof(ae_int32x2));
            __Pragma("no_reorder");
            d2=AE_MULZAAD32_HH_LL(Tk,Tk);
            {
                ae_int32x2 T0,T1,T2;
                AE_L32X2X2_I(T0,T1,(const ae_int32x4*)T,0);
                T2=AE_L32X2_I((ae_int32x2*)T,2*sizeof(ae_int32x2));
                n2=AE_MULZAAD32_HH_LL(T0,T0);
                n2=AE_MAX64(n2,AE_MULZAAD32_HH_LL(T1,T1));
                n2=AE_MAX64(n2,AE_MULZAAD32_HH_LL(T2,T2));
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
            pT=T;
            Tk=AE_MULFP32X2RAS(Tk,rden);
            AE_MOVT32X2(Tk,AE_NEG32S(Tk),bnegim);
            Tk=AE_SLAA32S(Tk,e_den+e_rden-1);
            pT[0]=AE_MULFC32RAS(AE_SLAA32S(pT[0],e_num),Tk);
            pT[1]=AE_MULFC32RAS(AE_SLAA32S(pT[1],e_num),Tk);
            pT[2]=AE_MULFC32RAS(AE_SLAA32S(pT[2],e_num),Tk);
            e=e_den-e_num+1;
        }
        // scale B and C (2-way shifts possible!)
        qB=qB+expB-e;
        qC=qC+expC-e;
        // Gauss-Jordan elimination: might be made using curcular addressing on B/C
        {
            ae_int32x2 esh_coef,expB_coef,expC_coef;
            ae_int32x2 Ti,Ckn;
            ae_int32x2 Bk0,Bk1,Bk2;
            pB=(ae_int32x2 *)B;
            pC=(ae_int32x2 *)C;
            expB_coef=AE_SLAA32S(1,expB);
            expC_coef=AE_SLAA32S(1,expC);
            {
                ae_int32x2 b0,b1,b2,b3,b4,b5,b6,b7,b8,c0,c1,c2;
                AE_L32X2X2_I (c0,c1,(const ae_int32x4*)pC,0*sizeof(ae_int32x4));
                c2=AE_L32X2_I(      (const ae_int32x2*)pC,1*sizeof(ae_int32x4));
                AE_L32X2X2_I (b0,b1,(const ae_int32x4*)pB,0*sizeof(ae_int32x4));
                b2=AE_L32X2_I(      (const ae_int32x2*)pB,1*sizeof(ae_int32x4));
                AE_L32X2X2_I (b3,b4,(const ae_int32x4*)pB,2*sizeof(ae_int32x4));
                b5=AE_L32X2_I(      (const ae_int32x2*)pB,3*sizeof(ae_int32x4));
                AE_L32X2X2_I (b6,b7,(const ae_int32x4*)pB,4*sizeof(ae_int32x4));
                b8=AE_L32X2_X(      (const ae_int32x2*)pB,5*sizeof(ae_int32x4));

                AE_MUL2P32X4(c0,c1,c0,c1,expC_coef,expC_coef);
                AE_MUL2P32X4(b0,b1,b0,b1,expB_coef,expB_coef);
                AE_MUL2P32X4(c2,b2,c2,b2,expC_coef,expB_coef);
                AE_MUL2P32X4(b3,b4,b3,b4,expB_coef,expB_coef);
                AE_MUL2P32X4(b6,b7,b6,b7,expB_coef,expB_coef);
                AE_MUL2P32X4(b5,b8,b5,b8,expB_coef,expB_coef);

                AE_S32X2X2_I (c0,c1,(ae_int32x4*)pC,0*sizeof(ae_int32x4));
                AE_S32X2_I   (c2,   (ae_int32x2*)pC,1*sizeof(ae_int32x4));
                AE_S32X2X2_I (b0,b1,(ae_int32x4*)pB,0*sizeof(ae_int32x4));
                AE_S32X2_I   (b2,   (ae_int32x2*)pB,1*sizeof(ae_int32x4));
                AE_S32X2X2_I (b3,b4,(ae_int32x4*)pB,2*sizeof(ae_int32x4));
                AE_S32X2_I   (b5,   (ae_int32x2*)pB,3*sizeof(ae_int32x4));
                AE_S32X2X2_I (b6,b7,(ae_int32x4*)pB,4*sizeof(ae_int32x4));
                AE_S32X2_X   (b8,   (ae_int32x2*)pB,5*sizeof(ae_int32x4));
            }

            pB=(B+k*N0);
            pT=(T+k);
            pC=(C+k);
            Bk2=AE_L32X2_I(pB,2*sizeof(ae_int32x2));
            AE_L32X2X2_XC2(Bk0,Bk1,castxcc(ae_int32x4,pB),N0*sizeof(ae_int32x2));
            AE_L32X2_XC (Ckn,pC,sizeof(ae_int32x2));
            AE_ADDCIRC32X2_XC1(pT,sizeof(ae_int32x2));
            pCwr=pC;
            pBwr=pB;
            esh_coef=AE_SLAA32S(0x40000000,-e+1);
            {
                ae_int32x2 Bin0,Bin1,Bin2;
                ae_int32x2 Bin3,Bin4,Bin5;
                ae_int32x2 Cin0,Cin1,Ti0,Ti1;
                AE_L32X2_XC1(Ti0,pT,sizeof(ae_int32x2));
                AE_L32X2_XC1(Ti1,pT,sizeof(ae_int32x2));
                AE_L32X2_XC (Cin0,pC,sizeof(ae_int32x2));
                AE_L32X2_XC (Cin1,pC,sizeof(ae_int32x2));
                Bin2=AE_L32X2_I(pB,2*sizeof(ae_int32x2));
                AE_L32X2X2_XC2 (Bin0,Bin1,castxcc(ae_int32x4,pB),N0*sizeof(ae_int32x2));
                Bin5=AE_L32X2_I(pB,2*sizeof(ae_int32x2));
                AE_L32X2X2_XC2 (Bin3,Bin4,castxcc(ae_int32x4,pB),N0*sizeof(ae_int32x2));
                Ti0=AE_NEG32S(Ti0);
                Ti1=AE_NEG32S(Ti1);
                AE_MULF2P32X4RAS(Cin0,Cin1,Cin0,Cin1,esh_coef,esh_coef);
                AE_MULAFC32RAS(Cin0,Ckn,Ti0);
                AE_MULAFC32RAS(Cin1,Ckn,Ti1);
                AE_MULF2P32X4RAS(Bin0,Bin1,Bin0,Bin1,esh_coef,esh_coef);
                AE_MULF2P32X4RAS(Bin3,Bin4,Bin3,Bin4,esh_coef,esh_coef);
                AE_MULF2P32X4RAS(Bin2,Bin5,Bin2,Bin5,esh_coef,esh_coef);
                AE_MULAFC32RAS(Bin0,Bk0,Ti0);
                AE_MULAFC32RAS(Bin1,Bk1,Ti0);
                AE_MULAFC32RAS(Bin2,Bk2,Ti0);
                AE_MULAFC32RAS(Bin3,Bk0,Ti1);
                AE_MULAFC32RAS(Bin4,Bk1,Ti1);
                AE_MULAFC32RAS(Bin5,Bk2,Ti1);

                AE_S32X2_XC   (Cin0,pCwr,sizeof(ae_int32x2));
                AE_S32X2_XC   (Cin1,pCwr,sizeof(ae_int32x2));
                AE_S32X2_I    (Bin2,pBwr,2*sizeof(ae_int32x2));
                AE_S32X2X2_XC2(Bin0,Bin1,castxcc(ae_int32x4,pBwr),N0*sizeof(ae_int32x2));
                AE_S32X2_I    (Bin5,pBwr,2*sizeof(ae_int32x2));
                AE_S32X2X2_XC2(Bin3,Bin4,castxcc(ae_int32x4,pBwr),N0*sizeof(ae_int32x2));
            }
            Ti=AE_L32X2_I(pT,0);
            Ckn=AE_MULFC32RAS(Ckn,Ti);
            AE_S32X2_XC   (Ckn,pCwr,sizeof(ae_int32x2));
            Bk0=AE_MULFC32RAS(Bk0,Ti);
            Bk1=AE_MULFC32RAS(Bk1,Ti);
            Bk2=AE_MULFC32RAS(Bk2,Ti);
            AE_S32X2_I(Bk2,pBwr,2*sizeof(ae_int32x2));
            AE_S32X2X2_XC2(Bk0,Bk1,castxcc(ae_int32x4,pBwr),N0*sizeof(ae_int32x2));
            pBwr=(B+k);
            AE_S32X2_XP(AE_ZERO32(),pBwr,N0*sizeof(ae_int32x2));
            AE_S32X2_XP(AE_ZERO32(),pBwr,N0*sizeof(ae_int32x2));
            AE_S32X2_XP(AE_ZERO32(),pBwr,N0*sizeof(ae_int32x2));
            __Pragma("no_reorder");
        }
    }
    // copy back to the output
    {
        ae_int32x2 c0,c1,c2;
        pY=(ae_int32x2 *)y;
        AE_L32X2X2_I(c0,c1,(const ae_int32x4*)C,0*sizeof(ae_int32x2));
        c2=AE_L32X2_I((const ae_int32x2*)C,2*sizeof(ae_int32x2));
        AE_S32X2_I(c0,pY,0*sizeof(ae_int32x2));
        AE_S32X2_I(c1,pY,1*sizeof(ae_int32x2));
        AE_S32X2_I(c2,pY,2*sizeof(ae_int32x2));
    }
    return qC;
}
// scratch allocation
size_t cmtx_gjelim3x3_32x32_getScratchSize   () 
{ 
    return  (3*4+2*4)*sizeof(complex_fract32);  
}
