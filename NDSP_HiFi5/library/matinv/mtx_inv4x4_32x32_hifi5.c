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
 * Real Matrix Inversion 4x4, 32-bit fixed point API
 * Optimized code for HiFi5
 */

/* Cross-platform data type definitions. */
#include "NatureDSP_types.h"
/* Common helper macros. */
#include "common.h"
/* Matrix functions */
#include "NatureDSP_Signal_matinv.h"

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
int  mtx_inv4x4_32x32  (void* pScr, int32_t *x, int qX) 
{
    const int N=4;
    int k,n;
    ae_valignx2 ax;
    int32_t *B; // [N][N]
    int32_t *C; // [N][N]
    int32_t *T; // [N]
    ae_int32x2 * restrict pB; // [N][N]
    ae_int32x2 * restrict pC; // [N][N]
    ae_int32x2 * restrict pT; // [N]
    ae_int32x2 * restrict pX; // [N]
    ae_int32x2 * restrict pBwr; // [N][N]
    ae_int32x2 * restrict pCwr; // [N][N]

    int qB,qC; // fixed point representations
    NASSERT_ALIGN16(pScr);
    // allocate on scratch
    B=(int32_t *)pScr;
    C=B+N*N;
    T=C+N*N;
    NASSERT_ALIGN16(B);
    NASSERT_ALIGN16(C);
    NASSERT_ALIGN16(T);
    // setup circular pointers for B and T
    WUR_AE_CBEGIN0((uintptr_t)B);    WUR_AE_CEND0((uintptr_t)(B+N*N));
    WUR_AE_CBEGIN1((uintptr_t)T);    WUR_AE_CEND1((uintptr_t)(T+N));

    // copy input
    {
        pB=(ae_int32x2*)B;
        pX=(ae_int32x2*)x;
        pC=(ae_int32x2*)C;
        ax=AE_LA128_PP(pX);
        for (k=0; k<(N>>1)*(N>>1); k++) 
        {
            ae_int32x2 b0,b1;
            AE_LA32X2X2_IP(b0,b1,ax,castxcc(ae_int32x4,pX));
            AE_S32X2X2_IP (b0,b1,castxcc(ae_int32x4,pB),sizeof(ae_int32x4));
            AE_S32X2X2_IP (AE_ZERO32(),AE_ZERO32(),castxcc(ae_int32x4,pC),sizeof(ae_int32x4));
        }
        pC=(ae_int32x2*)C;
        for (k=0; k<N; k++) 
        {
            AE_S32_L_XP(0x7fffffff,castxcc(ae_int32,pC),(N+1)*sizeof(int32_t));
        }
    }
    qB=31;
    qC=31+(31-qX); // representation of inverted matrix

    for (k=0; k<N; k++)
    {
        int i,imax;
        int e,expB,expC;
        // find matrix normalization scale
        ae_int16x4 minnsaB=31,minnsaC=31;
        pB=(ae_int32x2*)B;
        pC=(ae_int32x2*)C;
        for(n=0; n<((N*N)>>2); n++) 
        {
            ae_int32x2 b0,b1,c0,c1;
            AE_L32X2X2_IP(b0,b1,castxcc(ae_int32x4,pB),sizeof(ae_int32x4));
            AE_L32X2X2_IP(c0,c1,castxcc(ae_int32x4,pC),sizeof(ae_int32x4));
            minnsaB=AE_MIN16(minnsaB,AE_NSA32X4(b0,b1));
            minnsaC=AE_MIN16(minnsaC,AE_NSA32X4(c0,c1));
        }
        expB=AE_RMIN16X4(minnsaB);
        expC=AE_RMIN16X4(minnsaC);
        // pivoting
        {
            int off;
            ae_int32x2 imax2,n2;
            ae_int32x2 bmax=AE_ZERO32();
            pB=(ae_int32x2*)(&B[k*N+k]);
            imax2=n2=k;
            for (n=k; n<N; n++)
            {
                ae_int32x2 b;
                xtbool2 bbig;
                AE_L32_XP(b,castxcc(ae_int32,pB),N*sizeof(int32_t));
                b=AE_ABS32S(b);
                bbig=AE_LE32(bmax,b);
                AE_MOVT32X2(imax2,n2,bbig);
                bmax=AE_MAX32(bmax,b);
                n2=AE_ADD32(n2,1);
            }
            imax=AE_MOVAD32_L(imax2);
            off=(int)((imax-k)*sizeof(int32_t)*N);
            pB=(ae_int32x2*)&B[k*N];
            pC=(ae_int32x2*)&C[k*N];        
            for (n=0; n<(N/2); n++) 
            {
                ae_int32x2 bk,bi,ck,ci;
                bk=AE_L32X2_I(pB,0);
                bi=AE_L32X2_X(pB,off);
                ck=AE_L32X2_I(pC,0);
                ci=AE_L32X2_X(pC,off);
                AE_S32X2_X (bk,pB,off);
                AE_S32X2_IP(bi,pB,sizeof(ae_int32x2));
                AE_S32X2_X (ck,pC,off);
                AE_S32X2_IP(ci,pC,sizeof(ae_int32x2));
            }
        }
        // find normalization factor
        {
            ae_int32x2 rden,tk,maxden,maxnum;
            int e_den,e_num;
            pT=(ae_int32x2*)T;
            pB=(ae_int32x2*)&B[k];
            for (n=0; n<N; n++) 
            {
                ae_int32x2 b;
                AE_L32_XP(b,castxcc(ae_int32,pB),N*sizeof(int32_t));
                AE_S32_L_IP(b,castxcc(ae_int32,pT),sizeof(int32_t));
            }
            pT=(ae_int32x2*)T;
            tk=AE_L32_X((const ae_int32*)pT,k*sizeof(int32_t));
            AE_S32_L_X(AE_SLAA32S(1,qB),(ae_int32*)pT,k*sizeof(int32_t));
            {
                ae_int32x2 t0,t1;
                AE_L32X2X2_I (t0,t1,(const ae_int32x4*)pT,0*sizeof(ae_int32x4));
                maxnum=AE_MAXABS32S(t0,t1);
                maxnum=AE_MAX32(maxnum,AE_SEL32_LH(maxnum,maxnum));
            }
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
            {
                ae_int32x2 t0,t1;
                AE_L32X2X2_I (t0,t1,(const ae_int32x4*)pT,0*sizeof(ae_int32x4));
                t0=AE_SLAA32S(t0,e_num);
                t1=AE_SLAA32S(t1,e_num);
                AE_MULF2P32X4RAS(t0,t1,t0,t1,rden,rden);
                AE_S32X2X2_I (t0,t1,(ae_int32x4*)pT,0*sizeof(ae_int32x4));
            }
        }
        // scale B and C (2-way shifts possible!)
        qB=qB+expB-e;
        qC=qC+expC-e;
        // Gauss-Jordan elimination: might be made using curcular addressing on B/C
        {
            ae_int32x2 Ti;
            ae_int32x2 esh_coef;
            ae_int32x2 Bk0,Bk1;
            ae_int32x2 Ck0,Ck1;
            ae_int32x2 expB_coef,expC_coef;
            pB=(ae_int32x2 *)B;
            pC=(ae_int32x2 *)C;
            pBwr=(ae_int32x2 *)B;
            pCwr=(ae_int32x2 *)C;
            expB_coef=AE_SLAA32S(1,expB);
            expC_coef=AE_SLAA32S(1,expC);
            for (n=0; n<((N*N)>>2); n++)
            {
                ae_int32x2 b0,b1,c0,c1;
                AE_L32X2X2_IP (b0,b1,castxcc(ae_int32x4,pB),sizeof(ae_int32x4));
                AE_L32X2X2_IP (c0,c1,castxcc(ae_int32x4,pC),sizeof(ae_int32x4));
                AE_MUL2P32X4(b0,b1,b0,b1,expB_coef,expB_coef);
                AE_MUL2P32X4(c0,c1,c0,c1,expC_coef,expC_coef);
                AE_S32X2X2_IP (b0,b1,castxcc(ae_int32x4,pBwr),sizeof(ae_int32x4));
                AE_S32X2X2_IP (c0,c1,castxcc(ae_int32x4,pCwr),sizeof(ae_int32x4));
            }
            __Pragma("no_reorder")
            pB=(ae_int32x2*)(B+k*N);
            pC=(ae_int32x2*)(((uintptr_t)pB)+(int)((uintptr_t)C-(uintptr_t)B));
            AE_L32X2X2_I  (Ck0,Ck1,(const ae_int32x4*)pC,0*sizeof(ae_int32x2));
            AE_L32X2X2_XC (Bk0,Bk1,castxcc(ae_int32x4,pB),2*sizeof(ae_int32x2));
            pT=(ae_int32x2*)(T+k);
            AE_ADDCIRC32X2_XC1(pT,sizeof(int32_t));
            pBwr=pB;
//            NASSERT(e>0);
            esh_coef=AE_SLAA32S(0x40000000,-e+1);
            for (i=0; i<N-1; i++)
            {
                ae_int32x2 Bin0,Bin1;
                ae_int32x2 Cin0,Cin1;
                AE_L32_XC1(Ti,castxcc(ae_int32,pT),sizeof(int32_t));
                pC=(ae_int32x2*)(((uintptr_t)pB)+(int)((uintptr_t)C-(uintptr_t)B));
                pCwr=pC;
                AE_L32X2X2_I  (Cin0,Cin1,(const ae_int32x4*)pC,0*sizeof(ae_int32x2));
                AE_L32X2X2_XC (Bin0,Bin1,castxcc(ae_int32x4,pB),2*sizeof(ae_int32x2));

                AE_MULF2P32X4RAS(Bin0,Bin1,Bin0,Bin1,esh_coef,esh_coef);
                AE_MULF2P32X4RAS(Cin0,Cin1,Cin0,Cin1,esh_coef,esh_coef);
                AE_MULSF2P32X4RAS(Bin0,Bin1,Bk0,Bk1,Ti,Ti);
                AE_MULSF2P32X4RAS(Cin0,Cin1,Ck0,Ck1,Ti,Ti);

                AE_S32X2X2_I  (Cin0,Cin1,(      ae_int32x4*)pCwr,0*sizeof(ae_int32x2));
                AE_S32X2X2_XC (Bin0,Bin1,castxcc(ae_int32x4,pBwr),2*sizeof(ae_int32x2));
            }
            AE_L32_XC1(Ti,castxcc(ae_int32,pT),sizeof(int32_t));
            pCwr=(ae_int32x2*)(((uintptr_t)pBwr)+(int)((uintptr_t)C-(uintptr_t)B));
            AE_MULF2P32X4RAS(Bk0,Bk1,Bk0,Bk1,Ti,Ti);
            AE_MULF2P32X4RAS(Ck0,Ck1,Ck0,Ck1,Ti,Ti);
            AE_S32X2X2_I  (Ck0,Ck1,(      ae_int32x4*)pCwr,0*sizeof(ae_int32x2));
            AE_S32X2X2_XC (Bk0,Bk1,castxcc(ae_int32x4,pBwr),2*sizeof(ae_int32x2));
            __Pragma("no_reorder");
        }
        pB=(ae_int32x2*)(B+k);
        for (i=0; i<N; i++)  AE_S32_L_XP(AE_ZERO32(),castxcc(ae_int32,pB),N*sizeof(int32_t));
    }
    // copy back to the output
    {
        pX=(ae_int32x2*)x;
        pC=(ae_int32x2*)C;
        ax=AE_ZALIGN128();
        for (k=0; k<(N>>1)*(N>>1); k++) 
        {
            ae_int32x2 b0,b1;
            AE_L32X2X2_IP(b0,b1,castxcc(ae_int32x4,pC),sizeof(ae_int32x4));
            AE_SA32X2X2_IP(b0,b1,ax,castxcc(ae_int32x4,pX));
        }
        AE_SA128POS_FP(ax,pX);
    }
    return qC;
}

size_t mtx_inv4x4_32x32_getScratchSize   () 
{ 
    const int N=4; 
    return (2*N*N+N)*sizeof(int32_t); 
}
