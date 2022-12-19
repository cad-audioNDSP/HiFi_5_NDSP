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
  NatureDSP Signal Processing Library. Vector matematics
    Divide 64x64 
    Code optimized for HiFi5
  IntegrIT, 2006-2019
*/

#include "NatureDSP_Signal_math.h"
#include "NatureDSP_types.h"
#include "common.h"
#define SCR_SZ (MAX_ALLOCA_SZ/(sizeof(int64_t)*2))

/*
    internal functon called with N multiple of 8 and aligned scratch
*/
static void vec_divide64x64_internal
                (int64_t * restrict frac, 
                 int16_t *exp,
                 const int64_t * restrict x, 
                 const int64_t * restrict y, int N,
                 ae_int64* scratch)
{
    const ae_int64* restrict pX;
    const ae_int64* restrict pY;
    const ae_int64* restrict pScrRd;
          ae_int64* restrict pScrWr;
    const ae_int64* restrict pFractRd;
          ae_int64* restrict pFractWr;
          ae_int16* restrict pExp;
    int n,m,M;
    NASSERT(N>0 && N%8==0);

    for (m=0; m<N; m+=SCR_SZ,x+=SCR_SZ,y+=SCR_SZ,frac+=SCR_SZ,exp+=SCR_SZ)
    {
        ae_valignx2 aF;
        M=XT_MIN(N-m,SCR_SZ);
        NASSERT(M%8==0);
        /*
            normalization phase. Here we compute output exponent, 
            save the sign of numerator in the LSB of exp[], and normalized 
            numerator/denominator in the scratch memory
        */
        aF=AE_ZALIGN128();
        pExp    =(ae_int16*) exp;
        pFractWr=(ae_int64*)frac;
        pScrWr  =(ae_int64*)scratch;
        pX      =(const ae_int64*)x;
        pY      =(const ae_int64*)y;
        __Pragma("loop_count factor=4,min=1")
        for (n=0; n<M; n+=2)
        {
            ae_int64 x0,y0;
            ae_int64 x1,y1;
            int expx0,expy0,exponent0;
            int expx1,expy1,exponent1;
            xtbool sy0, bzero0;
            xtbool sy1, bzero1;

            AE_L64_IP(x0,pX,sizeof(int64_t));
            AE_L64_IP(y0,pY,sizeof(int64_t));
            AE_L64_IP(x1,pX,sizeof(int64_t));
            AE_L64_IP(y1,pY,sizeof(int64_t));

            sy0=AE_LT64(y0,AE_ZERO64());
            sy1=AE_LT64(y1,AE_ZERO64());
            bzero0=AE_EQ64(y0,AE_ZERO64());
            bzero1=AE_EQ64(y1,AE_ZERO64());
            y0=AE_ABS64S(y0);
            y1=AE_ABS64S(y1);
            expx0=AE_NSA64(x0);
            expx1=AE_NSA64(x1);
            expy0=AE_NSA64(y0);
            expy1=AE_NSA64(y1);
            x0=AE_SLAA64(x0,expx0);
            x1=AE_SLAA64(x1,expx1);
            y0=AE_SLAA64(y0,expy0);
            y1=AE_SLAA64(y1,expy1);
            exponent0=expy0-expx0+1;
            exponent1=expy1-expx1+1;
            XT_MOVNEZ(exponent0,0x40,AE_MOVAB(bzero0));
            XT_MOVNEZ(exponent1,0x40,AE_MOVAB(bzero1));
            // set numerator MAX_INT64/MIN_INT64 according to the 
            // sign if denominator is zero
            AE_MOVT64(x0,AE_SUB64(MAX_INT64,AE_SRAI64(x0,63)),bzero0);
            AE_MOVT64(x1,AE_SUB64(MAX_INT64,AE_SRAI64(x1,63)),bzero1);
            AE_MOVT64(y0,0x4000000000000000LL,bzero0);
            AE_MOVT64(y1,0x4000000000000000LL,bzero1);
            exp[n+0]=(exponent0<<1)|AE_MOVAB(sy0);
            exp[n+1]=(exponent1<<1)|AE_MOVAB(sy1);
            AE_SA32X2X2_IP(AE_MOVINT32X2_FROMINT64(x0),AE_MOVINT32X2_FROMINT64(x1),aF,castxcc(ae_int32x4,pFractWr));
            AE_S64X2_IP(y0,y1,castxcc(ae_int64x2,pScrWr),2*sizeof(ae_int64x2));
        }
        AE_SA128POS_FP(aF,pFractWr);
        /* at the next phase, we compute approximate reciprocal 
           in 32-bit accuracy using SIMD operations
        */
        __Pragma("no_reorder")
        pScrRd  =(const ae_int64*)scratch;
        pScrWr  =(      ae_int64*)(scratch+2);
        __Pragma("loop_count factor=2,min=1")
        for (n=0; n<M; n+=4)
        {
            ae_int32x2 e0,x0,y0;
            ae_int32x2 e1,x1,y1;
            ae_int64 z0,z1,z2,z3;
            AE_L64X2_IP(z0,z1,castxcc(ae_int64x2,pScrRd),2*sizeof(ae_int64x2));
            AE_L64X2_IP(z2,z3,castxcc(ae_int64x2,pScrRd),2*sizeof(ae_int64x2));
            x0=AE_SEL32_HH(AE_MOVINT32X2_FROMINT64(z0),AE_MOVINT32X2_FROMINT64(z1));
            x1=AE_SEL32_HH(AE_MOVINT32X2_FROMINT64(z2),AE_MOVINT32X2_FROMINT64(z3));
            y0=AE_SUB32((int32_t)0xBAEC0000UL,x0);
            y1=AE_SUB32((int32_t)0xBAEC0000UL,x1);
            AE_MOVD32X4(e0,e1,0x40000000,0x40000000); 
            AE_MULSF2P32X4RAS(e0,e1,x0,x1,y0,y1); AE_MULAF2P32X4RAS(y0,y1,y0,y1,AE_ADD32S(e0,e0),AE_ADD32S(e1,e1));
            AE_MOVD32X4(e0,e1,0x40000000,0x40000000); 
            AE_MULSF2P32X4RAS(e0,e1,x0,x1,y0,y1); AE_MULAF2P32X4RAS(y0,y1,y0,y1,AE_ADD32S(e0,e0),AE_ADD32S(e1,e1));
            AE_MOVD32X4(e0,e1,0x40000000,0x40000000); 
            AE_MULSF2P32X4RAS(e0,e1,x0,x1,y0,y1); AE_MULAF2P32X4RAS(y0,y1,y0,y1,AE_ADD32S(e0,e0),AE_ADD32S(e1,e1));
            z0=AE_MOVINT64_FROMINT32X2(AE_SEL32_HH(y0,AE_ZERO32()));
            z1=AE_MOVINT64_FROMINT32X2(AE_SEL32_LH(y0,AE_ZERO32()));
            z2=AE_MOVINT64_FROMINT32X2(AE_SEL32_HH(y1,AE_ZERO32()));
            z3=AE_MOVINT64_FROMINT32X2(AE_SEL32_LH(y1,AE_ZERO32()));
            AE_S64X2_IP(z0,z1,castxcc(ae_int64x2,pScrWr),2*sizeof(ae_int64x2));

            AE_S64X2_IP(z2,z3,castxcc(ae_int64x2,pScrWr),2*sizeof(ae_int64x2));
        }
        /* last 2 iterations are done in better accuracy: 
            first one use simpler computation of 1-x*y, next 
            one adds more 8 bits for precision of this matter 
            but able to use simpler last multiply
        */
        __Pragma("no_reorder")
        pScrRd  =(const ae_int64*)scratch;
        pScrWr  =(      ae_int64*)scratch+2;
        __Pragma("loop_count factor=4,min=1")
        for (n=0; n<M; n+=2)
        {
            ae_int64   x0,y0,z0,e0;
            ae_int64   x1,y1,z1,e1;
            ae_ep      Aep0,Aep1;
            AE_L64X2_IP(x0,x1,castxcc(ae_int64x2,pScrRd),sizeof(ae_int64x2));
            AE_L64X2_IP(y0,y1,castxcc(ae_int64x2,pScrRd),sizeof(ae_int64x2));

            AE_MULZAAD32USEP_HL_LH(Aep0,e0,AE_MOVINT32X2_FROMINT64(x0),AE_MOVINT32X2_FROMINT64(y0));
            AE_MULZAAD32USEP_HL_LH(Aep1,e1,AE_MOVINT32X2_FROMINT64(x1),AE_MOVINT32X2_FROMINT64(y1));
            e0=AE_SRAI72(Aep0,e0,32);
            e1=AE_SRAI72(Aep1,e1,32);
            AE_MULA32_HH(e0,AE_MOVINT32X2_FROMINT64(x0),AE_MOVINT32X2_FROMINT64(y0));
            AE_MULA32_HH(e1,AE_MOVINT32X2_FROMINT64(x1),AE_MOVINT32X2_FROMINT64(y1));
            e0=AE_SLAI64(e0,3);
            e1=AE_SLAI64(e1,3);

            AE_MULZAAD32USEP_HL_LH(Aep0,z0,AE_MOVINT32X2_FROMINT64(y0),AE_MOVINT32X2_FROMINT64(e0));
            AE_MULZAAD32USEP_HL_LH(Aep1,z1,AE_MOVINT32X2_FROMINT64(y1),AE_MOVINT32X2_FROMINT64(e1));
            z0=AE_SRAI72(Aep0,z0,32);
            z1=AE_SRAI72(Aep1,z1,32);
            AE_MULA32_HH(z0,AE_MOVINT32X2_FROMINT64(y0),AE_MOVINT32X2_FROMINT64(e0));
            AE_MULA32_HH(z1,AE_MOVINT32X2_FROMINT64(y1),AE_MOVINT32X2_FROMINT64(e1));
            y0=AE_SUB64S(y0,z0);
            y1=AE_SUB64S(y1,z1);

            AE_S64X2_IP(y0,y1,castxcc(ae_int64x2,pScrWr),2*sizeof(ae_int64x2));
        }
        __Pragma("no_reorder")
        pScrRd  =(const ae_int64*)scratch;
        pScrWr  =(      ae_int64*)scratch;
        __Pragma("loop_count factor=4,min=1")
        for (n=0; n<M; n+=2)
        {
            ae_int64   x0,y0,z0,e0;
            ae_int64   x1,y1,z1,e1;
            ae_ep      Aep0,Aep1;
            ae_int32x2 t0,t1;

            AE_L64X2_IP(x0,x1,castxcc(ae_int64x2,pScrRd),sizeof(ae_int64x2));
            AE_L64X2_IP(y0,y1,castxcc(ae_int64x2,pScrRd),sizeof(ae_int64x2));

            z0=AE_SRLI64(AE_MUL32U_LL(AE_MOVINT32X2_FROMINT64(x0),AE_MOVINT32X2_FROMINT64(y0)),32);
            z1=AE_SRLI64(AE_MUL32U_LL(AE_MOVINT32X2_FROMINT64(x1),AE_MOVINT32X2_FROMINT64(y1)),32);
            Aep0=AE_MOVEA(0);  AE_MULAAD32USEP_HL_LH(Aep0,z0,AE_MOVINT32X2_FROMINT64(x0),AE_MOVINT32X2_FROMINT64(y0));
            Aep1=AE_MOVEA(0);  AE_MULAAD32USEP_HL_LH(Aep1,z1,AE_MOVINT32X2_FROMINT64(x1),AE_MOVINT32X2_FROMINT64(y1));
            z0=AE_SLLI64(z0,10);
            z1=AE_SLLI64(z1,10);
            e0=AE_SLLI64(AE_MUL32_HH(AE_MOVINT32X2_FROMINT64(x0),AE_MOVINT32X2_FROMINT64(y0)),42);
            e1=AE_SLLI64(AE_MUL32_HH(AE_MOVINT32X2_FROMINT64(x1),AE_MOVINT32X2_FROMINT64(y1)),42);
            e0=AE_ADD64(z0,e0);
            e1=AE_ADD64(z1,e1);

            AE_MULF2P32X4RAS(t0,t1,
                AE_MOVINT32X2_FROMINT64(y0),AE_MOVINT32X2_FROMINT64(y1),
                AE_MOVINT32X2_FROMINT64(e0),AE_MOVINT32X2_FROMINT64(e1));

            AE_S64X2_IP  (y0,y1,castxcc(ae_int64x2,pScrWr),sizeof(ae_int64x2));
            AE_S32X2X2_IP(t0,t1,castxcc(ae_int32x4,pScrWr),sizeof(ae_int32x4));
        }
        /* last stage: very accurate multiply by numerator with right sign */
        __Pragma("no_reorder")
        pScrRd  =(const ae_int64*)scratch;
        pFractRd=(const ae_int64*)frac;
        pFractWr=(      ae_int64*)frac;
        aF=AE_LA128_PP(pFractRd);
        pExp    =(      ae_int16*)exp ;
        __Pragma("loop_count factor=4,min=1")
        for (n=0; n<M; n+=2)
        {
            ae_ep Aep0,Aep1;
            ae_int64   x0,r0,z0;
            ae_int64   x1,r1,z1;
            ae_int32x2 dr0,dr1;
            ae_int16x4 e0,e1;
            xtbool sy0,sy1;
            AE_L64X2_IP(r0,r1,castxcc(ae_int64x2,pScrRd),sizeof(ae_int64x2));
            AE_L32X2X2_IP(dr0,dr1,castxcc(ae_int32x4,pScrRd),sizeof(ae_int32x4));
            {
                ae_int32x2 t0,t1;
                AE_LA32X2X2_IP(t0,t1,aF,castxcc(ae_int32x4,pFractRd));
                x0=AE_MOVINT64_FROMINT32X2(t0);
                x1=AE_MOVINT64_FROMINT32X2(t1);
            }
            Aep0=AE_MOVEA(0); z0=AE_SRLI64(AE_MUL32U_LL(AE_MOVINT32X2_FROMINT64(r0),AE_MOVINT32X2_FROMINT64(x0)),32);
            Aep1=AE_MOVEA(0); z1=AE_SRLI64(AE_MUL32U_LL(AE_MOVINT32X2_FROMINT64(r1),AE_MOVINT32X2_FROMINT64(x1)),32);
            AE_MULAAD32USEP_HL_LH(Aep0,z0,AE_MOVINT32X2_FROMINT64(r0),AE_MOVINT32X2_FROMINT64(x0));
            AE_MULAAD32USEP_HL_LH(Aep1,z1,AE_MOVINT32X2_FROMINT64(r1),AE_MOVINT32X2_FROMINT64(x1));
            z0=AE_SRAI72(Aep0,z0,31);
            z1=AE_SRAI72(Aep1,z1,31);
            AE_MULAF32S_HH(z0,AE_MOVINT32X2_FROMINT64(r0),AE_MOVINT32X2_FROMINT64(x0));
            AE_MULAF32S_HH(z1,AE_MOVINT32X2_FROMINT64(r1),AE_MOVINT32X2_FROMINT64(x1));
            z0=AE_SUB64S(z0,AE_SRAI64(AE_MUL32_HH(dr0,AE_MOVINT32X2_FROMINT64(x0)),39));
            z1=AE_SUB64S(z1,AE_SRAI64(AE_MUL32_HH(dr1,AE_MOVINT32X2_FROMINT64(x1)),39));
            e0=AE_L16_I(pExp,0*sizeof(int16_t));
            e1=AE_L16_I(pExp,1*sizeof(int16_t));
            sy0=AE_MOVBA(AE_MOVAD16_0(e0)&1);
            sy1=AE_MOVBA(AE_MOVAD16_0(e1)&1);
            AE_S16_0_IP(AE_SRAI16(e0,1),pExp,sizeof(int16_t));
            AE_S16_0_IP(AE_SRAI16(e1,1),pExp,sizeof(int16_t));
            AE_MOVT64(z0,AE_NEG64S(z0),sy0);
            AE_MOVT64(z1,AE_NEG64S(z1),sy1);
            AE_S64_IP(z0,pFractWr,sizeof(int64_t));
            AE_S64_IP(z1,pFractWr,sizeof(int64_t));
        }
    }
}

/*-------------------------------------------------------------------------
  Division
  These routines perform pair-wise division of vectors written in Q63, Q31 or 
  Q15 format. They return the fractional and exponential portion of the division 
  result. Since the division may generate result greater than 1, it returns 
  fractional portion frac in Q(63-exp), Q(31-exp) or Q(15-exp) format and 
  exponent exp so true division result in the Q0.31 may be found by shifting 
  fractional part left by exponent value.
  Additional routine makes integer division of 64-bit number to 32-bit 
  denominator forming 32-bit result. If result is overflown, 0x7fffffff 
  or 0x80000000 is returned depending on the signs of inputs.
  For division to 0, the result is not defined.

  Two versions of routines are available: regular versions (vec_divide64x32i,
  vec_divide64x64, vec_divide32x32, vec_divide16x16) work 
  with arbitrary arguments, faster versions (vec_divide32x32_fast, 
  vec_divide16x16_fast) apply some restrictions.

  Accuracy is measured as accuracy of fractional part (mantissa):
  vec_divide64x32i, scl_divide64x32                      :  1 LSB   
  vec_divide64x64                                        :  2 LSB 
  vec_divide32x32, vec_divide32x32_fast                  :  2 LSB (1.8e-9) 
  scl_divide32x32                                        :  2 LSB (4.8e-7) 
  vec_divide16x16, scl_divide16x16, vec_divide16x16_fast :  2 LSB (1.2e-4)

  Precision: 
  64x32i integer division, 64-bit nominator, 32-bit denominator, 32-bit output. 
  64x64  fractional division, 64-bit inputs, 64-bit output. 
  32x32  fractional division, 32-bit inputs, 32-bit output. 
  16x16  fractional division, 16-bit inputs, 16-bit output. 

  Input:
  x[N]    nominator, 64-bit integer, Q63, Q31 or Q15
  y[N]    denominator, 32-bit integer, Q63, Q31 or Q15
  N       length of vectors
  Output:
  frac[N] fractional parts of result, Q(63-exp), Q(31-exp) or Q(15-exp)
  exp[N]  exponents of result 

  Restriction:
  For regular versions (vec_divide64x32i, vec_divide64x64, vec_divide32x32,
  vec_divide16x16) :
  x,y,frac,exp should not overlap

  For faster versions (vec_divide32x32_fast, vec_divide16x16_fast) :
  x,y,frac,exp  should not overlap
  x,y,frac      to be aligned by 16-byte boundary, N - multiple of 4.

  Scalar versions:
  ----------------
  scl_divide64x32(): integer remainder
  Return packed value: 
  scl_divide64x64():
  bits 55:0 fractional part
  bits 63:56 exponent
  scl_divide32x32():
  bits 23:0 fractional part
  bits 31:24 exponent
  scl_divide16x16():
  bits 15:0 fractional part
  bits 31:16 exponent
-------------------------------------------------------------------------*/
void vec_divide64x64 
                (int64_t * restrict frac, 
                 int16_t *exp,
                 const int64_t * restrict x, 
                 const int64_t * restrict y, int N)
{
          ae_int64* restrict pX;
          ae_int64* restrict pY;
    const ae_int64* restrict pFractRd;
    const ae_int16* restrict pExp;
    ae_int64 ALIGN(32) scratch[SCR_SZ*2];
    int n;
    if (N<=0) return;
    if (N&7) 
    {
        int64_t ALIGN(32) xbuf[8],ybuf[8],fbuf[8];
        int16_t ebuf[8];
        pX=(ae_int64*)xbuf;
        pY=(ae_int64*)ybuf;
        AE_S64X2_I(0, 0, (ae_int64x2*)pX, 0);
        AE_S64X2_I(0, 0, (ae_int64x2*)pX, 1 * sizeof(ae_int64x2));
        AE_S64X2_I(0, 0, (ae_int64x2*)pX, 2 * sizeof(ae_int64x2));
        AE_S64X2_I(0, 0, (ae_int64x2*)pX, 3 * sizeof(ae_int64x2));
        AE_S64X2_I(0, 0, (ae_int64x2*)pY, 0);
        AE_S64X2_I(0, 0, (ae_int64x2*)pY, 1 * sizeof(ae_int64x2));
        AE_S64X2_I(0, 0, (ae_int64x2*)pY, 2 * sizeof(ae_int64x2));
        AE_S64X2_I(0, 0, (ae_int64x2*)pY, 3 * sizeof(ae_int64x2));
        __Pragma("loop_count min=1, max=7")
        for (n=0; n<(N&7); n++) 
        {
            ae_int64 t;
            AE_L64_IP(t,castxcc(ae_int64,x),sizeof(int64_t));
            AE_S64_IP(t,pX,sizeof(int64_t));
            AE_L64_IP(t,castxcc(ae_int64,y),sizeof(int64_t));
            AE_S64_IP(t,pY,sizeof(int64_t));
        }
        vec_divide64x64_internal (fbuf, ebuf,xbuf, ybuf, 8, scratch);
        pFractRd=(const ae_int64* )fbuf;
        pExp    =(const ae_int16* )ebuf;
        __Pragma("loop_count min=1, max=7")
        for (n=0; n<(N&7); n++) 
        { 
            ae_int64 t;
            ae_int16x4 q;
            AE_L64_IP(t,pFractRd,sizeof(int64_t));
            AE_S64_IP(t,castxcc(ae_int64,frac),sizeof(int64_t));
            AE_L16_IP(q,pExp,sizeof(ae_int16));
            AE_S16_0_IP(q,castxcc(ae_int16,exp),sizeof(ae_int16));
        }
        N&=~7;
    }
    if ( N) vec_divide64x64_internal(frac,exp,x,y,N,scratch);
} /* vec_divide64x64() */
