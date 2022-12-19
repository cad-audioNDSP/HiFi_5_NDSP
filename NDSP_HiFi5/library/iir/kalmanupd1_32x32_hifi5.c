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
    32-bit data, 32-bit coefficients
    C code optimized for HiFi5
  IntegrIT, 2006-2019
*/

/* Portable data types. */
#include "NatureDSP_types.h"
/* Common utility and macros declarations. */
#include "common.h"
/* Signal Processing Library API. */
#include "NatureDSP_Signal_iir.h"

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
#if 0
DISCARD_FUN(void, kalmanupd1_32x32, ( void * pScr, 
                                      int32_t * K, 
                                const int32_t * U, 
                                const int32_t * H,
                                const int32_t * R,
                                int N, int qK, int qU, int qH, int qR ));
#else
void kalmanupd1_32x32( void * pScr, 
                       int32_t * K, 
                 const int32_t * U, 
                 const int32_t * H,
                 const int32_t * R,
                 int N, int qK, int qU, int qH, int qR )
{
    NASSERT_ALIGN(K, ALIGN_SIZE);
    NASSERT_ALIGN(U, ALIGN_SIZE);
    NASSERT_ALIGN(H, ALIGN_SIZE);
    NASSERT_ALIGN(R, ALIGN_SIZE);
    NASSERT(0==(N%32));
    const ae_int32x4 * restrict pX;
    const ae_int32x4 * restrict pY;
          ae_int32x4 * restrict pZ;

    ae_int64 f, g;
    ae_int32x2 h;
    int16_t eh;
    int eU, eH, eR, ef, qf, qg, qh, sft;
    int n;
    NASSERT_ALIGN(K, ALIGN_SIZE);
    NASSERT_ALIGN(U, ALIGN_SIZE);
    NASSERT_ALIGN(H, ALIGN_SIZE);
    NASSERT_ALIGN(R, ALIGN_SIZE);
    NASSERT(0==(N%32));
    if (N<=0) return;
    /* Retrieve NSAs for all input vectors. */
    {
        ae_int16x4 eU0,eU1;
        ae_int16x4 eH0,eH1;
        eU0=31; eU1=31;
        eH0=31; eH1=31;
        pX=(const ae_int32x4 *)U;
        pY=(const ae_int32x4 *)H;
        __Pragma("loop_count min=4")
        for (n=0; n<(N>>3); n++) 
        {
            ae_int32x2 x0,x1,x2,x3,y0,y1,y2,y3;
            AE_L32X2X2_IP(x0,x1,pX,sizeof(ae_int32x4));
            AE_L32X2X2_IP(x2,x3,pX,sizeof(ae_int32x4));
            AE_L32X2X2_IP(y0,y1,pY,sizeof(ae_int32x4));
            AE_L32X2X2_IP(y2,y3,pY,sizeof(ae_int32x4));
            eU0=AE_MIN16(eU0,AE_NSA32X4(x0,x1));
            eU1=AE_MIN16(eU1,AE_NSA32X4(x2,x3));
            eH0=AE_MIN16(eH0,AE_NSA32X4(y0,y1));
            eH1=AE_MIN16(eH1,AE_NSA32X4(y2,y3));
        }
        eU=AE_RMIN16X4(AE_MIN16(eU0,eU1));
        eH=AE_RMIN16X4(AE_MIN16(eH0,eH1));
    }

    eR = NSA(R[0]);
    /* f <- dot(H,U);     */
#if 0
    {
        f = 0; 
        for ( n=0; n<N; n++ ) {
            /* Q(qH+qU+eH+eU-15) <- Q(qH+eH)*Q(qU+eU) - 15 w/ asym. rounding */
            f += ((int64_t)(H[n]<<eH)*(U[n]<<eU)+(1LL<<14))>>15;
        }
    }
#else
    {
        ae_int32x2 sH,sU;
        ae_int64 a0=0,a1=0;
        pX=(const ae_int32x4 *)H;
        pY=(const ae_int32x4 *)U;
        eH=XT_MIN(eH,29);
        eU=XT_MIN(eU,29);
        sH=AE_SLAA32(1,eH);
        sU=AE_SLAA32(1,eU);
        __Pragma("loop_count min=4")
        for ( n=0; n<(N>>3); n++ ) 
        {
            ae_int32x2 x0,x1,y0,y1;
            AE_L32X2X2_IP(x0,x1,pX,sizeof(ae_int32x4));
            AE_L32X2X2_IP(y0,y1,pY,sizeof(ae_int32x4));
            AE_MUL2P32X4(x0,x1,x0,x1,sH,sH);
            AE_MUL2P32X4(y0,y1,y0,y1,sU,sU);
            AE_MULAAF2D32RA_HH_LL(a0,a1,x0,x1,y0,y1);

            AE_L32X2X2_IP(x0,x1,pX,sizeof(ae_int32x4));
            AE_L32X2X2_IP(y0,y1,pY,sizeof(ae_int32x4));
            AE_MUL2P32X4(x0,x1,x0,x1,sH,sH);
            AE_MUL2P32X4(y0,y1,y0,y1,sU,sU);
            AE_MULAAF2D32RA_HH_LL(a0,a1,x0,x1,y0,y1);
        }
        f=AE_ADD64(a0,a1);
    }
#endif
    ef = AE_NSA64(f);
    qf = qH+qU+eH+eU-15+ef;
    f = AE_SLAA64(f,ef); /* f: normal 64-bit Q(qf) */
    /*
     * g <- f+R; g: 64-bit Q(qg)
     * Fixed point position for the sum is selected from magnitudes of both summands.
     * Note that f and R should be non-negative from the essence of Kalman filter.
     */
    {
        xtbool qfbig;
        qfbig = AE_MOVBA (qf>(qR+eR+32)?1:0);
        ae_int64 gbig;
        gbig=AE_ADD64( AE_SRAA64(f,(qf-(qR+eR+32))) , AE_SLAA64(R[0],(eR+32)));
        g = AE_ADD64(f,AE_SLAA64S(R[0], qf-qR));
        AE_MOVT64(g,gbig,qfbig);
    }
    qg = XT_MIN(qf,qR+eR+32);
    /* Re-normalize the sum when it overflows into the sign bit. */
    {
        xtbool blezero;
        blezero=AE_LT64(g,0);
        AE_MOVT64(g,AE_SRLI64(g,1),blezero);
        qg=qg-AE_MOVAB(blezero);
    }
    /* h <- 1/g;     */ 
#if 0
    /* Q(118-qg-eh) <- Q(55+63)/[Q(qg)+eh] */
    w = (ae_int64)scl_recip64x64((int64_t)g); 
    eh = (int16_t)((int64_t)AE_SRAI64(w,56)); w=AE_SLAI64(w,8); /* Normal Q(126-qg-eh) */
    /* Q(94-qg-eh) <- Q(126-qg-eh) - 32 w/ asym. rounding */
    h = AE_ROUND32X2F64SASYM(w,w);
    qh = 94-qg-eh; /* h: 32-bit normal Q(qh) */
#else
    /* Q(118-qg-eh) <- Q(55+63)/[Q(qg)+eh] */
    {
        ae_int64  lx=g;
        ae_f32x2 xx,yy,ee;
        xtbool sx,bzero=AE_EQ64(lx,AE_ZERO64());
        ae_int64  y,e;
        ae_ep Aep; 
        ae_int64 A,B;
        int exp;
        sx=AE_LT64(lx,AE_ZERO64());
        lx=AE_ABS64S(lx);
        exp=AE_NSA64(lx);
        lx=AE_SLAA64(lx,exp);// x in 0.5..1
        exp++;
        /* first stage use 3 iterations computing in 32-bit accuracy */
        xx=AE_TRUNCI32F64S(lx,0);
        yy=AE_SUB32((int32_t)0xBAEC0000UL,xx);
        ee=0x40000000; AE_MULSFP32X2RAS(ee,xx,yy); AE_MULAFP32X2RAS(yy,yy,AE_SLLI32(ee,1));
        ee=0x40000000; AE_MULSFP32X2RAS(ee,xx,yy); AE_MULAFP32X2RAS(yy,yy,AE_SLLI32(ee,1));
        ee=0x40000000; AE_MULSFP32X2RAS(ee,xx,yy); AE_MULAFP32X2RAS(yy,yy,AE_SLLI32(ee,1));
        y=AE_MOVINT64_FROMF32X2(AE_SEL32_LL(yy,0));
        /* last 1 iteration is done in better accuracy */
        AE_MULZAAD32USEP_HL_LH(Aep,A,AE_MOVINT32X2_FROMINT64(lx),AE_MOVINT32X2_FROMINT64(y));
        A=AE_SRAI72(Aep,A,29);
        B=AE_MUL32_HH(AE_MOVINT32X2_FROMINT64(lx),AE_MOVINT32X2_FROMINT64(y));
        B=AE_SLLI64(B,3);
        e=AE_ADD64(B,A);

        AE_MULZAAD32USEP_HL_LH(Aep,A,AE_MOVINT32X2_FROMINT64(y),AE_MOVINT32X2_FROMINT64(e));
        A=AE_SRAI72(Aep,A,32);
        AE_MULA32_HH(A,AE_MOVINT32X2_FROMINT64(y),AE_MOVINT32X2_FROMINT64(e));
        y=AE_SUB64S(y,A);

        AE_MOVT64(y,AE_NEG64(y),sx);
        AE_MOVT64(y,0x7fffffffffffffffUL,bzero);
        h = AE_ROUND32X2F64SASYM(y,y);
        eh=exp;
    }
    /* Q(94-qg-eh) <- Q(126-qg-eh) - 32 w/ asym. rounding */
    qh = 94-qg-eh; /* h: 32-bit normal Q(qh) */
#endif
    /*
     * K <- h*U
     */
#if 0
    {
        int32_t r;
        sft = qK-qU-qh;
        r = sft<0 ? 1LL<<-(sft+1) : 0;
        for ( n=0; n<N; n++ ) 
        {
        /* Q(qK) = Q(qh)*Q(qU) + qK-qU-qh w/ asym. rounding */
            K[n] = L_sature_64(shls_64(add_64((int64_t)h*U[n], r), sft));
        }
    }
#else
    sft = XT_MIN(qK-qU-qh+32,63);
    pX=(const ae_int32x4 *)U;
    pZ=(      ae_int32x4 *)K;
    __Pragma("loop_count min=4")
    for ( n=0; n<(N>>3); n++ ) 
    {
        ae_int64 a0,a1,a2,a3,a4,a5,a6,a7;
        ae_int32x2 x0,x1,x2,x3;
        AE_L32X2X2_IP(x0,x1,pX,sizeof(ae_int32x4));
        AE_L32X2X2_IP(x2,x3,pX,sizeof(ae_int32x4));
        AE_MUL32X2S_HH_LL(a0,a1,x0,h);
        AE_MUL32X2S_HH_LL(a2,a3,x1,h);
        AE_MUL32X2S_HH_LL(a4,a5,x2,h);
        AE_MUL32X2S_HH_LL(a6,a7,x3,h);
        x0=AE_TRUNCA32X2F64S(a0,a1,sft);
        x1=AE_TRUNCA32X2F64S(a2,a3,sft);
        x2=AE_TRUNCA32X2F64S(a4,a5,sft);
        x3=AE_TRUNCA32X2F64S(a6,a7,sft);
        AE_S32X2X2_IP(x0,x1,pZ,sizeof(ae_int32x4));
        AE_S32X2X2_IP(x2,x3,pZ,sizeof(ae_int32x4));
    }
#endif
} /* kalmanupd1_32x32() */
#endif

/* Returns: size of scratch memory area, in bytes. */
size_t kalmanupd1_32x32_getScratchSize(int N)
{
    NASSERT(0==(N%32));
    (void)N;
    return 0;
}
