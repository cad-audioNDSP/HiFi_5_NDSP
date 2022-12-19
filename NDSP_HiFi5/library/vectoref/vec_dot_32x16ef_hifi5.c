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
    Vector dot product for Emulated Floating Point
    optimized code for HiFi5 core
*/
/* Library API */
#include "NatureDSP_Signal_vector.h"
#include "common.h"
// code optimized for HiFi5 core

/*-------------------------------------------------------------------------
  Vector Dot Product for Emulated Floating Point
  routines compute dot product of vectors represented in emulated floating 
  point format

  Input:
  xmant[N]  mantissa of input data
  ymant[N]  mantissa of input data
  xexp[N]   exponent of input data
  yexp[N]   exponent of input data
  N         length of vectors
  Output:
  zmant[1]  mantissa of output data
  zexp[1]   exponent of output data

  Restriction:
  xmant,ymant,xexp,yexp,zmant,zexp should not overlap
-------------------------------------------------------------------------*/
void   vec_dot_32x16ef (      int32_t  * restrict zmant,       int16_t  * restrict zexp, 
                        const int32_t  * restrict xmant, const int16_t  * restrict xexp, 
                        const int32_t  * restrict ymant, const int16_t  * restrict yexp, 
                        int N)
{
    ae_valignx2 aXexp,aXmant,aYexp,aYmant;
    const ae_int16x8 * restrict pXexp;
    const ae_int16x8 * restrict pYexp;
    const ae_int32x4 * restrict pXmant;
    const ae_int32x4 * restrict pYmant;
    ae_int64 a;
    int maxexp;
    int n;
    zmant[0]=0; zexp[0]=0;
    if (N<=0) return;
    /* compute common exponent taking into account zero mantissa words */
    pXexp =(const ae_int16x8 *)xexp ;
    pYexp =(const ae_int16x8 *)yexp ;
    pXmant=(const ae_int32x4 *)xmant;
    pYmant=(const ae_int32x4 *)ymant;
    if (N<=7)
    {
        ae_int16x4  vmaxexp=MIN_INT16;
        for(n=0; n<N; n++) 
        {
            ae_int32x2 xm0,ym0;
            ae_int16x4 xe0,ye0;
            AE_L32_IP(xm0,castxcc(ae_int32,pXmant),sizeof(int32_t));
            AE_L32_IP(ym0,castxcc(ae_int32,pYmant),sizeof(int32_t));
            AE_L16_IP(xe0,castxcc(ae_int16,pXexp ),sizeof(int16_t)) ;
            AE_L16_IP(ye0,castxcc(ae_int16,pYexp ),sizeof(int16_t)) ;
            AE_MOVT16X4(xe0,MIN_INT16,AE_EQ16(AE_SAT16X4(xm0,xm0),AE_ZERO16()));
            AE_MOVT16X4(ye0,MIN_INT16,AE_EQ16(AE_SAT16X4(ym0,ym0),AE_ZERO16()));
            vmaxexp=AE_MAX16(vmaxexp,AE_ADD16S(xe0,ye0));
        }
        maxexp=AE_RMAX16X4(vmaxexp);
    }
    else
    {
        ae_int32x2 xm0,xm1,xm2,xm3;
        ae_int32x2 ym0,ym1,ym2,ym3;
        ae_int16x4 xe0,xe1;
        ae_int16x4 ye0,ye1;
        ae_int16x4 vmaxexp=MIN_INT16;
        xtbool4 bxmzero0,bymzero0;
        xtbool4 bxmzero1,bymzero1;
        int N0;
        N0=((N-1)&7)+1;
        pXmant=(ae_int32x4 *)xmant;
        pYmant=(ae_int32x4 *)ymant;
        pXexp =(ae_int16x8 *)xexp ;
        pYexp =(ae_int16x8 *)yexp ;
        aXmant=AE_LA128_PP(pXmant);
        aYmant=AE_LA128_PP(pYmant);
        aXexp =AE_LA128_PP(pXexp );
        aYexp =AE_LA128_PP(pYexp );
        AE_LA32X2X2_IP(xm0,xm1,aXmant,pXmant);
        AE_LA32X2X2_IP(xm2,xm3,aXmant,pXmant);
        AE_LA32X2X2_IP(ym0,ym1,aYmant,pYmant);
        AE_LA32X2X2_IP(ym2,ym3,aYmant,pYmant);
        AE_LA16X4X2_IP(xe0,xe1,aXexp ,pXexp);
        AE_LA16X4X2_IP(ye0,ye1,aYexp ,pYexp);

        bxmzero0=AE_EQ16(AE_SAT16X4(xm0,xm1),AE_ZERO16());
        bxmzero1=AE_EQ16(AE_SAT16X4(xm2,xm3),AE_ZERO16());
        bymzero0=AE_EQ16(AE_SAT16X4(ym0,ym1),AE_ZERO16());
        bymzero1=AE_EQ16(AE_SAT16X4(ym2,ym3),AE_ZERO16());
        AE_MOVT16X4(xe0,MIN_INT16,bxmzero0);
        AE_MOVT16X4(xe1,MIN_INT16,bxmzero1);
        AE_MOVT16X4(ye0,MIN_INT16,bymzero0);
        AE_MOVT16X4(ye1,MIN_INT16,bymzero1);
        vmaxexp=AE_MAX16(vmaxexp,AE_ADD16S(xe0,ye0));
        vmaxexp=AE_MAX16(vmaxexp,AE_ADD16S(xe1,ye1));
        pXmant=(ae_int32x4 *)(xmant+N0);
        pYmant=(ae_int32x4 *)(ymant+N0);
        pXexp =(ae_int16x8 *)(xexp +N0);
        pYexp =(ae_int16x8 *)(yexp +N0);
        aXmant=AE_LA128_PP(pXmant);
        aYmant=AE_LA128_PP(pYmant);
        aXexp =AE_LA128_PP(pXexp );
        aYexp =AE_LA128_PP(pYexp );
        for (n=0; n<((N-N0)>>3); n++) 
        {
            AE_LA32X2X2_IP(xm0,xm1,aXmant,pXmant);
            AE_LA32X2X2_IP(xm2,xm3,aXmant,pXmant);
            AE_LA32X2X2_IP(ym0,ym1,aYmant,pYmant);
            AE_LA32X2X2_IP(ym2,ym3,aYmant,pYmant);
            AE_LA16X4X2_IP(xe0,xe1,aXexp ,pXexp);
            AE_LA16X4X2_IP(ye0,ye1,aYexp ,pYexp);
#if 0
            bxmzero0=AE_EQ16(AE_SAT16X4(xm0,xm1),AE_ZERO16());
            bxmzero1=AE_EQ16(AE_SAT16X4(xm2,xm3),AE_ZERO16());
            bymzero0=AE_EQ16(AE_SAT16X4(ym0,ym1),AE_ZERO16());
            bymzero1=AE_EQ16(AE_SAT16X4(ym2,ym3),AE_ZERO16());
            AE_MOVT16X4(xe0,MIN_INT16,bxmzero0);
            AE_MOVT16X4(xe1,MIN_INT16,bxmzero1);
            AE_MOVT16X4(ye0,MIN_INT16,bymzero0);
            AE_MOVT16X4(ye1,MIN_INT16,bymzero1);
#else
        AE_EXPADD16_H(xe0, 0, xm0);
        AE_EXPADD16_L(xe0, 0, xm1);
        AE_EXPADD16_H(xe1, 0, xm2);
        AE_EXPADD16_L(xe1, 0, xm3);
        AE_EXPADD16_H(ye0, 0, ym0);
        AE_EXPADD16_L(ye0, 0, ym1);
        AE_EXPADD16_H(ye1, 0, ym2);
        AE_EXPADD16_L(ye1, 0, ym3);
#endif
            vmaxexp=AE_MAX16(vmaxexp,AE_ADD16S(xe0,ye0));
            vmaxexp=AE_MAX16(vmaxexp,AE_ADD16S(xe1,ye1));
        }
        maxexp=AE_RMAX16X4(vmaxexp);
    }
    /* accumulate products with proper right shift */
    a=AE_ZERO64();
    pXmant=(ae_int32x4 *)xmant;
    pYmant=(ae_int32x4 *)ymant;
    pXexp =(ae_int16x8 *)xexp ;
    pYexp =(ae_int16x8 *)yexp ;
    __Pragma("loop_count min=0,max=3")
    __Pragma("no_unroll")
    for(n=0; n<(N&3); n++) 
    {
        ae_int64 t;
        ae_int16x4 xe0,ye0;
        ae_int32x2 xm0,ym0;
        int sh;
        AE_L16_IP(xe0,castxcc(ae_int16,pXexp),sizeof(int16_t));
        AE_L16_IP(ye0,castxcc(ae_int16,pYexp),sizeof(int16_t));
        AE_L32_IP(xm0,castxcc(ae_int32,pXmant),sizeof(int32_t));
        AE_L32_IP(ym0,castxcc(ae_int32,pYmant),sizeof(int32_t));
        sh=AE_MOVAD16_0(AE_SUB16S(AE_ADD16S(xe0,ye0),(ae_int16x4)(int16_t)(maxexp)));
        t=AE_MULF32RA_HH(xm0,ym0);    /* accumulate in Q47 */
        a=AE_ADD64S(a,AE_SLAA64S(t,sh));
    }
    N&=~3;
    /* process the most data block by 4-element portions */
    ae_int64 a0=AE_ZERO64();
    {
        ae_valign aXexp,aYexp;
        aXexp =AE_LA64_PP(castxcc(ae_int16x4,pXexp ));
        aYexp =AE_LA64_PP(castxcc(ae_int16x4,pYexp ));
        aXmant=AE_LA128_PP(castxcc(ae_int32x4,pXmant));
        aYmant=AE_LA128_PP(castxcc(ae_int32x2,pYmant));
        for(n=0; n<(N>>2); n++) 
        {
            ae_int16x4 sh,xsh,ysh;
            ae_int32x2 x0,x1,y0,y1;
            ae_int64 t0,t1,t2,t3;
            AE_LA16X4_IP(xsh,aXexp,castxcc(ae_int16x4,pXexp ));
            AE_LA16X4_IP(ysh,aYexp,castxcc(ae_int16x4,pYexp ));
            sh=AE_SUB16S(AE_ADD16S(xsh,ysh),(ae_int16x4)(int16_t)maxexp);
            AE_LA32X2X2_IP(x0,x1,aXmant,pXmant);
            AE_LA32X2X2_IP(y0,y1,aYmant,pYmant);
            AE_MULF32X2R_HH_LL(t0,t1,x0,y0);
            AE_MULF32X2R_HH_LL(t2,t3,x1,y1);
            t0=AE_SLAA64S(t0,AE_MOVAD16_3(sh));
            t1=AE_SLAA64S(t1,AE_MOVAD16_2(sh));
            t2=AE_SLAA64S(t2,AE_MOVAD16_1(sh));
            t3=AE_SLAA64S(t3,AE_MOVAD16_0(sh));
            a =AE_ADD64S(a,t0);
            a =AE_ADD64S(a,t1);
            a =AE_ADD64S(a,t2);
            a =AE_ADD64S(a,t3);
        }
    }
    a=AE_ADD64S(a,a0);


    /* normalize Q47 result */
    {
        int nsa0;
        ae_int32x2 nsa,x;
        nsa0=AE_NSA64(a);
        x=AE_TRUNCA32X2F64S(a,a,nsa0);
        nsa=AE_MOVDA32(maxexp-nsa0+16);
        /* underflow processing */
        AE_MOVT32X2(x,0,AE_LT32(nsa,MIN_INT16));
        /* overflow processing */
        AE_MOVT32X2(x,AE_SLLI32S(x,31),AE_LT32(MAX_INT16,nsa));
        AE_S32_L_IP(x,castxcc(ae_int32,zmant),sizeof(int32_t));
        AE_S16_0_IP(AE_SAT16X4(nsa,nsa),castxcc(ae_int16,zexp),sizeof(int16_t));
    }
}
