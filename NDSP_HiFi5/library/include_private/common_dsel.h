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

#ifndef __COMMON_DSEL_H

#include <common.h>

/* Emulate AE_DSEL32X2 by means of AE_DSEL16X4. Input argument e specifies the
 * selection pattern, analogous to AE_DSEL16X4: lower 16-bit word in each 32-bit
 * lane holds a 2-bit index for the respective lane of the output argument b, higher
 * 16-bit word of the same lane contains a 2-bit index for the output argument a. */
#if !defined(AE_DSEL32X2) && defined(AE_DSEL16X4)
#define AE_DSEL32X2(a,b,c,d,e) __AE_DSEL32X2(&a, &b, &c, &d, &e)
inline_ void __AE_DSEL32X2(ae_int32x2 * a, ae_int32x2 * b, const ae_int32x2 * cc, const ae_int32x2 * d, const ae_int32x2 * e)
{
    static const int16_t ALIGN(16) ishfl[3*4] = { 2,2,0,0, 3,3,1,1, 257,0,257,0 };
    ae_int16x4 vs, vse, vso, ofs;
    ae_int16x4 _a, _b, _c, _d, _e;
    _c = AE_MOVINT16X4_FROMINT32X2(*cc);
    _d = AE_MOVINT16X4_FROMINT32X2(*d);
    _e = AE_MOVINT16X4_FROMINT32X2(*e);
    AE_L16X4X2_I(vse, vso, (ae_int16x8*)ishfl, 0);
    ofs = AE_L16X4_I((ae_int16x4*)ishfl, 8*sizeof(int16_t));
    _e = AE_SLAI16(_e, 1);
    vse = AE_SHFL16X4(_e, vse);
    vso = AE_SLAI16(AE_SHFL16X4(_e, vso), 8);
    vs = AE_ADD16(AE_ADD16(vse, vso), ofs);
    AE_DSEL16X4(_a, _b, _c, _d, vs);
    *a = AE_MOVINT32X2_FROMINT16X4(_a);
    *b = AE_MOVINT32X2_FROMF16X4(_b);
}
#endif

#endif /* __COMMON_DSEL_H */
