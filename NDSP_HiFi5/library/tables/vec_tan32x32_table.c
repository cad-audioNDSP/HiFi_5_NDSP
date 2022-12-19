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
#include "vec_tan32x32_table.h"
#include "common.h"

/*
x=(-1:1e-5:1);
s=pow2(1,-16); x=(s:s:1/4); x=[-x(end:-1:1) x];
[pt]=polyfit(x,tan(pi*x)./x,9);
[pq]=polyfit(x,cot(pi*x).*x,10);
fprintf(1,'pt=');
fprintf(1,'%d ',round(pow2(pt(1:end),20)));
fprintf(1,'\n');
fprintf(1,'pq=');
fprintf(1,'%d ',round(pow2(pq(1:end),24)));
fprintf(1,'\n');
*/
const int32_t ALIGN(32) Pt[] = { 1321653322, 133294643, 43623826, 10831170, 3294206 };//Q20
const int32_t ALIGN(32) Qt[] = { -12847376, -10574258, -10870638, -11559912, -17569060, 5340354 };//Q24
