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
* Test module for testing cycle performance (Vectorized Math)
*/

#include "mips.h"
#include "config.h"
#include "packages.h"
#include LIBRARY_HEADER(math)

void mips_mathv1(int isFull, int isVerbose, FILE * fout)
{
	PROFILE_NORMALIZED(     1, isVerbose, vec_recip16x16,           (mips.out0.i16,mips.out1.i16,mips.inp0.i16,200              ),fout,"N=200",prf_cyclespts,200);
#if 0 //HiFi3/3z API
	PROFILE_NORMALIZED(     1, isVerbose, vec_recip24x24,           (mips.out0.i32,mips.out1.i16,mips.inp0.i32,200              ),fout,"N=200",prf_cyclespts,200);
#endif
	PROFILE_NORMALIZED(     1, isVerbose, vec_recip32x32,           (mips.out0.i32, mips.out1.i16, mips.inp0.i32, 200           ),fout,"N=200",prf_cyclespts,200);
	PROFILE_NORMALIZED(     1, isVerbose, vec_recip64x64,           (mips.out0.i64, mips.out1.i16, mips.inp0.i64, 200           ),fout,"N=200",prf_cyclespts,200);
	PROFILE_NORMALIZED(isFull, isVerbose, vec_divide16x16,          (mips.out0.i16,mips.out1.i16,mips.inp0.i16,mips.inp1.i16,200     ),fout,"N=200",prf_cyclespts,200);
#if 0 //HiFi3/3z API
	PROFILE_NORMALIZED(isFull, isVerbose, vec_divide24x24,          (mips.out0.i32, mips.out1.i16, mips.inp0.i32, mips.inp1.i32, 200 ),fout,"N=200",prf_cyclespts,200);
#endif
	PROFILE_NORMALIZED(isFull, isVerbose, vec_divide32x32,          (mips.out0.i32,mips.out1.i16,mips.inp0.i32,mips.inp1.i32,200     ),fout,"N=200",prf_cyclespts,200);
	PROFILE_NORMALIZED(     1, isVerbose, vec_divide64x32i,         (mips.out0.i32,         mips.inp0.i64,mips.inp1.i32,200     ),fout,"N=200",prf_cyclespts,200);
	PROFILE_NORMALIZED(     1, isVerbose, vec_divide64x64,          (mips.out0.i64,mips.out1.i16,mips.inp0.i64,mips.inp1.i64,200     ),fout,"N=200",prf_cyclespts,200);
	PROFILE_NORMALIZED(     1, isVerbose, vec_log2_32x32,           (mips.out0.i32, mips.inp0.i32, 200                     ),fout,"N=200",prf_cyclespts,200);
	PROFILE_NORMALIZED(     1, isVerbose, vec_logn_32x32,           (mips.out0.i32, mips.inp0.i32, 200                     ),fout,"N=200",prf_cyclespts,200);
	PROFILE_NORMALIZED(     1, isVerbose, vec_log10_32x32,          (mips.out0.i32, mips.inp0.i32, 200                     ),fout,"N=200",prf_cyclespts,200);
#if 0 //HiFi3/3z API
	PROFILE_NORMALIZED(     1, isVerbose, vec_log2_24x24,           (mips.out0.i32, mips.inp0.i32, 200                     ),fout,"N=200",prf_cyclespts,200);
	PROFILE_NORMALIZED(     1, isVerbose, vec_logn_24x24,           (mips.out0.i32, mips.inp0.i32, 200                     ),fout,"N=200",prf_cyclespts,200);
	PROFILE_NORMALIZED(     1, isVerbose, vec_log10_24x24,          (mips.out0.i32, mips.inp0.i32, 200                     ),fout,"N=200",prf_cyclespts,200);
	PROFILE_NORMALIZED(     1, isVerbose, vec_antilog2_24x24,       (mips.out0.i32, mips.inp0.i32, 200                     ),fout,"N=200",prf_cyclespts,200);
	PROFILE_NORMALIZED(     1, isVerbose, vec_antilogn_24x24,       (mips.out0.i32, mips.inp0.i32, 200                     ),fout,"N=200",prf_cyclespts,200);
	PROFILE_NORMALIZED(     1, isVerbose, vec_antilog10_24x24,      (mips.out0.i32, mips.inp0.i32, 200                     ),fout,"N=200",prf_cyclespts,200);
#endif
	PROFILE_NORMALIZED(     1, isVerbose, vec_antilog2_32x32,       (mips.out0.i32, mips.inp0.i32, 200                     ),fout,"N=200",prf_cyclespts,200);
	PROFILE_NORMALIZED(     1, isVerbose, vec_antilogn_32x32,       (mips.out0.i32, mips.inp0.i32, 200                     ),fout,"N=200",prf_cyclespts,200);
	PROFILE_NORMALIZED(     1, isVerbose, vec_antilog10_32x32,      (mips.out0.i32, mips.inp0.i32, 200                     ),fout,"N=200",prf_cyclespts,200);
	PROFILE_NORMALIZED(     1, isVerbose, vec_pow_32x32,            (mips.out0.i32, mips.out0.i16,mips.inp0.i32,mips.inp1.i32, 200   ),fout,"N=200",prf_cyclespts,200);
	PROFILE_NORMALIZED(isFull, isVerbose, vec_sine32x32,            (mips.out0.i32, mips.inp0.i32, 200                     ),fout,"N=200",prf_cyclespts,200);
	PROFILE_NORMALIZED(isFull, isVerbose, vec_cosine32x32,          (mips.out0.i32, mips.inp0.i32, 200                     ),fout,"N=200",prf_cyclespts,200);
#if 0 //HiFi3/3z API
	PROFILE_NORMALIZED(isFull, isVerbose, vec_sine24x24,            (mips.out0.i32, mips.inp0.i32, 200                     ),fout,"N=200",prf_cyclespts,200);
	PROFILE_NORMALIZED(isFull, isVerbose, vec_cosine24x24,          (mips.out0.i32, mips.inp0.i32, 200                     ),fout,"N=200",prf_cyclespts,200);
#endif
	PROFILE_NORMALIZED(     1, isVerbose, vec_tan32x32,             (mips.out0.i32, mips.inp0.i32, 200                     ),fout,"N=200",prf_cyclespts,200);
#if 0 //HiFi3/3z API
	PROFILE_NORMALIZED(     1, isVerbose, vec_tan24x24,             (mips.out0.i32, mips.inp0.i32, 200                     ),fout,"N=200",prf_cyclespts,200);
#endif
	PROFILE_NORMALIZED(     1, isVerbose, vec_atan32x32,            (mips.out0.i32, mips.inp0.i32, 200                     ),fout,"N=200",prf_cyclespts,200);
#if 0 //HiFi3/3z API
	PROFILE_NORMALIZED(     1, isVerbose, vec_atan24x24,            (mips.out0.i32, mips.inp0.i32, 200                     ),fout,"N=200",prf_cyclespts,200);
	PROFILE_NORMALIZED(     1, isVerbose, vec_atan2_24x24,          (mips.out0.i32, mips.inp0.i32, mips.inp1.i32, 200           ),fout,"N=200",prf_cyclespts,200);
#endif
	PROFILE_NORMALIZED(     1, isVerbose, vec_sqrt16x16,            (mips.out0.i16, mips.inp0.i16, 200                     ),fout,"N=200",prf_cyclespts,200);
#if 0 //HiFi3/3z API
	PROFILE_NORMALIZED(isFull, isVerbose, vec_sqrt24x24,            (mips.out0.i32, mips.inp0.i32, 200                     ),fout,"N=200",prf_cyclespts,200);
#endif
	PROFILE_NORMALIZED(     1, isVerbose, vec_sqrt32x16,            (mips.out0.i16, mips.inp0.i32, 200                     ),fout,"N=200",prf_cyclespts,200);
	PROFILE_NORMALIZED(isFull, isVerbose, vec_sqrt32x32,            (mips.out0.i32, mips.inp0.i32, 200                     ),fout,"N=200",prf_cyclespts,200);
	PROFILE_NORMALIZED(     1, isVerbose, vec_sqrt64x32,            (mips.out0.i32, mips.inp0.i64, 200                     ),fout,"N=200",prf_cyclespts,200);
	PROFILE_NORMALIZED(     1, isVerbose, vec_rsqrt16x16,           (mips.out0.i16, mips.out1.i16, mips.inp0.i16, 200           ),fout,"N=200",prf_cyclespts,200);
	PROFILE_NORMALIZED(     1, isVerbose, vec_rsqrt32x32,           (mips.out0.i32, mips.out1.i16, mips.inp0.i32, 200           ),fout,"N=200",prf_cyclespts,200);
	PROFILE_NORMALIZED(     1, isVerbose, vec_sigmoid32x32,         (mips.out0.i32, mips.inp0.i32, 200                     ),fout,"N=200",prf_cyclespts,200);
	PROFILE_NORMALIZED(     1, isVerbose, vec_softmax32x32,         (mips.out0.i32, mips.inp0.i32, 200                     ),fout,"N=200",prf_cyclespts,200);
	PROFILE_NORMALIZED(     1, isVerbose, vec_tanh32x32,            (mips.out0.i32, mips.inp0.i32, 200                     ),fout,"N=200",prf_cyclespts,200);
	PROFILE_NORMALIZED(     1, isVerbose, vec_relu32x32,            (mips.out0.i32, mips.inp0.i32, 10000,200               ),fout,"N=200",prf_cyclespts,200);
}
