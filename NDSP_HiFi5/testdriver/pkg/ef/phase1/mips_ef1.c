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
* Test module for testing cycle performance (Vector Operations, emulated float)
*/
#include "types.h"
#include "config.h"
#include "packages.h"
#include LIBRARY_HEADER(vector)
#include "mips.h"

void mips_ef1(int isFull, int isVerbose, FILE * fout)
{
	PROFILE_NORMALIZED(1, isVerbose, vec_add_32x16ef, (mips.out0.i32, mips.out1.i16, mips.inp0.i32, mips.inp1.i16, mips.inp0.i32, mips.inp1.i16, 200), fout, "N=200", prf_cyclespts, 200);
	PROFILE_NORMALIZED(1, isVerbose, vec_mul_32x16ef, (mips.out0.i32, mips.out1.i16, mips.inp0.i32, mips.inp1.i16, mips.inp0.i32, mips.inp1.i16, 200), fout, "N=200", prf_cyclespts, 200);
	PROFILE_NORMALIZED(1, isVerbose, vec_mac_32x16ef, (mips.out0.i32, mips.out1.i16, mips.inp0.i32, mips.inp1.i16, mips.inp0.i32[0], mips.inp1.i16[0], 200), fout, "N=200", prf_cyclespts, 200);
	PROFILE_NORMALIZED(1, isVerbose, vec_dot_32x16ef, (mips.out0.i32, mips.out1.i16, mips.inp0.i32, mips.inp1.i16, mips.inp0.i32, mips.inp1.i16, 200), fout, "N=200", prf_cyclespts, 200);
	PROFILE_SIMPLE(1, isVerbose, scl_add_32x16ef, (mips.out0.i32, mips.out1.i16, mips.inp0.i32[0], mips.inp1.i16[0], mips.inp0.i32[0], mips.inp1.i16[0]), fout, "", prf_cycle);
	PROFILE_SIMPLE(1, isVerbose, scl_mul_32x16ef, (mips.out0.i32, mips.out1.i16, mips.inp0.i32[0], mips.inp1.i16[0], mips.inp0.i32[0], mips.inp1.i16[0]), fout, "", prf_cycle);
	PROFILE_SIMPLE(1, isVerbose, scl_mac_32x16ef, (mips.out0.i32, mips.out1.i16, mips.inp0.i32[0], mips.inp1.i16[0], mips.inp0.i32[0], mips.inp1.i16[0]), fout, "", prf_cycle);
}
