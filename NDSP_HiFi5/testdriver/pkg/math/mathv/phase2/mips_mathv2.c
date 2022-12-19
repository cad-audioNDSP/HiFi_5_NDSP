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

void mips_mathv2(int isFull, int isVerbose, FILE * fout)
{
    PROFILE_NORMALIZED(1, isVerbose, vec_int2float,           (mips.out0.f32, mips.inp0.i32, 1,200                   ),fout,"N=200",prf_cyclespts,200);
    PROFILE_NORMALIZED(1, isVerbose, vec_float2int,           (mips.out0.i32, mips.inp0.f32, 1,200                   ),fout,"N=200",prf_cyclespts,200);
    PROFILE_NORMALIZED(1, isVerbose, vec_sinef     ,          (mips.out0.f32, mips.inp0.f32, 200                     ),fout,"N=200",prf_cyclespts,200);
    PROFILE_NORMALIZED(1, isVerbose, vec_cosinef   ,          (mips.out0.f32, mips.inp0.f32, 200                     ),fout,"N=200",prf_cyclespts,200);
    PROFILE_NORMALIZED(1, isVerbose, vec_tanf      ,          (mips.out0.f32, mips.inp0.f32, 200                     ),fout,"N=200",prf_cyclespts,200);
    PROFILE_NORMALIZED(1, isVerbose, vec_log2f     ,          (mips.out0.f32, mips.inp0.f32, 200                     ),fout,"N=200",prf_cyclespts,200);
    PROFILE_NORMALIZED(1, isVerbose, vec_log10f    ,          (mips.out0.f32, mips.inp0.f32, 200                     ),fout,"N=200",prf_cyclespts,200);
    PROFILE_NORMALIZED(1, isVerbose, vec_lognf     ,          (mips.out0.f32, mips.inp0.f32, 200                     ),fout,"N=200",prf_cyclespts,200);
    PROFILE_NORMALIZED(1, isVerbose, vec_antilog2f ,          (mips.out0.f32, mips.inp0.f32, 200                     ),fout,"N=200",prf_cyclespts,200);
    PROFILE_NORMALIZED(1, isVerbose, vec_antilognf ,          (mips.out0.f32, mips.inp0.f32, 200                     ),fout,"N=200",prf_cyclespts,200);
    PROFILE_NORMALIZED(1, isVerbose, vec_antilog10f,          (mips.out0.f32, mips.inp0.f32, 200                     ),fout,"N=200",prf_cyclespts,200);
    PROFILE_NORMALIZED(1, isVerbose, vec_powf,                (mips.out0.f32, mips.inp0.f32,mips.inp1.f32, 200       ),fout,"N=200",prf_cyclespts,200);
    PROFILE_NORMALIZED(1, isVerbose, vec_atanf     ,          (mips.out0.f32, mips.inp0.f32, 200                     ),fout,"N=200",prf_cyclespts,200);
    PROFILE_NORMALIZED(1, isVerbose, vec_atan2f    ,          (mips.out0.f32, mips.inp0.f32,mips.inp1.f32, 200       ),fout,"N=200",prf_cyclespts,200);
    PROFILE_NORMALIZED(1, isVerbose, vec_sigmoidf,            (mips.out0.f32, mips.inp0.f32, 200                     ),fout,"N=200",prf_cyclespts,200);
    PROFILE_NORMALIZED(1, isVerbose, vec_softmaxf,            (mips.out0.f32, mips.inp0.f32, 200                     ),fout,"N=200",prf_cyclespts,200);
    PROFILE_NORMALIZED(1, isVerbose, vec_tanhf,               (mips.out0.f32, mips.inp0.f32, 200                     ),fout,"N=200",prf_cyclespts,200);
    PROFILE_NORMALIZED(1, isVerbose, vec_reluf,               (mips.out0.f32, mips.inp0.f32, 1.f, 200                ),fout,"N=200",prf_cyclespts,200);
}
