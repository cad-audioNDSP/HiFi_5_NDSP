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
#include "NatureDSP_Signal_fft.h"
#include "common.h"
#include "fft_x16_common.h"
#include "fft_16x16_stages.h"

/* Twiddles tables for fft_cplx32x16, ifft_cplx32x16, fft_cplx16x16, ifft_cplx16x16, N=48 */
#define N 48

/********** Twiddles table N=48 stage 1 radix 3 ******************/
ALIGN(32) static const int16_t _fft48_tw1_[] =
{
    (int16_t)0x7FFF, (int16_t)0x0000, (int16_t)0x7FFF, (int16_t)0x0000, (int16_t)0x7EE8, (int16_t)0xEF4B, (int16_t)0x7BA3, (int16_t)0xDEDF,
    (int16_t)0x7BA3, (int16_t)0xDEDF, (int16_t)0x6EDA, (int16_t)0xC000, (int16_t)0x7642, (int16_t)0xCF04, (int16_t)0x5A82, (int16_t)0xA57E,
    (int16_t)0x6EDA, (int16_t)0xC000, (int16_t)0x4000, (int16_t)0x9126, (int16_t)0x658D, (int16_t)0xB214, (int16_t)0x2121, (int16_t)0x845D,
    (int16_t)0x5A82, (int16_t)0xA57E, (int16_t)0x0000, (int16_t)0x8000, (int16_t)0x4DEC, (int16_t)0x9A73, (int16_t)0xDEDF, (int16_t)0x845D,
    (int16_t)0x4000, (int16_t)0x9126, (int16_t)0xC000, (int16_t)0x9126, (int16_t)0x30FC, (int16_t)0x89BE, (int16_t)0xA57E, (int16_t)0xA57E,
    (int16_t)0x2121, (int16_t)0x845D, (int16_t)0x9126, (int16_t)0xC000, (int16_t)0x10B5, (int16_t)0x8118, (int16_t)0x845D, (int16_t)0xDEDF,
    (int16_t)0x0000, (int16_t)0x8000, (int16_t)0x8000, (int16_t)0x0000, (int16_t)0xEF4B, (int16_t)0x8118, (int16_t)0x845D, (int16_t)0x2121,
    (int16_t)0xDEDF, (int16_t)0x845D, (int16_t)0x9126, (int16_t)0x4000, (int16_t)0xCF04, (int16_t)0x89BE, (int16_t)0xA57E, (int16_t)0x5A82,
};

/********** Twiddles table N=48 stage 2 radix 4 ******************/
ALIGN(32) static const int16_t _fft48_tw2_[] =
{
    (int16_t)0x7FFF, (int16_t)0x0000, (int16_t)0x7FFF, (int16_t)0x0000, (int16_t)0x7FFF, (int16_t)0x0000, (int16_t)0x7642, (int16_t)0xCF04,
    (int16_t)0x5A82, (int16_t)0xA57E, (int16_t)0x30FC, (int16_t)0x89BE, (int16_t)0x5A82, (int16_t)0xA57E, (int16_t)0x0000, (int16_t)0x8000,
    (int16_t)0xA57E, (int16_t)0xA57E, (int16_t)0x30FC, (int16_t)0x89BE, (int16_t)0xA57E, (int16_t)0xA57E, (int16_t)0x89BE, (int16_t)0x30FC, 
};

static const int tw_step_tab[] =
{
    1, 1, 0
}; 
static const cint16ptr_t_fft tw_tab[] = 
{
    _fft48_tw1_, _fft48_tw2_, NULL
};

static const fn_fft_stage fft_stg_tab_s2[] = 
{
    (fn_fft_stage)fft_16x16_stage_first_scl2_DFT3,
    (fn_fft_stage)fft_16x16_stage_inner_scl2_DFT4_v3,
    (fn_fft_stage)fft_16x16_stage_last_scl2_DFT4
};
//static const fn_fft_stage fft_stg_tab_s3[] = 
//{
//    (fn_fft_stage)NULL, //(fn_fft_stage)fft_16x16_stage_first_scl3_DFT3,
//    (fn_fft_stage)NULL, // not implemented yet fft_16x16_stage_inner_scl3_DFT4,
//    (fn_fft_stage)NULL //(fn_fft_stage)fft_16x16_stage_last_scl3_DFT4
//};
static const fn_fft_stage ifft_stg_tab_s2[] =
{
    (fn_fft_stage)ifft_16x16_stage_first_scl2_DFT3,
    (fn_fft_stage)fft_16x16_stage_inner_scl2_DFT4_v3,
    (fn_fft_stage)fft_16x16_stage_last_scl2_DFT4
};
//static const fn_fft_stage ifft_stg_tab_s3[] =
//{
//    (fn_fft_stage)NULL,//(fn_fft_stage)ifft_16x16_stage_first_scl3_DFT3,
//    (fn_fft_stage)NULL, // not implemented yet fft_16x16_stage_inner_scl3_DFT4,
//    (fn_fft_stage)NULL,//(fn_fft_stage)fft_16x16_stage_last_scl3_DFT4
//};

const fft_cplx_x16_descr_t __cfft_x16_descr48 = 
{
    N, 
    tw_step_tab,
    tw_tab,
    NULL,
    fft_stg_tab_s2,
    NULL //, 
   // fft_stg_tab_s3
};     
const fft_cplx_x16_descr_t __cifft_x16_descr48 =
{
    N, 
    tw_step_tab,
    tw_tab,
    NULL,
    ifft_stg_tab_s2,
    NULL //, 
   // ifft_stg_tab_s3
};     

const fft_handle_t cnfft16_48  = (const fft_handle_t)&__cfft_x16_descr48;
const fft_handle_t cinfft16_48 = (const fft_handle_t)&__cifft_x16_descr48;
