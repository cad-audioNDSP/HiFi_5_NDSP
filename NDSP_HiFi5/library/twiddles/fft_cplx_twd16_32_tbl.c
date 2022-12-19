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
#include "fft_32x16_stages.h"

/* Twiddles tables for fft_cplx32x16, ifft_cplx32x16, fft_cplx16x16, ifft_cplx16x16, N=32 */
#define N 32

/****************** stage 1, stage 2 radix 4 ******************/
ALIGN(32) static const int16_t __fft16_tw1[] =
{
    (int16_t)0x7FFF, (int16_t)0x0000, (int16_t)0x7D8A, (int16_t)0xE707, (int16_t)0x7FFF, (int16_t)0x0000, (int16_t)0x7642, (int16_t)0xCF04,
    (int16_t)0x7FFF, (int16_t)0x0000, (int16_t)0x6A6E, (int16_t)0xB8E3, (int16_t)0x7642, (int16_t)0xCF04, (int16_t)0x6A6E, (int16_t)0xB8E3,
    (int16_t)0x5A82, (int16_t)0xA57E, (int16_t)0x30FC, (int16_t)0x89BE, (int16_t)0x30FC, (int16_t)0x89BE, (int16_t)0xE707, (int16_t)0x8276,
    (int16_t)0x5A82, (int16_t)0xA57E, (int16_t)0x471D, (int16_t)0x9592, (int16_t)0x0000, (int16_t)0x8000, (int16_t)0xCF04, (int16_t)0x89BE,
    (int16_t)0xA57E, (int16_t)0xA57E, (int16_t)0x8276, (int16_t)0xE707, (int16_t)0x30FC, (int16_t)0x89BE, (int16_t)0x18F9, (int16_t)0x8276,
    (int16_t)0xA57E, (int16_t)0xA57E, (int16_t)0x89BE, (int16_t)0xCF04, (int16_t)0x89BE, (int16_t)0x30FC, (int16_t)0xB8E3, (int16_t)0x6A6E
};
static const int tw_step_tab[] =
{
    1, 4
}; 
static const cint16ptr_t_fft tw_tab[] = 
{
    __fft16_tw1, NULL
};

static const fn_fft_stage fft32x16_stg_tab_s2[] = 
{
    (fn_fft_stage)fft_32x16_stage_first_scl2_DFT4,
    (fn_fft_stage)fft_32x16_stage_last_scl2_DFT8
};
static const fn_fft_stage fft32x16_stg_tab_s3[] = 
{
    (fn_fft_stage)fft_32x16_stage_first_scl3_DFT4,
    (fn_fft_stage)fft_32x16_stage_last_scl3_DFT8
};
static const fn_fft_stage ifft32x16_stg_tab_s2[] =
{
    (fn_fft_stage)ifft_32x16_stage_first_scl2_DFT4,
    (fn_fft_stage)ifft_32x16_stage_last_scl2_DFT8
};
static const fn_fft_stage ifft32x16_stg_tab_s3[] =
{
    (fn_fft_stage)ifft_32x16_stage_first_scl3_DFT4,
    (fn_fft_stage)ifft_32x16_stage_last_scl3_DFT8
};
static const fn_fft_stage fft_stg_tab_s2[] = 
{
    (fn_fft_stage)fft_16x16_stage_first_scl2_DFT4,
    (fn_fft_stage)fft_16x16_stage_last_scl2_DFT8
};
static const fn_fft_stage fft_stg_tab_s3[] = 
{
    (fn_fft_stage)fft_16x16_stage_first_scl3_DFT4,
    (fn_fft_stage)fft_16x16_stage_last_scl3_DFT8
};
static const fn_fft_stage ifft_stg_tab_s2[] =
{
    (fn_fft_stage)ifft_16x16_stage_first_scl2_DFT4,
    (fn_fft_stage)ifft_16x16_stage_last_scl2_DFT8
};
static const fn_fft_stage ifft_stg_tab_s3[] =
{
    (fn_fft_stage)ifft_16x16_stage_first_scl3_DFT4,
    (fn_fft_stage)ifft_16x16_stage_last_scl3_DFT8
};

const fft_cplx_x16_descr_t __cfft_x16_descr32 = 
{
    N, 
    tw_step_tab,
    tw_tab,
    fft32x16_stg_tab_s2,
    fft_stg_tab_s2,
    fft32x16_stg_tab_s3, 
    fft_stg_tab_s3
};     
const fft_cplx_x16_descr_t __cifft_x16_descr32 =
{
    N, 
    tw_step_tab,
    tw_tab,
    ifft32x16_stg_tab_s2,
    ifft_stg_tab_s2,
    ifft32x16_stg_tab_s3, 
    ifft_stg_tab_s3
};     
const fft_handle_t cfft16_32 = (const fft_handle_t)&__cfft_x16_descr32;
const fft_handle_t cifft16_32 = (const fft_handle_t)&__cifft_x16_descr32;
