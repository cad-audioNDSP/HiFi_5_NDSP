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

/* Twiddles tables for fft_cplx32x16, ifft_cplx32x16, fft_cplx16x16, ifft_cplx16x16, N=192 */
#define N 192


/********** Twiddles table N=192 stage 1 radix 4 ******************/
ALIGN(32) static const int16_t _fft192_tw1_[] =
{
    (int16_t)0x7FFF, (int16_t)0x0000, (int16_t)0x7FEE, (int16_t)0xFBD0, (int16_t)0x7FFF, (int16_t)0x0000, (int16_t)0x7FBA, (int16_t)0xF7A1,
    (int16_t)0x7FFF, (int16_t)0x0000, (int16_t)0x7F62, (int16_t)0xF374, (int16_t)0x7FBA, (int16_t)0xF7A1, (int16_t)0x7F62, (int16_t)0xF374,
    (int16_t)0x7EE8, (int16_t)0xEF4B, (int16_t)0x7D8A, (int16_t)0xE707, (int16_t)0x7D8A, (int16_t)0xE707, (int16_t)0x7A7D, (int16_t)0xDAD8,
    (int16_t)0x7EE8, (int16_t)0xEF4B, (int16_t)0x7E4A, (int16_t)0xEB26, (int16_t)0x7BA3, (int16_t)0xDEDF, (int16_t)0x7935, (int16_t)0xD6DB,
    (int16_t)0x7642, (int16_t)0xCF04, (int16_t)0x70E3, (int16_t)0xC3A9, (int16_t)0x7D8A, (int16_t)0xE707, (int16_t)0x7CA8, (int16_t)0xE2EF,
    (int16_t)0x7642, (int16_t)0xCF04, (int16_t)0x72CD, (int16_t)0xC763, (int16_t)0x6A6E, (int16_t)0xB8E3, (int16_t)0x62F2, (int16_t)0xAECC,
    (int16_t)0x7BA3, (int16_t)0xDEDF, (int16_t)0x7A7D, (int16_t)0xDAD8, (int16_t)0x6EDA, (int16_t)0xC000, (int16_t)0x6A6E, (int16_t)0xB8E3,
    (int16_t)0x5A82, (int16_t)0xA57E, (int16_t)0x5134, (int16_t)0x9D0E, (int16_t)0x7935, (int16_t)0xD6DB, (int16_t)0x77CC, (int16_t)0xD2E9,
    (int16_t)0x658D, (int16_t)0xB214, (int16_t)0x603C, (int16_t)0xAB9B, (int16_t)0x471D, (int16_t)0x9592, (int16_t)0x3C57, (int16_t)0x8F1D,
    (int16_t)0x7642, (int16_t)0xCF04, (int16_t)0x7497, (int16_t)0xCB2C, (int16_t)0x5A82, (int16_t)0xA57E, (int16_t)0x5465, (int16_t)0x9FC4,
    (int16_t)0x30FC, (int16_t)0x89BE, (int16_t)0x2528, (int16_t)0x8583, (int16_t)0x72CD, (int16_t)0xC763, (int16_t)0x70E3, (int16_t)0xC3A9,
    (int16_t)0x4DEC, (int16_t)0x9A73, (int16_t)0x471D, (int16_t)0x9592, (int16_t)0x18F9, (int16_t)0x8276, (int16_t)0x0C8C, (int16_t)0x809E,
    (int16_t)0x6EDA, (int16_t)0xC000, (int16_t)0x6CB3, (int16_t)0xBC68, (int16_t)0x4000, (int16_t)0x9126, (int16_t)0x389D, (int16_t)0x8D33,
    (int16_t)0x0000, (int16_t)0x8000, (int16_t)0xF374, (int16_t)0x809E, (int16_t)0x6A6E, (int16_t)0xB8E3, (int16_t)0x680B, (int16_t)0xB571,
    (int16_t)0x30FC, (int16_t)0x89BE, (int16_t)0x2925, (int16_t)0x86CB, (int16_t)0xE707, (int16_t)0x8276, (int16_t)0xDAD8, (int16_t)0x8583,
    (int16_t)0x658D, (int16_t)0xB214, (int16_t)0x62F2, (int16_t)0xAECC, (int16_t)0x2121, (int16_t)0x845D, (int16_t)0x18F9, (int16_t)0x8276,
    (int16_t)0xCF04, (int16_t)0x89BE, (int16_t)0xC3A9, (int16_t)0x8F1D, (int16_t)0x603C, (int16_t)0xAB9B, (int16_t)0x5D6C, (int16_t)0xA880,
    (int16_t)0x10B5, (int16_t)0x8118, (int16_t)0x085F, (int16_t)0x8046, (int16_t)0xB8E3, (int16_t)0x9592, (int16_t)0xAECC, (int16_t)0x9D0E,
    (int16_t)0x5A82, (int16_t)0xA57E, (int16_t)0x5780, (int16_t)0xA294, (int16_t)0x0000, (int16_t)0x8000, (int16_t)0xF7A1, (int16_t)0x8046,
    (int16_t)0xA57E, (int16_t)0xA57E, (int16_t)0x9D0E, (int16_t)0xAECC, (int16_t)0x5465, (int16_t)0x9FC4, (int16_t)0x5134, (int16_t)0x9D0E,
    (int16_t)0xEF4B, (int16_t)0x8118, (int16_t)0xE707, (int16_t)0x8276, (int16_t)0x9592, (int16_t)0xB8E3, (int16_t)0x8F1D, (int16_t)0xC3A9,
    (int16_t)0x4DEC, (int16_t)0x9A73, (int16_t)0x4A8F, (int16_t)0x97F5, (int16_t)0xDEDF, (int16_t)0x845D, (int16_t)0xD6DB, (int16_t)0x86CB,
    (int16_t)0x89BE, (int16_t)0xCF04, (int16_t)0x8583, (int16_t)0xDAD8, (int16_t)0x471D, (int16_t)0x9592, (int16_t)0x4398, (int16_t)0x934D,
    (int16_t)0xCF04, (int16_t)0x89BE, (int16_t)0xC763, (int16_t)0x8D33, (int16_t)0x8276, (int16_t)0xE707, (int16_t)0x809E, (int16_t)0xF374,
    (int16_t)0x4000, (int16_t)0x9126, (int16_t)0x3C57, (int16_t)0x8F1D, (int16_t)0xC000, (int16_t)0x9126, (int16_t)0xB8E3, (int16_t)0x9592,
    (int16_t)0x8000, (int16_t)0x0000, (int16_t)0x809E, (int16_t)0x0C8C, (int16_t)0x389D, (int16_t)0x8D33, (int16_t)0x34D4, (int16_t)0x8B69,
    (int16_t)0xB214, (int16_t)0x9A73, (int16_t)0xAB9B, (int16_t)0x9FC4, (int16_t)0x8276, (int16_t)0x18F9, (int16_t)0x8583, (int16_t)0x2528,
    (int16_t)0x30FC, (int16_t)0x89BE, (int16_t)0x2D17, (int16_t)0x8834, (int16_t)0xA57E, (int16_t)0xA57E, (int16_t)0x9FC4, (int16_t)0xAB9B,
    (int16_t)0x89BE, (int16_t)0x30FC, (int16_t)0x8F1D, (int16_t)0x3C57, (int16_t)0x2925, (int16_t)0x86CB, (int16_t)0x2528, (int16_t)0x8583,
    (int16_t)0x9A73, (int16_t)0xB214, (int16_t)0x9592, (int16_t)0xB8E3, (int16_t)0x9592, (int16_t)0x471D, (int16_t)0x9D0E, (int16_t)0x5134,
    (int16_t)0x2121, (int16_t)0x845D, (int16_t)0x1D11, (int16_t)0x8358, (int16_t)0x9126, (int16_t)0xC000, (int16_t)0x8D33, (int16_t)0xC763,
    (int16_t)0xA57E, (int16_t)0x5A82, (int16_t)0xAECC, (int16_t)0x62F2, (int16_t)0x18F9, (int16_t)0x8276, (int16_t)0x14DA, (int16_t)0x81B6,
    (int16_t)0x89BE, (int16_t)0xCF04, (int16_t)0x86CB, (int16_t)0xD6DB, (int16_t)0xB8E3, (int16_t)0x6A6E, (int16_t)0xC3A9, (int16_t)0x70E3,
    (int16_t)0x10B5, (int16_t)0x8118, (int16_t)0x0C8C, (int16_t)0x809E, (int16_t)0x845D, (int16_t)0xDEDF, (int16_t)0x8276, (int16_t)0xE707,
    (int16_t)0xCF04, (int16_t)0x7642, (int16_t)0xDAD8, (int16_t)0x7A7D, (int16_t)0x085F, (int16_t)0x8046, (int16_t)0x0430, (int16_t)0x8012,
    (int16_t)0x8118, (int16_t)0xEF4B, (int16_t)0x8046, (int16_t)0xF7A1, (int16_t)0xE707, (int16_t)0x7D8A, (int16_t)0xF374, (int16_t)0x7F62,
};

/********** Twiddles table N=192 stage 2 radix 6 ******************/
ALIGN(32) static const int16_t _fft192_tw2_[] =
{
    (int16_t)0x7FFF, (int16_t)0x0000, (int16_t)0x7FFF, (int16_t)0x0000, (int16_t)0x7FFF, (int16_t)0x0000, (int16_t)0x7FFF, (int16_t)0x0000,
    (int16_t)0x7FFF, (int16_t)0x0000, (int16_t)0x7EE8, (int16_t)0xEF4B, (int16_t)0x7BA3, (int16_t)0xDEDF, (int16_t)0x7642, (int16_t)0xCF04,
    (int16_t)0x6EDA, (int16_t)0xC000, (int16_t)0x658D, (int16_t)0xB214, (int16_t)0x7BA3, (int16_t)0xDEDF, (int16_t)0x6EDA, (int16_t)0xC000,
    (int16_t)0x5A82, (int16_t)0xA57E, (int16_t)0x4000, (int16_t)0x9126, (int16_t)0x2121, (int16_t)0x845D, (int16_t)0x7642, (int16_t)0xCF04,
    (int16_t)0x5A82, (int16_t)0xA57E, (int16_t)0x30FC, (int16_t)0x89BE, (int16_t)0x0000, (int16_t)0x8000, (int16_t)0xCF04, (int16_t)0x89BE,
    (int16_t)0x6EDA, (int16_t)0xC000, (int16_t)0x4000, (int16_t)0x9126, (int16_t)0x0000, (int16_t)0x8000, (int16_t)0xC000, (int16_t)0x9126,
    (int16_t)0x9126, (int16_t)0xC000, (int16_t)0x658D, (int16_t)0xB214, (int16_t)0x2121, (int16_t)0x845D, (int16_t)0xCF04, (int16_t)0x89BE,
    (int16_t)0x9126, (int16_t)0xC000, (int16_t)0x8118, (int16_t)0x10B5, (int16_t)0x5A82, (int16_t)0xA57E, (int16_t)0x0000, (int16_t)0x8000,
    (int16_t)0xA57E, (int16_t)0xA57E, (int16_t)0x8000, (int16_t)0x0000, (int16_t)0xA57E, (int16_t)0x5A82, (int16_t)0x4DEC, (int16_t)0x9A73,
    (int16_t)0xDEDF, (int16_t)0x845D, (int16_t)0x89BE, (int16_t)0xCF04, (int16_t)0x9126, (int16_t)0x4000, (int16_t)0xEF4B, (int16_t)0x7EE8,
};

static const int tw_step_tab[] =
{
    1, 1, 0
}; 

static const cint16ptr_t_fft tw_tab[] = 
{
    _fft192_tw1_, _fft192_tw2_,  NULL
};


static const fn_fft_stage fft_stg_tab_s2[] = 
{
    (fn_fft_stage)fft_16x16_stage_first_scl2_DFT4,
    (fn_fft_stage)fft_16x16_stage_inner_scl2_DFT6,
    (fn_fft_stage)fft_16x16_stage_last_scl2_DFT8
};
static const fn_fft_stage fft_stg_tab_s3[] = 
{
    (fn_fft_stage)fft_16x16_stage_first_scl3_DFT4,
    (fn_fft_stage)fft_16x16_stage_inner_scl3_DFT6,
    (fn_fft_stage)fft_16x16_stage_last_scl3_DFT8
};
static const fn_fft_stage ifft_stg_tab_s2[] =
{
    (fn_fft_stage)ifft_16x16_stage_first_scl2_DFT4,
    (fn_fft_stage)fft_16x16_stage_inner_scl2_DFT6,
    (fn_fft_stage)ifft_16x16_stage_last_scl2_DFT8
};
static const fn_fft_stage ifft_stg_tab_s3[] =
{
    (fn_fft_stage)ifft_16x16_stage_first_scl3_DFT4,
    (fn_fft_stage)fft_16x16_stage_inner_scl3_DFT6,
    (fn_fft_stage)ifft_16x16_stage_last_scl3_DFT8
};
const fft_cplx_x16_descr_t __cfft_x16_descr192 = 
{
    N, 
    tw_step_tab,
    tw_tab,
    NULL,
    fft_stg_tab_s2,
    NULL, 
    fft_stg_tab_s3
};     
const fft_cplx_x16_descr_t __cifft_x16_descr192 =
{
    N, 
    tw_step_tab,
    tw_tab,
    NULL,
    ifft_stg_tab_s2,
    NULL, 
    ifft_stg_tab_s3
};     
const fft_handle_t cnfft16_192 = (const fft_handle_t)&__cfft_x16_descr192;
const fft_handle_t cinfft16_192 = (const fft_handle_t)&__cifft_x16_descr192;
