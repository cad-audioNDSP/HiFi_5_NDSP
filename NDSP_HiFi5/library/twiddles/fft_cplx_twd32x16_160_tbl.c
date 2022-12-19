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
#include "fft_32x16_stages.h"
/* Twiddles tables for fft_cplx32x16, ifft_cplx32x16, fft_cplx16x16, ifft_cplx16x16, N=160 */
#define N 160

/********** Twiddles table N=160 stage 1 radix 4 ******************/
ALIGN(32) static const int16_t _fft160_tw1_[] =
{
    (int16_t)0x7FFF, (int16_t)0x0000, (int16_t)0x7FE7, (int16_t)0xFAFA, (int16_t)0x7FFF, (int16_t)0x0000, (int16_t)0x7F9B, (int16_t)0xF5F5,
    (int16_t)0x7FFF, (int16_t)0x0000, (int16_t)0x7F1D, (int16_t)0xF0F5, (int16_t)0x7F9B, (int16_t)0xF5F5, (int16_t)0x7F1D, (int16_t)0xF0F5,
    (int16_t)0x7E6D, (int16_t)0xEBFA, (int16_t)0x7C77, (int16_t)0xE21E, (int16_t)0x7C77, (int16_t)0xE21E, (int16_t)0x7817, (int16_t)0xD3B2,
    (int16_t)0x7E6D, (int16_t)0xEBFA, (int16_t)0x7D8A, (int16_t)0xE707, (int16_t)0x79BC, (int16_t)0xD872, (int16_t)0x7642, (int16_t)0xCF04,
    (int16_t)0x720D, (int16_t)0xC5E4, (int16_t)0x6A6E, (int16_t)0xB8E3, (int16_t)0x7C77, (int16_t)0xE21E, (int16_t)0x7B32, (int16_t)0xDD41,
    (int16_t)0x720D, (int16_t)0xC5E4, (int16_t)0x6D23, (int16_t)0xBD1F, (int16_t)0x6155, (int16_t)0xACDF, (int16_t)0x56E3, (int16_t)0xA202,
    (int16_t)0x79BC, (int16_t)0xD872, (int16_t)0x7817, (int16_t)0xD3B2, (int16_t)0x678E, (int16_t)0xB4C3, (int16_t)0x6155, (int16_t)0xACDF,
    (int16_t)0x4B3D, (int16_t)0x9872, (int16_t)0x3E8B, (int16_t)0x9052, (int16_t)0x7642, (int16_t)0xCF04, (int16_t)0x743E, (int16_t)0xCA69,
    (int16_t)0x5A82, (int16_t)0xA57E, (int16_t)0x5321, (int16_t)0x9EAB, (int16_t)0x30FC, (int16_t)0x89BE, (int16_t)0x22BF, (int16_t)0x84CE,
    (int16_t)0x720D, (int16_t)0xC5E4, (int16_t)0x6FAE, (int16_t)0xC175, (int16_t)0x4B3D, (int16_t)0x9872, (int16_t)0x42E1, (int16_t)0x92DD,
    (int16_t)0x1406, (int16_t)0x8193, (int16_t)0x0506, (int16_t)0x8019, (int16_t)0x6D23, (int16_t)0xBD1F, (int16_t)0x6A6E, (int16_t)0xB8E3,
    (int16_t)0x3A1C, (int16_t)0x8DF3, (int16_t)0x30FC, (int16_t)0x89BE, (int16_t)0xF5F5, (int16_t)0x8065, (int16_t)0xE707, (int16_t)0x8276,
    (int16_t)0x678E, (int16_t)0xB4C3, (int16_t)0x6485, (int16_t)0xB0C2, (int16_t)0x278E, (int16_t)0x8644, (int16_t)0x1DE2, (int16_t)0x8389,
    (int16_t)0xD872, (int16_t)0x8644, (int16_t)0xCA69, (int16_t)0x8BC2, (int16_t)0x6155, (int16_t)0xACDF, (int16_t)0x5DFE, (int16_t)0xA91D,
    (int16_t)0x1406, (int16_t)0x8193, (int16_t)0x0A0B, (int16_t)0x8065, (int16_t)0xBD1F, (int16_t)0x92DD, (int16_t)0xB0C2, (int16_t)0x9B7B,
    (int16_t)0x5A82, (int16_t)0xA57E, (int16_t)0x56E3, (int16_t)0xA202, (int16_t)0x0000, (int16_t)0x8000, (int16_t)0xF5F5, (int16_t)0x8065,
    (int16_t)0xA57E, (int16_t)0xA57E, (int16_t)0x9B7B, (int16_t)0xB0C2, (int16_t)0x5321, (int16_t)0x9EAB, (int16_t)0x4F3E, (int16_t)0x9B7B,
    (int16_t)0xEBFA, (int16_t)0x8193, (int16_t)0xE21E, (int16_t)0x8389, (int16_t)0x92DD, (int16_t)0xBD1F, (int16_t)0x8BC2, (int16_t)0xCA69,
    (int16_t)0x4B3D, (int16_t)0x9872, (int16_t)0x471D, (int16_t)0x9592, (int16_t)0xD872, (int16_t)0x8644, (int16_t)0xCF04, (int16_t)0x89BE,
    (int16_t)0x8644, (int16_t)0xD872, (int16_t)0x8276, (int16_t)0xE707, (int16_t)0x42E1, (int16_t)0x92DD, (int16_t)0x3E8B, (int16_t)0x9052,
    (int16_t)0xC5E4, (int16_t)0x8DF3, (int16_t)0xBD1F, (int16_t)0x92DD, (int16_t)0x8065, (int16_t)0xF5F5, (int16_t)0x8019, (int16_t)0x0506,
    (int16_t)0x3A1C, (int16_t)0x8DF3, (int16_t)0x3597, (int16_t)0x8BC2, (int16_t)0xB4C3, (int16_t)0x9872, (int16_t)0xACDF, (int16_t)0x9EAB,
    (int16_t)0x8193, (int16_t)0x1406, (int16_t)0x84CE, (int16_t)0x22BF, (int16_t)0x30FC, (int16_t)0x89BE, (int16_t)0x2C4E, (int16_t)0x87E9,
    (int16_t)0xA57E, (int16_t)0xA57E, (int16_t)0x9EAB, (int16_t)0xACDF, (int16_t)0x89BE, (int16_t)0x30FC, (int16_t)0x9052, (int16_t)0x3E8B,
    (int16_t)0x278E, (int16_t)0x8644, (int16_t)0x22BF, (int16_t)0x84CE, (int16_t)0x9872, (int16_t)0xB4C3, (int16_t)0x92DD, (int16_t)0xBD1F,
    (int16_t)0x9872, (int16_t)0x4B3D, (int16_t)0xA202, (int16_t)0x56E3, (int16_t)0x1DE2, (int16_t)0x8389, (int16_t)0x18F9, (int16_t)0x8276,
    (int16_t)0x8DF3, (int16_t)0xC5E4, (int16_t)0x89BE, (int16_t)0xCF04, (int16_t)0xACDF, (int16_t)0x6155, (int16_t)0xB8E3, (int16_t)0x6A6E,
    (int16_t)0x1406, (int16_t)0x8193, (int16_t)0x0F0B, (int16_t)0x80E3, (int16_t)0x8644, (int16_t)0xD872, (int16_t)0x8389, (int16_t)0xE21E,
    (int16_t)0xC5E4, (int16_t)0x720D, (int16_t)0xD3B2, (int16_t)0x7817, (int16_t)0x0A0B, (int16_t)0x8065, (int16_t)0x0506, (int16_t)0x8019,
    (int16_t)0x8193, (int16_t)0xEBFA, (int16_t)0x8065, (int16_t)0xF5F5, (int16_t)0xE21E, (int16_t)0x7C77, (int16_t)0xF0F5, (int16_t)0x7F1D,
};
/********** Twiddles table N=160 stage 2 radix 5 ******************/
ALIGN(32) static const int16_t _fft160_tw2_[] =
{
    (int16_t)0x7FFF, (int16_t)0x0000, (int16_t)0x7FFF, (int16_t)0x0000, (int16_t)0x7FFF, (int16_t)0x0000, (int16_t)0x7FFF, (int16_t)0x0000,
    (int16_t)0x7E6D, (int16_t)0xEBFA, (int16_t)0x79BC, (int16_t)0xD872, (int16_t)0x720D, (int16_t)0xC5E4, (int16_t)0x678E, (int16_t)0xB4C3,
    (int16_t)0x79BC, (int16_t)0xD872, (int16_t)0x678E, (int16_t)0xB4C3, (int16_t)0x4B3D, (int16_t)0x9872, (int16_t)0x278E, (int16_t)0x8644,
    (int16_t)0x720D, (int16_t)0xC5E4, (int16_t)0x4B3D, (int16_t)0x9872, (int16_t)0x1406, (int16_t)0x8193, (int16_t)0xD872, (int16_t)0x8644,
    (int16_t)0x678E, (int16_t)0xB4C3, (int16_t)0x278E, (int16_t)0x8644, (int16_t)0xD872, (int16_t)0x8644, (int16_t)0x9872, (int16_t)0xB4C3,
    (int16_t)0x5A82, (int16_t)0xA57E, (int16_t)0x0000, (int16_t)0x8000, (int16_t)0xA57E, (int16_t)0xA57E, (int16_t)0x8000, (int16_t)0x0000,
    (int16_t)0x4B3D, (int16_t)0x9872, (int16_t)0xD872, (int16_t)0x8644, (int16_t)0x8644, (int16_t)0xD872, (int16_t)0x9872, (int16_t)0x4B3D,
    (int16_t)0x3A1C, (int16_t)0x8DF3, (int16_t)0xB4C3, (int16_t)0x9872, (int16_t)0x8193, (int16_t)0x1406, (int16_t)0xD872, (int16_t)0x79BC,
};
static const int tw_step_tab[] =
{
    1, 1, 0, 
}; 
static const cint16ptr_t_fft tw_tab[] = 
{
    _fft160_tw1_, _fft160_tw2_, NULL
};

static const fn_fft_stage fft32x16_stg_tab_s2[] = 
{
    (fn_fft_stage)fft_32x16_stage_first_scl2_DFT4,
    (fn_fft_stage)fft_32x16_stage_inner_scl2_DFT5_v4,
    (fn_fft_stage)fft_32x16_stage_last_scl2_DFT8
};
static const fn_fft_stage fft32x16_stg_tab_s3[] = 
{
    (fn_fft_stage)fft_32x16_stage_first_scl3_DFT4,
    (fn_fft_stage)fft_32x16_stage_inner_scl3_DFT5_v4,
    (fn_fft_stage)fft_32x16_stage_last_scl3_DFT8
};
static const fn_fft_stage ifft32x16_stg_tab_s2[] =
{
    (fn_fft_stage)ifft_32x16_stage_first_scl2_DFT4,
    (fn_fft_stage)fft_32x16_stage_inner_scl2_DFT5_v4,
    (fn_fft_stage)ifft_32x16_stage_last_scl2_DFT8
};
static const fn_fft_stage ifft32x16_stg_tab_s3[] =
{
    (fn_fft_stage)ifft_32x16_stage_first_scl3_DFT4,
    (fn_fft_stage)fft_32x16_stage_inner_scl3_DFT5_v4,
    (fn_fft_stage)ifft_32x16_stage_last_scl3_DFT8
};

const fft_cplx_x16_descr_t __cfft_32x16_descr160 = 
{
    N, 
    tw_step_tab,
    tw_tab,
    fft32x16_stg_tab_s2,
    NULL, 
    fft32x16_stg_tab_s3,
    NULL
};     
const fft_cplx_x16_descr_t __cifft_32x16_descr160 = 
{
    N, 
    tw_step_tab,
    tw_tab,
    ifft32x16_stg_tab_s2,
    NULL, 
    ifft32x16_stg_tab_s3,
    NULL
};     
const fft_handle_t cnfft32x16_160 = (const fft_handle_t)&__cfft_32x16_descr160;
const fft_handle_t cinfft32x16_160 = (const fft_handle_t)&__cifft_32x16_descr160;
