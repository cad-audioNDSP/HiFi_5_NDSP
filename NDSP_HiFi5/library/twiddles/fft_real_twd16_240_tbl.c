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
/* Twiddles tables for fft_cplx32x16, ifft_cplx32x16, fft_cplx16x16, ifft_cplx16x16, N=120 */
#define N 120


/********** Twiddles table N=120 stage 1 radix 4 ******************/
ALIGN(32) static const int16_t _fft120_tw1_[] =
{
    (int16_t)0x7FFF, (int16_t)0x0000, (int16_t)0x7FD3, (int16_t)0xF94D, (int16_t)0x7FFF, (int16_t)0x0000,
    (int16_t)0x7F4C, (int16_t)0xF29F, (int16_t)0x7FFF, (int16_t)0x0000, (int16_t)0x7E6D, (int16_t)0xEBFA,
    (int16_t)0x7F4C, (int16_t)0xF29F, (int16_t)0x7E6D, (int16_t)0xEBFA, (int16_t)0x7D34, (int16_t)0xE563,
    (int16_t)0x79BC, (int16_t)0xD872, (int16_t)0x79BC, (int16_t)0xD872, (int16_t)0x720D, (int16_t)0xC5E4,
    (int16_t)0x7D34, (int16_t)0xE563, (int16_t)0x7BA3, (int16_t)0xDEDF, (int16_t)0x74EF, (int16_t)0xCBF0,
    (int16_t)0x6EDA, (int16_t)0xC000, (int16_t)0x678E, (int16_t)0xB4C3, (int16_t)0x5A82, (int16_t)0xA57E,
    (int16_t)0x79BC, (int16_t)0xD872, (int16_t)0x7780, (int16_t)0xD221, (int16_t)0x678E, (int16_t)0xB4C3,
    (int16_t)0x5F1F, (int16_t)0xAA5A, (int16_t)0x4B3D, (int16_t)0x9872, (int16_t)0x3A1C, (int16_t)0x8DF3,
    (int16_t)0x74EF, (int16_t)0xCBF0, (int16_t)0x720D, (int16_t)0xC5E4, (int16_t)0x55A6, (int16_t)0xA0E1,
    (int16_t)0x4B3D, (int16_t)0x9872, (int16_t)0x278E, (int16_t)0x8644, (int16_t)0x1406, (int16_t)0x8193,
    (int16_t)0x6EDA, (int16_t)0xC000, (int16_t)0x6B5A, (int16_t)0xBA49, (int16_t)0x4000, (int16_t)0x9126,
    (int16_t)0x3410, (int16_t)0x8B11, (int16_t)0x0000, (int16_t)0x8000, (int16_t)0xEBFA, (int16_t)0x8193,
    (int16_t)0x678E, (int16_t)0xB4C3, (int16_t)0x637A, (int16_t)0xAF72, (int16_t)0x278E, (int16_t)0x8644,
    (int16_t)0x1A9D, (int16_t)0x82CC, (int16_t)0xD872, (int16_t)0x8644, (int16_t)0xC5E4, (int16_t)0x8DF3,
    (int16_t)0x5F1F, (int16_t)0xAA5A, (int16_t)0x5A82, (int16_t)0xA57E, (int16_t)0x0D61, (int16_t)0x80B4,
    (int16_t)0x0000, (int16_t)0x8000, (int16_t)0xB4C3, (int16_t)0x9872, (int16_t)0xA57E, (int16_t)0xA57E,
    (int16_t)0x55A6, (int16_t)0xA0E1, (int16_t)0x508E, (int16_t)0x9C86, (int16_t)0xF29F, (int16_t)0x80B4,
    (int16_t)0xE563, (int16_t)0x82CC, (int16_t)0x9872, (int16_t)0xB4C3, (int16_t)0x8DF3, (int16_t)0xC5E4,
    (int16_t)0x4B3D, (int16_t)0x9872, (int16_t)0x45B7, (int16_t)0x94A6, (int16_t)0xD872, (int16_t)0x8644,
    (int16_t)0xCBF0, (int16_t)0x8B11, (int16_t)0x8644, (int16_t)0xD872, (int16_t)0x8193, (int16_t)0xEBFA,
    (int16_t)0x4000, (int16_t)0x9126, (int16_t)0x3A1C, (int16_t)0x8DF3, (int16_t)0xC000, (int16_t)0x9126,
    (int16_t)0xB4C3, (int16_t)0x9872, (int16_t)0x8000, (int16_t)0x0000, (int16_t)0x8193, (int16_t)0x1406,
    (int16_t)0x3410, (int16_t)0x8B11, (int16_t)0x2DDF, (int16_t)0x8880, (int16_t)0xAA5A, (int16_t)0xA0E1,
    (int16_t)0xA0E1, (int16_t)0xAA5A, (int16_t)0x8644, (int16_t)0x278E, (int16_t)0x8DF3, (int16_t)0x3A1C,
    (int16_t)0x278E, (int16_t)0x8644, (int16_t)0x2121, (int16_t)0x845D, (int16_t)0x9872, (int16_t)0xB4C3,
    (int16_t)0x9126, (int16_t)0xC000, (int16_t)0x9872, (int16_t)0x4B3D, (int16_t)0xA57E, (int16_t)0x5A82,
    (int16_t)0x1A9D, (int16_t)0x82CC, (int16_t)0x1406, (int16_t)0x8193, (int16_t)0x8B11, (int16_t)0xCBF0,
    (int16_t)0x8644, (int16_t)0xD872, (int16_t)0xB4C3, (int16_t)0x678E, (int16_t)0xC5E4, (int16_t)0x720D,
    (int16_t)0x0D61, (int16_t)0x80B4, (int16_t)0x06B3, (int16_t)0x802D, (int16_t)0x82CC, (int16_t)0xE563,
    (int16_t)0x80B4, (int16_t)0xF29F, (int16_t)0xD872, (int16_t)0x79BC, (int16_t)0xEBFA, (int16_t)0x7E6D,
};

/********** Twiddles table N=120 stage 2 radix 5 ******************/
ALIGN(32) static const int16_t _fft120_tw2_[] =
{
    (int16_t)0x7FFF, (int16_t)0x0000, (int16_t)0x7FFF, (int16_t)0x0000, (int16_t)0x7FFF, (int16_t)0x0000, (int16_t)0x7FFF, (int16_t)0x0000,
    (int16_t)0x7D34, (int16_t)0xE563, (int16_t)0x74EF, (int16_t)0xCBF0, (int16_t)0x678E, (int16_t)0xB4C3, (int16_t)0x55A6, (int16_t)0xA0E1,
    (int16_t)0x74EF, (int16_t)0xCBF0, (int16_t)0x55A6, (int16_t)0xA0E1, (int16_t)0x278E, (int16_t)0x8644, (int16_t)0xF29F, (int16_t)0x80B4,
    (int16_t)0x678E, (int16_t)0xB4C3, (int16_t)0x278E, (int16_t)0x8644, (int16_t)0xD872, (int16_t)0x8644, (int16_t)0x9872, (int16_t)0xB4C3,
    (int16_t)0x55A6, (int16_t)0xA0E1, (int16_t)0xF29F, (int16_t)0x80B4, (int16_t)0x9872, (int16_t)0xB4C3, (int16_t)0x82CC, (int16_t)0x1A9D,
    (int16_t)0x4000, (int16_t)0x9126, (int16_t)0xC000, (int16_t)0x9126, (int16_t)0x8000, (int16_t)0x0000, (int16_t)0xC000, (int16_t)0x6EDA,
};
static const int tw_step_tab[] =
{
    1, 1, 0
}; 
static const cint16ptr_t_fft tw_tab[] = 
{
    _fft120_tw1_, _fft120_tw2_, NULL
};

static const fn_fft_stage fft_stg_tab_s2[] = 
{
    (fn_fft_stage)fft_16x16_stage_first_scl2_DFT4,
    (fn_fft_stage)fft_16x16_stage_inner_scl2_DFT5,
    (fn_fft_stage)fft_16x16_stage_last_scl2_DFT6
};
static const fn_fft_stage fft_stg_tab_s3[] = 
{
    (fn_fft_stage)fft_16x16_stage_first_scl3_DFT4,
    (fn_fft_stage)fft_16x16_stage_inner_scl3_DFT5,
    (fn_fft_stage)fft_16x16_stage_last_scl3_DFT6
};
static const fn_fft_stage ifft_stg_tab_s2[] =
{
    (fn_fft_stage)ifft_16x16_stage_first_scl2_DFT4,
    (fn_fft_stage)fft_16x16_stage_inner_scl2_DFT5,
    (fn_fft_stage)ifft_16x16_stage_last_scl2_DFT6
};
static const fn_fft_stage ifft_stg_tab_s3[] =
{
    (fn_fft_stage)ifft_16x16_stage_first_scl3_DFT4,
    (fn_fft_stage)fft_16x16_stage_inner_scl3_DFT5,
    (fn_fft_stage)ifft_16x16_stage_last_scl3_DFT6
};
static const fft_cplx_x16_descr_t __cfft_descr = 
{
    N, 
    tw_step_tab,
    tw_tab,
    NULL,
    fft_stg_tab_s2,
    NULL, 
    fft_stg_tab_s3
};     
static const fft_cplx_x16_descr_t __cifft_descr =
{
    N, 
    tw_step_tab,
    tw_tab,
    NULL,
    ifft_stg_tab_s2,
    NULL, 
    ifft_stg_tab_s3
};     

/* Twiddles tables for fft_real32x16, ifft_real32x16, fft_real16x16, ifft_real16x16, N=240 */
ALIGN(32) static const int16_t __fft_real16_tw[] =
{
(int16_t)0x0000,(int16_t)0x7fff,
(int16_t)0x035a,(int16_t)0x7ff5,
(int16_t)0x06b3,(int16_t)0x7fd3,
(int16_t)0x0a0b,(int16_t)0x7f9b,
(int16_t)0x0d61,(int16_t)0x7f4c,
(int16_t)0x10b5,(int16_t)0x7ee8,
(int16_t)0x1406,(int16_t)0x7e6d,
(int16_t)0x1753,(int16_t)0x7ddb,
(int16_t)0x1a9d,(int16_t)0x7d34,
(int16_t)0x1de2,(int16_t)0x7c77,
(int16_t)0x2121,(int16_t)0x7ba3,
(int16_t)0x245b,(int16_t)0x7abb,
(int16_t)0x278e,(int16_t)0x79bc,
(int16_t)0x2aba,(int16_t)0x78a8,
(int16_t)0x2ddf,(int16_t)0x7780,
(int16_t)0x30fc,(int16_t)0x7642,
(int16_t)0x3410,(int16_t)0x74ef,
(int16_t)0x371b,(int16_t)0x7388,
(int16_t)0x3a1c,(int16_t)0x720d,
(int16_t)0x3d14,(int16_t)0x707d,
(int16_t)0x4000,(int16_t)0x6eda,
(int16_t)0x42e1,(int16_t)0x6d23,
(int16_t)0x45b7,(int16_t)0x6b5a,
(int16_t)0x4880,(int16_t)0x697d,
(int16_t)0x4b3d,(int16_t)0x678e,
(int16_t)0x4dec,(int16_t)0x658d,
(int16_t)0x508e,(int16_t)0x637a,
(int16_t)0x5321,(int16_t)0x6155,
(int16_t)0x55a6,(int16_t)0x5f1f,
(int16_t)0x581c,(int16_t)0x5cd9,
(int16_t)0x5a82,(int16_t)0x5a82,
(int16_t)0x5cd9,(int16_t)0x581c,
(int16_t)0x5f1f,(int16_t)0x55a6,
(int16_t)0x6155,(int16_t)0x5321,
(int16_t)0x637a,(int16_t)0x508e,
(int16_t)0x658d,(int16_t)0x4dec,
(int16_t)0x678e,(int16_t)0x4b3d,
(int16_t)0x697d,(int16_t)0x4880,
(int16_t)0x6b5a,(int16_t)0x45b7,
(int16_t)0x6d23,(int16_t)0x42e1,
(int16_t)0x6eda,(int16_t)0x4000,
(int16_t)0x707d,(int16_t)0x3d14,
(int16_t)0x720d,(int16_t)0x3a1c,
(int16_t)0x7388,(int16_t)0x371b,
(int16_t)0x74ef,(int16_t)0x3410,
(int16_t)0x7642,(int16_t)0x30fc,
(int16_t)0x7780,(int16_t)0x2ddf,
(int16_t)0x78a8,(int16_t)0x2aba,
(int16_t)0x79bc,(int16_t)0x278e,
(int16_t)0x7abb,(int16_t)0x245b,
(int16_t)0x7ba3,(int16_t)0x2121,
(int16_t)0x7c77,(int16_t)0x1de2,
(int16_t)0x7d34,(int16_t)0x1a9d,
(int16_t)0x7ddb,(int16_t)0x1753,
(int16_t)0x7e6d,(int16_t)0x1406,
(int16_t)0x7ee8,(int16_t)0x10b5,
(int16_t)0x7f4c,(int16_t)0x0d61,
(int16_t)0x7f9b,(int16_t)0x0a0b,
(int16_t)0x7fd3,(int16_t)0x06b3,
(int16_t)0x7ff5,(int16_t)0x035a
};

static const fft_real_x16_descr_t __rfft_descr =
{
    &__cfft_descr,
    __fft_real16_tw
};
static const fft_real_x16_descr_t __rifft_descr =
{
    &__cifft_descr,
    __fft_real16_tw
};
const fft_handle_t rnfft16_240 =  (const fft_handle_t)&__rfft_descr;
const fft_handle_t rinfft16_240 = (const fft_handle_t)&__rifft_descr;
