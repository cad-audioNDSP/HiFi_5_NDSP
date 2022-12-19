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

/* Twiddles tables for fft_real32x16, ifft_real32x16, fft_real16x16, ifft_real16x16, N=64 */
ALIGN(32) static const int16_t __fft_real16_tw[] =
{
    (int16_t)0x0000,(int16_t)0x7fff,
    (int16_t)0x0c8c,(int16_t)0x7f62,
    (int16_t)0x18f9,(int16_t)0x7d8a,
    (int16_t)0x2528,(int16_t)0x7a7d,
    (int16_t)0x30fc,(int16_t)0x7642,
    (int16_t)0x3c57,(int16_t)0x70e3,
    (int16_t)0x471d,(int16_t)0x6a6e,
    (int16_t)0x5134,(int16_t)0x62f2,
    (int16_t)0x5a82,(int16_t)0x5a82,
    (int16_t)0x62f2,(int16_t)0x5134,
    (int16_t)0x6a6e,(int16_t)0x471d,
    (int16_t)0x70e3,(int16_t)0x3c57,
    (int16_t)0x7642,(int16_t)0x30fc,
    (int16_t)0x7a7d,(int16_t)0x2528,
    (int16_t)0x7d8a,(int16_t)0x18f9,
    (int16_t)0x7f62,(int16_t)0x0c8c
};

static const fft_real_x16_descr_t __rfft_descr =
{
    &__cfft_x16_descr32,
    __fft_real16_tw
};
static const fft_real_x16_descr_t __rifft_descr =
{
    &__cifft_x16_descr32,
    __fft_real16_tw
};
const fft_handle_t rfft16_64 =  (const fft_handle_t)&__rfft_descr;
const fft_handle_t rifft16_64 = (const fft_handle_t)&__rifft_descr;
