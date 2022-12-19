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

/* Twiddles tables for fft_real32x16, ifft_real32x16, fft_real16x16, ifft_real16x16, N=480 */
ALIGN(32) static const int16_t  __fft_real16_tw288[144] =
{
    (int16_t)0x0000, (int16_t)0x7FFF, (int16_t)0x02CB, (int16_t)0x7FF8, (int16_t)0x0595, (int16_t)0x7FE1, (int16_t)0x085F, (int16_t)0x7FBA,
    (int16_t)0x0B28, (int16_t)0x7F83, (int16_t)0x0DEF, (int16_t)0x7F3D, (int16_t)0x10B5, (int16_t)0x7EE8, (int16_t)0x1379, (int16_t)0x7E83,
    (int16_t)0x163A, (int16_t)0x7E0E, (int16_t)0x18F9, (int16_t)0x7D8A, (int16_t)0x1BB4, (int16_t)0x7CF7, (int16_t)0x1E6C, (int16_t)0x7C55,
    (int16_t)0x2121, (int16_t)0x7BA3, (int16_t)0x23D1, (int16_t)0x7AE3, (int16_t)0x267E, (int16_t)0x7A13, (int16_t)0x2925, (int16_t)0x7935,
    (int16_t)0x2BC7, (int16_t)0x7848, (int16_t)0x2E64, (int16_t)0x774C, (int16_t)0x30FC, (int16_t)0x7642, (int16_t)0x338D, (int16_t)0x7529,
    (int16_t)0x3618, (int16_t)0x7402, (int16_t)0x389D, (int16_t)0x72CD, (int16_t)0x3B1B, (int16_t)0x718A, (int16_t)0x3D91, (int16_t)0x7039,
    (int16_t)0x4000, (int16_t)0x6EDA, (int16_t)0x4267, (int16_t)0x6D6E, (int16_t)0x44C6, (int16_t)0x6BF4, (int16_t)0x471D, (int16_t)0x6A6E,
    (int16_t)0x496B, (int16_t)0x68DA, (int16_t)0x4BB0, (int16_t)0x673A, (int16_t)0x4DEC, (int16_t)0x658D, (int16_t)0x501E, (int16_t)0x63D3,
    (int16_t)0x5247, (int16_t)0x620E, (int16_t)0x5465, (int16_t)0x603C, (int16_t)0x567A, (int16_t)0x5E5F, (int16_t)0x5883, (int16_t)0x5C76,
    (int16_t)0x5A82, (int16_t)0x5A82, (int16_t)0x5C76, (int16_t)0x5883, (int16_t)0x5E5F, (int16_t)0x567A, (int16_t)0x603C, (int16_t)0x5465,
    (int16_t)0x620E, (int16_t)0x5247, (int16_t)0x63D3, (int16_t)0x501E, (int16_t)0x658D, (int16_t)0x4DEC, (int16_t)0x673A, (int16_t)0x4BB0,
    (int16_t)0x68DA, (int16_t)0x496B, (int16_t)0x6A6E, (int16_t)0x471D, (int16_t)0x6BF4, (int16_t)0x44C6, (int16_t)0x6D6E, (int16_t)0x4267,
    (int16_t)0x6EDA, (int16_t)0x4000, (int16_t)0x7039, (int16_t)0x3D91, (int16_t)0x718A, (int16_t)0x3B1B, (int16_t)0x72CD, (int16_t)0x389D,
    (int16_t)0x7402, (int16_t)0x3618, (int16_t)0x7529, (int16_t)0x338D, (int16_t)0x7642, (int16_t)0x30FC, (int16_t)0x774C, (int16_t)0x2E64,
    (int16_t)0x7848, (int16_t)0x2BC7, (int16_t)0x7935, (int16_t)0x2925, (int16_t)0x7A13, (int16_t)0x267E, (int16_t)0x7AE3, (int16_t)0x23D1,
    (int16_t)0x7BA3, (int16_t)0x2121, (int16_t)0x7C55, (int16_t)0x1E6C, (int16_t)0x7CF7, (int16_t)0x1BB4, (int16_t)0x7D8A, (int16_t)0x18F9,
    (int16_t)0x7E0E, (int16_t)0x163A, (int16_t)0x7E83, (int16_t)0x1379, (int16_t)0x7EE8, (int16_t)0x10B5, (int16_t)0x7F3D, (int16_t)0x0DEF,
    (int16_t)0x7F83, (int16_t)0x0B28, (int16_t)0x7FBA, (int16_t)0x085F, (int16_t)0x7FE1, (int16_t)0x0595, (int16_t)0x7FF8, (int16_t)0x02CB,
};

static const fft_real_x16_descr_t __rfft_descr =
{
    &__cfft_x16_descr144,
    __fft_real16_tw288
};
static const fft_real_x16_descr_t __rifft_descr =
{
    &__cifft_x16_descr144,
    __fft_real16_tw288
};
const fft_handle_t rnfft16_288 =  (const fft_handle_t)&__rfft_descr;
const fft_handle_t rinfft16_288 = (const fft_handle_t)&__rifft_descr;
