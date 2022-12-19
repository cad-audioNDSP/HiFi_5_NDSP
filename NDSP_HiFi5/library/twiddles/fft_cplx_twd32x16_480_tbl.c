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

/* Twiddles tables for fft_cplx32x16, ifft_cplx32x16, fft_cplx16x16, ifft_cplx16x16, N=480 */
#define N 480


/********** Twiddles table N=480 stage 1 radix 4 ******************/
ALIGN(32) static const int16_t _fft480_tw1_[] =
{
    (int16_t)0x7FFF, (int16_t)0x0000, (int16_t)0x7FFD, (int16_t)0xFE53, (int16_t)0x7FFF, (int16_t)0x0000, (int16_t)0x7FF5, (int16_t)0xFCA6,
    (int16_t)0x7FFF, (int16_t)0x0000, (int16_t)0x7FE7, (int16_t)0xFAFA, (int16_t)0x7FF5, (int16_t)0xFCA6, (int16_t)0x7FE7, (int16_t)0xFAFA,
    (int16_t)0x7FD3, (int16_t)0xF94D, (int16_t)0x7F9B, (int16_t)0xF5F5, (int16_t)0x7F9B, (int16_t)0xF5F5, (int16_t)0x7F1D, (int16_t)0xF0F5,
    (int16_t)0x7FD3, (int16_t)0xF94D, (int16_t)0x7FBA, (int16_t)0xF7A1, (int16_t)0x7F4C, (int16_t)0xF29F, (int16_t)0x7EE8, (int16_t)0xEF4B,
    (int16_t)0x7E6D, (int16_t)0xEBFA, (int16_t)0x7D8A, (int16_t)0xE707, (int16_t)0x7F9B, (int16_t)0xF5F5, (int16_t)0x7F77, (int16_t)0xF44A,
    (int16_t)0x7E6D, (int16_t)0xEBFA, (int16_t)0x7DDB, (int16_t)0xE8AD, (int16_t)0x7C77, (int16_t)0xE21E, (int16_t)0x7B32, (int16_t)0xDD41,
    (int16_t)0x7F4C, (int16_t)0xF29F, (int16_t)0x7F1D, (int16_t)0xF0F5, (int16_t)0x7D34, (int16_t)0xE563, (int16_t)0x7C77, (int16_t)0xE21E,
    (int16_t)0x79BC, (int16_t)0xD872, (int16_t)0x7817, (int16_t)0xD3B2, (int16_t)0x7EE8, (int16_t)0xEF4B, (int16_t)0x7EAD, (int16_t)0xEDA2,
    (int16_t)0x7BA3, (int16_t)0xDEDF, (int16_t)0x7ABB, (int16_t)0xDBA5, (int16_t)0x7642, (int16_t)0xCF04, (int16_t)0x743E, (int16_t)0xCA69,
    (int16_t)0x7E6D, (int16_t)0xEBFA, (int16_t)0x7E27, (int16_t)0xEA53, (int16_t)0x79BC, (int16_t)0xD872, (int16_t)0x78A8, (int16_t)0xD546,
    (int16_t)0x720D, (int16_t)0xC5E4, (int16_t)0x6FAE, (int16_t)0xC175, (int16_t)0x7DDB, (int16_t)0xE8AD, (int16_t)0x7D8A, (int16_t)0xE707,
    (int16_t)0x7780, (int16_t)0xD221, (int16_t)0x7642, (int16_t)0xCF04, (int16_t)0x6D23, (int16_t)0xBD1F, (int16_t)0x6A6E, (int16_t)0xB8E3,
    (int16_t)0x7D34, (int16_t)0xE563, (int16_t)0x7CD8, (int16_t)0xE3C0, (int16_t)0x74EF, (int16_t)0xCBF0, (int16_t)0x7388, (int16_t)0xC8E5,
    (int16_t)0x678E, (int16_t)0xB4C3, (int16_t)0x6485, (int16_t)0xB0C2, (int16_t)0x7C77, (int16_t)0xE21E, (int16_t)0x7C10, (int16_t)0xE07E,
    (int16_t)0x720D, (int16_t)0xC5E4, (int16_t)0x707D, (int16_t)0xC2EC, (int16_t)0x6155, (int16_t)0xACDF, (int16_t)0x5DFE, (int16_t)0xA91D,
    (int16_t)0x7BA3, (int16_t)0xDEDF, (int16_t)0x7B32, (int16_t)0xDD41, (int16_t)0x6EDA, (int16_t)0xC000, (int16_t)0x6D23, (int16_t)0xBD1F,
    (int16_t)0x5A82, (int16_t)0xA57E, (int16_t)0x56E3, (int16_t)0xA202, (int16_t)0x7ABB, (int16_t)0xDBA5, (int16_t)0x7A3E, (int16_t)0xDA0B,
    (int16_t)0x6B5A, (int16_t)0xBA49, (int16_t)0x697D, (int16_t)0xB780, (int16_t)0x5321, (int16_t)0x9EAB, (int16_t)0x4F3E, (int16_t)0x9B7B,
    (int16_t)0x79BC, (int16_t)0xD872, (int16_t)0x7935, (int16_t)0xD6DB, (int16_t)0x678E, (int16_t)0xB4C3, (int16_t)0x658D, (int16_t)0xB214,
    (int16_t)0x4B3D, (int16_t)0x9872, (int16_t)0x471D, (int16_t)0x9592, (int16_t)0x78A8, (int16_t)0xD546, (int16_t)0x7817, (int16_t)0xD3B2,
    (int16_t)0x637A, (int16_t)0xAF72, (int16_t)0x6155, (int16_t)0xACDF, (int16_t)0x42E1, (int16_t)0x92DD, (int16_t)0x3E8B, (int16_t)0x9052,
    (int16_t)0x7780, (int16_t)0xD221, (int16_t)0x76E3, (int16_t)0xD092, (int16_t)0x5F1F, (int16_t)0xAA5A, (int16_t)0x5CD9, (int16_t)0xA7E4,
    (int16_t)0x3A1C, (int16_t)0x8DF3, (int16_t)0x3597, (int16_t)0x8BC2, (int16_t)0x7642, (int16_t)0xCF04, (int16_t)0x759B, (int16_t)0xCD79,
    (int16_t)0x5A82, (int16_t)0xA57E, (int16_t)0x581C, (int16_t)0xA327, (int16_t)0x30FC, (int16_t)0x89BE, (int16_t)0x2C4E, (int16_t)0x87E9,
    (int16_t)0x74EF, (int16_t)0xCBF0, (int16_t)0x743E, (int16_t)0xCA69, (int16_t)0x55A6, (int16_t)0xA0E1, (int16_t)0x5321, (int16_t)0x9EAB,
    (int16_t)0x278E, (int16_t)0x8644, (int16_t)0x22BF, (int16_t)0x84CE, (int16_t)0x7388, (int16_t)0xC8E5, (int16_t)0x72CD, (int16_t)0xC763,
    (int16_t)0x508E, (int16_t)0x9C86, (int16_t)0x4DEC, (int16_t)0x9A73, (int16_t)0x1DE2, (int16_t)0x8389, (int16_t)0x18F9, (int16_t)0x8276,
    (int16_t)0x720D, (int16_t)0xC5E4, (int16_t)0x7147, (int16_t)0xC467, (int16_t)0x4B3D, (int16_t)0x9872, (int16_t)0x4880, (int16_t)0x9683,
    (int16_t)0x1406, (int16_t)0x8193, (int16_t)0x0F0B, (int16_t)0x80E3, (int16_t)0x707D, (int16_t)0xC2EC, (int16_t)0x6FAE, (int16_t)0xC175,
    (int16_t)0x45B7, (int16_t)0x94A6, (int16_t)0x42E1, (int16_t)0x92DD, (int16_t)0x0A0B, (int16_t)0x8065, (int16_t)0x0506, (int16_t)0x8019,
    (int16_t)0x6EDA, (int16_t)0xC000, (int16_t)0x6E01, (int16_t)0xBE8E, (int16_t)0x4000, (int16_t)0x9126, (int16_t)0x3D14, (int16_t)0x8F83,
    (int16_t)0x0000, (int16_t)0x8000, (int16_t)0xFAFA, (int16_t)0x8019, (int16_t)0x6D23, (int16_t)0xBD1F, (int16_t)0x6C41, (int16_t)0xBBB3,
    (int16_t)0x3A1C, (int16_t)0x8DF3, (int16_t)0x371B, (int16_t)0x8C78, (int16_t)0xF5F5, (int16_t)0x8065, (int16_t)0xF0F5, (int16_t)0x80E3,
    (int16_t)0x6B5A, (int16_t)0xBA49, (int16_t)0x6A6E, (int16_t)0xB8E3, (int16_t)0x3410, (int16_t)0x8B11, (int16_t)0x30FC, (int16_t)0x89BE,
    (int16_t)0xEBFA, (int16_t)0x8193, (int16_t)0xE707, (int16_t)0x8276, (int16_t)0x697D, (int16_t)0xB780, (int16_t)0x6888, (int16_t)0xB620,
    (int16_t)0x2DDF, (int16_t)0x8880, (int16_t)0x2ABA, (int16_t)0x8758, (int16_t)0xE21E, (int16_t)0x8389, (int16_t)0xDD41, (int16_t)0x84CE,
    (int16_t)0x678E, (int16_t)0xB4C3, (int16_t)0x668F, (int16_t)0xB36A, (int16_t)0x278E, (int16_t)0x8644, (int16_t)0x245B, (int16_t)0x8545,
    (int16_t)0xD872, (int16_t)0x8644, (int16_t)0xD3B2, (int16_t)0x87E9, (int16_t)0x658D, (int16_t)0xB214, (int16_t)0x6485, (int16_t)0xB0C2,
    (int16_t)0x2121, (int16_t)0x845D, (int16_t)0x1DE2, (int16_t)0x8389, (int16_t)0xCF04, (int16_t)0x89BE, (int16_t)0xCA69, (int16_t)0x8BC2,
    (int16_t)0x637A, (int16_t)0xAF72, (int16_t)0x6269, (int16_t)0xAE27, (int16_t)0x1A9D, (int16_t)0x82CC, (int16_t)0x1753, (int16_t)0x8225,
    (int16_t)0xC5E4, (int16_t)0x8DF3, (int16_t)0xC175, (int16_t)0x9052, (int16_t)0x6155, (int16_t)0xACDF, (int16_t)0x603C, (int16_t)0xAB9B,
    (int16_t)0x1406, (int16_t)0x8193, (int16_t)0x10B5, (int16_t)0x8118, (int16_t)0xBD1F, (int16_t)0x92DD, (int16_t)0xB8E3, (int16_t)0x9592,
    (int16_t)0x5F1F, (int16_t)0xAA5A, (int16_t)0x5DFE, (int16_t)0xA91D, (int16_t)0x0D61, (int16_t)0x80B4, (int16_t)0x0A0B, (int16_t)0x8065,
    (int16_t)0xB4C3, (int16_t)0x9872, (int16_t)0xB0C2, (int16_t)0x9B7B, (int16_t)0x5CD9, (int16_t)0xA7E4, (int16_t)0x5BB0, (int16_t)0xA6AF,
    (int16_t)0x06B3, (int16_t)0x802D, (int16_t)0x035A, (int16_t)0x800B, (int16_t)0xACDF, (int16_t)0x9EAB, (int16_t)0xA91D, (int16_t)0xA202,
    (int16_t)0x5A82, (int16_t)0xA57E, (int16_t)0x5951, (int16_t)0xA450, (int16_t)0x0000, (int16_t)0x8000, (int16_t)0xFCA6, (int16_t)0x800B,
    (int16_t)0xA57E, (int16_t)0xA57E, (int16_t)0xA202, (int16_t)0xA91D, (int16_t)0x581C, (int16_t)0xA327, (int16_t)0x56E3, (int16_t)0xA202,
    (int16_t)0xF94D, (int16_t)0x802D, (int16_t)0xF5F5, (int16_t)0x8065, (int16_t)0x9EAB, (int16_t)0xACDF, (int16_t)0x9B7B, (int16_t)0xB0C2,
    (int16_t)0x55A6, (int16_t)0xA0E1, (int16_t)0x5465, (int16_t)0x9FC4, (int16_t)0xF29F, (int16_t)0x80B4, (int16_t)0xEF4B, (int16_t)0x8118,
    (int16_t)0x9872, (int16_t)0xB4C3, (int16_t)0x9592, (int16_t)0xB8E3, (int16_t)0x5321, (int16_t)0x9EAB, (int16_t)0x51D9, (int16_t)0x9D97,
    (int16_t)0xEBFA, (int16_t)0x8193, (int16_t)0xE8AD, (int16_t)0x8225, (int16_t)0x92DD, (int16_t)0xBD1F, (int16_t)0x9052, (int16_t)0xC175,
    (int16_t)0x508E, (int16_t)0x9C86, (int16_t)0x4F3E, (int16_t)0x9B7B, (int16_t)0xE563, (int16_t)0x82CC, (int16_t)0xE21E, (int16_t)0x8389,
    (int16_t)0x8DF3, (int16_t)0xC5E4, (int16_t)0x8BC2, (int16_t)0xCA69, (int16_t)0x4DEC, (int16_t)0x9A73, (int16_t)0x4C96, (int16_t)0x9971,
    (int16_t)0xDEDF, (int16_t)0x845D, (int16_t)0xDBA5, (int16_t)0x8545, (int16_t)0x89BE, (int16_t)0xCF04, (int16_t)0x87E9, (int16_t)0xD3B2,
    (int16_t)0x4B3D, (int16_t)0x9872, (int16_t)0x49E0, (int16_t)0x9778, (int16_t)0xD872, (int16_t)0x8644, (int16_t)0xD546, (int16_t)0x8758,
    (int16_t)0x8644, (int16_t)0xD872, (int16_t)0x84CE, (int16_t)0xDD41, (int16_t)0x4880, (int16_t)0x9683, (int16_t)0x471D, (int16_t)0x9592,
    (int16_t)0xD221, (int16_t)0x8880, (int16_t)0xCF04, (int16_t)0x89BE, (int16_t)0x8389, (int16_t)0xE21E, (int16_t)0x8276, (int16_t)0xE707,
    (int16_t)0x45B7, (int16_t)0x94A6, (int16_t)0x444D, (int16_t)0x93BF, (int16_t)0xCBF0, (int16_t)0x8B11, (int16_t)0xC8E5, (int16_t)0x8C78,
    (int16_t)0x8193, (int16_t)0xEBFA, (int16_t)0x80E3, (int16_t)0xF0F5, (int16_t)0x42E1, (int16_t)0x92DD, (int16_t)0x4172, (int16_t)0x91FF,
    (int16_t)0xC5E4, (int16_t)0x8DF3, (int16_t)0xC2EC, (int16_t)0x8F83, (int16_t)0x8065, (int16_t)0xF5F5, (int16_t)0x8019, (int16_t)0xFAFA,
    (int16_t)0x4000, (int16_t)0x9126, (int16_t)0x3E8B, (int16_t)0x9052, (int16_t)0xC000, (int16_t)0x9126, (int16_t)0xBD1F, (int16_t)0x92DD,
    (int16_t)0x8000, (int16_t)0x0000, (int16_t)0x8019, (int16_t)0x0506, (int16_t)0x3D14, (int16_t)0x8F83, (int16_t)0x3B99, (int16_t)0x8EB9,
    (int16_t)0xBA49, (int16_t)0x94A6, (int16_t)0xB780, (int16_t)0x9683, (int16_t)0x8065, (int16_t)0x0A0B, (int16_t)0x80E3, (int16_t)0x0F0B,
    (int16_t)0x3A1C, (int16_t)0x8DF3, (int16_t)0x389D, (int16_t)0x8D33, (int16_t)0xB4C3, (int16_t)0x9872, (int16_t)0xB214, (int16_t)0x9A73,
    (int16_t)0x8193, (int16_t)0x1406, (int16_t)0x8276, (int16_t)0x18F9, (int16_t)0x371B, (int16_t)0x8C78, (int16_t)0x3597, (int16_t)0x8BC2,
    (int16_t)0xAF72, (int16_t)0x9C86, (int16_t)0xACDF, (int16_t)0x9EAB, (int16_t)0x8389, (int16_t)0x1DE2, (int16_t)0x84CE, (int16_t)0x22BF,
    (int16_t)0x3410, (int16_t)0x8B11, (int16_t)0x3287, (int16_t)0x8A65, (int16_t)0xAA5A, (int16_t)0xA0E1, (int16_t)0xA7E4, (int16_t)0xA327,
    (int16_t)0x8644, (int16_t)0x278E, (int16_t)0x87E9, (int16_t)0x2C4E, (int16_t)0x30FC, (int16_t)0x89BE, (int16_t)0x2F6E, (int16_t)0x891D,
    (int16_t)0xA57E, (int16_t)0xA57E, (int16_t)0xA327, (int16_t)0xA7E4, (int16_t)0x89BE, (int16_t)0x30FC, (int16_t)0x8BC2, (int16_t)0x3597,
    (int16_t)0x2DDF, (int16_t)0x8880, (int16_t)0x2C4E, (int16_t)0x87E9, (int16_t)0xA0E1, (int16_t)0xAA5A, (int16_t)0x9EAB, (int16_t)0xACDF,
    (int16_t)0x8DF3, (int16_t)0x3A1C, (int16_t)0x9052, (int16_t)0x3E8B, (int16_t)0x2ABA, (int16_t)0x8758, (int16_t)0x2925, (int16_t)0x86CB,
    (int16_t)0x9C86, (int16_t)0xAF72, (int16_t)0x9A73, (int16_t)0xB214, (int16_t)0x92DD, (int16_t)0x42E1, (int16_t)0x9592, (int16_t)0x471D,
    (int16_t)0x278E, (int16_t)0x8644, (int16_t)0x25F5, (int16_t)0x85C2, (int16_t)0x9872, (int16_t)0xB4C3, (int16_t)0x9683, (int16_t)0xB780,
    (int16_t)0x9872, (int16_t)0x4B3D, (int16_t)0x9B7B, (int16_t)0x4F3E, (int16_t)0x245B, (int16_t)0x8545, (int16_t)0x22BF, (int16_t)0x84CE,
    (int16_t)0x94A6, (int16_t)0xBA49, (int16_t)0x92DD, (int16_t)0xBD1F, (int16_t)0x9EAB, (int16_t)0x5321, (int16_t)0xA202, (int16_t)0x56E3,
    (int16_t)0x2121, (int16_t)0x845D, (int16_t)0x1F82, (int16_t)0x83F0, (int16_t)0x9126, (int16_t)0xC000, (int16_t)0x8F83, (int16_t)0xC2EC,
    (int16_t)0xA57E, (int16_t)0x5A82, (int16_t)0xA91D, (int16_t)0x5DFE, (int16_t)0x1DE2, (int16_t)0x8389, (int16_t)0x1C40, (int16_t)0x8328,
    (int16_t)0x8DF3, (int16_t)0xC5E4, (int16_t)0x8C78, (int16_t)0xC8E5, (int16_t)0xACDF, (int16_t)0x6155, (int16_t)0xB0C2, (int16_t)0x6485,
    (int16_t)0x1A9D, (int16_t)0x82CC, (int16_t)0x18F9, (int16_t)0x8276, (int16_t)0x8B11, (int16_t)0xCBF0, (int16_t)0x89BE, (int16_t)0xCF04,
    (int16_t)0xB4C3, (int16_t)0x678E, (int16_t)0xB8E3, (int16_t)0x6A6E, (int16_t)0x1753, (int16_t)0x8225, (int16_t)0x15AD, (int16_t)0x81D9,
    (int16_t)0x8880, (int16_t)0xD221, (int16_t)0x8758, (int16_t)0xD546, (int16_t)0xBD1F, (int16_t)0x6D23, (int16_t)0xC175, (int16_t)0x6FAE,
    (int16_t)0x1406, (int16_t)0x8193, (int16_t)0x125E, (int16_t)0x8153, (int16_t)0x8644, (int16_t)0xD872, (int16_t)0x8545, (int16_t)0xDBA5,
    (int16_t)0xC5E4, (int16_t)0x720D, (int16_t)0xCA69, (int16_t)0x743E, (int16_t)0x10B5, (int16_t)0x8118, (int16_t)0x0F0B, (int16_t)0x80E3,
    (int16_t)0x845D, (int16_t)0xDEDF, (int16_t)0x8389, (int16_t)0xE21E, (int16_t)0xCF04, (int16_t)0x7642, (int16_t)0xD3B2, (int16_t)0x7817,
    (int16_t)0x0D61, (int16_t)0x80B4, (int16_t)0x0BB6, (int16_t)0x8089, (int16_t)0x82CC, (int16_t)0xE563, (int16_t)0x8225, (int16_t)0xE8AD,
    (int16_t)0xD872, (int16_t)0x79BC, (int16_t)0xDD41, (int16_t)0x7B32, (int16_t)0x0A0B, (int16_t)0x8065, (int16_t)0x085F, (int16_t)0x8046,
    (int16_t)0x8193, (int16_t)0xEBFA, (int16_t)0x8118, (int16_t)0xEF4B, (int16_t)0xE21E, (int16_t)0x7C77, (int16_t)0xE707, (int16_t)0x7D8A,
    (int16_t)0x06B3, (int16_t)0x802D, (int16_t)0x0506, (int16_t)0x8019, (int16_t)0x80B4, (int16_t)0xF29F, (int16_t)0x8065, (int16_t)0xF5F5,
    (int16_t)0xEBFA, (int16_t)0x7E6D, (int16_t)0xF0F5, (int16_t)0x7F1D, (int16_t)0x035A, (int16_t)0x800B, (int16_t)0x01AD, (int16_t)0x8003,
    (int16_t)0x802D, (int16_t)0xF94D, (int16_t)0x800B, (int16_t)0xFCA6, (int16_t)0xF5F5, (int16_t)0x7F9B, (int16_t)0xFAFA, (int16_t)0x7FE7,
};

/********** Twiddles table N=480 stage 2 radix 3 ******************/
ALIGN(32) static const int16_t _fft480_tw2_[] =
{
    (int16_t)0x7FFF, (int16_t)0x0000, (int16_t)0x7FFF, (int16_t)0x0000, (int16_t)0x7FD3, (int16_t)0xF94D, (int16_t)0x7F4C, (int16_t)0xF29F,
    (int16_t)0x7F4C, (int16_t)0xF29F, (int16_t)0x7D34, (int16_t)0xE563, (int16_t)0x7E6D, (int16_t)0xEBFA, (int16_t)0x79BC, (int16_t)0xD872,
    (int16_t)0x7D34, (int16_t)0xE563, (int16_t)0x74EF, (int16_t)0xCBF0, (int16_t)0x7BA3, (int16_t)0xDEDF, (int16_t)0x6EDA, (int16_t)0xC000,
    (int16_t)0x79BC, (int16_t)0xD872, (int16_t)0x678E, (int16_t)0xB4C3, (int16_t)0x7780, (int16_t)0xD221, (int16_t)0x5F1F, (int16_t)0xAA5A,
    (int16_t)0x74EF, (int16_t)0xCBF0, (int16_t)0x55A6, (int16_t)0xA0E1, (int16_t)0x720D, (int16_t)0xC5E4, (int16_t)0x4B3D, (int16_t)0x9872,
    (int16_t)0x6EDA, (int16_t)0xC000, (int16_t)0x4000, (int16_t)0x9126, (int16_t)0x6B5A, (int16_t)0xBA49, (int16_t)0x3410, (int16_t)0x8B11,
    (int16_t)0x678E, (int16_t)0xB4C3, (int16_t)0x278E, (int16_t)0x8644, (int16_t)0x637A, (int16_t)0xAF72, (int16_t)0x1A9D, (int16_t)0x82CC,
    (int16_t)0x5F1F, (int16_t)0xAA5A, (int16_t)0x0D61, (int16_t)0x80B4, (int16_t)0x5A82, (int16_t)0xA57E, (int16_t)0x0000, (int16_t)0x8000,
    (int16_t)0x55A6, (int16_t)0xA0E1, (int16_t)0xF29F, (int16_t)0x80B4, (int16_t)0x508E, (int16_t)0x9C86, (int16_t)0xE563, (int16_t)0x82CC,
    (int16_t)0x4B3D, (int16_t)0x9872, (int16_t)0xD872, (int16_t)0x8644, (int16_t)0x45B7, (int16_t)0x94A6, (int16_t)0xCBF0, (int16_t)0x8B11,
    (int16_t)0x4000, (int16_t)0x9126, (int16_t)0xC000, (int16_t)0x9126, (int16_t)0x3A1C, (int16_t)0x8DF3, (int16_t)0xB4C3, (int16_t)0x9872,
    (int16_t)0x3410, (int16_t)0x8B11, (int16_t)0xAA5A, (int16_t)0xA0E1, (int16_t)0x2DDF, (int16_t)0x8880, (int16_t)0xA0E1, (int16_t)0xAA5A,
    (int16_t)0x278E, (int16_t)0x8644, (int16_t)0x9872, (int16_t)0xB4C3, (int16_t)0x2121, (int16_t)0x845D, (int16_t)0x9126, (int16_t)0xC000,
    (int16_t)0x1A9D, (int16_t)0x82CC, (int16_t)0x8B11, (int16_t)0xCBF0, (int16_t)0x1406, (int16_t)0x8193, (int16_t)0x8644, (int16_t)0xD872,
    (int16_t)0x0D61, (int16_t)0x80B4, (int16_t)0x82CC, (int16_t)0xE563, (int16_t)0x06B3, (int16_t)0x802D, (int16_t)0x80B4, (int16_t)0xF29F,
    (int16_t)0x0000, (int16_t)0x8000, (int16_t)0x8000, (int16_t)0x0000, (int16_t)0xF94D, (int16_t)0x802D, (int16_t)0x80B4, (int16_t)0x0D61,
    (int16_t)0xF29F, (int16_t)0x80B4, (int16_t)0x82CC, (int16_t)0x1A9D, (int16_t)0xEBFA, (int16_t)0x8193, (int16_t)0x8644, (int16_t)0x278E,
    (int16_t)0xE563, (int16_t)0x82CC, (int16_t)0x8B11, (int16_t)0x3410, (int16_t)0xDEDF, (int16_t)0x845D, (int16_t)0x9126, (int16_t)0x4000,
    (int16_t)0xD872, (int16_t)0x8644, (int16_t)0x9872, (int16_t)0x4B3D, (int16_t)0xD221, (int16_t)0x8880, (int16_t)0xA0E1, (int16_t)0x55A6,
    (int16_t)0xCBF0, (int16_t)0x8B11, (int16_t)0xAA5A, (int16_t)0x5F1F, (int16_t)0xC5E4, (int16_t)0x8DF3, (int16_t)0xB4C3, (int16_t)0x678E,
};

/********** Twiddles table N=480 stage 3 radix 5 ******************/
ALIGN(32) static const int16_t _fft480_tw3_[] =
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
    1, 1, 1, 5
}; 
static const cint16ptr_t_fft tw_tab[] = 
{
    _fft480_tw1_, _fft480_tw2_, _fft480_tw3_, NULL
};

static const fn_fft_stage fft32x16_stg_tab_s2[] = 
{
    (fn_fft_stage)fft_32x16_stage_first_scl2_DFT4,
    (fn_fft_stage)fft_32x16_stage_inner_scl2_DFT3,
    (fn_fft_stage)fft_32x16_stage_inner_scl2_DFT5,
    (fn_fft_stage)fft_32x16_stage_last_scl2_DFT8
};
static const fn_fft_stage fft32x16_stg_tab_s3[] = 
{
    (fn_fft_stage)fft_32x16_stage_first_scl3_DFT4,
    (fn_fft_stage)fft_32x16_stage_inner_scl3_DFT3,
    (fn_fft_stage)fft_32x16_stage_inner_scl3_DFT5,
    (fn_fft_stage)fft_32x16_stage_last_scl3_DFT8
};
static const fn_fft_stage ifft32x16_stg_tab_s2[] =
{
    (fn_fft_stage)ifft_32x16_stage_first_scl2_DFT4,
    (fn_fft_stage)fft_32x16_stage_inner_scl2_DFT3,
    (fn_fft_stage)fft_32x16_stage_inner_scl2_DFT5,
    (fn_fft_stage)ifft_32x16_stage_last_scl2_DFT8
};
static const fn_fft_stage ifft32x16_stg_tab_s3[] =
{
    (fn_fft_stage)ifft_32x16_stage_first_scl3_DFT4,
    (fn_fft_stage)fft_32x16_stage_inner_scl3_DFT3,
    (fn_fft_stage)fft_32x16_stage_inner_scl3_DFT5,
    (fn_fft_stage)ifft_32x16_stage_last_scl3_DFT8
};

static const fft_cplx_x16_descr_t __cfft_descr = 
{
    N, 
    tw_step_tab,
    tw_tab,
    fft32x16_stg_tab_s2,
    NULL, 
    fft32x16_stg_tab_s3,
    NULL
};     
static const fft_cplx_x16_descr_t __cifft_descr =
{
    N, 
    tw_step_tab,
    tw_tab,
    ifft32x16_stg_tab_s2,
    NULL, 
    ifft32x16_stg_tab_s3,
    NULL
};     
const fft_handle_t cnfft32x16_480 = (const fft_handle_t)&__cfft_descr;
const fft_handle_t cinfft32x16_480 = (const fft_handle_t)&__cifft_descr;
