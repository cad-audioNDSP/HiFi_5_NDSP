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

/* Twiddles tables for fft_cplx32x16, ifft_cplx32x16, fft_cplx16x16, ifft_cplx16x16, N=320 */
#define N 320

/********** Twiddles table N=320 stage 1 radix 4 ******************/
ALIGN(32) static const int16_t __fft16_tw1[] = 
{
    (int16_t)0x7FFF, (int16_t)0x0000, (int16_t)0x7FFA, (int16_t)0xFD7D, (int16_t)0x7FFF, (int16_t)0x0000, (int16_t)0x7FE7, (int16_t)0xFAFA,
    (int16_t)0x7FFF, (int16_t)0x0000, (int16_t)0x7FC7, (int16_t)0xF877, (int16_t)0x7FE7, (int16_t)0xFAFA, (int16_t)0x7FC7, (int16_t)0xF877,
    (int16_t)0x7F9B, (int16_t)0xF5F5, (int16_t)0x7F1D, (int16_t)0xF0F5, (int16_t)0x7F1D, (int16_t)0xF0F5, (int16_t)0x7E02, (int16_t)0xE980,
    (int16_t)0x7F9B, (int16_t)0xF5F5, (int16_t)0x7F62, (int16_t)0xF374, (int16_t)0x7E6D, (int16_t)0xEBFA, (int16_t)0x7D8A, (int16_t)0xE707,
    (int16_t)0x7C77, (int16_t)0xE21E, (int16_t)0x7A7D, (int16_t)0xDAD8, (int16_t)0x7F1D, (int16_t)0xF0F5, (int16_t)0x7ECB, (int16_t)0xEE76,
    (int16_t)0x7C77, (int16_t)0xE21E, (int16_t)0x7B32, (int16_t)0xDD41, (int16_t)0x7817, (int16_t)0xD3B2, (int16_t)0x7546, (int16_t)0xCCB4,
    (int16_t)0x7E6D, (int16_t)0xEBFA, (int16_t)0x7E02, (int16_t)0xE980, (int16_t)0x79BC, (int16_t)0xD872, (int16_t)0x7817, (int16_t)0xD3B2,
    (int16_t)0x720D, (int16_t)0xC5E4, (int16_t)0x6E6E, (int16_t)0xBF47, (int16_t)0x7D8A, (int16_t)0xE707, (int16_t)0x7D07, (int16_t)0xE492,
    (int16_t)0x7642, (int16_t)0xCF04, (int16_t)0x743E, (int16_t)0xCA69, (int16_t)0x6A6E, (int16_t)0xB8E3, (int16_t)0x660F, (int16_t)0xB2BF,
    (int16_t)0x7C77, (int16_t)0xE21E, (int16_t)0x7BDA, (int16_t)0xDFAE, (int16_t)0x720D, (int16_t)0xC5E4, (int16_t)0x6FAE, (int16_t)0xC175,
    (int16_t)0x6155, (int16_t)0xACDF, (int16_t)0x5C45, (int16_t)0xA749, (int16_t)0x7B32, (int16_t)0xDD41, (int16_t)0x7A7D, (int16_t)0xDAD8,
    (int16_t)0x6D23, (int16_t)0xBD1F, (int16_t)0x6A6E, (int16_t)0xB8E3, (int16_t)0x56E3, (int16_t)0xA202, (int16_t)0x5134, (int16_t)0x9D0E,
    (int16_t)0x79BC, (int16_t)0xD872, (int16_t)0x78EF, (int16_t)0xD610, (int16_t)0x678E, (int16_t)0xB4C3, (int16_t)0x6485, (int16_t)0xB0C2,
    (int16_t)0x4B3D, (int16_t)0x9872, (int16_t)0x4502, (int16_t)0x9432, (int16_t)0x7817, (int16_t)0xD3B2, (int16_t)0x7732, (int16_t)0xD159,
    (int16_t)0x6155, (int16_t)0xACDF, (int16_t)0x5DFE, (int16_t)0xA91D, (int16_t)0x3E8B, (int16_t)0x9052, (int16_t)0x37DC, (int16_t)0x8CD5,
    (int16_t)0x7642, (int16_t)0xCF04, (int16_t)0x7546, (int16_t)0xCCB4, (int16_t)0x5A82, (int16_t)0xA57E, (int16_t)0x56E3, (int16_t)0xA202,
    (int16_t)0x30FC, (int16_t)0x89BE, (int16_t)0x29F0, (int16_t)0x8711, (int16_t)0x743E, (int16_t)0xCA69, (int16_t)0x732B, (int16_t)0xC824,
    (int16_t)0x5321, (int16_t)0x9EAB, (int16_t)0x4F3E, (int16_t)0x9B7B, (int16_t)0x22BF, (int16_t)0x84CE, (int16_t)0x1B6E, (int16_t)0x82F9,
    (int16_t)0x720D, (int16_t)0xC5E4, (int16_t)0x70E3, (int16_t)0xC3A9, (int16_t)0x4B3D, (int16_t)0x9872, (int16_t)0x471D, (int16_t)0x9592,
    (int16_t)0x1406, (int16_t)0x8193, (int16_t)0x0C8C, (int16_t)0x809E, (int16_t)0x6FAE, (int16_t)0xC175, (int16_t)0x6E6E, (int16_t)0xBF47,
    (int16_t)0x42E1, (int16_t)0x92DD, (int16_t)0x3E8B, (int16_t)0x9052, (int16_t)0x0506, (int16_t)0x8019, (int16_t)0xFD7D, (int16_t)0x8006,
    (int16_t)0x6D23, (int16_t)0xBD1F, (int16_t)0x6BCE, (int16_t)0xBAFE, (int16_t)0x3A1C, (int16_t)0x8DF3, (int16_t)0x3597, (int16_t)0x8BC2,
    (int16_t)0xF5F5, (int16_t)0x8065, (int16_t)0xEE76, (int16_t)0x8135, (int16_t)0x6A6E, (int16_t)0xB8E3, (int16_t)0x6903, (int16_t)0xB6D0,
    (int16_t)0x30FC, (int16_t)0x89BE, (int16_t)0x2C4E, (int16_t)0x87E9, (int16_t)0xE707, (int16_t)0x8276, (int16_t)0xDFAE, (int16_t)0x8426,
    (int16_t)0x678E, (int16_t)0xB4C3, (int16_t)0x660F, (int16_t)0xB2BF, (int16_t)0x278E, (int16_t)0x8644, (int16_t)0x22BF, (int16_t)0x84CE,
    (int16_t)0xD872, (int16_t)0x8644, (int16_t)0xD159, (int16_t)0x88CE, (int16_t)0x6485, (int16_t)0xB0C2, (int16_t)0x62F2, (int16_t)0xAECC,
    (int16_t)0x1DE2, (int16_t)0x8389, (int16_t)0x18F9, (int16_t)0x8276, (int16_t)0xCA69, (int16_t)0x8BC2, (int16_t)0xC3A9, (int16_t)0x8F1D,
    (int16_t)0x6155, (int16_t)0xACDF, (int16_t)0x5FAE, (int16_t)0xAAFA, (int16_t)0x1406, (int16_t)0x8193, (int16_t)0x0F0B, (int16_t)0x80E3,
    (int16_t)0xBD1F, (int16_t)0x92DD, (int16_t)0xB6D0, (int16_t)0x96FD, (int16_t)0x5DFE, (int16_t)0xA91D, (int16_t)0x5C45, (int16_t)0xA749,
    (int16_t)0x0A0B, (int16_t)0x8065, (int16_t)0x0506, (int16_t)0x8019, (int16_t)0xB0C2, (int16_t)0x9B7B, (int16_t)0xAAFA, (int16_t)0xA052,
    (int16_t)0x5A82, (int16_t)0xA57E, (int16_t)0x58B7, (int16_t)0xA3BB, (int16_t)0x0000, (int16_t)0x8000, (int16_t)0xFAFA, (int16_t)0x8019,
    (int16_t)0xA57E, (int16_t)0xA57E, (int16_t)0xA052, (int16_t)0xAAFA, (int16_t)0x56E3, (int16_t)0xA202, (int16_t)0x5506, (int16_t)0xA052,
    (int16_t)0xF5F5, (int16_t)0x8065, (int16_t)0xF0F5, (int16_t)0x80E3, (int16_t)0x9B7B, (int16_t)0xB0C2, (int16_t)0x96FD, (int16_t)0xB6D0,
    (int16_t)0x5321, (int16_t)0x9EAB, (int16_t)0x5134, (int16_t)0x9D0E, (int16_t)0xEBFA, (int16_t)0x8193, (int16_t)0xE707, (int16_t)0x8276,
    (int16_t)0x92DD, (int16_t)0xBD1F, (int16_t)0x8F1D, (int16_t)0xC3A9, (int16_t)0x4F3E, (int16_t)0x9B7B, (int16_t)0x4D41, (int16_t)0x99F1,
    (int16_t)0xE21E, (int16_t)0x8389, (int16_t)0xDD41, (int16_t)0x84CE, (int16_t)0x8BC2, (int16_t)0xCA69, (int16_t)0x88CE, (int16_t)0xD159,
    (int16_t)0x4B3D, (int16_t)0x9872, (int16_t)0x4930, (int16_t)0x96FD, (int16_t)0xD872, (int16_t)0x8644, (int16_t)0xD3B2, (int16_t)0x87E9,
    (int16_t)0x8644, (int16_t)0xD872, (int16_t)0x8426, (int16_t)0xDFAE, (int16_t)0x471D, (int16_t)0x9592, (int16_t)0x4502, (int16_t)0x9432,
    (int16_t)0xCF04, (int16_t)0x89BE, (int16_t)0xCA69, (int16_t)0x8BC2, (int16_t)0x8276, (int16_t)0xE707, (int16_t)0x8135, (int16_t)0xEE76,
    (int16_t)0x42E1, (int16_t)0x92DD, (int16_t)0x40B9, (int16_t)0x9192, (int16_t)0xC5E4, (int16_t)0x8DF3, (int16_t)0xC175, (int16_t)0x9052,
    (int16_t)0x8065, (int16_t)0xF5F5, (int16_t)0x8006, (int16_t)0xFD7D, (int16_t)0x3E8B, (int16_t)0x9052, (int16_t)0x3C57, (int16_t)0x8F1D,
    (int16_t)0xBD1F, (int16_t)0x92DD, (int16_t)0xB8E3, (int16_t)0x9592, (int16_t)0x8019, (int16_t)0x0506, (int16_t)0x809E, (int16_t)0x0C8C,
    (int16_t)0x3A1C, (int16_t)0x8DF3, (int16_t)0x37DC, (int16_t)0x8CD5, (int16_t)0xB4C3, (int16_t)0x9872, (int16_t)0xB0C2, (int16_t)0x9B7B,
    (int16_t)0x8193, (int16_t)0x1406, (int16_t)0x82F9, (int16_t)0x1B6E, (int16_t)0x3597, (int16_t)0x8BC2, (int16_t)0x334C, (int16_t)0x8ABA,
    (int16_t)0xACDF, (int16_t)0x9EAB, (int16_t)0xA91D, (int16_t)0xA202, (int16_t)0x84CE, (int16_t)0x22BF, (int16_t)0x8711, (int16_t)0x29F0,
    (int16_t)0x30FC, (int16_t)0x89BE, (int16_t)0x2EA7, (int16_t)0x88CE, (int16_t)0xA57E, (int16_t)0xA57E, (int16_t)0xA202, (int16_t)0xA91D,
    (int16_t)0x89BE, (int16_t)0x30FC, (int16_t)0x8CD5, (int16_t)0x37DC, (int16_t)0x2C4E, (int16_t)0x87E9, (int16_t)0x29F0, (int16_t)0x8711,
    (int16_t)0x9EAB, (int16_t)0xACDF, (int16_t)0x9B7B, (int16_t)0xB0C2, (int16_t)0x9052, (int16_t)0x3E8B, (int16_t)0x9432, (int16_t)0x4502,
    (int16_t)0x278E, (int16_t)0x8644, (int16_t)0x2528, (int16_t)0x8583, (int16_t)0x9872, (int16_t)0xB4C3, (int16_t)0x9592, (int16_t)0xB8E3,
    (int16_t)0x9872, (int16_t)0x4B3D, (int16_t)0x9D0E, (int16_t)0x5134, (int16_t)0x22BF, (int16_t)0x84CE, (int16_t)0x2052, (int16_t)0x8426,
    (int16_t)0x92DD, (int16_t)0xBD1F, (int16_t)0x9052, (int16_t)0xC175, (int16_t)0xA202, (int16_t)0x56E3, (int16_t)0xA749, (int16_t)0x5C45,
    (int16_t)0x1DE2, (int16_t)0x8389, (int16_t)0x1B6E, (int16_t)0x82F9, (int16_t)0x8DF3, (int16_t)0xC5E4, (int16_t)0x8BC2, (int16_t)0xCA69,
    (int16_t)0xACDF, (int16_t)0x6155, (int16_t)0xB2BF, (int16_t)0x660F, (int16_t)0x18F9, (int16_t)0x8276, (int16_t)0x1680, (int16_t)0x81FE,
    (int16_t)0x89BE, (int16_t)0xCF04, (int16_t)0x87E9, (int16_t)0xD3B2, (int16_t)0xB8E3, (int16_t)0x6A6E, (int16_t)0xBF47, (int16_t)0x6E6E,
    (int16_t)0x1406, (int16_t)0x8193, (int16_t)0x118A, (int16_t)0x8135, (int16_t)0x8644, (int16_t)0xD872, (int16_t)0x84CE, (int16_t)0xDD41,
    (int16_t)0xC5E4, (int16_t)0x720D, (int16_t)0xCCB4, (int16_t)0x7546, (int16_t)0x0F0B, (int16_t)0x80E3, (int16_t)0x0C8C, (int16_t)0x809E,
    (int16_t)0x8389, (int16_t)0xE21E, (int16_t)0x8276, (int16_t)0xE707, (int16_t)0xD3B2, (int16_t)0x7817, (int16_t)0xDAD8, (int16_t)0x7A7D,
    (int16_t)0x0A0B, (int16_t)0x8065, (int16_t)0x0789, (int16_t)0x8039, (int16_t)0x8193, (int16_t)0xEBFA, (int16_t)0x80E3, (int16_t)0xF0F5,
    (int16_t)0xE21E, (int16_t)0x7C77, (int16_t)0xE980, (int16_t)0x7E02, (int16_t)0x0506, (int16_t)0x8019, (int16_t)0x0283, (int16_t)0x8006,
    (int16_t)0x8065, (int16_t)0xF5F5, (int16_t)0x8019, (int16_t)0xFAFA, (int16_t)0xF0F5, (int16_t)0x7F1D, (int16_t)0xF877, (int16_t)0x7FC7
};

/********** Twiddles table N=320 stage 2 radix 4 ******************/
ALIGN(32) static const int16_t __fft16_tw2[] = 
{
    (int16_t)0x7FFF, (int16_t)0x0000, (int16_t)0x7FFF, (int16_t)0x0000, (int16_t)0x7FFF, (int16_t)0x0000, (int16_t)0x7F9B, (int16_t)0xF5F5,
    (int16_t)0x7E6D, (int16_t)0xEBFA, (int16_t)0x7C77, (int16_t)0xE21E, (int16_t)0x7E6D, (int16_t)0xEBFA, (int16_t)0x79BC, (int16_t)0xD872,
    (int16_t)0x720D, (int16_t)0xC5E4, (int16_t)0x7C77, (int16_t)0xE21E, (int16_t)0x720D, (int16_t)0xC5E4, (int16_t)0x6155, (int16_t)0xACDF,
    (int16_t)0x79BC, (int16_t)0xD872, (int16_t)0x678E, (int16_t)0xB4C3, (int16_t)0x4B3D, (int16_t)0x9872, (int16_t)0x7642, (int16_t)0xCF04,
    (int16_t)0x5A82, (int16_t)0xA57E, (int16_t)0x30FC, (int16_t)0x89BE, (int16_t)0x720D, (int16_t)0xC5E4, (int16_t)0x4B3D, (int16_t)0x9872,
    (int16_t)0x1406, (int16_t)0x8193, (int16_t)0x6D23, (int16_t)0xBD1F, (int16_t)0x3A1C, (int16_t)0x8DF3, (int16_t)0xF5F5, (int16_t)0x8065,
    (int16_t)0x678E, (int16_t)0xB4C3, (int16_t)0x278E, (int16_t)0x8644, (int16_t)0xD872, (int16_t)0x8644, (int16_t)0x6155, (int16_t)0xACDF,
    (int16_t)0x1406, (int16_t)0x8193, (int16_t)0xBD1F, (int16_t)0x92DD, (int16_t)0x5A82, (int16_t)0xA57E, (int16_t)0x0000, (int16_t)0x8000,
    (int16_t)0xA57E, (int16_t)0xA57E, (int16_t)0x5321, (int16_t)0x9EAB, (int16_t)0xEBFA, (int16_t)0x8193, (int16_t)0x92DD, (int16_t)0xBD1F,
    (int16_t)0x4B3D, (int16_t)0x9872, (int16_t)0xD872, (int16_t)0x8644, (int16_t)0x8644, (int16_t)0xD872, (int16_t)0x42E1, (int16_t)0x92DD,
    (int16_t)0xC5E4, (int16_t)0x8DF3, (int16_t)0x8065, (int16_t)0xF5F5, (int16_t)0x3A1C, (int16_t)0x8DF3, (int16_t)0xB4C3, (int16_t)0x9872,
    (int16_t)0x8193, (int16_t)0x1406, (int16_t)0x30FC, (int16_t)0x89BE, (int16_t)0xA57E, (int16_t)0xA57E, (int16_t)0x89BE, (int16_t)0x30FC,
    (int16_t)0x278E, (int16_t)0x8644, (int16_t)0x9872, (int16_t)0xB4C3, (int16_t)0x9872, (int16_t)0x4B3D, (int16_t)0x1DE2, (int16_t)0x8389,
    (int16_t)0x8DF3, (int16_t)0xC5E4, (int16_t)0xACDF, (int16_t)0x6155, (int16_t)0x1406, (int16_t)0x8193, (int16_t)0x8644, (int16_t)0xD872,
    (int16_t)0xC5E4, (int16_t)0x720D, (int16_t)0x0A0B, (int16_t)0x8065, (int16_t)0x8193, (int16_t)0xEBFA, (int16_t)0xE21E, (int16_t)0x7C77,
};

/********** Twiddles table N=320 stage 3 radix 4 ******************/
ALIGN(32) static const int16_t __fft16_tw3[] = 
{
    (int16_t)0x7FFF, (int16_t)0x0000, (int16_t)0x7FFF, (int16_t)0x0000, (int16_t)0x7FFF, (int16_t)0x0000,
    (int16_t)0x79BC, (int16_t)0xD872, (int16_t)0x678E, (int16_t)0xB4C3, (int16_t)0x4B3D, (int16_t)0x9872,
    (int16_t)0x678E, (int16_t)0xB4C3, (int16_t)0x278E, (int16_t)0x8644, (int16_t)0xD872, (int16_t)0x8644,
    (int16_t)0x4B3D, (int16_t)0x9872, (int16_t)0xD872, (int16_t)0x8644, (int16_t)0x8644, (int16_t)0xD872,
    (int16_t)0x278E, (int16_t)0x8644, (int16_t)0x9872, (int16_t)0xB4C3, (int16_t)0x9872, (int16_t)0x4B3D, 
};

static const int tw_step_tab[] =
{
    1, 1, 1, 0
}; 
static const cint16ptr_t_fft tw_tab[] = 
{
    __fft16_tw1, __fft16_tw2, __fft16_tw3, NULL
};

static const fn_fft_stage fft_stg_tab_s2[] = 
{
    (fn_fft_stage)fft_16x16_stage_first_scl2_DFT4,
    (fn_fft_stage)fft_16x16_stage_inner_scl2_DFT4x2_v4,
    (fn_fft_stage)fft_16x16_stage_inner_scl2_DFT4x2,
    (fn_fft_stage)fft_16x16_stage_last_scl2_DFT5
};
static const fn_fft_stage fft_stg_tab_s3[] = 
{
    (fn_fft_stage)fft_16x16_stage_first_scl3_DFT4,
    (fn_fft_stage)fft_16x16_stage_inner_scl3_DFT4x2_v4,
    (fn_fft_stage)fft_16x16_stage_inner_scl3_DFT4x2,
    (fn_fft_stage)fft_16x16_stage_last_scl3_DFT5
};
static const fn_fft_stage ifft_stg_tab_s2[] =
{
    (fn_fft_stage)ifft_16x16_stage_first_scl2_DFT4,
    (fn_fft_stage)fft_16x16_stage_inner_scl2_DFT4x2_v4,
    (fn_fft_stage)fft_16x16_stage_inner_scl2_DFT4x2,
    (fn_fft_stage)ifft_16x16_stage_last_scl2_DFT5
};
static const fn_fft_stage ifft_stg_tab_s3[] =
{
    (fn_fft_stage)ifft_16x16_stage_first_scl3_DFT4,
    (fn_fft_stage)fft_16x16_stage_inner_scl3_DFT4x2_v4,
    (fn_fft_stage)fft_16x16_stage_inner_scl3_DFT4x2,
    (fn_fft_stage)ifft_16x16_stage_last_scl3_DFT5
};

const fft_cplx_x16_descr_t __cfft_x16_descr320 =
{
    N, 
    tw_step_tab,
    tw_tab,
    NULL,
    fft_stg_tab_s2,
    NULL, 
    fft_stg_tab_s3
};     
const fft_cplx_x16_descr_t __cifft_x16_descr320 =
{
    N, 
    tw_step_tab,
    tw_tab,
    NULL,
    ifft_stg_tab_s2,
    NULL, 
    ifft_stg_tab_s3
};     
const fft_handle_t cnfft16_320 = (const fft_handle_t)&__cfft_x16_descr320;
const fft_handle_t cinfft16_320 = (const fft_handle_t)&__cifft_x16_descr320;
