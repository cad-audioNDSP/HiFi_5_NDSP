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
#include "fft_twiddles32x32.h"
#include "common.h"

/* Twiddles tables for fft_cplx32x32, ifft_cplx32x32 */

/****************** N=320 stage 1 radix 4 ******************/
ALIGN(32) static const int32_t __fft320_tw1[] =
{
    (int32_t)0x7FFFFFFF, (int32_t)0x00000000, (int32_t)0x7FFFFFFF, (int32_t)0x00000000, (int32_t)0x7FFFFFFF, (int32_t)0x00000000, (int32_t)0x7FF9AF04, (int32_t)0xFD7CA4A6,
    (int32_t)0x7FE6BCB0, (int32_t)0xFAF988CC, (int32_t)0x7FC72AE2, (int32_t)0xF876EBE8, (int32_t)0x7FE6BCB0, (int32_t)0xFAF988CC, (int32_t)0x7F9AFCB9, (int32_t)0xF5F50D67,
    (int32_t)0x7F1CDE01, (int32_t)0xF0F488D9, (int32_t)0x7FC72AE2, (int32_t)0xF876EBE8, (int32_t)0x7F1CDE01, (int32_t)0xF0F488D9, (int32_t)0x7E01B096, (int32_t)0xE97F81EB,
    (int32_t)0x7F9AFCB9, (int32_t)0xF5F50D67, (int32_t)0x7E6C9251, (int32_t)0xEBF9F498, (int32_t)0x7C769E18, (int32_t)0xE21E765A, (int32_t)0x7F62368F, (int32_t)0xF3742CA2,
    (int32_t)0x7D8A5F40, (int32_t)0xE70747C4, (int32_t)0x7A7D055B, (int32_t)0xDAD7F3A2, (int32_t)0x7F1CDE01, (int32_t)0xF0F488D9, (int32_t)0x7C769E18, (int32_t)0xE21E765A,
    (int32_t)0x7816A759, (int32_t)0xD3B26FB0, (int32_t)0x7ECAF9E5, (int32_t)0xEE76612D, (int32_t)0x7B31BBB2, (int32_t)0xDD417079, (int32_t)0x7545A5A0, (int32_t)0xCCB44322,
    (int32_t)0x7E6C9251, (int32_t)0xEBF9F498, (int32_t)0x79BC384D, (int32_t)0xD8722192, (int32_t)0x720C8075, (int32_t)0xC5E3A3A9, (int32_t)0x7E01B096, (int32_t)0xE97F81EB,
    (int32_t)0x7816A759, (int32_t)0xD3B26FB0, (int32_t)0x6E6E1492, (int32_t)0xBF469E83, (int32_t)0x7D8A5F40, (int32_t)0xE70747C4, (int32_t)0x7641AF3D, (int32_t)0xCF043AB3,
    (int32_t)0x6A6D98A4, (int32_t)0xB8E31319, (int32_t)0x7D06AA16, (int32_t)0xE4918486, (int32_t)0x743E0918, (int32_t)0xCA695B94, (int32_t)0x660E9A6A, (int32_t)0xB2BEADCC,
    (int32_t)0x7C769E18, (int32_t)0xE21E765A, (int32_t)0x720C8075, (int32_t)0xC5E3A3A9, (int32_t)0x6154FB91, (int32_t)0xACDEE2E8, (int32_t)0x7BDA497D, (int32_t)0xDFAE5B23,
    (int32_t)0x6FADF2FC, (int32_t)0xC174DBF2, (int32_t)0x5C44EE40, (int32_t)0xA748E9CE, (int32_t)0x7B31BBB2, (int32_t)0xDD417079, (int32_t)0x6D23501B, (int32_t)0xBD1EC45C,
    (int32_t)0x56E2F15D, (int32_t)0xA201B853, (int32_t)0x7A7D055B, (int32_t)0xDAD7F3A2, (int32_t)0x6A6D98A4, (int32_t)0xB8E31319, (int32_t)0x5133CC94, (int32_t)0x9D0DFE54,
    (int32_t)0x79BC384D, (int32_t)0xD8722192, (int32_t)0x678DDE6E, (int32_t)0xB4C373EE, (int32_t)0x4B3C8C12, (int32_t)0x98722192, (int32_t)0x78EF678F, (int32_t)0xD61036DB,
    (int32_t)0x648543E4, (int32_t)0xB0C1878B, (int32_t)0x45027C0C, (int32_t)0x943239C7, (int32_t)0x7816A759, (int32_t)0xD3B26FB0, (int32_t)0x6154FB91, (int32_t)0xACDEE2E8,
    (int32_t)0x3E8B240E, (int32_t)0x90520D04, (int32_t)0x77320D0D, (int32_t)0xD15907D9, (int32_t)0x5DFE47AD, (int32_t)0xA91D0EA3, (int32_t)0x37DC420C, (int32_t)0x8CD50C59,
    (int32_t)0x7641AF3D, (int32_t)0xCF043AB3, (int32_t)0x5A82799A, (int32_t)0xA57D8666, (int32_t)0x30FBC54D, (int32_t)0x89BE50C3, (int32_t)0x7545A5A0, (int32_t)0xCCB44322,
    (int32_t)0x56E2F15D, (int32_t)0xA201B853, (int32_t)0x29EFC925, (int32_t)0x87109871, (int32_t)0x743E0918, (int32_t)0xCA695B94, (int32_t)0x53211D18, (int32_t)0x9EAB046F,
    (int32_t)0x22BE8F87, (int32_t)0x84CE444E, (int32_t)0x732AF3A7, (int32_t)0xC823BDF4, (int32_t)0x4F3E7875, (int32_t)0x9B7ABC1C, (int32_t)0x1B6E7B7A, (int32_t)0x82F955EA,
    (int32_t)0x720C8075, (int32_t)0xC5E3A3A9, (int32_t)0x4B3C8C12, (int32_t)0x98722192, (int32_t)0x14060B68, (int32_t)0x81936DAF, (int32_t)0x70E2CBC6, (int32_t)0xC3A94590,
    (int32_t)0x471CECE7, (int32_t)0x9592675C, (int32_t)0x0C8BD35E, (int32_t)0x809DC971, (int32_t)0x6FADF2FC, (int32_t)0xC174DBF2, (int32_t)0x42E13BA4, (int32_t)0x92DCAFE5,
    (int32_t)0x05067734, (int32_t)0x80194350, (int32_t)0x6E6E1492, (int32_t)0xBF469E83, (int32_t)0x3E8B240E, (int32_t)0x90520D04, (int32_t)0xFD7CA4A6, (int32_t)0x800650FC,
    (int32_t)0x6D23501B, (int32_t)0xBD1EC45C, (int32_t)0x3A1C5C57, (int32_t)0x8DF37F8B, (int32_t)0xF5F50D67, (int32_t)0x80650347, (int32_t)0x6BCDC639, (int32_t)0xBAFD83F4,
    (int32_t)0x3596A46C, (int32_t)0x8BC1F6E8, (int32_t)0xEE76612D, (int32_t)0x8135061B, (int32_t)0x6A6D98A4, (int32_t)0xB8E31319, (int32_t)0x30FBC54D, (int32_t)0x89BE50C3,
    (int32_t)0xE70747C4, (int32_t)0x8275A0C0, (int32_t)0x6902EA1D, (int32_t)0xB6CFA6F1, (int32_t)0x2C4D9050, (int32_t)0x87E958A7, (int32_t)0xDFAE5B23, (int32_t)0x8425B683,
    (int32_t)0x678DDE6E, (int32_t)0xB4C373EE, (int32_t)0x278DDE6E, (int32_t)0x8643C7B3, (int32_t)0xD8722192, (int32_t)0x8643C7B3, (int32_t)0x660E9A6A, (int32_t)0xB2BEADCC,
    (int32_t)0x22BE8F87, (int32_t)0x84CE444E, (int32_t)0xD15907D9, (int32_t)0x88CDF2F3, (int32_t)0x648543E4, (int32_t)0xB0C1878B, (int32_t)0x1DE189A6, (int32_t)0x838961E8,
    (int32_t)0xCA695B94, (int32_t)0x8BC1F6E8, (int32_t)0x62F201AC, (int32_t)0xAECC336C, (int32_t)0x18F8B83C, (int32_t)0x8275A0C0, (int32_t)0xC3A94590, (int32_t)0x8F1D343A,
    (int32_t)0x6154FB91, (int32_t)0xACDEE2E8, (int32_t)0x14060B68, (int32_t)0x81936DAF, (int32_t)0xBD1EC45C, (int32_t)0x92DCAFE5, (int32_t)0x5FAE5A55, (int32_t)0xAAF9C6AF,
    (int32_t)0x0F0B7727, (int32_t)0x80E321FF, (int32_t)0xB6CFA6F1, (int32_t)0x96FD15E3, (int32_t)0x5DFE47AD, (int32_t)0xA91D0EA3, (int32_t)0x0A0AF299, (int32_t)0x80650347,
    (int32_t)0xB0C1878B, (int32_t)0x9B7ABC1C, (int32_t)0x5C44EE40, (int32_t)0xA748E9CE, (int32_t)0x05067734, (int32_t)0x80194350, (int32_t)0xAAF9C6AF, (int32_t)0xA051A5AB,
    (int32_t)0x5A82799A, (int32_t)0xA57D8666, (int32_t)0x00000000, (int32_t)0x80000000, (int32_t)0xA57D8666, (int32_t)0xA57D8666, (int32_t)0x58B71632, (int32_t)0xA3BB11C0,
    (int32_t)0xFAF988CC, (int32_t)0x80194350, (int32_t)0xA051A5AB, (int32_t)0xAAF9C6AF, (int32_t)0x56E2F15D, (int32_t)0xA201B853, (int32_t)0xF5F50D67, (int32_t)0x80650347,
    (int32_t)0x9B7ABC1C, (int32_t)0xB0C1878B, (int32_t)0x55063951, (int32_t)0xA051A5AB, (int32_t)0xF0F488D9, (int32_t)0x80E321FF, (int32_t)0x96FD15E3, (int32_t)0xB6CFA6F1,
    (int32_t)0x53211D18, (int32_t)0x9EAB046F, (int32_t)0xEBF9F498, (int32_t)0x81936DAF, (int32_t)0x92DCAFE5, (int32_t)0xBD1EC45C, (int32_t)0x5133CC94, (int32_t)0x9D0DFE54,
    (int32_t)0xE70747C4, (int32_t)0x8275A0C0, (int32_t)0x8F1D343A, (int32_t)0xC3A94590, (int32_t)0x4F3E7875, (int32_t)0x9B7ABC1C, (int32_t)0xE21E765A, (int32_t)0x838961E8,
    (int32_t)0x8BC1F6E8, (int32_t)0xCA695B94, (int32_t)0x4D415234, (int32_t)0x99F16596, (int32_t)0xDD417079, (int32_t)0x84CE444E, (int32_t)0x88CDF2F3, (int32_t)0xD15907D9,
    (int32_t)0x4B3C8C12, (int32_t)0x98722192, (int32_t)0xD8722192, (int32_t)0x8643C7B3, (int32_t)0x8643C7B3, (int32_t)0xD8722192, (int32_t)0x4930590F, (int32_t)0x96FD15E3,
    (int32_t)0xD3B26FB0, (int32_t)0x87E958A7, (int32_t)0x8425B683, (int32_t)0xDFAE5B23, (int32_t)0x471CECE7, (int32_t)0x9592675C, (int32_t)0xCF043AB3, (int32_t)0x89BE50C3,
    (int32_t)0x8275A0C0, (int32_t)0xE70747C4, (int32_t)0x45027C0C, (int32_t)0x943239C7, (int32_t)0xCA695B94, (int32_t)0x8BC1F6E8, (int32_t)0x8135061B, (int32_t)0xEE76612D,
    (int32_t)0x42E13BA4, (int32_t)0x92DCAFE5, (int32_t)0xC5E3A3A9, (int32_t)0x8DF37F8B, (int32_t)0x80650347, (int32_t)0xF5F50D67, (int32_t)0x40B9617D, (int32_t)0x9191EB6E,
    (int32_t)0xC174DBF2, (int32_t)0x90520D04, (int32_t)0x800650FC, (int32_t)0xFD7CA4A6, (int32_t)0x3E8B240E, (int32_t)0x90520D04, (int32_t)0xBD1EC45C, (int32_t)0x92DCAFE5,
    (int32_t)0x80194350, (int32_t)0x05067734, (int32_t)0x3C56BA70, (int32_t)0x8F1D343A, (int32_t)0xB8E31319, (int32_t)0x9592675C, (int32_t)0x809DC971, (int32_t)0x0C8BD35E,
    (int32_t)0x3A1C5C57, (int32_t)0x8DF37F8B, (int32_t)0xB4C373EE, (int32_t)0x98722192, (int32_t)0x81936DAF, (int32_t)0x14060B68, (int32_t)0x37DC420C, (int32_t)0x8CD50C59,
    (int32_t)0xB0C1878B, (int32_t)0x9B7ABC1C, (int32_t)0x82F955EA, (int32_t)0x1B6E7B7A, (int32_t)0x3596A46C, (int32_t)0x8BC1F6E8, (int32_t)0xACDEE2E8, (int32_t)0x9EAB046F,
    (int32_t)0x84CE444E, (int32_t)0x22BE8F87, (int32_t)0x334BBCDE, (int32_t)0x8ABA5A60, (int32_t)0xA91D0EA3, (int32_t)0xA201B853, (int32_t)0x87109871, (int32_t)0x29EFC925,
    (int32_t)0x30FBC54D, (int32_t)0x89BE50C3, (int32_t)0xA57D8666, (int32_t)0xA57D8666, (int32_t)0x89BE50C3, (int32_t)0x30FBC54D, (int32_t)0x2EA6F827, (int32_t)0x88CDF2F3,
    (int32_t)0xA201B853, (int32_t)0xA91D0EA3, (int32_t)0x8CD50C59, (int32_t)0x37DC420C, (int32_t)0x2C4D9050, (int32_t)0x87E958A7, (int32_t)0x9EAB046F, (int32_t)0xACDEE2E8,
    (int32_t)0x90520D04, (int32_t)0x3E8B240E, (int32_t)0x29EFC925, (int32_t)0x87109871, (int32_t)0x9B7ABC1C, (int32_t)0xB0C1878B, (int32_t)0x943239C7, (int32_t)0x45027C0C,
    (int32_t)0x278DDE6E, (int32_t)0x8643C7B3, (int32_t)0x98722192, (int32_t)0xB4C373EE, (int32_t)0x98722192, (int32_t)0x4B3C8C12, (int32_t)0x25280C5E, (int32_t)0x8582FAA5,
    (int32_t)0x9592675C, (int32_t)0xB8E31319, (int32_t)0x9D0DFE54, (int32_t)0x5133CC94, (int32_t)0x22BE8F87, (int32_t)0x84CE444E, (int32_t)0x92DCAFE5, (int32_t)0xBD1EC45C,
    (int32_t)0xA201B853, (int32_t)0x56E2F15D, (int32_t)0x2051A4DD, (int32_t)0x8425B683, (int32_t)0x90520D04, (int32_t)0xC174DBF2, (int32_t)0xA748E9CE, (int32_t)0x5C44EE40,
    (int32_t)0x1DE189A6, (int32_t)0x838961E8, (int32_t)0x8DF37F8B, (int32_t)0xC5E3A3A9, (int32_t)0xACDEE2E8, (int32_t)0x6154FB91, (int32_t)0x1B6E7B7A, (int32_t)0x82F955EA,
    (int32_t)0x8BC1F6E8, (int32_t)0xCA695B94, (int32_t)0xB2BEADCC, (int32_t)0x660E9A6A, (int32_t)0x18F8B83C, (int32_t)0x8275A0C0, (int32_t)0x89BE50C3, (int32_t)0xCF043AB3,
    (int32_t)0xB8E31319, (int32_t)0x6A6D98A4, (int32_t)0x16807E15, (int32_t)0x81FE4F6A, (int32_t)0x87E958A7, (int32_t)0xD3B26FB0, (int32_t)0xBF469E83, (int32_t)0x6E6E1492,
    (int32_t)0x14060B68, (int32_t)0x81936DAF, (int32_t)0x8643C7B3, (int32_t)0xD8722192, (int32_t)0xC5E3A3A9, (int32_t)0x720C8075, (int32_t)0x11899ED3, (int32_t)0x8135061B,
    (int32_t)0x84CE444E, (int32_t)0xDD417079, (int32_t)0xCCB44322, (int32_t)0x7545A5A0, (int32_t)0x0F0B7727, (int32_t)0x80E321FF, (int32_t)0x838961E8, (int32_t)0xE21E765A,
    (int32_t)0xD3B26FB0, (int32_t)0x7816A759, (int32_t)0x0C8BD35E, (int32_t)0x809DC971, (int32_t)0x8275A0C0, (int32_t)0xE70747C4, (int32_t)0xDAD7F3A2, (int32_t)0x7A7D055B,
    (int32_t)0x0A0AF299, (int32_t)0x80650347, (int32_t)0x81936DAF, (int32_t)0xEBF9F498, (int32_t)0xE21E765A, (int32_t)0x7C769E18, (int32_t)0x07891418, (int32_t)0x8038D51E,
    (int32_t)0x80E321FF, (int32_t)0xF0F488D9, (int32_t)0xE97F81EB, (int32_t)0x7E01B096, (int32_t)0x05067734, (int32_t)0x80194350, (int32_t)0x80650347, (int32_t)0xF5F50D67,
    (int32_t)0xF0F488D9, (int32_t)0x7F1CDE01, (int32_t)0x02835B5A, (int32_t)0x800650FC, (int32_t)0x80194350, (int32_t)0xFAF988CC, (int32_t)0xF876EBE8, (int32_t)0x7FC72AE2,
};

/****************** N=320 stage 2 radix 4 ******************/
ALIGN(32) static const int32_t __fft320_tw2[] =
{
    (int32_t)0x7FFFFFFF, (int32_t)0x00000000, (int32_t)0x7FFFFFFF, (int32_t)0x00000000, (int32_t)0x7FFFFFFF, (int32_t)0x00000000, (int32_t)0x7F9AFCB9, (int32_t)0xF5F50D67,
    (int32_t)0x7E6C9251, (int32_t)0xEBF9F498, (int32_t)0x7C769E18, (int32_t)0xE21E765A, (int32_t)0x7E6C9251, (int32_t)0xEBF9F498, (int32_t)0x79BC384D, (int32_t)0xD8722192,
    (int32_t)0x720C8075, (int32_t)0xC5E3A3A9, (int32_t)0x7C769E18, (int32_t)0xE21E765A, (int32_t)0x720C8075, (int32_t)0xC5E3A3A9, (int32_t)0x6154FB91, (int32_t)0xACDEE2E8,
    (int32_t)0x79BC384D, (int32_t)0xD8722192, (int32_t)0x678DDE6E, (int32_t)0xB4C373EE, (int32_t)0x4B3C8C12, (int32_t)0x98722192, (int32_t)0x7641AF3D, (int32_t)0xCF043AB3,
    (int32_t)0x5A82799A, (int32_t)0xA57D8666, (int32_t)0x30FBC54D, (int32_t)0x89BE50C3, (int32_t)0x720C8075, (int32_t)0xC5E3A3A9, (int32_t)0x4B3C8C12, (int32_t)0x98722192,
    (int32_t)0x14060B68, (int32_t)0x81936DAF, (int32_t)0x6D23501B, (int32_t)0xBD1EC45C, (int32_t)0x3A1C5C57, (int32_t)0x8DF37F8B, (int32_t)0xF5F50D67, (int32_t)0x80650347,
    (int32_t)0x678DDE6E, (int32_t)0xB4C373EE, (int32_t)0x278DDE6E, (int32_t)0x8643C7B3, (int32_t)0xD8722192, (int32_t)0x8643C7B3, (int32_t)0x6154FB91, (int32_t)0xACDEE2E8,
    (int32_t)0x14060B68, (int32_t)0x81936DAF, (int32_t)0xBD1EC45C, (int32_t)0x92DCAFE5, (int32_t)0x5A82799A, (int32_t)0xA57D8666, (int32_t)0x00000000, (int32_t)0x80000000,
    (int32_t)0xA57D8666, (int32_t)0xA57D8666, (int32_t)0x53211D18, (int32_t)0x9EAB046F, (int32_t)0xEBF9F498, (int32_t)0x81936DAF, (int32_t)0x92DCAFE5, (int32_t)0xBD1EC45C,
    (int32_t)0x4B3C8C12, (int32_t)0x98722192, (int32_t)0xD8722192, (int32_t)0x8643C7B3, (int32_t)0x8643C7B3, (int32_t)0xD8722192, (int32_t)0x42E13BA4, (int32_t)0x92DCAFE5,
    (int32_t)0xC5E3A3A9, (int32_t)0x8DF37F8B, (int32_t)0x80650347, (int32_t)0xF5F50D67, (int32_t)0x3A1C5C57, (int32_t)0x8DF37F8B, (int32_t)0xB4C373EE, (int32_t)0x98722192,
    (int32_t)0x81936DAF, (int32_t)0x14060B68, (int32_t)0x30FBC54D, (int32_t)0x89BE50C3, (int32_t)0xA57D8666, (int32_t)0xA57D8666, (int32_t)0x89BE50C3, (int32_t)0x30FBC54D,
    (int32_t)0x278DDE6E, (int32_t)0x8643C7B3, (int32_t)0x98722192, (int32_t)0xB4C373EE, (int32_t)0x98722192, (int32_t)0x4B3C8C12, (int32_t)0x1DE189A6, (int32_t)0x838961E8,
    (int32_t)0x8DF37F8B, (int32_t)0xC5E3A3A9, (int32_t)0xACDEE2E8, (int32_t)0x6154FB91, (int32_t)0x14060B68, (int32_t)0x81936DAF, (int32_t)0x8643C7B3, (int32_t)0xD8722192,
    (int32_t)0xC5E3A3A9, (int32_t)0x720C8075, (int32_t)0x0A0AF299, (int32_t)0x80650347, (int32_t)0x81936DAF, (int32_t)0xEBF9F498, (int32_t)0xE21E765A, (int32_t)0x7C769E18,
};

/****************** N=320 stage 3 radix 4 ******************/
ALIGN(32) static const int32_t __fft320_tw3[] =
{
    (int32_t)0x7FFFFFFF, (int32_t)0x00000000, (int32_t)0x7FFFFFFF, (int32_t)0x00000000, (int32_t)0x7FFFFFFF, (int32_t)0x00000000,
    (int32_t)0x79BC384D, (int32_t)0xD8722192, (int32_t)0x678DDE6E, (int32_t)0xB4C373EE, (int32_t)0x4B3C8C12, (int32_t)0x98722192,
    (int32_t)0x678DDE6E, (int32_t)0xB4C373EE, (int32_t)0x278DDE6E, (int32_t)0x8643C7B3, (int32_t)0xD8722192, (int32_t)0x8643C7B3,
    (int32_t)0x4B3C8C12, (int32_t)0x98722192, (int32_t)0xD8722192, (int32_t)0x8643C7B3, (int32_t)0x8643C7B3, (int32_t)0xD8722192,
    (int32_t)0x278DDE6E, (int32_t)0x8643C7B3, (int32_t)0x98722192, (int32_t)0xB4C373EE, (int32_t)0x98722192, (int32_t)0x4B3C8C12,
};
#define N 320
static const fft_cplx32x32_stage_t s2_tab[] = 
{
    fft_stageS2_DFT4_first_32x32,
    fft_stageS2_DFT4x2_32x32,
    fft_stageS2_DFT4x2_32x32,
    fft_stageS2_DFT5_last_32x32,
    NULL
};
static const fft_cplx32x32_stage_t s3_tab[] =
{
    fft_stageS3_DFT4_first_32x32,
    fft_stageS3_DFT4x2_32x32,
    fft_stageS3_DFT4x2_32x32,
    fft_stageS3_DFT5_last_32x32,
    NULL
};
static const fft_cplx32x32_stage_t is2_tab[] = 
{
    ifft_stageS2_DFT4_first_32x32,
    fft_stageS2_DFT4x2_32x32,
    fft_stageS2_DFT4x2_32x32,
    ifft_stageS2_DFT5_last_32x32,
    NULL
};
static const fft_cplx32x32_stage_t is3_tab[] =
{
    ifft_stageS3_DFT4_first_32x32,
    fft_stageS3_DFT4x2_32x32,
    fft_stageS3_DFT4x2_32x32,
    ifft_stageS3_DFT5_last_32x32,
    NULL
};
static const int tw_step_tab[] =
{
    1, 1, 1, 1, 
}; 
static const cint32ptr_t tw_tab[] = 
{
    __fft320_tw1, __fft320_tw2, __fft320_tw3, NULL
};
static const fft_cplx32x32_descr_t __cfft_descr320_32x32 =
{
    N, 
    s2_tab, 
    s3_tab, 
    tw_step_tab,
    tw_tab
};     
static const fft_cplx32x32_descr_t __cifft_descr320_32x32 =
{
    N, 
    is2_tab, 
    is3_tab, 
    tw_step_tab,
    tw_tab
};     
const fft_handle_t cnfft32_320 = (const fft_handle_t)&__cfft_descr320_32x32;
const fft_handle_t cinfft32_320 = (const fft_handle_t)&__cifft_descr320_32x32;
