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

/* Twiddles tables for fft_cplx32x16, ifft_cplx32x16, fft_cplx16x16, ifft_cplx16x16, N=512 */
#define N 512

/****************** stage 1 radix 4 ******************/
ALIGN(32) static const int16_t __fft16_tw1[] =
{
    (int16_t)0x7FFF, (int16_t)0x0000, (int16_t)0x7FFE, (int16_t)0xFE6E, (int16_t)0x7FFF, (int16_t)0x0000, (int16_t)0x7FF6, (int16_t)0xFCDC,
    (int16_t)0x7FFF, (int16_t)0x0000, (int16_t)0x7FEA, (int16_t)0xFB4A, (int16_t)0x7FF6, (int16_t)0xFCDC, (int16_t)0x7FEA, (int16_t)0xFB4A,
    (int16_t)0x7FD9, (int16_t)0xF9B8, (int16_t)0x7FA7, (int16_t)0xF695, (int16_t)0x7FA7, (int16_t)0xF695, (int16_t)0x7F38, (int16_t)0xF1E4,
    (int16_t)0x7FD9, (int16_t)0xF9B8, (int16_t)0x7FC2, (int16_t)0xF827, (int16_t)0x7F62, (int16_t)0xF374, (int16_t)0x7F0A, (int16_t)0xF055,
    (int16_t)0x7E9D, (int16_t)0xED38, (int16_t)0x7DD6, (int16_t)0xE892, (int16_t)0x7FA7, (int16_t)0xF695, (int16_t)0x7F87, (int16_t)0xF505,
    (int16_t)0x7E9D, (int16_t)0xED38, (int16_t)0x7E1E, (int16_t)0xEA1E, (int16_t)0x7CE4, (int16_t)0xE3F4, (int16_t)0x7BC6, (int16_t)0xDF61,
    (int16_t)0x7F62, (int16_t)0xF374, (int16_t)0x7F38, (int16_t)0xF1E4, (int16_t)0x7D8A, (int16_t)0xE707, (int16_t)0x7CE4, (int16_t)0xE3F4,
    (int16_t)0x7A7D, (int16_t)0xDAD8, (int16_t)0x790A, (int16_t)0xD65C, (int16_t)0x7F0A, (int16_t)0xF055, (int16_t)0x7ED6, (int16_t)0xEEC6,
    (int16_t)0x7C2A, (int16_t)0xE0E6, (int16_t)0x7B5D, (int16_t)0xDDDC, (int16_t)0x776C, (int16_t)0xD1EF, (int16_t)0x75A6, (int16_t)0xCD92,
    (int16_t)0x7E9D, (int16_t)0xED38, (int16_t)0x7E60, (int16_t)0xEBAB, (int16_t)0x7A7D, (int16_t)0xDAD8, (int16_t)0x798A, (int16_t)0xD7D9,
    (int16_t)0x73B6, (int16_t)0xC946, (int16_t)0x719E, (int16_t)0xC50D, (int16_t)0x7E1E, (int16_t)0xEA1E, (int16_t)0x7DD6, (int16_t)0xE892,
    (int16_t)0x7885, (int16_t)0xD4E1, (int16_t)0x776C, (int16_t)0xD1EF, (int16_t)0x6F5F, (int16_t)0xC0E9, (int16_t)0x6CF9, (int16_t)0xBCDA,
    (int16_t)0x7D8A, (int16_t)0xE707, (int16_t)0x7D3A, (int16_t)0xE57D, (int16_t)0x7642, (int16_t)0xCF04, (int16_t)0x7505, (int16_t)0xCC21,
    (int16_t)0x6A6E, (int16_t)0xB8E3, (int16_t)0x67BD, (int16_t)0xB505, (int16_t)0x7CE4, (int16_t)0xE3F4, (int16_t)0x7C89, (int16_t)0xE26D,
    (int16_t)0x73B6, (int16_t)0xC946, (int16_t)0x7255, (int16_t)0xC673, (int16_t)0x64E9, (int16_t)0xB140, (int16_t)0x61F1, (int16_t)0xAD97,
    (int16_t)0x7C2A, (int16_t)0xE0E6, (int16_t)0x7BC6, (int16_t)0xDF61, (int16_t)0x70E3, (int16_t)0xC3A9, (int16_t)0x6F5F, (int16_t)0xC0E9,
    (int16_t)0x5ED7, (int16_t)0xAA0A, (int16_t)0x5B9D, (int16_t)0xA69C, (int16_t)0x7B5D, (int16_t)0xDDDC, (int16_t)0x7AEF, (int16_t)0xDC59,
    (int16_t)0x6DCA, (int16_t)0xBE32, (int16_t)0x6C24, (int16_t)0xBB85, (int16_t)0x5843, (int16_t)0xA34C, (int16_t)0x54CA, (int16_t)0xA01C,
    (int16_t)0x7A7D, (int16_t)0xDAD8, (int16_t)0x7A06, (int16_t)0xD958, (int16_t)0x6A6E, (int16_t)0xB8E3, (int16_t)0x68A7, (int16_t)0xB64C,
    (int16_t)0x5134, (int16_t)0x9D0E, (int16_t)0x4D81, (int16_t)0x9A22, (int16_t)0x798A, (int16_t)0xD7D9, (int16_t)0x790A, (int16_t)0xD65C,
    (int16_t)0x66D0, (int16_t)0xB3C0, (int16_t)0x64E9, (int16_t)0xB140, (int16_t)0x49B4, (int16_t)0x9759, (int16_t)0x45CD, (int16_t)0x94B5,
    (int16_t)0x7885, (int16_t)0xD4E1, (int16_t)0x77FB, (int16_t)0xD367, (int16_t)0x62F2, (int16_t)0xAECC, (int16_t)0x60EC, (int16_t)0xAC65,
    (int16_t)0x41CE, (int16_t)0x9236, (int16_t)0x3DB8, (int16_t)0x8FDD, (int16_t)0x776C, (int16_t)0xD1EF, (int16_t)0x76D9, (int16_t)0xD079,
    (int16_t)0x5ED7, (int16_t)0xAA0A, (int16_t)0x5CB4, (int16_t)0xA7BD, (int16_t)0x398D, (int16_t)0x8DAB, (int16_t)0x354E, (int16_t)0x8BA0,
    (int16_t)0x7642, (int16_t)0xCF04, (int16_t)0x75A6, (int16_t)0xCD92, (int16_t)0x5A82, (int16_t)0xA57E, (int16_t)0x5843, (int16_t)0xA34C,
    (int16_t)0x30FC, (int16_t)0x89BE, (int16_t)0x2C99, (int16_t)0x8805, (int16_t)0x7505, (int16_t)0xCC21, (int16_t)0x7460, (int16_t)0xCAB2,
    (int16_t)0x55F6, (int16_t)0xA129, (int16_t)0x539B, (int16_t)0x9F14, (int16_t)0x2827, (int16_t)0x8676, (int16_t)0x23A7, (int16_t)0x8511,
    (int16_t)0x73B6, (int16_t)0xC946, (int16_t)0x7308, (int16_t)0xC7DB, (int16_t)0x5134, (int16_t)0x9D0E, (int16_t)0x4EC0, (int16_t)0x9B17,
    (int16_t)0x1F1A, (int16_t)0x83D6, (int16_t)0x1A83, (int16_t)0x82C6, (int16_t)0x7255, (int16_t)0xC673, (int16_t)0x719E, (int16_t)0xC50D,
    (int16_t)0x4C40, (int16_t)0x9930, (int16_t)0x49B4, (int16_t)0x9759, (int16_t)0x15E2, (int16_t)0x81E2, (int16_t)0x113A, (int16_t)0x812A,
    (int16_t)0x70E3, (int16_t)0xC3A9, (int16_t)0x7023, (int16_t)0xC248, (int16_t)0x471D, (int16_t)0x9592, (int16_t)0x447B, (int16_t)0x93DC,
    (int16_t)0x0C8C, (int16_t)0x809E, (int16_t)0x07D9, (int16_t)0x803E, (int16_t)0x6F5F, (int16_t)0xC0E9, (int16_t)0x6E97, (int16_t)0xBF8C,
    (int16_t)0x41CE, (int16_t)0x9236, (int16_t)0x3F17, (int16_t)0x90A1, (int16_t)0x0324, (int16_t)0x800A, (int16_t)0xFE6E, (int16_t)0x8002,
    (int16_t)0x6DCA, (int16_t)0xBE32, (int16_t)0x6CF9, (int16_t)0xBCDA, (int16_t)0x3C57, (int16_t)0x8F1D, (int16_t)0x398D, (int16_t)0x8DAB,
    (int16_t)0xF9B8, (int16_t)0x8027, (int16_t)0xF505, (int16_t)0x8079, (int16_t)0x6C24, (int16_t)0xBB85, (int16_t)0x6B4B, (int16_t)0xBA33,
    (int16_t)0x36BA, (int16_t)0x8C4A, (int16_t)0x33DF, (int16_t)0x8AFB, (int16_t)0xF055, (int16_t)0x80F6, (int16_t)0xEBAB, (int16_t)0x81A0,
    (int16_t)0x6A6E, (int16_t)0xB8E3, (int16_t)0x698C, (int16_t)0xB796, (int16_t)0x30FC, (int16_t)0x89BE, (int16_t)0x2E11, (int16_t)0x8894,
    (int16_t)0xE707, (int16_t)0x8276, (int16_t)0xE26D, (int16_t)0x8377, (int16_t)0x68A7, (int16_t)0xB64C, (int16_t)0x67BD, (int16_t)0xB505,
    (int16_t)0x2B1F, (int16_t)0x877B, (int16_t)0x2827, (int16_t)0x8676, (int16_t)0xDDDC, (int16_t)0x84A3, (int16_t)0xD958, (int16_t)0x85FA,
    (int16_t)0x66D0, (int16_t)0xB3C0, (int16_t)0x65DE, (int16_t)0xB27F, (int16_t)0x2528, (int16_t)0x8583, (int16_t)0x2224, (int16_t)0x84A3,
    (int16_t)0xD4E1, (int16_t)0x877B, (int16_t)0xD079, (int16_t)0x8927, (int16_t)0x64E9, (int16_t)0xB140, (int16_t)0x63EF, (int16_t)0xB005,
    (int16_t)0x1F1A, (int16_t)0x83D6, (int16_t)0x1C0C, (int16_t)0x831C, (int16_t)0xCC21, (int16_t)0x8AFB, (int16_t)0xC7DB, (int16_t)0x8CF8,
    (int16_t)0x62F2, (int16_t)0xAECC, (int16_t)0x61F1, (int16_t)0xAD97, (int16_t)0x18F9, (int16_t)0x8276, (int16_t)0x15E2, (int16_t)0x81E2,
    (int16_t)0xC3A9, (int16_t)0x8F1D, (int16_t)0xBF8C, (int16_t)0x9169, (int16_t)0x60EC, (int16_t)0xAC65, (int16_t)0x5FE4, (int16_t)0xAB36,
    (int16_t)0x12C8, (int16_t)0x8163, (int16_t)0x0FAB, (int16_t)0x80F6, (int16_t)0xBB85, (int16_t)0x93DC, (int16_t)0xB796, (int16_t)0x9674,
    (int16_t)0x5ED7, (int16_t)0xAA0A, (int16_t)0x5DC8, (int16_t)0xA8E2, (int16_t)0x0C8C, (int16_t)0x809E, (int16_t)0x096B, (int16_t)0x8059,
    (int16_t)0xB3C0, (int16_t)0x9930, (int16_t)0xB005, (int16_t)0x9C11, (int16_t)0x5CB4, (int16_t)0xA7BD, (int16_t)0x5B9D, (int16_t)0xA69C,
    (int16_t)0x0648, (int16_t)0x8027, (int16_t)0x0324, (int16_t)0x800A, (int16_t)0xAC65, (int16_t)0x9F14, (int16_t)0xA8E2, (int16_t)0xA238,
    (int16_t)0x5A82, (int16_t)0xA57E, (int16_t)0x5964, (int16_t)0xA463, (int16_t)0x0000, (int16_t)0x8000, (int16_t)0xFCDC, (int16_t)0x800A,
    (int16_t)0xA57E, (int16_t)0xA57E, (int16_t)0xA238, (int16_t)0xA8E2, (int16_t)0x5843, (int16_t)0xA34C, (int16_t)0x571E, (int16_t)0xA238,
    (int16_t)0xF9B8, (int16_t)0x8027, (int16_t)0xF695, (int16_t)0x8059, (int16_t)0x9F14, (int16_t)0xAC65, (int16_t)0x9C11, (int16_t)0xB005,
    (int16_t)0x55F6, (int16_t)0xA129, (int16_t)0x54CA, (int16_t)0xA01C, (int16_t)0xF374, (int16_t)0x809E, (int16_t)0xF055, (int16_t)0x80F6,
    (int16_t)0x9930, (int16_t)0xB3C0, (int16_t)0x9674, (int16_t)0xB796, (int16_t)0x539B, (int16_t)0x9F14, (int16_t)0x5269, (int16_t)0x9E0F,
    (int16_t)0xED38, (int16_t)0x8163, (int16_t)0xEA1E, (int16_t)0x81E2, (int16_t)0x93DC, (int16_t)0xBB85, (int16_t)0x9169, (int16_t)0xBF8C,
    (int16_t)0x5134, (int16_t)0x9D0E, (int16_t)0x4FFB, (int16_t)0x9C11, (int16_t)0xE707, (int16_t)0x8276, (int16_t)0xE3F4, (int16_t)0x831C,
    (int16_t)0x8F1D, (int16_t)0xC3A9, (int16_t)0x8CF8, (int16_t)0xC7DB, (int16_t)0x4EC0, (int16_t)0x9B17, (int16_t)0x4D81, (int16_t)0x9A22,
    (int16_t)0xE0E6, (int16_t)0x83D6, (int16_t)0xDDDC, (int16_t)0x84A3, (int16_t)0x8AFB, (int16_t)0xCC21, (int16_t)0x8927, (int16_t)0xD079,
    (int16_t)0x4C40, (int16_t)0x9930, (int16_t)0x4AFB, (int16_t)0x9843, (int16_t)0xDAD8, (int16_t)0x8583, (int16_t)0xD7D9, (int16_t)0x8676,
    (int16_t)0x877B, (int16_t)0xD4E1, (int16_t)0x85FA, (int16_t)0xD958, (int16_t)0x49B4, (int16_t)0x9759, (int16_t)0x486A, (int16_t)0x9674,
    (int16_t)0xD4E1, (int16_t)0x877B, (int16_t)0xD1EF, (int16_t)0x8894, (int16_t)0x84A3, (int16_t)0xDDDC, (int16_t)0x8377, (int16_t)0xE26D,
    (int16_t)0x471D, (int16_t)0x9592, (int16_t)0x45CD, (int16_t)0x94B5, (int16_t)0xCF04, (int16_t)0x89BE, (int16_t)0xCC21, (int16_t)0x8AFB,
    (int16_t)0x8276, (int16_t)0xE707, (int16_t)0x81A0, (int16_t)0xEBAB, (int16_t)0x447B, (int16_t)0x93DC, (int16_t)0x4326, (int16_t)0x9307,
    (int16_t)0xC946, (int16_t)0x8C4A, (int16_t)0xC673, (int16_t)0x8DAB, (int16_t)0x80F6, (int16_t)0xF055, (int16_t)0x8079, (int16_t)0xF505,
    (int16_t)0x41CE, (int16_t)0x9236, (int16_t)0x4074, (int16_t)0x9169, (int16_t)0xC3A9, (int16_t)0x8F1D, (int16_t)0xC0E9, (int16_t)0x90A1,
    (int16_t)0x8027, (int16_t)0xF9B8, (int16_t)0x8002, (int16_t)0xFE6E, (int16_t)0x3F17, (int16_t)0x90A1, (int16_t)0x3DB8, (int16_t)0x8FDD,
    (int16_t)0xBE32, (int16_t)0x9236, (int16_t)0xBB85, (int16_t)0x93DC, (int16_t)0x800A, (int16_t)0x0324, (int16_t)0x803E, (int16_t)0x07D9,
    (int16_t)0x3C57, (int16_t)0x8F1D, (int16_t)0x3AF3, (int16_t)0x8E62, (int16_t)0xB8E3, (int16_t)0x9592, (int16_t)0xB64C, (int16_t)0x9759,
    (int16_t)0x809E, (int16_t)0x0C8C, (int16_t)0x812A, (int16_t)0x113A, (int16_t)0x398D, (int16_t)0x8DAB, (int16_t)0x3825, (int16_t)0x8CF8,
    (int16_t)0xB3C0, (int16_t)0x9930, (int16_t)0xB140, (int16_t)0x9B17, (int16_t)0x81E2, (int16_t)0x15E2, (int16_t)0x82C6, (int16_t)0x1A83,
    (int16_t)0x36BA, (int16_t)0x8C4A, (int16_t)0x354E, (int16_t)0x8BA0, (int16_t)0xAECC, (int16_t)0x9D0E, (int16_t)0xAC65, (int16_t)0x9F14,
    (int16_t)0x83D6, (int16_t)0x1F1A, (int16_t)0x8511, (int16_t)0x23A7, (int16_t)0x33DF, (int16_t)0x8AFB, (int16_t)0x326E, (int16_t)0x8A5A,
    (int16_t)0xAA0A, (int16_t)0xA129, (int16_t)0xA7BD, (int16_t)0xA34C, (int16_t)0x8676, (int16_t)0x2827, (int16_t)0x8805, (int16_t)0x2C99,
    (int16_t)0x30FC, (int16_t)0x89BE, (int16_t)0x2F87, (int16_t)0x8927, (int16_t)0xA57E, (int16_t)0xA57E, (int16_t)0xA34C, (int16_t)0xA7BD,
    (int16_t)0x89BE, (int16_t)0x30FC, (int16_t)0x8BA0, (int16_t)0x354E, (int16_t)0x2E11, (int16_t)0x8894, (int16_t)0x2C99, (int16_t)0x8805,
    (int16_t)0xA129, (int16_t)0xAA0A, (int16_t)0x9F14, (int16_t)0xAC65, (int16_t)0x8DAB, (int16_t)0x398D, (int16_t)0x8FDD, (int16_t)0x3DB8,
    (int16_t)0x2B1F, (int16_t)0x877B, (int16_t)0x29A4, (int16_t)0x86F6, (int16_t)0x9D0E, (int16_t)0xAECC, (int16_t)0x9B17, (int16_t)0xB140,
    (int16_t)0x9236, (int16_t)0x41CE, (int16_t)0x94B5, (int16_t)0x45CD, (int16_t)0x2827, (int16_t)0x8676, (int16_t)0x26A8, (int16_t)0x85FA,
    (int16_t)0x9930, (int16_t)0xB3C0, (int16_t)0x9759, (int16_t)0xB64C, (int16_t)0x9759, (int16_t)0x49B4, (int16_t)0x9A22, (int16_t)0x4D81,
    (int16_t)0x2528, (int16_t)0x8583, (int16_t)0x23A7, (int16_t)0x8511, (int16_t)0x9592, (int16_t)0xB8E3, (int16_t)0x93DC, (int16_t)0xBB85,
    (int16_t)0x9D0E, (int16_t)0x5134, (int16_t)0xA01C, (int16_t)0x54CA, (int16_t)0x2224, (int16_t)0x84A3, (int16_t)0x209F, (int16_t)0x843A,
    (int16_t)0x9236, (int16_t)0xBE32, (int16_t)0x90A1, (int16_t)0xC0E9, (int16_t)0xA34C, (int16_t)0x5843, (int16_t)0xA69C, (int16_t)0x5B9D,
    (int16_t)0x1F1A, (int16_t)0x83D6, (int16_t)0x1D93, (int16_t)0x8377, (int16_t)0x8F1D, (int16_t)0xC3A9, (int16_t)0x8DAB, (int16_t)0xC673,
    (int16_t)0xAA0A, (int16_t)0x5ED7, (int16_t)0xAD97, (int16_t)0x61F1, (int16_t)0x1C0C, (int16_t)0x831C, (int16_t)0x1A83, (int16_t)0x82C6,
    (int16_t)0x8C4A, (int16_t)0xC946, (int16_t)0x8AFB, (int16_t)0xCC21, (int16_t)0xB140, (int16_t)0x64E9, (int16_t)0xB505, (int16_t)0x67BD,
    (int16_t)0x18F9, (int16_t)0x8276, (int16_t)0x176E, (int16_t)0x822A, (int16_t)0x89BE, (int16_t)0xCF04, (int16_t)0x8894, (int16_t)0xD1EF,
    (int16_t)0xB8E3, (int16_t)0x6A6E, (int16_t)0xBCDA, (int16_t)0x6CF9, (int16_t)0x15E2, (int16_t)0x81E2, (int16_t)0x1455, (int16_t)0x81A0,
    (int16_t)0x877B, (int16_t)0xD4E1, (int16_t)0x8676, (int16_t)0xD7D9, (int16_t)0xC0E9, (int16_t)0x6F5F, (int16_t)0xC50D, (int16_t)0x719E,
    (int16_t)0x12C8, (int16_t)0x8163, (int16_t)0x113A, (int16_t)0x812A, (int16_t)0x8583, (int16_t)0xDAD8, (int16_t)0x84A3, (int16_t)0xDDDC,
    (int16_t)0xC946, (int16_t)0x73B6, (int16_t)0xCD92, (int16_t)0x75A6, (int16_t)0x0FAB, (int16_t)0x80F6, (int16_t)0x0E1C, (int16_t)0x80C8,
    (int16_t)0x83D6, (int16_t)0xE0E6, (int16_t)0x831C, (int16_t)0xE3F4, (int16_t)0xD1EF, (int16_t)0x776C, (int16_t)0xD65C, (int16_t)0x790A,
    (int16_t)0x0C8C, (int16_t)0x809E, (int16_t)0x0AFB, (int16_t)0x8079, (int16_t)0x8276, (int16_t)0xE707, (int16_t)0x81E2, (int16_t)0xEA1E,
    (int16_t)0xDAD8, (int16_t)0x7A7D, (int16_t)0xDF61, (int16_t)0x7BC6, (int16_t)0x096B, (int16_t)0x8059, (int16_t)0x07D9, (int16_t)0x803E,
    (int16_t)0x8163, (int16_t)0xED38, (int16_t)0x80F6, (int16_t)0xF055, (int16_t)0xE3F4, (int16_t)0x7CE4, (int16_t)0xE892, (int16_t)0x7DD6,
    (int16_t)0x0648, (int16_t)0x8027, (int16_t)0x04B6, (int16_t)0x8016, (int16_t)0x809E, (int16_t)0xF374, (int16_t)0x8059, (int16_t)0xF695,
    (int16_t)0xED38, (int16_t)0x7E9D, (int16_t)0xF1E4, (int16_t)0x7F38, (int16_t)0x0324, (int16_t)0x800A, (int16_t)0x0192, (int16_t)0x8002,
    (int16_t)0x8027, (int16_t)0xF9B8, (int16_t)0x800A, (int16_t)0xFCDC, (int16_t)0xF695, (int16_t)0x7FA7, (int16_t)0xFB4A, (int16_t)0x7FEA,
};
/****************** stage 2 radix 4 ******************/
ALIGN(32) static const int16_t __fft16_tw2[] =
{
(int16_t)0x7fff,(int16_t)0x0000,(int16_t)0x7fff,(int16_t)0x0000,(int16_t)0x7fff,(int16_t)0x0000,
(int16_t)0x7fd9,(int16_t)0xf9b8,(int16_t)0x7f62,(int16_t)0xf374,(int16_t)0x7e9d,(int16_t)0xed38,
(int16_t)0x7f62,(int16_t)0xf374,(int16_t)0x7d8a,(int16_t)0xe707,(int16_t)0x7a7d,(int16_t)0xdad8,
(int16_t)0x7e9d,(int16_t)0xed38,(int16_t)0x7a7d,(int16_t)0xdad8,(int16_t)0x73b6,(int16_t)0xc946,
(int16_t)0x7d8a,(int16_t)0xe707,(int16_t)0x7642,(int16_t)0xcf04,(int16_t)0x6a6e,(int16_t)0xb8e3,
(int16_t)0x7c2a,(int16_t)0xe0e6,(int16_t)0x70e3,(int16_t)0xc3a9,(int16_t)0x5ed7,(int16_t)0xaa0a,
(int16_t)0x7a7d,(int16_t)0xdad8,(int16_t)0x6a6e,(int16_t)0xb8e3,(int16_t)0x5134,(int16_t)0x9d0e,
(int16_t)0x7885,(int16_t)0xd4e1,(int16_t)0x62f2,(int16_t)0xaecc,(int16_t)0x41ce,(int16_t)0x9236,
(int16_t)0x7642,(int16_t)0xcf04,(int16_t)0x5a82,(int16_t)0xa57e,(int16_t)0x30fc,(int16_t)0x89be,
(int16_t)0x73b6,(int16_t)0xc946,(int16_t)0x5134,(int16_t)0x9d0e,(int16_t)0x1f1a,(int16_t)0x83d6,
(int16_t)0x70e3,(int16_t)0xc3a9,(int16_t)0x471d,(int16_t)0x9592,(int16_t)0x0c8c,(int16_t)0x809e,
(int16_t)0x6dca,(int16_t)0xbe32,(int16_t)0x3c57,(int16_t)0x8f1d,(int16_t)0xf9b8,(int16_t)0x8027,
(int16_t)0x6a6e,(int16_t)0xb8e3,(int16_t)0x30fc,(int16_t)0x89be,(int16_t)0xe707,(int16_t)0x8276,
(int16_t)0x66d0,(int16_t)0xb3c0,(int16_t)0x2528,(int16_t)0x8583,(int16_t)0xd4e1,(int16_t)0x877b,
(int16_t)0x62f2,(int16_t)0xaecc,(int16_t)0x18f9,(int16_t)0x8276,(int16_t)0xc3a9,(int16_t)0x8f1d,
(int16_t)0x5ed7,(int16_t)0xaa0a,(int16_t)0x0c8c,(int16_t)0x809e,(int16_t)0xb3c0,(int16_t)0x9930,
(int16_t)0x5a82,(int16_t)0xa57e,(int16_t)0x0000,(int16_t)0x8000,(int16_t)0xa57e,(int16_t)0xa57e,
(int16_t)0x55f6,(int16_t)0xa129,(int16_t)0xf374,(int16_t)0x809e,(int16_t)0x9930,(int16_t)0xb3c0,
(int16_t)0x5134,(int16_t)0x9d0e,(int16_t)0xe707,(int16_t)0x8276,(int16_t)0x8f1d,(int16_t)0xc3a9,
(int16_t)0x4c40,(int16_t)0x9930,(int16_t)0xdad8,(int16_t)0x8583,(int16_t)0x877b,(int16_t)0xd4e1,
(int16_t)0x471d,(int16_t)0x9592,(int16_t)0xcf04,(int16_t)0x89be,(int16_t)0x8276,(int16_t)0xe707,
(int16_t)0x41ce,(int16_t)0x9236,(int16_t)0xc3a9,(int16_t)0x8f1d,(int16_t)0x8027,(int16_t)0xf9b8,
(int16_t)0x3c57,(int16_t)0x8f1d,(int16_t)0xb8e3,(int16_t)0x9592,(int16_t)0x809e,(int16_t)0x0c8c,
(int16_t)0x36ba,(int16_t)0x8c4a,(int16_t)0xaecc,(int16_t)0x9d0e,(int16_t)0x83d6,(int16_t)0x1f1a,
(int16_t)0x30fc,(int16_t)0x89be,(int16_t)0xa57e,(int16_t)0xa57e,(int16_t)0x89be,(int16_t)0x30fc,
(int16_t)0x2b1f,(int16_t)0x877b,(int16_t)0x9d0e,(int16_t)0xaecc,(int16_t)0x9236,(int16_t)0x41ce,
(int16_t)0x2528,(int16_t)0x8583,(int16_t)0x9592,(int16_t)0xb8e3,(int16_t)0x9d0e,(int16_t)0x5134,
(int16_t)0x1f1a,(int16_t)0x83d6,(int16_t)0x8f1d,(int16_t)0xc3a9,(int16_t)0xaa0a,(int16_t)0x5ed7,
(int16_t)0x18f9,(int16_t)0x8276,(int16_t)0x89be,(int16_t)0xcf04,(int16_t)0xb8e3,(int16_t)0x6a6e,
(int16_t)0x12c8,(int16_t)0x8163,(int16_t)0x8583,(int16_t)0xdad8,(int16_t)0xc946,(int16_t)0x73b6,
(int16_t)0x0c8c,(int16_t)0x809e,(int16_t)0x8276,(int16_t)0xe707,(int16_t)0xdad8,(int16_t)0x7a7d,
(int16_t)0x0648,(int16_t)0x8027,(int16_t)0x809e,(int16_t)0xf374,(int16_t)0xed38,(int16_t)0x7e9d
};
/****************** stage 3 radix 4 ******************/
ALIGN(32) static const int16_t __fft16_tw3[] =
{
(int16_t)0x7fff,(int16_t)0x0000,(int16_t)0x7fff,(int16_t)0x0000,(int16_t)0x7fff,(int16_t)0x0000,
(int16_t)0x7d8a,(int16_t)0xe707,(int16_t)0x7642,(int16_t)0xcf04,(int16_t)0x6a6e,(int16_t)0xb8e3,
(int16_t)0x7642,(int16_t)0xcf04,(int16_t)0x5a82,(int16_t)0xa57e,(int16_t)0x30fc,(int16_t)0x89be,
(int16_t)0x6a6e,(int16_t)0xb8e3,(int16_t)0x30fc,(int16_t)0x89be,(int16_t)0xe707,(int16_t)0x8276,
(int16_t)0x5a82,(int16_t)0xa57e,(int16_t)0x0000,(int16_t)0x8000,(int16_t)0xa57e,(int16_t)0xa57e,
(int16_t)0x471d,(int16_t)0x9592,(int16_t)0xcf04,(int16_t)0x89be,(int16_t)0x8276,(int16_t)0xe707,
(int16_t)0x30fc,(int16_t)0x89be,(int16_t)0xa57e,(int16_t)0xa57e,(int16_t)0x89be,(int16_t)0x30fc,
(int16_t)0x18f9,(int16_t)0x8276,(int16_t)0x89be,(int16_t)0xcf04,(int16_t)0xb8e3,(int16_t)0x6a6e
};
static const int tw_step_tab[] =
{
    1, 1, 1, 0
}; 
static const cint16ptr_t_fft tw_tab[] = 
{
    __fft16_tw1, __fft16_tw2, __fft16_tw3,  NULL
};

static const fn_fft_stage fft32x16_stg_tab_s2[] = 
{
    (fn_fft_stage)fft_32x16_stage_first_scl2_DFT4,
    (fn_fft_stage)fft_32x16_stage_inner_scl2_DFT4x2_v4,
    (fn_fft_stage)fft_32x16_stage_inner_scl2_DFT4x2,
    (fn_fft_stage)fft_32x16_stage_last_scl2_DFT8
};
static const fn_fft_stage fft32x16_stg_tab_s3[] = 
{
    (fn_fft_stage)fft_32x16_stage_first_scl3_DFT4,
    (fn_fft_stage)fft_32x16_stage_inner_scl3_DFT4x2_v4,
    (fn_fft_stage)fft_32x16_stage_inner_scl3_DFT4x2,
    (fn_fft_stage)fft_32x16_stage_last_scl3_DFT8
};
static const fn_fft_stage ifft32x16_stg_tab_s2[] =
{
    (fn_fft_stage)ifft_32x16_stage_first_scl2_DFT4,
    (fn_fft_stage)fft_32x16_stage_inner_scl2_DFT4x2_v4,
    (fn_fft_stage)fft_32x16_stage_inner_scl2_DFT4x2,
    (fn_fft_stage)ifft_32x16_stage_last_scl2_DFT8
};
static const fn_fft_stage ifft32x16_stg_tab_s3[] =
{
    (fn_fft_stage)ifft_32x16_stage_first_scl3_DFT4,
    (fn_fft_stage)fft_32x16_stage_inner_scl3_DFT4x2_v4,
    (fn_fft_stage)fft_32x16_stage_inner_scl3_DFT4x2,
    (fn_fft_stage)ifft_32x16_stage_last_scl3_DFT8
};

static const fn_fft_stage fft_stg_tab_s2[] = 
{
    (fn_fft_stage)fft_16x16_stage_first_scl2_DFT4,
    (fn_fft_stage)fft_16x16_stage_inner_scl2_DFT4x2_v4,
    (fn_fft_stage)fft_16x16_stage_inner_scl2_DFT4x2,
    (fn_fft_stage)fft_16x16_stage_last_scl2_DFT8
};
static const fn_fft_stage fft_stg_tab_s3[] = 
{
    (fn_fft_stage)fft_16x16_stage_first_scl3_DFT4,
    (fn_fft_stage)fft_16x16_stage_inner_scl3_DFT4x2_v4,
    (fn_fft_stage)fft_16x16_stage_inner_scl3_DFT4x2,
    (fn_fft_stage)fft_16x16_stage_last_scl3_DFT8
};
static const fn_fft_stage ifft_stg_tab_s2[] =
{
    (fn_fft_stage)ifft_16x16_stage_first_scl2_DFT4,
    (fn_fft_stage)fft_16x16_stage_inner_scl2_DFT4x2_v4,
    (fn_fft_stage)fft_16x16_stage_inner_scl2_DFT4x2,
    (fn_fft_stage)ifft_16x16_stage_last_scl2_DFT8
};
static const fn_fft_stage ifft_stg_tab_s3[] =
{
    (fn_fft_stage)ifft_16x16_stage_first_scl3_DFT4,
    (fn_fft_stage)fft_16x16_stage_inner_scl3_DFT4x2_v4,
    (fn_fft_stage)fft_16x16_stage_inner_scl3_DFT4x2,
    (fn_fft_stage)ifft_16x16_stage_last_scl3_DFT8
};

const fft_cplx_x16_descr_t __cfft_x16_descr512 = 
{
    N, 
    tw_step_tab,
    tw_tab,
    fft32x16_stg_tab_s2,
    fft_stg_tab_s2,
    fft32x16_stg_tab_s3, 
    fft_stg_tab_s3
};     
const fft_cplx_x16_descr_t __cifft_x16_descr512 =
{
    N, 
    tw_step_tab,
    tw_tab,
    ifft32x16_stg_tab_s2,
    ifft_stg_tab_s2,
    ifft32x16_stg_tab_s3, 
    ifft_stg_tab_s3
};     
const fft_handle_t cfft16_512 = (const fft_handle_t)&__cfft_x16_descr512;
const fft_handle_t cifft16_512 = (const fft_handle_t)&__cifft_x16_descr512;
