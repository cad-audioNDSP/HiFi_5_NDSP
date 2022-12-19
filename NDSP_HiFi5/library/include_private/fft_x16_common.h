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
/*
	NatureDSP_Signal library. FFT part
    Twiddle factor tables and descriptor structures for 32x16 & 16x16 FFTs
	IntegrIT, 2006-2018
*/

#ifndef __FFT_X16_COMMON_H__
#define __FFT_X16_COMMON_H__
#include "NatureDSP_types.h"
#include "NatureDSP_Signal_fft.h"


typedef const int16_t* cint16ptr_t_fft;

/*
 * 32x16-bit & 16x16-bit complex-valued FFT descriptor structure.
 */
typedef int (*fn_fft_stage)(const int16_t *tw, const void *x, void *y, int N, int *v, int tw_step, int *bexp);
typedef struct
{
    const int N;
    /* Twiddle factors (access step and tables) */
    const int *tw_step; 
    const cint16ptr_t_fft *twd;
    const fn_fft_stage* fnstages_32x16_s2;
    const fn_fft_stage* fnstages_16x16_s2;
    const fn_fft_stage* fnstages_32x16_s3;
    const fn_fft_stage* fnstages_16x16_s3;
} fft_cplx_x16_descr_t;


/*
 *   32x16-bit or 16x16-bit real-valued FFT descriptor structure.
 */
typedef struct
{
    /* Handle of half-size complex FFT */
    const fft_cplx_x16_descr_t *cfft_hdl; 
    const int16_t *twd;
} fft_real_x16_descr_t;

/* Descriptors that are used in real-valued FFT as half-sized complex FFT */
extern const fft_cplx_x16_descr_t __cfft_x16_descr16;
extern const fft_cplx_x16_descr_t __cfft_x16_descr32;
extern const fft_cplx_x16_descr_t __cfft_x16_descr64;
extern const fft_cplx_x16_descr_t __cfft_x16_descr128;
extern const fft_cplx_x16_descr_t __cfft_x16_descr256;
extern const fft_cplx_x16_descr_t __cfft_x16_descr512;
extern const fft_cplx_x16_descr_t __cfft_x16_descr1024;
extern const fft_cplx_x16_descr_t __cfft_x16_descr2048;
extern const fft_cplx_x16_descr_t __cfft_x16_descr4096;

extern const fft_cplx_x16_descr_t __cifft_x16_descr16;
extern const fft_cplx_x16_descr_t __cifft_x16_descr32;
extern const fft_cplx_x16_descr_t __cifft_x16_descr64;
extern const fft_cplx_x16_descr_t __cifft_x16_descr128;
extern const fft_cplx_x16_descr_t __cifft_x16_descr256;
extern const fft_cplx_x16_descr_t __cifft_x16_descr512;
extern const fft_cplx_x16_descr_t __cifft_x16_descr1024;
extern const fft_cplx_x16_descr_t __cifft_x16_descr2048;
extern const fft_cplx_x16_descr_t __cifft_x16_descr4096;

extern const fft_cplx_x16_descr_t __cfft_x16_descr160;
extern const fft_cplx_x16_descr_t __cfft_x16_descr192;
extern const fft_cplx_x16_descr_t __cfft_x16_descr240;
extern const fft_cplx_x16_descr_t __cfft_x16_descr160;
extern const fft_cplx_x16_descr_t __cfft_x16_descr192;
extern const fft_cplx_x16_descr_t __cfft_x16_descr240;
extern const fft_cplx_x16_descr_t __cfft_x16_descr320;
extern const fft_cplx_x16_descr_t __cifft_x16_descr144;
extern const fft_cplx_x16_descr_t __cfft_x16_descr144;


extern const fft_cplx_x16_descr_t __cfft_32x16_descr160;
extern const fft_cplx_x16_descr_t __cfft_32x16_descr192;
extern const fft_cplx_x16_descr_t __cfft_32x16_descr240;
extern const fft_cplx_x16_descr_t __cfft_32x16_descr160;
extern const fft_cplx_x16_descr_t __cfft_32x16_descr192;
extern const fft_cplx_x16_descr_t __cfft_32x16_descr240;

extern const fft_cplx_x16_descr_t __cifft_x16_descr160;
extern const fft_cplx_x16_descr_t __cifft_x16_descr192;
extern const fft_cplx_x16_descr_t __cifft_x16_descr240;
extern const fft_cplx_x16_descr_t __cifft_x16_descr160;
extern const fft_cplx_x16_descr_t __cifft_x16_descr192;
extern const fft_cplx_x16_descr_t __cifft_x16_descr240;
extern const fft_cplx_x16_descr_t __cifft_x16_descr320;
extern const fft_cplx_x16_descr_t __cfft_x16_descr176;
extern const fft_cplx_x16_descr_t __cifft_x16_descr176;

extern const fft_cplx_x16_descr_t __cifft_32x16_descr160;
extern const fft_cplx_x16_descr_t __cifft_32x16_descr192;
extern const fft_cplx_x16_descr_t __cifft_32x16_descr240;
extern const fft_cplx_x16_descr_t __cifft_32x16_descr160;
extern const fft_cplx_x16_descr_t __cifft_32x16_descr192;
extern const fft_cplx_x16_descr_t __cifft_32x16_descr240;

extern const fft_cplx_x16_descr_t __cfft_x16_descr288;
extern const fft_cplx_x16_descr_t __cfft_x16_descr48;
extern const fft_cplx_x16_descr_t __cifft_x16_descr48;
extern const fft_cplx_x16_descr_t __cifft_x16_descr288;

extern const  fft_handle_t cnfft16_640; 
extern const  fft_handle_t cnfft16_576;
extern const  fft_handle_t cnfft16_176;
extern const  fft_handle_t cnfft16_352;
extern const  fft_handle_t cnfft16_144;
extern const  fft_handle_t cnfft16_288;
extern const  fft_handle_t cinfft16_144;
extern const  fft_handle_t cnfft16_96;
extern const  fft_handle_t cinfft16_96;
extern const  fft_handle_t rnfft16_576;
extern const  fft_handle_t rinfft16_576;
extern const  fft_handle_t rnfft16_144;
extern const  fft_handle_t rnfft16_352; 
extern const  fft_handle_t rinfft16_144;
extern const  fft_handle_t rnfft16_96;
extern const  fft_handle_t rinfft16_96;
extern const  fft_handle_t rnfft16_176;
extern const  fft_handle_t rinfft16_176;
extern const  fft_handle_t cinfft16_288;
extern const  fft_handle_t cinfft16_176;
extern const  fft_handle_t cinfft16_352;
extern const  fft_handle_t cinfft16_576;
extern const  fft_handle_t cinfft16_640;
extern const  fft_handle_t rinfft16_288;
extern const  fft_handle_t rinfft16_352;
extern const  fft_handle_t rnfft16_288;
extern const  fft_handle_t rnfft16_640;
extern const  fft_handle_t rinfft16_640;


#endif // __FFT_X16_COMMON_H__
