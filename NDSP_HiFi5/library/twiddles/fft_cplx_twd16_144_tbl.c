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

/* Twiddles tables for fft_cplx32x16, ifft_cplx32x16, fft_cplx16x16, ifft_cplx16x16, N=144 */
#define N 144

/********** Twiddles table N=144 stage 1 radix 4 ******************/
ALIGN(32) static const int16_t _fft144_tw1_[] =
#if 0 // old HIFI4 packing
{
    (int16_t)0x7FFF, (int16_t)0x0000, (int16_t)0x7FFF, (int16_t)0x0000, (int16_t)0x7FFF, (int16_t)0x0000, (int16_t)0x7FE1, (int16_t)0xFA6B,
    (int16_t)0x7F83, (int16_t)0xF4D8, (int16_t)0x7EE8, (int16_t)0xEF4B, (int16_t)0x7F83, (int16_t)0xF4D8, (int16_t)0x7E0E, (int16_t)0xE9C6,
    (int16_t)0x7BA3, (int16_t)0xDEDF, (int16_t)0x7EE8, (int16_t)0xEF4B, (int16_t)0x7BA3, (int16_t)0xDEDF, (int16_t)0x7642, (int16_t)0xCF04,
    (int16_t)0x7E0E, (int16_t)0xE9C6, (int16_t)0x7848, (int16_t)0xD439, (int16_t)0x6EDA, (int16_t)0xC000, (int16_t)0x7CF7, (int16_t)0xE44C,
    (int16_t)0x7402, (int16_t)0xC9E8, (int16_t)0x658D, (int16_t)0xB214, (int16_t)0x7BA3, (int16_t)0xDEDF, (int16_t)0x6EDA, (int16_t)0xC000,
    (int16_t)0x5A82, (int16_t)0xA57E, (int16_t)0x7A13, (int16_t)0xD982, (int16_t)0x68DA, (int16_t)0xB695, (int16_t)0x4DEC, (int16_t)0x9A73,
    (int16_t)0x7848, (int16_t)0xD439, (int16_t)0x620E, (int16_t)0xADB9, (int16_t)0x4000, (int16_t)0x9126, (int16_t)0x7642, (int16_t)0xCF04,
    (int16_t)0x5A82, (int16_t)0xA57E, (int16_t)0x30FC, (int16_t)0x89BE, (int16_t)0x7402, (int16_t)0xC9E8, (int16_t)0x5247, (int16_t)0x9DF2,
    (int16_t)0x2121, (int16_t)0x845D, (int16_t)0x718A, (int16_t)0xC4E5, (int16_t)0x496B, (int16_t)0x9726, (int16_t)0x10B5, (int16_t)0x8118,
    (int16_t)0x6EDA, (int16_t)0xC000, (int16_t)0x4000, (int16_t)0x9126, (int16_t)0x0000, (int16_t)0x8000, (int16_t)0x6BF4, (int16_t)0xBB3A,
    (int16_t)0x3618, (int16_t)0x8BFE, (int16_t)0xEF4B, (int16_t)0x8118, (int16_t)0x68DA, (int16_t)0xB695, (int16_t)0x2BC7, (int16_t)0x87B8,
    (int16_t)0xDEDF, (int16_t)0x845D, (int16_t)0x658D, (int16_t)0xB214, (int16_t)0x2121, (int16_t)0x845D, (int16_t)0xCF04, (int16_t)0x89BE,
    (int16_t)0x620E, (int16_t)0xADB9, (int16_t)0x163A, (int16_t)0x81F2, (int16_t)0xC000, (int16_t)0x9126, (int16_t)0x5E5F, (int16_t)0xA986,
    (int16_t)0x0B28, (int16_t)0x807D, (int16_t)0xB214, (int16_t)0x9A73, (int16_t)0x5A82, (int16_t)0xA57E, (int16_t)0x0000, (int16_t)0x8000,
    (int16_t)0xA57E, (int16_t)0xA57E, (int16_t)0x567A, (int16_t)0xA1A1, (int16_t)0xF4D8, (int16_t)0x807D, (int16_t)0x9A73, (int16_t)0xB214,
    (int16_t)0x5247, (int16_t)0x9DF2, (int16_t)0xE9C6, (int16_t)0x81F2, (int16_t)0x9126, (int16_t)0xC000, (int16_t)0x4DEC, (int16_t)0x9A73,
    (int16_t)0xDEDF, (int16_t)0x845D, (int16_t)0x89BE, (int16_t)0xCF04, (int16_t)0x496B, (int16_t)0x9726, (int16_t)0xD439, (int16_t)0x87B8,
    (int16_t)0x845D, (int16_t)0xDEDF, (int16_t)0x44C6, (int16_t)0x940C, (int16_t)0xC9E8, (int16_t)0x8BFE, (int16_t)0x8118, (int16_t)0xEF4B,
    (int16_t)0x4000, (int16_t)0x9126, (int16_t)0xC000, (int16_t)0x9126, (int16_t)0x8000, (int16_t)0x0000, (int16_t)0x3B1B, (int16_t)0x8E76,
    (int16_t)0xB695, (int16_t)0x9726, (int16_t)0x8118, (int16_t)0x10B5, (int16_t)0x3618, (int16_t)0x8BFE, (int16_t)0xADB9, (int16_t)0x9DF2,
    (int16_t)0x845D, (int16_t)0x2121, (int16_t)0x30FC, (int16_t)0x89BE, (int16_t)0xA57E, (int16_t)0xA57E, (int16_t)0x89BE, (int16_t)0x30FC,
    (int16_t)0x2BC7, (int16_t)0x87B8, (int16_t)0x9DF2, (int16_t)0xADB9, (int16_t)0x9126, (int16_t)0x4000, (int16_t)0x267E, (int16_t)0x85ED,
    (int16_t)0x9726, (int16_t)0xB695, (int16_t)0x9A73, (int16_t)0x4DEC, (int16_t)0x2121, (int16_t)0x845D, (int16_t)0x9126, (int16_t)0xC000,
    (int16_t)0xA57E, (int16_t)0x5A82, (int16_t)0x1BB4, (int16_t)0x8309, (int16_t)0x8BFE, (int16_t)0xC9E8, (int16_t)0xB214, (int16_t)0x658D,
    (int16_t)0x163A, (int16_t)0x81F2, (int16_t)0x87B8, (int16_t)0xD439, (int16_t)0xC000, (int16_t)0x6EDA, (int16_t)0x10B5, (int16_t)0x8118,
    (int16_t)0x845D, (int16_t)0xDEDF, (int16_t)0xCF04, (int16_t)0x7642, (int16_t)0x0B28, (int16_t)0x807D, (int16_t)0x81F2, (int16_t)0xE9C6,
    (int16_t)0xDEDF, (int16_t)0x7BA3, (int16_t)0x0595, (int16_t)0x801F, (int16_t)0x807D, (int16_t)0xF4D8, (int16_t)0xEF4B, (int16_t)0x7EE8,
};
#else
/* new HiFi5 packing format for the first stage
 
    N=144;
    twd = exp(-2j*pi*[1;2;3]*(0:N/4-1)/N);
    twd=twd(:);
    for k=1:N/8
        t=twd((k-1)*6+1:k*6);
        t=reshape(t,3,2).';
        t=t(:);
        twd((k-1)*6+1:k*6)=t;
    end
    twd = reshape([real(twd(:).');imag(twd(:).')],1,2*numel(twd));
    twd = bitand(65535,(double(int16(round(twd*32768.)))+65536));
    for k=1:N/8
        t=twd((k-1)*12+1:k*12);
        fprintf(1,'\t');
        for n=1:12
        fprintf(1,'(int16_t)0x%4s, ',dec2hex(t(n)));
        end
        fprintf(1,'\n');
    end
*/
{
    (int16_t)0x7FFF, (int16_t)0x0000, (int16_t)0x7FE1, (int16_t)0xFA6B, (int16_t)0x7FFF, (int16_t)0x0000, (int16_t)0x7F83, (int16_t)0xF4D8, (int16_t)0x7FFF, (int16_t)0x0000, (int16_t)0x7EE8, (int16_t)0xEF4B, 
    (int16_t)0x7F83, (int16_t)0xF4D8, (int16_t)0x7EE8, (int16_t)0xEF4B, (int16_t)0x7E0E, (int16_t)0xE9C6, (int16_t)0x7BA3, (int16_t)0xDEDF, (int16_t)0x7BA3, (int16_t)0xDEDF, (int16_t)0x7642, (int16_t)0xCF04, 
    (int16_t)0x7E0E, (int16_t)0xE9C6, (int16_t)0x7CF7, (int16_t)0xE44C, (int16_t)0x7848, (int16_t)0xD439, (int16_t)0x7402, (int16_t)0xC9E8, (int16_t)0x6EDA, (int16_t)0xC000, (int16_t)0x658D, (int16_t)0xB214, 
    (int16_t)0x7BA3, (int16_t)0xDEDF, (int16_t)0x7A13, (int16_t)0xD982, (int16_t)0x6EDA, (int16_t)0xC000, (int16_t)0x68DA, (int16_t)0xB695, (int16_t)0x5A82, (int16_t)0xA57E, (int16_t)0x4DEC, (int16_t)0x9A73, 
    (int16_t)0x7848, (int16_t)0xD439, (int16_t)0x7642, (int16_t)0xCF04, (int16_t)0x620E, (int16_t)0xADB9, (int16_t)0x5A82, (int16_t)0xA57E, (int16_t)0x4000, (int16_t)0x9126, (int16_t)0x30FC, (int16_t)0x89BE, 
    (int16_t)0x7402, (int16_t)0xC9E8, (int16_t)0x718A, (int16_t)0xC4E5, (int16_t)0x5247, (int16_t)0x9DF2, (int16_t)0x496B, (int16_t)0x9726, (int16_t)0x2121, (int16_t)0x845D, (int16_t)0x10B5, (int16_t)0x8118, 
    (int16_t)0x6EDA, (int16_t)0xC000, (int16_t)0x6BF4, (int16_t)0xBB3A, (int16_t)0x4000, (int16_t)0x9126, (int16_t)0x3618, (int16_t)0x8BFE, (int16_t)0x0000, (int16_t)0x8000, (int16_t)0xEF4B, (int16_t)0x8118, 
    (int16_t)0x68DA, (int16_t)0xB695, (int16_t)0x658D, (int16_t)0xB214, (int16_t)0x2BC7, (int16_t)0x87B8, (int16_t)0x2121, (int16_t)0x845D, (int16_t)0xDEDF, (int16_t)0x845D, (int16_t)0xCF04, (int16_t)0x89BE, 
    (int16_t)0x620E, (int16_t)0xADB9, (int16_t)0x5E5F, (int16_t)0xA986, (int16_t)0x163A, (int16_t)0x81F2, (int16_t)0x0B28, (int16_t)0x807D, (int16_t)0xC000, (int16_t)0x9126, (int16_t)0xB214, (int16_t)0x9A73, 
    (int16_t)0x5A82, (int16_t)0xA57E, (int16_t)0x567A, (int16_t)0xA1A1, (int16_t)0x0000, (int16_t)0x8000, (int16_t)0xF4D8, (int16_t)0x807D, (int16_t)0xA57E, (int16_t)0xA57E, (int16_t)0x9A73, (int16_t)0xB214, 
    (int16_t)0x5247, (int16_t)0x9DF2, (int16_t)0x4DEC, (int16_t)0x9A73, (int16_t)0xE9C6, (int16_t)0x81F2, (int16_t)0xDEDF, (int16_t)0x845D, (int16_t)0x9126, (int16_t)0xC000, (int16_t)0x89BE, (int16_t)0xCF04, 
    (int16_t)0x496B, (int16_t)0x9726, (int16_t)0x44C6, (int16_t)0x940C, (int16_t)0xD439, (int16_t)0x87B8, (int16_t)0xC9E8, (int16_t)0x8BFE, (int16_t)0x845D, (int16_t)0xDEDF, (int16_t)0x8118, (int16_t)0xEF4B, 
    (int16_t)0x4000, (int16_t)0x9126, (int16_t)0x3B1B, (int16_t)0x8E76, (int16_t)0xC000, (int16_t)0x9126, (int16_t)0xB695, (int16_t)0x9726, (int16_t)0x8000, (int16_t)0x0000, (int16_t)0x8118, (int16_t)0x10B5, 
    (int16_t)0x3618, (int16_t)0x8BFE, (int16_t)0x30FC, (int16_t)0x89BE, (int16_t)0xADB9, (int16_t)0x9DF2, (int16_t)0xA57E, (int16_t)0xA57E, (int16_t)0x845D, (int16_t)0x2121, (int16_t)0x89BE, (int16_t)0x30FC, 
    (int16_t)0x2BC7, (int16_t)0x87B8, (int16_t)0x267E, (int16_t)0x85ED, (int16_t)0x9DF2, (int16_t)0xADB9, (int16_t)0x9726, (int16_t)0xB695, (int16_t)0x9126, (int16_t)0x4000, (int16_t)0x9A73, (int16_t)0x4DEC, 
    (int16_t)0x2121, (int16_t)0x845D, (int16_t)0x1BB4, (int16_t)0x8309, (int16_t)0x9126, (int16_t)0xC000, (int16_t)0x8BFE, (int16_t)0xC9E8, (int16_t)0xA57E, (int16_t)0x5A82, (int16_t)0xB214, (int16_t)0x658D, 
    (int16_t)0x163A, (int16_t)0x81F2, (int16_t)0x10B5, (int16_t)0x8118, (int16_t)0x87B8, (int16_t)0xD439, (int16_t)0x845D, (int16_t)0xDEDF, (int16_t)0xC000, (int16_t)0x6EDA, (int16_t)0xCF04, (int16_t)0x7642, 
    (int16_t)0x0B28, (int16_t)0x807D, (int16_t)0x0595, (int16_t)0x801F, (int16_t)0x81F2, (int16_t)0xE9C6, (int16_t)0x807D, (int16_t)0xF4D8, (int16_t)0xDEDF, (int16_t)0x7BA3, (int16_t)0xEF4B, (int16_t)0x7EE8
};
#endif

/********** Twiddles table N=144 stage 2 radix 3 ******************/
ALIGN(32) static const int16_t _fft144_tw2_[] =
{
    (int16_t)0x7FFF, (int16_t)0x0000, (int16_t)0x7FFF, (int16_t)0x0000, (int16_t)0x7E0E, (int16_t)0xE9C6, (int16_t)0x7848, (int16_t)0xD439,
    (int16_t)0x7848, (int16_t)0xD439, (int16_t)0x620E, (int16_t)0xADB9, (int16_t)0x6EDA, (int16_t)0xC000, (int16_t)0x4000, (int16_t)0x9126,
    (int16_t)0x620E, (int16_t)0xADB9, (int16_t)0x163A, (int16_t)0x81F2, (int16_t)0x5247, (int16_t)0x9DF2, (int16_t)0xE9C6, (int16_t)0x81F2,
    (int16_t)0x4000, (int16_t)0x9126, (int16_t)0xC000, (int16_t)0x9126, (int16_t)0x2BC7, (int16_t)0x87B8, (int16_t)0x9DF2, (int16_t)0xADB9,
    (int16_t)0x163A, (int16_t)0x81F2, (int16_t)0x87B8, (int16_t)0xD439, (int16_t)0x0000, (int16_t)0x8000, (int16_t)0x8000, (int16_t)0x0000,
    (int16_t)0xE9C6, (int16_t)0x81F2, (int16_t)0x87B8, (int16_t)0x2BC7, (int16_t)0xD439, (int16_t)0x87B8, (int16_t)0x9DF2, (int16_t)0x5247,
};

/********** Twiddles table N=144 stage 3 radix 3 ******************/
ALIGN(32) static const int16_t _fft144_tw3_[] =
{
    (int16_t)0x7FFF, (int16_t)0x0000, (int16_t)0x7FFF, (int16_t)0x0000, (int16_t)0x6EDA, (int16_t)0xC000, (int16_t)0x4000, (int16_t)0x9126,
    (int16_t)0x4000, (int16_t)0x9126, (int16_t)0xC000, (int16_t)0x9126, (int16_t)0x0000, (int16_t)0x8000, (int16_t)0x8000, (int16_t)0x0000,
};

static const int tw_step_tab[] =
{
    1, 1, 1, 0
}; 
static const cint16ptr_t_fft tw_tab[] = 
{
    _fft144_tw1_, _fft144_tw2_, _fft144_tw3_, NULL
};

static const fn_fft_stage fft_stg_tab_s2[] = 
{
    (fn_fft_stage)fft_16x16_stage_first_scl2_DFT4,
    (fn_fft_stage)fft_16x16_stage_inner_scl2_DFT3,
    (fn_fft_stage)fft_16x16_stage_inner_scl2_DFT3,
    (fn_fft_stage)fft_16x16_stage_last_scl2_DFT4
};
#if 0
static const fn_fft_stage fft_stg_tab_s3[] = 
{
    (fn_fft_stage)NULL,//(fn_fft_stage)fft_16x16_stage_first_scl3_DFT4,
    (fn_fft_stage)NULL,//(fn_fft_stage)fft_16x16_stage_inner_scl3_DFT3,
    (fn_fft_stage)NULL,//(fn_fft_stage)fft_16x16_stage_inner_scl3_DFT3,
    (fn_fft_stage)NULL //(fn_fft_stage)fft_16x16_stage_last_scl3_DFT4
};
#endif
static const fn_fft_stage ifft_stg_tab_s2[] =
{
    (fn_fft_stage)ifft_16x16_stage_first_scl2_DFT4,
    (fn_fft_stage)fft_16x16_stage_inner_scl2_DFT3,
    (fn_fft_stage)fft_16x16_stage_inner_scl2_DFT3,
    (fn_fft_stage)ifft_16x16_stage_last_scl2_DFT4
};
#if 0
static const fn_fft_stage ifft_stg_tab_s3[] =
{
    (fn_fft_stage)NULL,//(fn_fft_stage)ifft_16x16_stage_first_scl3_DFT4,
    (fn_fft_stage)NULL,//(fn_fft_stage)fft_16x16_stage_inner_scl3_DFT3,
    (fn_fft_stage)NULL,//(fn_fft_stage)fft_16x16_stage_inner_scl3_DFT3,
    (fn_fft_stage)NULL //(fn_fft_stage)ifft_16x16_stage_last_scl3_DFT4
};
#endif
const fft_cplx_x16_descr_t __cfft_x16_descr144 =
{
    N, 
    tw_step_tab,
    tw_tab,
    NULL,
    fft_stg_tab_s2,
    NULL 
    //fft_stg_tab_s3
};     
const fft_cplx_x16_descr_t __cifft_x16_descr144 =
{
    N, 
    tw_step_tab,
    tw_tab,
    NULL,
    ifft_stg_tab_s2,
    NULL 
   // ifft_stg_tab_s3
};     
const fft_handle_t cnfft16_144 = (const fft_handle_t)&__cfft_x16_descr144;
const fft_handle_t cinfft16_144 = (const fft_handle_t)&__cifft_x16_descr144;
