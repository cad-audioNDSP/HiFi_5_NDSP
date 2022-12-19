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

/********** Twiddles table N=88 stage 1 radix 4 ******************/
ALIGN(32) static const int16_t _fft88_tw1_[] =
/* older HiFi4 packing:
{
    (int16_t)0x7FFF, (int16_t)0x0000, (int16_t)0x7FFF, (int16_t)0x0000, (int16_t)0x7FFF, (int16_t)0x0000,
    (int16_t)0x7FAD, (int16_t)0xF6DE, (int16_t)0x7EB2, (int16_t)0xEDC9, (int16_t)0x7D13, (int16_t)0xE4CB,
    (int16_t)0x7EB2, (int16_t)0xEDC9, (int16_t)0x7AD1, (int16_t)0xDBF0, (int16_t)0x746F, (int16_t)0xCAD4,
    (int16_t)0x7D13, (int16_t)0xE4CB, (int16_t)0x746F, (int16_t)0xCAD4, (int16_t)0x6678, (int16_t)0xB34B,
    (int16_t)0x7AD1, (int16_t)0xDBF0, (int16_t)0x6BAE, (int16_t)0xBACC, (int16_t)0x53D2, (int16_t)0x9F44,
    (int16_t)0x77EE, (int16_t)0xD345, (int16_t)0x60BC, (int16_t)0xAC2E, (int16_t)0x3D58, (int16_t)0x8FA8,
    (int16_t)0x746F, (int16_t)0xCAD4, (int16_t)0x53D2, (int16_t)0x9F44, (int16_t)0x2410, (int16_t)0x852F,
    (int16_t)0x7058, (int16_t)0xC2A8, (int16_t)0x4534, (int16_t)0x9452, (int16_t)0x0922, (int16_t)0x8053,
    (int16_t)0x6BAE, (int16_t)0xBACC, (int16_t)0x352C, (int16_t)0x8B91, (int16_t)0xEDC9, (int16_t)0x814E,
    (int16_t)0x6678, (int16_t)0xB34B, (int16_t)0x2410, (int16_t)0x852F, (int16_t)0xD345, (int16_t)0x8812,
    (int16_t)0x60BC, (int16_t)0xAC2E, (int16_t)0x1237, (int16_t)0x814E, (int16_t)0xBACC, (int16_t)0x9452,
    (int16_t)0x5A82, (int16_t)0xA57E, (int16_t)0x0000, (int16_t)0x8000, (int16_t)0xA57E, (int16_t)0xA57E,
    (int16_t)0x53D2, (int16_t)0x9F44, (int16_t)0xEDC9, (int16_t)0x814E, (int16_t)0x9452, (int16_t)0xBACC,
    (int16_t)0x4CB5, (int16_t)0x9988, (int16_t)0xDBF0, (int16_t)0x852F, (int16_t)0x8812, (int16_t)0xD345,
    (int16_t)0x4534, (int16_t)0x9452, (int16_t)0xCAD4, (int16_t)0x8B91, (int16_t)0x814E, (int16_t)0xEDC9,
    (int16_t)0x3D58, (int16_t)0x8FA8, (int16_t)0xBACC, (int16_t)0x9452, (int16_t)0x8053, (int16_t)0x0922,
    (int16_t)0x352C, (int16_t)0x8B91, (int16_t)0xAC2E, (int16_t)0x9F44, (int16_t)0x852F, (int16_t)0x2410,
    (int16_t)0x2CBB, (int16_t)0x8812, (int16_t)0x9F44, (int16_t)0xAC2E, (int16_t)0x8FA8, (int16_t)0x3D58,
    (int16_t)0x2410, (int16_t)0x852F, (int16_t)0x9452, (int16_t)0xBACC, (int16_t)0x9F44, (int16_t)0x53D2,
    (int16_t)0x1B35, (int16_t)0x82ED, (int16_t)0x8B91, (int16_t)0xCAD4, (int16_t)0xB34B, (int16_t)0x6678,
    (int16_t)0x1237, (int16_t)0x814E, (int16_t)0x852F, (int16_t)0xDBF0, (int16_t)0xCAD4, (int16_t)0x746F,
    (int16_t)0x0922, (int16_t)0x8053, (int16_t)0x814E, (int16_t)0xEDC9, (int16_t)0xE4CB, (int16_t)0x7D13,
};*/
{
/* newer HiFi5 packing
    N=176;
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
    (int16_t)0x7FFF, (int16_t)0x0000, (int16_t)0x7FAD, (int16_t)0xF6DE, (int16_t)0x7FFF, (int16_t)0x0000, (int16_t)0x7EB2, (int16_t)0xEDC9, (int16_t)0x7FFF, (int16_t)0x0000, (int16_t)0x7D13, (int16_t)0xE4CB, 
    (int16_t)0x7EB2, (int16_t)0xEDC9, (int16_t)0x7D13, (int16_t)0xE4CB, (int16_t)0x7AD1, (int16_t)0xDBF0, (int16_t)0x746F, (int16_t)0xCAD4, (int16_t)0x746F, (int16_t)0xCAD4, (int16_t)0x6678, (int16_t)0xB34B, 
    (int16_t)0x7AD1, (int16_t)0xDBF0, (int16_t)0x77EE, (int16_t)0xD345, (int16_t)0x6BAE, (int16_t)0xBACC, (int16_t)0x60BC, (int16_t)0xAC2E, (int16_t)0x53D2, (int16_t)0x9F44, (int16_t)0x3D58, (int16_t)0x8FA8, 
    (int16_t)0x746F, (int16_t)0xCAD4, (int16_t)0x7058, (int16_t)0xC2A8, (int16_t)0x53D2, (int16_t)0x9F44, (int16_t)0x4534, (int16_t)0x9452, (int16_t)0x2410, (int16_t)0x852F, (int16_t)0x0922, (int16_t)0x8053, 
    (int16_t)0x6BAE, (int16_t)0xBACC, (int16_t)0x6678, (int16_t)0xB34B, (int16_t)0x352C, (int16_t)0x8B91, (int16_t)0x2410, (int16_t)0x852F, (int16_t)0xEDC9, (int16_t)0x814E, (int16_t)0xD345, (int16_t)0x8812, 
    (int16_t)0x60BC, (int16_t)0xAC2E, (int16_t)0x5A82, (int16_t)0xA57E, (int16_t)0x1237, (int16_t)0x814E, (int16_t)0x0000, (int16_t)0x8000, (int16_t)0xBACC, (int16_t)0x9452, (int16_t)0xA57E, (int16_t)0xA57E, 
    (int16_t)0x53D2, (int16_t)0x9F44, (int16_t)0x4CB5, (int16_t)0x9988, (int16_t)0xEDC9, (int16_t)0x814E, (int16_t)0xDBF0, (int16_t)0x852F, (int16_t)0x9452, (int16_t)0xBACC, (int16_t)0x8812, (int16_t)0xD345, 
    (int16_t)0x4534, (int16_t)0x9452, (int16_t)0x3D58, (int16_t)0x8FA8, (int16_t)0xCAD4, (int16_t)0x8B91, (int16_t)0xBACC, (int16_t)0x9452, (int16_t)0x814E, (int16_t)0xEDC9, (int16_t)0x8053, (int16_t)0x0922, 
    (int16_t)0x352C, (int16_t)0x8B91, (int16_t)0x2CBB, (int16_t)0x8812, (int16_t)0xAC2E, (int16_t)0x9F44, (int16_t)0x9F44, (int16_t)0xAC2E, (int16_t)0x852F, (int16_t)0x2410, (int16_t)0x8FA8, (int16_t)0x3D58, 
    (int16_t)0x2410, (int16_t)0x852F, (int16_t)0x1B35, (int16_t)0x82ED, (int16_t)0x9452, (int16_t)0xBACC, (int16_t)0x8B91, (int16_t)0xCAD4, (int16_t)0x9F44, (int16_t)0x53D2, (int16_t)0xB34B, (int16_t)0x6678, 
    (int16_t)0x1237, (int16_t)0x814E, (int16_t)0x0922, (int16_t)0x8053, (int16_t)0x852F, (int16_t)0xDBF0, (int16_t)0x814E, (int16_t)0xEDC9, (int16_t)0xCAD4, (int16_t)0x746F, (int16_t)0xE4CB, (int16_t)0x7D13
};

/********** Twiddles table N=88 stage 2 radix 2 ******************/
ALIGN(32) static const int16_t _fft88_tw2_[] =
{
    (int16_t)0x7FFF, (int16_t)0x0000,
    (int16_t)0x7AD1, (int16_t)0xDBF0,
    (int16_t)0x6BAE, (int16_t)0xBACC,
    (int16_t)0x53D2, (int16_t)0x9F44,
    (int16_t)0x352C, (int16_t)0x8B91,
    (int16_t)0x1237, (int16_t)0x814E,
    (int16_t)0xEDC9, (int16_t)0x814E,
    (int16_t)0xCAD4, (int16_t)0x8B91,
    (int16_t)0xAC2E, (int16_t)0x9F44,
    (int16_t)0x9452, (int16_t)0xBACC,
    (int16_t)0x852F, (int16_t)0xDBF0,
};

/* Twiddles tables for fft_real32x16, ifft_real32x16, fft_real16x16, ifft_real16x16, N=176 */
ALIGN(32) static const int16_t __fft_real16_tw176[88] =
{
    (int16_t)0x0000, (int16_t)0x7FFF, (int16_t)0x0492, (int16_t)0x7FEB, (int16_t)0x0922, (int16_t)0x7FAD, (int16_t)0x0DAF, (int16_t)0x7F44,
    (int16_t)0x1237, (int16_t)0x7EB2, (int16_t)0x16BA, (int16_t)0x7DF7, (int16_t)0x1B35, (int16_t)0x7D13, (int16_t)0x1FA8, (int16_t)0x7C06,
    (int16_t)0x2410, (int16_t)0x7AD1, (int16_t)0x286C, (int16_t)0x7973, (int16_t)0x2CBB, (int16_t)0x77EE, (int16_t)0x30FC, (int16_t)0x7642,
    (int16_t)0x352C, (int16_t)0x746F, (int16_t)0x394C, (int16_t)0x7276, (int16_t)0x3D58, (int16_t)0x7058, (int16_t)0x4151, (int16_t)0x6E15,
    (int16_t)0x4534, (int16_t)0x6BAE, (int16_t)0x4900, (int16_t)0x6924, (int16_t)0x4CB5, (int16_t)0x6678, (int16_t)0x5051, (int16_t)0x63AB,
    (int16_t)0x53D2, (int16_t)0x60BC, (int16_t)0x5739, (int16_t)0x5DAF, (int16_t)0x5A82, (int16_t)0x5A82, (int16_t)0x5DAF, (int16_t)0x5739,
    (int16_t)0x60BC, (int16_t)0x53D2, (int16_t)0x63AB, (int16_t)0x5051, (int16_t)0x6678, (int16_t)0x4CB5, (int16_t)0x6924, (int16_t)0x4900,
    (int16_t)0x6BAE, (int16_t)0x4534, (int16_t)0x6E15, (int16_t)0x4151, (int16_t)0x7058, (int16_t)0x3D58, (int16_t)0x7276, (int16_t)0x394C,
    (int16_t)0x746F, (int16_t)0x352C, (int16_t)0x7642, (int16_t)0x30FC, (int16_t)0x77EE, (int16_t)0x2CBB, (int16_t)0x7973, (int16_t)0x286C,
    (int16_t)0x7AD1, (int16_t)0x2410, (int16_t)0x7C06, (int16_t)0x1FA8, (int16_t)0x7D13, (int16_t)0x1B35, (int16_t)0x7DF7, (int16_t)0x16BA,
    (int16_t)0x7EB2, (int16_t)0x1237, (int16_t)0x7F44, (int16_t)0x0DAF, (int16_t)0x7FAD, (int16_t)0x0922, (int16_t)0x7FEB, (int16_t)0x0492,
};
static const int tw_step_tab[] =
{
    1, 1, 0
};
static const cint16ptr_t_fft tw_tab[] =
{
    _fft88_tw1_, _fft88_tw2_,  NULL
};


static const fn_fft_stage fft_stg_tab_s2[] = 
{
    (fn_fft_stage)fft_16x16_stage_first_scl2_DFT4,
    (fn_fft_stage)fft_16x16_stage_inner_scl2_DFT2_v4,
    (fn_fft_stage)fft_16x16_stage_last_scl2_DFT11
};
static const fn_fft_stage fft_stg_tab_s3[] = 
{
    (fn_fft_stage)fft_16x16_stage_first_scl3_DFT4,
    NULL, // not implemented (fn_fft_stage)fft_16x16_stage_inner_scl3_DFT2_v3,
    NULL  // not implemented (fn_fft_stage)fft_16x16_stage_last_scl3_DFT11
};
static const fn_fft_stage ifft_stg_tab_s2[] =
{
    (fn_fft_stage)ifft_16x16_stage_first_scl2_DFT4,
    (fn_fft_stage)fft_16x16_stage_inner_scl2_DFT2_v4,
    (fn_fft_stage)ifft_16x16_stage_last_scl2_DFT11
};
static const fn_fft_stage ifft_stg_tab_s3[] =
{
    (fn_fft_stage)ifft_16x16_stage_first_scl3_DFT4,
    NULL, // not implemented (fn_fft_stage)fft_16x16_stage_inner_scl3_DFT2_v3,
    NULL  // not implemented (fn_fft_stage)ifft_16x16_stage_last_scl3_DFT11
};
static const fft_cplx_x16_descr_t __cfft_x16_descr88 =
{
    88,
    tw_step_tab,
    tw_tab,
    NULL,
    fft_stg_tab_s2,
    NULL, 
    fft_stg_tab_s3
};
static const fft_cplx_x16_descr_t __cifft_x16_descr88 =
{
    88,
    tw_step_tab,
    tw_tab,
    NULL,
    ifft_stg_tab_s2,
    NULL, 
    ifft_stg_tab_s3
};


static const fft_real_x16_descr_t __rfft_descr =
{
    &__cfft_x16_descr88,
    __fft_real16_tw176
};
static const fft_real_x16_descr_t __rifft_descr =
{
    &__cifft_x16_descr88,
    __fft_real16_tw176
};
const fft_handle_t rnfft16_176 =  (const fft_handle_t)&__rfft_descr;
const fft_handle_t rinfft16_176 = (const fft_handle_t)&__rifft_descr;
