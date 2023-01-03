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
#include "NatureDSP_types.h"
#include "NatureDSP_Signal_img.h"
#include "imgfft_common.h"
#include "common.h"
#include "img_common.h"
#include "fft_x16_common.h"

/*-------------------------------------------------------------------------
  Image 2D-FFT makes special kind of Fourier transform over image data. 
  It supports a number of resolutions. For processing of arbitrary sizes, user 
  might pad original data to the closest available dimension.

  Important notes:
  1.  Accuracy of the fixed point 2D FFT suffers from the presence of unknown 
      overall image brightness so for the optimal results, mean pixel value 
      should be precomputed (i.e. by histogram function) and passed on the 
      input of FFT. 
  2.  FFT processes the data with autoscaling similarly as it is done in 
      regular FFT routines so it returns the scaled transform and number of 
      right shifts done during scaling procedure. Inverse FFT uses this scale 
      factor to recover original data range and saturates results to 8-bit 
      unsigned/signed or 16-bit signed data.
  3.  FFT accepts input in desired format and form output in form of 2D 
      array of complex 16-bit numbers
  4.  Supported dimensions: 64, 96, 128, 144, 176, 240, 256, 288, 320, 352, 
      384, 480, 512, 576, 640

  Image formats:
  gu8    8-bit unsigned grayscale data
  gs8    8-bit signed grayscale data
  gs16   16-bit grayscale data

  Input:
  img      input image
  sz[1]    size, see list of supported dimensions
  mean     mean value (average intensity) over the image, Q15.16
  Output:
  y[2*w*h] output spectrum. Real and imaginary data are interleaved and real 
           data goes first. w and h - width and height of original image

  Restrictions:
  img,y        should be aligned on a 16-bytes boundary
  image        dimensions should be selected from the list of supported 
               dimensions
  image stride should be a multiple of 8 for 8-bit images and 4 for 16-bit 
               images
-------------------------------------------------------------------------*/
/* forward transform */
int imgfft_gu8 (void *pScr, 
               complex_fract16* y, const void* img, 
               int32_t mean, const imgsize_t* sz)
{
    NASSERT(pScr);
    NASSERT(y);
    NASSERT(img);
    NASSERT_ALIGN(pScr,HIFI_SIMD_WIDTH);
    NASSERT_ALIGN(y   ,HIFI_SIMD_WIDTH);
    NASSERT_ALIGN(img ,HIFI_SIMD_WIDTH);
    imgsize_validate(sz,1,1);
    NASSERT(sz->width== 64 || sz->width== 96 || sz->width==128 || sz->width==144 || 
            sz->width==176 || sz->width==240 || sz->width==256 || sz->width==288 || 
            sz->width==320 || sz->width==352 || sz->width==384 || sz->width==480 || 
            sz->width==512 || sz->width==576 || sz->width==640);
    NASSERT(sz->height== 64 || sz->height== 96 || sz->height==128 || sz->height==144 || 
            sz->height==176 || sz->height==240 || sz->height==256 || sz->height==288 || 
            sz->height==320 || sz->height==352 || sz->height==384 || sz->height==480 || 
            sz->height==512 || sz->height==576 || sz->height==640);
    /*
        1. copy image row by row to the temporary buffer with 16-bit conversion, 
           remove mean value and make the row FFTs
        2. find common scale for all ffts
        3. rescale results
        4. copy columns by colums and make column FFTs
        5. find common scale for all ffts
        6. rescale results 
    */
    int M = sz->height;
    int N = sz->width; 
    int stride = sz->stride; 
    int16_t dcBias = (int16_t)(mean >> 9);
    int i, j, expRows, expCols;
    int16_t *tmpX2 = (int16_t*)y;
    uint8_t *x = (uint8_t*)img; 
    const signed char *pIn = (const signed char *)img;
    ae_int16x8 * restrict pOut;
    ae_int16x4 t0, t1; 
    /* Convert uint8_t samples to int16 and subtract the mean value */
    for (i = 0; i < M; i++)
    {
        pIn = (const signed char  *)(x + i*stride);
        pOut = (ae_int16x8 *)(tmpX2 + 2 * i*N); 
        for (j = 0; j < (N>>3); j++)
        {
            //tmpX2[2 * i*N + j] = ((int16_t)x[i*stride + j] << 7) - dcBias;
            AE_L8X4F_IP(t0, pIn, 4 * sizeof(int8_t));
            AE_L8X4F_IP(t1, pIn, 4 * sizeof(int8_t));
            t0 = AE_MOVINT16X4_FROMINT32X2(AE_SRLI32(AE_MOVINT32X2_FROMINT16X4(t0), 1));
            t1 = AE_MOVINT16X4_FROMINT32X2(AE_SRLI32(AE_MOVINT32X2_FROMINT16X4(t1), 1));
            t0 = AE_SUB16S(t0, dcBias); 
            t1 = AE_SUB16S(t1, dcBias); 
            AE_S16X4X2_IP(t0, t1, pOut, 2*sizeof(t0)); 
        }
    }

    expRows = ApplyRFFT_toRows((complex_fract16*)pScr, y, (int16_t*)y, M, N);
    expCols = ApplyFFT_toCols ((complex_fract16*)pScr, y, y, M, N);

    return expCols + expRows;
}

/* Request the scratch size:*/ 
size_t imgfft_gu8_getScratchSize   (const imgsize_t* sz)
{
    int M = sz->height;
    int N = sz->width;
    imgsize_validate(sz,1,1);

    return sizeof(complex_fract16) * XT_MAX(N+M, 4*M+N/2+1);
}
