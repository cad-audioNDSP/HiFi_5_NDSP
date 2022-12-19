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
  RMS calculation routines for 16-bit and 32-bit data.

  IntegrIT, 2006-2019
*/

#ifndef __RMS_H
#define __RMS_H

/* Portable data types. */
#include "types.h"

/*-------------------------------------------------------------------------
Calculate relative Root-Mean-Square (RMS) power of a real 16-bit or 32-bit
signal with a full-scale sine wave used as a reference. Result is a Q8.7
value measured in dB.

NOTE
  Functions may be also used to obtain the relative RMS power of a complex
  waveform against the full-scale complex wave.

Parameters:
  Input:
    x[N]          Real input signal
  Returned value:
                  Relative RMS power in dB, Q8.7
-------------------------------------------------------------------------*/

int16_t rms16( const int16_t * x, int N );
int16_t rms32( const int32_t * x, int N );

#endif // __RMS_H
