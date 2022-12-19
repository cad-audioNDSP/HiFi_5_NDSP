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
* Test module for testing cycle performance (Scalar Math)
*/

#include "mips.h"
#include "config.h"
#include "packages.h"
#include LIBRARY_HEADER(math)

void mips_maths1(int isFull, int isVerbose, FILE * fout)
{
  PROFILE_SIMPLE(isFull, isVerbose, scl_recip16x16,  (7002 ),                    fout,"",prf_cycle);
  PROFILE_SIMPLE(isFull, isVerbose, scl_recip32x32,  (7002 ),                    fout,"",prf_cycle);
#if 0 //HiFi3/3z API
  PROFILE_SIMPLE(isFull, isVerbose, scl_recip24x24,  (-1966),                    fout,"",prf_cycle);
#endif
  PROFILE_SIMPLE(isFull, isVerbose, scl_recip64x64,  (7002 ),                    fout,"",prf_cycle);
  PROFILE_SIMPLE(isFull, isVerbose, scl_divide16x16, (-17621, -29508),           fout,"",prf_cycle);
  PROFILE_SIMPLE(isFull, isVerbose, scl_divide32x32, (-1154751292, -1933789767), fout,"",prf_cycle);
#if 0 //HiFi3/3z API
  PROFILE_SIMPLE(isFull, isVerbose, scl_divide24x24, (1154751292, -1933789767),  fout,"",prf_cycle);
#endif
  PROFILE_SIMPLE(isFull, isVerbose, scl_divide64x32, (-1154751292, -1933789767), fout,"",prf_cycle);
  PROFILE_SIMPLE(isFull, isVerbose, scl_divide64x64, (-1154751292, -1933789767), fout,"",prf_cycle);
  PROFILE_SIMPLE(isFull, isVerbose, scl_log2_32x32, (496366179),                 fout,"",prf_cycle);
  PROFILE_SIMPLE(isFull, isVerbose, scl_logn_32x32, (496366179),                 fout,"",prf_cycle);
  PROFILE_SIMPLE(isFull, isVerbose, scl_log10_32x32, (496366179),                fout,"",prf_cycle);
#if 0 //HiFi3/3z API
  PROFILE_SIMPLE(isFull, isVerbose, scl_log2_24x24, (496366179),                 fout,"",prf_cycle);
  PROFILE_SIMPLE(isFull, isVerbose, scl_logn_24x24, (496366179),                 fout,"",prf_cycle);
  PROFILE_SIMPLE(isFull, isVerbose, scl_log10_24x24, (496366179),                fout,"",prf_cycle);
#endif
  PROFILE_SIMPLE(isFull, isVerbose, scl_antilog2_32x32, (-1010430329),           fout,"",prf_cycle);
  PROFILE_SIMPLE(isFull, isVerbose, scl_antilogn_32x32, (-1010430329),           fout,"",prf_cycle);
  PROFILE_SIMPLE(isFull, isVerbose, scl_antilog10_32x32, (-1010430329),          fout,"",prf_cycle);
#if 0 //HiFi3/3z API
  PROFILE_SIMPLE(isFull, isVerbose, scl_antilog2_24x24, (-1010430329),           fout,"",prf_cycle);
  PROFILE_SIMPLE(isFull, isVerbose, scl_antilogn_24x24, (-1010430329),           fout,"",prf_cycle);
  PROFILE_SIMPLE(isFull, isVerbose, scl_antilog10_24x24, (-1010430329),          fout,"",prf_cycle);
#endif
  PROFILE_SIMPLE(isFull, isVerbose, scl_sqrt16x16, (-29508),                     fout,"",prf_cycle);
  PROFILE_SIMPLE(isFull, isVerbose, scl_sqrt32x16, (-1154751292),                fout,"",prf_cycle);
  PROFILE_SIMPLE(isFull, isVerbose, scl_sqrt32x32, (-1154751292),                fout,"",prf_cycle);
#if 0 //HiFi3/3z API
  PROFILE_SIMPLE(isFull, isVerbose, scl_sqrt24x24, (-1154751292),                fout,"",prf_cycle);
#endif
  PROFILE_SIMPLE(isFull, isVerbose, scl_sqrt64x32, (-1105678961268363246LL),     fout,"",prf_cycle);
  PROFILE_SIMPLE(isFull, isVerbose, scl_sine32x32, (-1154751292),                fout,"",prf_cycle);
  PROFILE_SIMPLE(isFull, isVerbose, scl_cosine32x32, (-1154751292),              fout,"",prf_cycle);
#if 0 //HiFi3/3z API
  PROFILE_SIMPLE(isFull, isVerbose, scl_sine24x24, (-1154751292),                fout,"",prf_cycle);
  PROFILE_SIMPLE(isFull, isVerbose, scl_cosine24x24, (-1154751292),              fout,"",prf_cycle);
#endif
  PROFILE_SIMPLE(isFull, isVerbose, scl_tan32x32, (2147483640),                  fout,"",prf_cycle);
#if 0 //HiFi3/3z API
  PROFILE_SIMPLE(isFull, isVerbose, scl_tan24x24, (2147483640),                  fout,"",prf_cycle);
#endif
  PROFILE_SIMPLE(isFull, isVerbose, scl_atan32x32, (-1154751292),                fout,"",prf_cycle);
#if 0 //HiFi3/3z API
  PROFILE_SIMPLE(isFull, isVerbose, scl_atan24x24, (-1154751292),                fout,"",prf_cycle);
  PROFILE_SIMPLE(isFull, isVerbose, scl_atan2_24x24, (-1154751292, -1010430329), fout,"",prf_cycle);
#endif
  PROFILE_SIMPLE(isFull, isVerbose, scl_rsqrt16x16, (29508),                    fout,"",prf_cycle);
  PROFILE_SIMPLE(isFull, isVerbose, scl_rsqrt32x32, (1154751292),               fout,"",prf_cycle);
  PROFILE_SIMPLE(isFull, isVerbose, scl_sigmoid32x32, (-1154751292),             fout,"",prf_cycle);
  PROFILE_SIMPLE(isFull, isVerbose, scl_tanh32x32, (-1154751292),                fout,"",prf_cycle);
  PROFILE_SIMPLE(isFull, isVerbose, scl_relu32x32, (1154751292, 1010430329),     fout,"",prf_cycle);
}
