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
 * Test environamnt configuration.
 * Defines the abstraction layer between generic components of the test 
 * environamnet and the library, such that the test environment can refer
 * to general functions of the library (identification, diagnostics, etc) in
 * a unified way.
 */

#ifndef __CONFIG_H
#define __CONFIG_H

#define LIBRARY_NAME           "NatureDSP Signal"
/* Library API */
#define LIBRARY_HEADER(pkg)    LIBRARY_HEADER_##pkg
/* General APIs of the library */
#define LIBRARY_HEADER_types   "NatureDSP_types.h"          /* Portable data types                          */
#define LIBRARY_HEADER_id      "NatureDSP_Signal_id.h"      /* Library identification                       */
#define LIBRARY_HEADER_diag    "NatureDSP_Signal_diag.h"    /* Library diagnostics                          */
/* Functional APIs of the library */
#define LIBRARY_HEADER_fir     "NatureDSP_Signal_fir.h"     /* FIR Filters and Related Functions            */
#define LIBRARY_HEADER_iir     "NatureDSP_Signal_iir.h"     /* IIR Filters                                  */
#define LIBRARY_HEADER_fft     "NatureDSP_Signal_fft.h"     /* FFT Routines                                 */
#define LIBRARY_HEADER_matop   "NatureDSP_Signal_matop.h"   /* Matrix Operations                            */
#define LIBRARY_HEADER_matinv  "NatureDSP_Signal_matinv.h"  /* Matrix Decomposition and Inversion Functions */
#define LIBRARY_HEADER_vector  "NatureDSP_Signal_vector.h"  /* Vector Operations                            */
#define LIBRARY_HEADER_math    "NatureDSP_Signal_math.h"    /* Math Functions                               */
#define LIBRARY_HEADER_complex "NatureDSP_Signal_complex.h" /* Complex Math Functions                       */
#define LIBRARY_HEADER_fit     "NatureDSP_Signal_fit.h"     /* Fitting and Interpolation Routines           */
#define LIBRARY_HEADER_audio   "NatureDSP_Signal_audio.h"   /* Audio processing functions                   */
#define LIBRARY_HEADER_img     "NatureDSP_Signal_img.h"     /* Image processing functions                   */

/* Symbol name for a text annotation for a library function. */
#define ANNOTATE_FUN_REF(fun)   NatureDSP_Signal_annotation_##fun

#include LIBRARY_HEADER(id)
#include LIBRARY_HEADER(diag)

/* Library identification */
#define GET_LIBRARY_VERSION(version_string)       NatureDSP_Signal_get_library_version(version_string)
#define GET_LIBRARY_API_VERSION(version_string)   NatureDSP_Signal_get_library_api_version(version_string)
#define IS_PRESENT(fun)                           NatureDSP_Signal_isPresent(fun)
/* Hardware capabilities info */
#define GET_ISA_OPT(opt)                          NatureDSP_Signal_get_isa_opt(NATUREDSP_ISA_OPT_##opt)
/* Floating-point exceptions control */
#define FECLEAREXCEPT(excepts)                    NatureDSP_Signal_feclearexcept(excepts)
#define FERAISEEXCEPT(excepts)                    NatureDSP_Signal_feraiseexcept(excepts)
#define FETESTEXCEPT(excepts)                     NatureDSP_Signal_fetestexcept(excepts)

/*
 * Calculate the number of data elements to store a block-ordered matrix, taking into account
 * the optimal memory alignment featured by the target platform. 
 * For a library with single-matrix, replace the definitions below with this stub:
 *  #define MATALLOCMXNN(elemSize,M,N)    (M)*(N)
 */

#define MATALLOCMXNN(elemSize,M,N)    (M)*(N)

#endif /* __CONFIG_H */
