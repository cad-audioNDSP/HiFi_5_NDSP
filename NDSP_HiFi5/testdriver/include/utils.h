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
 * Utility finctions
 */

#ifndef __UTILS_H__
#define __UTILS_H__

#ifdef __cplusplus
extern "C" {
#endif

/* Portable data types */
#include "types.h"

#ifndef PI
    #define PI  3.1415926535897932384626433832795
#endif

#ifdef _MSC_VER
    #define ALIGN(x)    _declspec(align(x)) 
#else
    #define ALIGN(x)    __attribute__((aligned(x))) 
#endif

void dump_open (const char *pFileName);
void dump_close();
void dump_i16  (const void *x, int size, const char *info);
void dump_i24  (const void *x, int size, const char *info);
void dump_f24  (const void *x, int size, const char *info);
void dump_i32  (const void *x, int size, const char *info);

#define RAND_RESET_A 0
#define RAND_RESET_B 0
void Rand_reset(int a, int b);
int  Rand(void);
void Rand_i32(void *x, int n);
void Rand_i16(void *x, int n);
void Rand_f32(void *x, int n); /* [-1.f,1.f] */

/* 32-bit CRC update for a data block. */
uint32_t crc32( uint32_t crc, const uint8_t * buf, int len );

#ifdef __cplusplus
};
#endif
#endif/*__UTILS_H__*/
