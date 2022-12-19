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
 * Aligning wrapper over standard malloc/free.
 */

#ifndef MALLOC16_H__
#define MALLOC16_H__

/* Portable data types. */
#include "types.h"

/*
 Allocate and align a memory area.
 Input:
   sz     Area size, in bytes
   align  A power of 2 number specifying the alignment boundary, in bytes.
          Set it to zero to select the default alignment, that is the SIMD
          vector width for the target platform.
 Returns:
   Pointer to the aligned memory area if succeeded, or null pointer
   if failed to allocate the memory. 
*/
void * mallocAlign(size_t sz, size_t align);

/*
 Free an aligned memory area.
 Input:
   ptr  Pointer to a memory area allocated with mallocAlign()
 Output:
   None
*/
void freeAlign(void * ptr);

/*-------------------------------------------------------------------------------
additional functions allowing to allocate memory not from heap but from 
specified region
For allocation, user should map the region by mallocAlignInitRegion() and 
call mallocAlignRegion() instead of mallocAlign
-------------------------------------------------------------------------------*/
void mallocAlignInitRegion(void* region, size_t sz);
void * mallocAlignRegion(size_t sz, size_t align);


#endif
