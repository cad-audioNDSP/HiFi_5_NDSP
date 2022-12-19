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
 * Test data vectors tools and SEQ-file reader
 */

#ifndef __VECTOOLS_H
#define __VECTOOLS_H

#include <stdio.h>

/* Cross-platform data type definitions. */
#include "types.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Data vector format. Basic data type (e.g. FMT_FRACT16) may be ORed 
 * with FMT_REAL/FMT_CPLX. */
#define FMT_DTYPE_MASK   0x0f
#define FMT_UINT8        0x00
#define FMT_INT8         0x01
#define FMT_INT16        0x02
#define FMT_INT32        0x03
#define FMT_INT64        0x04
#define FMT_FRACT16      0x05
#define FMT_FRACT32      0x06
#define FMT_FLOAT16      0x07
#define FMT_FLOAT32      0x08
#define FMT_FLOAT64      0x09
#define FMT_EF32X16      0x0a      // emulated floating point (32b mantissa/16exponent) packed in 64b word
#define FMT_DOMAIN_MASK  0x30
#define FMT_REAL         0x10
#define FMT_CPLX         0x20

/* alignment/region allocation flags (isAligned in vecAlloc,vecsAlloc ) */
#define ALLOC_UNALIGNED  0x0000  /* unaligned alloc (alignment on element size) */
#define ALLOC_ALIGNED    0x0001  /* aligned on system natural alignment         */
#define ALLOC_REGIONMASK 0xf000  /* mask for region flags */
#define ALLOC_HEAP       0x0000  /* allocate in heap                            */
#define ALLOC_DRAM0      0x1000  /* allocate in dram0 (requires calling of 
                                   vecInit first)                              */

/* Basic data vector structure. */
typedef struct tagVec
{
  void *    ptr;    /* Memory block pointer                                     */
  int       fmt;    /* Vector format, see FMT_ flags                            */
  int       offset; /* First elem offset to emulate the unaligned data pointer  */
  uint32_t  gzCrc;  /* CRC for guard zones                                      */
  int       nElem;  /* Number of data items in the vector                       */
  size_t    szElem; /* Data item size, in bytes.                                */
  size_t    szBulk; /* Allocated memory block size, in bytes                    */
}
tVec;

/*------------------------------------------------------------------------------*
 *                     Vector allocation and deletion                           *
 *------------------------------------------------------------------------------*/

/* Return data element size for specified format, in bytes. */
size_t vecElemSize(int fmt);

/* Allocate memory for a single vector. Last argument, if not zero, is a pointer
 * to vector initialization data. Function returns zero if failed. 

    alignmentFlags see ALLOC_xxx 
*/
int vecAlloc( tVec * vec, int nElem, int alignmentFlags, int fmt, const void * initData );

/* Allocate memory for a few vectors of the same data format. Repeat last two
 * arguments for each new vector instance. The argument list must end up 
 * with zero pointer. Returns the number of successfully allocated vectors. 

   alignmentFlags see ALLOC_xxx 
 */
int vecsAlloc( int alignmentFlags, int fmt, tVec * vec, int nElem, ... );

/* Make a copy of the input vector. For an unaligned input vector, the replicated
 * vector will be also unaligned, possibly with different offset from the nearest
 * aligned address. Function return zero if failed. */
int vecClone( tVec * vecOut, const tVec * vecIn );

/* Free memory for a single vector. Return zero if guard memory areas check fails. */
int vecFree( tVec * vec );

/* Free memory for a list of vectors. The last argument must be zero pointer.
 * Return zero if guard memory areas check fails for any vector. */
int vecsFree( tVec * vec, ... );

/* 
    initialize allocation region

    alignmentFlags see ALLOC_xxx 

    NOTE:
    for allocation in heap this function is not neccesary
    for allocation in dram this function should be called before
*/
void vecInitRegion(int alignmentFlags);
/*------------------------------------------------------------------------------*
 *                       Vector data access functions                           *
 *------------------------------------------------------------------------------*/

/* Return the total size of vector data items, in bytes. */
size_t vecGetSize( const tVec * vec );

/* Return a pointer to element vec[n]. */
void * vecGetElem( const tVec * vec, int n );

/* Return a pointer to element vec[n] with appropriate type casting. */
uint8_t         * vecGetElem_u8   ( const tVec * vec, int n );
int8_t          * vecGetElem_i8   ( const tVec * vec, int n );
int16_t         * vecGetElem_i16  ( const tVec * vec, int n );
int32_t         * vecGetElem_i32  ( const tVec * vec, int n );
int64_t         * vecGetElem_i64  ( const tVec * vec, int n );
uint64_t        * vecGetElem_ef32x16  ( const tVec * vec, int n );
fract16         * vecGetElem_fr16 ( const tVec * vec, int n );
fract32         * vecGetElem_fr32 ( const tVec * vec, int n );
float16_t       * vecGetElem_fl16 ( const tVec * vec, int n );
float32_t       * vecGetElem_fl32 ( const tVec * vec, int n );
float64_t       * vecGetElem_fl64 ( const tVec * vec, int n );
complex_float16    * vecGetElem_fl16c( const tVec * vec, int n );
complex_float   * vecGetElem_fl32c( const tVec * vec, int n );
complex_double  * vecGetElem_fl64c( const tVec * vec, int n );
complex_fract16 * vecGetElem_fr16c( const tVec * vec, int n );
complex_fract32 * vecGetElem_fr32c( const tVec * vec, int n );

/* Print the velue of the specified element to a string buffer. */
void vecPrintElem( char * buf, const tVec * vec, int n );

/*------------------------------------------------------------------------------*
 *                           Data conversions                                   *
 *------------------------------------------------------------------------------*/

/* Convert 16-bit PCM samples to vector data format. When the target format is
 * fractional (either Q15 or Q31), data are normalized and shifted to the right
 * by bexp bit positions so that block exponent of output data equals bexp.
 * Return the overall shift amount, with positive values corresponding to the
 * left shift. */
int vecFromPcm16( tVec * vec, const fract16 * x, int bexp );
/* Blockwise variant. */
void bvecFromPcm16( tVec * vec, const fract16 * x, int N, int L, int bexp, int16_t shift[] );
/* Streamwise variant. */
void svecFromPcm16(tVec * vec, const fract16 * x, int N, int L, int bexp, int16_t shift[]);

/* Convert 64-bit floating point values to vector data format. When the target
 * format is fractional, we assume either Q15 or Q31. */
void vecFromFp64( tVec * vec, const float64_t * x );

/* Convert vector data to 64-bit floating point numbers. If the source vector
 * is of a fractional format (we assume Q15 or Q31), then data are additionally
 * scaled by 2^-shift. */
void vecToFp64( float64_t * z, const tVec * vec, int shift );
/* Blockwise variant. */
void bvecToFp64( float64_t * z, const tVec * vec, int N, int L, int16_t shift[] );
/* Streamwise variant. */
void svecToFp64(float64_t * z, const tVec * vec, int N, int L, int16_t shift[]);

/*------------------------------------------------------------------------------*
 *                  Vectors comparison and data validation                      *
 *------------------------------------------------------------------------------*/

/* Verify that every vector element belongs to the closed range defined by the lower
 * and upper bound vectors. Return zero if any vector element lies outside the range,
 * in which case *failIx may optionally receive the index of the first failed item.
 * Notes:
 *   1. All three vector arguments must have identical data format. The number of
 *      elements in vec must not exceed the length of lowerVec and upperVec.
 *   2. For complex data types the real and imaginary components are verified 
 *      independently against real and imaginary parts of lower and upper range
 *      bounds.
 *   3. For floating-point data formats: if either of lower and upper range bounds
 *      is NaN, then the test data item must be also NaN to pass the check. */
int vecCheckRange( const tVec * vec, const tVec * lowerVec, const tVec * upperVec, int * failIx );

/*------------------------------------------------------------------------------*
 *                        SEQ-file reader functions                             *
 *------------------------------------------------------------------------------*/

typedef struct
{
    char * filename;
    FILE * file;
}
tagSeqFile;

typedef tagSeqFile* tSeqFile;

/* Open a SEQ-file and return either the file descriptor or zero if failed. */
tSeqFile seqFileOpen( const char * fname );
/* Close a SEQ-file. */
void seqFileClose( tSeqFile seqFile );
/* fscanf-equivalent reading function for a SEQ-file. Returns the number of data
 * fields that are successfully assigned. */
int seqFileScanf( tSeqFile seqFile, const char * fmt, ... );
/* Read a data vector. Argument vec specifies vector length and data format. It
 * must also have a memory block allocated to accomodate the loaded data. The
 * function returns zero if failed to read vector data. */
int seqFileReadVec( tSeqFile seqFile, const tVec * vec );
/* Read a few data vectors. Repeat the last argument for each new vector. The
 * arguments list must end up with zero pointer. The function returns zero if
 * failed to read some or all of the specified vectors. */
int seqFileReadVecs( tSeqFile seqFile, const tVec * vec, ... );

#ifdef __cplusplus
};
#endif

#endif /* __VECTOOLS_H */
