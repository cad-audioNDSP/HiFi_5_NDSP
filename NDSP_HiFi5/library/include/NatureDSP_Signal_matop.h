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
#ifndef __NATUREDSP_SIGNAL_MATOP_H__
#define __NATUREDSP_SIGNAL_MATOP_H__

#include "NatureDSP_types.h"

#ifdef __cplusplus
extern "C" {
#endif

/*===========================================================================
  Matrix Operations:
  mtx_mpy              Matrix Multiply
  [c]mtx_vecmpy        Matrix by Vector Multiply
  [c]mtx_vecmpyt       Vector by Vector Multiply
  [c]mtx_lrmpy         Three Matrices Product
  mtx_transpose        Matrix Transpose
===========================================================================*/

/*-------------------------------------------------------------------------
  Matrix Multiply
  These functions compute the expression z = 2^lsh * x * y for the matrices 
  x and y. The columnar dimension of x must match the row dimension of y. 
  The resulting matrix has the same number of rows as x and the same number 
  of columns as y.
  Transposing API allows to interpret input yt as transposed matrix y.

  NOTE: lsh factor is not relevant for floating point routines.

  Functions require scratch memory for storing intermediate data. This 
  scratch memory area should be aligned on 16 byte boundary and its size is 
  calculated by dedicated scratch allocation functions.

  Two versions of functions available: regular version (mtx_mpy[t]32x32, 
  mtx_mpy[t]16x16, mtx_mpy[t]8x16, mtx_mpy[t]8x8, mtx[t]_mpyf) with 
  arbitrary arguments and faster version (mtx_mpy[t]32x32_fast, 
  mtx_mpy[t]16x16_fast, mtx_mpy[t]8x16_fast, mtx_mpy[t]8x8_fast, 
  mtx_mpy[t]f_fast, cntx_mpyt32x32_fast, cntx_mpytf_fast) that apply 
  some restrictions

  Precision:
  32x32 32-bit inputs, 32-bit output
  16x16 16-bit inputs, 16-bit output
  8x8   8-bit inputs, 8-bit output
  8x16  8/16-bit inputs, 16-bit output
  f     floating point

  Input:
  x[M*N]      input matrix x, Q7, Q15, Q31 or floating point
  y[N*P]      input matrix y, Q7, Q15, Q31 or floating point
  yt[P*N]     transposed input matrix y. Q31,Q15, Q7 floating point. (for 
              transposing API only)
  M           number of rows in matrix x and z
  N           number of columns in matrix x and number of rows in matrix y
  P           number of columns in matrices y and z
  lsh         left shift applied to the result (applied to the fixed-
              point functions only) 
  Output:
  z[M*P]      output matrix z, Q7, Q15, Q31 or floating point 
  Scratch:
  pScr        size in bytes defined by corresponding scratch allocation 
              functions

  Restrictions:
  For regular routines mpy[t]32x32, mtx_mpy[t]16x16, mtx_mpy[t]8x8, 
  mtx_mpy[t]8x16, mtx_mpy[t]f):
  pScr    aligned on 16-byte boundary
  x,y,z   should not overlap

  For faster routines (mtx_mpy[t]32x32_fast, mtx_mpy[t]16x16_fast, 
  mtx_mpy[t]8x8_fast, mtx_mpy[t]8x16_fast, 
  mtx_mpy[t]f_fast):
  x,y,z       should not overlap
  x,y,z,pScr  aligned on 16-byte boundary
  M,N,P       multiplies of 4 for mtx_mpy[t]32x32_fast, mtx_mpy[t]16x16_fast, 
              mtx_mpy[t]8x8_fast, mtx_mpy[t]8x16_fast, mtx_mpy[t]f_fast
              multiplies of 32 for cntx_mpyt32x32_fast, cntx_mpytf_fast
  lsh         should be in range:
              -31...31 for mtx_mpy32x32, mtx_mpy32x32_fast, cntx_mpyt32x32_fast, 
                       cntx_mpytf_fast
              -15...15 for mtx_mpy16x16, mtx_mpy16x16_fast, mtx_mpy[t]8x8, 
                       mtx_mpy[t]8x8_fast, mtx_mpy[t]8x16, 
                       mtx_mpy[t]8x16_fast 

-------------------------------------------------------------------------*/
void mtx_mpy32x32 ( void* pScr,
                    int32_t* z,
              const int32_t* x,
              const int32_t* y,
              int M, int N, int P, int lsh );
void mtx_mpy16x16 ( void* pScr,
                    int16_t* z,
              const int16_t* x,
              const int16_t* y,
              int M, int N, int P, int lsh );
void mtx_mpy8x8  (  void* pScr,
                     int8_t*  z, const int8_t*  x, const int8_t*  y,
                     int M, int N, int P, int lsh );
void mtx_mpy8x16 (  void* pScr,
                     int16_t*  z, const int8_t*  x, const int16_t*  y,
                     int M, int N, int P, int lsh );
void mtx_mpy32x32_fast ( void* pScr,int32_t* z,
                   const int32_t* x,
                   const int32_t* y,
                   int M, int N, int P, int lsh );
void mtx_mpy16x16_fast ( void* pScr, int16_t* z,
                   const int16_t* x,
                   const int16_t* y,
                   int M, int N, int P, int lsh );
void mtx_mpy8x8_fast(  void* pScr,
                     int8_t*  z, const int8_t*  x, const int8_t*  y,
                     int M, int N, int P, int lsh );
void mtx_mpy8x16_fast (  void* pScr,
                     int16_t*  z, const int8_t*  x, const int16_t*  y,
                     int M, int N, int P, int lsh );

void mtx_mpyf ( void* pScr, float32_t* z,
          const float32_t* x,
          const float32_t* y,
          int M, int N, int P);
void mtx_mpyf_fast ( void* pScr, float32_t* z,
               const float32_t* x,
               const float32_t* y,
               int M, int N, int P);

// Transposing API:
void mtx_mpyt16x16 ( void* pScr,
                     int16_t*  z, const int16_t*  x, const int16_t*  yt,
                     int M, int N, int P, int lsh );
void mtx_mpyt8x8 ( void* pScr,
                     int8_t*  z, const int8_t*  x, const int8_t*  yt,
                     int M, int N, int P, int lsh );
void mtx_mpyt8x16( void* pScr,
                     int16_t*  z, const int8_t*  x, const int16_t*  yt,
                     int M, int N, int P, int lsh );
void mtx_mpyt32x32 ( void* pScr,
                     int32_t*  z, const int32_t*  x, const int32_t*  yt,
                     int M, int N, int P, int lsh );
void mtx_mpytf (     void* pScr,
                     float32_t* z, const float32_t* x, const float32_t* yt,
                     int M, int N, int P);
void mtx_mpyt16x16_fast (  void* pScr,
                     int16_t*  z, const int16_t*  x, const int16_t*  yt,
                     int M, int N, int P, int lsh );
void mtx_mpyt8x8_fast (  void* pScr,
                     int8_t*  z, const int8_t*  x, const int8_t*  yt,
                     int M, int N, int P, int lsh );
void mtx_mpyt8x16_fast (  void* pScr,
                     int16_t*  z, const int8_t*  x, const int16_t*  yt,
                     int M, int N, int P, int lsh );
void mtx_mpyt32x32_fast (  void* pScr,
                     int32_t*  z, const int32_t*  x, const int32_t*  yt,
                     int M, int N, int P, int lsh );
void mtx_mpytf_fast( void* pScr,
                     float32_t* z, const float32_t* x, const float32_t* yt,
                     int M, int N, int P);
void cmtx_mpyt32x32_fast ( void* pScr,
                     complex_fract32*  z, 
                     const complex_fract32*  x, 
                     const complex_fract32*  yt,
                     int M, int N, int P, int lsh );
void cmtx_mpytf_fast ( void* pScr,
                     complex_float*  z, 
                     const complex_float *  x, 
                     const complex_float *  yt,
                     int M, int N, int P );

// scratch allocation functions 
size_t mtx_mpy16x16_getScratchSize      (int M, int N, int P);
size_t mtx_mpy8x8_getScratchSize        (int M, int N, int P);
size_t mtx_mpy8x16_getScratchSize       (int M, int N, int P);
size_t mtx_mpy32x32_getScratchSize      (int M, int N, int P);
size_t mtx_mpyf_getScratchSize          (int M, int N, int P);
size_t mtx_mpy16x16_fast_getScratchSize (int M, int N, int P);
size_t mtx_mpy8x8_fast_getScratchSize   (int M, int N, int P);
size_t mtx_mpy8x16_fast_getScratchSize  (int M, int N, int P);
size_t mtx_mpy32x32_fast_getScratchSize (int M, int N, int P);
size_t mtx_mpyf_fast_getScratchSize     (int M, int N, int P);

size_t mtx_mpyt16x16_getScratchSize       (int M, int N, int P);
size_t mtx_mpyt8x8_getScratchSize         (int M, int N, int P);
size_t mtx_mpyt8x16_getScratchSize        (int M, int N, int P);
size_t mtx_mpyt32x32_getScratchSize       (int M, int N, int P);
size_t mtx_mpytf_getScratchSize           (int M, int N, int P);
size_t mtx_mpyt16x16_fast_getScratchSize  (int M, int N, int P);
size_t mtx_mpyt8x8_fast_getScratchSize    (int M, int N, int P);
size_t mtx_mpyt8x16_fast_getScratchSize   (int M, int N, int P);
size_t mtx_mpyt32x32_fast_getScratchSize  (int M, int N, int P);
size_t mtx_mpytf_fast_getScratchSize      (int M, int N, int P);
size_t cmtx_mpyt32x32_fast_getScratchSize (int M, int N, int P);
size_t cmtx_mpytf_fast_getScratchSize     (int M, int N, int P);

/*-------------------------------------------------------------------------
  Matrix by Vector Multiply
  These functions compute the expression z = 2^lsh * x * y for the matrices 
  x and vector y. 
  NOTE: lsh factor is not relevant for floating point routines.

  Two versions of functions available: regular version (mtx_vecmpy32x32,  
  mtx_vecmpy16x16, mtx_vecmpy8x8, mtx_vecmpy8x16, mtx_vecmpyf) with arbitrary 
  arguments and faster version (mtx_vecmpy32x32_fast, mtx_vecmpy16x16_fast, 
  mtx_vecmpy8x8_fast, mtx_vecmpy8x16_fast,  mtx_vecmpyf_fast) that apply 
  some restrictions

  Precision: 
  32x32 32-bit input, 32-bit output
  16x16 16-bit input, 16-bit output
  8x8   8-bit inputs, 8-bit output
  8x16  8/16-bit inputs, 16-bit output
  f     floating point

  Input:
  x[M*N] input matrix,Q31,Q15 or floating point
  y[N]   input vector,Q31,Q15 or floating point
  M      number of rows in matrix x
  N      number of columns in matrix x
  lsh    additional left shift(applied to the fixed-
         point functions only) 
  Output:
  z[M]   output vector,Q31,Q15 or floating point

  Restriction:
  For regular routines (mtx_vecmpy32x32, mtx_vecmpy16x16, mtx_vecmpy8x8,
  mtx_vecmpy8x16,  mtx_vecmpyf)
  x,y,z should not overlap

  For faster routines  (mtx_vecmpy32x32_fast, mtx_vecmpy16x16_fast, 
  mtx_vecmpy8x8_fast, mtx_vecmpy8x16_fast, mtx_vecmpyf_fast)
  x,y,z   should not overlap
  x,y     aligned on 16-byte boundary
  N, M    multiples of 4
  lsh     should be in range:
          -31...31 for mtx_vecmpy32x32, mtx_vecmpy32x32_fast
          -15...15 for mtx_vecmpy16x16, mtx_vecmpy16x16_fast, 
                   mtx_vecmpy8x8_fast, mtx_vecmpy8x16_fast   
-------------------------------------------------------------------------*/
void mtx_vecmpy32x32 ( int32_t* z,
                 const int32_t* x,
                 const int32_t* y,
                 int M, int N, int lsh);
void mtx_vecmpy16x16 ( int16_t* z,
                 const int16_t* x,
                 const int16_t* y,
                 int M, int N, int lsh);
void mtx_vecmpy8x8 ( int8_t*  z,
               const int8_t*  x,
               const int8_t*  y,
               int M, int N, int lsh);
void mtx_vecmpy8x16( int16_t*  z,
               const int8_t *  x,
               const int16_t*  y,
               int M, int N, int lsh);
void mtx_vecmpy32x32_fast ( int32_t* z,
                      const int32_t* x,
                      const int32_t* y,
                      int M, int N, int lsh);
void mtx_vecmpy16x16_fast ( int16_t* z,
                      const int16_t* x,
                      const int16_t* y,
                      int M, int N, int lsh);
void mtx_vecmpy8x8_fast ( int8_t*  z,
               const int8_t*  x,
               const int8_t*  y,
               int M, int N, int lsh);
void mtx_vecmpy8x16_fast ( int16_t*  z,
               const int8_t *  x,
               const int16_t*  y,
               int M, int N, int lsh);
void mtx_vecmpyf ( float32_t* z,
             const float32_t* x,
             const float32_t* y,
             int M, int N);
void mtx_vecmpyf_fast ( float32_t* z,
                  const float32_t* x,
                  const float32_t* y,
                  int M, int N);

/*-------------------------------------------------------------------------
  Vector by Vector Multiply 
  These functions compute the expression y = 2^lsh*x*xt for the input column
  vector x and its Hermitian transpose xt.

  NOTE: lsh factor is not relevant for floating point routines.

  Precision:
  32x32   32-bit input, 32-bit output
  f       floating point

  Input:
  x[N]    input vector, Q31 or floating point
  N       size of vector x
  lsh     bidirectional left shift applied to the result (fixed point 
          functions only).
  Output:
  y[N*N]  output matrix, Q31 or floating point

  Restrictions:
  x,y     should not overlap
  x,y     aligned on 16-byte boundary
  N       multiple of 32
  lsh     -31...31
-------------------------------------------------------------------------*/
void cmtx_vecmpyt32x32_fast ( complex_fract32* y,
                        const complex_fract32* x,
                        int N, int lsh );
void cmtx_vecmpytf_fast (     complex_float* y,
                        const complex_float* x,
                        int N );

/*-------------------------------------------------------------------------
  Three Matrices Product
  These functions compute the expression z = 2^lsh*((2^-rsh*x*y)*xt) for 
  square input matrices x and y, where xt is the Hermitian transpose of the
  matrix x.

  NOTE: 2^lsh and 2^-rsh factors are not relevant for floating point routines.

  Functions require scratch memory for storing intermediate data. This scratch
  memory area should be aligned on 16-byte boundary; its size is calculated
  by the corresponding scratch allocation function.

  Precision:
  32x32   32-bit input, 32-bit output
  f       floating point

  Input:
  x[N*N]  input matrix x, Q31 or floating point
  y[N*N]  input matrix y, Q31 or floating point
  N       number of rows and columns in matrices x, y and z
  rsh     right shift applied to intermediate results to avoid overflow (fixed 
          point functions only)
  lsh     bidirectional left shift applied to the result (fixed point functions
          only)
  Output:
  z[N*N]  output matrix z, Q31 or floating point
  Scratch:
  pScr    scratch memory area with size in bytes defined by the corresponding
          scratch allocation function

  Restrictions:
  x,y,z   should not overlap
  x,y,z   aligned on 16-byte boundary
  N       multiple of 32
  rsh     >=0
  lsh     -31...31
-------------------------------------------------------------------------*/
void cmtx_lrmpy32x32_fast ( void* pScr,
                        complex_fract32* z,
                  const complex_fract32* x, 
                  const complex_fract32* y,
                  int N, int rsh, int lsh );
void cmtx_lrmpyf_fast ( void* pScr,
                        complex_float* z,
                  const complex_float* x, 
                  const complex_float* y,
                  int N );

// scratch allocation functions 
size_t cmtx_lrmpy32x32_fast_getScratchSize(int N);
size_t cmtx_lrmpyf_fast_getScratchSize    (int N);

/*-------------------------------------------------------------------------
  Matrix Transpose
  These functions transpose matrices.

  Precision: 
  32x32 32-bit input, 32-bit output
  16x16 16-bit input, 16-bit output
  8x8   8-bit inputs, 8-bit output
  f     floating point

  Input:
  x[M][N] input matrix,Q31,Q15,Q7 or floating point
  M       number of rows in matrix x
  N       number of columns in matrix x
  Output:
  y[N][M] output vector,Q31,Q15,Q7 or floating point

  Restriction:
  For regular routines (mtx_transpose_32x32, mtx_transpose_16x16, 
  mtx_transpose_8x8, mtx_transposef):
  x,y should not overlap

  For faster routines (mtx_transpose 32x32_fast, mtx_transpose 16x16_fast, 
  mtx_transpose_8x8_fast, mtx_transposef_fast)
  x,y   should not overlap
  x,y   aligned on 16-byte boundary
  N and M are multiples of 4
-------------------------------------------------------------------------*/
void mtx_transpose32x32 (int32_t*    y, const int32_t*     x, int M, int N);
void mtx_transpose16x16 (int16_t*    y, const int16_t*     x, int M, int N);
void mtx_transpose8x8   (int8_t *    y, const int8_t *     x, int M, int N);
void mtx_transposef     (float32_t*  y, const float32_t *  x, int M, int N);

void mtx_transpose32x32_fast (int32_t  *  y, const int32_t*     x, int M, int N);
void mtx_transpose16x16_fast (int16_t  *  y, const int16_t*     x, int M, int N);
void mtx_transpose8x8_fast   (int8_t   *  y, const int8_t *     x, int M, int N);
void mtx_transposef_fast     (float32_t*  y, const float32_t *  x, int M, int N);

#ifdef __cplusplus
}
#endif

#endif/* __NATUREDSP_SIGNAL_MATOP_H__ */
