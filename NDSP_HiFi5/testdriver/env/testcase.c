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
 * Test case textual representation
 */

#include "types.h"
#include "testcase.h"

/* Test case semantic types. */
const char * testCaseTypeStr[MAXTESTCASE+1] = {
  "ACCURACY"                    , /*   0 */
  "KNOWN DATA"                  , /*   1 */
  "BAD SIZE/PARAM"              , /*   2 */
  "SOME AMOUNT OF OUTLIERS"     , /*   3 */
  "EDGE CASES"                  , /*   4 */
  "SPECIAL NUMBERS"             , /*   5 */
  "FORWARD FFT TEST"            , /*   6 */
  "INVERSE FFT TEST"            , /*   7 */
  "FFT RECONSTRUCTION TEST"     , /*   8 */
  "2D FFT ROW-WISE OVERFLOW"    , /*   9 */
  "2D FFT COLUMN-WISE OVERFLOW" , /*  10 */
  "2D FFT WGN"                  , /*  11 */
  "2D FFT SPECTRAL SPOTS"       , /*  12 */
  "SMALL SAMPLE BLOCK SIZE"     , /*  13 */
  "SWEPT SINEWAVE TEST"         , /*  14 */
  "NOISE TEST"                  , /*  15 */
  "SPREAD SPECTRUM TEST"        , /*  16 */
    /* Cadence's specific test patterns */
  "INP_DATA_PAT_MPI_TO_PPI_F"            , /*  17 */
  "INP_DATA_PAT_M1_TO_P1_F"              , /*  18 */
  "INP_DATA_PAT_M10_TO_P10_F"            , /*  19 */
  "INP_DATA_PAT_M100_TO_P100_F"          , /*  20 */
  "INP_DATA_PAT_M1000_TO_P1000_F"        , /*  21 */
  "INP_DATA_PAT_M10000_TO_P10000_F"      , /*  22 */
  "INP_DATA_PAT_M100000_TO_P100000_F"    , /*  23 */
  "INP_DATA_PAT_RAND_F"                  , /*  24 */
  "INP_DATA_PAT_EXP_LIM_F"               , /*  25 */
  "INP_DATA_PAT_MANT_LIM_F"              , /*  26 */
  "INP_DATA_PAT_SUBNORMAL_F"             , /*  27 */
  "INP_DATA_PAT_SPECIAL_F"               , /*  28 */
  "INP_DATA_PAT_0_TO_P100000_F"          , /*  29 */
  "INP_DATA_PAT_M100000_TO_0_F"          , /*  30 */

    "INP_DATA_PAT_MPI_TO_PPI_F_MPI_TO_PPI_F" ,/* 31 */
    "INP_DATA_PAT_MPI_TO_PPI_F_M1_TO_P1_F" ,/* 32 */
    "INP_DATA_PAT_MPI_TO_PPI_F_M10_TO_P10_F" ,/* 33 */
    "INP_DATA_PAT_MPI_TO_PPI_F_M100_TO_P100_F" ,/* 34 */
    "INP_DATA_PAT_MPI_TO_PPI_F_M1000_TO_P1000_F" ,/* 35 */
    "INP_DATA_PAT_MPI_TO_PPI_F_M10000_TO_P10000_F" ,/* 36 */
    "INP_DATA_PAT_MPI_TO_PPI_F_M100000_TO_P100000_F" ,/* 37 */
    "INP_DATA_PAT_MPI_TO_PPI_F_RAND_F" ,/* 38 */
    "INP_DATA_PAT_MPI_TO_PPI_F_EXP_LIM_F" ,/* 39 */
    "INP_DATA_PAT_MPI_TO_PPI_F_MANT_LIM_F" ,/* 40 */
    "INP_DATA_PAT_MPI_TO_PPI_F_SUBNORMAL_F" ,/* 41 */
    "INP_DATA_PAT_MPI_TO_PPI_F_SPECIAL_F" ,/* 42 */
    "INP_DATA_PAT_MPI_TO_PPI_F_0_TO_P100000_F" ,/* 43 */
    "INP_DATA_PAT_MPI_TO_PPI_F_M100000_TO_0_F" ,/* 44 */
    "INP_DATA_PAT_M1_TO_P1_F_MPI_TO_PPI_F" ,/* 45 */
    "INP_DATA_PAT_M1_TO_P1_F_M1_TO_P1_F" ,/* 46 */
    "INP_DATA_PAT_M1_TO_P1_F_M10_TO_P10_F" ,/* 47 */
    "INP_DATA_PAT_M1_TO_P1_F_M100_TO_P100_F" ,/* 48 */
    "INP_DATA_PAT_M1_TO_P1_F_M1000_TO_P1000_F" ,/* 49 */
    "INP_DATA_PAT_M1_TO_P1_F_M10000_TO_P10000_F" ,/* 50 */
    "INP_DATA_PAT_M1_TO_P1_F_M100000_TO_P100000_F" ,/* 51 */
    "INP_DATA_PAT_M1_TO_P1_F_RAND_F" ,/* 52 */
    "INP_DATA_PAT_M1_TO_P1_F_EXP_LIM_F" ,/* 53 */
    "INP_DATA_PAT_M1_TO_P1_F_MANT_LIM_F" ,/* 54 */
    "INP_DATA_PAT_M1_TO_P1_F_SUBNORMAL_F" ,/* 55 */
    "INP_DATA_PAT_M1_TO_P1_F_SPECIAL_F" ,/* 56 */
    "INP_DATA_PAT_M1_TO_P1_F_0_TO_P100000_F" ,/* 57 */
    "INP_DATA_PAT_M1_TO_P1_F_M100000_TO_0_F" ,/* 58 */
    "INP_DATA_PAT_M10_TO_P10_F_MPI_TO_PPI_F" ,/* 59 */
    "INP_DATA_PAT_M10_TO_P10_F_M1_TO_P1_F" ,/* 60 */
    "INP_DATA_PAT_M10_TO_P10_F_M10_TO_P10_F" ,/* 61 */
    "INP_DATA_PAT_M10_TO_P10_F_M100_TO_P100_F" ,/* 62 */
    "INP_DATA_PAT_M10_TO_P10_F_M1000_TO_P1000_F" ,/* 63 */
    "INP_DATA_PAT_M10_TO_P10_F_M10000_TO_P10000_F" ,/* 64 */
    "INP_DATA_PAT_M10_TO_P10_F_M100000_TO_P100000_F" ,/* 65 */
    "INP_DATA_PAT_M10_TO_P10_F_RAND_F" ,/* 66 */
    "INP_DATA_PAT_M10_TO_P10_F_EXP_LIM_F" ,/* 67 */
    "INP_DATA_PAT_M10_TO_P10_F_MANT_LIM_F" ,/* 68 */
    "INP_DATA_PAT_M10_TO_P10_F_SUBNORMAL_F" ,/* 69 */
    "INP_DATA_PAT_M10_TO_P10_F_SPECIAL_F" ,/* 70 */
    "INP_DATA_PAT_M10_TO_P10_F_0_TO_P100000_F" ,/* 71 */
    "INP_DATA_PAT_M10_TO_P10_F_M100000_TO_0_F" ,/* 72 */
    "INP_DATA_PAT_M100_TO_P100_F_MPI_TO_PPI_F" ,/* 73 */
    "INP_DATA_PAT_M100_TO_P100_F_M1_TO_P1_F" ,/* 74 */
    "INP_DATA_PAT_M100_TO_P100_F_M10_TO_P10_F" ,/* 75 */
    "INP_DATA_PAT_M100_TO_P100_F_M100_TO_P100_F" ,/* 76 */
    "INP_DATA_PAT_M100_TO_P100_F_M1000_TO_P1000_F" ,/* 77 */
    "INP_DATA_PAT_M100_TO_P100_F_M10000_TO_P10000_F" ,/* 78 */
    "INP_DATA_PAT_M100_TO_P100_F_M100000_TO_P100000_F" ,/* 79 */
    "INP_DATA_PAT_M100_TO_P100_F_RAND_F" ,/* 80 */
    "INP_DATA_PAT_M100_TO_P100_F_EXP_LIM_F" ,/* 81 */
    "INP_DATA_PAT_M100_TO_P100_F_MANT_LIM_F" ,/* 82 */
    "INP_DATA_PAT_M100_TO_P100_F_SUBNORMAL_F" ,/* 83 */
    "INP_DATA_PAT_M100_TO_P100_F_SPECIAL_F" ,/* 84 */
    "INP_DATA_PAT_M100_TO_P100_F_0_TO_P100000_F" ,/* 85 */
    "INP_DATA_PAT_M100_TO_P100_F_M100000_TO_0_F" ,/* 86 */
    "INP_DATA_PAT_M1000_TO_P1000_F_MPI_TO_PPI_F" ,/* 87 */
    "INP_DATA_PAT_M1000_TO_P1000_F_M1_TO_P1_F" ,/* 88 */
    "INP_DATA_PAT_M1000_TO_P1000_F_M10_TO_P10_F" ,/* 89 */
    "INP_DATA_PAT_M1000_TO_P1000_F_M100_TO_P100_F" ,/* 90 */
    "INP_DATA_PAT_M1000_TO_P1000_F_M1000_TO_P1000_F" ,/* 91 */
    "INP_DATA_PAT_M1000_TO_P1000_F_M10000_TO_P10000_F" ,/* 92 */
    "INP_DATA_PAT_M1000_TO_P1000_F_M100000_TO_P100000_F" ,/* 93 */
    "INP_DATA_PAT_M1000_TO_P1000_F_RAND_F" ,/* 94 */
    "INP_DATA_PAT_M1000_TO_P1000_F_EXP_LIM_F" ,/* 95 */
    "INP_DATA_PAT_M1000_TO_P1000_F_MANT_LIM_F" ,/* 96 */
    "INP_DATA_PAT_M1000_TO_P1000_F_SUBNORMAL_F" ,/* 97 */
    "INP_DATA_PAT_M1000_TO_P1000_F_SPECIAL_F" ,/* 98 */
    "INP_DATA_PAT_M1000_TO_P1000_F_0_TO_P100000_F" ,/* 99 */
    "INP_DATA_PAT_M1000_TO_P1000_F_M100000_TO_0_F" ,/* 100 */
    "INP_DATA_PAT_M10000_TO_P10000_F_MPI_TO_PPI_F" ,/* 101 */
    "INP_DATA_PAT_M10000_TO_P10000_F_M1_TO_P1_F" ,/* 102 */
    "INP_DATA_PAT_M10000_TO_P10000_F_M10_TO_P10_F" ,/* 103 */
    "INP_DATA_PAT_M10000_TO_P10000_F_M100_TO_P100_F" ,/* 104 */
    "INP_DATA_PAT_M10000_TO_P10000_F_M1000_TO_P1000_F" ,/* 105 */
    "INP_DATA_PAT_M10000_TO_P10000_F_M10000_TO_P10000_F" ,/* 106 */
    "INP_DATA_PAT_M10000_TO_P10000_F_M100000_TO_P100000_F" ,/* 107 */
    "INP_DATA_PAT_M10000_TO_P10000_F_RAND_F" ,/* 108 */
    "INP_DATA_PAT_M10000_TO_P10000_F_EXP_LIM_F" ,/* 109 */
    "INP_DATA_PAT_M10000_TO_P10000_F_MANT_LIM_F" ,/* 110 */
    "INP_DATA_PAT_M10000_TO_P10000_F_SUBNORMAL_F" ,/* 111 */
    "INP_DATA_PAT_M10000_TO_P10000_F_SPECIAL_F" ,/* 112 */
    "INP_DATA_PAT_M10000_TO_P10000_F_0_TO_P100000_F" ,/* 113 */
    "INP_DATA_PAT_M10000_TO_P10000_F_M100000_TO_0_F" ,/* 114 */
    "INP_DATA_PAT_M100000_TO_P100000_F_MPI_TO_PPI_F" ,/* 115 */
    "INP_DATA_PAT_M100000_TO_P100000_F_M1_TO_P1_F" ,/* 116 */
    "INP_DATA_PAT_M100000_TO_P100000_F_M10_TO_P10_F" ,/* 117 */
    "INP_DATA_PAT_M100000_TO_P100000_F_M100_TO_P100_F" ,/* 118 */
    "INP_DATA_PAT_M100000_TO_P100000_F_M1000_TO_P1000_F" ,/* 119 */
    "INP_DATA_PAT_M100000_TO_P100000_F_M10000_TO_P10000_F" ,/* 120 */
    "INP_DATA_PAT_M100000_TO_P100000_F_M100000_TO_P100000_F" ,/* 121 */
    "INP_DATA_PAT_M100000_TO_P100000_F_RAND_F" ,/* 122 */
    "INP_DATA_PAT_M100000_TO_P100000_F_EXP_LIM_F" ,/* 123 */
    "INP_DATA_PAT_M100000_TO_P100000_F_MANT_LIM_F" ,/* 124 */
    "INP_DATA_PAT_M100000_TO_P100000_F_SUBNORMAL_F" ,/* 125 */
    "INP_DATA_PAT_M100000_TO_P100000_F_SPECIAL_F" ,/* 126 */
    "INP_DATA_PAT_M100000_TO_P100000_F_0_TO_P100000_F" ,/* 127 */
    "INP_DATA_PAT_M100000_TO_P100000_F_M100000_TO_0_F" ,/* 128 */
    "INP_DATA_PAT_RAND_F_MPI_TO_PPI_F" ,/* 129 */
    "INP_DATA_PAT_RAND_F_M1_TO_P1_F" ,/* 130 */
    "INP_DATA_PAT_RAND_F_M10_TO_P10_F" ,/* 131 */
    "INP_DATA_PAT_RAND_F_M100_TO_P100_F" ,/* 132 */
    "INP_DATA_PAT_RAND_F_M1000_TO_P1000_F" ,/* 133 */
    "INP_DATA_PAT_RAND_F_M10000_TO_P10000_F" ,/* 134 */
    "INP_DATA_PAT_RAND_F_M100000_TO_P100000_F" ,/* 135 */
    "INP_DATA_PAT_RAND_F_RAND_F" ,/* 136 */
    "INP_DATA_PAT_RAND_F_EXP_LIM_F" ,/* 137 */
    "INP_DATA_PAT_RAND_F_MANT_LIM_F" ,/* 138 */
    "INP_DATA_PAT_RAND_F_SUBNORMAL_F" ,/* 139 */
    "INP_DATA_PAT_RAND_F_SPECIAL_F" ,/* 140 */
    "INP_DATA_PAT_RAND_F_0_TO_P100000_F" ,/* 141 */
    "INP_DATA_PAT_RAND_F_M100000_TO_0_F" ,/* 142 */
    "INP_DATA_PAT_EXP_LIM_F_MPI_TO_PPI_F" ,/* 143 */
    "INP_DATA_PAT_EXP_LIM_F_M1_TO_P1_F" ,/* 144 */
    "INP_DATA_PAT_EXP_LIM_F_M10_TO_P10_F" ,/* 145 */
    "INP_DATA_PAT_EXP_LIM_F_M100_TO_P100_F" ,/* 146 */
    "INP_DATA_PAT_EXP_LIM_F_M1000_TO_P1000_F" ,/* 147 */
    "INP_DATA_PAT_EXP_LIM_F_M10000_TO_P10000_F" ,/* 148 */
    "INP_DATA_PAT_EXP_LIM_F_M100000_TO_P100000_F" ,/* 149 */
    "INP_DATA_PAT_EXP_LIM_F_RAND_F" ,/* 150 */
    "INP_DATA_PAT_EXP_LIM_F_EXP_LIM_F" ,/* 151 */
    "INP_DATA_PAT_EXP_LIM_F_MANT_LIM_F" ,/* 152 */
    "INP_DATA_PAT_EXP_LIM_F_SUBNORMAL_F" ,/* 153 */
    "INP_DATA_PAT_EXP_LIM_F_SPECIAL_F" ,/* 154 */
    "INP_DATA_PAT_EXP_LIM_F_0_TO_P100000_F" ,/* 155 */
    "INP_DATA_PAT_EXP_LIM_F_M100000_TO_0_F" ,/* 156 */
    "INP_DATA_PAT_MANT_LIM_F_MPI_TO_PPI_F" ,/* 157 */
    "INP_DATA_PAT_MANT_LIM_F_M1_TO_P1_F" ,/* 158 */
    "INP_DATA_PAT_MANT_LIM_F_M10_TO_P10_F" ,/* 159 */
    "INP_DATA_PAT_MANT_LIM_F_M100_TO_P100_F" ,/* 160 */
    "INP_DATA_PAT_MANT_LIM_F_M1000_TO_P1000_F" ,/* 161 */
    "INP_DATA_PAT_MANT_LIM_F_M10000_TO_P10000_F" ,/* 162 */
    "INP_DATA_PAT_MANT_LIM_F_M100000_TO_P100000_F" ,/* 163 */
    "INP_DATA_PAT_MANT_LIM_F_RAND_F" ,/* 164 */
    "INP_DATA_PAT_MANT_LIM_F_EXP_LIM_F" ,/* 165 */
    "INP_DATA_PAT_MANT_LIM_F_MANT_LIM_F" ,/* 166 */
    "INP_DATA_PAT_MANT_LIM_F_SUBNORMAL_F" ,/* 167 */
    "INP_DATA_PAT_MANT_LIM_F_SPECIAL_F" ,/* 168 */
    "INP_DATA_PAT_MANT_LIM_F_0_TO_P100000_F" ,/* 169 */
    "INP_DATA_PAT_MANT_LIM_F_M100000_TO_0_F" ,/* 170 */
    "INP_DATA_PAT_SUBNORMAL_F_MPI_TO_PPI_F" ,/* 171 */
    "INP_DATA_PAT_SUBNORMAL_F_M1_TO_P1_F" ,/* 172 */
    "INP_DATA_PAT_SUBNORMAL_F_M10_TO_P10_F" ,/* 173 */
    "INP_DATA_PAT_SUBNORMAL_F_M100_TO_P100_F" ,/* 174 */
    "INP_DATA_PAT_SUBNORMAL_F_M1000_TO_P1000_F" ,/* 175 */
    "INP_DATA_PAT_SUBNORMAL_F_M10000_TO_P10000_F" ,/* 176 */
    "INP_DATA_PAT_SUBNORMAL_F_M100000_TO_P100000_F" ,/* 177 */
    "INP_DATA_PAT_SUBNORMAL_F_RAND_F" ,/* 178 */
    "INP_DATA_PAT_SUBNORMAL_F_EXP_LIM_F" ,/* 179 */
    "INP_DATA_PAT_SUBNORMAL_F_MANT_LIM_F" ,/* 180 */
    "INP_DATA_PAT_SUBNORMAL_F_SUBNORMAL_F" ,/* 181 */
    "INP_DATA_PAT_SUBNORMAL_F_SPECIAL_F" ,/* 182 */
    "INP_DATA_PAT_SUBNORMAL_F_0_TO_P100000_F" ,/* 183 */
    "INP_DATA_PAT_SUBNORMAL_F_M100000_TO_0_F" ,/* 184 */
    "INP_DATA_PAT_SPECIAL_F_MPI_TO_PPI_F" ,/* 185 */
    "INP_DATA_PAT_SPECIAL_F_M1_TO_P1_F" ,/* 186 */
    "INP_DATA_PAT_SPECIAL_F_M10_TO_P10_F" ,/* 187 */
    "INP_DATA_PAT_SPECIAL_F_M100_TO_P100_F" ,/* 188 */
    "INP_DATA_PAT_SPECIAL_F_M1000_TO_P1000_F" ,/* 189 */
    "INP_DATA_PAT_SPECIAL_F_M10000_TO_P10000_F" ,/* 190 */
    "INP_DATA_PAT_SPECIAL_F_M100000_TO_P100000_F" ,/* 191 */
    "INP_DATA_PAT_SPECIAL_F_RAND_F" ,/* 192 */
    "INP_DATA_PAT_SPECIAL_F_EXP_LIM_F" ,/* 193 */
    "INP_DATA_PAT_SPECIAL_F_MANT_LIM_F" ,/* 194 */
    "INP_DATA_PAT_SPECIAL_F_SUBNORMAL_F" ,/* 195 */
    "INP_DATA_PAT_SPECIAL_F_SPECIAL_F" ,/* 196 */
    "INP_DATA_PAT_SPECIAL_F_0_TO_P100000_F" ,/* 197 */
    "INP_DATA_PAT_SPECIAL_F_M100000_TO_0_F" ,/* 198 */
    "INP_DATA_PAT_0_TO_P100000_F_MPI_TO_PPI_F" ,/* 199 */
    "INP_DATA_PAT_0_TO_P100000_F_M1_TO_P1_F" ,/* 200 */
    "INP_DATA_PAT_0_TO_P100000_F_M10_TO_P10_F" ,/* 201 */
    "INP_DATA_PAT_0_TO_P100000_F_M100_TO_P100_F" ,/* 202 */
    "INP_DATA_PAT_0_TO_P100000_F_M1000_TO_P1000_F" ,/* 203 */
    "INP_DATA_PAT_0_TO_P100000_F_M10000_TO_P10000_F" ,/* 204 */
    "INP_DATA_PAT_0_TO_P100000_F_M100000_TO_P100000_F" ,/* 205 */
    "INP_DATA_PAT_0_TO_P100000_F_RAND_F" ,/* 206 */
    "INP_DATA_PAT_0_TO_P100000_F_EXP_LIM_F" ,/* 207 */
    "INP_DATA_PAT_0_TO_P100000_F_MANT_LIM_F" ,/* 208 */
    "INP_DATA_PAT_0_TO_P100000_F_SUBNORMAL_F" ,/* 209 */
    "INP_DATA_PAT_0_TO_P100000_F_SPECIAL_F" ,/* 210 */
    "INP_DATA_PAT_0_TO_P100000_F_0_TO_P100000_F" ,/* 211 */
    "INP_DATA_PAT_0_TO_P100000_F_M100000_TO_0_F" ,/* 212 */
    "INP_DATA_PAT_M100000_TO_0_F_MPI_TO_PPI_F" ,/* 213 */
    "INP_DATA_PAT_M100000_TO_0_F_M1_TO_P1_F" ,/* 214 */
    "INP_DATA_PAT_M100000_TO_0_F_M10_TO_P10_F" ,/* 215 */
    "INP_DATA_PAT_M100000_TO_0_F_M100_TO_P100_F" ,/* 216 */
    "INP_DATA_PAT_M100000_TO_0_F_M1000_TO_P1000_F" ,/* 217 */
    "INP_DATA_PAT_M100000_TO_0_F_M10000_TO_P10000_F" ,/* 218 */
    "INP_DATA_PAT_M100000_TO_0_F_M100000_TO_P100000_F" ,/* 219 */
    "INP_DATA_PAT_M100000_TO_0_F_RAND_F" ,/* 220 */
    "INP_DATA_PAT_M100000_TO_0_F_EXP_LIM_F" ,/* 221 */
    "INP_DATA_PAT_M100000_TO_0_F_MANT_LIM_F" ,/* 222 */
    "INP_DATA_PAT_M100000_TO_0_F_SUBNORMAL_F" ,/* 223 */
    "INP_DATA_PAT_M100000_TO_0_F_SPECIAL_F" ,/* 224 */
    "INP_DATA_PAT_M100000_TO_0_F_0_TO_P100000_F" ,/* 225 */
    "INP_DATA_PAT_M100000_TO_0_F_M100000_TO_0_F" ,/* 226 */

  "UNKNOWN"                     , /* >226 */
};
