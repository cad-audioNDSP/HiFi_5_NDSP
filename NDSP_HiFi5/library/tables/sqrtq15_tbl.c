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
  NatureDSP Signal Processing Library. Matrix Operations
    Square Root tables
  IntegrIT, 2006-2019
*/

#include "NatureDSP_types.h"
#include "common.h"
#include "sqrtq15_tbl.h"

/*
    Polynomial coefficients for sqrt()
    taken by following matlab code:
    function calcsqrt16x16_coef()
    P=3;
    N=8;
    X=[];Y=[]; Z=[];
    for n=1:N
        x=1*((n-1)/N+(0:pow2(1,-16):1/N));
        z=sqrt(x);
        x0=1*(n-1)/N;
        p=polyfit(x-x0,z,P);
        y=polyval(p,x-x0);
        X=[X x];
        Y=[Y y];
        Z=[Z z];
        q=zeros(1,P+1); q=q+15;
        q(P+1)=15;
        for k=P:-1:1
            q(k)=q(k+1);
        end
        q(1)=14;
        fprintf(1,'%d ',sat16(round(pow2(p,q))));fprintf(1,'\n');
    end

    Finally, the order of columns is reversed to facilitate 
    using of HiFi5 SEL16X4 instruction
*/
const int16_t ALIGN(32) sqrtq15_tbl[] ={ 
     1206, 1726, 2623,  4341,  8169, 19131, 32767, 32767,
    -4970,-6249,-8185,-11369,-17296,-30877,-32768,-32768,
    17514,18917,20721, 23165, 26741, 32719, 32767, 32767,
    30652,28378,25905, 23171, 20066, 16384, 11587,  1470};
