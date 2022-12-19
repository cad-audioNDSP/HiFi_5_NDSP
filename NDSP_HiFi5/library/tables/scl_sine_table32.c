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
#include "scl_sine_table32.h"
#include "common.h"


/* -------------------------------------------------
	sine table: contains packed cos/sin 
	
	MATLAB CODE:
	phase=(0:512-1)/512*(2*pi);
	tbl=S_round_l(sin(phase)*32768.);
   	tbl = (max(min(tbl,32767),-32768));
	tbl=bitand((tbl+65536),65535);
	for k=1:16:length(tbl)
		fprintf(1,'0x%04x,0x%04x,0x%04x,0x%04x,0x%04x,0x%04x,0x%04x,0x%04x,0x%04x,0x%04x,0x%04x,0x%04x,0x%04x,0x%04x,0x%04x,0x%04x,\n',tbl(k:k+16));
	end


    ord=9;
    tbl=round(sin(2*pi*(0:2^ord-1)/(2^ord))*2^31);
    tbl=min(2^31-1,max(-2^31,tbl));
    for n=1:8:length(tbl)
        for k=0:7
            x=tbl(n+k);
            if (x<0) x=x+2^32; end
            fprintf(1,'0x%08x,',x);
        end
        fprintf(1,'\n');
    end

------------------------------------------------ */ 
#ifdef COMPILER_MSVC
const int32_t _declspec(align(32)) sine_table32[] =
#else
const int32_t                     sine_table32[] __attribute__((aligned(32))) =
#endif
{
  -1621691, 21637337, -160631214, 684490768, -1387189892, 843314609 };
