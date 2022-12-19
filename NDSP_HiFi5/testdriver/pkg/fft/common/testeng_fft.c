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
 * Test engine extension for FFT tests
 */

/* Test engine extension for FFT. */
#include "testeng_fft.h"

const void* GetFFT_handle(int N, int FFT_Id)
{
    int i, j;
    /* 32-bit handlers for real FFT */
    const  N_handle_pair_t  rfft32x32_hndls[NUM_RFFT32X32] =
    {
        { 32  , rfft32_32    }, { 64  , rfft32_64    }, { 128 , rfft32_128   }, { 256 , rfft32_256   },
        { 512 , rfft32_512   }, { 1024, rfft32_1024  }, { 2048, rfft32_2048  }, { 4096, rfft32_4096  },
        { 8192, rfft32_8192  }, { 12  , rnfft32_12   }, { 24  , rnfft32_24   }, { 36  , rnfft32_36   },
        { 48  , rnfft32_48   }, { 60  , rnfft32_60   }, { 72  , rnfft32_72   }, { 96  , rnfft32_96   },
        { 108 , rnfft32_108  }, { 120 , rnfft32_120  }, { 144 , rnfft32_144  }, { 180 , rnfft32_180  },
        { 192 , rnfft32_192  }, { 216 , rnfft32_216  }, { 240 , rnfft32_240  }, { 288 , rnfft32_288  },
        { 300 , rnfft32_300  }, { 324 , rnfft32_324  }, { 360 , rnfft32_360  }, { 432 , rnfft32_432  },
        { 480 , rnfft32_480  }, { 540 , rnfft32_540  }, { 576 , rnfft32_576  }, { 768 , rnfft32_768  },
        { 960 , rnfft32_960  }, { 30  , rnfft32_30   }, { 90  , rnfft32_90   }, { 384 , rnfft32_384  },
        { 720 , rnfft32_720  }, { 1152, rnfft32_1152 }, { 1440, rnfft32_1440 }, { 1536, rnfft32_1536 }, 
        { 1920, rnfft32_1920 }, { 160 , rnfft32_160  }, { 320 , rnfft32_320  }
    };
    const  N_handle_pair_t  rifft32x32_hndls[NUM_RFFT32X32] =
    {
        { 32  , rifft32_32    }, { 64  , rifft32_64    }, { 128 , rifft32_128   }, { 256 , rifft32_256   },
        { 512 , rifft32_512   }, { 1024, rifft32_1024  }, { 2048, rifft32_2048  }, { 4096, rifft32_4096  },
        { 8192, rifft32_8192  }, { 12  , rinfft32_12   }, { 24  , rinfft32_24   }, { 36  , rinfft32_36   },
        { 48  , rinfft32_48   }, { 60  , rinfft32_60   }, { 72  , rinfft32_72   }, { 96  , rinfft32_96   },
        { 108 , rinfft32_108  }, { 120 , rinfft32_120  }, { 144 , rinfft32_144  }, { 180 , rinfft32_180  },
        { 192 , rinfft32_192  }, { 216 , rinfft32_216  }, { 240 , rinfft32_240  }, { 288 , rinfft32_288  },
        { 300 , rinfft32_300  }, { 324 , rinfft32_324  }, { 360 , rinfft32_360  }, { 432 , rinfft32_432  },
        { 480 , rinfft32_480  }, { 540 , rinfft32_540  }, { 576 , rinfft32_576  }, { 768 , rinfft32_768  },
        { 960 , rinfft32_960  }, { 30  , rinfft32_30   }, { 90  , rinfft32_90   }, { 384 , rinfft32_384  },
        { 720 , rinfft32_720  }, { 1152, rinfft32_1152 }, { 1440, rinfft32_1440 }, { 1536, rinfft32_1536 },
        { 1920, rinfft32_1920 }, { 160 , rinfft32_160  }, { 320 , rinfft32_320  }
    };
    /* 32-bit handlers for complex FFT */
    const  N_handle_pair_t  cfft32x32_hndls[NUM_CFFT32X32] =
    {
        { 32 , cfft32_32   }, { 64  , cfft32_64   }, { 128 , cfft32_128  }, { 256 , cfft32_256  },
        { 512, cfft32_512  }, { 1024, cfft32_1024 }, { 2048, cfft32_2048 }, { 4096, cfft32_4096 },
        { 16 , cfft32_16   }, { 12  , cnfft32_12  }, { 24  , cnfft32_24  }, { 36  , cnfft32_36  },
        { 48 , cnfft32_48  }, { 60  , cnfft32_60  }, { 72  , cnfft32_72  }, { 96  , cnfft32_96  },
        { 108, cnfft32_108 }, { 120 , cnfft32_120 }, { 144 , cnfft32_144 }, { 180 , cnfft32_180 },
        { 192, cnfft32_192 }, { 216 , cnfft32_216 }, { 240 , cnfft32_240 }, { 288 , cnfft32_288 },
        { 300, cnfft32_300 }, { 324 , cnfft32_324 }, { 360 , cnfft32_360 }, { 432 , cnfft32_432 },
        { 480, cnfft32_480 }, { 540 , cnfft32_540 }, { 576 , cnfft32_576 }, { 768 , cnfft32_768 },
        { 960, cnfft32_960 }, { 600 , cnfft32_600 }, { 400 , cnfft32_400 }, { 384 , cnfft32_384 },
        { 200, cnfft32_200 }, { 160 , cnfft32_160 }, { 100 , cnfft32_100 }, { 80  , cnfft32_80  },
        { 320, cnfft32_320 }
    };
    const  N_handle_pair_t  cifft32x32_hndls[NUM_CFFT32X32] =
    {
        { 32 , cifft32_32   }, { 64  , cifft32_64   }, { 128 , cifft32_128  }, { 256 , cifft32_256  },
        { 512, cifft32_512  }, { 1024, cifft32_1024 }, { 2048, cifft32_2048 }, { 4096, cifft32_4096 },
        { 16 , cifft32_16   }, { 12  , cinfft32_12  }, { 24  , cinfft32_24  }, { 36  , cinfft32_36  },
        { 48 , cinfft32_48  }, { 60  , cinfft32_60  }, { 72  , cinfft32_72  }, { 96  , cinfft32_96  },
        { 108, cinfft32_108 }, { 120 , cinfft32_120 }, { 144 , cinfft32_144 }, { 180 , cinfft32_180 },
        { 192, cinfft32_192 }, { 216 , cinfft32_216 }, { 240 , cinfft32_240 }, { 288 , cinfft32_288 },
        { 300, cinfft32_300 }, { 324 , cinfft32_324 }, { 360 , cinfft32_360 }, { 432 , cinfft32_432 },
        { 480, cinfft32_480 }, { 540 , cinfft32_540 }, { 576 , cinfft32_576 }, { 768 , cinfft32_768 },
        { 960, cinfft32_960 }, { 600 , cinfft32_600 }, { 400 , cinfft32_400 }, { 384 , cinfft32_384 },
        { 200, cinfft32_200 }, { 160 , cinfft32_160 }, { 100 , cinfft32_100 }, { 80  , cinfft32_80  },
        { 320, cinfft32_320 }
    };
    /* 16-bit handlers for real FFT */
    const  N_handle_pair_t  rfft16x16_hndls[NUM_RFFT16X16] =
    {
        { 32  , rfft16_32    }, { 64  , rfft16_64    }, { 128 , rfft16_128   }, { 256 , rfft16_256   },
        { 512 , rfft16_512   }, { 1024, rfft16_1024  }, { 2048, rfft16_2048  }, { 4096, rfft16_4096  },
        { 8192, rfft16_8192  },
        { 160 , rnfft16_160  }, { 192 , rnfft16_192  }, { 240 , rnfft16_240  }, { 320 , rnfft16_320  },
        { 384 , rnfft16_384  }, { 480 , rnfft16_480  }
    };
    const  N_handle_pair_t  rifft16x16_hndls[NUM_RFFT16X16] =
    {
        { 32  , rifft16_32    }, { 64  , rifft16_64    }, { 128 , rifft16_128   }, { 256 , rifft16_256   },
        { 512 , rifft16_512   }, { 1024, rifft16_1024  }, { 2048, rifft16_2048  }, { 4096, rifft16_4096  },
        { 8192, rifft16_8192  },
        { 160 , rinfft16_160  }, { 192 , rinfft16_192  }, { 240 , rinfft16_240  }, { 320 , rinfft16_320  },
        { 384 , rinfft16_384  }, { 480 , rinfft16_480  }
    };
    /* 16-bit handlers for complex FFT */
    const  N_handle_pair_t  cfft16x16_hndls[NUM_CFFT16X16] =
    {
        { 16  , cfft16_16    }, { 32  , cfft16_32    }, { 64  , cfft16_64    }, { 128 , cfft16_128   },
        { 256 , cfft16_256   }, { 512 , cfft16_512   }, { 1024, cfft16_1024  }, { 2048, cfft16_2048  },
        { 4096, cfft16_4096  },
        { 160 , cnfft16_160  }, { 192 , cnfft16_192  }, { 240 , cnfft16_240  }, { 320 , cnfft16_320  },
        { 384 , cnfft16_384  }, { 480 , cnfft16_480  }
    };
    const  N_handle_pair_t  cifft16x16_hndls[NUM_CFFT16X16] =
    {
        { 16  , cifft16_16    }, { 32  , cifft16_32    }, { 64  , cifft16_64    }, { 128 , cifft16_128   },
        { 256 , cifft16_256   }, { 512 , cifft16_512   }, { 1024, cifft16_1024  }, { 2048, cifft16_2048  },
        { 4096, cifft16_4096  },
        { 160 , cinfft16_160  }, { 192 , cinfft16_192  }, { 240 , cinfft16_240  }, { 320 , cinfft16_320  },
        { 384 , cinfft16_384  }, { 480 , cinfft16_480  }
    };
    /* 32x16-bit handlers for real FFT */
    const  N_handle_pair_t  rfft32x16_hndls[NUM_RFFT32X16] =
    {
        { 32  , rfft16_32    }, { 64  , rfft16_64    }, { 128 , rfft16_128   }, { 256 , rfft16_256   },
        { 512 , rfft16_512   }, { 1024, rfft16_1024  }, { 2048, rfft16_2048  }, { 4096, rfft16_4096  },
        { 8192, rfft16_8192  },
        { 160 , rnfft32x16_160  }, { 192 , rnfft32x16_192  }, { 240 , rnfft32x16_240  }, { 320 , rnfft32x16_320  },
        { 384 , rnfft32x16_384  }, { 480 , rnfft32x16_480  }
    };
    const  N_handle_pair_t  rifft32x16_hndls[NUM_RFFT32X16] =
    {
        { 32  , rifft16_32    }, { 64  , rifft16_64    }, { 128 , rifft16_128   }, { 256 , rifft16_256   },
        { 512 , rifft16_512   }, { 1024, rifft16_1024  }, { 2048, rifft16_2048  }, { 4096, rifft16_4096  },
        { 8192, rifft16_8192  },
        { 160 , rinfft32x16_160  }, { 192 , rinfft32x16_192  }, { 240 , rinfft32x16_240  }, { 320 , rinfft32x16_320  },
        { 384 , rinfft32x16_384  }, { 480 , rinfft32x16_480  }
    };
    /* 32x16-bit handlers for complex FFT */
    const  N_handle_pair_t  cfft32x16_hndls[NUM_CFFT32X16] =
    {
        { 16  , cfft16_16    }, { 32  , cfft16_32    }, { 64  , cfft16_64    }, { 128 , cfft16_128   },
        { 256 , cfft16_256   }, { 512 , cfft16_512   }, { 1024, cfft16_1024  }, { 2048, cfft16_2048  },
        { 4096, cfft16_4096  },
        { 160 , cnfft32x16_160  }, { 192 , cnfft32x16_192  }, { 240 , cnfft32x16_240  }, { 320 , cnfft32x16_320  },
        { 384 , cnfft32x16_384  }, { 480 , cnfft32x16_480  }
    };
    const  N_handle_pair_t  cifft32x16_hndls[NUM_CFFT32X16] =
    {
        { 16  , cifft16_16    }, { 32  , cifft16_32    }, { 64  , cifft16_64    }, { 128 , cifft16_128   },
        { 256 , cifft16_256   }, { 512 , cifft16_512   }, { 1024, cifft16_1024  }, { 2048, cifft16_2048  },
        { 4096, cifft16_4096  },
        { 160 , cinfft32x16_160  }, { 192 , cinfft32x16_192  }, { 240 , cinfft32x16_240  }, { 320 , cinfft32x16_320  },
        { 384 , cinfft32x16_384  }, { 480 , cinfft32x16_480  }
    };
    /* 24-bit handlers for real FFT */
#if 0 //HiFi3/3z API
    const  N_handle_pair_t  rifft24x24_hndls[NUM_RFFT24X24] =
    {
        { 32  , rifft24_32   },
        { 64  , rifft24_64   },
        { 128 , rifft24_128  },
        { 256 , rifft24_256  },
        { 512 , rifft24_512  },
        { 1024, rifft24_1024 },
        { 2048, rifft24_2048 },
        { 4096, rifft24_4096 },
        { 8192, rifft24_8192 }
    };
    const  N_handle_pair_t  rfft24x24_hndls[NUM_RFFT24X24] =
    {
        { 32  , rfft24_32   },
        { 64  , rfft24_64   },
        { 128 , rfft24_128  },
        { 256 , rfft24_256  },
        { 512 , rfft24_512  },
        { 1024, rfft24_1024 },
        { 2048, rfft24_2048 },
        { 4096, rfft24_4096 },
        { 8192, rfft24_8192 }
    };
    /* 24-bit handlers for complex FFT */
    const  N_handle_pair_t  cfft24x24_hndls[NUM_CFFT24X24] =
    {
        { 16  , cfft24_16   },
        { 32  , cfft24_32   },
        { 64  , cfft24_64   },
        { 128 , cfft24_128  },
        { 256 , cfft24_256  },
        { 512 , cfft24_512  },
        { 1024, cfft24_1024 },
        { 2048, cfft24_2048 },
        { 4096, cfft24_4096 }
    };
    const  N_handle_pair_t  cifft24x24_hndls[NUM_CFFT24X24] =
    {
        { 16  , cifft24_16   },
        { 32  , cifft24_32   },
        { 64  , cifft24_64   },
        { 128 , cifft24_128  },
        { 256 , cifft24_256  },
        { 512 , cifft24_512  },
        { 1024, cifft24_1024 },
        { 2048, cifft24_2048 },
        { 4096, cifft24_4096 }
    };
#endif
    const FFT_handle_tab_t h_tab[] =
    {
        {
            MAKE_FFT_ID(TE_FFT_TESTCASE_FORWARD, TE_FFT_REAL, FMT_FRACT32), NUM_RFFT32X32, rfft32x32_hndls
        }, 
        {
            MAKE_FFT_ID(TE_FFT_TESTCASE_INVERSE, TE_FFT_REAL, FMT_FRACT32), NUM_RFFT32X32, rifft32x32_hndls
        },
        {
            MAKE_FFT_ID(TE_FFT_TESTCASE_FORWARD, TE_FFT_REAL, FMT_FRACT16), NUM_RFFT16X16, rfft16x16_hndls
        },
        {
            MAKE_FFT_ID(TE_FFT_TESTCASE_INVERSE, TE_FFT_REAL, FMT_FRACT16), NUM_RFFT16X16, rifft16x16_hndls
        },
        {
            MAKE_FFT_ID(TE_FFT_TESTCASE_FORWARD, TE_FFT_REAL | TE_FFT_32X16, FMT_FRACT32), NUM_RFFT32X16, rfft32x16_hndls
        },
        {
            MAKE_FFT_ID(TE_FFT_TESTCASE_INVERSE, TE_FFT_REAL | TE_FFT_32X16, FMT_FRACT32), NUM_RFFT32X16, rifft32x16_hndls
        },
#if 0 //HiFi3/3z API
        {
            MAKE_FFT_ID(TE_FFT_TESTCASE_FORWARD, TE_FFT_REAL | TE_FFT_UNPACKED24, FMT_FRACT32), NUM_RFFT24X24, rfft24x24_hndls
        },
        {
            MAKE_FFT_ID(TE_FFT_TESTCASE_INVERSE, TE_FFT_REAL | TE_FFT_UNPACKED24, FMT_FRACT32), NUM_RFFT24X24, rifft24x24_hndls             
        },
#endif
        {
            MAKE_FFT_ID(TE_FFT_TESTCASE_FORWARD, TE_FFT_CPLX, FMT_FRACT32), NUM_CFFT32X32, cfft32x32_hndls
        },
        {
            MAKE_FFT_ID(TE_FFT_TESTCASE_INVERSE, TE_FFT_CPLX, FMT_FRACT32), NUM_CFFT32X32, cifft32x32_hndls
        },
        {
            MAKE_FFT_ID(TE_FFT_TESTCASE_FORWARD, TE_FFT_CPLX, FMT_FRACT16), NUM_CFFT16X16, cfft16x16_hndls
        },
        {
            MAKE_FFT_ID(TE_FFT_TESTCASE_INVERSE, TE_FFT_CPLX, FMT_FRACT16), NUM_CFFT16X16, cifft16x16_hndls
        },        
        {
            MAKE_FFT_ID(TE_FFT_TESTCASE_FORWARD, TE_FFT_CPLX | TE_FFT_32X16, FMT_FRACT32), NUM_CFFT32X16, cfft32x16_hndls
        },
        {
            MAKE_FFT_ID(TE_FFT_TESTCASE_INVERSE, TE_FFT_CPLX | TE_FFT_32X16, FMT_FRACT32), NUM_CFFT32X16, cifft32x16_hndls
        },
#if 0//HiFi3/3z API
        {
            MAKE_FFT_ID(TE_FFT_TESTCASE_FORWARD, (TE_FFT_CPLX | TE_FFT_UNPACKED24), FMT_FRACT32), NUM_CFFT24X24, cfft24x24_hndls
        },
        {
            MAKE_FFT_ID(TE_FFT_TESTCASE_INVERSE, (TE_FFT_CPLX | TE_FFT_UNPACKED24), FMT_FRACT32), NUM_CFFT24X24, cifft24x24_hndls
        },
#endif
    };

    for (i=0; i<H_TAB_SZ; i++)
    {
        if (h_tab[i].FFT_Id == FFT_Id)
        {
            for (j = 0; j<h_tab[i].numN; j++)
            {
                if (h_tab[i].pNh[j].N == N)
                {
                    return h_tab[i].pNh[j].h; 
                }
            }
        }
    }
    
    return NULL;
}/* GetFFT_handle() */


/* unpack 24-bit stream to 32-bit data */
static void unpack24(uint8_t *x,int N)
{
    int n;
    ASSERT(N%4==0);
    for (n=N/4-1; n>=0; n--)
    {
        uint8_t b0,b1,b2,b3;
        b0=0;
        b1=x[3*n+0];
        b2=x[3*n+1];
        b3=x[3*n+2];
        x[4*n+0]=b0;
        x[4*n+1]=b1;
        x[4*n+2]=b2;
        x[4*n+3]=b3;
    }
}

/* pack 32-bit data to 24-bit stream */
static void pack24(uint8_t *x,int N)
{
    int n;
    ASSERT(N%4==0);
    for (n=0; n<N/4; n++)
    {
        uint8_t b0,b1,b2;
        b0=x[4*n+1];
        b1=x[4*n+2];
        b2=x[4*n+3];
        x[3*n+0]=b0;
        x[3*n+1]=b1;
        x[3*n+2]=b2;
    }
}

/* Create a target algorithm instance and set tTestEngContext::target fields.
 * Return zero if failed. */
int te_create_fft( tTestEngContext * context )
{
  tTestEngContext_fft_int * context_fft;

  tVec * tblSet = 0;
  int baseSize, tblNum, res;

  /*
   * Allocate and initialize the context structure.
   */

  context_fft = (tTestEngContext_fft_int * )malloc( sizeof(*context_fft) );

  if ( !( res = ( 0 != context_fft ) ) )
  {
    printf( "te_create_fft(): malloc() failed\n" );
  }
  else
  {
    memset( context_fft, 0, sizeof(*context_fft) );
  }

  /*
   * Load twiddle factor tables.
   */

  if ( res )
  {
    int tblIx;
    tVec tbl_fp64c;
	const char *dir = getVectorsDir(context->isFull);
    char *fname;

    memset( &tbl_fp64c, 0, sizeof(tbl_fp64c) );

    res = 0;

    if ( 2 != seqFileScanf( context->seqFile, "%d %d", &baseSize, &tblNum ) )
    {
      printf( "te_create_fft(): bad SEQ-file format (1)\n" );
    }
    /* Allocate 64-bit vector of max FFT size. */
    else if ( !vecAlloc( &tbl_fp64c, baseSize<<(tblNum-1), 
                         TE_ALIGN_YES, FMT_FLOAT64|FMT_CPLX, 0 ) )
    {
      printf( "te_create_fft(): failed to allocate tbl_fp64c\n" );
    }
    /* Allocate vectors for all tables. */
    else if ( !( tblSet = (tVec*)calloc( tblNum, sizeof(tVec) ) ) )
    {
      printf( "te_create_fft(): failed to allocate the tables set\n" );
    }
    else res = 1;

    /* Read and convert twiddle tables. */
    fname=(char*)malloc(strlen(dir)+64+2+6);
    for ( tblIx=0; tblIx<tblNum && res; tblIx++ )
    {
      int tblSize = ( baseSize << tblIx )*3/4;
      int r;
      FILE * f = 0;
      int fmtTwd;
      char twdFname[64];
      fmtTwd = context->desc->fmt & FMT_DTYPE_MASK;
      if(context->desc->extraParam & TE_FFT_TWD16) fmtTwd = FMT_FRACT16;    /* load 16-bit twiddles instead of the same as input/output */

      res = 0;
      /* Read the twiddle factor table filename. */
      twdFname[0]=0;
      r=seqFileScanf( context->seqFile, "%63s", twdFname );
      sprintf(fname,"%s/fft_twd/%s",dir,twdFname);
      if ( !r )
      {
        printf( "te_create_fft(): bad SEQ-file format (2), tblIx=%d\n", tblIx );
      }
      /* Open the table file. */
      else if ( !( f = fopen( fname, "rb" ) ) )
      {
        printf( "te_create_fft(): failed to open %s for reading\n", fname );
      }
      /* Read twiddle factors in 64-bit FP format. */
      else if ( tblSize != (int)fread( vecGetElem( &tbl_fp64c, 0 ), sz_fp64c, tblSize, f ) )
      {
        printf( "te_create_fft(): failed to read %s\n", fname );
      }
      /* Allocate the twiddle table vector. */
      else if ( !vecAlloc( &tblSet[tblIx], tblSize, 
                           context->desc->isAligned, 
                           fmtTwd | FMT_CPLX, 0 ) )
      {
        printf( "te_create_fft(): failed to allocate twiddle table, tblIx=%d", tblIx );
      }
      else
      {
        res = 1;
        /* Convert twiddle factors to the target FFT format. */
        vecFromFp64( &tblSet[tblIx], (float64_t*)vecGetElem( &tbl_fp64c, 0 ) );
      }

      if ( f ) fclose( f );
    }
    free(fname);

    if ( !res && tblSet )
    {
      for ( tblIx=0; tblIx<tblNum; tblIx++ )
      {
        if ( tblSet[tblIx].szBulk > 0 ) vecFree( &tblSet[tblIx] );
      }

      free( tblSet );
    }

    if ( tbl_fp64c.szBulk > 0 ) vecFree( &tbl_fp64c );
  }

  if ( res )
  {
    context_fft->twd.baseSize = baseSize;
    context_fft->twd.tblNum   = tblNum;
    context_fft->twd.tblSet   = tblSet;

    context->target.handle = &context_fft->ext;
  }

  {
    typedef void (*tFxn)( const void * x, void * y, const void * twdTbl, int twdStep, int N );
    tFxn fxn_fwd,fxn_inv ;
    fxn_fwd = (tFxn)((const tTestEngDesc_fft *)context->desc)->frwTransFxn;
    fxn_inv = (tFxn)((const tTestEngDesc_fft *)context->desc)->invTransFxn;
    if ((NULL!=fxn_fwd && !IS_PRESENT(fxn_fwd)) || 
        (NULL!=fxn_inv && !IS_PRESENT(fxn_inv)) ||
        (NULL==fxn_fwd && NULL==fxn_inv))
    {
        // FUT is not defined
        return -1;
    }
  }
  return (res);

} /* te_create_fft() */

/* Destroy the target algorithm instance and free memory block(s) allocated
 * for the target object. Return zero whenever any violation of memory area
 * bounds is detected. */
int te_destroy_fft( tTestEngContext * context )
{
  tTestEngContext_fft_int * context_fft;

  ASSERT( context );

  if ( 0 != ( context_fft = (tTestEngContext_fft_int *)context->target.handle ) )
  {
    int tblIx;
    tVec * tblSet;

    if ( 0 != ( tblSet = context_fft->twd.tblSet ) )
    {
      for ( tblIx=0; tblIx<context_fft->twd.tblNum; tblIx++ )
      {
        if ( tblSet[tblIx].szBulk > 0 ) vecFree( &tblSet[tblIx] );
      }

      free( tblSet );
    }

    free( context_fft );
  }

  return (1);

} /* te_destroy_fft() */

/* Allocate in/out vectors for the next test case, and load the data set
 * from the SEQ-file. Return zero if failed. */
int te_load_fft( tTestEngContext * context )
{
  tTestEngContext_fft_int * context_fft = (tTestEngContext_fft_int *)context->target.handle;
  tVec BEXP, Z, Zlo, Zhi;
  int res = 0;
  int isFract = ( ( FMT_FRACT16 == ( context->desc->fmt & FMT_DTYPE_MASK ) ) ||
                  ( FMT_FRACT32 == ( context->desc->fmt & FMT_DTYPE_MASK ) ) );

  NASSERT( context_fft );

  memset( &context_fft->ext, 0, sizeof(context_fft->ext) );

  memset( &BEXP, 0, sizeof(BEXP) );
  memset( &Z   , 0, sizeof(Z   ) );
  memset( &Zlo , 0, sizeof(Zlo ) );
  memset( &Zhi , 0, sizeof(Zhi ) );

  /* If FFT supports the scaling option, read the scaling method from the SEQ-file. */
  if ( ( context->desc->extraParam & TE_FFT_OPT_SCALE_METH ) &&
       ( 1 != seqFileScanf( context->seqFile, "%d", &context_fft->ext.scale_method ) ) )
  {
    printf( "te_load_fft(): bad SEQ-file format (a)\n" );
  }
  /* For a fixed point blockwise FFT, allocate a vector for temporal storage of block exponent. */
  else if ( 0 != ( context->desc->extraParam & TE_FFT_BLOCKWISE ) && isFract &&
            !vecAlloc( &BEXP, context->args.L_DIM_, TE_ALIGN_NO, FMT_INT16, 0 ) )
  {
      printf("te_load_fft(): failed to allocate BEXP, L=%d\n", context->args.L_DIM_);
  }
  /* Read input data filename. */
  else if ( 1 != seqFileScanf( context->seqFile, "%63s", &context_fft->ext.fInName ) )
  {
    printf( "te_load_fft(): bad SEQ-file format (b)\n" );
  }
  /* Read reference data filename. */
  else if ( 1 != seqFileScanf( context->seqFile, "%63s", &context_fft->ext.fRefName ) )
  {
    printf( "te_load_fft(): bad SEQ-file format (c)\n" );
  }
  /* Allocate vectors for SINAD verification. */
  else if ( 3 != vecsAlloc( TE_ALIGN_NO, FMT_FLOAT32, &Z, 1, &Zlo, 1, &Zhi, 1, 0 ) )
  {
    printf( "te_load_fft(): failed to allocate vectors Z/Zlo/Zhi\n" );
  }
  /* Read the minimum SINAD value from the SEQ-file. */
  else if ( 1 != seqFileScanf( context->seqFile, "%f", vecGetElem_fl32( &Zlo, 0 ) ) )
  {
    printf( "te_load_fft(): bad SEQ-file format (d)\n" );
  }
  else
  {
    /* Set SINAD upper limit to infinity. */
    *vecGetElem_fl32( &Zhi, 0 ) = INFINITY;

    memset( &context->dataSet, 0, sizeof(context->dataSet) );

    context->dataSet.X   = BEXP;
    context->dataSet.Z   = Z;
    context->dataSet.Zlo = Zlo;
    context->dataSet.Zhi = Zhi;

    res = 1;
  }

  if ( !res )
  {
    if ( BEXP.szBulk ) vecFree( &BEXP );
    if ( Z   .szBulk ) vecFree( &Z    );
    if ( Zlo .szBulk ) vecFree( &Zlo  );
    if ( Zhi .szBulk ) vecFree( &Zhi  );
  }

  return (res);

} /* te_load_fft() */

/* Return a pointer to twiddle factor table. If step parameter
 * is non-zero, then the table is selected from the set of available
 * tables in dependence of the test frame counter, with appropriate
 * stride amount returned through step. If step is zero, return the
 * twiddle table such that stride is 1. Return zero if found no table
 * for the requested FFT size. */
void * te_get_twd_tbl( tTestEngContext * context, int fftSize, int * step )
{
  tTestEngContext_fft_int * context_fft;
  int baseSize, tblNum, cnt;
  int tblIx, stride = 0;

  ASSERT( context && context->target.handle );
  //ASSERT( !( fftSize & (fftSize-1) ) );

  context_fft = (tTestEngContext_fft_int *)context->target.handle;
  baseSize = context_fft->twd.baseSize;
  tblNum   = context_fft->twd.tblNum;

  if ( step )
  {
    tblIx = S_exp0_l( baseSize ) - S_exp0_l( fftSize );

    if ( tblIx < tblNum )
    {
      /* Select the "randomizer" among the frame counter (plain FFT test) or
       * test case number (2D FFT test). */
      cnt = ( context_fft->ext.frameCnt > 0 ? context_fft->ext.frameCnt : context->args.caseNum );
      /* Choose a table for size >= fftSize */
      tblIx += ( stride = ( cnt % ( tblNum - tblIx ) ) );
    }
  }
  else
  {
    tblIx = S_exp0_l( baseSize ) - S_exp0_l( fftSize );
  }

  if ( 0 <= tblIx && tblIx < tblNum )
  {
    if ( step ) *step = ( 1 << stride );

    return ( vecGetElem( &context_fft->twd.tblSet[tblIx], 0 ) );
  }

  return (0);

} /* te_get_twd_tbl() */

/* Apply the FFT function to a single frame of test data, any FFT routine excepting feature
 * rich fixed point FFTs and blockwise FFTs. */
int te_frameProc_stnd_fft( tTestEngContext * context, 
                     const int16_t         * in,
                           float64_t       * out,
                     tVec * xVec, tVec * yVec )
{
  typedef void (*tFxn)( const void * x, void * y, const void * twdTbl, int twdStep, int N, tTestEngContext_fft *);

  tFxn fxn = NULL;
  tTestEngContext_fft * context_fft;
  void *px, *py, *twdTbl;

  int bexp, shift;
  int N, logN, twdStep;
  int noReuse, doInplace, isRealForward;

  uint32_t crcSum = 0;

  NASSERT( context && context->desc && context->target.handle );
  NASSERT( in && out && xVec && yVec );
  NASSERT(0 == (context->args.N_DIM_ & (context->args.N_DIM_ - 1)));
  //NASSERT( 0 == ( context->desc->extraParam & TE_FFT_OPT_SCALE_METH ) );
  
  context_fft = (tTestEngContext_fft *)context->target.handle;
  N           = context->args.N_DIM_;
  logN        = 30 - S_exp0_l( N );

  /* If the FFT routine supports inplace operation, try it every second frame. */
  doInplace = ( context->desc->extraParam & TE_FFT_OPT_INPLACE ) && ( context_fft->frameCnt & 1 );
  /* Also check if the target FFT is allowed to reuse the input buffer for intermediate data. */
  noReuse = !( context->desc->extraParam & TE_FFT_OPT_REUSE_INBUF ) && !doInplace;
  /* Real-valued forward FFT requires special block exponent on input. */
  isRealForward = ( context->desc->extraParam & TE_FFT_REAL ) &&
                  ( context->args.caseType == TE_FFT_TESTCASE_FORWARD );

  /* For all fixed point FFTs, block exponent of input data should be at least 1, but for
   * real-valued forward FFT zero is allowed. */
  bexp = ( isRealForward ? 0 : 1 );

  /* Convert 16-bit PCM input data to target FFT format. */
  shift = vecFromPcm16( xVec, (int16_t*)in, bexp );

  /* Select in/out buffers for the FFT, and wipe the output buffer. */
  if ( doInplace )
  {
    if ( vecGetSize( xVec ) < vecGetSize( yVec ) )
    {
      memcpy( vecGetElem( yVec, 0 ), vecGetElem( xVec, 0 ), vecGetSize( xVec ) );
      px = py = vecGetElem( yVec, 0 );
    }
    else
    {
      px = py = vecGetElem( xVec, 0 );
    }
  }
  else
  {
    memset( vecGetElem( yVec, 0 ), 0, vecGetSize( yVec ) );
    px = vecGetElem( xVec, 0 );
    py = vecGetElem( yVec, 0 );
  }

  /* Select the target FFT routine (either forward or inverse). */
  if ( context->args.caseType == TE_FFT_TESTCASE_FORWARD )
  {
    fxn = (tFxn)((const tTestEngDesc_fft *)context->desc)->frwTransFxn;
    /* Compensate for scaling shift performed by fixed point FFTs. */
    shift -= logN;
  }
  else if ( context->args.caseType == TE_FFT_TESTCASE_INVERSE )
  {
    fxn = (tFxn)((const tTestEngDesc_fft *)context->desc)->invTransFxn;
    /* For fixed point inverse FFT, we have to divide output signal by FFT size.
     * Just don't compensate for the scaling shift performed by the FFT routine. */
  }
  else
  {
    NASSERT( !"Bad test case type!" );
  }

  /* If reuse is not allowed, make sure the buffer stays intact. */
  if ( noReuse ) crcSum = crc32( 0, (uint8_t*)px, vecGetSize( xVec ) );

  if (((tTestEngContext_fft_int *)context->target.handle)->twd.baseSize == 0)
  { 
      /* Some FFTs functions not use external twiddle tables, in this case baseSize
         (first twiddle table fft size) set to 0 in the *.seq file  */
      twdStep = 0;
      twdTbl = NULL;
  }
  else
  {    
      /* Select a twiddle factor table for FFT size >= N. */
      if ( !( twdTbl = te_get_twd_tbl( context, N, &twdStep ) ) )
      {
         printf( "te_frameProc_stnd_fft(): no twiddle factor table for N=%d\n", N );
        return (0);
      }
  }

  { /* verification reporting */
    char str[30];
    sprintf(str,"N=%d,scalingOption=%d",N,context_fft->scale_method);
    vReportAdd((tReportFUT*)&fxn,1,str,context_fft->fInName,context->args.caseType,context_fft->dataSize);
  }
  /* Apply the target FFT routine. */
  fxn(py, px, twdTbl, twdStep, N, context_fft);

  if ( doInplace && vecGetSize( xVec ) >= vecGetSize( yVec ) )
  {
    memcpy( vecGetElem( yVec, 0 ), py, vecGetSize( yVec ) );
  }
  
  if ( noReuse && crcSum != crc32( 0, (uint8_t*)px, vecGetSize( xVec ) ) )
  {
    printf( "te_frameProc_stnd_fft(): target FFT has corrupted the input buffer\n" );
    return ( 0 );
  }

  /* Convert output data to complex 64-bit floating point and rescale them. */
  vecToFp64( (float64_t*)out, yVec, shift );

  return (1);

} /* te_frameProc_stnd_fft() */


/* Apply the FFT function to a single frame of test data, feature rich fixed point
 * FFT (with scaling method option). */
int te_frameProc_stnd_scl_fft( tTestEngContext * context, 
                         const int16_t         * in,
                               float64_t       * out,
                         tVec * xVec, tVec * yVec)
{

  typedef int(*tFxn)(void * y, void * x, const void * twdTbl, ...);
  tFxn fxn = NULL;
  tTestEngContext_fft * context_fft;
  void *px, *py, *twdTbl;
  int bexp=0, shift;
  int N, logN, twdStep;
  int noReuse, doInplace;

  uint32_t crcSum = 0;

  NASSERT( context && context->desc && context->target.handle );
  NASSERT( in && out && xVec && yVec );
  NASSERT( 0 != ( context->desc->extraParam & TE_FFT_OPT_SCALE_METH ) );
  
  context_fft = (tTestEngContext_fft *)context->target.handle;
  N = context->args.N_DIM_;
  logN        = 30 - S_exp0_l( N );

  /* If the FFT routine supports inplace operation, try it every second frame. */
  doInplace = ( context->desc->extraParam & TE_FFT_OPT_INPLACE ) && ( context_fft->frameCnt & 1 );
  /* Also check if the target FFT is allowed to reuse the input buffer for intermediate data. */
  noReuse = !( context->desc->extraParam & TE_FFT_OPT_REUSE_INBUF ) && !doInplace;

    /* Select the target FFT routine. */
    if ( context->args.caseType == TE_FFT_TESTCASE_FORWARD )
    {
        fxn = (tFxn)((const tTestEngDesc_fft *)context->desc)->frwTransFxn;
    }
    else if ( context->args.caseType == TE_FFT_TESTCASE_INVERSE )
    {
        fxn = (tFxn)((const tTestEngDesc_fft *)context->desc)->invTransFxn;
    }
    else   
    {
    NASSERT( !"Bad test case type!" );
    }

    if (context_fft->scale_method == 0)
    {
        // Set prescaling if scaling option is no-scaling
        bexp = logN + 1;
    }
    else
    {
        bexp = 0;
    }
 

    if ( (xVec->fmt & FMT_DTYPE_MASK) == FMT_FRACT16  &&  context_fft->scale_method == 2)
    {
        /* For fract16 FFT format with dynamic scaling copy input data as is */
        int num_int16 = (xVec->fmt & FMT_CPLX) ? 2 * xVec->nElem : xVec->nElem; 
        memcpy(vecGetElem(xVec, 0), in, sizeof(int16_t) * num_int16);
        shift = 0; 
    }
    else if (   (xVec->fmt & FMT_DTYPE_MASK) == FMT_FRACT32  &&  context_fft->scale_method == 2 )
    {
        int k; 
        int32_t *pdst = (int32_t*)vecGetElem(xVec, 0);     
        int num_int32 = (xVec->fmt & FMT_CPLX) ? 2 * xVec->nElem : xVec->nElem;

        /* For FFTs with 32-bit dynamic scaling - extend dynamic range to 32 bits */
        shift = -16 * (context_fft->frameCnt & 1);
        for (k = 0; k < num_int32; k++)
        {
            pdst[k] = ((int32_t)in[k]) << (16 + shift);
        }        
    }
#if 1
    else if ((xVec->fmt & FMT_DTYPE_MASK) == FMT_FRACT32  &&  context_fft->scale_method == 1)
    {
        int k;
        int32_t *pdst = (int32_t*)vecGetElem(xVec, 0);
        int num_int32 = (xVec->fmt & FMT_CPLX) ? 2 * xVec->nElem : xVec->nElem;

        /* For FFTs with 24-bit dynamic scaling, put data into high 16 bits */
        shift = 0;
        for (k = 0; k < num_int32; k++)
        {
            pdst[k] = ((int32_t)in[k]) << (16 + shift);
        }
    }
#endif
    else
    {
        /* Convert 16-bit PCM input data to target FFT format. */
        shift = vecFromPcm16( xVec, (int16_t*)in, bexp );
    }

  /* Select in/out buffers for the FFT, and wipe the output buffer. */
  if ( doInplace )
  {
    if ( vecGetSize( xVec ) < vecGetSize( yVec ) )
    {
      memcpy( vecGetElem( yVec, 0 ), vecGetElem( xVec, 0 ), vecGetSize( xVec ) );
      px = py = vecGetElem( yVec, 0 );
    }
    else
    {
      px = py = vecGetElem( xVec, 0 );
    }
  }
  else
  {
    memset( vecGetElem( yVec, 0 ), 0, vecGetSize( yVec ) );
    px = vecGetElem( xVec, 0 );
    py = vecGetElem( yVec, 0 );
  }

  /* If reuse is not allowed, make sure the buffer stays intact. */
  if ( noReuse ) crcSum = crc32( 0, (uint8_t*)px, vecGetSize( xVec ) );

  /* Select a twiddle factor table for FFT size >= N. */
  if ( !( twdTbl = te_get_twd_tbl( context, N, &twdStep ) ) )
  {
  //  printf( "te_frameProc_stnd_scl_fft(): no twiddle factor table for N=%d\n", N );
   // return (0);
  }

  /* Apply the target FFT routine. */
  if (context->desc->extraParam & TE_FFT_PACKED24)
  {
      /* reallocate xVec, yVec to fit largest size */
      int res;
      tVec XVec,YVec;
      void* pX,*pY;
      size_t szX=vecGetSize( xVec ),szY=vecGetSize( yVec );
      size_t maxsz=MAX(szX,szY);
      res =vecAlloc(&XVec, maxsz/xVec->szElem, ((uintptr_t)px)%8==0, xVec->fmt, NULL);
      res&=vecAlloc(&YVec, maxsz/yVec->szElem, ((uintptr_t)py)%8==0, yVec->fmt, NULL);
      if(!res)
      {
        printf( "te_frameProc_stnd_scl_fft(): does not able to allocate memory\n");
        return (0);
      }
      pX=vecGetElem( &XVec, 0 );
      pY=vecGetElem( &YVec, 0 );
      memset(pX,0xcc,maxsz);
      memset(pY,0xcc,maxsz);
      memcpy( pX, px, szX );
      pack24((uint8_t*)pX, szX);
      bexp = fxn(pY, pX, twdTbl, twdStep, N, context_fft->scale_method, context_fft);
      unpack24((uint8_t*)pY, szY);
      memcpy( py, pY, szY );
      vecsFree(&XVec,&YVec,NULL);
  }
  else
  {
      { /* verification reporting */
        char str[30];
        sprintf(str,"N=%d,scalingOption=%d",N,context_fft->scale_method);
        vReportAdd((tReportFUT*)&fxn,1,str,context_fft->fInName,context->args.caseType,context_fft->dataSize);
      }
      if (twdTbl != NULL)
      {
         bexp = fxn( py, px, twdTbl, twdStep, N, context_fft->scale_method, context_fft);
      }
      else
      {
          const void *h = GetFFT_handle(N, MAKE_FFT_ID(context->args.caseType, context->desc->extraParam, context->desc->fmt));
          if (h==NULL)
          {
              printf( "te_frameProc_stnd_scl_fft(): corresponding FFT handle not found\n" );
              return ( 0 );
          }
          bexp = fxn(py, px, h, context_fft->scale_method);
      }
  }

  if ( doInplace && vecGetSize( xVec ) >= vecGetSize( yVec ) )
  {
    memcpy( vecGetElem( yVec, 0 ), py, vecGetSize( yVec ) );
  }
  
  if ( noReuse && crcSum != crc32( 0, (uint8_t*)px, vecGetSize( xVec ) ) )
  {
    printf( "te_frameProc_stnd_scl_fft(): target FFT has corrupted the input buffer\n" );
    return ( 0 );
  }

  /* Compensate for scaling shift performed by fixed point forward FFTs. */
  shift -= bexp;
  /* For inverse FFTs,  we also have to divide output signal by FFT size. */
  if (context->args.caseType == TE_FFT_TESTCASE_INVERSE && !(N&(N - 1)))
  {
      // Doesn't work  for inverse mixed radix FFT
      shift += logN;
  }
  /* Convert output data to complex 64-bit floating point and rescale them. */
  vecToFp64( (float64_t*)out, yVec, shift );

  if (context->args.caseType == TE_FFT_TESTCASE_INVERSE && (N&(N - 1)) > 0 )
  {
      // Scaling by 1/N  for inverse mixed radix FFT
      int k; 
      float64_t f = 1.0/N; 
      int num_out = (yVec->fmt&FMT_CPLX)? 2 * yVec->nElem: yVec->nElem; 
      for (k = 0; k < num_out; k++) out[k] *= f;
  }



  return (1);

} /* te_frameProc_stnd_scl_fft() */

/* Apply the FFT function to a single frame of test data, any blockwise FFT routine. */
int te_frameProc_stnd_blkfft( tTestEngContext * context, 
                     const int16_t         * in,
                           float64_t       * out,
                     tVec * xVec, tVec * yVec) 
{
  tTestEngContext_fft* context_fft = (tTestEngContext_fft *)context->target.handle;
  typedef void (*tFxn)( const void * x, void * y, const void * twdTbl, int twdStep, int L, int N , tTestEngContext_fft*);

  tFxn fxn = NULL;
  /*tTestEngContext_fft * context_fft; */
  void *px, *py, *twdTbl;
  int16_t * shift = 0;

  int bexp;
  int N, logN, twdStep;
  int l, L;
  int isRealForward, isRealInverse;

  NASSERT( context && context->desc && context->target.handle );
  NASSERT( in && out && xVec && yVec );
  NASSERT(0 == (context->args.N_DIM_ & (context->args.N_DIM_ - 1)));
  NASSERT( 0 != ( context->desc->extraParam & TE_FFT_BLOCKWISE ) );
  
  /*context_fft = (tTestEngContext_fft *)context->target.handle;*/
  N = context->args.N_DIM_;
  logN        = 30 - S_exp0_l( N );

  isRealForward = ( context->desc->extraParam & TE_FFT_REAL ) &&
                  ( context->args.caseType == TE_FFT_TESTCASE_FORWARD );
  isRealInverse = ( context->desc->extraParam & TE_FFT_REAL ) &&
                  ( context->args.caseType == TE_FFT_TESTCASE_INVERSE );

  L = ( xVec->nElem >> ( isRealInverse ? logN-1 : logN ) );

  /* For all fixed point FFTs, block exponent of input data should be at least 1, but for
   * real-valued forward FFT zero is allowed. */
  bexp = ( isRealForward ? 0 : 1 );
  /* Prepare shift amount array for fixed point data. */
  if ( context->dataSet.X.szBulk > 0 ) shift = vecGetElem_i16( &context->dataSet.X, 0 );

  /* Convert 16-bit PCM input data to target FFT format. */
  bvecFromPcm16( xVec, (int16_t*)in, ( isRealInverse ? N/2 : N ), L, bexp, shift );

  /* Select in/out buffers for the FFT, and wipe the output buffer. */
  px = vecGetElem( xVec, 0 );
  py = vecGetElem( yVec, 0 );
  memset( py, 0, vecGetSize( yVec ) );

  /* Select the target FFT routine (either forward or inverse). */
  if ( context->args.caseType == TE_FFT_TESTCASE_FORWARD )
  {
      fxn = (tFxn)((const tTestEngDesc_fft *)context->desc)->frwTransFxn;
    /* Compensate for scaling shift performed by fixed point FFTs. */
    if ( shift ) for ( l=0; l<L; l++ ) shift[l] -= logN;
  }
  else if ( context->args.caseType == TE_FFT_TESTCASE_INVERSE )
  {
      fxn = (tFxn)((const tTestEngDesc_fft *)context->desc)->invTransFxn;
    /* For fixed point inverse FFT, we have to divide output signal by FFT size.
     * Just don't compensate for the scaling shift performed by the FFT routine. */
  }
  else
  {
    NASSERT( !"Bad test case type!" );
  }

  if (((tTestEngContext_fft_int *)context->target.handle)->twd.baseSize == 0)
  {
      /* Some FFTs functions not use external twiddle tables, in this case baseSize
      (first twiddle table fft size) set to 0 in the *.seq file  */
      twdStep = 0;
      twdTbl = NULL;
  }
  else
  {
      /* Select a twiddle factor table for FFT size >= N. */
      if ( !( twdTbl = te_get_twd_tbl( context, N, &twdStep ) ) )
      {
        printf( "te_frameProc_stnd_blkfft(): no twiddle factor table for N=%d\n", N );
        return (0);
      }
  }
  /* Apply the target FFT routine. */
  fxn( px, py, twdTbl, twdStep, L, N , context_fft);

  /* Convert output data to complex 64-bit floating point and rescale them. */
  bvecToFp64( (float64_t*)out, yVec, ( isRealForward ? N/2 : N ), L, shift );

  return (1);

} /* te_frameProc_stnd_blkfft() */

/* Test engine methods for 2D FFT tests:
 *   Allocate vectors for the next in turn test case, and load the data
 *   set from a SEQ-file. Return zero if failed. */
int te_load_fft2d( tTestEngContext * context )
{
  int N, nElem, res = 0;
  int isCplx, baseFmt, fmtX, fmtY;

  ASSERT( context );

  N = context->args.N_DIM_;
  /* Don't allocate data vectors for an out-of-domain FFT size, or not a power of two. */
  nElem = ( ( N>=4 && N<=128 && !(N&(N-1)) ) ? N*N : 0 );

  isCplx = ( 0 != ( context->desc->fmt & FMT_CPLX ) );
  baseFmt  = ( context->desc->fmt & FMT_DTYPE_MASK );

  fmtX = ( ( isCplx || context->desc->extraParam == TE_FFT_TESTCASE_INVERSE ) ? FMT_CPLX : FMT_REAL );
  fmtY = ( ( isCplx || context->desc->extraParam == TE_FFT_TESTCASE_FORWARD ) ? FMT_CPLX : FMT_REAL );

  memset( &context->dataSet, 0, sizeof(context->dataSet) );

  /* Allocate input data vector memory, use target data format. */
  if ( !vecAlloc( &context->dataSet.X, nElem, context->desc->isAligned, baseFmt|fmtX, 0 ) )
  {
    printf( "te_load_fft2d(): failed to allocate vector X; "
            "fmt = 0x%02x, nElem = %d\n", (unsigned)(baseFmt|fmtX), nElem );
  }
  /* Allocate reference data vector memory. Always double precision FP! */
  else if ( !vecAlloc( &context->dataSet.Y, nElem, TE_ALIGN_YES,
                       FMT_FLOAT64 | fmtY, 0 ) )
  {
    printf( "te_load_fft2d(): failed to allocate vector Y; "
            "fmt = 0x%02x, nElem = %d\n", FMT_FLOAT64 | fmtY, nElem );
  }
  /* Allocate a single element vectors for SINAD values, single precision FP. */
  else if ( 3 != vecsAlloc( TE_ALIGN_YES, FMT_FLOAT32,
                            &context->dataSet.Z, 1,
                            &context->dataSet.Zlo, 1,
                            &context->dataSet.Zhi, 1, 0 ) )
  {
    printf( "te_load_fft2d(): failed to allocate vectors Z/Zlo/Zhi; "
            "fmt = 0x%02x, nElem = 1\n", FMT_FLOAT32 );
  }
  /* Load input and reference data vectors from the SEQ-file. */
  else if ( !seqFileReadVecs( context->seqFile, &context->dataSet.X,
                                                &context->dataSet.Y, 0 ) )
  {
    printf( "te_load_fft2d(): failed to read input/reference vectors data; "
            "fmt = 0x%02x, nElem = %d\n", (unsigned)context->desc->fmt, nElem );
  }
  /* Read the SINAD threshold. */
  else if ( !seqFileReadVec( context->seqFile, &context->dataSet.Zlo ) )
  {
    printf( "te_load_fft2d(): failed to read SINAD threshold\n" );
  }
  else res = 1;

  if ( !res )
  {
    if ( context->dataSet.X  .szBulk ) vecFree( &context->dataSet.X   );
    if ( context->dataSet.Y  .szBulk ) vecFree( &context->dataSet.Y   );
    if ( context->dataSet.Z  .szBulk ) vecFree( &context->dataSet.Z   );
    if ( context->dataSet.Zlo.szBulk ) vecFree( &context->dataSet.Zlo );
    if ( context->dataSet.Zhi.szBulk ) vecFree( &context->dataSet.Zhi );
  }
  else
  {
    /* We allow any SINAD value greater than or equal to the threshold value:
     * upper limit is set to infinity. */
    *vecGetElem_fl32( &context->dataSet.Zhi, 0 ) = INFINITY;
  }

  return (res);

} /* te_load_fft2d() */

/* Test engine methods for 2D FFT tests:
 *   Apply the target function to test case data. */
void te_process_fft2d( tTestEngContext * context )
{
  typedef void (*tFxn)( const void * x, void * temp, void * y,
                        const void * twiddle_table, int twiddle_stride, int N );

  tFxn fxn = NULL;
  void * px = 0, * py = 0, * ptemp = 0;
  tVec tempVec, resVec, resVecFp64;
  int isCplx, baseFmt, outFmt;
  int doInplace;
  int n, N, nElem;
  const void * twdTbl;
  int twdStep;
  uint32_t crc;

  ASSERT( context );

  memset( &tempVec   , 0, sizeof(tempVec   ) );
  memset( &resVec    , 0, sizeof(resVec    ) );
  memset( &resVecFp64, 0, sizeof(resVecFp64) );

  N = context->args.N_DIM_;
  /* Don't allocate data vectors for an out-of-domain FFT size, or not a power of two. */
  nElem = ( ( N>=4 && N<=128 && !(N&(N-1)) ) ? N*N : 0 );

  /* Set invalid SINAD result that will indicate an abnormal completion of the
   * test procedure. */
  *vecGetElem_fl32( &context->dataSet.Z, 0 ) = -INFINITY;

  /* Select a twiddle factor table for FFT size >= N. */
  if ( !( twdTbl = te_get_twd_tbl( context, ( nElem>0 ? N : 4 ), &twdStep ) ) )
  {
    printf( "te_process_fft2d(): no twiddle factor table for N=%d\n", N );
    return;
  }

  /* Allocate the scratch vector. */
  if ( !vecAlloc( &tempVec, nElem, TE_ALIGN_NO, context->desc->fmt | FMT_CPLX, 0 ) )
  {
    printf( "te_process_fft2d(): failed to allocate temp vector memory; "
            "fmt = 0x%02x, nElem = %d\n", (unsigned)( context->desc->fmt | FMT_CPLX ), nElem );
    return;
  }
  else ptemp = vecGetElem( &tempVec, 0 );

  /* Complex-valued 2D FFT? */
  isCplx = ( 0 != ( context->desc->fmt & FMT_CPLX ) );
  /* Base data type and width (fract<16|32>, float<32|64>_t). */
  baseFmt = ( context->desc->fmt & FMT_DTYPE_MASK );
  /* Output data format: complex or real. */
  outFmt = ( ( isCplx || context->desc->extraParam == TE_FFT_TESTCASE_FORWARD ) ? FMT_CPLX : FMT_REAL );
  /* For even-numbered test cases, try in-place operation. */
  doInplace = !( context->args.caseNum & 1 );
  /* Target FFT routine. */
  fxn = (tFxn)context->target.fut;

  /* Out-of-place operation. */
  if ( !doInplace )
  {
    /* Allocate the results vector. */
    if ( !vecAlloc( &resVec, nElem, TE_ALIGN_NO, baseFmt|outFmt, 0 ) )
    {
      printf( "te_process_fft2d(): failed to allocate output vector memory; "
              "fmt = 0x%02x, nElem = %d\n", (unsigned)(baseFmt|outFmt), nElem );
      return;
    }

    px = vecGetElem( &context->dataSet.X, 0 );
    py = vecGetElem( &resVec, 0 );
    /* Protect input data of corruption with a CRC checksum. */
    crc = crc32( 0, (uint8_t*)px, nElem*context->dataSet.X.szElem );
    /* Wipe the output buffer. */
    memset( py, 0, nElem*resVec.szElem );

    /* Invoke the target function. */
    fxn( px, ptemp, py, twdTbl, twdStep, N );

    /* Check input buffer integrity. */
    if ( crc != crc32( 0, (uint8_t*)px, nElem*context->dataSet.X.szElem ) )
    {
      printf( "te_process_fft2d(): target FFT routine spoiled the input buffer\n" );
      return;
    }
  }
  /* In-place operation; output data size greater than or equal to input data size. */
  else if ( outFmt == FMT_CPLX )
  {
    /* Allocate the results vector. */
    if ( !vecAlloc( &resVec, nElem, TE_ALIGN_NO, baseFmt|outFmt, 0 ) )
    {
      printf( "te_process_fft2d(): failed to allocate output vector memory; "
              "fmt = 0x%02x, nElem = %d\n", (unsigned)(baseFmt|outFmt), nElem );
      return;
    }

    px = vecGetElem( &context->dataSet.X, 0 );
    py = vecGetElem( &resVec, 0 );
    /* Copy input data to the output buffer. */
    memcpy( py, px, nElem*context->dataSet.X.szElem );

    /* Invoke the target function. */
    fxn( py, ptemp, py, twdTbl, twdStep, N );
  }
  /* In-place operation; output data size less than input data size. */
  else
  {
    px = vecGetElem( &context->dataSet.X, 0 );
    /* Invoke the target function. */
    fxn( px, ptemp, px, twdTbl, twdStep, N );

    /* Allocate the results vector and initialize it with FFT image disposed in
     * the input buffer. */
    if ( !vecAlloc( &resVec, nElem, TE_ALIGN_NO, baseFmt|outFmt, px ) )
    {
      printf( "te_process_fft2d(): failed to allocate output vector memory; "
              "fmt = 0x%02x, nElem = %d\n", (unsigned)(baseFmt|outFmt), nElem );
      return;
    }
  }

  /* Allocate vector for double precision FP Fourier image. */
  if ( !vecAlloc( &resVecFp64, nElem, TE_ALIGN_YES, FMT_FLOAT64 | outFmt, 0 ) )
  {
    printf( "te_process_fft2d(): failed to allocate DP FP result vector memory; "
            "fmt = 0x%02x, nElem = %d\n", (unsigned)( FMT_FLOAT64 | outFmt ), nElem );
    return;
  }

  /* Convert the results to double precision FP. */
  vecToFp64( (float64_t*)vecGetElem( &resVecFp64, 0 ), &resVec, 0 );

  /*
   * Compare results to reference data and estimate the SINAD.
   */

  {
    const float64_t *pref, *pres;

    float64_t err, p, q;
    float32_t sinad;

    pref = (float64_t*)vecGetElem( &context->dataSet.Y, 0 );
    pres = (float64_t*)vecGetElem( &resVecFp64        , 0 );

    N = ( outFmt == FMT_CPLX ? 2*nElem : nElem );

    for ( p=q=0., n=0; n<N; n++ )
    {
      err = pres[n] - pref[n];

      /* |signal|^2 */
      p += pref[n]*pref[n];
      /* |noise+distortion|^2 */
      q += err*err;
    }

    /* Take care of 0/0 and NaNs! */ 
    sinad = 10.f*log10f( q == 0. ? INFINITY : (float32_t)(p/q) );
    /* SINAD is to be validated by the Test Engine. */
    *vecGetElem_fl32( &context->dataSet.Z, 0 ) = sinad;
  }

  /* Free vectors. */
  if ( tempVec   .szBulk > 0 ) vecFree( &tempVec    );
  if ( resVec    .szBulk > 0 ) vecFree( &resVec     );
  if ( resVecFp64.szBulk > 0 ) vecFree( &resVecFp64 );

}
