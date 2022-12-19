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
 * Test procedures for matrix inversion and related functions
 */

#include <string.h>
#include <stdlib.h>

/* Cross-platform data type definitions. */
#include "types.h"
/* Test environment configuration. */
#include "config.h"
#include "packages.h"
/* DSP Library API. */
#include LIBRARY_HEADER(matinv)
/* Test engine API. */
#include "../../common/testeng_matinv.h"

static const tMtxInvApi mtx_inv2x2_32x32_Api    = {0,mtx_inv2x2_32x32_getScratchSize  };
static const tMtxInvApi mtx_inv3x3_32x32_Api    = {0,mtx_inv3x3_32x32_getScratchSize  };
static const tMtxInvApi mtx_inv4x4_32x32_Api    = {0,mtx_inv4x4_32x32_getScratchSize  };
static const tMtxInvApi mtx_inv6x6_32x32_Api    = {0,mtx_inv6x6_32x32_getScratchSize  };
static const tMtxInvApi mtx_inv8x8_32x32_Api    = {0,mtx_inv8x8_32x32_getScratchSize  };
static const tMtxInvApi mtx_inv10x10_32x32_Api  = {0,mtx_inv10x10_32x32_getScratchSize};
static const tMtxInvApi cmtx_inv2x2_32x32_Api   = {0,cmtx_inv2x2_32x32_getScratchSize  };
static const tMtxInvApi cmtx_inv3x3_32x32_Api   = {0,cmtx_inv3x3_32x32_getScratchSize  };
static const tMtxInvApi cmtx_inv4x4_32x32_Api   = {0,cmtx_inv4x4_32x32_getScratchSize  };
static const tMtxInvApi cmtx_inv6x4_32x32_Api   = {0,cmtx_inv6x6_32x32_getScratchSize  };
static const tMtxInvApi cmtx_inv8x4_32x32_Api   = {0,cmtx_inv8x8_32x32_getScratchSize  };
static const tMtxInvApi cmtx_inv10x10_32x32_Api = {0,cmtx_inv10x10_32x32_getScratchSize};
static const tMtxInvApi mtx_gjelim2x2_32x32_Api    = {0,mtx_gjelim2x2_32x32_getScratchSize  };
static const tMtxInvApi mtx_gjelim3x3_32x32_Api    = {0,mtx_gjelim3x3_32x32_getScratchSize  };
static const tMtxInvApi mtx_gjelim4x4_32x32_Api    = {0,mtx_gjelim4x4_32x32_getScratchSize  };
static const tMtxInvApi mtx_gjelim6x6_32x32_Api    = {0,mtx_gjelim6x6_32x32_getScratchSize  };
static const tMtxInvApi mtx_gjelim8x8_32x32_Api    = {0,mtx_gjelim8x8_32x32_getScratchSize  };
static const tMtxInvApi mtx_gjelim10x10_32x32_Api  = {0,mtx_gjelim10x10_32x32_getScratchSize};
static const tMtxInvApi cmtx_gjelim2x2_32x32_Api   = {0,cmtx_gjelim2x2_32x32_getScratchSize  };
static const tMtxInvApi cmtx_gjelim3x3_32x32_Api   = {0,cmtx_gjelim3x3_32x32_getScratchSize  };
static const tMtxInvApi cmtx_gjelim4x4_32x32_Api   = {0,cmtx_gjelim4x4_32x32_getScratchSize  };
static const tMtxInvApi cmtx_gjelim6x6_32x32_Api   = {0,cmtx_gjelim6x6_32x32_getScratchSize  };
static const tMtxInvApi cmtx_gjelim8x8_32x32_Api   = {0,cmtx_gjelim8x8_32x32_getScratchSize  };
static const tMtxInvApi cmtx_gjelim10x10_32x32_Api = {0,cmtx_gjelim10x10_32x32_getScratchSize};

static const tTestEngDesc descr_mtx_gjelim2x2_32x32 =  { FMT_REAL | FMT_FLOAT32, MTXINV_PLAIN,&mtx_gjelim2x2_32x32_Api  ,TE_DIM_NUM_1, TE_ALIGN_NO, NULL, NULL, &te_loadFxn_mtxgjelim, &te_processFxn_gjelim };
static const tTestEngDesc descr_mtx_gjelim3x3_32x32 =  { FMT_REAL | FMT_FLOAT32, MTXINV_PLAIN,&mtx_gjelim3x3_32x32_Api  ,TE_DIM_NUM_1, TE_ALIGN_NO, NULL, NULL, &te_loadFxn_mtxgjelim, &te_processFxn_gjelim };
static const tTestEngDesc descr_mtx_gjelim4x4_32x32 =  { FMT_REAL | FMT_FLOAT32, MTXINV_PLAIN,&mtx_gjelim4x4_32x32_Api  ,TE_DIM_NUM_1, TE_ALIGN_NO, NULL, NULL, &te_loadFxn_mtxgjelim, &te_processFxn_gjelim };
static const tTestEngDesc descr_mtx_gjelim6x6_32x32 =  { FMT_REAL | FMT_FLOAT32, MTXINV_PLAIN,&mtx_gjelim6x6_32x32_Api  ,TE_DIM_NUM_1, TE_ALIGN_NO, NULL, NULL, &te_loadFxn_mtxgjelim, &te_processFxn_gjelim };
static const tTestEngDesc descr_mtx_gjelim8x8_32x32 =  { FMT_REAL | FMT_FLOAT32, MTXINV_PLAIN,&mtx_gjelim8x8_32x32_Api  ,TE_DIM_NUM_1, TE_ALIGN_NO, NULL, NULL, &te_loadFxn_mtxgjelim, &te_processFxn_gjelim };
static const tTestEngDesc descr_mtx_gjelim10x10_32x32 ={ FMT_REAL | FMT_FLOAT32, MTXINV_PLAIN,&mtx_gjelim10x10_32x32_Api,TE_DIM_NUM_1, TE_ALIGN_NO, NULL, NULL, &te_loadFxn_mtxgjelim, &te_processFxn_gjelim };

static const tTestEngDesc descr_cmtx_gjelim2x2_32x32 =  { FMT_CPLX | FMT_FLOAT32, MTXINV_PLAIN,&cmtx_gjelim2x2_32x32_Api  ,TE_DIM_NUM_1, TE_ALIGN_NO, NULL, NULL, &te_loadFxn_mtxgjelim, &te_processFxn_gjelim };
static const tTestEngDesc descr_cmtx_gjelim3x3_32x32 =  { FMT_CPLX | FMT_FLOAT32, MTXINV_PLAIN,&cmtx_gjelim3x3_32x32_Api  ,TE_DIM_NUM_1, TE_ALIGN_NO, NULL, NULL, &te_loadFxn_mtxgjelim, &te_processFxn_gjelim };
static const tTestEngDesc descr_cmtx_gjelim4x4_32x32 =  { FMT_CPLX | FMT_FLOAT32, MTXINV_PLAIN,&cmtx_gjelim4x4_32x32_Api  ,TE_DIM_NUM_1, TE_ALIGN_NO, NULL, NULL, &te_loadFxn_mtxgjelim, &te_processFxn_gjelim };
static const tTestEngDesc descr_cmtx_gjelim6x6_32x32 =  { FMT_CPLX | FMT_FLOAT32, MTXINV_PLAIN,&cmtx_gjelim6x6_32x32_Api  ,TE_DIM_NUM_1, TE_ALIGN_NO, NULL, NULL, &te_loadFxn_mtxgjelim, &te_processFxn_gjelim };
static const tTestEngDesc descr_cmtx_gjelim8x8_32x32 =  { FMT_CPLX | FMT_FLOAT32, MTXINV_PLAIN,&cmtx_gjelim8x8_32x32_Api  ,TE_DIM_NUM_1, TE_ALIGN_NO, NULL, NULL, &te_loadFxn_mtxgjelim, &te_processFxn_gjelim };
static const tTestEngDesc descr_cmtx_gjelim10x10_32x32 ={ FMT_CPLX | FMT_FLOAT32, MTXINV_PLAIN,&cmtx_gjelim10x10_32x32_Api,TE_DIM_NUM_1, TE_ALIGN_NO, NULL, NULL, &te_loadFxn_mtxgjelim, &te_processFxn_gjelim };

static const tTestEngDesc descr_mtx_inv2x2_32x32 =  { FMT_REAL | FMT_FLOAT32, MTXINV_PLAIN,&mtx_inv2x2_32x32_Api  ,TE_DIM_NUM_1, TE_ALIGN_NO, NULL, NULL, &te_loadFxn_mtxinv, &te_processFxn_matinv };
static const tTestEngDesc descr_mtx_inv3x3_32x32 =  { FMT_REAL | FMT_FLOAT32, MTXINV_PLAIN,&mtx_inv3x3_32x32_Api  ,TE_DIM_NUM_1, TE_ALIGN_NO, NULL, NULL, &te_loadFxn_mtxinv, &te_processFxn_matinv };
static const tTestEngDesc descr_mtx_inv4x4_32x32 =  { FMT_REAL | FMT_FLOAT32, MTXINV_PLAIN,&mtx_inv4x4_32x32_Api  ,TE_DIM_NUM_1, TE_ALIGN_NO, NULL, NULL, &te_loadFxn_mtxinv, &te_processFxn_matinv };
static const tTestEngDesc descr_mtx_inv6x6_32x32 =  { FMT_REAL | FMT_FLOAT32, MTXINV_PLAIN,&mtx_inv6x6_32x32_Api  ,TE_DIM_NUM_1, TE_ALIGN_NO, NULL, NULL, &te_loadFxn_mtxinv, &te_processFxn_matinv };
static const tTestEngDesc descr_mtx_inv8x8_32x32 =  { FMT_REAL | FMT_FLOAT32, MTXINV_PLAIN,&mtx_inv8x8_32x32_Api  ,TE_DIM_NUM_1, TE_ALIGN_NO, NULL, NULL, &te_loadFxn_mtxinv, &te_processFxn_matinv };
static const tTestEngDesc descr_mtx_inv10x10_32x32 ={ FMT_REAL | FMT_FLOAT32, MTXINV_PLAIN,&mtx_inv10x10_32x32_Api,TE_DIM_NUM_1, TE_ALIGN_NO, NULL, NULL, &te_loadFxn_mtxinv, &te_processFxn_matinv };
static const tTestEngDesc descr_cmtx_inv2x2_32x32 =  { FMT_CPLX | FMT_FLOAT32, MTXINV_PLAIN,&cmtx_inv2x2_32x32_Api  ,TE_DIM_NUM_1, TE_ALIGN_NO, NULL, NULL, &te_loadFxn_mtxinv, &te_processFxn_matinv };
static const tTestEngDesc descr_cmtx_inv3x3_32x32 =  { FMT_CPLX | FMT_FLOAT32, MTXINV_PLAIN,&cmtx_inv3x3_32x32_Api  ,TE_DIM_NUM_1, TE_ALIGN_NO, NULL, NULL, &te_loadFxn_mtxinv, &te_processFxn_matinv };
static const tTestEngDesc descr_cmtx_inv4x4_32x32 =  { FMT_CPLX | FMT_FLOAT32, MTXINV_PLAIN,&cmtx_inv4x4_32x32_Api  ,TE_DIM_NUM_1, TE_ALIGN_NO, NULL, NULL, &te_loadFxn_mtxinv, &te_processFxn_matinv };
static const tTestEngDesc descr_cmtx_inv6x6_32x32 =  { FMT_CPLX | FMT_FLOAT32, MTXINV_PLAIN,&cmtx_inv6x4_32x32_Api  ,TE_DIM_NUM_1, TE_ALIGN_NO, NULL, NULL, &te_loadFxn_mtxinv, &te_processFxn_matinv };
static const tTestEngDesc descr_cmtx_inv8x8_32x32 =  { FMT_CPLX | FMT_FLOAT32, MTXINV_PLAIN,&cmtx_inv8x4_32x32_Api  ,TE_DIM_NUM_1, TE_ALIGN_NO, NULL, NULL, &te_loadFxn_mtxinv, &te_processFxn_matinv };
static const tTestEngDesc descr_cmtx_inv10x10_32x32 ={ FMT_CPLX | FMT_FLOAT32, MTXINV_PLAIN,&cmtx_inv10x10_32x32_Api,TE_DIM_NUM_1, TE_ALIGN_NO, NULL, NULL, &te_loadFxn_mtxinv, &te_processFxn_matinv };

static const tTestMtxinv tests_32x32[] = 
{
  { &descr_cmtx_inv2x2_32x32,   (tTestEngTarget)&cmtx_inv2x2_32x32,   1,"gj1/cmtx_inv2x2_32x32_cond10.seq"    },
  { &descr_cmtx_inv2x2_32x32,   (tTestEngTarget)&cmtx_inv2x2_32x32,   1,"gj1/cmtx_inv2x2_32x32_cond100.seq"   },
  { &descr_cmtx_inv2x2_32x32,   (tTestEngTarget)&cmtx_inv2x2_32x32,   1,"gj1/cmtx_inv2x2_32x32_cond1000.seq"  },
  { &descr_cmtx_inv3x3_32x32,   (tTestEngTarget)&cmtx_inv3x3_32x32,   1,"gj1/cmtx_inv3x3_32x32_cond10.seq"    },
  { &descr_cmtx_inv3x3_32x32,   (tTestEngTarget)&cmtx_inv3x3_32x32,   1,"gj1/cmtx_inv3x3_32x32_cond100.seq"   },
  { &descr_cmtx_inv3x3_32x32,   (tTestEngTarget)&cmtx_inv3x3_32x32,   1,"gj1/cmtx_inv3x3_32x32_cond1000.seq"  },
  { &descr_cmtx_inv4x4_32x32,   (tTestEngTarget)&cmtx_inv4x4_32x32,   1,"gj1/cmtx_inv4x4_32x32_cond10.seq"    },
  { &descr_cmtx_inv4x4_32x32,   (tTestEngTarget)&cmtx_inv4x4_32x32,   1,"gj1/cmtx_inv4x4_32x32_cond100.seq"   },
  { &descr_cmtx_inv4x4_32x32,   (tTestEngTarget)&cmtx_inv4x4_32x32,   1,"gj1/cmtx_inv4x4_32x32_cond1000.seq"  },
  { &descr_cmtx_inv6x6_32x32,   (tTestEngTarget)&cmtx_inv6x6_32x32,   1,"gj1/cmtx_inv6x6_32x32_cond10.seq"    },
  { &descr_cmtx_inv6x6_32x32,   (tTestEngTarget)&cmtx_inv6x6_32x32,   1,"gj1/cmtx_inv6x6_32x32_cond100.seq"   },
  { &descr_cmtx_inv6x6_32x32,   (tTestEngTarget)&cmtx_inv6x6_32x32,   1,"gj1/cmtx_inv6x6_32x32_cond1000.seq"  },
  { &descr_cmtx_inv8x8_32x32,   (tTestEngTarget)&cmtx_inv8x8_32x32,   1,"gj1/cmtx_inv8x8_32x32_cond10.seq"    },
  { &descr_cmtx_inv8x8_32x32,   (tTestEngTarget)&cmtx_inv8x8_32x32,   1,"gj1/cmtx_inv8x8_32x32_cond100.seq"   },
  { &descr_cmtx_inv8x8_32x32,   (tTestEngTarget)&cmtx_inv8x8_32x32,   1,"gj1/cmtx_inv8x8_32x32_cond1000.seq"  },
  { &descr_cmtx_inv10x10_32x32, (tTestEngTarget)&cmtx_inv10x10_32x32, 1,"gj1/cmtx_inv10x10_32x32_cond10.seq"    },
  { &descr_cmtx_inv10x10_32x32, (tTestEngTarget)&cmtx_inv10x10_32x32, 1,"gj1/cmtx_inv10x10_32x32_cond100.seq"   },
  { &descr_cmtx_inv10x10_32x32, (tTestEngTarget)&cmtx_inv10x10_32x32, 1,"gj1/cmtx_inv10x10_32x32_cond1000.seq"  },

  { &descr_mtx_inv2x2_32x32,   (tTestEngTarget)&mtx_inv2x2_32x32,   1,"gj1/mtx_inv2x2_32x32_cond10.seq"    },
  { &descr_mtx_inv2x2_32x32,   (tTestEngTarget)&mtx_inv2x2_32x32,   1,"gj1/mtx_inv2x2_32x32_cond100.seq"   },
  { &descr_mtx_inv2x2_32x32,   (tTestEngTarget)&mtx_inv2x2_32x32,   1,"gj1/mtx_inv2x2_32x32_cond1000.seq"  },
  { &descr_mtx_inv3x3_32x32,   (tTestEngTarget)&mtx_inv3x3_32x32,   1,"gj1/mtx_inv3x3_32x32_cond10.seq"    },
  { &descr_mtx_inv3x3_32x32,   (tTestEngTarget)&mtx_inv3x3_32x32,   1,"gj1/mtx_inv3x3_32x32_cond100.seq"   },
  { &descr_mtx_inv3x3_32x32,   (tTestEngTarget)&mtx_inv3x3_32x32,   1,"gj1/mtx_inv3x3_32x32_cond1000.seq"  },
  { &descr_mtx_inv4x4_32x32,   (tTestEngTarget)&mtx_inv4x4_32x32,   1,"gj1/mtx_inv4x4_32x32_cond10.seq"    },
  { &descr_mtx_inv4x4_32x32,   (tTestEngTarget)&mtx_inv4x4_32x32,   1,"gj1/mtx_inv4x4_32x32_cond100.seq"   },
  { &descr_mtx_inv4x4_32x32,   (tTestEngTarget)&mtx_inv4x4_32x32,   1,"gj1/mtx_inv4x4_32x32_cond1000.seq"  },
  { &descr_mtx_inv6x6_32x32,   (tTestEngTarget)&mtx_inv6x6_32x32,   1,"gj1/mtx_inv6x6_32x32_cond10.seq"    },
  { &descr_mtx_inv6x6_32x32,   (tTestEngTarget)&mtx_inv6x6_32x32,   1,"gj1/mtx_inv6x6_32x32_cond100.seq"   },
  { &descr_mtx_inv6x6_32x32,   (tTestEngTarget)&mtx_inv6x6_32x32,   1,"gj1/mtx_inv6x6_32x32_cond1000.seq"  },
  { &descr_mtx_inv8x8_32x32,   (tTestEngTarget)&mtx_inv8x8_32x32,   1,"gj1/mtx_inv8x8_32x32_cond10.seq"    },
  { &descr_mtx_inv8x8_32x32,   (tTestEngTarget)&mtx_inv8x8_32x32,   1,"gj1/mtx_inv8x8_32x32_cond100.seq"   },
  { &descr_mtx_inv8x8_32x32,   (tTestEngTarget)&mtx_inv8x8_32x32,   1,"gj1/mtx_inv8x8_32x32_cond1000.seq"  },
  { &descr_mtx_inv10x10_32x32, (tTestEngTarget)&mtx_inv10x10_32x32, 1,"gj1/mtx_inv10x10_32x32_cond10.seq"    },
  { &descr_mtx_inv10x10_32x32, (tTestEngTarget)&mtx_inv10x10_32x32, 1,"gj1/mtx_inv10x10_32x32_cond100.seq"   },
  { &descr_mtx_inv10x10_32x32, (tTestEngTarget)&mtx_inv10x10_32x32, 1,"gj1/mtx_inv10x10_32x32_cond1000.seq"  },

  { &descr_cmtx_gjelim2x2_32x32,   (tTestEngTarget)&cmtx_gjelim2x2_32x32,   1,"gj1/cmtx_gjelim2x2_32x32_cond10.seq"    },
  { &descr_cmtx_gjelim2x2_32x32,   (tTestEngTarget)&cmtx_gjelim2x2_32x32,   1,"gj1/cmtx_gjelim2x2_32x32_cond100.seq"   },
  { &descr_cmtx_gjelim2x2_32x32,   (tTestEngTarget)&cmtx_gjelim2x2_32x32,   1,"gj1/cmtx_gjelim2x2_32x32_cond1000.seq"  },
  { &descr_cmtx_gjelim3x3_32x32,   (tTestEngTarget)&cmtx_gjelim3x3_32x32,   1,"gj1/cmtx_gjelim3x3_32x32_cond10.seq"    },
  { &descr_cmtx_gjelim3x3_32x32,   (tTestEngTarget)&cmtx_gjelim3x3_32x32,   1,"gj1/cmtx_gjelim3x3_32x32_cond100.seq"   },
  { &descr_cmtx_gjelim3x3_32x32,   (tTestEngTarget)&cmtx_gjelim3x3_32x32,   1,"gj1/cmtx_gjelim3x3_32x32_cond1000.seq"  },
  { &descr_cmtx_gjelim4x4_32x32,   (tTestEngTarget)&cmtx_gjelim4x4_32x32,   1,"gj1/cmtx_gjelim4x4_32x32_cond10.seq"    },
  { &descr_cmtx_gjelim4x4_32x32,   (tTestEngTarget)&cmtx_gjelim4x4_32x32,   1,"gj1/cmtx_gjelim4x4_32x32_cond100.seq"   },
  { &descr_cmtx_gjelim4x4_32x32,   (tTestEngTarget)&cmtx_gjelim4x4_32x32,   1,"gj1/cmtx_gjelim4x4_32x32_cond1000.seq"  },
  { &descr_cmtx_gjelim6x6_32x32,   (tTestEngTarget)&cmtx_gjelim6x6_32x32,   1,"gj1/cmtx_gjelim6x6_32x32_cond10.seq"    },
  { &descr_cmtx_gjelim6x6_32x32,   (tTestEngTarget)&cmtx_gjelim6x6_32x32,   1,"gj1/cmtx_gjelim6x6_32x32_cond100.seq"   },
  { &descr_cmtx_gjelim6x6_32x32,   (tTestEngTarget)&cmtx_gjelim6x6_32x32,   1,"gj1/cmtx_gjelim6x6_32x32_cond1000.seq"  },
  { &descr_cmtx_gjelim8x8_32x32,   (tTestEngTarget)&cmtx_gjelim8x8_32x32,   1,"gj1/cmtx_gjelim8x8_32x32_cond10.seq"    },
  { &descr_cmtx_gjelim8x8_32x32,   (tTestEngTarget)&cmtx_gjelim8x8_32x32,   1,"gj1/cmtx_gjelim8x8_32x32_cond100.seq"   },
  { &descr_cmtx_gjelim8x8_32x32,   (tTestEngTarget)&cmtx_gjelim8x8_32x32,   1,"gj1/cmtx_gjelim8x8_32x32_cond1000.seq"  },
  { &descr_cmtx_gjelim10x10_32x32, (tTestEngTarget)&cmtx_gjelim10x10_32x32, 1,"gj1/cmtx_gjelim10x10_32x32_cond10.seq"    },
  { &descr_cmtx_gjelim10x10_32x32, (tTestEngTarget)&cmtx_gjelim10x10_32x32, 1,"gj1/cmtx_gjelim10x10_32x32_cond100.seq"   },
  { &descr_cmtx_gjelim10x10_32x32, (tTestEngTarget)&cmtx_gjelim10x10_32x32, 1,"gj1/cmtx_gjelim10x10_32x32_cond1000.seq"  },

  { &descr_mtx_gjelim2x2_32x32,   (tTestEngTarget)&mtx_gjelim2x2_32x32,   1,"gj1/mtx_gjelim2x2_32x32_cond10.seq"    },
  { &descr_mtx_gjelim2x2_32x32,   (tTestEngTarget)&mtx_gjelim2x2_32x32,   1,"gj1/mtx_gjelim2x2_32x32_cond100.seq"   },
  { &descr_mtx_gjelim2x2_32x32,   (tTestEngTarget)&mtx_gjelim2x2_32x32,   1,"gj1/mtx_gjelim2x2_32x32_cond1000.seq"  },
  { &descr_mtx_gjelim3x3_32x32,   (tTestEngTarget)&mtx_gjelim3x3_32x32,   1,"gj1/mtx_gjelim3x3_32x32_cond10.seq"    },
  { &descr_mtx_gjelim3x3_32x32,   (tTestEngTarget)&mtx_gjelim3x3_32x32,   1,"gj1/mtx_gjelim3x3_32x32_cond100.seq"   },
  { &descr_mtx_gjelim3x3_32x32,   (tTestEngTarget)&mtx_gjelim3x3_32x32,   1,"gj1/mtx_gjelim3x3_32x32_cond1000.seq"  },
  { &descr_mtx_gjelim4x4_32x32,   (tTestEngTarget)&mtx_gjelim4x4_32x32,   1,"gj1/mtx_gjelim4x4_32x32_cond10.seq"    },
  { &descr_mtx_gjelim4x4_32x32,   (tTestEngTarget)&mtx_gjelim4x4_32x32,   1,"gj1/mtx_gjelim4x4_32x32_cond100.seq"   },
  { &descr_mtx_gjelim4x4_32x32,   (tTestEngTarget)&mtx_gjelim4x4_32x32,   1,"gj1/mtx_gjelim4x4_32x32_cond1000.seq"  },
  { &descr_mtx_gjelim6x6_32x32,   (tTestEngTarget)&mtx_gjelim6x6_32x32,   1,"gj1/mtx_gjelim6x6_32x32_cond10.seq"    },
  { &descr_mtx_gjelim6x6_32x32,   (tTestEngTarget)&mtx_gjelim6x6_32x32,   1,"gj1/mtx_gjelim6x6_32x32_cond100.seq"   },
  { &descr_mtx_gjelim6x6_32x32,   (tTestEngTarget)&mtx_gjelim6x6_32x32,   1,"gj1/mtx_gjelim6x6_32x32_cond1000.seq"  },
  { &descr_mtx_gjelim8x8_32x32,   (tTestEngTarget)&mtx_gjelim8x8_32x32,   1,"gj1/mtx_gjelim8x8_32x32_cond10.seq"    },
  { &descr_mtx_gjelim8x8_32x32,   (tTestEngTarget)&mtx_gjelim8x8_32x32,   1,"gj1/mtx_gjelim8x8_32x32_cond100.seq"   },
  { &descr_mtx_gjelim8x8_32x32,   (tTestEngTarget)&mtx_gjelim8x8_32x32,   1,"gj1/mtx_gjelim8x8_32x32_cond1000.seq"  },
  { &descr_mtx_gjelim10x10_32x32, (tTestEngTarget)&mtx_gjelim10x10_32x32, 1,"gj1/mtx_gjelim10x10_32x32_cond10.seq"    },
  { &descr_mtx_gjelim10x10_32x32, (tTestEngTarget)&mtx_gjelim10x10_32x32, 1,"gj1/mtx_gjelim10x10_32x32_cond100.seq"   },
  { &descr_mtx_gjelim10x10_32x32, (tTestEngTarget)&mtx_gjelim10x10_32x32, 1,"gj1/mtx_gjelim10x10_32x32_cond1000.seq"  },
  
  { NULL }  /* end of list */ 
};

/* Perform all tests for matrix inversion API functions. */
int func_gj1(int isFull, int isVerbose, int breakOnError)
{
  int res = 1;
  res &= (0!=te_ExecMtxInv(tests_32x32,isFull,isVerbose,breakOnError));
  return res;
}
