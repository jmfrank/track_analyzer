/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * diff_gradients_data.cpp
 *
 * Code generation for function 'diff_gradients_data'
 *
 */

/* Include files */
#include "diff_gradients_data.h"
#include "diff_gradients.h"
#include "rt_nonfinite.h"

/* Variable Definitions */
emlrtCTX emlrtRootTLSGlobal = NULL;
const volatile char_T *emlrtBreakCheckR2012bFlagVar = NULL;
emlrtContext emlrtContextGlobal = { true,/* bFirstTime */
  false,                               /* bInitialized */
  131483U,                             /* fVersionInfo */
  NULL,                                /* fErrorFunction */
  "diff_gradients",                    /* fFunctionName */
  NULL,                                /* fRTCallStack */
  false,                               /* bDebugMode */
  { 1858410525U, 2505464270U, 328108647U, 1256672073U },/* fSigWrd */
  NULL                                 /* fSigMem */
};

emlrtRSInfo q_emlrtRSI = { 76,         /* lineNo */
  "validateattributes",                /* fcnName */
  "/usr/local/MATLAB/R2019b/toolbox/eml/lib/matlab/lang/validateattributes.m"/* pathName */
};

emlrtRSInfo ab_emlrtRSI = { 15,        /* lineNo */
  "sum",                               /* fcnName */
  "/usr/local/MATLAB/R2019b/toolbox/eml/lib/matlab/datafun/sum.m"/* pathName */
};

emlrtRSInfo bb_emlrtRSI = { 96,        /* lineNo */
  "sumprod",                           /* fcnName */
  "/usr/local/MATLAB/R2019b/toolbox/eml/lib/matlab/datafun/private/sumprod.m"/* pathName */
};

emlrtRSInfo cb_emlrtRSI = { 124,       /* lineNo */
  "combineVectorElements",             /* fcnName */
  "/usr/local/MATLAB/R2019b/toolbox/eml/lib/matlab/datafun/private/combineVectorElements.m"/* pathName */
};

emlrtRSInfo db_emlrtRSI = { 184,       /* lineNo */
  "colMajorFlatIter",                  /* fcnName */
  "/usr/local/MATLAB/R2019b/toolbox/eml/lib/matlab/datafun/private/combineVectorElements.m"/* pathName */
};

emlrtRSInfo eb_emlrtRSI = { 21,        /* lineNo */
  "eml_int_forloop_overflow_check",    /* fcnName */
  "/usr/local/MATLAB/R2019b/toolbox/eml/lib/matlab/eml/eml_int_forloop_overflow_check.m"/* pathName */
};

emlrtRSInfo ib_emlrtRSI = { 14,        /* lineNo */
  "max",                               /* fcnName */
  "/usr/local/MATLAB/R2019b/toolbox/eml/lib/matlab/datafun/max.m"/* pathName */
};

emlrtRSInfo jb_emlrtRSI = { 20,        /* lineNo */
  "minOrMax",                          /* fcnName */
  "/usr/local/MATLAB/R2019b/toolbox/eml/eml/+coder/+internal/minOrMax.m"/* pathName */
};

emlrtRSInfo kb_emlrtRSI = { 45,        /* lineNo */
  "unaryOrBinaryDispatch",             /* fcnName */
  "/usr/local/MATLAB/R2019b/toolbox/eml/eml/+coder/+internal/minOrMax.m"/* pathName */
};

emlrtRSInfo lb_emlrtRSI = { 145,       /* lineNo */
  "unaryMinOrMax",                     /* fcnName */
  "/usr/local/MATLAB/R2019b/toolbox/eml/eml/+coder/+internal/unaryMinOrMax.m"/* pathName */
};

emlrtRSInfo mb_emlrtRSI = { 1019,      /* lineNo */
  "maxRealVectorOmitNaN",              /* fcnName */
  "/usr/local/MATLAB/R2019b/toolbox/eml/eml/+coder/+internal/unaryMinOrMax.m"/* pathName */
};

emlrtRSInfo nb_emlrtRSI = { 932,       /* lineNo */
  "minOrMaxRealVector",                /* fcnName */
  "/usr/local/MATLAB/R2019b/toolbox/eml/eml/+coder/+internal/unaryMinOrMax.m"/* pathName */
};

emlrtRSInfo ob_emlrtRSI = { 924,       /* lineNo */
  "minOrMaxRealVector",                /* fcnName */
  "/usr/local/MATLAB/R2019b/toolbox/eml/eml/+coder/+internal/unaryMinOrMax.m"/* pathName */
};

emlrtRSInfo pb_emlrtRSI = { 975,       /* lineNo */
  "findFirst",                         /* fcnName */
  "/usr/local/MATLAB/R2019b/toolbox/eml/eml/+coder/+internal/unaryMinOrMax.m"/* pathName */
};

emlrtRSInfo qb_emlrtRSI = { 992,       /* lineNo */
  "minOrMaxRealVectorKernel",          /* fcnName */
  "/usr/local/MATLAB/R2019b/toolbox/eml/eml/+coder/+internal/unaryMinOrMax.m"/* pathName */
};

emlrtRSInfo rb_emlrtRSI = { 20,        /* lineNo */
  "padarray",                          /* fcnName */
  "/usr/local/MATLAB/R2019b/toolbox/images/images/eml/padarray.m"/* pathName */
};

emlrtRSInfo sb_emlrtRSI = { 66,        /* lineNo */
  "padarray",                          /* fcnName */
  "/usr/local/MATLAB/R2019b/toolbox/images/images/eml/padarray.m"/* pathName */
};

emlrtRSInfo dc_emlrtRSI = { 189,       /* lineNo */
  "lincombSharedLibrary",              /* fcnName */
  "/usr/local/MATLAB/R2019b/toolbox/images/images/eml/imlincomb.m"/* pathName */
};

emlrtRSInfo fc_emlrtRSI = { 55,        /* lineNo */
  "power",                             /* fcnName */
  "/usr/local/MATLAB/R2019b/toolbox/eml/lib/matlab/ops/power.m"/* pathName */
};

emlrtRSInfo wd_emlrtRSI = { 106,       /* lineNo */
  "imfilter",                          /* fcnName */
  "/usr/local/MATLAB/R2019b/toolbox/images/images/eml/imfilter.m"/* pathName */
};

emlrtRSInfo xd_emlrtRSI = { 110,       /* lineNo */
  "imfilter",                          /* fcnName */
  "/usr/local/MATLAB/R2019b/toolbox/images/images/eml/imfilter.m"/* pathName */
};

emlrtRSInfo yd_emlrtRSI = { 773,       /* lineNo */
  "padImage",                          /* fcnName */
  "/usr/local/MATLAB/R2019b/toolbox/images/images/eml/imfilter.m"/* pathName */
};

emlrtRSInfo ae_emlrtRSI = { 843,       /* lineNo */
  "filterPartOrWhole",                 /* fcnName */
  "/usr/local/MATLAB/R2019b/toolbox/images/images/eml/imfilter.m"/* pathName */
};

emlrtRSInfo be_emlrtRSI = { 917,       /* lineNo */
  "imfiltercore",                      /* fcnName */
  "/usr/local/MATLAB/R2019b/toolbox/images/images/eml/imfilter.m"/* pathName */
};

emlrtRSInfo fe_emlrtRSI = { 66,        /* lineNo */
  "applyBinaryScalarFunction",         /* fcnName */
  "/usr/local/MATLAB/R2019b/toolbox/eml/eml/+coder/+internal/applyBinaryScalarFunction.m"/* pathName */
};

emlrtRSInfo ge_emlrtRSI = { 188,       /* lineNo */
  "flatIter",                          /* fcnName */
  "/usr/local/MATLAB/R2019b/toolbox/eml/eml/+coder/+internal/applyBinaryScalarFunction.m"/* pathName */
};

emlrtRSInfo ie_emlrtRSI = { 33,        /* lineNo */
  "applyScalarFunctionInPlace",        /* fcnName */
  "/usr/local/MATLAB/R2019b/toolbox/eml/eml/+coder/+internal/applyScalarFunctionInPlace.m"/* pathName */
};

emlrtRSInfo xe_emlrtRSI = { 32,        /* lineNo */
  "gauss",                             /* fcnName */
  "/media/data/Dropbox/code_bank/matlab/3Dnucleisegmentation/gaussgradient.m"/* pathName */
};

emlrtRSInfo gg_emlrtRSI = { 9,         /* lineNo */
  "int",                               /* fcnName */
  "/usr/local/MATLAB/R2019b/toolbox/eml/eml/+coder/+internal/+lapack/int.m"/* pathName */
};

emlrtRSInfo hg_emlrtRSI = { 8,         /* lineNo */
  "majority",                          /* fcnName */
  "/usr/local/MATLAB/R2019b/toolbox/eml/eml/+coder/+internal/+lapack/majority.m"/* pathName */
};

emlrtRSInfo ig_emlrtRSI = { 31,        /* lineNo */
  "infocheck",                         /* fcnName */
  "/usr/local/MATLAB/R2019b/toolbox/eml/eml/+coder/+internal/+lapack/infocheck.m"/* pathName */
};

emlrtRSInfo lg_emlrtRSI = { 174,       /* lineNo */
  "ceval_xgesvd",                      /* fcnName */
  "/usr/local/MATLAB/R2019b/toolbox/eml/eml/+coder/+internal/+lapack/xgesvd.m"/* pathName */
};

emlrtRSInfo mg_emlrtRSI = { 117,       /* lineNo */
  "ceval_xgesvd",                      /* fcnName */
  "/usr/local/MATLAB/R2019b/toolbox/eml/eml/+coder/+internal/+lapack/xgesvd.m"/* pathName */
};

emlrtRSInfo ng_emlrtRSI = { 110,       /* lineNo */
  "ceval_xgesvd",                      /* fcnName */
  "/usr/local/MATLAB/R2019b/toolbox/eml/eml/+coder/+internal/+lapack/xgesvd.m"/* pathName */
};

emlrtRSInfo og_emlrtRSI = { 59,        /* lineNo */
  "ceval_xgesvd",                      /* fcnName */
  "/usr/local/MATLAB/R2019b/toolbox/eml/eml/+coder/+internal/+lapack/xgesvd.m"/* pathName */
};

emlrtRTEInfo c_emlrtRTEI = { 13,       /* lineNo */
  37,                                  /* colNo */
  "validateinteger",                   /* fName */
  "/usr/local/MATLAB/R2019b/toolbox/eml/eml/+coder/+internal/+valattr/validateinteger.m"/* pName */
};

emlrtRTEInfo l_emlrtRTEI = { 26,       /* lineNo */
  27,                                  /* colNo */
  "unaryMinOrMax",                     /* fName */
  "/usr/local/MATLAB/R2019b/toolbox/eml/eml/+coder/+internal/unaryMinOrMax.m"/* pName */
};

emlrtRTEInfo m_emlrtRTEI = { 95,       /* lineNo */
  27,                                  /* colNo */
  "unaryMinOrMax",                     /* fName */
  "/usr/local/MATLAB/R2019b/toolbox/eml/eml/+coder/+internal/unaryMinOrMax.m"/* pName */
};

emlrtRTEInfo q_emlrtRTEI = { 13,       /* lineNo */
  9,                                   /* colNo */
  "sqrt",                              /* fName */
  "/usr/local/MATLAB/R2019b/toolbox/eml/lib/matlab/elfun/sqrt.m"/* pName */
};

emlrtRTEInfo t_emlrtRTEI = { 47,       /* lineNo */
  19,                                  /* colNo */
  "allOrAny",                          /* fName */
  "/usr/local/MATLAB/R2019b/toolbox/eml/eml/+coder/+internal/allOrAny.m"/* pName */
};

emlrtRTEInfo yd_emlrtRTEI = { 917,     /* lineNo */
  11,                                  /* colNo */
  "imfilter",                          /* fName */
  "/usr/local/MATLAB/R2019b/toolbox/images/images/eml/imfilter.m"/* pName */
};

emlrtRTEInfo vg_emlrtRTEI = { 59,      /* lineNo */
  9,                                   /* colNo */
  "imfilter",                          /* fName */
  "/usr/local/MATLAB/R2019b/toolbox/images/images/eml/imfilter.m"/* pName */
};

/* End of code generation (diff_gradients_data.cpp) */
