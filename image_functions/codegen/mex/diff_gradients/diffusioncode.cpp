/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * diffusioncode.cpp
 *
 * Code generation for function 'diffusioncode'
 *
 */

/* Include files */
#include "diffusioncode.h"
#include "diff_gradients.h"
#include "diff_gradients_data.h"
#include "diff_gradients_emxutil.h"
#include "exp.h"
#include "imfilter.h"
#include "libmwimfilter.h"
#include "libmwippfilter.h"
#include "padarray.h"
#include "power.h"
#include "rt_nonfinite.h"

/* Variable Definitions */
static emlrtRSInfo gc_emlrtRSI = { 82, /* lineNo */
  "diffusioncode",                     /* fcnName */
  "/media/data/Dropbox/code_bank/matlab/3Dnucleisegmentation/diffusioncode.m"/* pathName */
};

static emlrtRSInfo hc_emlrtRSI = { 83, /* lineNo */
  "diffusioncode",                     /* fcnName */
  "/media/data/Dropbox/code_bank/matlab/3Dnucleisegmentation/diffusioncode.m"/* pathName */
};

static emlrtRSInfo ic_emlrtRSI = { 84, /* lineNo */
  "diffusioncode",                     /* fcnName */
  "/media/data/Dropbox/code_bank/matlab/3Dnucleisegmentation/diffusioncode.m"/* pathName */
};

static emlrtRSInfo jc_emlrtRSI = { 85, /* lineNo */
  "diffusioncode",                     /* fcnName */
  "/media/data/Dropbox/code_bank/matlab/3Dnucleisegmentation/diffusioncode.m"/* pathName */
};

static emlrtRSInfo kc_emlrtRSI = { 86, /* lineNo */
  "diffusioncode",                     /* fcnName */
  "/media/data/Dropbox/code_bank/matlab/3Dnucleisegmentation/diffusioncode.m"/* pathName */
};

static emlrtRSInfo lc_emlrtRSI = { 87, /* lineNo */
  "diffusioncode",                     /* fcnName */
  "/media/data/Dropbox/code_bank/matlab/3Dnucleisegmentation/diffusioncode.m"/* pathName */
};

static emlrtRSInfo mc_emlrtRSI = { 88, /* lineNo */
  "diffusioncode",                     /* fcnName */
  "/media/data/Dropbox/code_bank/matlab/3Dnucleisegmentation/diffusioncode.m"/* pathName */
};

static emlrtRSInfo nc_emlrtRSI = { 89, /* lineNo */
  "diffusioncode",                     /* fcnName */
  "/media/data/Dropbox/code_bank/matlab/3Dnucleisegmentation/diffusioncode.m"/* pathName */
};

static emlrtRSInfo oc_emlrtRSI = { 94, /* lineNo */
  "diffusioncode",                     /* fcnName */
  "/media/data/Dropbox/code_bank/matlab/3Dnucleisegmentation/diffusioncode.m"/* pathName */
};

static emlrtRSInfo pc_emlrtRSI = { 95, /* lineNo */
  "diffusioncode",                     /* fcnName */
  "/media/data/Dropbox/code_bank/matlab/3Dnucleisegmentation/diffusioncode.m"/* pathName */
};

static emlrtRSInfo qc_emlrtRSI = { 96, /* lineNo */
  "diffusioncode",                     /* fcnName */
  "/media/data/Dropbox/code_bank/matlab/3Dnucleisegmentation/diffusioncode.m"/* pathName */
};

static emlrtRSInfo rc_emlrtRSI = { 97, /* lineNo */
  "diffusioncode",                     /* fcnName */
  "/media/data/Dropbox/code_bank/matlab/3Dnucleisegmentation/diffusioncode.m"/* pathName */
};

static emlrtRSInfo sc_emlrtRSI = { 98, /* lineNo */
  "diffusioncode",                     /* fcnName */
  "/media/data/Dropbox/code_bank/matlab/3Dnucleisegmentation/diffusioncode.m"/* pathName */
};

static emlrtRSInfo tc_emlrtRSI = { 99, /* lineNo */
  "diffusioncode",                     /* fcnName */
  "/media/data/Dropbox/code_bank/matlab/3Dnucleisegmentation/diffusioncode.m"/* pathName */
};

static emlrtRSInfo uc_emlrtRSI = { 100,/* lineNo */
  "diffusioncode",                     /* fcnName */
  "/media/data/Dropbox/code_bank/matlab/3Dnucleisegmentation/diffusioncode.m"/* pathName */
};

static emlrtRSInfo vc_emlrtRSI = { 101,/* lineNo */
  "diffusioncode",                     /* fcnName */
  "/media/data/Dropbox/code_bank/matlab/3Dnucleisegmentation/diffusioncode.m"/* pathName */
};

static emlrtRSInfo wc_emlrtRSI = { 104,/* lineNo */
  "diffusioncode",                     /* fcnName */
  "/media/data/Dropbox/code_bank/matlab/3Dnucleisegmentation/diffusioncode.m"/* pathName */
};

static emlrtRSInfo xc_emlrtRSI = { 105,/* lineNo */
  "diffusioncode",                     /* fcnName */
  "/media/data/Dropbox/code_bank/matlab/3Dnucleisegmentation/diffusioncode.m"/* pathName */
};

static emlrtRSInfo yc_emlrtRSI = { 106,/* lineNo */
  "diffusioncode",                     /* fcnName */
  "/media/data/Dropbox/code_bank/matlab/3Dnucleisegmentation/diffusioncode.m"/* pathName */
};

static emlrtRSInfo ad_emlrtRSI = { 107,/* lineNo */
  "diffusioncode",                     /* fcnName */
  "/media/data/Dropbox/code_bank/matlab/3Dnucleisegmentation/diffusioncode.m"/* pathName */
};

static emlrtRSInfo bd_emlrtRSI = { 108,/* lineNo */
  "diffusioncode",                     /* fcnName */
  "/media/data/Dropbox/code_bank/matlab/3Dnucleisegmentation/diffusioncode.m"/* pathName */
};

static emlrtRSInfo cd_emlrtRSI = { 109,/* lineNo */
  "diffusioncode",                     /* fcnName */
  "/media/data/Dropbox/code_bank/matlab/3Dnucleisegmentation/diffusioncode.m"/* pathName */
};

static emlrtRSInfo dd_emlrtRSI = { 110,/* lineNo */
  "diffusioncode",                     /* fcnName */
  "/media/data/Dropbox/code_bank/matlab/3Dnucleisegmentation/diffusioncode.m"/* pathName */
};

static emlrtRSInfo ed_emlrtRSI = { 111,/* lineNo */
  "diffusioncode",                     /* fcnName */
  "/media/data/Dropbox/code_bank/matlab/3Dnucleisegmentation/diffusioncode.m"/* pathName */
};

static emlrtRSInfo fd_emlrtRSI = { 113,/* lineNo */
  "diffusioncode",                     /* fcnName */
  "/media/data/Dropbox/code_bank/matlab/3Dnucleisegmentation/diffusioncode.m"/* pathName */
};

static emlrtRSInfo gd_emlrtRSI = { 114,/* lineNo */
  "diffusioncode",                     /* fcnName */
  "/media/data/Dropbox/code_bank/matlab/3Dnucleisegmentation/diffusioncode.m"/* pathName */
};

static emlrtRSInfo hd_emlrtRSI = { 115,/* lineNo */
  "diffusioncode",                     /* fcnName */
  "/media/data/Dropbox/code_bank/matlab/3Dnucleisegmentation/diffusioncode.m"/* pathName */
};

static emlrtRSInfo id_emlrtRSI = { 116,/* lineNo */
  "diffusioncode",                     /* fcnName */
  "/media/data/Dropbox/code_bank/matlab/3Dnucleisegmentation/diffusioncode.m"/* pathName */
};

static emlrtRSInfo jd_emlrtRSI = { 117,/* lineNo */
  "diffusioncode",                     /* fcnName */
  "/media/data/Dropbox/code_bank/matlab/3Dnucleisegmentation/diffusioncode.m"/* pathName */
};

static emlrtRSInfo kd_emlrtRSI = { 118,/* lineNo */
  "diffusioncode",                     /* fcnName */
  "/media/data/Dropbox/code_bank/matlab/3Dnucleisegmentation/diffusioncode.m"/* pathName */
};

static emlrtRSInfo ld_emlrtRSI = { 119,/* lineNo */
  "diffusioncode",                     /* fcnName */
  "/media/data/Dropbox/code_bank/matlab/3Dnucleisegmentation/diffusioncode.m"/* pathName */
};

static emlrtRSInfo md_emlrtRSI = { 120,/* lineNo */
  "diffusioncode",                     /* fcnName */
  "/media/data/Dropbox/code_bank/matlab/3Dnucleisegmentation/diffusioncode.m"/* pathName */
};

static emlrtRSInfo nd_emlrtRSI = { 122,/* lineNo */
  "diffusioncode",                     /* fcnName */
  "/media/data/Dropbox/code_bank/matlab/3Dnucleisegmentation/diffusioncode.m"/* pathName */
};

static emlrtRSInfo od_emlrtRSI = { 123,/* lineNo */
  "diffusioncode",                     /* fcnName */
  "/media/data/Dropbox/code_bank/matlab/3Dnucleisegmentation/diffusioncode.m"/* pathName */
};

static emlrtRSInfo pd_emlrtRSI = { 124,/* lineNo */
  "diffusioncode",                     /* fcnName */
  "/media/data/Dropbox/code_bank/matlab/3Dnucleisegmentation/diffusioncode.m"/* pathName */
};

static emlrtRSInfo qd_emlrtRSI = { 125,/* lineNo */
  "diffusioncode",                     /* fcnName */
  "/media/data/Dropbox/code_bank/matlab/3Dnucleisegmentation/diffusioncode.m"/* pathName */
};

static emlrtRSInfo rd_emlrtRSI = { 126,/* lineNo */
  "diffusioncode",                     /* fcnName */
  "/media/data/Dropbox/code_bank/matlab/3Dnucleisegmentation/diffusioncode.m"/* pathName */
};

static emlrtRSInfo sd_emlrtRSI = { 127,/* lineNo */
  "diffusioncode",                     /* fcnName */
  "/media/data/Dropbox/code_bank/matlab/3Dnucleisegmentation/diffusioncode.m"/* pathName */
};

static emlrtRSInfo td_emlrtRSI = { 128,/* lineNo */
  "diffusioncode",                     /* fcnName */
  "/media/data/Dropbox/code_bank/matlab/3Dnucleisegmentation/diffusioncode.m"/* pathName */
};

static emlrtRSInfo ud_emlrtRSI = { 129,/* lineNo */
  "diffusioncode",                     /* fcnName */
  "/media/data/Dropbox/code_bank/matlab/3Dnucleisegmentation/diffusioncode.m"/* pathName */
};

static emlrtRSInfo vd_emlrtRSI = { 131,/* lineNo */
  "diffusioncode",                     /* fcnName */
  "/media/data/Dropbox/code_bank/matlab/3Dnucleisegmentation/diffusioncode.m"/* pathName */
};

static emlrtMCInfo emlrtMCI = { 27,    /* lineNo */
  5,                                   /* colNo */
  "error",                             /* fName */
  "/usr/local/MATLAB/R2019b/toolbox/eml/lib/matlab/lang/error.m"/* pName */
};

static emlrtRTEInfo cb_emlrtRTEI = { 79,/* lineNo */
  13,                                  /* colNo */
  "diffusioncode",                     /* fName */
  "/media/data/Dropbox/code_bank/matlab/3Dnucleisegmentation/diffusioncode.m"/* pName */
};

static emlrtECInfo k_emlrtECI = { 2,   /* nDims */
  113,                                 /* lineNo */
  14,                                  /* colNo */
  "diffusioncode",                     /* fName */
  "/media/data/Dropbox/code_bank/matlab/3Dnucleisegmentation/diffusioncode.m"/* pName */
};

static emlrtECInfo l_emlrtECI = { 2,   /* nDims */
  114,                                 /* lineNo */
  14,                                  /* colNo */
  "diffusioncode",                     /* fName */
  "/media/data/Dropbox/code_bank/matlab/3Dnucleisegmentation/diffusioncode.m"/* pName */
};

static emlrtECInfo m_emlrtECI = { 2,   /* nDims */
  115,                                 /* lineNo */
  14,                                  /* colNo */
  "diffusioncode",                     /* fName */
  "/media/data/Dropbox/code_bank/matlab/3Dnucleisegmentation/diffusioncode.m"/* pName */
};

static emlrtECInfo n_emlrtECI = { 2,   /* nDims */
  116,                                 /* lineNo */
  14,                                  /* colNo */
  "diffusioncode",                     /* fName */
  "/media/data/Dropbox/code_bank/matlab/3Dnucleisegmentation/diffusioncode.m"/* pName */
};

static emlrtECInfo o_emlrtECI = { 2,   /* nDims */
  117,                                 /* lineNo */
  15,                                  /* colNo */
  "diffusioncode",                     /* fName */
  "/media/data/Dropbox/code_bank/matlab/3Dnucleisegmentation/diffusioncode.m"/* pName */
};

static emlrtECInfo p_emlrtECI = { 2,   /* nDims */
  118,                                 /* lineNo */
  15,                                  /* colNo */
  "diffusioncode",                     /* fName */
  "/media/data/Dropbox/code_bank/matlab/3Dnucleisegmentation/diffusioncode.m"/* pName */
};

static emlrtECInfo q_emlrtECI = { 2,   /* nDims */
  119,                                 /* lineNo */
  15,                                  /* colNo */
  "diffusioncode",                     /* fName */
  "/media/data/Dropbox/code_bank/matlab/3Dnucleisegmentation/diffusioncode.m"/* pName */
};

static emlrtECInfo r_emlrtECI = { 2,   /* nDims */
  120,                                 /* lineNo */
  15,                                  /* colNo */
  "diffusioncode",                     /* fName */
  "/media/data/Dropbox/code_bank/matlab/3Dnucleisegmentation/diffusioncode.m"/* pName */
};

static emlrtECInfo s_emlrtECI = { 2,   /* nDims */
  122,                                 /* lineNo */
  14,                                  /* colNo */
  "diffusioncode",                     /* fName */
  "/media/data/Dropbox/code_bank/matlab/3Dnucleisegmentation/diffusioncode.m"/* pName */
};

static emlrtECInfo t_emlrtECI = { 2,   /* nDims */
  123,                                 /* lineNo */
  14,                                  /* colNo */
  "diffusioncode",                     /* fName */
  "/media/data/Dropbox/code_bank/matlab/3Dnucleisegmentation/diffusioncode.m"/* pName */
};

static emlrtECInfo u_emlrtECI = { 2,   /* nDims */
  124,                                 /* lineNo */
  14,                                  /* colNo */
  "diffusioncode",                     /* fName */
  "/media/data/Dropbox/code_bank/matlab/3Dnucleisegmentation/diffusioncode.m"/* pName */
};

static emlrtECInfo v_emlrtECI = { 2,   /* nDims */
  125,                                 /* lineNo */
  14,                                  /* colNo */
  "diffusioncode",                     /* fName */
  "/media/data/Dropbox/code_bank/matlab/3Dnucleisegmentation/diffusioncode.m"/* pName */
};

static emlrtECInfo w_emlrtECI = { 2,   /* nDims */
  126,                                 /* lineNo */
  15,                                  /* colNo */
  "diffusioncode",                     /* fName */
  "/media/data/Dropbox/code_bank/matlab/3Dnucleisegmentation/diffusioncode.m"/* pName */
};

static emlrtECInfo x_emlrtECI = { 2,   /* nDims */
  127,                                 /* lineNo */
  15,                                  /* colNo */
  "diffusioncode",                     /* fName */
  "/media/data/Dropbox/code_bank/matlab/3Dnucleisegmentation/diffusioncode.m"/* pName */
};

static emlrtECInfo y_emlrtECI = { 2,   /* nDims */
  128,                                 /* lineNo */
  15,                                  /* colNo */
  "diffusioncode",                     /* fName */
  "/media/data/Dropbox/code_bank/matlab/3Dnucleisegmentation/diffusioncode.m"/* pName */
};

static emlrtECInfo ab_emlrtECI = { 2,  /* nDims */
  129,                                 /* lineNo */
  15,                                  /* colNo */
  "diffusioncode",                     /* fName */
  "/media/data/Dropbox/code_bank/matlab/3Dnucleisegmentation/diffusioncode.m"/* pName */
};

static emlrtECInfo bb_emlrtECI = { 2,  /* nDims */
  137,                                 /* lineNo */
  9,                                   /* colNo */
  "diffusioncode",                     /* fName */
  "/media/data/Dropbox/code_bank/matlab/3Dnucleisegmentation/diffusioncode.m"/* pName */
};

static emlrtECInfo cb_emlrtECI = { 2,  /* nDims */
  137,                                 /* lineNo */
  33,                                  /* colNo */
  "diffusioncode",                     /* fName */
  "/media/data/Dropbox/code_bank/matlab/3Dnucleisegmentation/diffusioncode.m"/* pName */
};

static emlrtECInfo db_emlrtECI = { 2,  /* nDims */
  138,                                 /* lineNo */
  9,                                   /* colNo */
  "diffusioncode",                     /* fName */
  "/media/data/Dropbox/code_bank/matlab/3Dnucleisegmentation/diffusioncode.m"/* pName */
};

static emlrtECInfo eb_emlrtECI = { 2,  /* nDims */
  138,                                 /* lineNo */
  33,                                  /* colNo */
  "diffusioncode",                     /* fName */
  "/media/data/Dropbox/code_bank/matlab/3Dnucleisegmentation/diffusioncode.m"/* pName */
};

static emlrtECInfo fb_emlrtECI = { 2,  /* nDims */
  139,                                 /* lineNo */
  9,                                   /* colNo */
  "diffusioncode",                     /* fName */
  "/media/data/Dropbox/code_bank/matlab/3Dnucleisegmentation/diffusioncode.m"/* pName */
};

static emlrtECInfo gb_emlrtECI = { 2,  /* nDims */
  139,                                 /* lineNo */
  35,                                  /* colNo */
  "diffusioncode",                     /* fName */
  "/media/data/Dropbox/code_bank/matlab/3Dnucleisegmentation/diffusioncode.m"/* pName */
};

static emlrtECInfo hb_emlrtECI = { 2,  /* nDims */
  140,                                 /* lineNo */
  9,                                   /* colNo */
  "diffusioncode",                     /* fName */
  "/media/data/Dropbox/code_bank/matlab/3Dnucleisegmentation/diffusioncode.m"/* pName */
};

static emlrtECInfo ib_emlrtECI = { 2,  /* nDims */
  140,                                 /* lineNo */
  35,                                  /* colNo */
  "diffusioncode",                     /* fName */
  "/media/data/Dropbox/code_bank/matlab/3Dnucleisegmentation/diffusioncode.m"/* pName */
};

static emlrtECInfo jb_emlrtECI = { 2,  /* nDims */
  135,                                 /* lineNo */
  15,                                  /* colNo */
  "diffusioncode",                     /* fName */
  "/media/data/Dropbox/code_bank/matlab/3Dnucleisegmentation/diffusioncode.m"/* pName */
};

static emlrtRTEInfo wd_emlrtRTEI = { 82,/* lineNo */
  5,                                   /* colNo */
  "diffusioncode",                     /* fName */
  "/media/data/Dropbox/code_bank/matlab/3Dnucleisegmentation/diffusioncode.m"/* pName */
};

static emlrtRTEInfo xd_emlrtRTEI = { 83,/* lineNo */
  5,                                   /* colNo */
  "diffusioncode",                     /* fName */
  "/media/data/Dropbox/code_bank/matlab/3Dnucleisegmentation/diffusioncode.m"/* pName */
};

static emlrtRTEInfo ae_emlrtRTEI = { 84,/* lineNo */
  5,                                   /* colNo */
  "diffusioncode",                     /* fName */
  "/media/data/Dropbox/code_bank/matlab/3Dnucleisegmentation/diffusioncode.m"/* pName */
};

static emlrtRTEInfo be_emlrtRTEI = { 85,/* lineNo */
  5,                                   /* colNo */
  "diffusioncode",                     /* fName */
  "/media/data/Dropbox/code_bank/matlab/3Dnucleisegmentation/diffusioncode.m"/* pName */
};

static emlrtRTEInfo ce_emlrtRTEI = { 86,/* lineNo */
  5,                                   /* colNo */
  "diffusioncode",                     /* fName */
  "/media/data/Dropbox/code_bank/matlab/3Dnucleisegmentation/diffusioncode.m"/* pName */
};

static emlrtRTEInfo de_emlrtRTEI = { 87,/* lineNo */
  5,                                   /* colNo */
  "diffusioncode",                     /* fName */
  "/media/data/Dropbox/code_bank/matlab/3Dnucleisegmentation/diffusioncode.m"/* pName */
};

static emlrtRTEInfo ee_emlrtRTEI = { 88,/* lineNo */
  5,                                   /* colNo */
  "diffusioncode",                     /* fName */
  "/media/data/Dropbox/code_bank/matlab/3Dnucleisegmentation/diffusioncode.m"/* pName */
};

static emlrtRTEInfo fe_emlrtRTEI = { 89,/* lineNo */
  5,                                   /* colNo */
  "diffusioncode",                     /* fName */
  "/media/data/Dropbox/code_bank/matlab/3Dnucleisegmentation/diffusioncode.m"/* pName */
};

static emlrtRTEInfo ge_emlrtRTEI = { 94,/* lineNo */
  20,                                  /* colNo */
  "diffusioncode",                     /* fName */
  "/media/data/Dropbox/code_bank/matlab/3Dnucleisegmentation/diffusioncode.m"/* pName */
};

static emlrtRTEInfo he_emlrtRTEI = { 113,/* lineNo */
  24,                                  /* colNo */
  "diffusioncode",                     /* fName */
  "/media/data/Dropbox/code_bank/matlab/3Dnucleisegmentation/diffusioncode.m"/* pName */
};

static emlrtRTEInfo ie_emlrtRTEI = { 104,/* lineNo */
  23,                                  /* colNo */
  "diffusioncode",                     /* fName */
  "/media/data/Dropbox/code_bank/matlab/3Dnucleisegmentation/diffusioncode.m"/* pName */
};

static emlrtRTEInfo je_emlrtRTEI = { 122,/* lineNo */
  23,                                  /* colNo */
  "diffusioncode",                     /* fName */
  "/media/data/Dropbox/code_bank/matlab/3Dnucleisegmentation/diffusioncode.m"/* pName */
};

static emlrtRTEInfo ke_emlrtRTEI = { 95,/* lineNo */
  20,                                  /* colNo */
  "diffusioncode",                     /* fName */
  "/media/data/Dropbox/code_bank/matlab/3Dnucleisegmentation/diffusioncode.m"/* pName */
};

static emlrtRTEInfo le_emlrtRTEI = { 137,/* lineNo */
  9,                                   /* colNo */
  "diffusioncode",                     /* fName */
  "/media/data/Dropbox/code_bank/matlab/3Dnucleisegmentation/diffusioncode.m"/* pName */
};

static emlrtRTEInfo me_emlrtRTEI = { 105,/* lineNo */
  23,                                  /* colNo */
  "diffusioncode",                     /* fName */
  "/media/data/Dropbox/code_bank/matlab/3Dnucleisegmentation/diffusioncode.m"/* pName */
};

static emlrtRTEInfo ne_emlrtRTEI = { 122,/* lineNo */
  51,                                  /* colNo */
  "diffusioncode",                     /* fName */
  "/media/data/Dropbox/code_bank/matlab/3Dnucleisegmentation/diffusioncode.m"/* pName */
};

static emlrtRTEInfo oe_emlrtRTEI = { 113,/* lineNo */
  50,                                  /* colNo */
  "diffusioncode",                     /* fName */
  "/media/data/Dropbox/code_bank/matlab/3Dnucleisegmentation/diffusioncode.m"/* pName */
};

static emlrtRTEInfo pe_emlrtRTEI = { 137,/* lineNo */
  33,                                  /* colNo */
  "diffusioncode",                     /* fName */
  "/media/data/Dropbox/code_bank/matlab/3Dnucleisegmentation/diffusioncode.m"/* pName */
};

static emlrtRTEInfo qe_emlrtRTEI = { 96,/* lineNo */
  20,                                  /* colNo */
  "diffusioncode",                     /* fName */
  "/media/data/Dropbox/code_bank/matlab/3Dnucleisegmentation/diffusioncode.m"/* pName */
};

static emlrtRTEInfo re_emlrtRTEI = { 138,/* lineNo */
  9,                                   /* colNo */
  "diffusioncode",                     /* fName */
  "/media/data/Dropbox/code_bank/matlab/3Dnucleisegmentation/diffusioncode.m"/* pName */
};

static emlrtRTEInfo se_emlrtRTEI = { 106,/* lineNo */
  23,                                  /* colNo */
  "diffusioncode",                     /* fName */
  "/media/data/Dropbox/code_bank/matlab/3Dnucleisegmentation/diffusioncode.m"/* pName */
};

static emlrtRTEInfo te_emlrtRTEI = { 114,/* lineNo */
  24,                                  /* colNo */
  "diffusioncode",                     /* fName */
  "/media/data/Dropbox/code_bank/matlab/3Dnucleisegmentation/diffusioncode.m"/* pName */
};

static emlrtRTEInfo ue_emlrtRTEI = { 97,/* lineNo */
  20,                                  /* colNo */
  "diffusioncode",                     /* fName */
  "/media/data/Dropbox/code_bank/matlab/3Dnucleisegmentation/diffusioncode.m"/* pName */
};

static emlrtRTEInfo ve_emlrtRTEI = { 138,/* lineNo */
  33,                                  /* colNo */
  "diffusioncode",                     /* fName */
  "/media/data/Dropbox/code_bank/matlab/3Dnucleisegmentation/diffusioncode.m"/* pName */
};

static emlrtRTEInfo we_emlrtRTEI = { 123,/* lineNo */
  23,                                  /* colNo */
  "diffusioncode",                     /* fName */
  "/media/data/Dropbox/code_bank/matlab/3Dnucleisegmentation/diffusioncode.m"/* pName */
};

static emlrtRTEInfo xe_emlrtRTEI = { 139,/* lineNo */
  9,                                   /* colNo */
  "diffusioncode",                     /* fName */
  "/media/data/Dropbox/code_bank/matlab/3Dnucleisegmentation/diffusioncode.m"/* pName */
};

static emlrtRTEInfo ye_emlrtRTEI = { 107,/* lineNo */
  23,                                  /* colNo */
  "diffusioncode",                     /* fName */
  "/media/data/Dropbox/code_bank/matlab/3Dnucleisegmentation/diffusioncode.m"/* pName */
};

static emlrtRTEInfo af_emlrtRTEI = { 98,/* lineNo */
  21,                                  /* colNo */
  "diffusioncode",                     /* fName */
  "/media/data/Dropbox/code_bank/matlab/3Dnucleisegmentation/diffusioncode.m"/* pName */
};

static emlrtRTEInfo bf_emlrtRTEI = { 114,/* lineNo */
  50,                                  /* colNo */
  "diffusioncode",                     /* fName */
  "/media/data/Dropbox/code_bank/matlab/3Dnucleisegmentation/diffusioncode.m"/* pName */
};

static emlrtRTEInfo cf_emlrtRTEI = { 139,/* lineNo */
  35,                                  /* colNo */
  "diffusioncode",                     /* fName */
  "/media/data/Dropbox/code_bank/matlab/3Dnucleisegmentation/diffusioncode.m"/* pName */
};

static emlrtRTEInfo df_emlrtRTEI = { 123,/* lineNo */
  51,                                  /* colNo */
  "diffusioncode",                     /* fName */
  "/media/data/Dropbox/code_bank/matlab/3Dnucleisegmentation/diffusioncode.m"/* pName */
};

static emlrtRTEInfo ef_emlrtRTEI = { 99,/* lineNo */
  21,                                  /* colNo */
  "diffusioncode",                     /* fName */
  "/media/data/Dropbox/code_bank/matlab/3Dnucleisegmentation/diffusioncode.m"/* pName */
};

static emlrtRTEInfo ff_emlrtRTEI = { 108,/* lineNo */
  24,                                  /* colNo */
  "diffusioncode",                     /* fName */
  "/media/data/Dropbox/code_bank/matlab/3Dnucleisegmentation/diffusioncode.m"/* pName */
};

static emlrtRTEInfo gf_emlrtRTEI = { 140,/* lineNo */
  9,                                   /* colNo */
  "diffusioncode",                     /* fName */
  "/media/data/Dropbox/code_bank/matlab/3Dnucleisegmentation/diffusioncode.m"/* pName */
};

static emlrtRTEInfo hf_emlrtRTEI = { 115,/* lineNo */
  24,                                  /* colNo */
  "diffusioncode",                     /* fName */
  "/media/data/Dropbox/code_bank/matlab/3Dnucleisegmentation/diffusioncode.m"/* pName */
};

static emlrtRTEInfo if_emlrtRTEI = { 100,/* lineNo */
  21,                                  /* colNo */
  "diffusioncode",                     /* fName */
  "/media/data/Dropbox/code_bank/matlab/3Dnucleisegmentation/diffusioncode.m"/* pName */
};

static emlrtRTEInfo jf_emlrtRTEI = { 140,/* lineNo */
  35,                                  /* colNo */
  "diffusioncode",                     /* fName */
  "/media/data/Dropbox/code_bank/matlab/3Dnucleisegmentation/diffusioncode.m"/* pName */
};

static emlrtRTEInfo kf_emlrtRTEI = { 124,/* lineNo */
  23,                                  /* colNo */
  "diffusioncode",                     /* fName */
  "/media/data/Dropbox/code_bank/matlab/3Dnucleisegmentation/diffusioncode.m"/* pName */
};

static emlrtRTEInfo lf_emlrtRTEI = { 109,/* lineNo */
  24,                                  /* colNo */
  "diffusioncode",                     /* fName */
  "/media/data/Dropbox/code_bank/matlab/3Dnucleisegmentation/diffusioncode.m"/* pName */
};

static emlrtRTEInfo mf_emlrtRTEI = { 101,/* lineNo */
  21,                                  /* colNo */
  "diffusioncode",                     /* fName */
  "/media/data/Dropbox/code_bank/matlab/3Dnucleisegmentation/diffusioncode.m"/* pName */
};

static emlrtRTEInfo nf_emlrtRTEI = { 115,/* lineNo */
  50,                                  /* colNo */
  "diffusioncode",                     /* fName */
  "/media/data/Dropbox/code_bank/matlab/3Dnucleisegmentation/diffusioncode.m"/* pName */
};

static emlrtRTEInfo of_emlrtRTEI = { 124,/* lineNo */
  51,                                  /* colNo */
  "diffusioncode",                     /* fName */
  "/media/data/Dropbox/code_bank/matlab/3Dnucleisegmentation/diffusioncode.m"/* pName */
};

static emlrtRTEInfo pf_emlrtRTEI = { 110,/* lineNo */
  24,                                  /* colNo */
  "diffusioncode",                     /* fName */
  "/media/data/Dropbox/code_bank/matlab/3Dnucleisegmentation/diffusioncode.m"/* pName */
};

static emlrtRTEInfo qf_emlrtRTEI = { 116,/* lineNo */
  24,                                  /* colNo */
  "diffusioncode",                     /* fName */
  "/media/data/Dropbox/code_bank/matlab/3Dnucleisegmentation/diffusioncode.m"/* pName */
};

static emlrtRTEInfo rf_emlrtRTEI = { 111,/* lineNo */
  24,                                  /* colNo */
  "diffusioncode",                     /* fName */
  "/media/data/Dropbox/code_bank/matlab/3Dnucleisegmentation/diffusioncode.m"/* pName */
};

static emlrtRTEInfo sf_emlrtRTEI = { 125,/* lineNo */
  23,                                  /* colNo */
  "diffusioncode",                     /* fName */
  "/media/data/Dropbox/code_bank/matlab/3Dnucleisegmentation/diffusioncode.m"/* pName */
};

static emlrtRTEInfo tf_emlrtRTEI = { 116,/* lineNo */
  50,                                  /* colNo */
  "diffusioncode",                     /* fName */
  "/media/data/Dropbox/code_bank/matlab/3Dnucleisegmentation/diffusioncode.m"/* pName */
};

static emlrtRTEInfo uf_emlrtRTEI = { 125,/* lineNo */
  51,                                  /* colNo */
  "diffusioncode",                     /* fName */
  "/media/data/Dropbox/code_bank/matlab/3Dnucleisegmentation/diffusioncode.m"/* pName */
};

static emlrtRTEInfo vf_emlrtRTEI = { 117,/* lineNo */
  25,                                  /* colNo */
  "diffusioncode",                     /* fName */
  "/media/data/Dropbox/code_bank/matlab/3Dnucleisegmentation/diffusioncode.m"/* pName */
};

static emlrtRTEInfo wf_emlrtRTEI = { 126,/* lineNo */
  24,                                  /* colNo */
  "diffusioncode",                     /* fName */
  "/media/data/Dropbox/code_bank/matlab/3Dnucleisegmentation/diffusioncode.m"/* pName */
};

static emlrtRTEInfo xf_emlrtRTEI = { 117,/* lineNo */
  52,                                  /* colNo */
  "diffusioncode",                     /* fName */
  "/media/data/Dropbox/code_bank/matlab/3Dnucleisegmentation/diffusioncode.m"/* pName */
};

static emlrtRTEInfo yf_emlrtRTEI = { 126,/* lineNo */
  53,                                  /* colNo */
  "diffusioncode",                     /* fName */
  "/media/data/Dropbox/code_bank/matlab/3Dnucleisegmentation/diffusioncode.m"/* pName */
};

static emlrtRTEInfo ag_emlrtRTEI = { 118,/* lineNo */
  25,                                  /* colNo */
  "diffusioncode",                     /* fName */
  "/media/data/Dropbox/code_bank/matlab/3Dnucleisegmentation/diffusioncode.m"/* pName */
};

static emlrtRTEInfo bg_emlrtRTEI = { 127,/* lineNo */
  24,                                  /* colNo */
  "diffusioncode",                     /* fName */
  "/media/data/Dropbox/code_bank/matlab/3Dnucleisegmentation/diffusioncode.m"/* pName */
};

static emlrtRTEInfo cg_emlrtRTEI = { 118,/* lineNo */
  52,                                  /* colNo */
  "diffusioncode",                     /* fName */
  "/media/data/Dropbox/code_bank/matlab/3Dnucleisegmentation/diffusioncode.m"/* pName */
};

static emlrtRTEInfo dg_emlrtRTEI = { 127,/* lineNo */
  53,                                  /* colNo */
  "diffusioncode",                     /* fName */
  "/media/data/Dropbox/code_bank/matlab/3Dnucleisegmentation/diffusioncode.m"/* pName */
};

static emlrtRTEInfo eg_emlrtRTEI = { 119,/* lineNo */
  25,                                  /* colNo */
  "diffusioncode",                     /* fName */
  "/media/data/Dropbox/code_bank/matlab/3Dnucleisegmentation/diffusioncode.m"/* pName */
};

static emlrtRTEInfo fg_emlrtRTEI = { 119,/* lineNo */
  52,                                  /* colNo */
  "diffusioncode",                     /* fName */
  "/media/data/Dropbox/code_bank/matlab/3Dnucleisegmentation/diffusioncode.m"/* pName */
};

static emlrtRTEInfo gg_emlrtRTEI = { 128,/* lineNo */
  24,                                  /* colNo */
  "diffusioncode",                     /* fName */
  "/media/data/Dropbox/code_bank/matlab/3Dnucleisegmentation/diffusioncode.m"/* pName */
};

static emlrtRTEInfo hg_emlrtRTEI = { 128,/* lineNo */
  53,                                  /* colNo */
  "diffusioncode",                     /* fName */
  "/media/data/Dropbox/code_bank/matlab/3Dnucleisegmentation/diffusioncode.m"/* pName */
};

static emlrtRTEInfo ig_emlrtRTEI = { 120,/* lineNo */
  25,                                  /* colNo */
  "diffusioncode",                     /* fName */
  "/media/data/Dropbox/code_bank/matlab/3Dnucleisegmentation/diffusioncode.m"/* pName */
};

static emlrtRTEInfo jg_emlrtRTEI = { 120,/* lineNo */
  52,                                  /* colNo */
  "diffusioncode",                     /* fName */
  "/media/data/Dropbox/code_bank/matlab/3Dnucleisegmentation/diffusioncode.m"/* pName */
};

static emlrtRTEInfo kg_emlrtRTEI = { 129,/* lineNo */
  24,                                  /* colNo */
  "diffusioncode",                     /* fName */
  "/media/data/Dropbox/code_bank/matlab/3Dnucleisegmentation/diffusioncode.m"/* pName */
};

static emlrtRTEInfo lg_emlrtRTEI = { 129,/* lineNo */
  53,                                  /* colNo */
  "diffusioncode",                     /* fName */
  "/media/data/Dropbox/code_bank/matlab/3Dnucleisegmentation/diffusioncode.m"/* pName */
};

static emlrtRTEInfo mg_emlrtRTEI = { 94,/* lineNo */
  9,                                   /* colNo */
  "diffusioncode",                     /* fName */
  "/media/data/Dropbox/code_bank/matlab/3Dnucleisegmentation/diffusioncode.m"/* pName */
};

static emlrtRTEInfo ng_emlrtRTEI = { 95,/* lineNo */
  9,                                   /* colNo */
  "diffusioncode",                     /* fName */
  "/media/data/Dropbox/code_bank/matlab/3Dnucleisegmentation/diffusioncode.m"/* pName */
};

static emlrtRTEInfo og_emlrtRTEI = { 96,/* lineNo */
  9,                                   /* colNo */
  "diffusioncode",                     /* fName */
  "/media/data/Dropbox/code_bank/matlab/3Dnucleisegmentation/diffusioncode.m"/* pName */
};

static emlrtRTEInfo pg_emlrtRTEI = { 97,/* lineNo */
  9,                                   /* colNo */
  "diffusioncode",                     /* fName */
  "/media/data/Dropbox/code_bank/matlab/3Dnucleisegmentation/diffusioncode.m"/* pName */
};

static emlrtRTEInfo qg_emlrtRTEI = { 98,/* lineNo */
  9,                                   /* colNo */
  "diffusioncode",                     /* fName */
  "/media/data/Dropbox/code_bank/matlab/3Dnucleisegmentation/diffusioncode.m"/* pName */
};

static emlrtRTEInfo rg_emlrtRTEI = { 99,/* lineNo */
  9,                                   /* colNo */
  "diffusioncode",                     /* fName */
  "/media/data/Dropbox/code_bank/matlab/3Dnucleisegmentation/diffusioncode.m"/* pName */
};

static emlrtRTEInfo sg_emlrtRTEI = { 100,/* lineNo */
  9,                                   /* colNo */
  "diffusioncode",                     /* fName */
  "/media/data/Dropbox/code_bank/matlab/3Dnucleisegmentation/diffusioncode.m"/* pName */
};

static emlrtRTEInfo tg_emlrtRTEI = { 101,/* lineNo */
  9,                                   /* colNo */
  "diffusioncode",                     /* fName */
  "/media/data/Dropbox/code_bank/matlab/3Dnucleisegmentation/diffusioncode.m"/* pName */
};

static emlrtRTEInfo ug_emlrtRTEI = { 25,/* lineNo */
  20,                                  /* colNo */
  "diffusioncode",                     /* fName */
  "/media/data/Dropbox/code_bank/matlab/3Dnucleisegmentation/diffusioncode.m"/* pName */
};

static emlrtRTEInfo wg_emlrtRTEI = { 113,/* lineNo */
  49,                                  /* colNo */
  "diffusioncode",                     /* fName */
  "/media/data/Dropbox/code_bank/matlab/3Dnucleisegmentation/diffusioncode.m"/* pName */
};

static emlrtRSInfo ah_emlrtRSI = { 27, /* lineNo */
  "error",                             /* fcnName */
  "/usr/local/MATLAB/R2019b/toolbox/eml/lib/matlab/lang/error.m"/* pathName */
};

/* Function Declarations */
static void error(const emlrtStack *sp, const mxArray *b, emlrtMCInfo *location);

/* Function Definitions */
static void error(const emlrtStack *sp, const mxArray *b, emlrtMCInfo *location)
{
  const mxArray *pArray;
  pArray = b;
  emlrtCallMATLABR2012b(sp, 0, NULL, 1, &pArray, "error", true, location);
}

void diffusioncode(const emlrtStack *sp, emxArray_real_T *im, real_T num_iter,
                   real_T kappa1, real_T kappa2, real_T option)
{
  int32_T i;
  emxArray_real_T *nablaN;
  emxArray_real_T *nablaS;
  emxArray_real_T *nablaW;
  emxArray_real_T *nablaE;
  emxArray_real_T *nablaNE;
  emxArray_real_T *nablaSE;
  emxArray_real_T *nablaSW;
  emxArray_real_T *nablaNW;
  emxArray_real_T *cN;
  emxArray_real_T *cS;
  emxArray_real_T *cW;
  emxArray_real_T *cE;
  emxArray_real_T *cNE;
  emxArray_real_T *cSE;
  emxArray_real_T *cSW;
  emxArray_real_T *cNW;
  emxArray_real_T *r;
  emxArray_real_T *a;
  emxArray_real_T *r1;
  emxArray_real_T *b_nablaN;
  int32_T t;
  real_T outSizeT[2];
  real_T startT[2];
  int32_T i1;
  int32_T loop_ub;
  boolean_T tooBig;
  int32_T i2;
  real_T padSizeT[2];
  real_T nonZeroKernel[2];
  real_T connDimsT[2];
  static const real_T kernel[9] = { 0.0, 0.0, 0.0, 1.0, -1.0, 0.0, 0.0, 0.0, 0.0
  };

  static const boolean_T conn[9] = { false, false, false, true, true, false,
    false, false, false };

  const mxArray *y;
  const mxArray *m;
  static const int32_T iv[2] = { 1, 16 };

  static const char_T varargin_1[16] = { 'I', 'n', 'c', 'o', 'r', 'r', 'e', 'c',
    't', ' ', 'o', 'p', 't', 'i', 'o', 'n' };

  int32_T b_loop_ub;
  emlrtStack st;
  emlrtStack b_st;
  emlrtStack c_st;
  emlrtStack d_st;
  st.prev = sp;
  st.tls = sp->tls;
  b_st.prev = &st;
  b_st.tls = st.tls;
  c_st.prev = &b_st;
  c_st.tls = b_st.tls;
  d_st.prev = &c_st;
  d_st.tls = c_st.tls;
  emlrtHeapReferenceStackEnterFcnR2012b(sp);

  /*  Adapated by Bhavna Rajasekaran */
  /*  References: */
  /*    P. Perona and J. Malik. */
  /*    Scale-Space and Edge Detection Using Anisotropic Diffusion. */
  /*    IEEE Transactions on Pattern Analysis and Machine Intelligence, */
  /*    12(7):629-639, July 1990. */
  /*    P. D. Kovesi. MATLAB and Octave Functions for Computer Vision and Image Processing. */
  /*    School of Computer Science & Software Engineering, */
  /*    The University of Western Australia. Available from: */
  /*    <http://www.csse.uwa.edu.au/~pk/research/matlabfns/>. */
  /*   */
  /*  Original Code: anisodiff3D.m */
  /*  Daniel Simoes Lopes */
  /*  ICIST */
  /*  Instituto Superior Tecnico - Universidade Tecnica de Lisboa */
  /*  danlopes (at) civil ist utl pt */
  /*  http://www.civil.ist.utl.pt/~danlopes */
  /*  Code modified to include forward and backward anisotropic diffusion */
  /* On the combined forward and backward anisotropic diffusion */
  /* scheme for the multispectral image enhancement */
  /* Bogdan Smolka Rastislav Lukac */
  /* ANISODIFF2D Conventional anisotropic diffusion */
  /*    DIFF_IM = ANISODIFF2D(IM, NUM_ITER, DELTA_T, KAPPA, OPTION) perfoms */
  /*    conventional anisotropic diffusion (Perona & Malik) upon a gray scale */
  /*    image. A 2D network structure of 8 neighboring nodes is considered for */
  /*    diffusion conduction. */
  /*  */
  /*        ARGUMENT DESCRIPTION: */
  /*                IM       - gray scale image (MxN). */
  /*                NUM_ITER - number of iterations. */
  /*                DELTA_T  - integration constant (0 <= delta_t <= 1/7). */
  /*                           Usually, due to numerical stability this */
  /*                           parameter is set to its maximum value. */
  /*                KAPPA    - gradient modulus threshold that controls the conduction. */
  /*                OPTION   - conduction coefficient functions proposed by Perona & Malik: */
  /*                           1 - c(x,y,t) = exp(-(nablaI/kappa).^2), */
  /*                               privileges high-contrast edges over low-contrast ones. */
  /*                           2 - c(x,y,t) = 1./(1 + (nablaI/kappa).^2), */
  /*                               privileges wide regions over smaller ones. */
  /* %% Extra options for forward-backward diffusions%%%%% */
  /*                           3- Forward-backward diffusion c(x,y,t) = */
  /*                           2.*exp(-(nablaI/kappa1).^2)-exp(-(nablaI/kappa2).^2), */
  /*                           4- Forward-backward diffusion c(x,y,t) =2./(1 +(nablaI/kappa1).^2)-1./(1 + (nablaI/kappa2).^2) */
  /*        OUTPUT DESCRIPTION: */
  /*                 DIFF_IM - (diffused) image with the largest scale-space parameter. */
  /*  */
  /*  */
  /*  Convert input image to double. */
  /*  PDE (partial differential equation) initial condition. */
  /*  Center pixel distances. */
  /* dd = sqrt(2); */
  /* dx = 0.5; */
  /* dy = 0.5; */
  /*  2D convolution masks - finite differences. */
  /*  Anisotropic diffusion. */
  i = static_cast<int32_T>(num_iter);
  emlrtForLoopVectorCheckR2012b(1.0, 1.0, num_iter, mxDOUBLE_CLASS,
    static_cast<int32_T>(num_iter), &cb_emlrtRTEI, sp);
  emxInit_real_T(sp, &nablaN, 2, &wd_emlrtRTEI, true);
  emxInit_real_T(sp, &nablaS, 2, &xd_emlrtRTEI, true);
  emxInit_real_T(sp, &nablaW, 2, &ae_emlrtRTEI, true);
  emxInit_real_T(sp, &nablaE, 2, &be_emlrtRTEI, true);
  emxInit_real_T(sp, &nablaNE, 2, &ce_emlrtRTEI, true);
  emxInit_real_T(sp, &nablaSE, 2, &de_emlrtRTEI, true);
  emxInit_real_T(sp, &nablaSW, 2, &ee_emlrtRTEI, true);
  emxInit_real_T(sp, &nablaNW, 2, &fe_emlrtRTEI, true);
  emxInit_real_T(sp, &cN, 2, &mg_emlrtRTEI, true);
  emxInit_real_T(sp, &cS, 2, &ng_emlrtRTEI, true);
  emxInit_real_T(sp, &cW, 2, &og_emlrtRTEI, true);
  emxInit_real_T(sp, &cE, 2, &pg_emlrtRTEI, true);
  emxInit_real_T(sp, &cNE, 2, &qg_emlrtRTEI, true);
  emxInit_real_T(sp, &cSE, 2, &rg_emlrtRTEI, true);
  emxInit_real_T(sp, &cSW, 2, &sg_emlrtRTEI, true);
  emxInit_real_T(sp, &cNW, 2, &tg_emlrtRTEI, true);
  emxInit_real_T(sp, &r, 2, &ug_emlrtRTEI, true);
  emxInit_real_T(sp, &a, 2, &vg_emlrtRTEI, true);
  emxInit_real_T(sp, &r1, 2, &wg_emlrtRTEI, true);
  emxInit_real_T(sp, &b_nablaN, 2, &je_emlrtRTEI, true);
  for (t = 0; t < i; t++) {
    /*  Finite differences. [imfilter(.,.,'conv') can be replaced by conv2(.,.,'same')] */
    st.site = &gc_emlrtRSI;
    outSizeT[0] = im->size[0];
    startT[0] = 1.0;
    outSizeT[1] = im->size[1];
    startT[1] = 1.0;
    if ((im->size[0] == 0) || (im->size[1] == 0)) {
      i1 = nablaN->size[0] * nablaN->size[1];
      nablaN->size[0] = im->size[0];
      nablaN->size[1] = im->size[1];
      emxEnsureCapacity_real_T(&st, nablaN, i1, &wd_emlrtRTEI);
      loop_ub = im->size[0] * im->size[1];
      for (i1 = 0; i1 < loop_ub; i1++) {
        nablaN->data[i1] = im->data[i1];
      }
    } else {
      b_st.site = &wd_emlrtRSI;
      c_st.site = &yd_emlrtRSI;
      padarray(&c_st, im, startT, a);
      b_st.site = &xd_emlrtRSI;
      c_st.site = &ae_emlrtRSI;
      tooBig = (outSizeT[0] > 65500.0);
      if ((!tooBig) || (!(outSizeT[1] > 65500.0))) {
        tooBig = false;
      }

      tooBig = !tooBig;
      d_st.site = &be_emlrtRSI;
      i1 = nablaN->size[0] * nablaN->size[1];
      nablaN->size[0] = static_cast<int32_T>(outSizeT[0]);
      nablaN->size[1] = static_cast<int32_T>(outSizeT[1]);
      emxEnsureCapacity_real_T(&d_st, nablaN, i1, &yd_emlrtRTEI);
      if (tooBig) {
        padSizeT[0] = a->size[0];
        startT[0] = 3.0;
        padSizeT[1] = a->size[1];
        startT[1] = 3.0;
        ippfilter_real64(&a->data[0], &nablaN->data[0], outSizeT, 2.0, padSizeT,
                         kernel, startT, true);
      } else {
        padSizeT[0] = a->size[0];
        nonZeroKernel[0] = 1.0;
        connDimsT[0] = 3.0;
        padSizeT[1] = a->size[1];
        nonZeroKernel[1] = -1.0;
        connDimsT[1] = 3.0;
        imfilter_real64(&a->data[0], &nablaN->data[0], 2.0, outSizeT, 2.0,
                        padSizeT, nonZeroKernel, 2.0, conn, 2.0, connDimsT,
                        startT, 2.0, true, true);
      }
    }

    i1 = im->size[0];
    i2 = nablaS->size[0] * nablaS->size[1];
    nablaS->size[0] = i1;
    i1 = im->size[1];
    nablaS->size[1] = i1;
    emxEnsureCapacity_real_T(sp, nablaS, i2, &xd_emlrtRTEI);
    loop_ub = im->size[0] * im->size[1];
    for (i1 = 0; i1 < loop_ub; i1++) {
      nablaS->data[i1] = im->data[i1];
    }

    st.site = &hc_emlrtRSI;
    imfilter(&st, nablaS);
    i1 = im->size[0];
    i2 = nablaW->size[0] * nablaW->size[1];
    nablaW->size[0] = i1;
    i1 = im->size[1];
    nablaW->size[1] = i1;
    emxEnsureCapacity_real_T(sp, nablaW, i2, &ae_emlrtRTEI);
    loop_ub = im->size[0] * im->size[1];
    for (i1 = 0; i1 < loop_ub; i1++) {
      nablaW->data[i1] = im->data[i1];
    }

    st.site = &ic_emlrtRSI;
    b_imfilter(&st, nablaW);
    i1 = im->size[0];
    i2 = nablaE->size[0] * nablaE->size[1];
    nablaE->size[0] = i1;
    i1 = im->size[1];
    nablaE->size[1] = i1;
    emxEnsureCapacity_real_T(sp, nablaE, i2, &be_emlrtRTEI);
    loop_ub = im->size[0] * im->size[1];
    for (i1 = 0; i1 < loop_ub; i1++) {
      nablaE->data[i1] = im->data[i1];
    }

    st.site = &jc_emlrtRSI;
    c_imfilter(&st, nablaE);
    i1 = im->size[0];
    i2 = nablaNE->size[0] * nablaNE->size[1];
    nablaNE->size[0] = i1;
    i1 = im->size[1];
    nablaNE->size[1] = i1;
    emxEnsureCapacity_real_T(sp, nablaNE, i2, &ce_emlrtRTEI);
    loop_ub = im->size[0] * im->size[1];
    for (i1 = 0; i1 < loop_ub; i1++) {
      nablaNE->data[i1] = im->data[i1];
    }

    st.site = &kc_emlrtRSI;
    d_imfilter(&st, nablaNE);
    i1 = im->size[0];
    i2 = nablaSE->size[0] * nablaSE->size[1];
    nablaSE->size[0] = i1;
    i1 = im->size[1];
    nablaSE->size[1] = i1;
    emxEnsureCapacity_real_T(sp, nablaSE, i2, &de_emlrtRTEI);
    loop_ub = im->size[0] * im->size[1];
    for (i1 = 0; i1 < loop_ub; i1++) {
      nablaSE->data[i1] = im->data[i1];
    }

    st.site = &lc_emlrtRSI;
    e_imfilter(&st, nablaSE);
    i1 = im->size[0];
    i2 = nablaSW->size[0] * nablaSW->size[1];
    nablaSW->size[0] = i1;
    i1 = im->size[1];
    nablaSW->size[1] = i1;
    emxEnsureCapacity_real_T(sp, nablaSW, i2, &ee_emlrtRTEI);
    loop_ub = im->size[0] * im->size[1];
    for (i1 = 0; i1 < loop_ub; i1++) {
      nablaSW->data[i1] = im->data[i1];
    }

    st.site = &mc_emlrtRSI;
    f_imfilter(&st, nablaSW);
    i1 = im->size[0];
    i2 = nablaNW->size[0] * nablaNW->size[1];
    nablaNW->size[0] = i1;
    i1 = im->size[1];
    nablaNW->size[1] = i1;
    emxEnsureCapacity_real_T(sp, nablaNW, i2, &fe_emlrtRTEI);
    loop_ub = im->size[0] * im->size[1];
    for (i1 = 0; i1 < loop_ub; i1++) {
      nablaNW->data[i1] = im->data[i1];
    }

    st.site = &nc_emlrtRSI;
    g_imfilter(&st, nablaNW);

    /*  Diffusion function. */
    if (option == 1.0) {
      kappa2 = 0.0;
      i1 = b_nablaN->size[0] * b_nablaN->size[1];
      b_nablaN->size[0] = nablaN->size[0];
      b_nablaN->size[1] = nablaN->size[1];
      emxEnsureCapacity_real_T(sp, b_nablaN, i1, &ge_emlrtRTEI);
      loop_ub = nablaN->size[0] * nablaN->size[1];
      for (i1 = 0; i1 < loop_ub; i1++) {
        b_nablaN->data[i1] = nablaN->data[i1] / kappa1;
      }

      st.site = &oc_emlrtRSI;
      power(&st, b_nablaN, cN);
      loop_ub = cN->size[0] * cN->size[1];
      for (i1 = 0; i1 < loop_ub; i1++) {
        cN->data[i1] = -cN->data[i1];
      }

      st.site = &oc_emlrtRSI;
      b_exp(&st, cN);
      i1 = b_nablaN->size[0] * b_nablaN->size[1];
      b_nablaN->size[0] = nablaS->size[0];
      b_nablaN->size[1] = nablaS->size[1];
      emxEnsureCapacity_real_T(sp, b_nablaN, i1, &ke_emlrtRTEI);
      loop_ub = nablaS->size[0] * nablaS->size[1];
      for (i1 = 0; i1 < loop_ub; i1++) {
        b_nablaN->data[i1] = nablaS->data[i1] / kappa1;
      }

      st.site = &pc_emlrtRSI;
      power(&st, b_nablaN, cS);
      loop_ub = cS->size[0] * cS->size[1];
      for (i1 = 0; i1 < loop_ub; i1++) {
        cS->data[i1] = -cS->data[i1];
      }

      st.site = &pc_emlrtRSI;
      b_exp(&st, cS);
      i1 = b_nablaN->size[0] * b_nablaN->size[1];
      b_nablaN->size[0] = nablaW->size[0];
      b_nablaN->size[1] = nablaW->size[1];
      emxEnsureCapacity_real_T(sp, b_nablaN, i1, &qe_emlrtRTEI);
      loop_ub = nablaW->size[0] * nablaW->size[1];
      for (i1 = 0; i1 < loop_ub; i1++) {
        b_nablaN->data[i1] = nablaW->data[i1] / kappa1;
      }

      st.site = &qc_emlrtRSI;
      power(&st, b_nablaN, cW);
      loop_ub = cW->size[0] * cW->size[1];
      for (i1 = 0; i1 < loop_ub; i1++) {
        cW->data[i1] = -cW->data[i1];
      }

      st.site = &qc_emlrtRSI;
      b_exp(&st, cW);
      i1 = b_nablaN->size[0] * b_nablaN->size[1];
      b_nablaN->size[0] = nablaE->size[0];
      b_nablaN->size[1] = nablaE->size[1];
      emxEnsureCapacity_real_T(sp, b_nablaN, i1, &ue_emlrtRTEI);
      loop_ub = nablaE->size[0] * nablaE->size[1];
      for (i1 = 0; i1 < loop_ub; i1++) {
        b_nablaN->data[i1] = nablaE->data[i1] / kappa1;
      }

      st.site = &rc_emlrtRSI;
      power(&st, b_nablaN, cE);
      loop_ub = cE->size[0] * cE->size[1];
      for (i1 = 0; i1 < loop_ub; i1++) {
        cE->data[i1] = -cE->data[i1];
      }

      st.site = &rc_emlrtRSI;
      b_exp(&st, cE);
      i1 = b_nablaN->size[0] * b_nablaN->size[1];
      b_nablaN->size[0] = nablaNE->size[0];
      b_nablaN->size[1] = nablaNE->size[1];
      emxEnsureCapacity_real_T(sp, b_nablaN, i1, &af_emlrtRTEI);
      loop_ub = nablaNE->size[0] * nablaNE->size[1];
      for (i1 = 0; i1 < loop_ub; i1++) {
        b_nablaN->data[i1] = nablaNE->data[i1] / kappa1;
      }

      st.site = &sc_emlrtRSI;
      power(&st, b_nablaN, cNE);
      loop_ub = cNE->size[0] * cNE->size[1];
      for (i1 = 0; i1 < loop_ub; i1++) {
        cNE->data[i1] = -cNE->data[i1];
      }

      st.site = &sc_emlrtRSI;
      b_exp(&st, cNE);
      i1 = b_nablaN->size[0] * b_nablaN->size[1];
      b_nablaN->size[0] = nablaSE->size[0];
      b_nablaN->size[1] = nablaSE->size[1];
      emxEnsureCapacity_real_T(sp, b_nablaN, i1, &ef_emlrtRTEI);
      loop_ub = nablaSE->size[0] * nablaSE->size[1];
      for (i1 = 0; i1 < loop_ub; i1++) {
        b_nablaN->data[i1] = nablaSE->data[i1] / kappa1;
      }

      st.site = &tc_emlrtRSI;
      power(&st, b_nablaN, cSE);
      loop_ub = cSE->size[0] * cSE->size[1];
      for (i1 = 0; i1 < loop_ub; i1++) {
        cSE->data[i1] = -cSE->data[i1];
      }

      st.site = &tc_emlrtRSI;
      b_exp(&st, cSE);
      i1 = b_nablaN->size[0] * b_nablaN->size[1];
      b_nablaN->size[0] = nablaSW->size[0];
      b_nablaN->size[1] = nablaSW->size[1];
      emxEnsureCapacity_real_T(sp, b_nablaN, i1, &if_emlrtRTEI);
      loop_ub = nablaSW->size[0] * nablaSW->size[1];
      for (i1 = 0; i1 < loop_ub; i1++) {
        b_nablaN->data[i1] = nablaSW->data[i1] / kappa1;
      }

      st.site = &uc_emlrtRSI;
      power(&st, b_nablaN, cSW);
      loop_ub = cSW->size[0] * cSW->size[1];
      for (i1 = 0; i1 < loop_ub; i1++) {
        cSW->data[i1] = -cSW->data[i1];
      }

      st.site = &uc_emlrtRSI;
      b_exp(&st, cSW);
      i1 = b_nablaN->size[0] * b_nablaN->size[1];
      b_nablaN->size[0] = nablaNW->size[0];
      b_nablaN->size[1] = nablaNW->size[1];
      emxEnsureCapacity_real_T(sp, b_nablaN, i1, &mf_emlrtRTEI);
      loop_ub = nablaNW->size[0] * nablaNW->size[1];
      for (i1 = 0; i1 < loop_ub; i1++) {
        b_nablaN->data[i1] = nablaNW->data[i1] / kappa1;
      }

      st.site = &vc_emlrtRSI;
      power(&st, b_nablaN, cNW);
      loop_ub = cNW->size[0] * cNW->size[1];
      for (i1 = 0; i1 < loop_ub; i1++) {
        cNW->data[i1] = -cNW->data[i1];
      }

      st.site = &vc_emlrtRSI;
      b_exp(&st, cNW);
    } else if (option == 2.0) {
      kappa2 = 0.0;
      i1 = b_nablaN->size[0] * b_nablaN->size[1];
      b_nablaN->size[0] = nablaN->size[0];
      b_nablaN->size[1] = nablaN->size[1];
      emxEnsureCapacity_real_T(sp, b_nablaN, i1, &ie_emlrtRTEI);
      loop_ub = nablaN->size[0] * nablaN->size[1];
      for (i1 = 0; i1 < loop_ub; i1++) {
        b_nablaN->data[i1] = nablaN->data[i1] / kappa1;
      }

      st.site = &wc_emlrtRSI;
      power(&st, b_nablaN, cN);
      loop_ub = cN->size[0] * cN->size[1];
      for (i1 = 0; i1 < loop_ub; i1++) {
        cN->data[i1]++;
      }

      loop_ub = cN->size[0] * cN->size[1];
      for (i1 = 0; i1 < loop_ub; i1++) {
        cN->data[i1] = 1.0 / cN->data[i1];
      }

      i1 = b_nablaN->size[0] * b_nablaN->size[1];
      b_nablaN->size[0] = nablaS->size[0];
      b_nablaN->size[1] = nablaS->size[1];
      emxEnsureCapacity_real_T(sp, b_nablaN, i1, &me_emlrtRTEI);
      loop_ub = nablaS->size[0] * nablaS->size[1];
      for (i1 = 0; i1 < loop_ub; i1++) {
        b_nablaN->data[i1] = nablaS->data[i1] / kappa1;
      }

      st.site = &xc_emlrtRSI;
      power(&st, b_nablaN, cS);
      loop_ub = cS->size[0] * cS->size[1];
      for (i1 = 0; i1 < loop_ub; i1++) {
        cS->data[i1]++;
      }

      loop_ub = cS->size[0] * cS->size[1];
      for (i1 = 0; i1 < loop_ub; i1++) {
        cS->data[i1] = 1.0 / cS->data[i1];
      }

      i1 = b_nablaN->size[0] * b_nablaN->size[1];
      b_nablaN->size[0] = nablaW->size[0];
      b_nablaN->size[1] = nablaW->size[1];
      emxEnsureCapacity_real_T(sp, b_nablaN, i1, &se_emlrtRTEI);
      loop_ub = nablaW->size[0] * nablaW->size[1];
      for (i1 = 0; i1 < loop_ub; i1++) {
        b_nablaN->data[i1] = nablaW->data[i1] / kappa1;
      }

      st.site = &yc_emlrtRSI;
      power(&st, b_nablaN, cW);
      loop_ub = cW->size[0] * cW->size[1];
      for (i1 = 0; i1 < loop_ub; i1++) {
        cW->data[i1]++;
      }

      loop_ub = cW->size[0] * cW->size[1];
      for (i1 = 0; i1 < loop_ub; i1++) {
        cW->data[i1] = 1.0 / cW->data[i1];
      }

      i1 = b_nablaN->size[0] * b_nablaN->size[1];
      b_nablaN->size[0] = nablaE->size[0];
      b_nablaN->size[1] = nablaE->size[1];
      emxEnsureCapacity_real_T(sp, b_nablaN, i1, &ye_emlrtRTEI);
      loop_ub = nablaE->size[0] * nablaE->size[1];
      for (i1 = 0; i1 < loop_ub; i1++) {
        b_nablaN->data[i1] = nablaE->data[i1] / kappa1;
      }

      st.site = &ad_emlrtRSI;
      power(&st, b_nablaN, cE);
      loop_ub = cE->size[0] * cE->size[1];
      for (i1 = 0; i1 < loop_ub; i1++) {
        cE->data[i1]++;
      }

      loop_ub = cE->size[0] * cE->size[1];
      for (i1 = 0; i1 < loop_ub; i1++) {
        cE->data[i1] = 1.0 / cE->data[i1];
      }

      i1 = b_nablaN->size[0] * b_nablaN->size[1];
      b_nablaN->size[0] = nablaNE->size[0];
      b_nablaN->size[1] = nablaNE->size[1];
      emxEnsureCapacity_real_T(sp, b_nablaN, i1, &ff_emlrtRTEI);
      loop_ub = nablaNE->size[0] * nablaNE->size[1];
      for (i1 = 0; i1 < loop_ub; i1++) {
        b_nablaN->data[i1] = nablaNE->data[i1] / kappa1;
      }

      st.site = &bd_emlrtRSI;
      power(&st, b_nablaN, cNE);
      loop_ub = cNE->size[0] * cNE->size[1];
      for (i1 = 0; i1 < loop_ub; i1++) {
        cNE->data[i1]++;
      }

      loop_ub = cNE->size[0] * cNE->size[1];
      for (i1 = 0; i1 < loop_ub; i1++) {
        cNE->data[i1] = 1.0 / cNE->data[i1];
      }

      i1 = b_nablaN->size[0] * b_nablaN->size[1];
      b_nablaN->size[0] = nablaSE->size[0];
      b_nablaN->size[1] = nablaSE->size[1];
      emxEnsureCapacity_real_T(sp, b_nablaN, i1, &lf_emlrtRTEI);
      loop_ub = nablaSE->size[0] * nablaSE->size[1];
      for (i1 = 0; i1 < loop_ub; i1++) {
        b_nablaN->data[i1] = nablaSE->data[i1] / kappa1;
      }

      st.site = &cd_emlrtRSI;
      power(&st, b_nablaN, cSE);
      loop_ub = cSE->size[0] * cSE->size[1];
      for (i1 = 0; i1 < loop_ub; i1++) {
        cSE->data[i1]++;
      }

      loop_ub = cSE->size[0] * cSE->size[1];
      for (i1 = 0; i1 < loop_ub; i1++) {
        cSE->data[i1] = 1.0 / cSE->data[i1];
      }

      i1 = b_nablaN->size[0] * b_nablaN->size[1];
      b_nablaN->size[0] = nablaSW->size[0];
      b_nablaN->size[1] = nablaSW->size[1];
      emxEnsureCapacity_real_T(sp, b_nablaN, i1, &pf_emlrtRTEI);
      loop_ub = nablaSW->size[0] * nablaSW->size[1];
      for (i1 = 0; i1 < loop_ub; i1++) {
        b_nablaN->data[i1] = nablaSW->data[i1] / kappa1;
      }

      st.site = &dd_emlrtRSI;
      power(&st, b_nablaN, cSW);
      loop_ub = cSW->size[0] * cSW->size[1];
      for (i1 = 0; i1 < loop_ub; i1++) {
        cSW->data[i1]++;
      }

      loop_ub = cSW->size[0] * cSW->size[1];
      for (i1 = 0; i1 < loop_ub; i1++) {
        cSW->data[i1] = 1.0 / cSW->data[i1];
      }

      i1 = b_nablaN->size[0] * b_nablaN->size[1];
      b_nablaN->size[0] = nablaNW->size[0];
      b_nablaN->size[1] = nablaNW->size[1];
      emxEnsureCapacity_real_T(sp, b_nablaN, i1, &rf_emlrtRTEI);
      loop_ub = nablaNW->size[0] * nablaNW->size[1];
      for (i1 = 0; i1 < loop_ub; i1++) {
        b_nablaN->data[i1] = nablaNW->data[i1] / kappa1;
      }

      st.site = &ed_emlrtRSI;
      power(&st, b_nablaN, cNW);
      loop_ub = cNW->size[0] * cNW->size[1];
      for (i1 = 0; i1 < loop_ub; i1++) {
        cNW->data[i1]++;
      }

      loop_ub = cNW->size[0] * cNW->size[1];
      for (i1 = 0; i1 < loop_ub; i1++) {
        cNW->data[i1] = 1.0 / cNW->data[i1];
      }
    } else if (option == 3.0) {
      i1 = b_nablaN->size[0] * b_nablaN->size[1];
      b_nablaN->size[0] = nablaN->size[0];
      b_nablaN->size[1] = nablaN->size[1];
      emxEnsureCapacity_real_T(sp, b_nablaN, i1, &he_emlrtRTEI);
      loop_ub = nablaN->size[0] * nablaN->size[1];
      for (i1 = 0; i1 < loop_ub; i1++) {
        b_nablaN->data[i1] = nablaN->data[i1] / kappa1;
      }

      st.site = &fd_emlrtRSI;
      power(&st, b_nablaN, cN);
      loop_ub = cN->size[0] * cN->size[1];
      for (i1 = 0; i1 < loop_ub; i1++) {
        cN->data[i1] = -cN->data[i1];
      }

      st.site = &fd_emlrtRSI;
      b_exp(&st, cN);
      loop_ub = cN->size[0] * cN->size[1];
      for (i1 = 0; i1 < loop_ub; i1++) {
        cN->data[i1] *= 2.0;
      }

      i1 = b_nablaN->size[0] * b_nablaN->size[1];
      b_nablaN->size[0] = nablaN->size[0];
      b_nablaN->size[1] = nablaN->size[1];
      emxEnsureCapacity_real_T(sp, b_nablaN, i1, &oe_emlrtRTEI);
      loop_ub = nablaN->size[0] * nablaN->size[1];
      for (i1 = 0; i1 < loop_ub; i1++) {
        b_nablaN->data[i1] = nablaN->data[i1] / kappa2;
      }

      st.site = &fd_emlrtRSI;
      power(&st, b_nablaN, r1);
      loop_ub = r1->size[0] * r1->size[1];
      for (i1 = 0; i1 < loop_ub; i1++) {
        r1->data[i1] = -r1->data[i1];
      }

      st.site = &fd_emlrtRSI;
      b_exp(&st, r1);
      emlrtSizeEqCheckNDR2012b(*(int32_T (*)[2])cN->size, *(int32_T (*)[2])
        r1->size, &k_emlrtECI, sp);
      loop_ub = cN->size[0] * cN->size[1];
      for (i1 = 0; i1 < loop_ub; i1++) {
        cN->data[i1] -= r1->data[i1];
      }

      i1 = b_nablaN->size[0] * b_nablaN->size[1];
      b_nablaN->size[0] = nablaS->size[0];
      b_nablaN->size[1] = nablaS->size[1];
      emxEnsureCapacity_real_T(sp, b_nablaN, i1, &te_emlrtRTEI);
      loop_ub = nablaS->size[0] * nablaS->size[1];
      for (i1 = 0; i1 < loop_ub; i1++) {
        b_nablaN->data[i1] = nablaS->data[i1] / kappa1;
      }

      st.site = &gd_emlrtRSI;
      power(&st, b_nablaN, cS);
      loop_ub = cS->size[0] * cS->size[1];
      for (i1 = 0; i1 < loop_ub; i1++) {
        cS->data[i1] = -cS->data[i1];
      }

      st.site = &gd_emlrtRSI;
      b_exp(&st, cS);
      loop_ub = cS->size[0] * cS->size[1];
      for (i1 = 0; i1 < loop_ub; i1++) {
        cS->data[i1] *= 2.0;
      }

      i1 = b_nablaN->size[0] * b_nablaN->size[1];
      b_nablaN->size[0] = nablaS->size[0];
      b_nablaN->size[1] = nablaS->size[1];
      emxEnsureCapacity_real_T(sp, b_nablaN, i1, &bf_emlrtRTEI);
      loop_ub = nablaS->size[0] * nablaS->size[1];
      for (i1 = 0; i1 < loop_ub; i1++) {
        b_nablaN->data[i1] = nablaS->data[i1] / kappa2;
      }

      st.site = &gd_emlrtRSI;
      power(&st, b_nablaN, r1);
      loop_ub = r1->size[0] * r1->size[1];
      for (i1 = 0; i1 < loop_ub; i1++) {
        r1->data[i1] = -r1->data[i1];
      }

      st.site = &gd_emlrtRSI;
      b_exp(&st, r1);
      emlrtSizeEqCheckNDR2012b(*(int32_T (*)[2])cS->size, *(int32_T (*)[2])
        r1->size, &l_emlrtECI, sp);
      loop_ub = cS->size[0] * cS->size[1];
      for (i1 = 0; i1 < loop_ub; i1++) {
        cS->data[i1] -= r1->data[i1];
      }

      i1 = b_nablaN->size[0] * b_nablaN->size[1];
      b_nablaN->size[0] = nablaW->size[0];
      b_nablaN->size[1] = nablaW->size[1];
      emxEnsureCapacity_real_T(sp, b_nablaN, i1, &hf_emlrtRTEI);
      loop_ub = nablaW->size[0] * nablaW->size[1];
      for (i1 = 0; i1 < loop_ub; i1++) {
        b_nablaN->data[i1] = nablaW->data[i1] / kappa1;
      }

      st.site = &hd_emlrtRSI;
      power(&st, b_nablaN, cW);
      loop_ub = cW->size[0] * cW->size[1];
      for (i1 = 0; i1 < loop_ub; i1++) {
        cW->data[i1] = -cW->data[i1];
      }

      st.site = &hd_emlrtRSI;
      b_exp(&st, cW);
      loop_ub = cW->size[0] * cW->size[1];
      for (i1 = 0; i1 < loop_ub; i1++) {
        cW->data[i1] *= 2.0;
      }

      i1 = b_nablaN->size[0] * b_nablaN->size[1];
      b_nablaN->size[0] = nablaW->size[0];
      b_nablaN->size[1] = nablaW->size[1];
      emxEnsureCapacity_real_T(sp, b_nablaN, i1, &nf_emlrtRTEI);
      loop_ub = nablaW->size[0] * nablaW->size[1];
      for (i1 = 0; i1 < loop_ub; i1++) {
        b_nablaN->data[i1] = nablaW->data[i1] / kappa2;
      }

      st.site = &hd_emlrtRSI;
      power(&st, b_nablaN, r1);
      loop_ub = r1->size[0] * r1->size[1];
      for (i1 = 0; i1 < loop_ub; i1++) {
        r1->data[i1] = -r1->data[i1];
      }

      st.site = &hd_emlrtRSI;
      b_exp(&st, r1);
      emlrtSizeEqCheckNDR2012b(*(int32_T (*)[2])cW->size, *(int32_T (*)[2])
        r1->size, &m_emlrtECI, sp);
      loop_ub = cW->size[0] * cW->size[1];
      for (i1 = 0; i1 < loop_ub; i1++) {
        cW->data[i1] -= r1->data[i1];
      }

      i1 = b_nablaN->size[0] * b_nablaN->size[1];
      b_nablaN->size[0] = nablaE->size[0];
      b_nablaN->size[1] = nablaE->size[1];
      emxEnsureCapacity_real_T(sp, b_nablaN, i1, &qf_emlrtRTEI);
      loop_ub = nablaE->size[0] * nablaE->size[1];
      for (i1 = 0; i1 < loop_ub; i1++) {
        b_nablaN->data[i1] = nablaE->data[i1] / kappa1;
      }

      st.site = &id_emlrtRSI;
      power(&st, b_nablaN, cE);
      loop_ub = cE->size[0] * cE->size[1];
      for (i1 = 0; i1 < loop_ub; i1++) {
        cE->data[i1] = -cE->data[i1];
      }

      st.site = &id_emlrtRSI;
      b_exp(&st, cE);
      loop_ub = cE->size[0] * cE->size[1];
      for (i1 = 0; i1 < loop_ub; i1++) {
        cE->data[i1] *= 2.0;
      }

      i1 = b_nablaN->size[0] * b_nablaN->size[1];
      b_nablaN->size[0] = nablaE->size[0];
      b_nablaN->size[1] = nablaE->size[1];
      emxEnsureCapacity_real_T(sp, b_nablaN, i1, &tf_emlrtRTEI);
      loop_ub = nablaE->size[0] * nablaE->size[1];
      for (i1 = 0; i1 < loop_ub; i1++) {
        b_nablaN->data[i1] = nablaE->data[i1] / kappa2;
      }

      st.site = &id_emlrtRSI;
      power(&st, b_nablaN, r1);
      loop_ub = r1->size[0] * r1->size[1];
      for (i1 = 0; i1 < loop_ub; i1++) {
        r1->data[i1] = -r1->data[i1];
      }

      st.site = &id_emlrtRSI;
      b_exp(&st, r1);
      emlrtSizeEqCheckNDR2012b(*(int32_T (*)[2])cE->size, *(int32_T (*)[2])
        r1->size, &n_emlrtECI, sp);
      loop_ub = cE->size[0] * cE->size[1];
      for (i1 = 0; i1 < loop_ub; i1++) {
        cE->data[i1] -= r1->data[i1];
      }

      i1 = b_nablaN->size[0] * b_nablaN->size[1];
      b_nablaN->size[0] = nablaNE->size[0];
      b_nablaN->size[1] = nablaNE->size[1];
      emxEnsureCapacity_real_T(sp, b_nablaN, i1, &vf_emlrtRTEI);
      loop_ub = nablaNE->size[0] * nablaNE->size[1];
      for (i1 = 0; i1 < loop_ub; i1++) {
        b_nablaN->data[i1] = nablaNE->data[i1] / kappa1;
      }

      st.site = &jd_emlrtRSI;
      power(&st, b_nablaN, cNE);
      loop_ub = cNE->size[0] * cNE->size[1];
      for (i1 = 0; i1 < loop_ub; i1++) {
        cNE->data[i1] = -cNE->data[i1];
      }

      st.site = &jd_emlrtRSI;
      b_exp(&st, cNE);
      loop_ub = cNE->size[0] * cNE->size[1];
      for (i1 = 0; i1 < loop_ub; i1++) {
        cNE->data[i1] *= 2.0;
      }

      i1 = b_nablaN->size[0] * b_nablaN->size[1];
      b_nablaN->size[0] = nablaNE->size[0];
      b_nablaN->size[1] = nablaNE->size[1];
      emxEnsureCapacity_real_T(sp, b_nablaN, i1, &xf_emlrtRTEI);
      loop_ub = nablaNE->size[0] * nablaNE->size[1];
      for (i1 = 0; i1 < loop_ub; i1++) {
        b_nablaN->data[i1] = nablaNE->data[i1] / kappa2;
      }

      st.site = &jd_emlrtRSI;
      power(&st, b_nablaN, r1);
      loop_ub = r1->size[0] * r1->size[1];
      for (i1 = 0; i1 < loop_ub; i1++) {
        r1->data[i1] = -r1->data[i1];
      }

      st.site = &jd_emlrtRSI;
      b_exp(&st, r1);
      emlrtSizeEqCheckNDR2012b(*(int32_T (*)[2])cNE->size, *(int32_T (*)[2])
        r1->size, &o_emlrtECI, sp);
      loop_ub = cNE->size[0] * cNE->size[1];
      for (i1 = 0; i1 < loop_ub; i1++) {
        cNE->data[i1] -= r1->data[i1];
      }

      i1 = b_nablaN->size[0] * b_nablaN->size[1];
      b_nablaN->size[0] = nablaSE->size[0];
      b_nablaN->size[1] = nablaSE->size[1];
      emxEnsureCapacity_real_T(sp, b_nablaN, i1, &ag_emlrtRTEI);
      loop_ub = nablaSE->size[0] * nablaSE->size[1];
      for (i1 = 0; i1 < loop_ub; i1++) {
        b_nablaN->data[i1] = nablaSE->data[i1] / kappa1;
      }

      st.site = &kd_emlrtRSI;
      power(&st, b_nablaN, cSE);
      loop_ub = cSE->size[0] * cSE->size[1];
      for (i1 = 0; i1 < loop_ub; i1++) {
        cSE->data[i1] = -cSE->data[i1];
      }

      st.site = &kd_emlrtRSI;
      b_exp(&st, cSE);
      loop_ub = cSE->size[0] * cSE->size[1];
      for (i1 = 0; i1 < loop_ub; i1++) {
        cSE->data[i1] *= 2.0;
      }

      i1 = b_nablaN->size[0] * b_nablaN->size[1];
      b_nablaN->size[0] = nablaSE->size[0];
      b_nablaN->size[1] = nablaSE->size[1];
      emxEnsureCapacity_real_T(sp, b_nablaN, i1, &cg_emlrtRTEI);
      loop_ub = nablaSE->size[0] * nablaSE->size[1];
      for (i1 = 0; i1 < loop_ub; i1++) {
        b_nablaN->data[i1] = nablaSE->data[i1] / kappa2;
      }

      st.site = &kd_emlrtRSI;
      power(&st, b_nablaN, r1);
      loop_ub = r1->size[0] * r1->size[1];
      for (i1 = 0; i1 < loop_ub; i1++) {
        r1->data[i1] = -r1->data[i1];
      }

      st.site = &kd_emlrtRSI;
      b_exp(&st, r1);
      emlrtSizeEqCheckNDR2012b(*(int32_T (*)[2])cSE->size, *(int32_T (*)[2])
        r1->size, &p_emlrtECI, sp);
      loop_ub = cSE->size[0] * cSE->size[1];
      for (i1 = 0; i1 < loop_ub; i1++) {
        cSE->data[i1] -= r1->data[i1];
      }

      i1 = b_nablaN->size[0] * b_nablaN->size[1];
      b_nablaN->size[0] = nablaSW->size[0];
      b_nablaN->size[1] = nablaSW->size[1];
      emxEnsureCapacity_real_T(sp, b_nablaN, i1, &eg_emlrtRTEI);
      loop_ub = nablaSW->size[0] * nablaSW->size[1];
      for (i1 = 0; i1 < loop_ub; i1++) {
        b_nablaN->data[i1] = nablaSW->data[i1] / kappa1;
      }

      st.site = &ld_emlrtRSI;
      power(&st, b_nablaN, cSW);
      loop_ub = cSW->size[0] * cSW->size[1];
      for (i1 = 0; i1 < loop_ub; i1++) {
        cSW->data[i1] = -cSW->data[i1];
      }

      st.site = &ld_emlrtRSI;
      b_exp(&st, cSW);
      loop_ub = cSW->size[0] * cSW->size[1];
      for (i1 = 0; i1 < loop_ub; i1++) {
        cSW->data[i1] *= 2.0;
      }

      i1 = b_nablaN->size[0] * b_nablaN->size[1];
      b_nablaN->size[0] = nablaSW->size[0];
      b_nablaN->size[1] = nablaSW->size[1];
      emxEnsureCapacity_real_T(sp, b_nablaN, i1, &fg_emlrtRTEI);
      loop_ub = nablaSW->size[0] * nablaSW->size[1];
      for (i1 = 0; i1 < loop_ub; i1++) {
        b_nablaN->data[i1] = nablaSW->data[i1] / kappa2;
      }

      st.site = &ld_emlrtRSI;
      power(&st, b_nablaN, r1);
      loop_ub = r1->size[0] * r1->size[1];
      for (i1 = 0; i1 < loop_ub; i1++) {
        r1->data[i1] = -r1->data[i1];
      }

      st.site = &ld_emlrtRSI;
      b_exp(&st, r1);
      emlrtSizeEqCheckNDR2012b(*(int32_T (*)[2])cSW->size, *(int32_T (*)[2])
        r1->size, &q_emlrtECI, sp);
      loop_ub = cSW->size[0] * cSW->size[1];
      for (i1 = 0; i1 < loop_ub; i1++) {
        cSW->data[i1] -= r1->data[i1];
      }

      i1 = b_nablaN->size[0] * b_nablaN->size[1];
      b_nablaN->size[0] = nablaNW->size[0];
      b_nablaN->size[1] = nablaNW->size[1];
      emxEnsureCapacity_real_T(sp, b_nablaN, i1, &ig_emlrtRTEI);
      loop_ub = nablaNW->size[0] * nablaNW->size[1];
      for (i1 = 0; i1 < loop_ub; i1++) {
        b_nablaN->data[i1] = nablaNW->data[i1] / kappa1;
      }

      st.site = &md_emlrtRSI;
      power(&st, b_nablaN, cNW);
      loop_ub = cNW->size[0] * cNW->size[1];
      for (i1 = 0; i1 < loop_ub; i1++) {
        cNW->data[i1] = -cNW->data[i1];
      }

      st.site = &md_emlrtRSI;
      b_exp(&st, cNW);
      loop_ub = cNW->size[0] * cNW->size[1];
      for (i1 = 0; i1 < loop_ub; i1++) {
        cNW->data[i1] *= 2.0;
      }

      i1 = b_nablaN->size[0] * b_nablaN->size[1];
      b_nablaN->size[0] = nablaNW->size[0];
      b_nablaN->size[1] = nablaNW->size[1];
      emxEnsureCapacity_real_T(sp, b_nablaN, i1, &jg_emlrtRTEI);
      loop_ub = nablaNW->size[0] * nablaNW->size[1];
      for (i1 = 0; i1 < loop_ub; i1++) {
        b_nablaN->data[i1] = nablaNW->data[i1] / kappa2;
      }

      st.site = &md_emlrtRSI;
      power(&st, b_nablaN, r1);
      loop_ub = r1->size[0] * r1->size[1];
      for (i1 = 0; i1 < loop_ub; i1++) {
        r1->data[i1] = -r1->data[i1];
      }

      st.site = &md_emlrtRSI;
      b_exp(&st, r1);
      emlrtSizeEqCheckNDR2012b(*(int32_T (*)[2])cNW->size, *(int32_T (*)[2])
        r1->size, &r_emlrtECI, sp);
      loop_ub = cNW->size[0] * cNW->size[1];
      for (i1 = 0; i1 < loop_ub; i1++) {
        cNW->data[i1] -= r1->data[i1];
      }
    } else if (option == 4.0) {
      i1 = b_nablaN->size[0] * b_nablaN->size[1];
      b_nablaN->size[0] = nablaN->size[0];
      b_nablaN->size[1] = nablaN->size[1];
      emxEnsureCapacity_real_T(sp, b_nablaN, i1, &je_emlrtRTEI);
      loop_ub = nablaN->size[0] * nablaN->size[1];
      for (i1 = 0; i1 < loop_ub; i1++) {
        b_nablaN->data[i1] = nablaN->data[i1] / kappa1;
      }

      st.site = &nd_emlrtRSI;
      power(&st, b_nablaN, cN);
      loop_ub = cN->size[0] * cN->size[1];
      for (i1 = 0; i1 < loop_ub; i1++) {
        cN->data[i1]++;
      }

      loop_ub = cN->size[0] * cN->size[1];
      for (i1 = 0; i1 < loop_ub; i1++) {
        cN->data[i1] = 2.0 / cN->data[i1];
      }

      i1 = b_nablaN->size[0] * b_nablaN->size[1];
      b_nablaN->size[0] = nablaN->size[0];
      b_nablaN->size[1] = nablaN->size[1];
      emxEnsureCapacity_real_T(sp, b_nablaN, i1, &ne_emlrtRTEI);
      loop_ub = nablaN->size[0] * nablaN->size[1];
      for (i1 = 0; i1 < loop_ub; i1++) {
        b_nablaN->data[i1] = nablaN->data[i1] / kappa2;
      }

      st.site = &nd_emlrtRSI;
      power(&st, b_nablaN, a);
      loop_ub = a->size[0] * a->size[1];
      for (i1 = 0; i1 < loop_ub; i1++) {
        a->data[i1]++;
      }

      loop_ub = a->size[0] * a->size[1];
      for (i1 = 0; i1 < loop_ub; i1++) {
        a->data[i1] = 1.0 / a->data[i1];
      }

      emlrtSizeEqCheckNDR2012b(*(int32_T (*)[2])cN->size, *(int32_T (*)[2])
        a->size, &s_emlrtECI, sp);
      loop_ub = cN->size[0] * cN->size[1];
      for (i1 = 0; i1 < loop_ub; i1++) {
        cN->data[i1] -= a->data[i1];
      }

      i1 = b_nablaN->size[0] * b_nablaN->size[1];
      b_nablaN->size[0] = nablaS->size[0];
      b_nablaN->size[1] = nablaS->size[1];
      emxEnsureCapacity_real_T(sp, b_nablaN, i1, &we_emlrtRTEI);
      loop_ub = nablaS->size[0] * nablaS->size[1];
      for (i1 = 0; i1 < loop_ub; i1++) {
        b_nablaN->data[i1] = nablaS->data[i1] / kappa1;
      }

      st.site = &od_emlrtRSI;
      power(&st, b_nablaN, cS);
      loop_ub = cS->size[0] * cS->size[1];
      for (i1 = 0; i1 < loop_ub; i1++) {
        cS->data[i1]++;
      }

      loop_ub = cS->size[0] * cS->size[1];
      for (i1 = 0; i1 < loop_ub; i1++) {
        cS->data[i1] = 2.0 / cS->data[i1];
      }

      i1 = b_nablaN->size[0] * b_nablaN->size[1];
      b_nablaN->size[0] = nablaS->size[0];
      b_nablaN->size[1] = nablaS->size[1];
      emxEnsureCapacity_real_T(sp, b_nablaN, i1, &df_emlrtRTEI);
      loop_ub = nablaS->size[0] * nablaS->size[1];
      for (i1 = 0; i1 < loop_ub; i1++) {
        b_nablaN->data[i1] = nablaS->data[i1] / kappa2;
      }

      st.site = &od_emlrtRSI;
      power(&st, b_nablaN, a);
      loop_ub = a->size[0] * a->size[1];
      for (i1 = 0; i1 < loop_ub; i1++) {
        a->data[i1]++;
      }

      loop_ub = a->size[0] * a->size[1];
      for (i1 = 0; i1 < loop_ub; i1++) {
        a->data[i1] = 1.0 / a->data[i1];
      }

      emlrtSizeEqCheckNDR2012b(*(int32_T (*)[2])cS->size, *(int32_T (*)[2])
        a->size, &t_emlrtECI, sp);
      loop_ub = cS->size[0] * cS->size[1];
      for (i1 = 0; i1 < loop_ub; i1++) {
        cS->data[i1] -= a->data[i1];
      }

      i1 = b_nablaN->size[0] * b_nablaN->size[1];
      b_nablaN->size[0] = nablaW->size[0];
      b_nablaN->size[1] = nablaW->size[1];
      emxEnsureCapacity_real_T(sp, b_nablaN, i1, &kf_emlrtRTEI);
      loop_ub = nablaW->size[0] * nablaW->size[1];
      for (i1 = 0; i1 < loop_ub; i1++) {
        b_nablaN->data[i1] = nablaW->data[i1] / kappa1;
      }

      st.site = &pd_emlrtRSI;
      power(&st, b_nablaN, cW);
      loop_ub = cW->size[0] * cW->size[1];
      for (i1 = 0; i1 < loop_ub; i1++) {
        cW->data[i1]++;
      }

      loop_ub = cW->size[0] * cW->size[1];
      for (i1 = 0; i1 < loop_ub; i1++) {
        cW->data[i1] = 2.0 / cW->data[i1];
      }

      i1 = b_nablaN->size[0] * b_nablaN->size[1];
      b_nablaN->size[0] = nablaW->size[0];
      b_nablaN->size[1] = nablaW->size[1];
      emxEnsureCapacity_real_T(sp, b_nablaN, i1, &of_emlrtRTEI);
      loop_ub = nablaW->size[0] * nablaW->size[1];
      for (i1 = 0; i1 < loop_ub; i1++) {
        b_nablaN->data[i1] = nablaW->data[i1] / kappa2;
      }

      st.site = &pd_emlrtRSI;
      power(&st, b_nablaN, a);
      loop_ub = a->size[0] * a->size[1];
      for (i1 = 0; i1 < loop_ub; i1++) {
        a->data[i1]++;
      }

      loop_ub = a->size[0] * a->size[1];
      for (i1 = 0; i1 < loop_ub; i1++) {
        a->data[i1] = 1.0 / a->data[i1];
      }

      emlrtSizeEqCheckNDR2012b(*(int32_T (*)[2])cW->size, *(int32_T (*)[2])
        a->size, &u_emlrtECI, sp);
      loop_ub = cW->size[0] * cW->size[1];
      for (i1 = 0; i1 < loop_ub; i1++) {
        cW->data[i1] -= a->data[i1];
      }

      i1 = b_nablaN->size[0] * b_nablaN->size[1];
      b_nablaN->size[0] = nablaE->size[0];
      b_nablaN->size[1] = nablaE->size[1];
      emxEnsureCapacity_real_T(sp, b_nablaN, i1, &sf_emlrtRTEI);
      loop_ub = nablaE->size[0] * nablaE->size[1];
      for (i1 = 0; i1 < loop_ub; i1++) {
        b_nablaN->data[i1] = nablaE->data[i1] / kappa1;
      }

      st.site = &qd_emlrtRSI;
      power(&st, b_nablaN, cE);
      loop_ub = cE->size[0] * cE->size[1];
      for (i1 = 0; i1 < loop_ub; i1++) {
        cE->data[i1]++;
      }

      loop_ub = cE->size[0] * cE->size[1];
      for (i1 = 0; i1 < loop_ub; i1++) {
        cE->data[i1] = 2.0 / cE->data[i1];
      }

      i1 = b_nablaN->size[0] * b_nablaN->size[1];
      b_nablaN->size[0] = nablaE->size[0];
      b_nablaN->size[1] = nablaE->size[1];
      emxEnsureCapacity_real_T(sp, b_nablaN, i1, &uf_emlrtRTEI);
      loop_ub = nablaE->size[0] * nablaE->size[1];
      for (i1 = 0; i1 < loop_ub; i1++) {
        b_nablaN->data[i1] = nablaE->data[i1] / kappa2;
      }

      st.site = &qd_emlrtRSI;
      power(&st, b_nablaN, a);
      loop_ub = a->size[0] * a->size[1];
      for (i1 = 0; i1 < loop_ub; i1++) {
        a->data[i1]++;
      }

      loop_ub = a->size[0] * a->size[1];
      for (i1 = 0; i1 < loop_ub; i1++) {
        a->data[i1] = 1.0 / a->data[i1];
      }

      emlrtSizeEqCheckNDR2012b(*(int32_T (*)[2])cE->size, *(int32_T (*)[2])
        a->size, &v_emlrtECI, sp);
      loop_ub = cE->size[0] * cE->size[1];
      for (i1 = 0; i1 < loop_ub; i1++) {
        cE->data[i1] -= a->data[i1];
      }

      i1 = b_nablaN->size[0] * b_nablaN->size[1];
      b_nablaN->size[0] = nablaNE->size[0];
      b_nablaN->size[1] = nablaNE->size[1];
      emxEnsureCapacity_real_T(sp, b_nablaN, i1, &wf_emlrtRTEI);
      loop_ub = nablaNE->size[0] * nablaNE->size[1];
      for (i1 = 0; i1 < loop_ub; i1++) {
        b_nablaN->data[i1] = nablaNE->data[i1] / kappa1;
      }

      st.site = &rd_emlrtRSI;
      power(&st, b_nablaN, cNE);
      loop_ub = cNE->size[0] * cNE->size[1];
      for (i1 = 0; i1 < loop_ub; i1++) {
        cNE->data[i1]++;
      }

      loop_ub = cNE->size[0] * cNE->size[1];
      for (i1 = 0; i1 < loop_ub; i1++) {
        cNE->data[i1] = 2.0 / cNE->data[i1];
      }

      i1 = b_nablaN->size[0] * b_nablaN->size[1];
      b_nablaN->size[0] = nablaNE->size[0];
      b_nablaN->size[1] = nablaNE->size[1];
      emxEnsureCapacity_real_T(sp, b_nablaN, i1, &yf_emlrtRTEI);
      loop_ub = nablaNE->size[0] * nablaNE->size[1];
      for (i1 = 0; i1 < loop_ub; i1++) {
        b_nablaN->data[i1] = nablaNE->data[i1] / kappa2;
      }

      st.site = &rd_emlrtRSI;
      power(&st, b_nablaN, a);
      loop_ub = a->size[0] * a->size[1];
      for (i1 = 0; i1 < loop_ub; i1++) {
        a->data[i1]++;
      }

      loop_ub = a->size[0] * a->size[1];
      for (i1 = 0; i1 < loop_ub; i1++) {
        a->data[i1] = 1.0 / a->data[i1];
      }

      emlrtSizeEqCheckNDR2012b(*(int32_T (*)[2])cNE->size, *(int32_T (*)[2])
        a->size, &w_emlrtECI, sp);
      loop_ub = cNE->size[0] * cNE->size[1];
      for (i1 = 0; i1 < loop_ub; i1++) {
        cNE->data[i1] -= a->data[i1];
      }

      i1 = b_nablaN->size[0] * b_nablaN->size[1];
      b_nablaN->size[0] = nablaSE->size[0];
      b_nablaN->size[1] = nablaSE->size[1];
      emxEnsureCapacity_real_T(sp, b_nablaN, i1, &bg_emlrtRTEI);
      loop_ub = nablaSE->size[0] * nablaSE->size[1];
      for (i1 = 0; i1 < loop_ub; i1++) {
        b_nablaN->data[i1] = nablaSE->data[i1] / kappa1;
      }

      st.site = &sd_emlrtRSI;
      power(&st, b_nablaN, cSE);
      loop_ub = cSE->size[0] * cSE->size[1];
      for (i1 = 0; i1 < loop_ub; i1++) {
        cSE->data[i1]++;
      }

      loop_ub = cSE->size[0] * cSE->size[1];
      for (i1 = 0; i1 < loop_ub; i1++) {
        cSE->data[i1] = 2.0 / cSE->data[i1];
      }

      i1 = b_nablaN->size[0] * b_nablaN->size[1];
      b_nablaN->size[0] = nablaSE->size[0];
      b_nablaN->size[1] = nablaSE->size[1];
      emxEnsureCapacity_real_T(sp, b_nablaN, i1, &dg_emlrtRTEI);
      loop_ub = nablaSE->size[0] * nablaSE->size[1];
      for (i1 = 0; i1 < loop_ub; i1++) {
        b_nablaN->data[i1] = nablaSE->data[i1] / kappa2;
      }

      st.site = &sd_emlrtRSI;
      power(&st, b_nablaN, a);
      loop_ub = a->size[0] * a->size[1];
      for (i1 = 0; i1 < loop_ub; i1++) {
        a->data[i1]++;
      }

      loop_ub = a->size[0] * a->size[1];
      for (i1 = 0; i1 < loop_ub; i1++) {
        a->data[i1] = 1.0 / a->data[i1];
      }

      emlrtSizeEqCheckNDR2012b(*(int32_T (*)[2])cSE->size, *(int32_T (*)[2])
        a->size, &x_emlrtECI, sp);
      loop_ub = cSE->size[0] * cSE->size[1];
      for (i1 = 0; i1 < loop_ub; i1++) {
        cSE->data[i1] -= a->data[i1];
      }

      i1 = b_nablaN->size[0] * b_nablaN->size[1];
      b_nablaN->size[0] = nablaSW->size[0];
      b_nablaN->size[1] = nablaSW->size[1];
      emxEnsureCapacity_real_T(sp, b_nablaN, i1, &gg_emlrtRTEI);
      loop_ub = nablaSW->size[0] * nablaSW->size[1];
      for (i1 = 0; i1 < loop_ub; i1++) {
        b_nablaN->data[i1] = nablaSW->data[i1] / kappa1;
      }

      st.site = &td_emlrtRSI;
      power(&st, b_nablaN, cSW);
      loop_ub = cSW->size[0] * cSW->size[1];
      for (i1 = 0; i1 < loop_ub; i1++) {
        cSW->data[i1]++;
      }

      loop_ub = cSW->size[0] * cSW->size[1];
      for (i1 = 0; i1 < loop_ub; i1++) {
        cSW->data[i1] = 2.0 / cSW->data[i1];
      }

      i1 = b_nablaN->size[0] * b_nablaN->size[1];
      b_nablaN->size[0] = nablaSW->size[0];
      b_nablaN->size[1] = nablaSW->size[1];
      emxEnsureCapacity_real_T(sp, b_nablaN, i1, &hg_emlrtRTEI);
      loop_ub = nablaSW->size[0] * nablaSW->size[1];
      for (i1 = 0; i1 < loop_ub; i1++) {
        b_nablaN->data[i1] = nablaSW->data[i1] / kappa2;
      }

      st.site = &td_emlrtRSI;
      power(&st, b_nablaN, a);
      loop_ub = a->size[0] * a->size[1];
      for (i1 = 0; i1 < loop_ub; i1++) {
        a->data[i1]++;
      }

      loop_ub = a->size[0] * a->size[1];
      for (i1 = 0; i1 < loop_ub; i1++) {
        a->data[i1] = 1.0 / a->data[i1];
      }

      emlrtSizeEqCheckNDR2012b(*(int32_T (*)[2])cSW->size, *(int32_T (*)[2])
        a->size, &y_emlrtECI, sp);
      loop_ub = cSW->size[0] * cSW->size[1];
      for (i1 = 0; i1 < loop_ub; i1++) {
        cSW->data[i1] -= a->data[i1];
      }

      i1 = b_nablaN->size[0] * b_nablaN->size[1];
      b_nablaN->size[0] = nablaNW->size[0];
      b_nablaN->size[1] = nablaNW->size[1];
      emxEnsureCapacity_real_T(sp, b_nablaN, i1, &kg_emlrtRTEI);
      loop_ub = nablaNW->size[0] * nablaNW->size[1];
      for (i1 = 0; i1 < loop_ub; i1++) {
        b_nablaN->data[i1] = nablaNW->data[i1] / kappa1;
      }

      st.site = &ud_emlrtRSI;
      power(&st, b_nablaN, cNW);
      loop_ub = cNW->size[0] * cNW->size[1];
      for (i1 = 0; i1 < loop_ub; i1++) {
        cNW->data[i1]++;
      }

      loop_ub = cNW->size[0] * cNW->size[1];
      for (i1 = 0; i1 < loop_ub; i1++) {
        cNW->data[i1] = 2.0 / cNW->data[i1];
      }

      i1 = b_nablaN->size[0] * b_nablaN->size[1];
      b_nablaN->size[0] = nablaNW->size[0];
      b_nablaN->size[1] = nablaNW->size[1];
      emxEnsureCapacity_real_T(sp, b_nablaN, i1, &lg_emlrtRTEI);
      loop_ub = nablaNW->size[0] * nablaNW->size[1];
      for (i1 = 0; i1 < loop_ub; i1++) {
        b_nablaN->data[i1] = nablaNW->data[i1] / kappa2;
      }

      st.site = &ud_emlrtRSI;
      power(&st, b_nablaN, a);
      loop_ub = a->size[0] * a->size[1];
      for (i1 = 0; i1 < loop_ub; i1++) {
        a->data[i1]++;
      }

      loop_ub = a->size[0] * a->size[1];
      for (i1 = 0; i1 < loop_ub; i1++) {
        a->data[i1] = 1.0 / a->data[i1];
      }

      emlrtSizeEqCheckNDR2012b(*(int32_T (*)[2])cNW->size, *(int32_T (*)[2])
        a->size, &ab_emlrtECI, sp);
      loop_ub = cNW->size[0] * cNW->size[1];
      for (i1 = 0; i1 < loop_ub; i1++) {
        cNW->data[i1] -= a->data[i1];
      }
    } else {
      st.site = &vd_emlrtRSI;
      y = NULL;
      m = emlrtCreateCharArray(2, iv);
      emlrtInitCharArrayR2013a(&st, 16, m, &varargin_1[0]);
      emlrtAssign(&y, m);
      b_st.site = &ah_emlrtRSI;
      error(&b_st, y, &emlrtMCI);
    }

    /*  Discrete PDE solution. */
    emlrtSizeEqCheckNDR2012b(*(int32_T (*)[2])cN->size, *(int32_T (*)[2])
      nablaN->size, &bb_emlrtECI, sp);
    emlrtSizeEqCheckNDR2012b(*(int32_T (*)[2])cS->size, *(int32_T (*)[2])
      nablaS->size, &cb_emlrtECI, sp);
    loop_ub = cN->size[0] * cN->size[1];
    i1 = nablaN->size[0] * nablaN->size[1];
    nablaN->size[0] = cN->size[0];
    nablaN->size[1] = cN->size[1];
    emxEnsureCapacity_real_T(sp, nablaN, i1, &le_emlrtRTEI);
    for (i1 = 0; i1 < loop_ub; i1++) {
      nablaN->data[i1] *= cN->data[i1];
    }

    loop_ub = cS->size[0] * cS->size[1];
    i1 = nablaS->size[0] * nablaS->size[1];
    nablaS->size[0] = cS->size[0];
    nablaS->size[1] = cS->size[1];
    emxEnsureCapacity_real_T(sp, nablaS, i1, &pe_emlrtRTEI);
    for (i1 = 0; i1 < loop_ub; i1++) {
      nablaS->data[i1] *= cS->data[i1];
    }

    emlrtSizeEqCheckNDR2012b(*(int32_T (*)[2])nablaN->size, *(int32_T (*)[2])
      nablaS->size, &bb_emlrtECI, sp);
    emlrtSizeEqCheckNDR2012b(*(int32_T (*)[2])cW->size, *(int32_T (*)[2])
      nablaW->size, &db_emlrtECI, sp);
    loop_ub = cW->size[0] * cW->size[1];
    i1 = nablaW->size[0] * nablaW->size[1];
    nablaW->size[0] = cW->size[0];
    nablaW->size[1] = cW->size[1];
    emxEnsureCapacity_real_T(sp, nablaW, i1, &re_emlrtRTEI);
    for (i1 = 0; i1 < loop_ub; i1++) {
      nablaW->data[i1] *= cW->data[i1];
    }

    emlrtSizeEqCheckNDR2012b(*(int32_T (*)[2])nablaN->size, *(int32_T (*)[2])
      nablaW->size, &bb_emlrtECI, sp);
    emlrtSizeEqCheckNDR2012b(*(int32_T (*)[2])cE->size, *(int32_T (*)[2])
      nablaE->size, &eb_emlrtECI, sp);
    loop_ub = cE->size[0] * cE->size[1];
    i1 = nablaE->size[0] * nablaE->size[1];
    nablaE->size[0] = cE->size[0];
    nablaE->size[1] = cE->size[1];
    emxEnsureCapacity_real_T(sp, nablaE, i1, &ve_emlrtRTEI);
    for (i1 = 0; i1 < loop_ub; i1++) {
      nablaE->data[i1] *= cE->data[i1];
    }

    emlrtSizeEqCheckNDR2012b(*(int32_T (*)[2])nablaN->size, *(int32_T (*)[2])
      nablaE->size, &bb_emlrtECI, sp);
    i1 = r1->size[0] * r1->size[1];
    r1->size[0] = cNE->size[0];
    r1->size[1] = cNE->size[1];
    emxEnsureCapacity_real_T(sp, r1, i1, &xe_emlrtRTEI);
    loop_ub = cNE->size[0] * cNE->size[1];
    for (i1 = 0; i1 < loop_ub; i1++) {
      r1->data[i1] = 0.49999999999999989 * cNE->data[i1];
    }

    emlrtSizeEqCheckNDR2012b(*(int32_T (*)[2])r1->size, *(int32_T (*)[2])
      nablaNE->size, &fb_emlrtECI, sp);
    loop_ub = r1->size[0] * r1->size[1];
    for (i1 = 0; i1 < loop_ub; i1++) {
      r1->data[i1] *= nablaNE->data[i1];
    }

    emlrtSizeEqCheckNDR2012b(*(int32_T (*)[2])nablaN->size, *(int32_T (*)[2])
      r1->size, &bb_emlrtECI, sp);
    i1 = r->size[0] * r->size[1];
    r->size[0] = cSE->size[0];
    r->size[1] = cSE->size[1];
    emxEnsureCapacity_real_T(sp, r, i1, &cf_emlrtRTEI);
    loop_ub = cSE->size[0] * cSE->size[1];
    for (i1 = 0; i1 < loop_ub; i1++) {
      r->data[i1] = 0.49999999999999989 * cSE->data[i1];
    }

    emlrtSizeEqCheckNDR2012b(*(int32_T (*)[2])r->size, *(int32_T (*)[2])
      nablaSE->size, &gb_emlrtECI, sp);
    loop_ub = r->size[0] * r->size[1];
    for (i1 = 0; i1 < loop_ub; i1++) {
      r->data[i1] *= nablaSE->data[i1];
    }

    emlrtSizeEqCheckNDR2012b(*(int32_T (*)[2])nablaN->size, *(int32_T (*)[2])
      r->size, &bb_emlrtECI, sp);
    i1 = nablaNE->size[0] * nablaNE->size[1];
    nablaNE->size[0] = cSW->size[0];
    nablaNE->size[1] = cSW->size[1];
    emxEnsureCapacity_real_T(sp, nablaNE, i1, &gf_emlrtRTEI);
    loop_ub = cSW->size[0] * cSW->size[1];
    for (i1 = 0; i1 < loop_ub; i1++) {
      nablaNE->data[i1] = 0.49999999999999989 * cSW->data[i1];
    }

    emlrtSizeEqCheckNDR2012b(*(int32_T (*)[2])nablaNE->size, *(int32_T (*)[2])
      nablaSW->size, &hb_emlrtECI, sp);
    loop_ub = nablaNE->size[0] * nablaNE->size[1];
    for (i1 = 0; i1 < loop_ub; i1++) {
      nablaNE->data[i1] *= nablaSW->data[i1];
    }

    emlrtSizeEqCheckNDR2012b(*(int32_T (*)[2])nablaN->size, *(int32_T (*)[2])
      nablaNE->size, &bb_emlrtECI, sp);
    i1 = a->size[0] * a->size[1];
    a->size[0] = cNW->size[0];
    a->size[1] = cNW->size[1];
    emxEnsureCapacity_real_T(sp, a, i1, &jf_emlrtRTEI);
    loop_ub = cNW->size[0] * cNW->size[1];
    for (i1 = 0; i1 < loop_ub; i1++) {
      a->data[i1] = 0.49999999999999989 * cNW->data[i1];
    }

    emlrtSizeEqCheckNDR2012b(*(int32_T (*)[2])a->size, *(int32_T (*)[2])
      nablaNW->size, &ib_emlrtECI, sp);
    loop_ub = a->size[0] * a->size[1];
    for (i1 = 0; i1 < loop_ub; i1++) {
      a->data[i1] *= nablaNW->data[i1];
    }

    emlrtSizeEqCheckNDR2012b(*(int32_T (*)[2])nablaN->size, *(int32_T (*)[2])
      a->size, &bb_emlrtECI, sp);
    loop_ub = nablaN->size[0] * nablaN->size[1];
    for (i1 = 0; i1 < loop_ub; i1++) {
      nablaN->data[i1] = 0.1429 * (((((((nablaN->data[i1] + nablaS->data[i1]) +
        nablaW->data[i1]) + nablaE->data[i1]) + r1->data[i1]) + r->data[i1]) +
        nablaNE->data[i1]) + a->data[i1]);
    }

    emlrtSizeEqCheckNDR2012b(*(int32_T (*)[2])im->size, *(int32_T (*)[2])
      nablaN->size, &jb_emlrtECI, sp);
    loop_ub = im->size[1];
    for (i1 = 0; i1 < loop_ub; i1++) {
      b_loop_ub = im->size[0];
      for (i2 = 0; i2 < b_loop_ub; i2++) {
        im->data[i2 + im->size[0] * i1] += nablaN->data[i2 + nablaN->size[0] *
          i1];
      }
    }

    if (*emlrtBreakCheckR2012bFlagVar != 0) {
      emlrtBreakCheckR2012b(sp);
    }
  }

  emxFree_real_T(&b_nablaN);
  emxFree_real_T(&r1);
  emxFree_real_T(&a);
  emxFree_real_T(&r);
  emxFree_real_T(&cNW);
  emxFree_real_T(&cSW);
  emxFree_real_T(&cSE);
  emxFree_real_T(&cNE);
  emxFree_real_T(&cE);
  emxFree_real_T(&cW);
  emxFree_real_T(&cS);
  emxFree_real_T(&cN);
  emxFree_real_T(&nablaNW);
  emxFree_real_T(&nablaSW);
  emxFree_real_T(&nablaSE);
  emxFree_real_T(&nablaNE);
  emxFree_real_T(&nablaE);
  emxFree_real_T(&nablaW);
  emxFree_real_T(&nablaS);
  emxFree_real_T(&nablaN);
  emlrtHeapReferenceStackLeaveFcnR2012b(sp);
}

/* End of code generation (diffusioncode.cpp) */
