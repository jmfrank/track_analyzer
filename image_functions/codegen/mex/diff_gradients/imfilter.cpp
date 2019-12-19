/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * imfilter.cpp
 *
 * Code generation for function 'imfilter'
 *
 */

/* Include files */
#include "imfilter.h"
#include "all.h"
#include "combineVectorElements.h"
#include "diag.h"
#include "diff_gradients.h"
#include "diff_gradients_data.h"
#include "diff_gradients_emxutil.h"
#include "eml_int_forloop_overflow_check.h"
#include "indexShapeCheck.h"
#include "libmwimfilter.h"
#include "libmwippfilter.h"
#include "mwmathutil.h"
#include "padarray.h"
#include "repmat.h"
#include "rt_nonfinite.h"
#include "sqrt.h"
#include "svd.h"
#include "validateattributes.h"

/* Variable Definitions */
static emlrtRSInfo ce_emlrtRSI = { 945,/* lineNo */
  "imfiltercoreAlgo",                  /* fcnName */
  "/usr/local/MATLAB/R2019b/toolbox/images/images/eml/imfilter.m"/* pathName */
};

static emlrtRSInfo de_emlrtRSI = { 957,/* lineNo */
  "imfiltercoreAlgo",                  /* fcnName */
  "/usr/local/MATLAB/R2019b/toolbox/images/images/eml/imfilter.m"/* pathName */
};

static emlrtRSInfo cf_emlrtRSI = { 55, /* lineNo */
  "imfilter",                          /* fcnName */
  "/usr/local/MATLAB/R2019b/toolbox/images/images/eml/imfilter.m"/* pathName */
};

static emlrtRSInfo df_emlrtRSI = { 59, /* lineNo */
  "imfilter",                          /* fcnName */
  "/usr/local/MATLAB/R2019b/toolbox/images/images/eml/imfilter.m"/* pathName */
};

static emlrtRSInfo ef_emlrtRSI = { 64, /* lineNo */
  "imfilter",                          /* fcnName */
  "/usr/local/MATLAB/R2019b/toolbox/images/images/eml/imfilter.m"/* pathName */
};

static emlrtRSInfo ff_emlrtRSI = { 66, /* lineNo */
  "imfilter",                          /* fcnName */
  "/usr/local/MATLAB/R2019b/toolbox/images/images/eml/imfilter.m"/* pathName */
};

static emlrtRSInfo gf_emlrtRSI = { 67, /* lineNo */
  "imfilter",                          /* fcnName */
  "/usr/local/MATLAB/R2019b/toolbox/images/images/eml/imfilter.m"/* pathName */
};

static emlrtRSInfo hf_emlrtRSI = { 68, /* lineNo */
  "imfilter",                          /* fcnName */
  "/usr/local/MATLAB/R2019b/toolbox/images/images/eml/imfilter.m"/* pathName */
};

static emlrtRSInfo if_emlrtRSI = { 84, /* lineNo */
  "imfilter",                          /* fcnName */
  "/usr/local/MATLAB/R2019b/toolbox/images/images/eml/imfilter.m"/* pathName */
};

static emlrtRSInfo jf_emlrtRSI = { 88, /* lineNo */
  "imfilter",                          /* fcnName */
  "/usr/local/MATLAB/R2019b/toolbox/images/images/eml/imfilter.m"/* pathName */
};

static emlrtRSInfo kf_emlrtRSI = { 600,/* lineNo */
  "isSeparable",                       /* fcnName */
  "/usr/local/MATLAB/R2019b/toolbox/images/images/eml/imfilter.m"/* pathName */
};

static emlrtRSInfo lf_emlrtRSI = { 603,/* lineNo */
  "isSeparable",                       /* fcnName */
  "/usr/local/MATLAB/R2019b/toolbox/images/images/eml/imfilter.m"/* pathName */
};

static emlrtRSInfo mf_emlrtRSI = { 606,/* lineNo */
  "isSeparable",                       /* fcnName */
  "/usr/local/MATLAB/R2019b/toolbox/images/images/eml/imfilter.m"/* pathName */
};

static emlrtRSInfo nf_emlrtRSI = { 607,/* lineNo */
  "isSeparable",                       /* fcnName */
  "/usr/local/MATLAB/R2019b/toolbox/images/images/eml/imfilter.m"/* pathName */
};

static emlrtRSInfo of_emlrtRSI = { 608,/* lineNo */
  "isSeparable",                       /* fcnName */
  "/usr/local/MATLAB/R2019b/toolbox/images/images/eml/imfilter.m"/* pathName */
};

static emlrtRSInfo pg_emlrtRSI = { 769,/* lineNo */
  "padImage",                          /* fcnName */
  "/usr/local/MATLAB/R2019b/toolbox/images/images/eml/imfilter.m"/* pathName */
};

static emlrtRSInfo qg_emlrtRSI = { 80, /* lineNo */
  "padarray",                          /* fcnName */
  "/usr/local/MATLAB/R2019b/toolbox/images/images/eml/padarray.m"/* pathName */
};

static emlrtRSInfo rg_emlrtRSI = { 736,/* lineNo */
  "getPaddingIndices",                 /* fcnName */
  "/usr/local/MATLAB/R2019b/toolbox/images/images/eml/padarray.m"/* pathName */
};

static emlrtRSInfo sg_emlrtRSI = { 823,/* lineNo */
  "filterPartOrWhole",                 /* fcnName */
  "/usr/local/MATLAB/R2019b/toolbox/images/images/eml/imfilter.m"/* pathName */
};

static emlrtECInfo i_emlrtECI = { -1,  /* nDims */
  846,                                 /* lineNo */
  9,                                   /* colNo */
  "ReplicatePad",                      /* fName */
  "/usr/local/MATLAB/R2019b/toolbox/images/images/eml/padarray.m"/* pName */
};

static emlrtBCInfo k_emlrtBCI = { -1,  /* iFirst */
  -1,                                  /* iLast */
  846,                                 /* lineNo */
  16,                                  /* colNo */
  "",                                  /* aName */
  "ReplicatePad",                      /* fName */
  "/usr/local/MATLAB/R2019b/toolbox/images/images/eml/padarray.m",/* pName */
  0                                    /* checkKind */
};

static emlrtBCInfo l_emlrtBCI = { -1,  /* iFirst */
  -1,                                  /* iLast */
  846,                                 /* lineNo */
  14,                                  /* colNo */
  "",                                  /* aName */
  "ReplicatePad",                      /* fName */
  "/usr/local/MATLAB/R2019b/toolbox/images/images/eml/padarray.m",/* pName */
  0                                    /* checkKind */
};

static emlrtDCInfo h_emlrtDCI = { 836, /* lineNo */
  32,                                  /* colNo */
  "ReplicatePad",                      /* fName */
  "/usr/local/MATLAB/R2019b/toolbox/images/images/eml/padarray.m",/* pName */
  4                                    /* checkKind */
};

static emlrtDCInfo i_emlrtDCI = { 830, /* lineNo */
  33,                                  /* colNo */
  "ReplicatePad",                      /* fName */
  "/usr/local/MATLAB/R2019b/toolbox/images/images/eml/padarray.m",/* pName */
  1                                    /* checkKind */
};

static emlrtDCInfo j_emlrtDCI = { 830, /* lineNo */
  33,                                  /* colNo */
  "ReplicatePad",                      /* fName */
  "/usr/local/MATLAB/R2019b/toolbox/images/images/eml/padarray.m",/* pName */
  4                                    /* checkKind */
};

static emlrtBCInfo m_emlrtBCI = { -1,  /* iFirst */
  -1,                                  /* iLast */
  100,                                 /* lineNo */
  30,                                  /* colNo */
  "",                                  /* aName */
  "padarray",                          /* fName */
  "/usr/local/MATLAB/R2019b/toolbox/images/images/eml/padarray.m",/* pName */
  0                                    /* checkKind */
};

static emlrtBCInfo n_emlrtBCI = { -1,  /* iFirst */
  -1,                                  /* iLast */
  100,                                 /* lineNo */
  21,                                  /* colNo */
  "",                                  /* aName */
  "padarray",                          /* fName */
  "/usr/local/MATLAB/R2019b/toolbox/images/images/eml/padarray.m",/* pName */
  0                                    /* checkKind */
};

static emlrtDCInfo k_emlrtDCI = { 83,  /* lineNo */
  56,                                  /* colNo */
  "padarray",                          /* fName */
  "/usr/local/MATLAB/R2019b/toolbox/images/images/eml/padarray.m",/* pName */
  1                                    /* checkKind */
};

static emlrtDCInfo l_emlrtDCI = { 83,  /* lineNo */
  56,                                  /* colNo */
  "padarray",                          /* fName */
  "/usr/local/MATLAB/R2019b/toolbox/images/images/eml/padarray.m",/* pName */
  4                                    /* checkKind */
};

static emlrtBCInfo r_emlrtBCI = { -1,  /* iFirst */
  -1,                                  /* iLast */
  820,                                 /* lineNo */
  21,                                  /* colNo */
  "",                                  /* aName */
  "filterPartOrWhole",                 /* fName */
  "/usr/local/MATLAB/R2019b/toolbox/images/images/eml/imfilter.m",/* pName */
  0                                    /* checkKind */
};

static emlrtBCInfo s_emlrtBCI = { -1,  /* iFirst */
  -1,                                  /* iLast */
  823,                                 /* lineNo */
  23,                                  /* colNo */
  "",                                  /* aName */
  "filterPartOrWhole",                 /* fName */
  "/usr/local/MATLAB/R2019b/toolbox/images/images/eml/imfilter.m",/* pName */
  0                                    /* checkKind */
};

static emlrtBCInfo t_emlrtBCI = { -1,  /* iFirst */
  -1,                                  /* iLast */
  67,                                  /* lineNo */
  20,                                  /* colNo */
  "",                                  /* aName */
  "imfilter",                          /* fName */
  "/usr/local/MATLAB/R2019b/toolbox/images/images/eml/imfilter.m",/* pName */
  0                                    /* checkKind */
};

static emlrtBCInfo u_emlrtBCI = { -1,  /* iFirst */
  -1,                                  /* iLast */
  67,                                  /* lineNo */
  32,                                  /* colNo */
  "",                                  /* aName */
  "imfilter",                          /* fName */
  "/usr/local/MATLAB/R2019b/toolbox/images/images/eml/imfilter.m",/* pName */
  0                                    /* checkKind */
};

static emlrtBCInfo v_emlrtBCI = { -1,  /* iFirst */
  -1,                                  /* iLast */
  68,                                  /* lineNo */
  20,                                  /* colNo */
  "",                                  /* aName */
  "imfilter",                          /* fName */
  "/usr/local/MATLAB/R2019b/toolbox/images/images/eml/imfilter.m",/* pName */
  0                                    /* checkKind */
};

static emlrtBCInfo w_emlrtBCI = { -1,  /* iFirst */
  -1,                                  /* iLast */
  68,                                  /* lineNo */
  33,                                  /* colNo */
  "",                                  /* aName */
  "imfilter",                          /* fName */
  "/usr/local/MATLAB/R2019b/toolbox/images/images/eml/imfilter.m",/* pName */
  0                                    /* checkKind */
};

static emlrtBCInfo x_emlrtBCI = { -1,  /* iFirst */
  -1,                                  /* iLast */
  820,                                 /* lineNo */
  27,                                  /* colNo */
  "",                                  /* aName */
  "filterPartOrWhole",                 /* fName */
  "/usr/local/MATLAB/R2019b/toolbox/images/images/eml/imfilter.m",/* pName */
  0                                    /* checkKind */
};

static emlrtRTEInfo ed_emlrtRTEI = { 736,/* lineNo */
  12,                                  /* colNo */
  "padarray",                          /* fName */
  "/usr/local/MATLAB/R2019b/toolbox/images/images/eml/padarray.m"/* pName */
};

static emlrtRTEInfo fd_emlrtRTEI = { 836,/* lineNo */
  5,                                   /* colNo */
  "padarray",                          /* fName */
  "/usr/local/MATLAB/R2019b/toolbox/images/images/eml/padarray.m"/* pName */
};

static emlrtRTEInfo gd_emlrtRTEI = { 28,/* lineNo */
  9,                                   /* colNo */
  "colon",                             /* fName */
  "/usr/local/MATLAB/R2019b/toolbox/eml/lib/matlab/ops/colon.m"/* pName */
};

static emlrtRTEInfo hd_emlrtRTEI = { 845,/* lineNo */
  9,                                   /* colNo */
  "padarray",                          /* fName */
  "/usr/local/MATLAB/R2019b/toolbox/images/images/eml/padarray.m"/* pName */
};

static emlrtRTEInfo id_emlrtRTEI = { 769,/* lineNo */
  9,                                   /* colNo */
  "imfilter",                          /* fName */
  "/usr/local/MATLAB/R2019b/toolbox/images/images/eml/imfilter.m"/* pName */
};

static emlrtRTEInfo jd_emlrtRTEI = { 80,/* lineNo */
  5,                                   /* colNo */
  "padarray",                          /* fName */
  "/usr/local/MATLAB/R2019b/toolbox/images/images/eml/padarray.m"/* pName */
};

static emlrtRTEInfo kd_emlrtRTEI = { 839,/* lineNo */
  9,                                   /* colNo */
  "padarray",                          /* fName */
  "/usr/local/MATLAB/R2019b/toolbox/images/images/eml/padarray.m"/* pName */
};

static emlrtRTEInfo ld_emlrtRTEI = { 845,/* lineNo */
  30,                                  /* colNo */
  "padarray",                          /* fName */
  "/usr/local/MATLAB/R2019b/toolbox/images/images/eml/padarray.m"/* pName */
};

static emlrtRTEInfo xg_emlrtRTEI = { 37,/* lineNo */
  5,                                   /* colNo */
  "imfilter",                          /* fName */
  "/usr/local/MATLAB/R2019b/toolbox/images/images/eml/imfilter.m"/* pName */
};

static emlrtRTEInfo yg_emlrtRTEI = { 17,/* lineNo */
  13,                                  /* colNo */
  "isinf",                             /* fName */
  "/usr/local/MATLAB/R2019b/toolbox/eml/lib/matlab/elmat/isinf.m"/* pName */
};

static emlrtRTEInfo ah_emlrtRTEI = { 17,/* lineNo */
  13,                                  /* colNo */
  "isnan",                             /* fName */
  "/usr/local/MATLAB/R2019b/toolbox/eml/lib/matlab/elmat/isnan.m"/* pName */
};

static emlrtRTEInfo bh_emlrtRTEI = { 820,/* lineNo */
  27,                                  /* colNo */
  "imfilter",                          /* fName */
  "/usr/local/MATLAB/R2019b/toolbox/images/images/eml/imfilter.m"/* pName */
};

static emlrtRTEInfo ch_emlrtRTEI = { 66,/* lineNo */
  18,                                  /* colNo */
  "imfilter",                          /* fName */
  "/usr/local/MATLAB/R2019b/toolbox/images/images/eml/imfilter.m"/* pName */
};

static emlrtRTEInfo dh_emlrtRTEI = { 110,/* lineNo */
  17,                                  /* colNo */
  "imfilter",                          /* fName */
  "/usr/local/MATLAB/R2019b/toolbox/images/images/eml/imfilter.m"/* pName */
};

static emlrtRTEInfo eh_emlrtRTEI = { 823,/* lineNo */
  9,                                   /* colNo */
  "imfilter",                          /* fName */
  "/usr/local/MATLAB/R2019b/toolbox/images/images/eml/imfilter.m"/* pName */
};

static emlrtRTEInfo fh_emlrtRTEI = { 606,/* lineNo */
  14,                                  /* colNo */
  "imfilter",                          /* fName */
  "/usr/local/MATLAB/R2019b/toolbox/images/images/eml/imfilter.m"/* pName */
};

static emlrtRTEInfo gh_emlrtRTEI = { 814,/* lineNo */
  8,                                   /* colNo */
  "imfilter",                          /* fName */
  "/usr/local/MATLAB/R2019b/toolbox/images/images/eml/imfilter.m"/* pName */
};

static emlrtRTEInfo hh_emlrtRTEI = { 67,/* lineNo */
  9,                                   /* colNo */
  "imfilter",                          /* fName */
  "/usr/local/MATLAB/R2019b/toolbox/images/images/eml/imfilter.m"/* pName */
};

static emlrtRTEInfo ih_emlrtRTEI = { 68,/* lineNo */
  9,                                   /* colNo */
  "imfilter",                          /* fName */
  "/usr/local/MATLAB/R2019b/toolbox/images/images/eml/imfilter.m"/* pName */
};

static emlrtRTEInfo jh_emlrtRTEI = { 608,/* lineNo */
  16,                                  /* colNo */
  "imfilter",                          /* fName */
  "/usr/local/MATLAB/R2019b/toolbox/images/images/eml/imfilter.m"/* pName */
};

static emlrtRTEInfo kh_emlrtRTEI = { 84,/* lineNo */
  13,                                  /* colNo */
  "imfilter",                          /* fName */
  "/usr/local/MATLAB/R2019b/toolbox/images/images/eml/imfilter.m"/* pName */
};

static emlrtRTEInfo lh_emlrtRTEI = { 820,/* lineNo */
  21,                                  /* colNo */
  "imfilter",                          /* fName */
  "/usr/local/MATLAB/R2019b/toolbox/images/images/eml/imfilter.m"/* pName */
};

static emlrtRTEInfo mh_emlrtRTEI = { 88,/* lineNo */
  13,                                  /* colNo */
  "imfilter",                          /* fName */
  "/usr/local/MATLAB/R2019b/toolbox/images/images/eml/imfilter.m"/* pName */
};

static emlrtRTEInfo nh_emlrtRTEI = { 917,/* lineNo */
  28,                                  /* colNo */
  "imfilter",                          /* fName */
  "/usr/local/MATLAB/R2019b/toolbox/images/images/eml/imfilter.m"/* pName */
};

static emlrtRTEInfo oh_emlrtRTEI = { 1,/* lineNo */
  14,                                  /* colNo */
  "imfilter",                          /* fName */
  "/usr/local/MATLAB/R2019b/toolbox/images/images/eml/imfilter.m"/* pName */
};

static emlrtRTEInfo ph_emlrtRTEI = { 603,/* lineNo */
  8,                                   /* colNo */
  "imfilter",                          /* fName */
  "/usr/local/MATLAB/R2019b/toolbox/images/images/eml/imfilter.m"/* pName */
};

static emlrtRTEInfo qh_emlrtRTEI = { 820,/* lineNo */
  9,                                   /* colNo */
  "imfilter",                          /* fName */
  "/usr/local/MATLAB/R2019b/toolbox/images/images/eml/imfilter.m"/* pName */
};

static emlrtRTEInfo rh_emlrtRTEI = { 823,/* lineNo */
  23,                                  /* colNo */
  "imfilter",                          /* fName */
  "/usr/local/MATLAB/R2019b/toolbox/images/images/eml/imfilter.m"/* pName */
};

/* Function Declarations */
static void padImage(const emlrtStack *sp, const emxArray_real_T *a_tmp, const
                     real_T pad[2], emxArray_real_T *a);

/* Function Definitions */
static void padImage(const emlrtStack *sp, const emxArray_real_T *a_tmp, const
                     real_T pad[2], emxArray_real_T *a)
{
  real_T y[2];
  uint32_T varargin_1_idx_0_tmp_tmp;
  real_T varargin_1_idx_0;
  uint32_T varargin_1_idx_1_tmp_tmp;
  real_T varargin_1_idx_1;
  emxArray_int32_T *idxA;
  emxArray_int8_T *onesVector;
  int32_T i;
  int32_T i1;
  int32_T loop_ub;
  emxArray_real_T *b_y;
  emxArray_uint32_T *idxDir;
  int32_T iv[1];
  int32_T b_i;
  int32_T i2;
  int32_T i3;
  int32_T i4;
  int32_T i5;
  emlrtStack st;
  emlrtStack b_st;
  emlrtStack c_st;
  st.prev = sp;
  st.tls = sp->tls;
  b_st.prev = &st;
  b_st.tls = st.tls;
  c_st.prev = &b_st;
  c_st.tls = b_st.tls;
  emlrtHeapReferenceStackEnterFcnR2012b(sp);
  st.site = &pg_emlrtRSI;
  b_st.site = &rb_emlrtRSI;
  b_validateattributes(&b_st, pad);
  if ((a_tmp->size[0] == 0) || (a_tmp->size[1] == 0)) {
    y[0] = static_cast<real_T>(a_tmp->size[0]) + 2.0 * pad[0];
    y[1] = static_cast<real_T>(a_tmp->size[1]) + 2.0 * pad[1];
    b_st.site = &sb_emlrtRSI;
    repmat(&b_st, y, a);
  } else {
    b_st.site = &qg_emlrtRSI;
    c_st.site = &rg_emlrtRSI;
    varargin_1_idx_0_tmp_tmp = static_cast<uint32_T>(a_tmp->size[0]);
    varargin_1_idx_0 = 2.0 * pad[0] + static_cast<real_T>
      (varargin_1_idx_0_tmp_tmp);
    varargin_1_idx_1_tmp_tmp = static_cast<uint32_T>(a_tmp->size[1]);
    varargin_1_idx_1 = 2.0 * pad[1] + static_cast<real_T>
      (varargin_1_idx_1_tmp_tmp);
    if ((varargin_1_idx_0 < varargin_1_idx_1) || (muDoubleScalarIsNaN
         (varargin_1_idx_0) && (!muDoubleScalarIsNaN(varargin_1_idx_1)))) {
      varargin_1_idx_0 = varargin_1_idx_1;
    }

    if (!(varargin_1_idx_0 >= 0.0)) {
      emlrtNonNegativeCheckR2012b(varargin_1_idx_0, &j_emlrtDCI, &c_st);
    }

    if (varargin_1_idx_0 != static_cast<int32_T>(muDoubleScalarFloor
         (varargin_1_idx_0))) {
      emlrtIntegerCheckR2012b(varargin_1_idx_0, &i_emlrtDCI, &c_st);
    }

    emxInit_int32_T(&c_st, &idxA, 2, &jd_emlrtRTEI, true);
    emxInit_int8_T(&c_st, &onesVector, 2, &fd_emlrtRTEI, true);
    i = idxA->size[0] * idxA->size[1];
    i1 = static_cast<int32_T>(varargin_1_idx_0);
    idxA->size[0] = i1;
    idxA->size[1] = 2;
    emxEnsureCapacity_int32_T(&c_st, idxA, i, &ed_emlrtRTEI);
    i = onesVector->size[0] * onesVector->size[1];
    onesVector->size[0] = 1;
    emxEnsureCapacity_int8_T(&c_st, onesVector, i, &fd_emlrtRTEI);
    if (!(pad[0] >= 0.0)) {
      emlrtNonNegativeCheckR2012b(pad[0], &h_emlrtDCI, &c_st);
    }

    loop_ub = static_cast<int32_T>(pad[0]);
    i = onesVector->size[0] * onesVector->size[1];
    onesVector->size[1] = loop_ub;
    emxEnsureCapacity_int8_T(&c_st, onesVector, i, &fd_emlrtRTEI);
    for (i = 0; i < loop_ub; i++) {
      onesVector->data[i] = 1;
    }

    emxInit_real_T(&c_st, &b_y, 2, &ld_emlrtRTEI, true);
    i = b_y->size[0] * b_y->size[1];
    b_y->size[0] = 1;
    loop_ub = static_cast<int32_T>((static_cast<real_T>(varargin_1_idx_0_tmp_tmp)
      - 1.0));
    b_y->size[1] = loop_ub + 1;
    emxEnsureCapacity_real_T(&c_st, b_y, i, &gd_emlrtRTEI);
    for (i = 0; i <= loop_ub; i++) {
      b_y->data[i] = static_cast<real_T>(i) + 1.0;
    }

    emxInit_uint32_T(&c_st, &idxDir, 2, &kd_emlrtRTEI, true);
    i = idxDir->size[0] * idxDir->size[1];
    idxDir->size[0] = 1;
    idxDir->size[1] = (onesVector->size[1] + b_y->size[1]) + onesVector->size[1];
    emxEnsureCapacity_uint32_T(&c_st, idxDir, i, &hd_emlrtRTEI);
    loop_ub = onesVector->size[1];
    for (i = 0; i < loop_ub; i++) {
      idxDir->data[i] = static_cast<uint32_T>(onesVector->data[i]);
    }

    loop_ub = b_y->size[1];
    for (i = 0; i < loop_ub; i++) {
      idxDir->data[i + onesVector->size[1]] = static_cast<uint32_T>
        (muDoubleScalarRound(b_y->data[i]));
    }

    loop_ub = onesVector->size[1];
    for (i = 0; i < loop_ub; i++) {
      idxDir->data[(i + onesVector->size[1]) + b_y->size[1]] =
        varargin_1_idx_0_tmp_tmp * onesVector->data[i];
    }

    if (1 > i1) {
      emlrtDynamicBoundsCheckR2012b(1, 1, i1, &l_emlrtBCI, &c_st);
    }

    if ((idxDir->size[1] < 1) || (idxDir->size[1] > i1)) {
      emlrtDynamicBoundsCheckR2012b(idxDir->size[1], 1, i1, &k_emlrtBCI, &c_st);
    }

    iv[0] = idxDir->size[1];
    emlrtSubAssignSizeCheckR2012b(&iv[0], 1, &idxDir->size[0], 2, &i_emlrtECI,
      &c_st);
    loop_ub = idxDir->size[1];
    for (i = 0; i < loop_ub; i++) {
      idxA->data[i] = static_cast<int32_T>(idxDir->data[i]);
    }

    i = onesVector->size[0] * onesVector->size[1];
    onesVector->size[0] = 1;
    emxEnsureCapacity_int8_T(&c_st, onesVector, i, &fd_emlrtRTEI);
    if (!(pad[1] >= 0.0)) {
      emlrtNonNegativeCheckR2012b(pad[1], &h_emlrtDCI, &c_st);
    }

    loop_ub = static_cast<int32_T>(pad[1]);
    i = onesVector->size[0] * onesVector->size[1];
    onesVector->size[1] = loop_ub;
    emxEnsureCapacity_int8_T(&c_st, onesVector, i, &fd_emlrtRTEI);
    for (i = 0; i < loop_ub; i++) {
      onesVector->data[i] = 1;
    }

    i = b_y->size[0] * b_y->size[1];
    b_y->size[0] = 1;
    loop_ub = static_cast<int32_T>((static_cast<real_T>(varargin_1_idx_1_tmp_tmp)
      - 1.0));
    b_y->size[1] = loop_ub + 1;
    emxEnsureCapacity_real_T(&c_st, b_y, i, &gd_emlrtRTEI);
    for (i = 0; i <= loop_ub; i++) {
      b_y->data[i] = static_cast<real_T>(i) + 1.0;
    }

    i = idxDir->size[0] * idxDir->size[1];
    idxDir->size[0] = 1;
    idxDir->size[1] = (onesVector->size[1] + b_y->size[1]) + onesVector->size[1];
    emxEnsureCapacity_uint32_T(&c_st, idxDir, i, &hd_emlrtRTEI);
    loop_ub = onesVector->size[1];
    for (i = 0; i < loop_ub; i++) {
      idxDir->data[i] = static_cast<uint32_T>(onesVector->data[i]);
    }

    loop_ub = b_y->size[1];
    for (i = 0; i < loop_ub; i++) {
      idxDir->data[i + onesVector->size[1]] = static_cast<uint32_T>
        (muDoubleScalarRound(b_y->data[i]));
    }

    loop_ub = onesVector->size[1];
    for (i = 0; i < loop_ub; i++) {
      idxDir->data[(i + onesVector->size[1]) + b_y->size[1]] =
        varargin_1_idx_1_tmp_tmp * onesVector->data[i];
    }

    emxFree_real_T(&b_y);
    emxFree_int8_T(&onesVector);
    if (1 > idxA->size[0]) {
      emlrtDynamicBoundsCheckR2012b(1, 1, idxA->size[0], &l_emlrtBCI, &c_st);
    }

    if ((idxDir->size[1] < 1) || (idxDir->size[1] > idxA->size[0])) {
      emlrtDynamicBoundsCheckR2012b(idxDir->size[1], 1, idxA->size[0],
        &k_emlrtBCI, &c_st);
    }

    iv[0] = idxDir->size[1];
    emlrtSubAssignSizeCheckR2012b(&iv[0], 1, &idxDir->size[0], 2, &i_emlrtECI,
      &c_st);
    loop_ub = idxDir->size[1];
    for (i = 0; i < loop_ub; i++) {
      idxA->data[i + idxA->size[0]] = static_cast<int32_T>(idxDir->data[i]);
    }

    emxFree_uint32_T(&idxDir);
    varargin_1_idx_0 = static_cast<real_T>(a_tmp->size[0]) + 2.0 * pad[0];
    if (!(varargin_1_idx_0 >= 0.0)) {
      emlrtNonNegativeCheckR2012b(varargin_1_idx_0, &l_emlrtDCI, &st);
    }

    if (varargin_1_idx_0 != static_cast<int32_T>(muDoubleScalarFloor
         (varargin_1_idx_0))) {
      emlrtIntegerCheckR2012b(varargin_1_idx_0, &k_emlrtDCI, &st);
    }

    varargin_1_idx_1 = static_cast<real_T>(a_tmp->size[1]) + 2.0 * pad[1];
    if (!(varargin_1_idx_1 >= 0.0)) {
      emlrtNonNegativeCheckR2012b(varargin_1_idx_1, &l_emlrtDCI, &st);
    }

    if (varargin_1_idx_1 != static_cast<int32_T>(muDoubleScalarFloor
         (varargin_1_idx_1))) {
      emlrtIntegerCheckR2012b(varargin_1_idx_1, &k_emlrtDCI, &st);
    }

    i = a->size[0] * a->size[1];
    a->size[0] = static_cast<int32_T>(varargin_1_idx_0);
    i1 = static_cast<int32_T>(varargin_1_idx_1);
    a->size[1] = i1;
    emxEnsureCapacity_real_T(&st, a, i, &id_emlrtRTEI);
    for (loop_ub = 0; loop_ub < i1; loop_ub++) {
      i = a->size[0];
      for (b_i = 0; b_i < i; b_i++) {
        i2 = b_i + 1;
        if ((i2 < 1) || (i2 > idxA->size[0])) {
          emlrtDynamicBoundsCheckR2012b(i2, 1, idxA->size[0], &m_emlrtBCI, &st);
        }

        i2 = idxA->data[i2 - 1];
        if ((i2 < 1) || (i2 > a_tmp->size[0])) {
          emlrtDynamicBoundsCheckR2012b(i2, 1, a_tmp->size[0], &m_emlrtBCI, &st);
        }

        i3 = loop_ub + 1;
        if ((i3 < 1) || (i3 > idxA->size[0])) {
          emlrtDynamicBoundsCheckR2012b(i3, 1, idxA->size[0], &m_emlrtBCI, &st);
        }

        i3 = idxA->data[(i3 + idxA->size[0]) - 1];
        if ((i3 < 1) || (i3 > a_tmp->size[1])) {
          emlrtDynamicBoundsCheckR2012b(i3, 1, a_tmp->size[1], &m_emlrtBCI, &st);
        }

        i4 = b_i + 1;
        if ((i4 < 1) || (i4 > a->size[0])) {
          emlrtDynamicBoundsCheckR2012b(i4, 1, a->size[0], &n_emlrtBCI, &st);
        }

        i5 = loop_ub + 1;
        if ((i5 < 1) || (i5 > a->size[1])) {
          emlrtDynamicBoundsCheckR2012b(i5, 1, a->size[1], &n_emlrtBCI, &st);
        }

        a->data[(i4 + a->size[0] * (i5 - 1)) - 1] = a_tmp->data[(i2 +
          a_tmp->size[0] * (i3 - 1)) - 1];
      }
    }

    emxFree_int32_T(&idxA);
  }

  emlrtHeapReferenceStackLeaveFcnR2012b(sp);
}

void b_imfilter(const emlrtStack *sp, emxArray_real_T *varargin_1)
{
  real_T outSizeT[2];
  real_T startT[2];
  emxArray_real_T *a;
  boolean_T tooBig;
  int32_T i;
  real_T padSizeT[2];
  real_T nonZeroKernel[2];
  real_T connDimsT[2];
  static const real_T kernel[9] = { 0.0, 1.0, 0.0, 0.0, -1.0, 0.0, 0.0, 0.0, 0.0
  };

  static const boolean_T conn[9] = { false, true, false, false, true, false,
    false, false, false };

  emlrtStack st;
  emlrtStack b_st;
  emlrtStack c_st;
  st.prev = sp;
  st.tls = sp->tls;
  b_st.prev = &st;
  b_st.tls = st.tls;
  c_st.prev = &b_st;
  c_st.tls = b_st.tls;
  emlrtHeapReferenceStackEnterFcnR2012b(sp);
  outSizeT[0] = varargin_1->size[0];
  startT[0] = 1.0;
  outSizeT[1] = varargin_1->size[1];
  startT[1] = 1.0;
  if ((varargin_1->size[0] != 0) && (varargin_1->size[1] != 0)) {
    emxInit_real_T(sp, &a, 2, &vg_emlrtRTEI, true);
    st.site = &wd_emlrtRSI;
    b_st.site = &yd_emlrtRSI;
    padarray(&b_st, varargin_1, startT, a);
    st.site = &xd_emlrtRSI;
    b_st.site = &ae_emlrtRSI;
    tooBig = (outSizeT[0] > 65500.0);
    if ((!tooBig) || (!(outSizeT[1] > 65500.0))) {
      tooBig = false;
    }

    tooBig = !tooBig;
    c_st.site = &be_emlrtRSI;
    i = varargin_1->size[0] * varargin_1->size[1];
    varargin_1->size[0] = static_cast<int32_T>(outSizeT[0]);
    varargin_1->size[1] = static_cast<int32_T>(outSizeT[1]);
    emxEnsureCapacity_real_T(&c_st, varargin_1, i, &yd_emlrtRTEI);
    if (tooBig) {
      padSizeT[0] = a->size[0];
      startT[0] = 3.0;
      padSizeT[1] = a->size[1];
      startT[1] = 3.0;
      ippfilter_real64(&a->data[0], &varargin_1->data[0], outSizeT, 2.0,
                       padSizeT, kernel, startT, true);
    } else {
      padSizeT[0] = a->size[0];
      nonZeroKernel[0] = 1.0;
      connDimsT[0] = 3.0;
      padSizeT[1] = a->size[1];
      nonZeroKernel[1] = -1.0;
      connDimsT[1] = 3.0;
      imfilter_real64(&a->data[0], &varargin_1->data[0], 2.0, outSizeT, 2.0,
                      padSizeT, nonZeroKernel, 2.0, conn, 2.0, connDimsT, startT,
                      2.0, true, true);
    }

    emxFree_real_T(&a);
  }

  emlrtHeapReferenceStackLeaveFcnR2012b(sp);
}

void c_imfilter(const emlrtStack *sp, emxArray_real_T *varargin_1)
{
  real_T outSizeT[2];
  real_T startT[2];
  emxArray_real_T *a;
  boolean_T tooBig;
  int32_T i;
  real_T padSizeT[2];
  real_T nonZeroKernel[2];
  real_T connDimsT[2];
  static const real_T kernel[9] = { 0.0, 0.0, 0.0, 0.0, -1.0, 0.0, 0.0, 1.0, 0.0
  };

  static const boolean_T conn[9] = { false, false, false, false, true, false,
    false, true, false };

  emlrtStack st;
  emlrtStack b_st;
  emlrtStack c_st;
  st.prev = sp;
  st.tls = sp->tls;
  b_st.prev = &st;
  b_st.tls = st.tls;
  c_st.prev = &b_st;
  c_st.tls = b_st.tls;
  emlrtHeapReferenceStackEnterFcnR2012b(sp);
  outSizeT[0] = varargin_1->size[0];
  startT[0] = 1.0;
  outSizeT[1] = varargin_1->size[1];
  startT[1] = 1.0;
  if ((varargin_1->size[0] != 0) && (varargin_1->size[1] != 0)) {
    emxInit_real_T(sp, &a, 2, &vg_emlrtRTEI, true);
    st.site = &wd_emlrtRSI;
    b_st.site = &yd_emlrtRSI;
    padarray(&b_st, varargin_1, startT, a);
    st.site = &xd_emlrtRSI;
    b_st.site = &ae_emlrtRSI;
    tooBig = (outSizeT[0] > 65500.0);
    if ((!tooBig) || (!(outSizeT[1] > 65500.0))) {
      tooBig = false;
    }

    tooBig = !tooBig;
    c_st.site = &be_emlrtRSI;
    i = varargin_1->size[0] * varargin_1->size[1];
    varargin_1->size[0] = static_cast<int32_T>(outSizeT[0]);
    varargin_1->size[1] = static_cast<int32_T>(outSizeT[1]);
    emxEnsureCapacity_real_T(&c_st, varargin_1, i, &yd_emlrtRTEI);
    if (tooBig) {
      padSizeT[0] = a->size[0];
      startT[0] = 3.0;
      padSizeT[1] = a->size[1];
      startT[1] = 3.0;
      ippfilter_real64(&a->data[0], &varargin_1->data[0], outSizeT, 2.0,
                       padSizeT, kernel, startT, true);
    } else {
      padSizeT[0] = a->size[0];
      nonZeroKernel[0] = -1.0;
      connDimsT[0] = 3.0;
      padSizeT[1] = a->size[1];
      nonZeroKernel[1] = 1.0;
      connDimsT[1] = 3.0;
      imfilter_real64(&a->data[0], &varargin_1->data[0], 2.0, outSizeT, 2.0,
                      padSizeT, nonZeroKernel, 2.0, conn, 2.0, connDimsT, startT,
                      2.0, true, true);
    }

    emxFree_real_T(&a);
  }

  emlrtHeapReferenceStackLeaveFcnR2012b(sp);
}

void d_imfilter(const emlrtStack *sp, emxArray_real_T *varargin_1)
{
  real_T outSizeT[2];
  real_T startT[2];
  emxArray_real_T *a;
  boolean_T tooBig;
  int32_T i;
  real_T padSizeT[2];
  real_T nonZeroKernel[2];
  real_T connDimsT[2];
  static const real_T kernel[9] = { 0.0, 0.0, 0.0, 0.0, -1.0, 0.0, 1.0, 0.0, 0.0
  };

  static const boolean_T conn[9] = { false, false, false, false, true, false,
    true, false, false };

  emlrtStack st;
  emlrtStack b_st;
  emlrtStack c_st;
  st.prev = sp;
  st.tls = sp->tls;
  b_st.prev = &st;
  b_st.tls = st.tls;
  c_st.prev = &b_st;
  c_st.tls = b_st.tls;
  emlrtHeapReferenceStackEnterFcnR2012b(sp);
  outSizeT[0] = varargin_1->size[0];
  startT[0] = 1.0;
  outSizeT[1] = varargin_1->size[1];
  startT[1] = 1.0;
  if ((varargin_1->size[0] != 0) && (varargin_1->size[1] != 0)) {
    emxInit_real_T(sp, &a, 2, &vg_emlrtRTEI, true);
    st.site = &wd_emlrtRSI;
    b_st.site = &yd_emlrtRSI;
    padarray(&b_st, varargin_1, startT, a);
    st.site = &xd_emlrtRSI;
    b_st.site = &ae_emlrtRSI;
    tooBig = (outSizeT[0] > 65500.0);
    if ((!tooBig) || (!(outSizeT[1] > 65500.0))) {
      tooBig = false;
    }

    tooBig = !tooBig;
    c_st.site = &be_emlrtRSI;
    i = varargin_1->size[0] * varargin_1->size[1];
    varargin_1->size[0] = static_cast<int32_T>(outSizeT[0]);
    varargin_1->size[1] = static_cast<int32_T>(outSizeT[1]);
    emxEnsureCapacity_real_T(&c_st, varargin_1, i, &yd_emlrtRTEI);
    if (tooBig) {
      padSizeT[0] = a->size[0];
      startT[0] = 3.0;
      padSizeT[1] = a->size[1];
      startT[1] = 3.0;
      ippfilter_real64(&a->data[0], &varargin_1->data[0], outSizeT, 2.0,
                       padSizeT, kernel, startT, true);
    } else {
      padSizeT[0] = a->size[0];
      nonZeroKernel[0] = -1.0;
      connDimsT[0] = 3.0;
      padSizeT[1] = a->size[1];
      nonZeroKernel[1] = 1.0;
      connDimsT[1] = 3.0;
      imfilter_real64(&a->data[0], &varargin_1->data[0], 2.0, outSizeT, 2.0,
                      padSizeT, nonZeroKernel, 2.0, conn, 2.0, connDimsT, startT,
                      2.0, true, true);
    }

    emxFree_real_T(&a);
  }

  emlrtHeapReferenceStackLeaveFcnR2012b(sp);
}

void e_imfilter(const emlrtStack *sp, emxArray_real_T *varargin_1)
{
  real_T outSizeT[2];
  real_T startT[2];
  emxArray_real_T *a;
  boolean_T tooBig;
  int32_T i;
  real_T padSizeT[2];
  real_T nonZeroKernel[2];
  real_T connDimsT[2];
  static const real_T kernel[9] = { 0.0, 0.0, 0.0, 0.0, -1.0, 0.0, 0.0, 0.0, 1.0
  };

  static const boolean_T conn[9] = { false, false, false, false, true, false,
    false, false, true };

  emlrtStack st;
  emlrtStack b_st;
  emlrtStack c_st;
  st.prev = sp;
  st.tls = sp->tls;
  b_st.prev = &st;
  b_st.tls = st.tls;
  c_st.prev = &b_st;
  c_st.tls = b_st.tls;
  emlrtHeapReferenceStackEnterFcnR2012b(sp);
  outSizeT[0] = varargin_1->size[0];
  startT[0] = 1.0;
  outSizeT[1] = varargin_1->size[1];
  startT[1] = 1.0;
  if ((varargin_1->size[0] != 0) && (varargin_1->size[1] != 0)) {
    emxInit_real_T(sp, &a, 2, &vg_emlrtRTEI, true);
    st.site = &wd_emlrtRSI;
    b_st.site = &yd_emlrtRSI;
    padarray(&b_st, varargin_1, startT, a);
    st.site = &xd_emlrtRSI;
    b_st.site = &ae_emlrtRSI;
    tooBig = (outSizeT[0] > 65500.0);
    if ((!tooBig) || (!(outSizeT[1] > 65500.0))) {
      tooBig = false;
    }

    tooBig = !tooBig;
    c_st.site = &be_emlrtRSI;
    i = varargin_1->size[0] * varargin_1->size[1];
    varargin_1->size[0] = static_cast<int32_T>(outSizeT[0]);
    varargin_1->size[1] = static_cast<int32_T>(outSizeT[1]);
    emxEnsureCapacity_real_T(&c_st, varargin_1, i, &yd_emlrtRTEI);
    if (tooBig) {
      padSizeT[0] = a->size[0];
      startT[0] = 3.0;
      padSizeT[1] = a->size[1];
      startT[1] = 3.0;
      ippfilter_real64(&a->data[0], &varargin_1->data[0], outSizeT, 2.0,
                       padSizeT, kernel, startT, true);
    } else {
      padSizeT[0] = a->size[0];
      nonZeroKernel[0] = -1.0;
      connDimsT[0] = 3.0;
      padSizeT[1] = a->size[1];
      nonZeroKernel[1] = 1.0;
      connDimsT[1] = 3.0;
      imfilter_real64(&a->data[0], &varargin_1->data[0], 2.0, outSizeT, 2.0,
                      padSizeT, nonZeroKernel, 2.0, conn, 2.0, connDimsT, startT,
                      2.0, true, true);
    }

    emxFree_real_T(&a);
  }

  emlrtHeapReferenceStackLeaveFcnR2012b(sp);
}

void f_imfilter(const emlrtStack *sp, emxArray_real_T *varargin_1)
{
  real_T outSizeT[2];
  real_T startT[2];
  emxArray_real_T *a;
  boolean_T tooBig;
  int32_T i;
  real_T padSizeT[2];
  real_T nonZeroKernel[2];
  real_T connDimsT[2];
  static const real_T kernel[9] = { 0.0, 0.0, 1.0, 0.0, -1.0, 0.0, 0.0, 0.0, 0.0
  };

  static const boolean_T conn[9] = { false, false, true, false, true, false,
    false, false, false };

  emlrtStack st;
  emlrtStack b_st;
  emlrtStack c_st;
  st.prev = sp;
  st.tls = sp->tls;
  b_st.prev = &st;
  b_st.tls = st.tls;
  c_st.prev = &b_st;
  c_st.tls = b_st.tls;
  emlrtHeapReferenceStackEnterFcnR2012b(sp);
  outSizeT[0] = varargin_1->size[0];
  startT[0] = 1.0;
  outSizeT[1] = varargin_1->size[1];
  startT[1] = 1.0;
  if ((varargin_1->size[0] != 0) && (varargin_1->size[1] != 0)) {
    emxInit_real_T(sp, &a, 2, &vg_emlrtRTEI, true);
    st.site = &wd_emlrtRSI;
    b_st.site = &yd_emlrtRSI;
    padarray(&b_st, varargin_1, startT, a);
    st.site = &xd_emlrtRSI;
    b_st.site = &ae_emlrtRSI;
    tooBig = (outSizeT[0] > 65500.0);
    if ((!tooBig) || (!(outSizeT[1] > 65500.0))) {
      tooBig = false;
    }

    tooBig = !tooBig;
    c_st.site = &be_emlrtRSI;
    i = varargin_1->size[0] * varargin_1->size[1];
    varargin_1->size[0] = static_cast<int32_T>(outSizeT[0]);
    varargin_1->size[1] = static_cast<int32_T>(outSizeT[1]);
    emxEnsureCapacity_real_T(&c_st, varargin_1, i, &yd_emlrtRTEI);
    if (tooBig) {
      padSizeT[0] = a->size[0];
      startT[0] = 3.0;
      padSizeT[1] = a->size[1];
      startT[1] = 3.0;
      ippfilter_real64(&a->data[0], &varargin_1->data[0], outSizeT, 2.0,
                       padSizeT, kernel, startT, true);
    } else {
      padSizeT[0] = a->size[0];
      nonZeroKernel[0] = 1.0;
      connDimsT[0] = 3.0;
      padSizeT[1] = a->size[1];
      nonZeroKernel[1] = -1.0;
      connDimsT[1] = 3.0;
      imfilter_real64(&a->data[0], &varargin_1->data[0], 2.0, outSizeT, 2.0,
                      padSizeT, nonZeroKernel, 2.0, conn, 2.0, connDimsT, startT,
                      2.0, true, true);
    }

    emxFree_real_T(&a);
  }

  emlrtHeapReferenceStackLeaveFcnR2012b(sp);
}

void g_imfilter(const emlrtStack *sp, emxArray_real_T *varargin_1)
{
  real_T outSizeT[2];
  real_T startT[2];
  emxArray_real_T *a;
  boolean_T tooBig;
  int32_T i;
  real_T padSizeT[2];
  real_T nonZeroKernel[2];
  real_T connDimsT[2];
  static const real_T kernel[9] = { 1.0, 0.0, 0.0, 0.0, -1.0, 0.0, 0.0, 0.0, 0.0
  };

  static const boolean_T conn[9] = { true, false, false, false, true, false,
    false, false, false };

  emlrtStack st;
  emlrtStack b_st;
  emlrtStack c_st;
  st.prev = sp;
  st.tls = sp->tls;
  b_st.prev = &st;
  b_st.tls = st.tls;
  c_st.prev = &b_st;
  c_st.tls = b_st.tls;
  emlrtHeapReferenceStackEnterFcnR2012b(sp);
  outSizeT[0] = varargin_1->size[0];
  startT[0] = 1.0;
  outSizeT[1] = varargin_1->size[1];
  startT[1] = 1.0;
  if ((varargin_1->size[0] != 0) && (varargin_1->size[1] != 0)) {
    emxInit_real_T(sp, &a, 2, &vg_emlrtRTEI, true);
    st.site = &wd_emlrtRSI;
    b_st.site = &yd_emlrtRSI;
    padarray(&b_st, varargin_1, startT, a);
    st.site = &xd_emlrtRSI;
    b_st.site = &ae_emlrtRSI;
    tooBig = (outSizeT[0] > 65500.0);
    if ((!tooBig) || (!(outSizeT[1] > 65500.0))) {
      tooBig = false;
    }

    tooBig = !tooBig;
    c_st.site = &be_emlrtRSI;
    i = varargin_1->size[0] * varargin_1->size[1];
    varargin_1->size[0] = static_cast<int32_T>(outSizeT[0]);
    varargin_1->size[1] = static_cast<int32_T>(outSizeT[1]);
    emxEnsureCapacity_real_T(&c_st, varargin_1, i, &yd_emlrtRTEI);
    if (tooBig) {
      padSizeT[0] = a->size[0];
      startT[0] = 3.0;
      padSizeT[1] = a->size[1];
      startT[1] = 3.0;
      ippfilter_real64(&a->data[0], &varargin_1->data[0], outSizeT, 2.0,
                       padSizeT, kernel, startT, true);
    } else {
      padSizeT[0] = a->size[0];
      nonZeroKernel[0] = 1.0;
      connDimsT[0] = 3.0;
      padSizeT[1] = a->size[1];
      nonZeroKernel[1] = -1.0;
      connDimsT[1] = 3.0;
      imfilter_real64(&a->data[0], &varargin_1->data[0], 2.0, outSizeT, 2.0,
                      padSizeT, nonZeroKernel, 2.0, conn, 2.0, connDimsT, startT,
                      2.0, true, true);
    }

    emxFree_real_T(&a);
  }

  emlrtHeapReferenceStackLeaveFcnR2012b(sp);
}

void h_imfilter(const emlrtStack *sp, emxArray_real_T *varargin_1, const
                emxArray_real_T *varargin_2)
{
  real_T outSizeT[2];
  real_T startT[2];
  real_T filter_center[2];
  emxArray_real_T *a;
  int32_T n;
  emxArray_real_T *s;
  int32_T b_a;
  emxArray_real_T *unusedU1;
  emxArray_real_T *b_s;
  emxArray_boolean_T *x;
  int32_T nz;
  emxArray_boolean_T *r;
  emxArray_real_T *c_s;
  boolean_T overflow;
  boolean_T b_varargin_2[2];
  int32_T k;
  emxArray_real_T *nonzero_h;
  emxArray_int32_T *r1;
  emxArray_real_T *v;
  int32_T idx;
  emxArray_boolean_T *connb;
  emxArray_real_T *hcol;
  real_T tol;
  boolean_T tooBig;
  emxArray_real_T *hrow;
  real_T padSizeT[2];
  real_T out_size_row[2];
  real_T start[2];
  boolean_T exitg1;
  real_T d;
  emxArray_boolean_T *b_connb;
  emxArray_int32_T *r2;
  emlrtStack st;
  emlrtStack b_st;
  emlrtStack c_st;
  emlrtStack d_st;
  emlrtStack e_st;
  emlrtStack f_st;
  emlrtStack g_st;
  emlrtStack h_st;
  emlrtStack i_st;
  emlrtStack j_st;
  st.prev = sp;
  st.tls = sp->tls;
  b_st.prev = &st;
  b_st.tls = st.tls;
  c_st.prev = &b_st;
  c_st.tls = b_st.tls;
  d_st.prev = &c_st;
  d_st.tls = c_st.tls;
  e_st.prev = &d_st;
  e_st.tls = d_st.tls;
  f_st.prev = &e_st;
  f_st.tls = e_st.tls;
  g_st.prev = &f_st;
  g_st.tls = f_st.tls;
  h_st.prev = &g_st;
  h_st.tls = g_st.tls;
  i_st.prev = &h_st;
  i_st.tls = h_st.tls;
  j_st.prev = &i_st;
  j_st.tls = i_st.tls;
  emlrtHeapReferenceStackEnterFcnR2012b(sp);
  outSizeT[0] = varargin_1->size[0];
  startT[0] = static_cast<real_T>(varargin_2->size[0]) - muDoubleScalarFloor((
    static_cast<real_T>(varargin_2->size[0]) + 1.0) / 2.0);
  outSizeT[1] = varargin_1->size[1];
  startT[1] = static_cast<real_T>(varargin_2->size[1]) - muDoubleScalarFloor((
    static_cast<real_T>(varargin_2->size[1]) + 1.0) / 2.0);
  if ((varargin_1->size[0] != 0) && (varargin_1->size[1] != 0)) {
    if ((varargin_2->size[0] == 0) || (varargin_2->size[1] == 0)) {
      filter_center[1] = varargin_1->size[1];
      n = varargin_1->size[0];
      b_a = varargin_1->size[0] * varargin_1->size[1];
      varargin_1->size[0] = n;
      emxEnsureCapacity_real_T(sp, varargin_1, b_a, &xg_emlrtRTEI);
      nz = static_cast<int32_T>(filter_center[1]);
      b_a = varargin_1->size[0] * varargin_1->size[1];
      varargin_1->size[1] = static_cast<int32_T>(filter_center[1]);
      emxEnsureCapacity_real_T(sp, varargin_1, b_a, &xg_emlrtRTEI);
      for (b_a = 0; b_a < nz; b_a++) {
        for (k = 0; k < n; k++) {
          varargin_1->data[k + varargin_1->size[0] * b_a] = 0.0;
        }
      }
    } else {
      st.site = &cf_emlrtRSI;
      emxInit_real_T(&st, &a, 2, &vg_emlrtRTEI, true);
      emxInit_real_T(&st, &s, 2, &fh_emlrtRTEI, true);
      emxInit_real_T(&st, &unusedU1, 2, &oh_emlrtRTEI, true);
      emxInit_real_T(&st, &b_s, 1, &ph_emlrtRTEI, true);
      emxInit_boolean_T(&st, &x, 1, &jh_emlrtRTEI, true);
      emxInit_boolean_T(&st, &r, 1, &ah_emlrtRTEI, true);
      emxInit_real_T(&st, &c_s, 2, &fh_emlrtRTEI, true);
      if (varargin_2->size[0] * varargin_2->size[1] >= 49) {
        b_varargin_2[0] = (varargin_2->size[0] != 1);
        b_varargin_2[1] = (varargin_2->size[1] != 1);
        if (all(b_varargin_2)) {
          b_a = x->size[0];
          x->size[0] = varargin_2->size[0] * varargin_2->size[1];
          emxEnsureCapacity_boolean_T(&st, x, b_a, &yg_emlrtRTEI);
          n = varargin_2->size[0] * varargin_2->size[1];
          for (b_a = 0; b_a < n; b_a++) {
            x->data[b_a] = muDoubleScalarIsInf(varargin_2->data[b_a]);
          }

          b_a = r->size[0];
          r->size[0] = varargin_2->size[0] * varargin_2->size[1];
          emxEnsureCapacity_boolean_T(&st, r, b_a, &ah_emlrtRTEI);
          n = varargin_2->size[0] * varargin_2->size[1];
          for (b_a = 0; b_a < n; b_a++) {
            r->data[b_a] = muDoubleScalarIsNaN(varargin_2->data[b_a]);
          }

          n = x->size[0];
          for (b_a = 0; b_a < n; b_a++) {
            x->data[b_a] = ((!x->data[b_a]) && (!r->data[b_a]));
          }

          b_st.site = &kf_emlrtRSI;
          if (b_all(&b_st, x)) {
            b_st.site = &lf_emlrtRSI;
            svd(&b_st, varargin_2, a, s, unusedU1);
            n = s->size[0] - 1;
            nz = s->size[1] - 1;
            b_a = c_s->size[0] * c_s->size[1];
            c_s->size[0] = s->size[0];
            c_s->size[1] = s->size[1];
            emxEnsureCapacity_real_T(&st, c_s, b_a, &fh_emlrtRTEI);
            for (b_a = 0; b_a <= nz; b_a++) {
              for (k = 0; k <= n; k++) {
                c_s->data[k + c_s->size[0] * b_a] = s->data[k + s->size[0] * b_a];
              }
            }

            b_a = s->size[0] * s->size[1];
            s->size[0] = c_s->size[0];
            s->size[1] = c_s->size[1];
            emxEnsureCapacity_real_T(&st, s, b_a, &fh_emlrtRTEI);
            n = c_s->size[1];
            for (b_a = 0; b_a < n; b_a++) {
              nz = c_s->size[0];
              for (k = 0; k < nz; k++) {
                s->data[k + s->size[0] * b_a] = c_s->data[k + c_s->size[0] * b_a];
              }
            }

            b_st.site = &mf_emlrtRSI;
            diag(&b_st, s, b_s);
            if (varargin_2->size[0] > varargin_2->size[1]) {
              nz = varargin_2->size[0];
            } else {
              nz = varargin_2->size[1];
            }

            b_st.site = &nf_emlrtRSI;
            c_st.site = &ib_emlrtRSI;
            d_st.site = &jb_emlrtRSI;
            e_st.site = &kb_emlrtRSI;
            if (b_s->size[0] < 1) {
              emlrtErrorWithMessageIdR2018a(&e_st, &m_emlrtRTEI,
                "Coder:toolbox:eml_min_or_max_varDimZero",
                "Coder:toolbox:eml_min_or_max_varDimZero", 0);
            }

            f_st.site = &lb_emlrtRSI;
            g_st.site = &mb_emlrtRSI;
            n = b_s->size[0];
            if (b_s->size[0] <= 2) {
              if (b_s->size[0] == 1) {
                tol = b_s->data[0];
              } else if ((b_s->data[0] < b_s->data[1]) || (muDoubleScalarIsNaN
                          (b_s->data[0]) && (!muDoubleScalarIsNaN(b_s->data[1]))))
              {
                tol = b_s->data[1];
              } else {
                tol = b_s->data[0];
              }
            } else {
              h_st.site = &ob_emlrtRSI;
              if (!muDoubleScalarIsNaN(b_s->data[0])) {
                idx = 1;
              } else {
                idx = 0;
                i_st.site = &pb_emlrtRSI;
                if (b_s->size[0] > 2147483646) {
                  j_st.site = &eb_emlrtRSI;
                  check_forloop_overflow_error(&j_st);
                }

                k = 2;
                exitg1 = false;
                while ((!exitg1) && (k <= b_s->size[0])) {
                  if (!muDoubleScalarIsNaN(b_s->data[k - 1])) {
                    idx = k;
                    exitg1 = true;
                  } else {
                    k++;
                  }
                }
              }

              if (idx == 0) {
                tol = b_s->data[0];
              } else {
                h_st.site = &nb_emlrtRSI;
                tol = b_s->data[idx - 1];
                b_a = idx + 1;
                i_st.site = &qb_emlrtRSI;
                overflow = ((idx + 1 <= b_s->size[0]) && (b_s->size[0] >
                  2147483646));
                if (overflow) {
                  j_st.site = &eb_emlrtRSI;
                  check_forloop_overflow_error(&j_st);
                }

                for (k = b_a; k <= n; k++) {
                  d = b_s->data[k - 1];
                  if (tol < d) {
                    tol = d;
                  }
                }
              }
            }

            tol = static_cast<real_T>(nz) * tol * 2.2204460492503131E-16;
            b_st.site = &of_emlrtRSI;
            b_a = x->size[0];
            x->size[0] = b_s->size[0];
            emxEnsureCapacity_boolean_T(&b_st, x, b_a, &jh_emlrtRTEI);
            n = b_s->size[0];
            for (b_a = 0; b_a < n; b_a++) {
              x->data[b_a] = (b_s->data[b_a] > tol);
            }

            c_st.site = &ab_emlrtRSI;
            d_st.site = &bb_emlrtRSI;
            nz = combineVectorElements(&d_st, x);
            overflow = (nz == 1);
          } else {
            overflow = false;
          }
        } else {
          overflow = false;
        }
      } else {
        overflow = false;
      }

      emxFree_boolean_T(&r);
      emxInit_real_T(sp, &nonzero_h, 1, &qh_emlrtRTEI, true);
      emxInit_int32_T(sp, &r1, 1, &rh_emlrtRTEI, true);
      if (overflow) {
        emxInit_real_T(sp, &v, 2, &oh_emlrtRTEI, true);
        st.site = &df_emlrtRSI;
        padImage(&st, varargin_1, startT, a);
        st.site = &ef_emlrtRSI;
        svd(&st, varargin_2, unusedU1, s, v);
        n = s->size[0] - 1;
        nz = s->size[1] - 1;
        b_a = c_s->size[0] * c_s->size[1];
        c_s->size[0] = s->size[0];
        c_s->size[1] = s->size[1];
        emxEnsureCapacity_real_T(sp, c_s, b_a, &ch_emlrtRTEI);
        for (b_a = 0; b_a <= nz; b_a++) {
          for (k = 0; k <= n; k++) {
            c_s->data[k + c_s->size[0] * b_a] = s->data[k + s->size[0] * b_a];
          }
        }

        b_a = s->size[0] * s->size[1];
        s->size[0] = c_s->size[0];
        s->size[1] = c_s->size[1];
        emxEnsureCapacity_real_T(sp, s, b_a, &ch_emlrtRTEI);
        n = c_s->size[1];
        for (b_a = 0; b_a < n; b_a++) {
          nz = c_s->size[0];
          for (k = 0; k < nz; k++) {
            s->data[k + s->size[0] * b_a] = c_s->data[k + c_s->size[0] * b_a];
          }
        }

        st.site = &ff_emlrtRSI;
        diag(&st, s, b_s);
        if (1 > unusedU1->size[1]) {
          emlrtDynamicBoundsCheckR2012b(1, 1, unusedU1->size[1], &t_emlrtBCI, sp);
        }

        if (1 > b_s->size[0]) {
          emlrtDynamicBoundsCheckR2012b(1, 1, b_s->size[0], &u_emlrtBCI, sp);
        }

        emxInit_real_T(sp, &hcol, 1, &hh_emlrtRTEI, true);
        tol = b_s->data[0];
        st.site = &gf_emlrtRSI;
        b_sqrt(&st, &tol);
        n = unusedU1->size[0];
        b_a = hcol->size[0];
        hcol->size[0] = unusedU1->size[0];
        emxEnsureCapacity_real_T(sp, hcol, b_a, &hh_emlrtRTEI);
        for (b_a = 0; b_a < n; b_a++) {
          hcol->data[b_a] = unusedU1->data[b_a] * tol;
        }

        if (1 > v->size[1]) {
          emlrtDynamicBoundsCheckR2012b(1, 1, v->size[1], &v_emlrtBCI, sp);
        }

        if (1 > b_s->size[0]) {
          emlrtDynamicBoundsCheckR2012b(1, 1, b_s->size[0], &w_emlrtBCI, sp);
        }

        emxInit_real_T(sp, &hrow, 2, &ih_emlrtRTEI, true);
        tol = b_s->data[0];
        st.site = &hf_emlrtRSI;
        b_sqrt(&st, &tol);
        n = v->size[0];
        b_a = hrow->size[0] * hrow->size[1];
        hrow->size[0] = 1;
        hrow->size[1] = v->size[0];
        emxEnsureCapacity_real_T(sp, hrow, b_a, &ih_emlrtRTEI);
        for (b_a = 0; b_a < n; b_a++) {
          hrow->data[b_a] = v->data[b_a] * tol;
        }

        emxFree_real_T(&v);
        out_size_row[0] = a->size[0];
        out_size_row[1] = static_cast<int32_T>(outSizeT[1]);
        start[0] = 0.0;
        start[1] = static_cast<int32_T>(startT[1]);
        st.site = &if_emlrtRSI;
        b_a = x->size[0];
        x->size[0] = hrow->size[1];
        emxEnsureCapacity_boolean_T(&st, x, b_a, &bh_emlrtRTEI);
        n = hrow->size[1];
        for (b_a = 0; b_a < n; b_a++) {
          x->data[b_a] = (hrow->data[b_a] != 0.0);
        }

        n = x->size[0];
        for (idx = 0; idx < n; idx++) {
          if (x->data[idx]) {
            b_a = idx + 1;
            if ((b_a < 1) || (b_a > hrow->size[1])) {
              emlrtDynamicBoundsCheckR2012b(b_a, 1, hrow->size[1], &x_emlrtBCI,
                &st);
            }
          }
        }

        b_a = x->size[0];
        x->size[0] = hrow->size[1];
        emxEnsureCapacity_boolean_T(&st, x, b_a, &bh_emlrtRTEI);
        n = hrow->size[1];
        for (b_a = 0; b_a < n; b_a++) {
          x->data[b_a] = (hrow->data[b_a] != 0.0);
        }

        n = x->size[0] - 1;
        nz = 0;
        for (idx = 0; idx <= n; idx++) {
          if (x->data[idx]) {
            nz++;
          }
        }

        b_a = r1->size[0];
        r1->size[0] = nz;
        emxEnsureCapacity_int32_T(&st, r1, b_a, &kh_emlrtRTEI);
        nz = 0;
        for (idx = 0; idx <= n; idx++) {
          if (x->data[idx]) {
            r1->data[nz] = idx + 1;
            nz++;
          }
        }

        b_st.site = &ae_emlrtRSI;
        b_a = nonzero_h->size[0];
        nonzero_h->size[0] = r1->size[0];
        emxEnsureCapacity_real_T(&b_st, nonzero_h, b_a, &lh_emlrtRTEI);
        n = r1->size[0];
        for (b_a = 0; b_a < n; b_a++) {
          nonzero_h->data[b_a] = hrow->data[r1->data[b_a] - 1];
        }

        emxInit_boolean_T(&b_st, &b_connb, 2, &gh_emlrtRTEI, true);
        b_a = b_connb->size[0] * b_connb->size[1];
        b_connb->size[0] = 1;
        b_connb->size[1] = hrow->size[1];
        emxEnsureCapacity_boolean_T(&b_st, b_connb, b_a, &gh_emlrtRTEI);
        n = hrow->size[0] * hrow->size[1];
        for (b_a = 0; b_a < n; b_a++) {
          b_connb->data[b_a] = (hrow->data[b_a] != 0.0);
        }

        overflow = false;
        if (static_cast<real_T>(r1->size[0]) / static_cast<real_T>(hrow->size[1])
            > 0.05) {
          overflow = true;
        }

        tooBig = (a->size[0] > 65500);
        if ((!tooBig) || (static_cast<int32_T>(outSizeT[1]) <= 65500)) {
          tooBig = false;
        }

        if (overflow && (!tooBig)) {
          overflow = true;
        } else {
          overflow = false;
        }

        c_st.site = &be_emlrtRSI;
        b_a = varargin_1->size[0] * varargin_1->size[1];
        varargin_1->size[0] = a->size[0];
        varargin_1->size[1] = static_cast<int32_T>(outSizeT[1]);
        emxEnsureCapacity_real_T(&c_st, varargin_1, b_a, &yd_emlrtRTEI);
        if (overflow) {
          d_st.site = &ce_emlrtRSI;
          padSizeT[0] = a->size[0];
          filter_center[0] = hrow->size[0];
          padSizeT[1] = a->size[1];
          filter_center[1] = hrow->size[1];
          ippfilter_real64(&a->data[0], &varargin_1->data[0], out_size_row, 2.0,
                           padSizeT, &hrow->data[0], filter_center, true);
        } else {
          d_st.site = &de_emlrtRSI;
          padSizeT[0] = a->size[0];
          filter_center[0] = b_connb->size[0];
          padSizeT[1] = a->size[1];
          filter_center[1] = b_connb->size[1];
          imfilter_real64(&a->data[0], &varargin_1->data[0], 2.0, out_size_row,
                          2.0, padSizeT, &nonzero_h->data[0], static_cast<real_T>
                          (r1->size[0]), &b_connb->data[0], 2.0, filter_center,
                          start, 2.0, true, true);
        }

        emxFree_boolean_T(&b_connb);
        emxFree_real_T(&hrow);
        start[0] = static_cast<int32_T>(startT[0]);
        start[1] = 0.0;
        st.site = &jf_emlrtRSI;
        b_a = x->size[0];
        x->size[0] = hcol->size[0];
        emxEnsureCapacity_boolean_T(&st, x, b_a, &bh_emlrtRTEI);
        n = hcol->size[0];
        for (b_a = 0; b_a < n; b_a++) {
          x->data[b_a] = (hcol->data[b_a] != 0.0);
        }

        n = x->size[0];
        for (idx = 0; idx < n; idx++) {
          if (x->data[idx]) {
            b_a = idx + 1;
            if ((b_a < 1) || (b_a > hcol->size[0])) {
              emlrtDynamicBoundsCheckR2012b(b_a, 1, hcol->size[0], &x_emlrtBCI,
                &st);
            }
          }
        }

        b_a = x->size[0];
        x->size[0] = hcol->size[0];
        emxEnsureCapacity_boolean_T(&st, x, b_a, &bh_emlrtRTEI);
        n = hcol->size[0];
        for (b_a = 0; b_a < n; b_a++) {
          x->data[b_a] = (hcol->data[b_a] != 0.0);
        }

        n = x->size[0] - 1;
        nz = 0;
        for (idx = 0; idx <= n; idx++) {
          if (x->data[idx]) {
            nz++;
          }
        }

        emxInit_int32_T(&st, &r2, 1, &bh_emlrtRTEI, true);
        b_a = r2->size[0];
        r2->size[0] = nz;
        emxEnsureCapacity_int32_T(&st, r2, b_a, &mh_emlrtRTEI);
        nz = 0;
        for (idx = 0; idx <= n; idx++) {
          if (x->data[idx]) {
            r2->data[nz] = idx + 1;
            nz++;
          }
        }

        b_st.site = &ae_emlrtRSI;
        b_a = nonzero_h->size[0];
        nonzero_h->size[0] = r2->size[0];
        emxEnsureCapacity_real_T(&b_st, nonzero_h, b_a, &lh_emlrtRTEI);
        n = r2->size[0];
        for (b_a = 0; b_a < n; b_a++) {
          nonzero_h->data[b_a] = hcol->data[r2->data[b_a] - 1];
        }

        b_a = x->size[0];
        x->size[0] = hcol->size[0];
        emxEnsureCapacity_boolean_T(&b_st, x, b_a, &gh_emlrtRTEI);
        n = hcol->size[0];
        for (b_a = 0; b_a < n; b_a++) {
          x->data[b_a] = (hcol->data[b_a] != 0.0);
        }

        overflow = false;
        if (static_cast<real_T>(r2->size[0]) / static_cast<real_T>(hcol->size[0])
            > 0.05) {
          overflow = true;
        }

        tooBig = (outSizeT[0] > 65500.0);
        if ((!tooBig) || (static_cast<int32_T>(outSizeT[1]) <= 65500)) {
          tooBig = false;
        }

        if (overflow && (!tooBig)) {
          overflow = true;
        } else {
          overflow = false;
        }

        c_st.site = &be_emlrtRSI;
        b_a = a->size[0] * a->size[1];
        a->size[0] = varargin_1->size[0];
        a->size[1] = varargin_1->size[1];
        emxEnsureCapacity_real_T(&c_st, a, b_a, &nh_emlrtRTEI);
        n = varargin_1->size[0] * varargin_1->size[1];
        for (b_a = 0; b_a < n; b_a++) {
          a->data[b_a] = varargin_1->data[b_a];
        }

        b_a = varargin_1->size[0] * varargin_1->size[1];
        varargin_1->size[0] = static_cast<int32_T>(outSizeT[0]);
        varargin_1->size[1] = static_cast<int32_T>(outSizeT[1]);
        emxEnsureCapacity_real_T(&c_st, varargin_1, b_a, &yd_emlrtRTEI);
        if (overflow) {
          d_st.site = &ce_emlrtRSI;
          padSizeT[0] = a->size[0];
          padSizeT[1] = a->size[1];
          filter_center[0] = hcol->size[0];
          filter_center[1] = 1.0;
          ippfilter_real64(&a->data[0], &varargin_1->data[0], outSizeT, 2.0,
                           padSizeT, &hcol->data[0], filter_center, true);
        } else {
          d_st.site = &de_emlrtRSI;
          padSizeT[0] = a->size[0];
          padSizeT[1] = a->size[1];
          filter_center[0] = x->size[0];
          filter_center[1] = 1.0;
          imfilter_real64(&a->data[0], &varargin_1->data[0], 2.0, outSizeT, 2.0,
                          padSizeT, &nonzero_h->data[0], static_cast<real_T>
                          (r2->size[0]), &x->data[0], 2.0, filter_center, start,
                          2.0, true, true);
        }

        emxFree_int32_T(&r2);
        emxFree_real_T(&hcol);
      } else {
        st.site = &wd_emlrtRSI;
        padImage(&st, varargin_1, startT, a);
        st.site = &xd_emlrtRSI;
        if ((varargin_2->size[0] == 1) || (varargin_2->size[1] == 1)) {
          b_a = x->size[0];
          x->size[0] = varargin_2->size[0] * varargin_2->size[1];
          emxEnsureCapacity_boolean_T(&st, x, b_a, &bh_emlrtRTEI);
          n = varargin_2->size[0] * varargin_2->size[1];
          for (b_a = 0; b_a < n; b_a++) {
            x->data[b_a] = (varargin_2->data[b_a] != 0.0);
          }

          n = x->size[0] - 1;
          nz = 0;
          for (idx = 0; idx <= n; idx++) {
            if (x->data[idx]) {
              nz++;
            }
          }

          b_a = nonzero_h->size[0];
          nonzero_h->size[0] = nz;
          emxEnsureCapacity_real_T(&st, nonzero_h, b_a, &dh_emlrtRTEI);
          nz = 0;
          for (idx = 0; idx <= n; idx++) {
            if (x->data[idx]) {
              b_a = varargin_2->size[0] * varargin_2->size[1];
              k = idx + 1;
              if ((k < 1) || (k > b_a)) {
                emlrtDynamicBoundsCheckR2012b(k, 1, b_a, &r_emlrtBCI, &st);
              }

              nonzero_h->data[nz] = varargin_2->data[k - 1];
              nz++;
            }
          }
        } else {
          b_st.site = &sg_emlrtRSI;
          indexShapeCheck(&b_st, *(int32_T (*)[2])varargin_2->size, *(int32_T (*)
            [2])varargin_2->size);
          n = varargin_2->size[0] * varargin_2->size[1] - 1;
          nz = 0;
          for (idx = 0; idx <= n; idx++) {
            if (varargin_2->data[idx] != 0.0) {
              nz++;
            }
          }

          b_a = r1->size[0];
          r1->size[0] = nz;
          emxEnsureCapacity_int32_T(&st, r1, b_a, &dh_emlrtRTEI);
          nz = 0;
          for (idx = 0; idx <= n; idx++) {
            if (varargin_2->data[idx] != 0.0) {
              r1->data[nz] = idx + 1;
              nz++;
            }
          }

          nz = varargin_2->size[0] * varargin_2->size[1];
          b_a = nonzero_h->size[0];
          nonzero_h->size[0] = r1->size[0];
          emxEnsureCapacity_real_T(&st, nonzero_h, b_a, &eh_emlrtRTEI);
          n = r1->size[0];
          for (b_a = 0; b_a < n; b_a++) {
            if ((r1->data[b_a] < 1) || (r1->data[b_a] > nz)) {
              emlrtDynamicBoundsCheckR2012b(r1->data[b_a], 1, nz, &s_emlrtBCI,
                &st);
            }

            nonzero_h->data[b_a] = varargin_2->data[r1->data[b_a] - 1];
          }
        }

        emxInit_boolean_T(&st, &connb, 2, &gh_emlrtRTEI, true);
        b_st.site = &ae_emlrtRSI;
        b_a = connb->size[0] * connb->size[1];
        connb->size[0] = varargin_2->size[0];
        connb->size[1] = varargin_2->size[1];
        emxEnsureCapacity_boolean_T(&b_st, connb, b_a, &gh_emlrtRTEI);
        n = varargin_2->size[0] * varargin_2->size[1];
        for (b_a = 0; b_a < n; b_a++) {
          connb->data[b_a] = (varargin_2->data[b_a] != 0.0);
        }

        overflow = false;
        if (static_cast<real_T>(nonzero_h->size[0]) / static_cast<real_T>
            ((varargin_2->size[0] * varargin_2->size[1])) > 0.05) {
          overflow = true;
        }

        tooBig = (outSizeT[0] > 65500.0);
        if ((!tooBig) || (!(outSizeT[1] > 65500.0))) {
          tooBig = false;
        }

        if (overflow && (!tooBig)) {
          overflow = true;
        } else {
          overflow = false;
        }

        c_st.site = &be_emlrtRSI;
        b_a = varargin_1->size[0] * varargin_1->size[1];
        varargin_1->size[0] = static_cast<int32_T>(outSizeT[0]);
        varargin_1->size[1] = static_cast<int32_T>(outSizeT[1]);
        emxEnsureCapacity_real_T(&c_st, varargin_1, b_a, &yd_emlrtRTEI);
        if (overflow) {
          d_st.site = &ce_emlrtRSI;
          padSizeT[0] = a->size[0];
          filter_center[0] = varargin_2->size[0];
          padSizeT[1] = a->size[1];
          filter_center[1] = varargin_2->size[1];
          ippfilter_real64(&a->data[0], &varargin_1->data[0], outSizeT, 2.0,
                           padSizeT, &varargin_2->data[0], filter_center, true);
        } else {
          d_st.site = &de_emlrtRSI;
          padSizeT[0] = a->size[0];
          filter_center[0] = connb->size[0];
          padSizeT[1] = a->size[1];
          filter_center[1] = connb->size[1];
          imfilter_real64(&a->data[0], &varargin_1->data[0], 2.0, outSizeT, 2.0,
                          padSizeT, &nonzero_h->data[0], static_cast<real_T>
                          (nonzero_h->size[0]), &connb->data[0], 2.0,
                          filter_center, startT, 2.0, true, true);
        }

        emxFree_boolean_T(&connb);
      }

      emxFree_real_T(&c_s);
      emxFree_boolean_T(&x);
      emxFree_int32_T(&r1);
      emxFree_real_T(&nonzero_h);
      emxFree_real_T(&b_s);
      emxFree_real_T(&unusedU1);
      emxFree_real_T(&s);
      emxFree_real_T(&a);
    }
  }

  emlrtHeapReferenceStackLeaveFcnR2012b(sp);
}

void imfilter(const emlrtStack *sp, emxArray_real_T *varargin_1)
{
  real_T outSizeT[2];
  real_T startT[2];
  emxArray_real_T *a;
  boolean_T tooBig;
  int32_T i;
  real_T padSizeT[2];
  real_T nonZeroKernel[2];
  real_T connDimsT[2];
  static const real_T kernel[9] = { 0.0, 0.0, 0.0, 0.0, -1.0, 1.0, 0.0, 0.0, 0.0
  };

  static const boolean_T conn[9] = { false, false, false, false, true, true,
    false, false, false };

  emlrtStack st;
  emlrtStack b_st;
  emlrtStack c_st;
  st.prev = sp;
  st.tls = sp->tls;
  b_st.prev = &st;
  b_st.tls = st.tls;
  c_st.prev = &b_st;
  c_st.tls = b_st.tls;
  emlrtHeapReferenceStackEnterFcnR2012b(sp);
  outSizeT[0] = varargin_1->size[0];
  startT[0] = 1.0;
  outSizeT[1] = varargin_1->size[1];
  startT[1] = 1.0;
  if ((varargin_1->size[0] != 0) && (varargin_1->size[1] != 0)) {
    emxInit_real_T(sp, &a, 2, &vg_emlrtRTEI, true);
    st.site = &wd_emlrtRSI;
    b_st.site = &yd_emlrtRSI;
    padarray(&b_st, varargin_1, startT, a);
    st.site = &xd_emlrtRSI;
    b_st.site = &ae_emlrtRSI;
    tooBig = (outSizeT[0] > 65500.0);
    if ((!tooBig) || (!(outSizeT[1] > 65500.0))) {
      tooBig = false;
    }

    tooBig = !tooBig;
    c_st.site = &be_emlrtRSI;
    i = varargin_1->size[0] * varargin_1->size[1];
    varargin_1->size[0] = static_cast<int32_T>(outSizeT[0]);
    varargin_1->size[1] = static_cast<int32_T>(outSizeT[1]);
    emxEnsureCapacity_real_T(&c_st, varargin_1, i, &yd_emlrtRTEI);
    if (tooBig) {
      padSizeT[0] = a->size[0];
      startT[0] = 3.0;
      padSizeT[1] = a->size[1];
      startT[1] = 3.0;
      ippfilter_real64(&a->data[0], &varargin_1->data[0], outSizeT, 2.0,
                       padSizeT, kernel, startT, true);
    } else {
      padSizeT[0] = a->size[0];
      nonZeroKernel[0] = -1.0;
      connDimsT[0] = 3.0;
      padSizeT[1] = a->size[1];
      nonZeroKernel[1] = 1.0;
      connDimsT[1] = 3.0;
      imfilter_real64(&a->data[0], &varargin_1->data[0], 2.0, outSizeT, 2.0,
                      padSizeT, nonZeroKernel, 2.0, conn, 2.0, connDimsT, startT,
                      2.0, true, true);
    }

    emxFree_real_T(&a);
  }

  emlrtHeapReferenceStackLeaveFcnR2012b(sp);
}

/* End of code generation (imfilter.cpp) */
