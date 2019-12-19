/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * svd1.cpp
 *
 * Code generation for function 'svd1'
 *
 */

/* Include files */
#include "svd1.h"
#include "diff_gradients.h"
#include "diff_gradients_emxutil.h"
#include "lapacke.h"
#include "mwmathutil.h"
#include "rt_nonfinite.h"

/* Variable Definitions */
static emlrtRSInfo xf_emlrtRSI = { 53, /* lineNo */
  "svd",                               /* fcnName */
  "/usr/local/MATLAB/R2019b/toolbox/eml/eml/+coder/+internal/svd.m"/* pathName */
};

static emlrtRSInfo yf_emlrtRSI = { 83, /* lineNo */
  "callLAPACK",                        /* fcnName */
  "/usr/local/MATLAB/R2019b/toolbox/eml/eml/+coder/+internal/svd.m"/* pathName */
};

static emlrtRSInfo ag_emlrtRSI = { 75, /* lineNo */
  "callLAPACK",                        /* fcnName */
  "/usr/local/MATLAB/R2019b/toolbox/eml/eml/+coder/+internal/svd.m"/* pathName */
};

static emlrtRSInfo bg_emlrtRSI = { 209,/* lineNo */
  "xgesdd",                            /* fcnName */
  "/usr/local/MATLAB/R2019b/toolbox/eml/eml/+coder/+internal/+lapack/xgesdd.m"/* pathName */
};

static emlrtRSInfo cg_emlrtRSI = { 179,/* lineNo */
  "xgesdd",                            /* fcnName */
  "/usr/local/MATLAB/R2019b/toolbox/eml/eml/+coder/+internal/+lapack/xgesdd.m"/* pathName */
};

static emlrtRSInfo dg_emlrtRSI = { 68, /* lineNo */
  "xgesdd",                            /* fcnName */
  "/usr/local/MATLAB/R2019b/toolbox/eml/eml/+coder/+internal/+lapack/xgesdd.m"/* pathName */
};

static emlrtRSInfo eg_emlrtRSI = { 61, /* lineNo */
  "xgesdd",                            /* fcnName */
  "/usr/local/MATLAB/R2019b/toolbox/eml/eml/+coder/+internal/+lapack/xgesdd.m"/* pathName */
};

static emlrtRSInfo fg_emlrtRSI = { 58, /* lineNo */
  "xgesdd",                            /* fcnName */
  "/usr/local/MATLAB/R2019b/toolbox/eml/eml/+coder/+internal/+lapack/xgesdd.m"/* pathName */
};

static emlrtRSInfo jg_emlrtRSI = { 31, /* lineNo */
  "xgesvd",                            /* fcnName */
  "/usr/local/MATLAB/R2019b/toolbox/eml/eml/+coder/+internal/+lapack/xgesvd.m"/* pathName */
};

static emlrtRSInfo kg_emlrtRSI = { 196,/* lineNo */
  "ceval_xgesvd",                      /* fcnName */
  "/usr/local/MATLAB/R2019b/toolbox/eml/eml/+coder/+internal/+lapack/xgesvd.m"/* pathName */
};

static emlrtRTEInfo u_emlrtRTEI = { 105,/* lineNo */
  5,                                   /* colNo */
  "callLAPACK",                        /* fName */
  "/usr/local/MATLAB/R2019b/toolbox/eml/eml/+coder/+internal/svd.m"/* pName */
};

static emlrtRTEInfo v_emlrtRTEI = { 45,/* lineNo */
  13,                                  /* colNo */
  "infocheck",                         /* fName */
  "/usr/local/MATLAB/R2019b/toolbox/eml/eml/+coder/+internal/+lapack/infocheck.m"/* pName */
};

static emlrtRTEInfo w_emlrtRTEI = { 48,/* lineNo */
  13,                                  /* colNo */
  "infocheck",                         /* fName */
  "/usr/local/MATLAB/R2019b/toolbox/eml/eml/+coder/+internal/+lapack/infocheck.m"/* pName */
};

static emlrtRTEInfo uc_emlrtRTEI = { 75,/* lineNo */
  63,                                  /* colNo */
  "svd",                               /* fName */
  "/usr/local/MATLAB/R2019b/toolbox/eml/eml/+coder/+internal/svd.m"/* pName */
};

static emlrtRTEInfo vc_emlrtRTEI = { 75,/* lineNo */
  9,                                   /* colNo */
  "svd",                               /* fName */
  "/usr/local/MATLAB/R2019b/toolbox/eml/eml/+coder/+internal/svd.m"/* pName */
};

static emlrtRTEInfo wc_emlrtRTEI = { 217,/* lineNo */
  5,                                   /* colNo */
  "xgesvd",                            /* fName */
  "/usr/local/MATLAB/R2019b/toolbox/eml/eml/+coder/+internal/+lapack/xgesvd.m"/* pName */
};

static emlrtRTEInfo xc_emlrtRTEI = { 31,/* lineNo */
  33,                                  /* colNo */
  "xgesvd",                            /* fName */
  "/usr/local/MATLAB/R2019b/toolbox/eml/eml/+coder/+internal/+lapack/xgesvd.m"/* pName */
};

static emlrtRTEInfo yc_emlrtRTEI = { 31,/* lineNo */
  5,                                   /* colNo */
  "xgesvd",                            /* fName */
  "/usr/local/MATLAB/R2019b/toolbox/eml/eml/+coder/+internal/+lapack/xgesvd.m"/* pName */
};

static emlrtRTEInfo ad_emlrtRTEI = { 75,/* lineNo */
  14,                                  /* colNo */
  "svd",                               /* fName */
  "/usr/local/MATLAB/R2019b/toolbox/eml/eml/+coder/+internal/svd.m"/* pName */
};

static emlrtRTEInfo bd_emlrtRTEI = { 120,/* lineNo */
  9,                                   /* colNo */
  "xgesvd",                            /* fName */
  "/usr/local/MATLAB/R2019b/toolbox/eml/eml/+coder/+internal/+lapack/xgesvd.m"/* pName */
};

/* Function Definitions */
void b_svd(const emlrtStack *sp, const emxArray_real_T *A, emxArray_real_T *U,
           emxArray_real_T *s, emxArray_real_T *V)
{
  emxArray_real_T *b_A;
  int32_T m;
  int32_T n;
  int32_T i;
  int32_T loop_ub;
  emxArray_real_T *Vt;
  int32_T i1;
  ptrdiff_t info_t;
  static const char_T fname[14] = { 'L', 'A', 'P', 'A', 'C', 'K', 'E', '_', 'd',
    'g', 'e', 's', 'd', 'd' };

  emxArray_real_T *superb;
  static const char_T b_fname[14] = { 'L', 'A', 'P', 'A', 'C', 'K', 'E', '_',
    'd', 'g', 'e', 's', 'v', 'd' };

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
  emxInit_real_T(sp, &b_A, 2, &uc_emlrtRTEI, true);
  st.site = &xf_emlrtRSI;
  m = A->size[0];
  n = A->size[1];
  b_st.site = &ag_emlrtRSI;
  i = b_A->size[0] * b_A->size[1];
  b_A->size[0] = A->size[0];
  b_A->size[1] = A->size[1];
  emxEnsureCapacity_real_T(&b_st, b_A, i, &uc_emlrtRTEI);
  loop_ub = A->size[0] * A->size[1];
  for (i = 0; i < loop_ub; i++) {
    b_A->data[i] = A->data[i];
  }

  emxInit_real_T(&b_st, &Vt, 2, &ad_emlrtRTEI, true);
  i = U->size[0] * U->size[1];
  U->size[0] = A->size[0];
  U->size[1] = A->size[0];
  emxEnsureCapacity_real_T(&b_st, U, i, &vc_emlrtRTEI);
  i = Vt->size[0] * Vt->size[1];
  Vt->size[0] = A->size[1];
  Vt->size[1] = A->size[1];
  emxEnsureCapacity_real_T(&b_st, Vt, i, &vc_emlrtRTEI);
  i = s->size[0];
  s->size[0] = muIntScalarMin_sint32(n, m);
  emxEnsureCapacity_real_T(&b_st, s, i, &vc_emlrtRTEI);
  c_st.site = &fg_emlrtRSI;
  if ((A->size[0] != 0) && (A->size[1] != 0)) {
    c_st.site = &eg_emlrtRSI;
    c_st.site = &dg_emlrtRSI;
    c_st.site = &cg_emlrtRSI;
    info_t = LAPACKE_dgesdd(102, 'A', (ptrdiff_t)A->size[0], (ptrdiff_t)A->size
      [1], &b_A->data[0], (ptrdiff_t)A->size[0], &s->data[0], &U->data[0],
      (ptrdiff_t)A->size[0], &Vt->data[0], (ptrdiff_t)A->size[1]);
    m = (int32_T)info_t;
    c_st.site = &bg_emlrtRSI;
    if (m < 0) {
      if (m == -1010) {
        emlrtErrorWithMessageIdR2018a(&c_st, &v_emlrtRTEI, "MATLAB:nomem",
          "MATLAB:nomem", 0);
      } else {
        emlrtErrorWithMessageIdR2018a(&c_st, &w_emlrtRTEI,
          "Coder:toolbox:LAPACKCallErrorInfo",
          "Coder:toolbox:LAPACKCallErrorInfo", 5, 4, 14, fname, 12, m);
      }
    }
  } else {
    m = 0;
  }

  if (m > 0) {
    b_st.site = &yf_emlrtRSI;
    c_st.site = &jg_emlrtRSI;
    i = b_A->size[0] * b_A->size[1];
    b_A->size[0] = A->size[0];
    b_A->size[1] = A->size[1];
    emxEnsureCapacity_real_T(&c_st, b_A, i, &xc_emlrtRTEI);
    loop_ub = A->size[0] * A->size[1];
    for (i = 0; i < loop_ub; i++) {
      b_A->data[i] = A->data[i];
    }

    m = A->size[0];
    n = A->size[1];
    m = muIntScalarMin_sint32(n, m);
    i = U->size[0] * U->size[1];
    U->size[0] = A->size[0];
    U->size[1] = A->size[0];
    emxEnsureCapacity_real_T(&c_st, U, i, &yc_emlrtRTEI);
    i = Vt->size[0] * Vt->size[1];
    Vt->size[0] = A->size[1];
    Vt->size[1] = A->size[1];
    emxEnsureCapacity_real_T(&c_st, Vt, i, &yc_emlrtRTEI);
    i = s->size[0];
    s->size[0] = m;
    emxEnsureCapacity_real_T(&c_st, s, i, &yc_emlrtRTEI);
    if ((A->size[0] != 0) && (A->size[1] != 0)) {
      emxInit_real_T(&c_st, &superb, 1, &bd_emlrtRTEI, true);
      if (m > 1) {
        i = superb->size[0];
        superb->size[0] = m - 1;
        emxEnsureCapacity_real_T(&c_st, superb, i, &yc_emlrtRTEI);
      } else {
        i = superb->size[0];
        superb->size[0] = 1;
        emxEnsureCapacity_real_T(&c_st, superb, i, &yc_emlrtRTEI);
      }

      info_t = LAPACKE_dgesvd(102, 'A', 'A', (ptrdiff_t)A->size[0], (ptrdiff_t)
        A->size[1], &b_A->data[0], (ptrdiff_t)A->size[0], &s->data[0], &U->data
        [0], (ptrdiff_t)A->size[0], &Vt->data[0], (ptrdiff_t)A->size[1],
        &superb->data[0]);
      m = (int32_T)info_t;
      emxFree_real_T(&superb);
    } else {
      m = 0;
    }

    i = V->size[0] * V->size[1];
    V->size[0] = Vt->size[1];
    V->size[1] = Vt->size[0];
    emxEnsureCapacity_real_T(&c_st, V, i, &wc_emlrtRTEI);
    loop_ub = Vt->size[0];
    for (i = 0; i < loop_ub; i++) {
      n = Vt->size[1];
      for (i1 = 0; i1 < n; i1++) {
        V->data[i1 + V->size[0] * i] = Vt->data[i + Vt->size[0] * i1];
      }
    }

    d_st.site = &kg_emlrtRSI;
    if (m < 0) {
      if (m == -1010) {
        emlrtErrorWithMessageIdR2018a(&d_st, &v_emlrtRTEI, "MATLAB:nomem",
          "MATLAB:nomem", 0);
      } else {
        emlrtErrorWithMessageIdR2018a(&d_st, &w_emlrtRTEI,
          "Coder:toolbox:LAPACKCallErrorInfo",
          "Coder:toolbox:LAPACKCallErrorInfo", 5, 4, 14, b_fname, 12, m);
      }
    }
  } else {
    i = V->size[0] * V->size[1];
    V->size[0] = Vt->size[1];
    V->size[1] = Vt->size[0];
    emxEnsureCapacity_real_T(&st, V, i, &wc_emlrtRTEI);
    loop_ub = Vt->size[0];
    for (i = 0; i < loop_ub; i++) {
      n = Vt->size[1];
      for (i1 = 0; i1 < n; i1++) {
        V->data[i1 + V->size[0] * i] = Vt->data[i + Vt->size[0] * i1];
      }
    }
  }

  emxFree_real_T(&b_A);
  emxFree_real_T(&Vt);
  if (m > 0) {
    emlrtErrorWithMessageIdR2018a(&st, &u_emlrtRTEI,
      "Coder:MATLAB:svd_NoConvergence", "Coder:MATLAB:svd_NoConvergence", 0);
  }

  emlrtHeapReferenceStackLeaveFcnR2012b(sp);
}

/* End of code generation (svd1.cpp) */
