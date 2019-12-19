/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * gaussgradient.cpp
 *
 * Code generation for function 'gaussgradient'
 *
 */

/* Include files */
#include "gaussgradient.h"
#include "abs.h"
#include "diff_gradients.h"
#include "diff_gradients_data.h"
#include "diff_gradients_emxutil.h"
#include "eml_int_forloop_overflow_check.h"
#include "imfilter.h"
#include "mwmathutil.h"
#include "rt_nonfinite.h"

/* Variable Definitions */
static emlrtRSInfo ec_emlrtRSI = { 45, /* lineNo */
  "mpower",                            /* fcnName */
  "/usr/local/MATLAB/R2019b/toolbox/eml/lib/matlab/ops/mpower.m"/* pathName */
};

static emlrtRSInfo se_emlrtRSI = { 13, /* lineNo */
  "gaussgradient",                     /* fcnName */
  "/media/data/Dropbox/code_bank/matlab/3Dnucleisegmentation/gaussgradient.m"/* pathName */
};

static emlrtRSInfo te_emlrtRSI = { 20, /* lineNo */
  "gaussgradient",                     /* fcnName */
  "/media/data/Dropbox/code_bank/matlab/3Dnucleisegmentation/gaussgradient.m"/* pathName */
};

static emlrtRSInfo ue_emlrtRSI = { 23, /* lineNo */
  "gaussgradient",                     /* fcnName */
  "/media/data/Dropbox/code_bank/matlab/3Dnucleisegmentation/gaussgradient.m"/* pathName */
};

static emlrtRSInfo ve_emlrtRSI = { 27, /* lineNo */
  "gaussgradient",                     /* fcnName */
  "/media/data/Dropbox/code_bank/matlab/3Dnucleisegmentation/gaussgradient.m"/* pathName */
};

static emlrtRSInfo we_emlrtRSI = { 28, /* lineNo */
  "gaussgradient",                     /* fcnName */
  "/media/data/Dropbox/code_bank/matlab/3Dnucleisegmentation/gaussgradient.m"/* pathName */
};

static emlrtRSInfo ye_emlrtRSI = { 35, /* lineNo */
  "dgauss",                            /* fcnName */
  "/media/data/Dropbox/code_bank/matlab/3Dnucleisegmentation/gaussgradient.m"/* pathName */
};

static emlrtRSInfo af_emlrtRSI = { 20, /* lineNo */
  "sum",                               /* fcnName */
  "/usr/local/MATLAB/R2019b/toolbox/eml/lib/matlab/datafun/sum.m"/* pathName */
};

static emlrtRSInfo bf_emlrtRSI = { 163,/* lineNo */
  "colMajorFlatIter",                  /* fcnName */
  "/usr/local/MATLAB/R2019b/toolbox/eml/lib/matlab/datafun/private/combineVectorElements.m"/* pathName */
};

static emlrtRTEInfo n_emlrtRTEI = { 17,/* lineNo */
  7,                                   /* colNo */
  "gaussgradient",                     /* fName */
  "/media/data/Dropbox/code_bank/matlab/3Dnucleisegmentation/gaussgradient.m"/* pName */
};

static emlrtRTEInfo o_emlrtRTEI = { 18,/* lineNo */
  11,                                  /* colNo */
  "gaussgradient",                     /* fName */
  "/media/data/Dropbox/code_bank/matlab/3Dnucleisegmentation/gaussgradient.m"/* pName */
};

static emlrtECInfo h_emlrtECI = { 2,   /* nDims */
  23,                                  /* lineNo */
  20,                                  /* colNo */
  "gaussgradient",                     /* fName */
  "/media/data/Dropbox/code_bank/matlab/3Dnucleisegmentation/gaussgradient.m"/* pName */
};

static emlrtRTEInfo p_emlrtRTEI = { 14,/* lineNo */
  9,                                   /* colNo */
  "log",                               /* fName */
  "/usr/local/MATLAB/R2019b/toolbox/eml/lib/matlab/elfun/log.m"/* pName */
};

static emlrtRTEInfo r_emlrtRTEI = { 43,/* lineNo */
  23,                                  /* colNo */
  "sumprod",                           /* fName */
  "/usr/local/MATLAB/R2019b/toolbox/eml/lib/matlab/datafun/private/sumprod.m"/* pName */
};

static emlrtRTEInfo s_emlrtRTEI = { 73,/* lineNo */
  9,                                   /* colNo */
  "sumprod",                           /* fName */
  "/usr/local/MATLAB/R2019b/toolbox/eml/lib/matlab/datafun/private/sumprod.m"/* pName */
};

static emlrtDCInfo e_emlrtDCI = { 16,  /* lineNo */
  10,                                  /* colNo */
  "gaussgradient",                     /* fName */
  "/media/data/Dropbox/code_bank/matlab/3Dnucleisegmentation/gaussgradient.m",/* pName */
  1                                    /* checkKind */
};

static emlrtDCInfo f_emlrtDCI = { 16,  /* lineNo */
  10,                                  /* colNo */
  "gaussgradient",                     /* fName */
  "/media/data/Dropbox/code_bank/matlab/3Dnucleisegmentation/gaussgradient.m",/* pName */
  4                                    /* checkKind */
};

static emlrtDCInfo g_emlrtDCI = { 16,  /* lineNo */
  15,                                  /* colNo */
  "gaussgradient",                     /* fName */
  "/media/data/Dropbox/code_bank/matlab/3Dnucleisegmentation/gaussgradient.m",/* pName */
  1                                    /* checkKind */
};

static emlrtBCInfo j_emlrtBCI = { -1,  /* iFirst */
  -1,                                  /* iLast */
  20,                                  /* lineNo */
  9,                                   /* colNo */
  "hx",                                /* aName */
  "gaussgradient",                     /* fName */
  "/media/data/Dropbox/code_bank/matlab/3Dnucleisegmentation/gaussgradient.m",/* pName */
  0                                    /* checkKind */
};

static emlrtRTEInfo ec_emlrtRTEI = { 16,/* lineNo */
  4,                                   /* colNo */
  "gaussgradient",                     /* fName */
  "/media/data/Dropbox/code_bank/matlab/3Dnucleisegmentation/gaussgradient.m"/* pName */
};

static emlrtRTEInfo fc_emlrtRTEI = { 20,/* lineNo */
  1,                                   /* colNo */
  "sum",                               /* fName */
  "/usr/local/MATLAB/R2019b/toolbox/eml/lib/matlab/datafun/sum.m"/* pName */
};

static emlrtRTEInfo gc_emlrtRTEI = { 124,/* lineNo */
  13,                                  /* colNo */
  "combineVectorElements",             /* fName */
  "/usr/local/MATLAB/R2019b/toolbox/eml/lib/matlab/datafun/private/combineVectorElements.m"/* pName */
};

static emlrtRTEInfo hc_emlrtRTEI = { 27,/* lineNo */
  1,                                   /* colNo */
  "gaussgradient",                     /* fName */
  "/media/data/Dropbox/code_bank/matlab/3Dnucleisegmentation/gaussgradient.m"/* pName */
};

static emlrtRTEInfo ic_emlrtRTEI = { 28,/* lineNo */
  1,                                   /* colNo */
  "gaussgradient",                     /* fName */
  "/media/data/Dropbox/code_bank/matlab/3Dnucleisegmentation/gaussgradient.m"/* pName */
};

static emlrtRTEInfo jc_emlrtRTEI = { 25,/* lineNo */
  4,                                   /* colNo */
  "gaussgradient",                     /* fName */
  "/media/data/Dropbox/code_bank/matlab/3Dnucleisegmentation/gaussgradient.m"/* pName */
};

static emlrtRTEInfo kc_emlrtRTEI = { 16,/* lineNo */
  1,                                   /* colNo */
  "gaussgradient",                     /* fName */
  "/media/data/Dropbox/code_bank/matlab/3Dnucleisegmentation/gaussgradient.m"/* pName */
};

static emlrtRTEInfo lc_emlrtRTEI = { 1,/* lineNo */
  18,                                  /* colNo */
  "gaussgradient",                     /* fName */
  "/media/data/Dropbox/code_bank/matlab/3Dnucleisegmentation/gaussgradient.m"/* pName */
};

static emlrtRTEInfo mc_emlrtRTEI = { 23,/* lineNo */
  16,                                  /* colNo */
  "gaussgradient",                     /* fName */
  "/media/data/Dropbox/code_bank/matlab/3Dnucleisegmentation/gaussgradient.m"/* pName */
};

static emlrtRTEInfo nc_emlrtRTEI = { 23,/* lineNo */
  20,                                  /* colNo */
  "gaussgradient",                     /* fName */
  "/media/data/Dropbox/code_bank/matlab/3Dnucleisegmentation/gaussgradient.m"/* pName */
};

/* Function Declarations */
static real_T gauss(real_T x, real_T sigma);

/* Function Definitions */
static real_T gauss(real_T x, real_T sigma)
{
  /* Gaussian */
  return muDoubleScalarExp(-(x * x) / (2.0 * (sigma * sigma))) / (sigma *
    2.5066282746310002);
}

void gaussgradient(const emlrtStack *sp, const emxArray_real_T *IM, real_T sigma,
                   emxArray_real_T *gx, emxArray_real_T *gy)
{
  real_T x;
  real_T halfsize;
  real_T size;
  int32_T i;
  emxArray_real_T *hx;
  int32_T vlen;
  int32_T b_i;
  emxArray_real_T *r;
  emxArray_real_T *y;
  emxArray_real_T *b_x;
  real_T u_idx_0;
  int32_T j;
  boolean_T p;
  int32_T xpageoffset;
  int32_T npages;
  boolean_T overflow;
  emlrtStack st;
  emlrtStack b_st;
  emlrtStack c_st;
  emlrtStack d_st;
  emlrtStack e_st;
  emlrtStack f_st;
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
  emlrtHeapReferenceStackEnterFcnR2012b(sp);

  /* GAUSSGRADIENT Gradient using first order derivative of Gaussian. */
  /*   [gx,gy]=gaussgradient(IM,sigma) outputs the gradient image gx and gy of */
  /*   image IM using a 2-D Gaussian kernel. Sigma is the standard deviation of */
  /*   this kernel along both directions. */
  /*  */
  /*   Contributed by Guanglei Xiong (xgl99@mails.tsinghua.edu.cn) */
  /*   at Tsinghua University, Beijing, China. */
  /* determine the appropriate size of kernel. The smaller epsilon, the larger */
  /* size. */
  st.site = &se_emlrtRSI;
  st.site = &se_emlrtRSI;
  x = 2.5066282746310002 * sigma * 0.01;
  if (x < 0.0) {
    emlrtErrorWithMessageIdR2018a(&st, &p_emlrtRTEI,
      "Coder:toolbox:ElFunDomainError", "Coder:toolbox:ElFunDomainError", 3, 4,
      3, "log");
  }

  x = muDoubleScalarLog(x);
  st.site = &se_emlrtRSI;
  x *= -2.0;
  if (x < 0.0) {
    emlrtErrorWithMessageIdR2018a(&st, &q_emlrtRTEI,
      "Coder:toolbox:ElFunDomainError", "Coder:toolbox:ElFunDomainError", 3, 4,
      4, "sqrt");
  }

  x = muDoubleScalarSqrt(x);
  halfsize = muDoubleScalarCeil(sigma * x);
  size = 2.0 * halfsize + 1.0;

  /* generate a 2-D Gaussian kernel along x direction */
  if (!(size >= 0.0)) {
    emlrtNonNegativeCheckR2012b(size, &f_emlrtDCI, sp);
  }

  i = static_cast<int32_T>(size);
  if (size != i) {
    emlrtIntegerCheckR2012b(size, &e_emlrtDCI, sp);
  }

  emxInit_real_T(sp, &hx, 2, &kc_emlrtRTEI, true);
  vlen = hx->size[0] * hx->size[1];
  hx->size[0] = i;
  emxEnsureCapacity_real_T(sp, hx, vlen, &ec_emlrtRTEI);
  if (size != static_cast<int32_T>(size)) {
    emlrtIntegerCheckR2012b(size, &g_emlrtDCI, sp);
  }

  vlen = hx->size[0] * hx->size[1];
  hx->size[1] = i;
  emxEnsureCapacity_real_T(sp, hx, vlen, &ec_emlrtRTEI);
  emlrtForLoopVectorCheckR2012b(1.0, 1.0, size, mxDOUBLE_CLASS,
    static_cast<int32_T>(size), &n_emlrtRTEI, sp);
  for (b_i = 0; b_i < i; b_i++) {
    emlrtForLoopVectorCheckR2012b(1.0, 1.0, size, mxDOUBLE_CLASS,
      static_cast<int32_T>(size), &o_emlrtRTEI, sp);
    if (0 <= i - 1) {
      u_idx_0 = ((static_cast<real_T>(b_i) + 1.0) - halfsize) - 1.0;
    }

    for (j = 0; j < i; j++) {
      x = ((static_cast<real_T>(j) + 1.0) - halfsize) - 1.0;
      st.site = &te_emlrtRSI;

      /* first order derivative of Gaussian */
      b_st.site = &ye_emlrtRSI;
      c_st.site = &ec_emlrtRSI;
      d_st.site = &fc_emlrtRSI;
      b_st.site = &ye_emlrtRSI;
      x = -x * gauss(x, sigma) / (sigma * sigma);
      vlen = static_cast<int32_T>((b_i + 1U));
      if ((vlen < 1) || (vlen > hx->size[0])) {
        emlrtDynamicBoundsCheckR2012b(vlen, 1, hx->size[0], &j_emlrtBCI, sp);
      }

      xpageoffset = static_cast<int32_T>((j + 1U));
      if ((xpageoffset < 1) || (xpageoffset > hx->size[1])) {
        emlrtDynamicBoundsCheckR2012b(xpageoffset, 1, hx->size[1], &j_emlrtBCI,
          sp);
      }

      st.site = &te_emlrtRSI;
      hx->data[(vlen + hx->size[0] * (xpageoffset - 1)) - 1] = gauss(u_idx_0,
        sigma) * x;
      if (*emlrtBreakCheckR2012bFlagVar != 0) {
        emlrtBreakCheckR2012b(sp);
      }
    }

    if (*emlrtBreakCheckR2012bFlagVar != 0) {
      emlrtBreakCheckR2012b(sp);
    }
  }

  emxInit_real_T(sp, &r, 2, &lc_emlrtRTEI, true);
  emxInit_real_T(sp, &y, 2, &mc_emlrtRTEI, true);
  emxInit_real_T(sp, &b_x, 2, &nc_emlrtRTEI, true);
  st.site = &ue_emlrtRSI;
  c_abs(&st, hx, b_x);
  st.site = &ue_emlrtRSI;
  c_abs(&st, hx, r);
  emlrtSizeEqCheckNDR2012b(*(int32_T (*)[2])b_x->size, *(int32_T (*)[2])r->size,
    &h_emlrtECI, sp);
  st.site = &ue_emlrtRSI;
  j = b_x->size[0] * b_x->size[1];
  for (i = 0; i < j; i++) {
    b_x->data[i] *= r->data[i];
  }

  emxFree_real_T(&r);
  b_st.site = &af_emlrtRSI;
  if ((b_x->size[0] == 1) && (b_x->size[1] != 1)) {
    emlrtErrorWithMessageIdR2018a(&b_st, &r_emlrtRTEI,
      "Coder:toolbox:autoDimIncompatibility",
      "Coder:toolbox:autoDimIncompatibility", 0);
  }

  p = false;
  if ((b_x->size[0] == 0) && (b_x->size[1] == 0)) {
    p = true;
  }

  if (p) {
    emlrtErrorWithMessageIdR2018a(&b_st, &s_emlrtRTEI,
      "Coder:toolbox:UnsupportedSpecialEmpty",
      "Coder:toolbox:UnsupportedSpecialEmpty", 0);
  }

  c_st.site = &bb_emlrtRSI;
  vlen = b_x->size[0];
  if ((b_x->size[0] == 0) || (b_x->size[1] == 0)) {
    i = y->size[0] * y->size[1];
    y->size[0] = 1;
    j = b_x->size[1];
    y->size[1] = j;
    emxEnsureCapacity_real_T(&c_st, y, i, &fc_emlrtRTEI);
    for (i = 0; i < j; i++) {
      y->data[i] = 0.0;
    }
  } else {
    d_st.site = &cb_emlrtRSI;
    npages = b_x->size[1];
    i = y->size[0] * y->size[1];
    y->size[0] = 1;
    y->size[1] = b_x->size[1];
    emxEnsureCapacity_real_T(&d_st, y, i, &gc_emlrtRTEI);
    e_st.site = &bf_emlrtRSI;
    if (b_x->size[1] > 2147483646) {
      f_st.site = &eb_emlrtRSI;
      check_forloop_overflow_error(&f_st);
    }

    if (0 <= b_x->size[1] - 1) {
      overflow = ((2 <= b_x->size[0]) && (b_x->size[0] > 2147483646));
    }

    for (b_i = 0; b_i < npages; b_i++) {
      xpageoffset = b_i * b_x->size[0];
      y->data[b_i] = b_x->data[xpageoffset];
      e_st.site = &db_emlrtRSI;
      if (overflow) {
        f_st.site = &eb_emlrtRSI;
        check_forloop_overflow_error(&f_st);
      }

      for (j = 2; j <= vlen; j++) {
        y->data[b_i] += b_x->data[(xpageoffset + j) - 1];
      }
    }
  }

  st.site = &ue_emlrtRSI;
  b_st.site = &af_emlrtRSI;
  c_st.site = &bb_emlrtRSI;
  vlen = y->size[1];
  if (y->size[1] == 0) {
    x = 0.0;
  } else {
    d_st.site = &cb_emlrtRSI;
    x = y->data[0];
    e_st.site = &db_emlrtRSI;
    overflow = ((2 <= y->size[1]) && (y->size[1] > 2147483646));
    if (overflow) {
      f_st.site = &eb_emlrtRSI;
      check_forloop_overflow_error(&f_st);
    }

    for (j = 2; j <= vlen; j++) {
      x += y->data[j - 1];
    }
  }

  emxFree_real_T(&y);
  st.site = &ue_emlrtRSI;
  if (x < 0.0) {
    emlrtErrorWithMessageIdR2018a(&st, &q_emlrtRTEI,
      "Coder:toolbox:ElFunDomainError", "Coder:toolbox:ElFunDomainError", 3, 4,
      4, "sqrt");
  }

  x = muDoubleScalarSqrt(x);
  j = hx->size[0] * hx->size[1];
  for (i = 0; i < j; i++) {
    hx->data[i] /= x;
  }

  /* generate a 2-D Gaussian kernel along y direction */
  /* 2-D filtering */
  i = gx->size[0] * gx->size[1];
  gx->size[0] = IM->size[0];
  gx->size[1] = IM->size[1];
  emxEnsureCapacity_real_T(sp, gx, i, &hc_emlrtRTEI);
  j = IM->size[0] * IM->size[1];
  for (i = 0; i < j; i++) {
    gx->data[i] = IM->data[i];
  }

  st.site = &ve_emlrtRSI;
  h_imfilter(&st, gx, hx);
  i = gy->size[0] * gy->size[1];
  gy->size[0] = IM->size[0];
  gy->size[1] = IM->size[1];
  emxEnsureCapacity_real_T(sp, gy, i, &ic_emlrtRTEI);
  j = IM->size[0] * IM->size[1];
  for (i = 0; i < j; i++) {
    gy->data[i] = IM->data[i];
  }

  i = b_x->size[0] * b_x->size[1];
  b_x->size[0] = hx->size[1];
  b_x->size[1] = hx->size[0];
  emxEnsureCapacity_real_T(sp, b_x, i, &jc_emlrtRTEI);
  j = hx->size[0];
  for (i = 0; i < j; i++) {
    xpageoffset = hx->size[1];
    for (vlen = 0; vlen < xpageoffset; vlen++) {
      b_x->data[vlen + b_x->size[0] * i] = hx->data[i + hx->size[0] * vlen];
    }
  }

  emxFree_real_T(&hx);
  st.site = &we_emlrtRSI;
  h_imfilter(&st, gy, b_x);
  emxFree_real_T(&b_x);
  emlrtHeapReferenceStackLeaveFcnR2012b(sp);
}

/* End of code generation (gaussgradient.cpp) */
