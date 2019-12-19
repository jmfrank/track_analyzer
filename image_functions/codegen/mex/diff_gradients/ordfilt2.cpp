/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * ordfilt2.cpp
 *
 * Code generation for function 'ordfilt2'
 *
 */

/* Include files */
#include "ordfilt2.h"
#include "abs.h"
#include "combineVectorElements.h"
#include "diff_gradients.h"
#include "diff_gradients_data.h"
#include "diff_gradients_emxutil.h"
#include "eml_int_forloop_overflow_check.h"
#include "libmwordfilt2.h"
#include "mwmathutil.h"
#include "padarray.h"
#include "rt_nonfinite.h"
#include "validateattributes.h"

/* Variable Definitions */
static emlrtRSInfo r_emlrtRSI = { 155, /* lineNo */
  "ordfilt2",                          /* fcnName */
  "/usr/local/MATLAB/R2019b/toolbox/images/images/eml/ordfilt2.m"/* pathName */
};

static emlrtRSInfo s_emlrtRSI = { 137, /* lineNo */
  "ordfilt2",                          /* fcnName */
  "/usr/local/MATLAB/R2019b/toolbox/images/images/eml/ordfilt2.m"/* pathName */
};

static emlrtRSInfo t_emlrtRSI = { 134, /* lineNo */
  "ordfilt2",                          /* fcnName */
  "/usr/local/MATLAB/R2019b/toolbox/images/images/eml/ordfilt2.m"/* pathName */
};

static emlrtRSInfo u_emlrtRSI = { 107, /* lineNo */
  "ordfilt2",                          /* fcnName */
  "/usr/local/MATLAB/R2019b/toolbox/images/images/eml/ordfilt2.m"/* pathName */
};

static emlrtRSInfo v_emlrtRSI = { 101, /* lineNo */
  "ordfilt2",                          /* fcnName */
  "/usr/local/MATLAB/R2019b/toolbox/images/images/eml/ordfilt2.m"/* pathName */
};

static emlrtRSInfo w_emlrtRSI = { 28,  /* lineNo */
  "ordfilt2",                          /* fcnName */
  "/usr/local/MATLAB/R2019b/toolbox/images/images/eml/ordfilt2.m"/* pathName */
};

static emlrtRSInfo x_emlrtRSI = { 25,  /* lineNo */
  "ordfilt2",                          /* fcnName */
  "/usr/local/MATLAB/R2019b/toolbox/images/images/eml/ordfilt2.m"/* pathName */
};

static emlrtRSInfo y_emlrtRSI = { 19,  /* lineNo */
  "ordfilt2",                          /* fcnName */
  "/usr/local/MATLAB/R2019b/toolbox/images/images/eml/ordfilt2.m"/* pathName */
};

static emlrtRSInfo fb_emlrtRSI = { 19, /* lineNo */
  "ind2sub",                           /* fcnName */
  "/usr/local/MATLAB/R2019b/toolbox/eml/lib/matlab/elmat/ind2sub.m"/* pathName */
};

static emlrtRSInfo xb_emlrtRSI = { 179,/* lineNo */
  "ordfilt2SharedLibrary",             /* fcnName */
  "/usr/local/MATLAB/R2019b/toolbox/images/images/eml/ordfilt2.m"/* pathName */
};

static emlrtECInfo j_emlrtECI = { -1,  /* nDims */
  153,                                 /* lineNo */
  25,                                  /* colNo */
  "ordfilt2",                          /* fName */
  "/usr/local/MATLAB/R2019b/toolbox/images/images/eml/ordfilt2.m"/* pName */
};

static emlrtBCInfo o_emlrtBCI = { -1,  /* iFirst */
  -1,                                  /* iLast */
  106,                                 /* lineNo */
  23,                                  /* colNo */
  "",                                  /* aName */
  "ordfilt2",                          /* fName */
  "/usr/local/MATLAB/R2019b/toolbox/images/images/eml/ordfilt2.m",/* pName */
  0                                    /* checkKind */
};

static emlrtRTEInfo ab_emlrtRTEI = { 28,/* lineNo */
  1,                                   /* colNo */
  "ordfilt2",                          /* fName */
  "/usr/local/MATLAB/R2019b/toolbox/images/images/eml/ordfilt2.m"/* pName */
};

static emlrtRTEInfo bb_emlrtRTEI = { 38,/* lineNo */
  15,                                  /* colNo */
  "ind2sub_indexClass",                /* fName */
  "/usr/local/MATLAB/R2019b/toolbox/eml/lib/matlab/elmat/ind2sub.m"/* pName */
};

static emlrtDCInfo m_emlrtDCI = { 102, /* lineNo */
  22,                                  /* colNo */
  "ordfilt2",                          /* fName */
  "/usr/local/MATLAB/R2019b/toolbox/images/images/eml/ordfilt2.m",/* pName */
  4                                    /* checkKind */
};

static emlrtBCInfo p_emlrtBCI = { -1,  /* iFirst */
  -1,                                  /* iLast */
  107,                                 /* lineNo */
  18,                                  /* colNo */
  "",                                  /* aName */
  "ordfilt2",                          /* fName */
  "/usr/local/MATLAB/R2019b/toolbox/images/images/eml/ordfilt2.m",/* pName */
  0                                    /* checkKind */
};

static emlrtBCInfo q_emlrtBCI = { -1,  /* iFirst */
  -1,                                  /* iLast */
  107,                                 /* lineNo */
  31,                                  /* colNo */
  "",                                  /* aName */
  "ordfilt2",                          /* fName */
  "/usr/local/MATLAB/R2019b/toolbox/images/images/eml/ordfilt2.m",/* pName */
  0                                    /* checkKind */
};

static emlrtRTEInfo nd_emlrtRTEI = { 25,/* lineNo */
  1,                                   /* colNo */
  "ordfilt2",                          /* fName */
  "/usr/local/MATLAB/R2019b/toolbox/images/images/eml/ordfilt2.m"/* pName */
};

static emlrtRTEInfo od_emlrtRTEI = { 102,/* lineNo */
  9,                                   /* colNo */
  "ordfilt2",                          /* fName */
  "/usr/local/MATLAB/R2019b/toolbox/images/images/eml/ordfilt2.m"/* pName */
};

static emlrtRTEInfo pd_emlrtRTEI = { 103,/* lineNo */
  9,                                   /* colNo */
  "ordfilt2",                          /* fName */
  "/usr/local/MATLAB/R2019b/toolbox/images/images/eml/ordfilt2.m"/* pName */
};

static emlrtRTEInfo qd_emlrtRTEI = { 153,/* lineNo */
  9,                                   /* colNo */
  "ordfilt2",                          /* fName */
  "/usr/local/MATLAB/R2019b/toolbox/images/images/eml/ordfilt2.m"/* pName */
};

static emlrtRTEInfo rd_emlrtRTEI = { 155,/* lineNo */
  13,                                  /* colNo */
  "ordfilt2",                          /* fName */
  "/usr/local/MATLAB/R2019b/toolbox/images/images/eml/ordfilt2.m"/* pName */
};

static emlrtRTEInfo sd_emlrtRTEI = { 99,/* lineNo */
  10,                                  /* colNo */
  "ordfilt2",                          /* fName */
  "/usr/local/MATLAB/R2019b/toolbox/images/images/eml/ordfilt2.m"/* pName */
};

static emlrtRTEInfo td_emlrtRTEI = { 99,/* lineNo */
  16,                                  /* colNo */
  "ordfilt2",                          /* fName */
  "/usr/local/MATLAB/R2019b/toolbox/images/images/eml/ordfilt2.m"/* pName */
};

static emlrtRTEInfo ud_emlrtRTEI = { 137,/* lineNo */
  5,                                   /* colNo */
  "ordfilt2",                          /* fName */
  "/usr/local/MATLAB/R2019b/toolbox/images/images/eml/ordfilt2.m"/* pName */
};

static emlrtRTEInfo vd_emlrtRTEI = { 134,/* lineNo */
  19,                                  /* colNo */
  "ordfilt2",                          /* fName */
  "/usr/local/MATLAB/R2019b/toolbox/images/images/eml/ordfilt2.m"/* pName */
};

static emlrtRSInfo bh_emlrtRSI = { 18, /* lineNo */
  "indexDivide",                       /* fcnName */
  "/usr/local/MATLAB/R2019b/toolbox/eml/eml/+coder/+internal/indexDivide.m"/* pathName */
};

/* Function Declarations */
static int32_T div_s32(const emlrtStack *sp, int32_T numerator, int32_T
  denominator);

/* Function Definitions */
static int32_T div_s32(const emlrtStack *sp, int32_T numerator, int32_T
  denominator)
{
  int32_T quotient;
  uint32_T b_numerator;
  uint32_T b_denominator;
  if (denominator == 0) {
    emlrtDivisionByZeroErrorR2012b(NULL, sp);
  } else {
    if (numerator < 0) {
      b_numerator = ~static_cast<uint32_T>(numerator) + 1U;
    } else {
      b_numerator = static_cast<uint32_T>(numerator);
    }

    if (denominator < 0) {
      b_denominator = ~static_cast<uint32_T>(denominator) + 1U;
    } else {
      b_denominator = static_cast<uint32_T>(denominator);
    }

    b_numerator /= b_denominator;
    if ((numerator < 0) != (denominator < 0)) {
      quotient = -static_cast<int32_T>(b_numerator);
    } else {
      quotient = static_cast<int32_T>(b_numerator);
    }
  }

  return quotient;
}

void ordfilt2(const emlrtStack *sp, emxArray_real_T *varargin_1, real_T
              varargin_2, const emxArray_real_T *varargin_3)
{
  emxArray_boolean_T *domain;
  int32_T i;
  int32_T nz;
  boolean_T guard1 = false;
  real_T domainSize[2];
  emxArray_boolean_T b_domain;
  real_T center[2];
  int32_T c_domain[1];
  int32_T d_domain[1];
  emxArray_real_T *cols;
  emxArray_real_T *rows;
  uint32_T b_index;
  real_T startIdx[2];
  int32_T b_i;
  int32_T a;
  emxArray_real_T *b_varargin_1;
  int32_T idx;
  real_T padSize;
  boolean_T exitg1;
  boolean_T overflow;
  real_T d;
  real_T ex;
  emxArray_real_T *Apad;
  emxArray_int32_T *indices;
  emlrtStack st;
  emlrtStack b_st;
  emlrtStack c_st;
  emlrtStack d_st;
  emlrtStack e_st;
  emlrtStack f_st;
  emlrtStack g_st;
  emlrtStack h_st;
  emlrtStack i_st;
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
  emlrtHeapReferenceStackEnterFcnR2012b(sp);
  emxInit_boolean_T(sp, &domain, 2, &nd_emlrtRTEI, true);
  st.site = &y_emlrtRSI;
  validateattributes(&st, varargin_2);
  st.site = &x_emlrtRSI;
  i = domain->size[0] * domain->size[1];
  domain->size[0] = varargin_3->size[0];
  domain->size[1] = varargin_3->size[1];
  emxEnsureCapacity_boolean_T(sp, domain, i, &nd_emlrtRTEI);
  nz = varargin_3->size[0] * varargin_3->size[1];
  for (i = 0; i < nz; i++) {
    domain->data[i] = (varargin_3->data[i] != 0.0);
  }

  guard1 = false;
  if (varargin_2 < 1.0) {
    guard1 = true;
  } else {
    st.site = &w_emlrtRSI;
    b_st.site = &ab_emlrtRSI;
    nz = domain->size[0] * domain->size[1];
    b_domain = *domain;
    c_domain[0] = nz;
    b_domain.size = &c_domain[0];
    b_domain.numDimensions = 1;
    c_st.site = &bb_emlrtRSI;
    nz = combineVectorElements(&c_st, &b_domain);
    if (varargin_2 > nz) {
      guard1 = true;
    }
  }

  if (guard1) {
    emlrtErrorWithMessageIdR2018a(sp, &ab_emlrtRTEI,
      "images:ordfilt2:orderNotValid", "images:ordfilt2:orderNotValid", 0);
  }

  domainSize[0] = domain->size[0];
  center[0] = muDoubleScalarFloor((static_cast<real_T>(domain->size[0]) + 1.0) /
    2.0);
  domainSize[1] = domain->size[1];
  center[1] = muDoubleScalarFloor((static_cast<real_T>(domain->size[1]) + 1.0) /
    2.0);
  st.site = &v_emlrtRSI;
  b_st.site = &ab_emlrtRSI;
  nz = domain->size[0] * domain->size[1];
  b_domain = *domain;
  d_domain[0] = nz;
  b_domain.size = &d_domain[0];
  b_domain.numDimensions = 1;
  c_st.site = &bb_emlrtRSI;
  nz = combineVectorElements(&c_st, &b_domain);
  if (nz < 0) {
    emlrtNonNegativeCheckR2012b(static_cast<real_T>(nz), &m_emlrtDCI, sp);
  }

  emxInit_real_T(sp, &cols, 1, &td_emlrtRTEI, true);
  i = cols->size[0];
  cols->size[0] = nz;
  emxEnsureCapacity_real_T(sp, cols, i, &od_emlrtRTEI);
  for (i = 0; i < nz; i++) {
    cols->data[i] = 0.0;
  }

  emxInit_real_T(sp, &rows, 1, &sd_emlrtRTEI, true);
  i = rows->size[0];
  rows->size[0] = nz;
  emxEnsureCapacity_real_T(sp, rows, i, &pd_emlrtRTEI);
  for (i = 0; i < nz; i++) {
    rows->data[i] = 0.0;
  }

  b_index = 1U;
  i = domain->size[0] * domain->size[1];
  if (0 <= i - 1) {
    startIdx[0] = domain->size[0];
    startIdx[1] = domain->size[1];
  }

  for (b_i = 0; b_i < i; b_i++) {
    a = domain->size[0] * domain->size[1];
    nz = b_i + 1;
    if ((nz < 1) || (nz > a)) {
      emlrtDynamicBoundsCheckR2012b(nz, 1, a, &o_emlrtBCI, sp);
    }

    st.site = &u_emlrtRSI;
    b_st.site = &fb_emlrtRSI;
    a = static_cast<int32_T>(static_cast<uint32_T>(startIdx[0]));
    if (b_i + 1 > a * static_cast<int32_T>(static_cast<uint32_T>(startIdx[1])))
    {
      emlrtErrorWithMessageIdR2018a(&b_st, &bb_emlrtRTEI,
        "Coder:MATLAB:ind2sub_IndexOutOfRange",
        "Coder:MATLAB:ind2sub_IndexOutOfRange", 0);
    }

    c_st.site = &bh_emlrtRSI;
    nz = div_s32(&c_st, b_i, a);
    if ((static_cast<int32_T>(b_index) < 1) || (static_cast<int32_T>(b_index) >
         rows->size[0])) {
      emlrtDynamicBoundsCheckR2012b(static_cast<int32_T>(b_index), 1, rows->
        size[0], &p_emlrtBCI, &st);
    }

    idx = static_cast<int32_T>(b_index) - 1;
    rows->data[idx] = (b_i - nz * a) + 1;
    if ((static_cast<int32_T>(b_index) < 1) || (static_cast<int32_T>(b_index) >
         cols->size[0])) {
      emlrtDynamicBoundsCheckR2012b(static_cast<int32_T>(b_index), 1, cols->
        size[0], &q_emlrtBCI, &st);
    }

    cols->data[idx] = nz + 1;
    b_index++;
  }

  emxFree_boolean_T(&domain);
  nz = rows->size[0];
  for (i = 0; i < nz; i++) {
    rows->data[i] -= center[0];
  }

  nz = cols->size[0];
  for (i = 0; i < nz; i++) {
    cols->data[i] -= center[1];
  }

  emxInit_real_T(sp, &b_varargin_1, 1, &vd_emlrtRTEI, true);
  st.site = &t_emlrtRSI;
  b_st.site = &t_emlrtRSI;
  b_abs(&b_st, rows, b_varargin_1);
  b_st.site = &ib_emlrtRSI;
  c_st.site = &jb_emlrtRSI;
  d_st.site = &kb_emlrtRSI;
  if (b_varargin_1->size[0] < 1) {
    emlrtErrorWithMessageIdR2018a(&d_st, &m_emlrtRTEI,
      "Coder:toolbox:eml_min_or_max_varDimZero",
      "Coder:toolbox:eml_min_or_max_varDimZero", 0);
  }

  e_st.site = &lb_emlrtRSI;
  f_st.site = &mb_emlrtRSI;
  nz = b_varargin_1->size[0];
  if (b_varargin_1->size[0] <= 2) {
    if (b_varargin_1->size[0] == 1) {
      padSize = b_varargin_1->data[0];
    } else if ((b_varargin_1->data[0] < b_varargin_1->data[1]) ||
               (muDoubleScalarIsNaN(b_varargin_1->data[0]) &&
                (!muDoubleScalarIsNaN(b_varargin_1->data[1])))) {
      padSize = b_varargin_1->data[1];
    } else {
      padSize = b_varargin_1->data[0];
    }
  } else {
    g_st.site = &ob_emlrtRSI;
    if (!muDoubleScalarIsNaN(b_varargin_1->data[0])) {
      idx = 1;
    } else {
      idx = 0;
      h_st.site = &pb_emlrtRSI;
      if (b_varargin_1->size[0] > 2147483646) {
        i_st.site = &eb_emlrtRSI;
        check_forloop_overflow_error(&i_st);
      }

      b_i = 2;
      exitg1 = false;
      while ((!exitg1) && (b_i <= b_varargin_1->size[0])) {
        if (!muDoubleScalarIsNaN(b_varargin_1->data[b_i - 1])) {
          idx = b_i;
          exitg1 = true;
        } else {
          b_i++;
        }
      }
    }

    if (idx == 0) {
      padSize = b_varargin_1->data[0];
    } else {
      g_st.site = &nb_emlrtRSI;
      padSize = b_varargin_1->data[idx - 1];
      a = idx + 1;
      h_st.site = &qb_emlrtRSI;
      overflow = ((idx + 1 <= b_varargin_1->size[0]) && (b_varargin_1->size[0] >
        2147483646));
      if (overflow) {
        i_st.site = &eb_emlrtRSI;
        check_forloop_overflow_error(&i_st);
      }

      for (b_i = a; b_i <= nz; b_i++) {
        d = b_varargin_1->data[b_i - 1];
        if (padSize < d) {
          padSize = d;
        }
      }
    }
  }

  st.site = &t_emlrtRSI;
  b_st.site = &t_emlrtRSI;
  b_abs(&b_st, cols, b_varargin_1);
  b_st.site = &ib_emlrtRSI;
  c_st.site = &jb_emlrtRSI;
  d_st.site = &kb_emlrtRSI;
  if (b_varargin_1->size[0] < 1) {
    emlrtErrorWithMessageIdR2018a(&d_st, &m_emlrtRTEI,
      "Coder:toolbox:eml_min_or_max_varDimZero",
      "Coder:toolbox:eml_min_or_max_varDimZero", 0);
  }

  e_st.site = &lb_emlrtRSI;
  f_st.site = &mb_emlrtRSI;
  nz = b_varargin_1->size[0];
  if (b_varargin_1->size[0] <= 2) {
    if (b_varargin_1->size[0] == 1) {
      ex = b_varargin_1->data[0];
    } else if ((b_varargin_1->data[0] < b_varargin_1->data[1]) ||
               (muDoubleScalarIsNaN(b_varargin_1->data[0]) &&
                (!muDoubleScalarIsNaN(b_varargin_1->data[1])))) {
      ex = b_varargin_1->data[1];
    } else {
      ex = b_varargin_1->data[0];
    }
  } else {
    g_st.site = &ob_emlrtRSI;
    if (!muDoubleScalarIsNaN(b_varargin_1->data[0])) {
      idx = 1;
    } else {
      idx = 0;
      h_st.site = &pb_emlrtRSI;
      if (b_varargin_1->size[0] > 2147483646) {
        i_st.site = &eb_emlrtRSI;
        check_forloop_overflow_error(&i_st);
      }

      b_i = 2;
      exitg1 = false;
      while ((!exitg1) && (b_i <= b_varargin_1->size[0])) {
        if (!muDoubleScalarIsNaN(b_varargin_1->data[b_i - 1])) {
          idx = b_i;
          exitg1 = true;
        } else {
          b_i++;
        }
      }
    }

    if (idx == 0) {
      ex = b_varargin_1->data[0];
    } else {
      g_st.site = &nb_emlrtRSI;
      ex = b_varargin_1->data[idx - 1];
      a = idx + 1;
      h_st.site = &qb_emlrtRSI;
      overflow = ((idx + 1 <= b_varargin_1->size[0]) && (b_varargin_1->size[0] >
        2147483646));
      if (overflow) {
        i_st.site = &eb_emlrtRSI;
        check_forloop_overflow_error(&i_st);
      }

      for (b_i = a; b_i <= nz; b_i++) {
        d = b_varargin_1->data[b_i - 1];
        if (ex < d) {
          ex = d;
        }
      }
    }
  }

  emxFree_real_T(&b_varargin_1);
  emxInit_real_T(&f_st, &Apad, 2, &ud_emlrtRTEI, true);
  padSize = muDoubleScalarMax(padSize, ex);
  center[0] = padSize;
  center[1] = padSize;
  st.site = &s_emlrtRSI;
  padarray(&st, varargin_1, center, Apad);
  if ((varargin_1->size[0] != 0) && (varargin_1->size[1] != 0)) {
    nz = cols->size[0];
    for (i = 0; i < nz; i++) {
      cols->data[i] *= static_cast<real_T>(Apad->size[0]);
    }

    if (cols->size[0] != rows->size[0]) {
      emlrtSizeEqCheck1DR2012b(cols->size[0], rows->size[0], &j_emlrtECI, sp);
    }

    emxInit_int32_T(sp, &indices, 1, &qd_emlrtRTEI, true);
    i = indices->size[0];
    indices->size[0] = cols->size[0];
    emxEnsureCapacity_int32_T(sp, indices, i, &qd_emlrtRTEI);
    nz = cols->size[0];
    for (i = 0; i < nz; i++) {
      d = muDoubleScalarRound(cols->data[i] + rows->data[i]);
      if (d < 2.147483648E+9) {
        if (d >= -2.147483648E+9) {
          a = static_cast<int32_T>(d);
        } else {
          a = MIN_int32_T;
        }
      } else if (d >= 2.147483648E+9) {
        a = MAX_int32_T;
      } else {
        a = 0;
      }

      indices->data[i] = a;
    }

    startIdx[0] = padSize;
    startIdx[1] = padSize;
    st.site = &r_emlrtRSI;
    nz = varargin_1->size[1];
    i = varargin_1->size[0] * varargin_1->size[1];
    varargin_1->size[1] = nz;
    emxEnsureCapacity_real_T(&st, varargin_1, i, &rd_emlrtRTEI);
    b_st.site = &xb_emlrtRSI;
    center[0] = varargin_1->size[0];
    center[1] = varargin_1->size[1];
    ordfilt2_real64(&Apad->data[0], static_cast<real_T>(Apad->size[0]), startIdx,
                    &indices->data[0], static_cast<real_T>(indices->size[0]),
                    domainSize, varargin_2 - 1.0, &varargin_1->data[0], center,
                    true);
    emxFree_int32_T(&indices);
  }

  emxFree_real_T(&Apad);
  emxFree_real_T(&cols);
  emxFree_real_T(&rows);
  emlrtHeapReferenceStackLeaveFcnR2012b(sp);
}

/* End of code generation (ordfilt2.cpp) */
