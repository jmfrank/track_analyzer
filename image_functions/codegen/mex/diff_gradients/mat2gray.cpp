/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * mat2gray.cpp
 *
 * Code generation for function 'mat2gray'
 *
 */

/* Include files */
#include "mat2gray.h"
#include "diff_gradients.h"
#include "diff_gradients_data.h"
#include "diff_gradients_emxutil.h"
#include "eml_int_forloop_overflow_check.h"
#include "imlincomb.h"
#include "mwmathutil.h"
#include "rt_nonfinite.h"

/* Variable Definitions */
static emlrtRSInfo je_emlrtRSI = { 33, /* lineNo */
  "mat2gray",                          /* fcnName */
  "/usr/local/MATLAB/R2019b/toolbox/images/images/mat2gray.m"/* pathName */
};

static emlrtRSInfo ke_emlrtRSI = { 42, /* lineNo */
  "mat2gray",                          /* fcnName */
  "/usr/local/MATLAB/R2019b/toolbox/images/images/mat2gray.m"/* pathName */
};

static emlrtRSInfo le_emlrtRSI = { 46, /* lineNo */
  "mat2gray",                          /* fcnName */
  "/usr/local/MATLAB/R2019b/toolbox/images/images/mat2gray.m"/* pathName */
};

static emlrtRSInfo me_emlrtRSI = { 14, /* lineNo */
  "min",                               /* fcnName */
  "/usr/local/MATLAB/R2019b/toolbox/eml/lib/matlab/datafun/min.m"/* pathName */
};

static emlrtRSInfo ne_emlrtRSI = { 155,/* lineNo */
  "unaryMinOrMax",                     /* fcnName */
  "/usr/local/MATLAB/R2019b/toolbox/eml/eml/+coder/+internal/unaryMinOrMax.m"/* pathName */
};

static emlrtRSInfo oe_emlrtRSI = { 1015,/* lineNo */
  "minRealVectorOmitNaN",              /* fcnName */
  "/usr/local/MATLAB/R2019b/toolbox/eml/eml/+coder/+internal/unaryMinOrMax.m"/* pathName */
};

static emlrtRSInfo pe_emlrtRSI = { 38, /* lineNo */
  "unaryOrBinaryDispatch",             /* fcnName */
  "/usr/local/MATLAB/R2019b/toolbox/eml/eml/+coder/+internal/minOrMax.m"/* pathName */
};

static emlrtRSInfo qe_emlrtRSI = { 62, /* lineNo */
  "binaryMinOrMax",                    /* fcnName */
  "/usr/local/MATLAB/R2019b/toolbox/eml/eml/+coder/+internal/binaryMinOrMax.m"/* pathName */
};

static emlrtRSInfo re_emlrtRSI = { 174,/* lineNo */
  "flatIter",                          /* fcnName */
  "/usr/local/MATLAB/R2019b/toolbox/eml/eml/+coder/+internal/applyBinaryScalarFunction.m"/* pathName */
};

static emlrtRTEInfo ac_emlrtRTEI = { 39,/* lineNo */
  4,                                   /* colNo */
  "mat2gray",                          /* fName */
  "/usr/local/MATLAB/R2019b/toolbox/images/images/mat2gray.m"/* pName */
};

static emlrtRTEInfo bc_emlrtRTEI = { 62,/* lineNo */
  10,                                  /* colNo */
  "binaryMinOrMax",                    /* fName */
  "/usr/local/MATLAB/R2019b/toolbox/eml/eml/+coder/+internal/binaryMinOrMax.m"/* pName */
};

static emlrtRTEInfo cc_emlrtRTEI = { 46,/* lineNo */
  11,                                  /* colNo */
  "mat2gray",                          /* fName */
  "/usr/local/MATLAB/R2019b/toolbox/images/images/mat2gray.m"/* pName */
};

/* Function Definitions */
void mat2gray(const emlrtStack *sp, const emxArray_real_T *A, emxArray_real_T
              *b_I)
{
  int32_T n;
  real_T ex;
  int32_T idx;
  boolean_T overflow;
  int32_T k;
  int32_T a;
  boolean_T exitg1;
  real_T delta;
  emxArray_real_T *minval;
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
  st.site = &je_emlrtRSI;
  b_st.site = &me_emlrtRSI;
  c_st.site = &jb_emlrtRSI;
  d_st.site = &kb_emlrtRSI;
  if (A->size[0] * A->size[1] < 1) {
    emlrtErrorWithMessageIdR2018a(&d_st, &m_emlrtRTEI,
      "Coder:toolbox:eml_min_or_max_varDimZero",
      "Coder:toolbox:eml_min_or_max_varDimZero", 0);
  }

  e_st.site = &ne_emlrtRSI;
  f_st.site = &oe_emlrtRSI;
  n = A->size[0] * A->size[1];
  if (A->size[0] * A->size[1] <= 2) {
    if (A->size[0] * A->size[1] == 1) {
      ex = A->data[0];
    } else if ((A->data[0] > A->data[1]) || (muDoubleScalarIsNaN(A->data[0]) &&
                (!muDoubleScalarIsNaN(A->data[1])))) {
      ex = A->data[1];
    } else {
      ex = A->data[0];
    }
  } else {
    g_st.site = &ob_emlrtRSI;
    if (!muDoubleScalarIsNaN(A->data[0])) {
      idx = 1;
    } else {
      idx = 0;
      h_st.site = &pb_emlrtRSI;
      overflow = ((2 <= A->size[0] * A->size[1]) && (A->size[0] * A->size[1] >
        2147483646));
      if (overflow) {
        i_st.site = &eb_emlrtRSI;
        check_forloop_overflow_error(&i_st);
      }

      k = 2;
      exitg1 = false;
      while ((!exitg1) && (k <= A->size[0] * A->size[1])) {
        if (!muDoubleScalarIsNaN(A->data[k - 1])) {
          idx = k;
          exitg1 = true;
        } else {
          k++;
        }
      }
    }

    if (idx == 0) {
      ex = A->data[0];
    } else {
      g_st.site = &nb_emlrtRSI;
      ex = A->data[idx - 1];
      a = idx + 1;
      h_st.site = &qb_emlrtRSI;
      overflow = ((idx + 1 <= A->size[0] * A->size[1]) && (A->size[0] * A->size
        [1] > 2147483646));
      if (overflow) {
        i_st.site = &eb_emlrtRSI;
        check_forloop_overflow_error(&i_st);
      }

      for (k = a; k <= n; k++) {
        if (ex > A->data[k - 1]) {
          ex = A->data[k - 1];
        }
      }
    }
  }

  st.site = &je_emlrtRSI;
  b_st.site = &ib_emlrtRSI;
  c_st.site = &jb_emlrtRSI;
  d_st.site = &kb_emlrtRSI;
  if (A->size[0] * A->size[1] < 1) {
    emlrtErrorWithMessageIdR2018a(&d_st, &m_emlrtRTEI,
      "Coder:toolbox:eml_min_or_max_varDimZero",
      "Coder:toolbox:eml_min_or_max_varDimZero", 0);
  }

  e_st.site = &lb_emlrtRSI;
  f_st.site = &mb_emlrtRSI;
  n = A->size[0] * A->size[1];
  if (A->size[0] * A->size[1] <= 2) {
    if (A->size[0] * A->size[1] == 1) {
      delta = A->data[0];
    } else if ((A->data[0] < A->data[1]) || (muDoubleScalarIsNaN(A->data[0]) &&
                (!muDoubleScalarIsNaN(A->data[1])))) {
      delta = A->data[1];
    } else {
      delta = A->data[0];
    }
  } else {
    g_st.site = &ob_emlrtRSI;
    if (!muDoubleScalarIsNaN(A->data[0])) {
      idx = 1;
    } else {
      idx = 0;
      h_st.site = &pb_emlrtRSI;
      overflow = ((2 <= A->size[0] * A->size[1]) && (A->size[0] * A->size[1] >
        2147483646));
      if (overflow) {
        i_st.site = &eb_emlrtRSI;
        check_forloop_overflow_error(&i_st);
      }

      k = 2;
      exitg1 = false;
      while ((!exitg1) && (k <= A->size[0] * A->size[1])) {
        if (!muDoubleScalarIsNaN(A->data[k - 1])) {
          idx = k;
          exitg1 = true;
        } else {
          k++;
        }
      }
    }

    if (idx == 0) {
      delta = A->data[0];
    } else {
      g_st.site = &nb_emlrtRSI;
      delta = A->data[idx - 1];
      a = idx + 1;
      h_st.site = &qb_emlrtRSI;
      overflow = ((idx + 1 <= A->size[0] * A->size[1]) && (A->size[0] * A->size
        [1] > 2147483646));
      if (overflow) {
        i_st.site = &eb_emlrtRSI;
        check_forloop_overflow_error(&i_st);
      }

      for (k = a; k <= n; k++) {
        if (delta < A->data[k - 1]) {
          delta = A->data[k - 1];
        }
      }
    }
  }

  if (delta == ex) {
    a = b_I->size[0] * b_I->size[1];
    b_I->size[0] = A->size[0];
    b_I->size[1] = A->size[1];
    emxEnsureCapacity_real_T(sp, b_I, a, &ac_emlrtRTEI);
    n = A->size[0] * A->size[1];
    for (a = 0; a < n; a++) {
      b_I->data[a] = A->data[a];
    }
  } else {
    delta = 1.0 / (delta - ex);
    st.site = &ke_emlrtRSI;
    b_imlincomb(&st, delta, A, -ex * delta, b_I);
  }

  emxInit_real_T(sp, &minval, 2, &cc_emlrtRTEI, true);
  st.site = &le_emlrtRSI;
  b_st.site = &me_emlrtRSI;
  c_st.site = &jb_emlrtRSI;
  d_st.site = &pe_emlrtRSI;
  e_st.site = &qe_emlrtRSI;
  a = b_I->size[0];
  n = minval->size[0] * minval->size[1];
  minval->size[0] = a;
  idx = b_I->size[1];
  minval->size[1] = idx;
  emxEnsureCapacity_real_T(&e_st, minval, n, &bc_emlrtRTEI);
  f_st.site = &fe_emlrtRSI;
  n = a * idx;
  g_st.site = &ge_emlrtRSI;
  overflow = ((1 <= n) && (n > 2147483646));
  if (overflow) {
    h_st.site = &eb_emlrtRSI;
    check_forloop_overflow_error(&h_st);
  }

  for (k = 0; k < n; k++) {
    minval->data[k] = muDoubleScalarMin(b_I->data[k], 1.0);
  }

  st.site = &le_emlrtRSI;
  b_st.site = &ib_emlrtRSI;
  c_st.site = &jb_emlrtRSI;
  d_st.site = &pe_emlrtRSI;
  e_st.site = &qe_emlrtRSI;
  a = minval->size[0];
  n = b_I->size[0] * b_I->size[1];
  b_I->size[0] = a;
  idx = minval->size[1];
  b_I->size[1] = idx;
  emxEnsureCapacity_real_T(&e_st, b_I, n, &bc_emlrtRTEI);
  f_st.site = &fe_emlrtRSI;
  n = a * idx;
  g_st.site = &re_emlrtRSI;
  overflow = ((1 <= n) && (n > 2147483646));
  if (overflow) {
    h_st.site = &eb_emlrtRSI;
    check_forloop_overflow_error(&h_st);
  }

  for (k = 0; k < n; k++) {
    b_I->data[k] = muDoubleScalarMax(0.0, minval->data[k]);
  }

  emxFree_real_T(&minval);
  emlrtHeapReferenceStackLeaveFcnR2012b(sp);
}

/* End of code generation (mat2gray.cpp) */
