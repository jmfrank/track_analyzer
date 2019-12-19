/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * power.cpp
 *
 * Code generation for function 'power'
 *
 */

/* Include files */
#include "power.h"
#include "diff_gradients.h"
#include "diff_gradients_data.h"
#include "diff_gradients_emxutil.h"
#include "eml_int_forloop_overflow_check.h"
#include "rt_nonfinite.h"

/* Variable Definitions */
static emlrtRSInfo ee_emlrtRSI = { 64, /* lineNo */
  "fltpower",                          /* fcnName */
  "/usr/local/MATLAB/R2019b/toolbox/eml/lib/matlab/ops/power.m"/* pathName */
};

static emlrtRTEInfo yb_emlrtRTEI = { 64,/* lineNo */
  5,                                   /* colNo */
  "power",                             /* fName */
  "/usr/local/MATLAB/R2019b/toolbox/eml/lib/matlab/ops/power.m"/* pName */
};

/* Function Definitions */
void power(const emlrtStack *sp, const emxArray_real_T *a, emxArray_real_T *y)
{
  uint32_T unnamed_idx_0;
  uint32_T unnamed_idx_1;
  int32_T nx;
  boolean_T overflow;
  int32_T k;
  emlrtStack st;
  emlrtStack b_st;
  emlrtStack c_st;
  emlrtStack d_st;
  emlrtStack e_st;
  st.prev = sp;
  st.tls = sp->tls;
  st.site = &fc_emlrtRSI;
  b_st.prev = &st;
  b_st.tls = st.tls;
  c_st.prev = &b_st;
  c_st.tls = b_st.tls;
  d_st.prev = &c_st;
  d_st.tls = c_st.tls;
  e_st.prev = &d_st;
  e_st.tls = d_st.tls;
  b_st.site = &ee_emlrtRSI;
  unnamed_idx_0 = static_cast<uint32_T>(a->size[0]);
  unnamed_idx_1 = static_cast<uint32_T>(a->size[1]);
  nx = y->size[0] * y->size[1];
  y->size[0] = static_cast<int32_T>(unnamed_idx_0);
  y->size[1] = static_cast<int32_T>(unnamed_idx_1);
  emxEnsureCapacity_real_T(&b_st, y, nx, &yb_emlrtRTEI);
  c_st.site = &fe_emlrtRSI;
  nx = static_cast<int32_T>(unnamed_idx_0) * static_cast<int32_T>(unnamed_idx_1);
  d_st.site = &ge_emlrtRSI;
  overflow = ((1 <= nx) && (nx > 2147483646));
  if (overflow) {
    e_st.site = &eb_emlrtRSI;
    check_forloop_overflow_error(&e_st);
  }

  for (k = 0; k < nx; k++) {
    y->data[k] = a->data[k] * a->data[k];
  }
}

/* End of code generation (power.cpp) */
