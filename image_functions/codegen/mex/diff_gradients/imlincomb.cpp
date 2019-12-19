/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * imlincomb.cpp
 *
 * Code generation for function 'imlincomb'
 *
 */

/* Include files */
#include "imlincomb.h"
#include "diff_gradients.h"
#include "diff_gradients_emxutil.h"
#include "libmwimlincomb_tbb.h"
#include "rt_nonfinite.h"

/* Variable Definitions */
static emlrtRSInfo yb_emlrtRSI = { 11, /* lineNo */
  "imlincomb",                         /* fcnName */
  "/usr/local/MATLAB/R2019b/toolbox/images/images/eml/imlincomb.m"/* pathName */
};

static emlrtRSInfo ac_emlrtRSI = { 27, /* lineNo */
  "imlincomb",                         /* fcnName */
  "/usr/local/MATLAB/R2019b/toolbox/images/images/eml/imlincomb.m"/* pathName */
};

static emlrtRSInfo bc_emlrtRSI = { 137,/* lineNo */
  "lincombGeneric",                    /* fcnName */
  "/usr/local/MATLAB/R2019b/toolbox/images/images/eml/imlincomb.m"/* pathName */
};

static emlrtRSInfo cc_emlrtRSI = { 139,/* lineNo */
  "lincombGeneric",                    /* fcnName */
  "/usr/local/MATLAB/R2019b/toolbox/images/images/eml/imlincomb.m"/* pathName */
};

static emlrtRTEInfo k_emlrtRTEI = { 95,/* lineNo */
  5,                                   /* colNo */
  "parseInputs",                       /* fName */
  "/usr/local/MATLAB/R2019b/toolbox/images/images/eml/imlincomb.m"/* pName */
};

static emlrtBCInfo g_emlrtBCI = { -1,  /* iFirst */
  -1,                                  /* iLast */
  218,                                 /* lineNo */
  55,                                  /* colNo */
  "",                                  /* aName */
  "lincombPortableCode",               /* fName */
  "/usr/local/MATLAB/R2019b/toolbox/images/images/eml/imlincomb.m",/* pName */
  0                                    /* checkKind */
};

static emlrtBCInfo h_emlrtBCI = { -1,  /* iFirst */
  -1,                                  /* iLast */
  218,                                 /* lineNo */
  39,                                  /* colNo */
  "",                                  /* aName */
  "lincombPortableCode",               /* fName */
  "/usr/local/MATLAB/R2019b/toolbox/images/images/eml/imlincomb.m",/* pName */
  0                                    /* checkKind */
};

static emlrtBCInfo i_emlrtBCI = { -1,  /* iFirst */
  -1,                                  /* iLast */
  223,                                 /* lineNo */
  9,                                   /* colNo */
  "",                                  /* aName */
  "lincombPortableCode",               /* fName */
  "/usr/local/MATLAB/R2019b/toolbox/images/images/eml/imlincomb.m",/* pName */
  0                                    /* checkKind */
};

static emlrtRTEInfo ub_emlrtRTEI = { 11,/* lineNo */
  1,                                   /* colNo */
  "imlincomb",                         /* fName */
  "/usr/local/MATLAB/R2019b/toolbox/images/images/eml/imlincomb.m"/* pName */
};

static emlrtRTEInfo vb_emlrtRTEI = { 139,/* lineNo */
  9,                                   /* colNo */
  "imlincomb",                         /* fName */
  "/usr/local/MATLAB/R2019b/toolbox/images/images/eml/imlincomb.m"/* pName */
};

static emlrtRTEInfo wb_emlrtRTEI = { 137,/* lineNo */
  9,                                   /* colNo */
  "imlincomb",                         /* fName */
  "/usr/local/MATLAB/R2019b/toolbox/images/images/eml/imlincomb.m"/* pName */
};

static emlrtRTEInfo xb_emlrtRTEI = { 1,/* lineNo */
  14,                                  /* colNo */
  "imlincomb",                         /* fName */
  "/usr/local/MATLAB/R2019b/toolbox/images/images/eml/imlincomb.m"/* pName */
};

static emlrtRTEInfo dc_emlrtRTEI = { 25,/* lineNo */
  9,                                   /* colNo */
  "imlincomb",                         /* fName */
  "/usr/local/MATLAB/R2019b/toolbox/images/images/eml/imlincomb.m"/* pName */
};

/* Function Definitions */
void b_imlincomb(const emlrtStack *sp, real_T varargin_1, const emxArray_real_T *
                 varargin_2, real_T varargin_3, emxArray_real_T *Z)
{
  int32_T i;
  real_T multipliers[2];
  int32_T k;
  int32_T i1;
  int32_T i2;
  emlrtStack st;
  emlrtStack b_st;
  st.prev = sp;
  st.tls = sp->tls;
  b_st.prev = &st;
  b_st.tls = st.tls;
  if (varargin_1 == 1.0) {
    i = Z->size[0] * Z->size[1];
    Z->size[0] = varargin_2->size[0];
    Z->size[1] = varargin_2->size[1];
    emxEnsureCapacity_real_T(sp, Z, i, &dc_emlrtRTEI);
    k = varargin_2->size[0] * varargin_2->size[1];
    for (i = 0; i < k; i++) {
      Z->data[i] = varargin_2->data[i] + varargin_3;
    }
  } else {
    st.site = &ac_emlrtRSI;
    if (varargin_2->size[0] * varargin_2->size[1] > 500000) {
      b_st.site = &bc_emlrtRSI;
      multipliers[0] = varargin_1;
      multipliers[1] = varargin_3;
      i = Z->size[0] * Z->size[1];
      Z->size[0] = varargin_2->size[0];
      Z->size[1] = varargin_2->size[1];
      emxEnsureCapacity_real_T(&b_st, Z, i, &wb_emlrtRTEI);
      imlincomb_tbb_real64(multipliers, 2.0, &Z->data[0], 0, static_cast<real_T>
                           ((varargin_2->size[0] * varargin_2->size[1])), 1.0,
                           &varargin_2->data[0]);
    } else {
      b_st.site = &cc_emlrtRSI;
      i = Z->size[0] * Z->size[1];
      Z->size[0] = varargin_2->size[0];
      Z->size[1] = varargin_2->size[1];
      emxEnsureCapacity_real_T(&b_st, Z, i, &vb_emlrtRTEI);
      i = varargin_2->size[0] * varargin_2->size[1];
      for (k = 0; k < i; k++) {
        i1 = varargin_2->size[0] * varargin_2->size[1];
        i2 = k + 1;
        if ((i2 < 1) || (i2 > i1)) {
          emlrtDynamicBoundsCheckR2012b(i2, 1, i1, &g_emlrtBCI, &b_st);
        }

        i1 = Z->size[0] * Z->size[1];
        i2 = static_cast<int32_T>((k + 1U));
        if ((i2 < 1) || (i2 > i1)) {
          emlrtDynamicBoundsCheckR2012b(i2, 1, i1, &i_emlrtBCI, &b_st);
        }

        Z->data[i2 - 1] = varargin_1 * varargin_2->data[k] + varargin_3;
      }
    }
  }
}

void imlincomb(const emlrtStack *sp, const emxArray_real_T *varargin_2, const
               emxArray_real_T *varargin_4, emxArray_real_T *Z)
{
  cell_wrap_14 r;
  int32_T i;
  int32_T loop_ub;
  cell_wrap_14 r1;
  real_T inputSize[2];
  uint32_T varargin_1[2];
  boolean_T p;
  boolean_T exitg1;
  boolean_T b_p;
  int32_T i1;
  int32_T i2;
  emlrtStack st;
  emlrtStack b_st;
  st.prev = sp;
  st.tls = sp->tls;
  b_st.prev = &st;
  b_st.tls = st.tls;
  emlrtHeapReferenceStackEnterFcnR2012b(sp);
  emxInitStruct_cell_wrap_14(sp, &r, &xb_emlrtRTEI, true);
  st.site = &yb_emlrtRSI;
  i = r.f1->size[0] * r.f1->size[1];
  r.f1->size[0] = varargin_2->size[0];
  r.f1->size[1] = varargin_2->size[1];
  emxEnsureCapacity_real_T(&st, r.f1, i, &ub_emlrtRTEI);
  loop_ub = varargin_2->size[0] * varargin_2->size[1];
  for (i = 0; i < loop_ub; i++) {
    r.f1->data[i] = varargin_2->data[i];
  }

  emxInitStruct_cell_wrap_14(&st, &r1, &xb_emlrtRTEI, true);
  i = r1.f1->size[0] * r1.f1->size[1];
  r1.f1->size[0] = varargin_4->size[0];
  r1.f1->size[1] = varargin_4->size[1];
  emxEnsureCapacity_real_T(&st, r1.f1, i, &ub_emlrtRTEI);
  loop_ub = varargin_4->size[0] * varargin_4->size[1];
  for (i = 0; i < loop_ub; i++) {
    r1.f1->data[i] = varargin_4->data[i];
  }

  inputSize[0] = r.f1->size[0];
  varargin_1[0] = static_cast<uint32_T>(r1.f1->size[0]);
  inputSize[1] = r.f1->size[1];
  varargin_1[1] = static_cast<uint32_T>(r1.f1->size[1]);
  p = true;
  loop_ub = 0;
  emxFreeStruct_cell_wrap_14(&r1);
  emxFreeStruct_cell_wrap_14(&r);
  exitg1 = false;
  while ((!exitg1) && (loop_ub < 2)) {
    if (static_cast<int32_T>(varargin_1[loop_ub]) != static_cast<int32_T>(
         static_cast<uint32_T>(inputSize[loop_ub]))) {
      p = false;
      exitg1 = true;
    } else {
      loop_ub++;
    }
  }

  b_p = p;
  if (!b_p) {
    emlrtErrorWithMessageIdR2018a(&st, &k_emlrtRTEI,
      "images:imlincomb:mismatchedArraySize",
      "images:imlincomb:mismatchedArraySize", 0);
  }

  st.site = &ac_emlrtRSI;
  if (varargin_2->size[0] * varargin_2->size[1] > 500000) {
    b_st.site = &bc_emlrtRSI;
    i = Z->size[0] * Z->size[1];
    Z->size[0] = varargin_2->size[0];
    Z->size[1] = varargin_2->size[1];
    emxEnsureCapacity_real_T(&b_st, Z, i, &wb_emlrtRTEI);
    inputSize[0] = 0.5;
    inputSize[1] = 0.5;
    imlincomb_tbb_real64(inputSize, 2.0, &Z->data[0], 0, static_cast<real_T>
                         ((varargin_2->size[0] * varargin_2->size[1])), 2.0,
                         &varargin_2->data[0], &varargin_4->data[0]);
  } else {
    b_st.site = &cc_emlrtRSI;
    i = Z->size[0] * Z->size[1];
    Z->size[0] = varargin_2->size[0];
    Z->size[1] = varargin_2->size[1];
    emxEnsureCapacity_real_T(&b_st, Z, i, &vb_emlrtRTEI);
    i = varargin_2->size[0] * varargin_2->size[1];
    for (loop_ub = 0; loop_ub < i; loop_ub++) {
      i1 = varargin_2->size[0] * varargin_2->size[1];
      i2 = loop_ub + 1;
      if ((i2 < 1) || (i2 > i1)) {
        emlrtDynamicBoundsCheckR2012b(i2, 1, i1, &g_emlrtBCI, &b_st);
      }

      i1 = varargin_4->size[0] * varargin_4->size[1];
      i2 = static_cast<int32_T>((loop_ub + 1U));
      if ((i2 < 1) || (i2 > i1)) {
        emlrtDynamicBoundsCheckR2012b(i2, 1, i1, &h_emlrtBCI, &b_st);
      }

      i1 = Z->size[0] * Z->size[1];
      if (i2 > i1) {
        emlrtDynamicBoundsCheckR2012b(i2, 1, i1, &i_emlrtBCI, &b_st);
      }

      Z->data[i2 - 1] = 0.5 * varargin_2->data[loop_ub] + 0.5 * varargin_4->
        data[i2 - 1];
    }
  }

  emlrtHeapReferenceStackLeaveFcnR2012b(sp);
}

/* End of code generation (imlincomb.cpp) */
