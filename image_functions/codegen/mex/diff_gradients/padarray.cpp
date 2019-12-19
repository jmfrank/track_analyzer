/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * padarray.cpp
 *
 * Code generation for function 'padarray'
 *
 */

/* Include files */
#include "padarray.h"
#include "diff_gradients.h"
#include "diff_gradients_data.h"
#include "diff_gradients_emxutil.h"
#include "eml_int_forloop_overflow_check.h"
#include "mwmathutil.h"
#include "repmat.h"
#include "rt_nonfinite.h"
#include "validateattributes.h"

/* Variable Definitions */
static emlrtRSInfo tb_emlrtRSI = { 72, /* lineNo */
  "padarray",                          /* fcnName */
  "/usr/local/MATLAB/R2019b/toolbox/images/images/eml/padarray.m"/* pathName */
};

static emlrtRSInfo vb_emlrtRSI = { 405,/* lineNo */
  "ConstantPad",                       /* fcnName */
  "/usr/local/MATLAB/R2019b/toolbox/images/images/eml/padarray.m"/* pathName */
};

static emlrtRSInfo wb_emlrtRSI = { 420,/* lineNo */
  "ConstantPad",                       /* fcnName */
  "/usr/local/MATLAB/R2019b/toolbox/images/images/eml/padarray.m"/* pathName */
};

static emlrtRTEInfo e_emlrtRTEI = { 398,/* lineNo */
  17,                                  /* colNo */
  "ConstantPad",                       /* fName */
  "/usr/local/MATLAB/R2019b/toolbox/images/images/eml/padarray.m"/* pName */
};

static emlrtRTEInfo f_emlrtRTEI = { 413,/* lineNo */
  21,                                  /* colNo */
  "ConstantPad",                       /* fName */
  "/usr/local/MATLAB/R2019b/toolbox/images/images/eml/padarray.m"/* pName */
};

static emlrtBCInfo emlrtBCI = { -1,    /* iFirst */
  -1,                                  /* iLast */
  400,                                 /* lineNo */
  17,                                  /* colNo */
  "",                                  /* aName */
  "ConstantPad",                       /* fName */
  "/usr/local/MATLAB/R2019b/toolbox/images/images/eml/padarray.m",/* pName */
  0                                    /* checkKind */
};

static emlrtBCInfo b_emlrtBCI = { -1,  /* iFirst */
  -1,                                  /* iLast */
  407,                                 /* lineNo */
  17,                                  /* colNo */
  "",                                  /* aName */
  "ConstantPad",                       /* fName */
  "/usr/local/MATLAB/R2019b/toolbox/images/images/eml/padarray.m",/* pName */
  0                                    /* checkKind */
};

static emlrtBCInfo c_emlrtBCI = { -1,  /* iFirst */
  -1,                                  /* iLast */
  414,                                 /* lineNo */
  17,                                  /* colNo */
  "",                                  /* aName */
  "ConstantPad",                       /* fName */
  "/usr/local/MATLAB/R2019b/toolbox/images/images/eml/padarray.m",/* pName */
  0                                    /* checkKind */
};

static emlrtBCInfo d_emlrtBCI = { -1,  /* iFirst */
  -1,                                  /* iLast */
  444,                                 /* lineNo */
  100,                                 /* colNo */
  "",                                  /* aName */
  "ConstantPad",                       /* fName */
  "/usr/local/MATLAB/R2019b/toolbox/images/images/eml/padarray.m",/* pName */
  0                                    /* checkKind */
};

static emlrtBCInfo e_emlrtBCI = { -1,  /* iFirst */
  -1,                                  /* iLast */
  444,                                 /* lineNo */
  17,                                  /* colNo */
  "",                                  /* aName */
  "ConstantPad",                       /* fName */
  "/usr/local/MATLAB/R2019b/toolbox/images/images/eml/padarray.m",/* pName */
  0                                    /* checkKind */
};

static emlrtBCInfo f_emlrtBCI = { -1,  /* iFirst */
  -1,                                  /* iLast */
  421,                                 /* lineNo */
  17,                                  /* colNo */
  "",                                  /* aName */
  "ConstantPad",                       /* fName */
  "/usr/local/MATLAB/R2019b/toolbox/images/images/eml/padarray.m",/* pName */
  0                                    /* checkKind */
};

static emlrtDCInfo c_emlrtDCI = { 253, /* lineNo */
  35,                                  /* colNo */
  "ConstantPad",                       /* fName */
  "/usr/local/MATLAB/R2019b/toolbox/images/images/eml/padarray.m",/* pName */
  1                                    /* checkKind */
};

static emlrtRTEInfo sb_emlrtRTEI = { 72,/* lineNo */
  13,                                  /* colNo */
  "padarray",                          /* fName */
  "/usr/local/MATLAB/R2019b/toolbox/images/images/eml/padarray.m"/* pName */
};

/* Function Definitions */
void padarray(const emlrtStack *sp, const emxArray_real_T *varargin_1, const
              real_T varargin_2[2], emxArray_real_T *b)
{
  real_T sizeB[2];
  real_T d;
  real_T d1;
  int32_T i;
  int32_T j;
  int32_T a;
  int32_T i1;
  int32_T b_b;
  int32_T b_i;
  boolean_T overflow;
  int32_T i2;
  int32_T a_tmp;
  int32_T i3;
  emlrtStack st;
  emlrtStack b_st;
  emlrtStack c_st;
  st.prev = sp;
  st.tls = sp->tls;
  st.site = &rb_emlrtRSI;
  b_st.prev = &st;
  b_st.tls = st.tls;
  c_st.prev = &b_st;
  c_st.tls = b_st.tls;
  b_validateattributes(&st, varargin_2);
  if ((varargin_1->size[0] == 0) || (varargin_1->size[1] == 0)) {
    sizeB[0] = static_cast<real_T>(varargin_1->size[0]) + 2.0 * varargin_2[0];
    sizeB[1] = static_cast<real_T>(varargin_1->size[1]) + 2.0 * varargin_2[1];
    st.site = &sb_emlrtRSI;
    repmat(&st, sizeB, b);
  } else {
    st.site = &tb_emlrtRSI;
    d = static_cast<real_T>(varargin_1->size[0]) + 2.0 * varargin_2[0];
    if (d != static_cast<int32_T>(muDoubleScalarFloor(d))) {
      emlrtIntegerCheckR2012b(d, &c_emlrtDCI, &st);
    }

    d1 = static_cast<real_T>(varargin_1->size[1]) + 2.0 * varargin_2[1];
    if (d1 != static_cast<int32_T>(muDoubleScalarFloor(d1))) {
      emlrtIntegerCheckR2012b(d1, &c_emlrtDCI, &st);
    }

    i = b->size[0] * b->size[1];
    b->size[0] = static_cast<int32_T>(d);
    b->size[1] = static_cast<int32_T>(d1);
    emxEnsureCapacity_real_T(&st, b, i, &sb_emlrtRTEI);
    i = static_cast<int32_T>(varargin_2[1]);
    emlrtForLoopVectorCheckR2012b(1.0, 1.0, varargin_2[1], mxDOUBLE_CLASS,
      static_cast<int32_T>(varargin_2[1]), &e_emlrtRTEI, &st);
    for (j = 0; j < i; j++) {
      i1 = b->size[0];
      for (b_i = 0; b_i < i1; b_i++) {
        i2 = b_i + 1;
        if ((i2 < 1) || (i2 > b->size[0])) {
          emlrtDynamicBoundsCheckR2012b(i2, 1, b->size[0], &emlrtBCI, &st);
        }

        a_tmp = static_cast<int32_T>((j + 1U));
        if ((a_tmp < 1) || (a_tmp > b->size[1])) {
          emlrtDynamicBoundsCheckR2012b(a_tmp, 1, b->size[1], &emlrtBCI, &st);
        }

        b->data[(i2 + b->size[0] * (a_tmp - 1)) - 1] = 0.0;
      }
    }

    a = (varargin_1->size[1] + i) + 1;
    b_b = b->size[1];
    b_st.site = &vb_emlrtRSI;
    overflow = (((varargin_1->size[1] + i) + 1 <= b->size[1]) && (b->size[1] >
      2147483646));
    if (overflow) {
      c_st.site = &eb_emlrtRSI;
      check_forloop_overflow_error(&c_st);
    }

    for (j = a; j <= b_b; j++) {
      i1 = b->size[0];
      for (b_i = 0; b_i < i1; b_i++) {
        i2 = b_i + 1;
        if ((i2 < 1) || (i2 > b->size[0])) {
          emlrtDynamicBoundsCheckR2012b(i2, 1, b->size[0], &b_emlrtBCI, &st);
        }

        if ((j < 1) || (j > b->size[1])) {
          emlrtDynamicBoundsCheckR2012b(j, 1, b->size[1], &b_emlrtBCI, &st);
        }

        b->data[(i2 + b->size[0] * (j - 1)) - 1] = 0.0;
      }
    }

    i1 = varargin_1->size[1];
    for (j = 0; j < i1; j++) {
      i2 = static_cast<int32_T>(varargin_2[0]);
      emlrtForLoopVectorCheckR2012b(1.0, 1.0, varargin_2[0], mxDOUBLE_CLASS,
        static_cast<int32_T>(varargin_2[0]), &f_emlrtRTEI, &st);
      for (b_i = 0; b_i < i2; b_i++) {
        a_tmp = static_cast<int32_T>((b_i + 1U));
        if ((a_tmp < 1) || (a_tmp > b->size[0])) {
          emlrtDynamicBoundsCheckR2012b(a_tmp, 1, b->size[0], &c_emlrtBCI, &st);
        }

        a = (j + i) + 1;
        if ((a < 1) || (a > b->size[1])) {
          emlrtDynamicBoundsCheckR2012b(a, 1, b->size[1], &c_emlrtBCI, &st);
        }

        b->data[(a_tmp + b->size[0] * (a - 1)) - 1] = 0.0;
      }
    }

    i1 = varargin_1->size[1];
    for (j = 0; j < i1; j++) {
      a_tmp = static_cast<int32_T>(varargin_2[0]);
      a = (a_tmp + varargin_1->size[0]) + 1;
      b_b = b->size[0];
      b_st.site = &wb_emlrtRSI;
      overflow = (((a_tmp + varargin_1->size[0]) + 1 <= b->size[0]) && (b->size
        [0] > 2147483646));
      if (overflow) {
        c_st.site = &eb_emlrtRSI;
        check_forloop_overflow_error(&c_st);
      }

      for (b_i = a; b_i <= b_b; b_i++) {
        if ((b_i < 1) || (b_i > b->size[0])) {
          emlrtDynamicBoundsCheckR2012b(b_i, 1, b->size[0], &f_emlrtBCI, &st);
        }

        i2 = (j + i) + 1;
        if ((i2 < 1) || (i2 > b->size[1])) {
          emlrtDynamicBoundsCheckR2012b(i2, 1, b->size[1], &f_emlrtBCI, &st);
        }

        b->data[(b_i + b->size[0] * (i2 - 1)) - 1] = 0.0;
      }
    }

    i1 = varargin_1->size[1];
    for (j = 0; j < i1; j++) {
      i2 = varargin_1->size[0];
      for (b_i = 0; b_i < i2; b_i++) {
        a_tmp = b_i + 1;
        if ((a_tmp < 1) || (a_tmp > varargin_1->size[0])) {
          emlrtDynamicBoundsCheckR2012b(a_tmp, 1, varargin_1->size[0],
            &d_emlrtBCI, &st);
        }

        a = j + 1;
        if ((a < 1) || (a > varargin_1->size[1])) {
          emlrtDynamicBoundsCheckR2012b(a, 1, varargin_1->size[1], &d_emlrtBCI,
            &st);
        }

        b_b = (b_i + static_cast<int32_T>(varargin_2[0])) + 1;
        if ((b_b < 1) || (b_b > b->size[0])) {
          emlrtDynamicBoundsCheckR2012b(b_b, 1, b->size[0], &e_emlrtBCI, &st);
        }

        i3 = (j + i) + 1;
        if ((i3 < 1) || (i3 > b->size[1])) {
          emlrtDynamicBoundsCheckR2012b(i3, 1, b->size[1], &e_emlrtBCI, &st);
        }

        b->data[(b_b + b->size[0] * (i3 - 1)) - 1] = varargin_1->data[(a_tmp +
          varargin_1->size[0] * (a - 1)) - 1];
      }
    }
  }
}

/* End of code generation (padarray.cpp) */
