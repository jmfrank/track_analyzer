/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * _coder_diff_gradients_api.cpp
 *
 * Code generation for function '_coder_diff_gradients_api'
 *
 */

/* Include files */
#include "_coder_diff_gradients_api.h"
#include "diff_gradients.h"
#include "diff_gradients_data.h"
#include "diff_gradients_emxutil.h"
#include "rt_nonfinite.h"
#include <string.h>

/* Variable Definitions */
static emlrtRTEInfo md_emlrtRTEI = { 1,/* lineNo */
  1,                                   /* colNo */
  "_coder_diff_gradients_api",         /* fName */
  ""                                   /* pName */
};

/* Function Declarations */
static void b_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId, emxArray_real_T *y);
static void c_emlrt_marshallIn(const emlrtStack *sp, const mxArray *params,
  const char_T *identifier, struct0_T *y);
static void d_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId, struct0_T *y);
static real_T e_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId);
static void emlrt_marshallIn(const emlrtStack *sp, const mxArray *b_I, const
  char_T *identifier, emxArray_real_T *y);
static const mxArray *emlrt_marshallOut(const emxArray_real_T *u);
static void f_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId, real_T y[3]);
static void g_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId, real_T y[100]);
static void h_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src, const
  emlrtMsgIdentifier *msgId, emxArray_real_T *ret);
static real_T i_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src, const
  emlrtMsgIdentifier *msgId);
static void j_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src, const
  emlrtMsgIdentifier *msgId, real_T ret[3]);
static void k_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src, const
  emlrtMsgIdentifier *msgId, real_T ret[100]);

/* Function Definitions */
static void b_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId, emxArray_real_T *y)
{
  h_emlrt_marshallIn(sp, emlrtAlias(u), parentId, y);
  emlrtDestroyArray(&u);
}

static void c_emlrt_marshallIn(const emlrtStack *sp, const mxArray *params,
  const char_T *identifier, struct0_T *y)
{
  emlrtMsgIdentifier thisId;
  thisId.fIdentifier = const_cast<const char *>(identifier);
  thisId.fParent = NULL;
  thisId.bParentIsCell = false;
  d_emlrt_marshallIn(sp, emlrtAlias(params), &thisId, y);
  emlrtDestroyArray(&params);
}

static void d_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId, struct0_T *y)
{
  emlrtMsgIdentifier thisId;
  static const char * fieldNames[32] = { "seg_channel", "nuclear_channel",
    "AbsMinVol", "AbsMaxVol", "alpha", "beta", "gamma", "epsilon", "delta", "bg",
    "sigmagradient", "noise_filter", "median_filter", "diffuse_iterations",
    "kappa1", "kappa2", "option", "thrshlevel", "percentile", "I_sm_sigma",
    "Hdepth", "nconn_BW", "nconn_BW2", "z_effect", "h_min_depth", "h_min_conn",
    "MeanFilterSensitivity", "MeanFilterNeighborhood", "OutlierThreshold",
    "StrelXYZ", "imclose_r", "rg_threshold" };

  static const int32_T dims = 0;
  thisId.fParent = parentId;
  thisId.bParentIsCell = false;
  emlrtCheckStructR2012b(sp, parentId, u, 32, fieldNames, 0U, &dims);
  thisId.fIdentifier = "seg_channel";
  y->seg_channel = e_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u,
    0, 0, "seg_channel")), &thisId);
  thisId.fIdentifier = "nuclear_channel";
  y->nuclear_channel = e_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp,
    u, 0, 1, "nuclear_channel")), &thisId);
  thisId.fIdentifier = "AbsMinVol";
  y->AbsMinVol = e_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0,
    2, "AbsMinVol")), &thisId);
  thisId.fIdentifier = "AbsMaxVol";
  y->AbsMaxVol = e_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0,
    3, "AbsMaxVol")), &thisId);
  thisId.fIdentifier = "alpha";
  y->alpha = e_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 4,
    "alpha")), &thisId);
  thisId.fIdentifier = "beta";
  y->beta = e_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 5,
    "beta")), &thisId);
  thisId.fIdentifier = "gamma";
  y->gamma = e_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 6,
    "gamma")), &thisId);
  thisId.fIdentifier = "epsilon";
  y->epsilon = e_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 7,
    "epsilon")), &thisId);
  thisId.fIdentifier = "delta";
  y->delta = e_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 8,
    "delta")), &thisId);
  thisId.fIdentifier = "bg";
  y->bg = e_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 9,
    "bg")), &thisId);
  thisId.fIdentifier = "sigmagradient";
  f_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 10,
    "sigmagradient")), &thisId, y->sigmagradient);
  thisId.fIdentifier = "noise_filter";
  g_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 11,
    "noise_filter")), &thisId, y->noise_filter);
  thisId.fIdentifier = "median_filter";
  f_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 12,
    "median_filter")), &thisId, y->median_filter);
  thisId.fIdentifier = "diffuse_iterations";
  y->diffuse_iterations = e_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b
    (sp, u, 0, 13, "diffuse_iterations")), &thisId);
  thisId.fIdentifier = "kappa1";
  y->kappa1 = e_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 14,
    "kappa1")), &thisId);
  thisId.fIdentifier = "kappa2";
  y->kappa2 = e_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 15,
    "kappa2")), &thisId);
  thisId.fIdentifier = "option";
  y->option = e_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 16,
    "option")), &thisId);
  thisId.fIdentifier = "thrshlevel";
  y->thrshlevel = e_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0,
    17, "thrshlevel")), &thisId);
  thisId.fIdentifier = "percentile";
  y->percentile = e_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0,
    18, "percentile")), &thisId);
  thisId.fIdentifier = "I_sm_sigma";
  y->I_sm_sigma = e_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0,
    19, "I_sm_sigma")), &thisId);
  thisId.fIdentifier = "Hdepth";
  y->Hdepth = e_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 20,
    "Hdepth")), &thisId);
  thisId.fIdentifier = "nconn_BW";
  y->nconn_BW = e_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0,
    21, "nconn_BW")), &thisId);
  thisId.fIdentifier = "nconn_BW2";
  y->nconn_BW2 = e_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0,
    22, "nconn_BW2")), &thisId);
  thisId.fIdentifier = "z_effect";
  y->z_effect = e_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0,
    23, "z_effect")), &thisId);
  thisId.fIdentifier = "h_min_depth";
  y->h_min_depth = e_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u,
    0, 24, "h_min_depth")), &thisId);
  thisId.fIdentifier = "h_min_conn";
  y->h_min_conn = e_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0,
    25, "h_min_conn")), &thisId);
  thisId.fIdentifier = "MeanFilterSensitivity";
  y->MeanFilterSensitivity = e_emlrt_marshallIn(sp, emlrtAlias
    (emlrtGetFieldR2017b(sp, u, 0, 26, "MeanFilterSensitivity")), &thisId);
  thisId.fIdentifier = "MeanFilterNeighborhood";
  f_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 27,
    "MeanFilterNeighborhood")), &thisId, y->MeanFilterNeighborhood);
  thisId.fIdentifier = "OutlierThreshold";
  y->OutlierThreshold = e_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp,
    u, 0, 28, "OutlierThreshold")), &thisId);
  thisId.fIdentifier = "StrelXYZ";
  f_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 29, "StrelXYZ")),
                     &thisId, y->StrelXYZ);
  thisId.fIdentifier = "imclose_r";
  y->imclose_r = e_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0,
    30, "imclose_r")), &thisId);
  thisId.fIdentifier = "rg_threshold";
  y->rg_threshold = e_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u,
    0, 31, "rg_threshold")), &thisId);
  emlrtDestroyArray(&u);
}

static real_T e_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId)
{
  real_T y;
  y = i_emlrt_marshallIn(sp, emlrtAlias(u), parentId);
  emlrtDestroyArray(&u);
  return y;
}

static void emlrt_marshallIn(const emlrtStack *sp, const mxArray *b_I, const
  char_T *identifier, emxArray_real_T *y)
{
  emlrtMsgIdentifier thisId;
  thisId.fIdentifier = const_cast<const char *>(identifier);
  thisId.fParent = NULL;
  thisId.bParentIsCell = false;
  b_emlrt_marshallIn(sp, emlrtAlias(b_I), &thisId, y);
  emlrtDestroyArray(&b_I);
}

static const mxArray *emlrt_marshallOut(const emxArray_real_T *u)
{
  const mxArray *y;
  const mxArray *m;
  static const int32_T iv[2] = { 0, 0 };

  y = NULL;
  m = emlrtCreateNumericArray(2, iv, mxDOUBLE_CLASS, mxREAL);
  emlrtMxSetData((mxArray *)m, &u->data[0]);
  emlrtSetDimensions((mxArray *)m, u->size, 2);
  emlrtAssign(&y, m);
  return y;
}

static void f_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId, real_T y[3])
{
  j_emlrt_marshallIn(sp, emlrtAlias(u), parentId, y);
  emlrtDestroyArray(&u);
}

static void g_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId, real_T y[100])
{
  k_emlrt_marshallIn(sp, emlrtAlias(u), parentId, y);
  emlrtDestroyArray(&u);
}

static void h_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src, const
  emlrtMsgIdentifier *msgId, emxArray_real_T *ret)
{
  static const int32_T dims[2] = { -1, -1 };

  const boolean_T bv[2] = { true, true };

  int32_T iv[2];
  int32_T i;
  emlrtCheckVsBuiltInR2012b(sp, msgId, src, "double", false, 2U, dims, &bv[0],
    iv);
  ret->allocatedSize = iv[0] * iv[1];
  i = ret->size[0] * ret->size[1];
  ret->size[0] = iv[0];
  ret->size[1] = iv[1];
  emxEnsureCapacity_real_T(sp, ret, i, (emlrtRTEInfo *)NULL);
  ret->data = (real_T *)emlrtMxGetData(src);
  ret->canFreeData = false;
  emlrtDestroyArray(&src);
}

static real_T i_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src, const
  emlrtMsgIdentifier *msgId)
{
  real_T ret;
  static const int32_T dims = 0;
  emlrtCheckBuiltInR2012b(sp, msgId, src, "double", false, 0U, &dims);
  ret = *(real_T *)emlrtMxGetData(src);
  emlrtDestroyArray(&src);
  return ret;
}

static void j_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src, const
  emlrtMsgIdentifier *msgId, real_T ret[3])
{
  static const int32_T dims[2] = { 1, 3 };

  real_T (*r)[3];
  emlrtCheckBuiltInR2012b(sp, msgId, src, "double", false, 2U, dims);
  r = (real_T (*)[3])emlrtMxGetData(src);
  ret[0] = (*r)[0];
  ret[1] = (*r)[1];
  ret[2] = (*r)[2];
  emlrtDestroyArray(&src);
}

static void k_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src, const
  emlrtMsgIdentifier *msgId, real_T ret[100])
{
  static const int32_T dims[2] = { 10, 10 };

  real_T (*r)[100];
  emlrtCheckBuiltInR2012b(sp, msgId, src, "double", false, 2U, dims);
  r = (real_T (*)[100])emlrtMxGetData(src);
  memcpy(&ret[0], &(*r)[0], 100U * sizeof(real_T));
  emlrtDestroyArray(&src);
}

void diff_gradients_api(const mxArray * const prhs[2], int32_T, const mxArray
  *plhs[1])
{
  emxArray_real_T *b_I;
  emxArray_real_T *J;
  struct0_T params;
  emlrtStack st = { NULL,              /* site */
    NULL,                              /* tls */
    NULL                               /* prev */
  };

  st.tls = emlrtRootTLSGlobal;
  emlrtHeapReferenceStackEnterFcnR2012b(&st);
  emxInit_real_T(&st, &b_I, 2, &md_emlrtRTEI, true);
  emxInit_real_T(&st, &J, 2, &md_emlrtRTEI, true);

  /* Marshall function inputs */
  b_I->canFreeData = false;
  emlrt_marshallIn(&st, emlrtAlias(prhs[0]), "I", b_I);
  c_emlrt_marshallIn(&st, emlrtAliasP(prhs[1]), "params", &params);

  /* Invoke the target function */
  diff_gradients(&st, b_I, &params, J);

  /* Marshall function outputs */
  J->canFreeData = false;
  plhs[0] = emlrt_marshallOut(J);
  emxFree_real_T(&J);
  emxFree_real_T(&b_I);
  emlrtHeapReferenceStackLeaveFcnR2012b(&st);
}

/* End of code generation (_coder_diff_gradients_api.cpp) */
