/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * diff_gradients.h
 *
 * Code generation for function 'diff_gradients'
 *
 */

#pragma once

/* Include files */
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include "mex.h"
#include "emlrt.h"
#include "rtwtypes.h"
#include "diff_gradients_types.h"

/* Function Declarations */
CODEGEN_EXPORT_SYM void diff_gradients(const emlrtStack *sp, const
  emxArray_real_T *b_I, const struct0_T *params, emxArray_real_T *J);

/* End of code generation (diff_gradients.h) */
