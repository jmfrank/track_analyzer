/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * abs.h
 *
 * Code generation for function 'abs'
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
void b_abs(const emlrtStack *sp, const emxArray_real_T *x, emxArray_real_T *y);
void c_abs(const emlrtStack *sp, const emxArray_real_T *x, emxArray_real_T *y);

/* End of code generation (abs.h) */
