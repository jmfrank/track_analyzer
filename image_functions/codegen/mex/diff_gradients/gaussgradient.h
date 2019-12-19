/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * gaussgradient.h
 *
 * Code generation for function 'gaussgradient'
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
void gaussgradient(const emlrtStack *sp, const emxArray_real_T *IM, real_T sigma,
                   emxArray_real_T *gx, emxArray_real_T *gy);

/* End of code generation (gaussgradient.h) */
