/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * diffusioncode.h
 *
 * Code generation for function 'diffusioncode'
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
void diffusioncode(const emlrtStack *sp, emxArray_real_T *im, real_T num_iter,
                   real_T kappa1, real_T kappa2, real_T option);

/* End of code generation (diffusioncode.h) */
