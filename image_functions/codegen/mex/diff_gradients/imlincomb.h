/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * imlincomb.h
 *
 * Code generation for function 'imlincomb'
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
void b_imlincomb(const emlrtStack *sp, real_T varargin_1, const emxArray_real_T *
                 varargin_2, real_T varargin_3, emxArray_real_T *Z);
void imlincomb(const emlrtStack *sp, const emxArray_real_T *varargin_2, const
               emxArray_real_T *varargin_4, emxArray_real_T *Z);

/* End of code generation (imlincomb.h) */
