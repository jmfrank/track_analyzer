/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * diff_gradients_initialize.cpp
 *
 * Code generation for function 'diff_gradients_initialize'
 *
 */

/* Include files */
#include "diff_gradients_initialize.h"
#include "_coder_diff_gradients_mex.h"
#include "diff_gradients.h"
#include "diff_gradients_data.h"
#include "rt_nonfinite.h"

/* Function Definitions */
void diff_gradients_initialize()
{
  emlrtStack st = { NULL,              /* site */
    NULL,                              /* tls */
    NULL                               /* prev */
  };

  mex_InitInfAndNan();
  mexFunctionCreateRootTLS();
  emlrtBreakCheckR2012bFlagVar = emlrtGetBreakCheckFlagAddressR2012b();
  st.tls = emlrtRootTLSGlobal;
  emlrtClearAllocCountR2012b(&st, false, 0U, 0);
  emlrtEnterRtStackR2012b(&st);
  emlrtLicenseCheckR2012b(&st, "Image_Toolbox", 2);
  emlrtFirstTimeR2012b(emlrtRootTLSGlobal);
}

/* End of code generation (diff_gradients_initialize.cpp) */
