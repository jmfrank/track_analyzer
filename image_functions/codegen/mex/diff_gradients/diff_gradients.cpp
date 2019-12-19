/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * diff_gradients.cpp
 *
 * Code generation for function 'diff_gradients'
 *
 */

/* Include files */
#include "diff_gradients.h"
#include "abs.h"
#include "diff_gradients_data.h"
#include "diff_gradients_emxutil.h"
#include "diffusioncode.h"
#include "eml_int_forloop_overflow_check.h"
#include "gaussgradient.h"
#include "imlincomb.h"
#include "mat2gray.h"
#include "mwmathutil.h"
#include "ordfilt2.h"
#include "rt_nonfinite.h"
#include "scalexpAlloc.h"

/* Variable Definitions */
static emlrtRSInfo emlrtRSI = { 7,     /* lineNo */
  "diff_gradients",                    /* fcnName */
  "/media/data/Dropbox/code_bank/matlab/sandbox_guan/cell_tracking/track_analyzer/image_functions/diff_gradients.m"/* pathName */
};

static emlrtRSInfo b_emlrtRSI = { 12,  /* lineNo */
  "diff_gradients",                    /* fcnName */
  "/media/data/Dropbox/code_bank/matlab/sandbox_guan/cell_tracking/track_analyzer/image_functions/diff_gradients.m"/* pathName */
};

static emlrtRSInfo c_emlrtRSI = { 15,  /* lineNo */
  "diff_gradients",                    /* fcnName */
  "/media/data/Dropbox/code_bank/matlab/sandbox_guan/cell_tracking/track_analyzer/image_functions/diff_gradients.m"/* pathName */
};

static emlrtRSInfo d_emlrtRSI = { 16,  /* lineNo */
  "diff_gradients",                    /* fcnName */
  "/media/data/Dropbox/code_bank/matlab/sandbox_guan/cell_tracking/track_analyzer/image_functions/diff_gradients.m"/* pathName */
};

static emlrtRSInfo e_emlrtRSI = { 17,  /* lineNo */
  "diff_gradients",                    /* fcnName */
  "/media/data/Dropbox/code_bank/matlab/sandbox_guan/cell_tracking/track_analyzer/image_functions/diff_gradients.m"/* pathName */
};

static emlrtRSInfo f_emlrtRSI = { 20,  /* lineNo */
  "diff_gradients",                    /* fcnName */
  "/media/data/Dropbox/code_bank/matlab/sandbox_guan/cell_tracking/track_analyzer/image_functions/diff_gradients.m"/* pathName */
};

static emlrtRSInfo g_emlrtRSI = { 21,  /* lineNo */
  "diff_gradients",                    /* fcnName */
  "/media/data/Dropbox/code_bank/matlab/sandbox_guan/cell_tracking/track_analyzer/image_functions/diff_gradients.m"/* pathName */
};

static emlrtRSInfo h_emlrtRSI = { 29,  /* lineNo */
  "diff_gradients",                    /* fcnName */
  "/media/data/Dropbox/code_bank/matlab/sandbox_guan/cell_tracking/track_analyzer/image_functions/diff_gradients.m"/* pathName */
};

static emlrtRSInfo i_emlrtRSI = { 30,  /* lineNo */
  "diff_gradients",                    /* fcnName */
  "/media/data/Dropbox/code_bank/matlab/sandbox_guan/cell_tracking/track_analyzer/image_functions/diff_gradients.m"/* pathName */
};

static emlrtRSInfo j_emlrtRSI = { 14,  /* lineNo */
  "medfilt2",                          /* fcnName */
  "/usr/local/MATLAB/R2019b/toolbox/images/images/eml/medfilt2.m"/* pathName */
};

static emlrtRSInfo k_emlrtRSI = { 36,  /* lineNo */
  "medfilt2",                          /* fcnName */
  "/usr/local/MATLAB/R2019b/toolbox/images/images/eml/medfilt2.m"/* pathName */
};

static emlrtRSInfo l_emlrtRSI = { 60,  /* lineNo */
  "medfilt2",                          /* fcnName */
  "/usr/local/MATLAB/R2019b/toolbox/images/images/eml/medfilt2.m"/* pathName */
};

static emlrtRSInfo m_emlrtRSI = { 61,  /* lineNo */
  "medfilt2",                          /* fcnName */
  "/usr/local/MATLAB/R2019b/toolbox/images/images/eml/medfilt2.m"/* pathName */
};

static emlrtRSInfo n_emlrtRSI = { 74,  /* lineNo */
  "medfilt2",                          /* fcnName */
  "/usr/local/MATLAB/R2019b/toolbox/images/images/eml/medfilt2.m"/* pathName */
};

static emlrtRSInfo o_emlrtRSI = { 151, /* lineNo */
  "parseInputs",                       /* fcnName */
  "/usr/local/MATLAB/R2019b/toolbox/images/images/eml/medfilt2.m"/* pathName */
};

static emlrtRSInfo p_emlrtRSI = { 204, /* lineNo */
  "parseMN",                           /* fcnName */
  "/usr/local/MATLAB/R2019b/toolbox/images/images/eml/medfilt2.m"/* pathName */
};

static emlrtRSInfo vg_emlrtRSI = { 13, /* lineNo */
  "hypot",                             /* fcnName */
  "/usr/local/MATLAB/R2019b/toolbox/eml/lib/matlab/elfun/hypot.m"/* pathName */
};

static emlrtRSInfo wg_emlrtRSI = { 46, /* lineNo */
  "applyBinaryScalarFunction",         /* fcnName */
  "/usr/local/MATLAB/R2019b/toolbox/eml/eml/+coder/+internal/applyBinaryScalarFunction.m"/* pathName */
};

static emlrtRSInfo xg_emlrtRSI = { 202,/* lineNo */
  "flatIter",                          /* fcnName */
  "/usr/local/MATLAB/R2019b/toolbox/eml/eml/+coder/+internal/applyBinaryScalarFunction.m"/* pathName */
};

static emlrtRSInfo yg_emlrtRSI = { 11, /* lineNo */
  "tanh",                              /* fcnName */
  "/usr/local/MATLAB/R2019b/toolbox/eml/lib/matlab/elfun/tanh.m"/* pathName */
};

static emlrtECInfo emlrtECI = { 2,     /* nDims */
  22,                                  /* lineNo */
  15,                                  /* colNo */
  "diff_gradients",                    /* fName */
  "/media/data/Dropbox/code_bank/matlab/sandbox_guan/cell_tracking/track_analyzer/image_functions/diff_gradients.m"/* pName */
};

static emlrtECInfo b_emlrtECI = { 2,   /* nDims */
  23,                                  /* lineNo */
  13,                                  /* colNo */
  "diff_gradients",                    /* fName */
  "/media/data/Dropbox/code_bank/matlab/sandbox_guan/cell_tracking/track_analyzer/image_functions/diff_gradients.m"/* pName */
};

static emlrtECInfo c_emlrtECI = { 2,   /* nDims */
  26,                                  /* lineNo */
  17,                                  /* colNo */
  "diff_gradients",                    /* fName */
  "/media/data/Dropbox/code_bank/matlab/sandbox_guan/cell_tracking/track_analyzer/image_functions/diff_gradients.m"/* pName */
};

static emlrtECInfo d_emlrtECI = { 2,   /* nDims */
  26,                                  /* lineNo */
  32,                                  /* colNo */
  "diff_gradients",                    /* fName */
  "/media/data/Dropbox/code_bank/matlab/sandbox_guan/cell_tracking/track_analyzer/image_functions/diff_gradients.m"/* pName */
};

static emlrtECInfo e_emlrtECI = { 2,   /* nDims */
  27,                                  /* lineNo */
  17,                                  /* colNo */
  "diff_gradients",                    /* fName */
  "/media/data/Dropbox/code_bank/matlab/sandbox_guan/cell_tracking/track_analyzer/image_functions/diff_gradients.m"/* pName */
};

static emlrtECInfo f_emlrtECI = { 2,   /* nDims */
  29,                                  /* lineNo */
  22,                                  /* colNo */
  "diff_gradients",                    /* fName */
  "/media/data/Dropbox/code_bank/matlab/sandbox_guan/cell_tracking/track_analyzer/image_functions/diff_gradients.m"/* pName */
};

static emlrtECInfo g_emlrtECI = { 2,   /* nDims */
  31,                                  /* lineNo */
  8,                                   /* colNo */
  "diff_gradients",                    /* fName */
  "/media/data/Dropbox/code_bank/matlab/sandbox_guan/cell_tracking/track_analyzer/image_functions/diff_gradients.m"/* pName */
};

static emlrtRTEInfo emlrtRTEI = { 19,  /* lineNo */
  23,                                  /* colNo */
  "scalexpAlloc",                      /* fName */
  "/usr/local/MATLAB/R2019b/toolbox/eml/eml/+coder/+internal/scalexpAlloc.m"/* pName */
};

static emlrtRTEInfo b_emlrtRTEI = { 14,/* lineNo */
  37,                                  /* colNo */
  "validatepositive",                  /* fName */
  "/usr/local/MATLAB/R2019b/toolbox/eml/eml/+coder/+internal/+valattr/validatepositive.m"/* pName */
};

static emlrtDCInfo emlrtDCI = { 29,    /* lineNo */
  19,                                  /* colNo */
  "medfilt2",                          /* fName */
  "/usr/local/MATLAB/R2019b/toolbox/images/images/eml/medfilt2.m",/* pName */
  1                                    /* checkKind */
};

static emlrtDCInfo b_emlrtDCI = { 29,  /* lineNo */
  19,                                  /* colNo */
  "medfilt2",                          /* fName */
  "/usr/local/MATLAB/R2019b/toolbox/images/images/eml/medfilt2.m",/* pName */
  4                                    /* checkKind */
};

static emlrtRTEInfo db_emlrtRTEI = { 7,/* lineNo */
  5,                                   /* colNo */
  "diff_gradients",                    /* fName */
  "/media/data/Dropbox/code_bank/matlab/sandbox_guan/cell_tracking/track_analyzer/image_functions/diff_gradients.m"/* pName */
};

static emlrtRTEInfo eb_emlrtRTEI = { 29,/* lineNo */
  5,                                   /* colNo */
  "medfilt2",                          /* fName */
  "/usr/local/MATLAB/R2019b/toolbox/images/images/eml/medfilt2.m"/* pName */
};

static emlrtRTEInfo fb_emlrtRTEI = { 12,/* lineNo */
  5,                                   /* colNo */
  "diff_gradients",                    /* fName */
  "/media/data/Dropbox/code_bank/matlab/sandbox_guan/cell_tracking/track_analyzer/image_functions/diff_gradients.m"/* pName */
};

static emlrtRTEInfo gb_emlrtRTEI = { 60,/* lineNo */
  17,                                  /* colNo */
  "medfilt2",                          /* fName */
  "/usr/local/MATLAB/R2019b/toolbox/images/images/eml/medfilt2.m"/* pName */
};

static emlrtRTEInfo hb_emlrtRTEI = { 61,/* lineNo */
  17,                                  /* colNo */
  "medfilt2",                          /* fName */
  "/usr/local/MATLAB/R2019b/toolbox/images/images/eml/medfilt2.m"/* pName */
};

static emlrtRTEInfo ib_emlrtRTEI = { 46,/* lineNo */
  6,                                   /* colNo */
  "applyBinaryScalarFunction",         /* fName */
  "/usr/local/MATLAB/R2019b/toolbox/eml/eml/+coder/+internal/applyBinaryScalarFunction.m"/* pName */
};

static emlrtRTEInfo jb_emlrtRTEI = { 13,/* lineNo */
  5,                                   /* colNo */
  "hypot",                             /* fName */
  "/usr/local/MATLAB/R2019b/toolbox/eml/lib/matlab/elfun/hypot.m"/* pName */
};

static emlrtRTEInfo kb_emlrtRTEI = { 22,/* lineNo */
  5,                                   /* colNo */
  "diff_gradients",                    /* fName */
  "/media/data/Dropbox/code_bank/matlab/sandbox_guan/cell_tracking/track_analyzer/image_functions/diff_gradients.m"/* pName */
};

static emlrtRTEInfo lb_emlrtRTEI = { 23,/* lineNo */
  23,                                  /* colNo */
  "diff_gradients",                    /* fName */
  "/media/data/Dropbox/code_bank/matlab/sandbox_guan/cell_tracking/track_analyzer/image_functions/diff_gradients.m"/* pName */
};

static emlrtRTEInfo mb_emlrtRTEI = { 27,/* lineNo */
  31,                                  /* colNo */
  "diff_gradients",                    /* fName */
  "/media/data/Dropbox/code_bank/matlab/sandbox_guan/cell_tracking/track_analyzer/image_functions/diff_gradients.m"/* pName */
};

static emlrtRTEInfo nb_emlrtRTEI = { 17,/* lineNo */
  5,                                   /* colNo */
  "diff_gradients",                    /* fName */
  "/media/data/Dropbox/code_bank/matlab/sandbox_guan/cell_tracking/track_analyzer/image_functions/diff_gradients.m"/* pName */
};

static emlrtRTEInfo ob_emlrtRTEI = { 26,/* lineNo */
  5,                                   /* colNo */
  "diff_gradients",                    /* fName */
  "/media/data/Dropbox/code_bank/matlab/sandbox_guan/cell_tracking/track_analyzer/image_functions/diff_gradients.m"/* pName */
};

static emlrtRTEInfo pb_emlrtRTEI = { 16,/* lineNo */
  10,                                  /* colNo */
  "diff_gradients",                    /* fName */
  "/media/data/Dropbox/code_bank/matlab/sandbox_guan/cell_tracking/track_analyzer/image_functions/diff_gradients.m"/* pName */
};

static emlrtRTEInfo qb_emlrtRTEI = { 2,/* lineNo */
  14,                                  /* colNo */
  "diff_gradients",                    /* fName */
  "/media/data/Dropbox/code_bank/matlab/sandbox_guan/cell_tracking/track_analyzer/image_functions/diff_gradients.m"/* pName */
};

/* Function Definitions */
void diff_gradients(const emlrtStack *sp, const emxArray_real_T *b_I, const
                    struct0_T *params, emxArray_real_T *J)
{
  boolean_T overflow;
  int32_T csz_idx_0;
  boolean_T exitg1;
  static const int8_T iv[2] = { 0, 1 };

  emxArray_real_T *domain;
  emxArray_real_T *b1;
  emxArray_real_T *b2;
  int32_T i;
  int32_T csz_idx_1;
  real_T order1;
  emxArray_real_T *imy;
  emxArray_real_T *Mag_fim;
  emxArray_real_T *Det_hessian;
  emxArray_real_T *L_imxy;
  emxArray_boolean_T *r;
  emlrtStack st;
  emlrtStack b_st;
  emlrtStack c_st;
  emlrtStack d_st;
  emlrtStack e_st;
  st.prev = sp;
  st.tls = sp->tls;
  b_st.prev = &st;
  b_st.tls = st.tls;
  c_st.prev = &b_st;
  c_st.tls = b_st.tls;
  d_st.prev = &c_st;
  d_st.tls = c_st.tls;
  e_st.prev = &d_st;
  e_st.tls = d_st.tls;
  emlrtHeapReferenceStackEnterFcnR2012b(sp);

  /*  diffusion gradients */
  /*     %%  Step1: Denoising Filters */
  /* G(:,:,zslice) = imfilter(I(:,:,zslice), h,'replicate'); */
  /* Optional filters median or deconvlucy */
  st.site = &emlrtRSI;
  b_st.site = &j_emlrtRSI;
  c_st.site = &o_emlrtRSI;
  d_st.site = &p_emlrtRSI;
  e_st.site = &q_emlrtRSI;
  overflow = true;
  csz_idx_0 = 0;
  exitg1 = false;
  while ((!exitg1) && (csz_idx_0 < 2)) {
    if (!(params->median_filter[iv[csz_idx_0]] <= 0.0)) {
      csz_idx_0++;
    } else {
      overflow = false;
      exitg1 = true;
    }
  }

  if (!overflow) {
    emlrtErrorWithMessageIdR2018a(&e_st, &b_emlrtRTEI,
      "Coder:toolbox:ValidateattributesexpectedPositive",
      "MATLAB:medfilt2:expectedPositive", 3, 4, 22, "input number 2, [M N],");
  }

  e_st.site = &q_emlrtRSI;
  overflow = true;
  csz_idx_0 = 0;
  exitg1 = false;
  while ((!exitg1) && (csz_idx_0 < 2)) {
    if ((!muDoubleScalarIsInf(params->median_filter[csz_idx_0])) &&
        (!muDoubleScalarIsNaN(params->median_filter[csz_idx_0])) &&
        (muDoubleScalarFloor(params->median_filter[csz_idx_0]) ==
         params->median_filter[csz_idx_0])) {
      csz_idx_0++;
    } else {
      overflow = false;
      exitg1 = true;
    }
  }

  if (!overflow) {
    emlrtErrorWithMessageIdR2018a(&e_st, &c_emlrtRTEI,
      "Coder:toolbox:ValidateattributesexpectedInteger",
      "MATLAB:medfilt2:expectedInteger", 3, 4, 22, "input number 2, [M N],");
  }

  emxInit_real_T(&st, &domain, 2, &eb_emlrtRTEI, true);
  emxInit_real_T(&st, &b1, 2, &gb_emlrtRTEI, true);
  emxInit_real_T(&st, &b2, 2, &hb_emlrtRTEI, true);
  if ((b_I->size[0] == 0) || (b_I->size[1] == 0)) {
    i = J->size[0] * J->size[1];
    J->size[0] = b_I->size[0];
    J->size[1] = b_I->size[1];
    emxEnsureCapacity_real_T(&st, J, i, &db_emlrtRTEI);
    csz_idx_0 = b_I->size[0] * b_I->size[1];
    for (i = 0; i < csz_idx_0; i++) {
      J->data[i] = b_I->data[i];
    }
  } else {
    if (!(params->median_filter[0] >= 0.0)) {
      emlrtNonNegativeCheckR2012b(params->median_filter[0], &b_emlrtDCI, &st);
    }

    if (params->median_filter[0] != static_cast<int32_T>(muDoubleScalarFloor
         (params->median_filter[0]))) {
      emlrtIntegerCheckR2012b(params->median_filter[0], &emlrtDCI, &st);
    }

    if (!(params->median_filter[1] >= 0.0)) {
      emlrtNonNegativeCheckR2012b(params->median_filter[1], &b_emlrtDCI, &st);
    }

    if (params->median_filter[1] != static_cast<int32_T>(muDoubleScalarFloor
         (params->median_filter[1]))) {
      emlrtIntegerCheckR2012b(params->median_filter[1], &emlrtDCI, &st);
    }

    i = static_cast<int32_T>(params->median_filter[0]);
    csz_idx_0 = domain->size[0] * domain->size[1];
    domain->size[0] = i;
    csz_idx_1 = static_cast<int32_T>(params->median_filter[1]);
    domain->size[1] = csz_idx_1;
    emxEnsureCapacity_real_T(&st, domain, csz_idx_0, &eb_emlrtRTEI);
    csz_idx_0 = i * csz_idx_1;
    for (i = 0; i < csz_idx_0; i++) {
      domain->data[i] = 1.0;
    }

    order1 = params->median_filter[0] * params->median_filter[1];
    if (muDoubleScalarRem(order1, 2.0) == 1.0) {
      i = J->size[0] * J->size[1];
      J->size[0] = b_I->size[0];
      J->size[1] = b_I->size[1];
      emxEnsureCapacity_real_T(&st, J, i, &db_emlrtRTEI);
      csz_idx_0 = b_I->size[0] * b_I->size[1];
      for (i = 0; i < csz_idx_0; i++) {
        J->data[i] = b_I->data[i];
      }

      b_st.site = &k_emlrtRSI;
      ordfilt2(&b_st, J, (order1 + 1.0) / 2.0, domain);
    } else {
      order1 /= 2.0;
      i = b1->size[0] * b1->size[1];
      b1->size[0] = b_I->size[0];
      b1->size[1] = b_I->size[1];
      emxEnsureCapacity_real_T(&st, b1, i, &gb_emlrtRTEI);
      csz_idx_0 = b_I->size[0] * b_I->size[1];
      for (i = 0; i < csz_idx_0; i++) {
        b1->data[i] = b_I->data[i];
      }

      b_st.site = &l_emlrtRSI;
      ordfilt2(&b_st, b1, order1, domain);
      i = b2->size[0] * b2->size[1];
      b2->size[0] = b_I->size[0];
      b2->size[1] = b_I->size[1];
      emxEnsureCapacity_real_T(&st, b2, i, &hb_emlrtRTEI);
      csz_idx_0 = b_I->size[0] * b_I->size[1];
      for (i = 0; i < csz_idx_0; i++) {
        b2->data[i] = b_I->data[i];
      }

      b_st.site = &m_emlrtRSI;
      ordfilt2(&b_st, b2, order1 + 1.0, domain);
      b_st.site = &n_emlrtRSI;
      imlincomb(&b_st, b1, b2, J);
    }
  }

  /* G(:,:,zslice) = deconvlucy(G(:,:,zslice),h); %%use the same 'h' as for Gaussian or define new variable */
  /* Step2: Perona & Malik non-linear isotropic diffusion filter, refer to */
  /* diffusioncode.m for details, here alter only diffuse_iterations & kappa1        */
  i = domain->size[0] * domain->size[1];
  domain->size[0] = b_I->size[0];
  domain->size[1] = b_I->size[1];
  emxEnsureCapacity_real_T(sp, domain, i, &fb_emlrtRTEI);
  csz_idx_0 = b_I->size[0] * b_I->size[1];
  for (i = 0; i < csz_idx_0; i++) {
    domain->data[i] = b_I->data[i];
  }

  emxInit_real_T(sp, &imy, 2, &pb_emlrtRTEI, true);
  st.site = &b_emlrtRSI;
  diffusioncode(&st, domain, params->diffuse_iterations, params->kappa1,
                params->kappa2, params->option);

  /* Step3a: Gauss gradient, 1st derivative of the image, refer to gaussgradient.m for details */
  st.site = &c_emlrtRSI;
  mat2gray(&st, domain, b1);
  st.site = &d_emlrtRSI;
  gaussgradient(&st, b1, params->sigmagradient[0], b2, imy);
  st.site = &e_emlrtRSI;
  b_st.site = &vg_emlrtRSI;
  c_st.site = &wg_emlrtRSI;
  if (b2->size[0] <= imy->size[0]) {
    csz_idx_0 = b2->size[0];
  } else {
    csz_idx_0 = imy->size[0];
  }

  if (b2->size[1] <= imy->size[1]) {
    csz_idx_1 = b2->size[1];
  } else {
    csz_idx_1 = imy->size[1];
  }

  i = domain->size[0] * domain->size[1];
  domain->size[0] = csz_idx_0;
  domain->size[1] = csz_idx_1;
  emxEnsureCapacity_real_T(&c_st, domain, i, &ib_emlrtRTEI);
  if (!dimagree(domain, b2, imy)) {
    emlrtErrorWithMessageIdR2018a(&c_st, &emlrtRTEI, "MATLAB:dimagree",
      "MATLAB:dimagree", 0);
  }

  emxInit_real_T(&c_st, &Mag_fim, 2, &nb_emlrtRTEI, true);
  i = Mag_fim->size[0] * Mag_fim->size[1];
  Mag_fim->size[0] = csz_idx_0;
  Mag_fim->size[1] = csz_idx_1;
  emxEnsureCapacity_real_T(&b_st, Mag_fim, i, &jb_emlrtRTEI);
  c_st.site = &fe_emlrtRSI;
  csz_idx_1 *= csz_idx_0;
  d_st.site = &xg_emlrtRSI;
  overflow = ((1 <= csz_idx_1) && (csz_idx_1 > 2147483646));
  if (overflow) {
    e_st.site = &eb_emlrtRSI;
    check_forloop_overflow_error(&e_st);
  }

  for (csz_idx_0 = 0; csz_idx_0 < csz_idx_1; csz_idx_0++) {
    Mag_fim->data[csz_idx_0] = muDoubleScalarHypot(b2->data[csz_idx_0],
      imy->data[csz_idx_0]);
  }

  emxInit_real_T(&c_st, &Det_hessian, 2, &ob_emlrtRTEI, true);
  emxInit_real_T(&c_st, &L_imxy, 2, &qb_emlrtRTEI, true);

  /* Step3b: Laplacian */
  st.site = &f_emlrtRSI;
  gaussgradient(&st, b2, params->sigmagradient[0], Det_hessian, L_imxy);
  st.site = &g_emlrtRSI;
  gaussgradient(&st, imy, params->sigmagradient[0], domain, b1);
  emlrtSizeEqCheckNDR2012b(*(int32_T (*)[2])Det_hessian->size, *(int32_T (*)[2])
    b1->size, &emlrtECI, sp);
  i = b2->size[0] * b2->size[1];
  b2->size[0] = Det_hessian->size[0];
  b2->size[1] = Det_hessian->size[1];
  emxEnsureCapacity_real_T(sp, b2, i, &kb_emlrtRTEI);
  csz_idx_0 = Det_hessian->size[0] * Det_hessian->size[1];
  emxFree_real_T(&imy);
  for (i = 0; i < csz_idx_0; i++) {
    b2->data[i] = Det_hessian->data[i] + b1->data[i];
  }

  emxInit_boolean_T(sp, &r, 2, &qb_emlrtRTEI, true);

  /* Laplacian (trace of the 2D matrix) */
  i = r->size[0] * r->size[1];
  r->size[0] = b2->size[0];
  r->size[1] = b2->size[1];
  emxEnsureCapacity_boolean_T(sp, r, i, &lb_emlrtRTEI);
  csz_idx_0 = b2->size[0] * b2->size[1];
  for (i = 0; i < csz_idx_0; i++) {
    r->data[i] = (b2->data[i] > 0.0);
  }

  emlrtSizeEqCheckNDR2012b(*(int32_T (*)[2])b2->size, *(int32_T (*)[2])r->size,
    &b_emlrtECI, sp);
  csz_idx_0 = b2->size[0] * b2->size[1];
  for (i = 0; i < csz_idx_0; i++) {
    b2->data[i] *= static_cast<real_T>(r->data[i]);
  }

  /* Step3c: Hessian Determniant */
  emlrtSizeEqCheckNDR2012b(*(int32_T (*)[2])Det_hessian->size, *(int32_T (*)[2])
    b1->size, &c_emlrtECI, sp);
  emlrtSizeEqCheckNDR2012b(*(int32_T (*)[2])L_imxy->size, *(int32_T (*)[2])
    domain->size, &d_emlrtECI, sp);
  csz_idx_0 = Det_hessian->size[0] * Det_hessian->size[1];
  for (i = 0; i < csz_idx_0; i++) {
    Det_hessian->data[i] *= b1->data[i];
  }

  emxFree_real_T(&b1);
  csz_idx_0 = L_imxy->size[0] * L_imxy->size[1];
  for (i = 0; i < csz_idx_0; i++) {
    L_imxy->data[i] *= domain->data[i];
  }

  emlrtSizeEqCheckNDR2012b(*(int32_T (*)[2])Det_hessian->size, *(int32_T (*)[2])
    L_imxy->size, &c_emlrtECI, sp);
  csz_idx_0 = Det_hessian->size[0] * Det_hessian->size[1];
  for (i = 0; i < csz_idx_0; i++) {
    Det_hessian->data[i] -= L_imxy->data[i];
  }

  emxFree_real_T(&L_imxy);
  i = r->size[0] * r->size[1];
  r->size[0] = Det_hessian->size[0];
  r->size[1] = Det_hessian->size[1];
  emxEnsureCapacity_boolean_T(sp, r, i, &mb_emlrtRTEI);
  csz_idx_0 = Det_hessian->size[0] * Det_hessian->size[1];
  for (i = 0; i < csz_idx_0; i++) {
    r->data[i] = (Det_hessian->data[i] < 0.0);
  }

  emlrtSizeEqCheckNDR2012b(*(int32_T (*)[2])Det_hessian->size, *(int32_T (*)[2])
    r->size, &e_emlrtECI, sp);
  csz_idx_0 = Det_hessian->size[0] * Det_hessian->size[1];
  for (i = 0; i < csz_idx_0; i++) {
    Det_hessian->data[i] *= static_cast<real_T>(r->data[i]);
  }

  emxFree_boolean_T(&r);

  /* Step4: Masking function using tanh on the Summed Derivatives from Step3 */
  csz_idx_0 = Mag_fim->size[0] * Mag_fim->size[1];
  for (i = 0; i < csz_idx_0; i++) {
    Mag_fim->data[i] *= params->alpha;
  }

  csz_idx_0 = b2->size[0] * b2->size[1];
  for (i = 0; i < csz_idx_0; i++) {
    b2->data[i] *= params->beta;
  }

  emlrtSizeEqCheckNDR2012b(*(int32_T (*)[2])Mag_fim->size, *(int32_T (*)[2])
    b2->size, &f_emlrtECI, sp);
  st.site = &h_emlrtRSI;
  c_abs(&st, Det_hessian, domain);
  csz_idx_0 = domain->size[0] * domain->size[1];
  emxFree_real_T(&Det_hessian);
  for (i = 0; i < csz_idx_0; i++) {
    domain->data[i] *= params->epsilon;
  }

  emlrtSizeEqCheckNDR2012b(*(int32_T (*)[2])Mag_fim->size, *(int32_T (*)[2])
    domain->size, &f_emlrtECI, sp);
  csz_idx_0 = Mag_fim->size[0] * Mag_fim->size[1];
  for (i = 0; i < csz_idx_0; i++) {
    Mag_fim->data[i] = params->gamma - ((Mag_fim->data[i] + b2->data[i]) +
      domain->data[i]) / params->delta;
  }

  emxFree_real_T(&b2);
  emxFree_real_T(&domain);
  st.site = &i_emlrtRSI;
  b_st.site = &yg_emlrtRSI;
  csz_idx_1 = Mag_fim->size[0] * Mag_fim->size[1];
  c_st.site = &ie_emlrtRSI;
  overflow = ((1 <= csz_idx_1) && (csz_idx_1 > 2147483646));
  if (overflow) {
    d_st.site = &eb_emlrtRSI;
    check_forloop_overflow_error(&d_st);
  }

  for (csz_idx_0 = 0; csz_idx_0 < csz_idx_1; csz_idx_0++) {
    Mag_fim->data[csz_idx_0] = muDoubleScalarTanh(Mag_fim->data[csz_idx_0]);
  }

  csz_idx_0 = Mag_fim->size[0] * Mag_fim->size[1];
  for (i = 0; i < csz_idx_0; i++) {
    Mag_fim->data[i] = 0.5 * (Mag_fim->data[i] + 1.0);
  }

  /*  masking function */
  emlrtSizeEqCheckNDR2012b(*(int32_T (*)[2])J->size, *(int32_T (*)[2])
    Mag_fim->size, &g_emlrtECI, sp);
  csz_idx_0 = J->size[0] * J->size[1];
  for (i = 0; i < csz_idx_0; i++) {
    J->data[i] *= Mag_fim->data[i];
  }

  emxFree_real_T(&Mag_fim);

  /*  masked image applied on smoothened image */
  emlrtHeapReferenceStackLeaveFcnR2012b(sp);
}

/* End of code generation (diff_gradients.cpp) */
