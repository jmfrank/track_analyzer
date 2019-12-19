/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * diff_gradients_types.h
 *
 * Code generation for function 'diff_gradients_types'
 *
 */

#pragma once

/* Include files */
#include "rtwtypes.h"

/* Type Definitions */
struct emxArray_int32_T
{
  int32_T *data;
  int32_T *size;
  int32_T allocatedSize;
  int32_T numDimensions;
  boolean_T canFreeData;
};

struct struct0_T
{
  real_T seg_channel;
  real_T nuclear_channel;
  real_T AbsMinVol;
  real_T AbsMaxVol;
  real_T alpha;
  real_T beta;
  real_T gamma;
  real_T epsilon;
  real_T delta;
  real_T bg;
  real_T sigmagradient[3];
  real_T noise_filter[100];
  real_T median_filter[3];
  real_T diffuse_iterations;
  real_T kappa1;
  real_T kappa2;
  real_T option;
  real_T thrshlevel;
  real_T percentile;
  real_T I_sm_sigma;
  real_T Hdepth;
  real_T nconn_BW;
  real_T nconn_BW2;
  real_T z_effect;
  real_T h_min_depth;
  real_T h_min_conn;
  real_T MeanFilterSensitivity;
  real_T MeanFilterNeighborhood[3];
  real_T OutlierThreshold;
  real_T StrelXYZ[3];
  real_T imclose_r;
  real_T rg_threshold;
};

struct emxArray_real_T
{
  real_T *data;
  int32_T *size;
  int32_T allocatedSize;
  int32_T numDimensions;
  boolean_T canFreeData;
};

struct emxArray_boolean_T
{
  boolean_T *data;
  int32_T *size;
  int32_T allocatedSize;
  int32_T numDimensions;
  boolean_T canFreeData;
};

struct cell_wrap_14
{
  emxArray_real_T *f1;
};

struct emxArray_int8_T
{
  int8_T *data;
  int32_T *size;
  int32_T allocatedSize;
  int32_T numDimensions;
  boolean_T canFreeData;
};

struct emxArray_uint32_T
{
  uint32_T *data;
  int32_T *size;
  int32_T allocatedSize;
  int32_T numDimensions;
  boolean_T canFreeData;
};

/* End of code generation (diff_gradients_types.h) */
