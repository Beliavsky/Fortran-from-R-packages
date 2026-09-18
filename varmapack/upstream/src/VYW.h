// VYW.h -- internal vector Yule-Walker setup utilities
#ifndef VYW_H
#define VYW_H

#include <stdbool.h>
#include "varmapack_config.h"

HIDDEN bool VYWFactorizeSolve(
  double A[],   // in   r×r×p, autoregressive parameter matrices
  double B[],   // in   r×r×q, moving average parameter matrices (ignored when q=0)
  double Sig[], // in   r×r, innovation covariance
  int p,        // in   number of autoregressive terms
  int q,        // in   number of moving-average terms
  int r,        // in   dimension of each x(t)
  double S[],   // out  r×r×(p+1), covariance sequence S0…Sp
  double C[],   // out  r×r×(q+1) buffer for Ci
  double G[]);  // out  r×r×(q+1) buffer for Gi

#endif
