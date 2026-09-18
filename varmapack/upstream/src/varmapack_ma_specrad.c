// varmapack_ma_specrad  Spectral radius of companion matrix of a VMA process
#include "error.h"
#include "VarmaPackUtil.h"

double varmapack_ma_specrad(double B[], int r, int q) {
  clear_error();
  return CompanionSpecrad(B, r, q, -1);
}
