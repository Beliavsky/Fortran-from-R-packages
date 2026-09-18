#ifndef LYAPUNOV_H
#define LYAPUNOV_H

#include <stdbool.h>
#include "varmapack_config.h"

HIDDEN bool LyapunovFactorizeSolve( double A[], double B[], double Sig[],
  int p, int q, int r, double S[], double C[], double G[]);

#endif
