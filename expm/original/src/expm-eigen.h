/*  ===== File part of R package expm =====
 *
 *  expm-eigen.h
 *
 *  Created by Christophe Dutang on 27/02/08.
 *
 */
#include "expm.h"

SEXP do_expm_eigen(SEXP x, SEXP tolin);
void expm_eigen(double *x, int n, double *z, double tol);


