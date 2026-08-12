/*  ===== File part of R package expm =====
 *
 *  logm-eigen.h
 *
 *  Created by Christophe Dutang on 13/05/08.
 *
 */
#include "expm.h"

SEXP do_logm_eigen(SEXP x, SEXP tolin);
void logm_eigen(double *x, int n, double *z, double tol);


