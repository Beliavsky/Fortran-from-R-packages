// VarmaPackUtil.h
// Utility functions for VARMA simulation and covariance computation
// (Part of varmapack / varmasim)

#include <stdio.h>
#include <stdbool.h>
#include "varmapack.h"
#include "varmapack_config.h"

HIDDEN void FindC( // Compute Ci = cov(x(t), eps(t-i))
  double A[],   // in   r×r×p  autoregressive parameter matrices
  double B[],   // in   r×r×q  moving average parameter matrices
  double Sig[], // in   r×r    covariance of the shock terms eps(t)
  int p,        // in   number of autoregressive terms
  int q,        // in   number of moving average terms
  int r,        // in   dimension of each x(t)
  double C[]    // out  r×r×(q+1)  Ci = cov(x(t), eps(t−i))
);

HIDDEN void FindG( // Compute Gi = cov(y(t), x(t-i))
  double B[],   // in   r×r×q  moving average parameter matrices
  double C[],   // in   r×r×(q+1)  Ci = cov(x(t), eps(t−i))
  int q,        // in   number of moving average terms
  int r,        // in   dimension of each x(t)
  double G[]    // out  r×r×(q+1)  Gi = cov(y(t), x(t−i))
);

HIDDEN bool FindS( // Compute S0,...,Sp and the corresponding C and G blocks
  double A[], double B[], double Sig[], int p, int q, int r,
  double S[], double C[], double G[]);

HIDDEN bool SBuild(char *uplo, double S[], double A[], double G[], int p, int q,
  int r, int n, double SS[]);

HIDDEN void CCBuild(double A[], double C[], int p, int q, int r, int n, double CC[]);

HIDDEN bool FindW(
  double B[],   // in   r×r×q  moving average parameter matrices
  double Sig[], // in   r×r    covariance of the shock terms eps(t)
  int q,        // in   number of moving average terms
  int r,        // in   dimension of each x(t)
  double W[]    // out  r×r×(q+1)  Wi = cov(y(t), y(t−i))
);

HIDDEN void WBuild(double Wlag[], int q, int r, int nW, double W[]);

HIDDEN void postmultiplySigmaPrime(double Y[], int ldY, int m, int h, double Sig[],
                                   int r, double Work[]);

HIDDEN bool drawConditionalStartupShocks(
  double B[], double Sig[], double residual[], int q, int r, int h, int t0,
  int firstShockTime, int firstActiveShock, int M, double E[], int ldE,
  randompack_rng *rng);

HIDDEN bool requirePositiveDefinite(double Sig[], int r);

// Build the lower block-triangular map from startup shocks to startup states.
HIDDEN void FindPsi(double *A, double *B, double *Psi, int p, int q, int r);

// Multiply each shock block of Psi by a positive-semidefinite factor of Sig.
HIDDEN bool FindPsiHat(double *Psi, double *Psi_hat, double *Sig, int r, int h);

// Factor positive-semidefinite Sig as L*L^T and report whether L is triangular.
HIDDEN bool psdFactor(double Sig[], int r, double L[], bool *triangular);

// Return the spectral radius of a companion matrix with first block row sign*P.
HIDDEN double CompanionSpecrad(double P[], int r, int k, double sign);
