#include "Lyapunov.h"

#include "BlasGateway.h"
#include "sb03mdGateway.h"
#include "VarmaPackUtil.h"
#include "VarmaUtilities.h"
#include "error.h"

static double *block(double X[], int ld, int r, int i, int j) {
  return X + i*r + ld*j*r;
}

static void add_identity_block(double X[], int ld, int r, int i, int j) {
  double *Xij = block(X, ld, r, i, j);
  for (int k=0; k<r; k++) Xij[k + ld*k] += 1;
}

static void build_state_matrices(double A[], double B[], double Sig[], int p, int q,
                                 int r, int nx, int ne, int nstate,
                                 double F[], double Q[]) {
  setzero(nstate*nstate, F);
  setzero(nstate*nstate, Q);
  for (int i=0; i<p; i++) {
    lacpy("All", r, r, A + i*r*r, r, block(F, nstate, r, 0, i), nstate);
  }
  for (int i=0; i<q; i++) {
    lacpy("All", r, r, B + i*r*r, r, block(F, nstate, r, 0, nx+i), nstate);
  }
  for (int i=1; i<nx; i++) add_identity_block(F, nstate, r, i, i-1);
  for (int i=1; i<ne; i++) add_identity_block(F, nstate, r, nx+i, nx+i-1);
  lacpy("All", r, r, Sig, r, block(Q, nstate, r, 0, 0), nstate);
  lacpy("All", r, r, Sig, r, block(Q, nstate, r, 0, nx), nstate);
  lacpy("All", r, r, Sig, r, block(Q, nstate, r, nx, 0), nstate);
  lacpy("All", r, r, Sig, r, block(Q, nstate, r, nx, nx), nstate);
}

HIDDEN bool LyapunovFactorizeSolve(double A[], double B[], double Sig[],
                                   int p, int q, int r,
                                   double S[], double C[], double G[]) {
  xAssert((A != 0 || p == 0) && Sig != 0 && S != 0);
  xAssert(B != 0 || q == 0);
  xAssert(p >= 0 && q >= 0 && r > 0);
  xAssert(C != 0 && G != 0);
  int nx = p > 0 ? p + 1 : 1;
  int ne = q + 1;
  int nb = nx + ne;
  int nstate = r*nb;
  int rr = r*r;
  double *F = 0;
  double *Q = 0;
  double *P = 0;
  bool ok = false;
  if (!ALLOC(F, nstate*nstate)) goto done;
  if (!ALLOC(Q, nstate*nstate)) goto done;
  if (!ALLOC(P, nstate*nstate)) goto done;
  build_state_matrices(A, B, Sig, p, q, r, nx, ne, nstate, F, Q);
  if (!sb03mdGateway(F, Q, nstate, P)) goto done;
  for (int i=0; i<=p; i++) {
    lacpy("All", r, r, block(P, nstate, r, 0, i), nstate, S + i*rr, r);
  }
  for (int i=0; i<=q; i++) {
    lacpy("All", r, r, block(P, nstate, r, 0, nx+i), nstate, C + i*rr, r);
  }
  FindG(B, C, q, r, G);
  ok = true;
done:
  FREE(P);
  FREE(Q);
  FREE(F);
  return ok;
}
