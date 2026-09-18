#include <string.h>
#include <stdbool.h>
#include <math.h>
#include "BlasGateway.h"
#include "error.h"
#include "VarmaUtilities.h"
#include "RandompackGateway.h"
#include "varmapack.h"

static void flipud(int m, int n, double src[], double dst[]);
static void hilb(double A[], int m);
static int find_named_case(const char *namev[], const char *name, int Ncase);
static void construct_rho_A(double A[], int p, int r);
static bool scale_to_rho(double A[], int p, int r, double rho);
static double scaled_specrad(double A[], double scale, int p, int r);

// Construct or inquire about built-in VARMA testcases.
bool varmapack_testcase( char *name, int *index, double rho, int *pp, int *qp, int *rp,
  double A[], double B[], double Sig[], randompack_rng *rng) {
  // With A, B, and Sig all zero, inquire by named case or positive index.
  // The name "max" returns the named-case count and maximum p, q, r.
  // The name "rho" is construction-only. Otherwise construct a named/indexed case,
  // an unnamed deterministic case (index -1), an unnamed random case (index 0),
  // or a case with target spectral radius (name "rho").
#define p15 5
#define r15 7
#define seed15 42
#define c15 0.05
#define n15 (r15*r15*p15)
#define p16 3
#define q16 3
#define r16 7
#define rho16 .98
  const char *namev[] = {// nr   p  q  r
    "tinyAR",      // 1    1  0  1
    "tinyMA",      // 2    0  1  1
    "tinyARMA",    // 3    1  1  1
    "smallAR1",    // 4    1  0  2
    "smallAR2",    // 5    2  0  2
    "smallMA1",    // 6    0  1  2
    "smallMA2",    // 7    0  2  2
    "smallARMA1",  // 8    1  1  2
    "smallARMA2",  // 9    1  2  2
    "mediumAR",    // 10   1  0  3
    "mediumMA1",   // 11   0  1  3
    "mediumARMA1", // 12   3  3  3
    "mediumARMA2", // 13   3  3  3
    "mediumMA2",   // 14   0  2  3
    "largeAR",     // 15   5  0  7
    "largeARMA"    // 16   3  3  7
  };
  int    pv[]      = { 1,  0,  1,  1,  2,  0,  0,  1,  1,  1,  0,  3,  3,  0, p15, p16};
  int    qv[]      = { 0,  1,  1,  0,  0,  1,  2,  1,  2,  0,  1,  3,  3,  2, 0, q16};
  int    rv[]      = { 1,  1,  1,  2,  2,  2,  2,  2,  2,  3,  3,  3,  3,  3, r15, r16};
  int *icase = index;
  int p, q, r, Ncase = sizeof(namev)/sizeof(char*);
  clear_error();
  // SANITY CHECKS:
  if (icase == 0) return fail_error("invalid argument");
  int nnull = !A + !B + !Sig;
  if (nnull != 0 && nnull !=3)
    return fail_error("invalid argument");
  bool MAX = name && !strcmp(name, "max");
  bool RHO = name && !strcmp(name, "rho");
  bool INQUIRY = nnull == 3;
  bool NAMED = name && strlen(name) > 0 && !RHO;
  if (INQUIRY && RHO) return fail_error("invalid argument");
  if (!INQUIRY) {
    if (NAMED) {
      *icase = find_named_case(namev, name, Ncase);
      if (*icase == 0) return fail_error("unknown testcase");
    }
    if (*icase <= 0 || RHO) {
      if (!pp || !qp || !rp) return fail_error("invalid argument");
      if (*pp < 0 || *qp < 0 || *rp <= 0) {
        return fail_error("invalid argument");
      }
      if (*icase == 0 && !RHO && rng == 0) return fail_error("invalid argument");
    }
    if (!RHO && (*icase < -1 || *icase > Ncase)) return fail_error("unknown testcase");
  }
  else if (!MAX && NAMED) {
    *icase = find_named_case(namev, name, Ncase);
    if (*icase == 0) return fail_error("unknown testcase");
  }
  else if (INQUIRY && !NAMED) {
    if (*icase == 0 || *icase == -1) {
      return fail_error("invalid argument");
    }
    if (*icase < -1 || *icase > Ncase) {
      return fail_error("unknown testcase");
    }
  }
  bool LARGE_ARMA = !RHO && *icase == 16;
  // MAX INQUIRY
  if (MAX) {
    int pmax = 0, qmax = 0, rmax = 0;
    *icase = Ncase;
    for (int k=0; k<Ncase; k++) {
      pmax = imax(pmax, pv[k]);
      qmax = imax(qmax, qv[k]);
      rmax = imax(rmax, rv[k]);
    }
    if (pp) *pp = pmax;
    if (qp) *qp = qmax;
    if (rp) *rp = rmax;
    return true;
  }
  // CASE INQUIRY
  else if (INQUIRY) {
    if (NAMED) {
      *icase = find_named_case(namev, name, Ncase);
    }
    else if (name) {
      snprintf(name, VARMAPACK_TESTCASE_NAME_LEN, "%s", namev[*icase - 1]);
    }
    int k = *icase - 1;
    if (pp) *pp = pv[k];
    if (qp) *qp = qv[k];
    if (rp) *rp = rv[k];
    return true;
  }
  // CONSTRUCT TESTCASE
  double
    A1[] = {0.5}, A3[] = {0.4}, A4[] = {0.1, 0.1,   0.1, 0.100000001},
    A5[] = {0.60, 0.05, 0.30, 0.02,   0.04, 0.60, 0.02, 0.30},
    A8[] = {0.3, 0.1,   0.4, 0.2}, A9[] = {0.2, 0.3,   0.2, 0.3}, A10[] = {
      0.35, 0.15, 0.15,//  0.11, 0.14, 0.07,
      0.25, 0.15, 0.05,//  0.12, 0.15, 0.08,
      0.15, 0.05, 0.01,//  0.13, 0.16, 0.09
    }, A12[9*3], A13[9*3], A33[] = {
      0.15, 0.10, 0.05,  0.11, 0.14, 0.17,  0.01, 0.04, 0.06,
      0.16, 0.11, 0.06,  0.12, 0.15, 0.18,  0.02, 0.05, 0.08,
      0.17, 0.12, 0.07,  0.13, 0.16, 0.19,  0.03, 0.06, 0.09
    }, B2[] = {0.5}, B3[] = {0.4},
    B6[] = {0.3, 0.1,   0.1, 0.2}, B7[] = {0.3, 0.1, 0.1, 0.1,   0.3, 0.1, 0.1, 0.1},
    B8[] = {0.2, 0.3,   0.2, 0.3}, B9[] = {0.4, 0.1, 0.2, 0.1,   0.2, 0.1, 0.3, 0.1},
    B11[] = {
      0.35, 0.25, 0.15, 0.25, 0.15, 0.05, 0.15, 0.05, 0.01
    }, B12[9*3], B13[9*3], B14[] =  {
      0.35, 0.25, 0.15,  0.11, 0.14, 0.17, 0.25, 0.15, 0.05,  0.12, 0.15, 0.18,
      0.15, 0.05, 0.01,  0.13, 0.16, 0.19
    }, S1[] = {0.8}, S2[] = {2, 1,   1, 3}, S3[] = {2, 1,   1, 2},
    S4[] = {2.0, 0.5, 0.0,   0.5, 2.0, 0.5,   0.0, 0.5, 1.0};
  copy(9*3, A33, 1, A12, 1);
  flipud(9, 3, A33, B12);
  flipud(9, 3, A33, A13);
  for (int j=0; j<3; j++)
    for (int i=3; i<9; i++) A13[i + j*9] += 0.01;
  copy(9*3, A33, 1, B13, 1);
  double *Av[] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  double *Bv[] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  double *Sigv[] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  int i, j;
  // Define the named cases with fixed coefficients.
  Av[0]=A1; Av[2]=A3; Av[3]=A4; Av[4]=A5; Av[7]=A8; Av[8]=A9;
  Av[9]=A10; Av[11]=A12; Av[12]=A13;
  Bv[1]=B2; Bv[2]=B3; Bv[5]=B6; Bv[6]=B7; Bv[7]=B8; Bv[8]=B9;
  Bv[10]=B11; Bv[11]=B12; Bv[12]=B13; Bv[13]=B14;
  Sigv[0]=Sigv[1]=Sigv[2]=S1; Sigv[3]=Sigv[4]=S2;
  Sigv[5]=Sigv[6]=Sigv[7]=Sigv[8]=S3;
  Sigv[9]=Sigv[10]=Sigv[11]=Sigv[12]=Sigv[13]=S4;
  if (RHO || *icase <= 0) {
    p = *pp;
    q = *qp;
    r = *rp;
  }
  else {
    p = pv[*icase-1];
    q = qv[*icase-1];
    r = rv[*icase-1];
  }
  if (RHO || LARGE_ARMA || *icase == -1 || *icase == 0) {
    if (Sig) {
      hilb(Sig, r);
      for (i = 0; i < r; i++) Sig[i + i*r] += 0.2; // add 0.2 to diagonal
    }
  }
  if (RHO || LARGE_ARMA) {
    if (A && p > 0) {
      construct_rho_A(A, p, r);
      if (!scale_to_rho(A, p, r, RHO ? rho : rho16)) return false;
    }
    if (A && p == 0 && rho != 0) return fail_error("invalid argument");
    if (B && q>0) for (i=0; i<r*r*q; i++) B[i] = 1.0/(q*r);
    if (LARGE_ARMA && name && !NAMED) {
      snprintf(name, VARMAPACK_TESTCASE_NAME_LEN, "%s", namev[*icase - 1]);
    }
  }
  else if (*icase == -1) {
    if (A && p>0) for (i=0; i<r*r*p; i++) A[i] = 0.5/(p*r);
    if (B && q>0) for (i=0; i<r*r*q; i++) B[i] = 1.0/(q*r);
  }
  else if (*icase == 0) {
    if (A && p>0) {
      if (!randompack_u01(A, r*r*p, rng)) return fail_error(randompack_last_error(rng));
      scal(r*r*p, 0.5/(p*r), A, 1);
      j = 0;
      double specrad = varmapack_specrad(A, r, p);
      if (isnan(specrad)) return false;
      while (specrad >= 1) {
        scal(r*r*p, 0.5, A, 1);
        specrad = varmapack_specrad(A, r, p);
        if (isnan(specrad)) return false;
        j++;
        if (j >= 10) return fail_error("internal error");
      }
    }
    if (B && q>0) {
      if (!randompack_u01(B, r*r*q, rng)) return fail_error(randompack_last_error(rng));
      scal(r*r*q, 1.0/(q*r), B, 1);
    }
  }
  else if (1 <= *icase && *icase <= Ncase) {
    if (*icase == 15) {
      double A15[n15], S6[r15*r15];
      randompack_rng *rng15 = randompack_create(0);
      if (rng15 == 0) return fail_error("allocation failed");
      if (!randompack_seed(seed15, 0, 0, rng15)) {
        fail_error(randompack_last_error(rng15));
        randompack_free(rng15);
        return false;
      }
      if (!randompack_u01(A15, n15, rng15)) {
        fail_error(randompack_last_error(rng15));
        randompack_free(rng15);
        return false;
      }
      randompack_free(rng15);
      scal(n15, c15, A15, 1);
      hilb(S6, r15);
      for (i=0; i<r15; i++) S6[i + i*r] += 1;
      lacpy("All", r, p*r, A15, r, A, r);
      if (Sig) copy(r*r, S6, 1, Sig, 1);
    }
    else {
      if (A && Av[*icase-1]) copytranspose(p*r, r, Av[*icase-1], p*r, A, r);
      if (B && Bv[*icase-1]) copytranspose(q*r, r, Bv[*icase-1], q*r, B, r);
      if (Sig && Sigv[*icase-1]) copy(r*r, Sigv[*icase-1], 1, Sig, 1);
    }
    if (name && !NAMED) {
      snprintf(name, VARMAPACK_TESTCASE_NAME_LEN, "%s", namev[*icase-1]);
    }
  }
  else return fail_error("internal error");
  if (pp) *pp = p;
  if (qp) *qp = q;
  if (rp) *rp = r;
  return true;
}

static void flipud(int m, int n, double *src, double *dst) {
  int i, j;
  for (j=0; j<n; j++)
    for (i=0; i<m; i++)
      dst[i + j*m] = src[m-1-i + j*m];
}

static void hilb(double A[], int n) {
  // Return n by n Hilbert matrix
  int i, j;
  for (j=0; j<n; j++)
    for (i=0; i<n; i++)
      A[i + j*n] = 1.0/(i + j + 1);
}

static int find_named_case(const char *namev[], const char *name, int Ncase) {
  for (int i = 0; i < Ncase; i++) {
    if (strcmp(namev[i], name) == 0) return i + 1;
  }
  return 0; // not found
}

static void construct_rho_A(double A[], int p, int r) {
  double scale = 2/(r*pow(p, 1.0/3.0));
  for (int j=0; j<p*r; j++) {
    int jj = j + 1;
    for (int i=0; i<r; i++) {
      int ii = i + 1;
      int mn = ii < jj ? ii : jj;
      int mx = ii > jj ? ii : jj;
      A[i + j*r] = scale*mn/mx;
    }
  }
}

static bool scale_to_rho(double A[], int p, int r, double rho) {
  double lo = 0, hi = 1, mid = 0, midrho = 0, hirho;
  if (rho < 0 || !isfinite(rho)) return fail_error("invalid argument");
  if (p == 0) {
    if (rho != 0) return fail_error("invalid argument");
    return true;
  }
  if (rho == 0) {
    setzero(r*r*p, A);
    return true;
  }
  hirho = scaled_specrad(A, hi, p, r);
  if (isnan(hirho)) return false;
  for (int iter=0; hirho <= rho; iter++) {
    hi = 2*hi;
    if (iter >= 60) return fail_error("internal error");
    hirho = scaled_specrad(A, hi, p, r);
    if (isnan(hirho)) return false;
  }
  for (int iter=0; iter<100; iter++) {
    mid = (lo + hi)/2;
    midrho = scaled_specrad(A, mid, p, r);
    if (isnan(midrho)) return false;
    if (fabs(midrho - rho) < 1e-6) {
      scal(r*r*p, mid, A, 1);
      return true;
    }
    if (midrho < rho) lo = mid;
    else hi = mid;
  }
  scal(r*r*p, (lo + hi)/2, A, 1);
  return true;
}

static double scaled_specrad(double A[], double scale, int p, int r) {
  double *As = 0, rho;
  if (!ALLOC(As, r*r*p)) {
    fail_error("allocation failed");
    return NAN;
  }
  copy(r*r*p, A, 1, As, 1);
  scal(r*r*p, scale, As, 1);
  rho = varmapack_specrad(As, r, p);
  FREE(As);
  return rho;
}
