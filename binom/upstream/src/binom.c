#include <Rmath.h>
#include <R.h>

#define NCP 0
#define LOWER_TAIL 1
#define UPPER_TAIL 0
#define LOG_P 0

double dbeta_shift(double x, double *p) {
  double y = p[0];
  double a = p[1];
  double b = p[2];
  return dbeta(x, a, b, NCP) - y;
}

double zeroin(double (*f)(double x, double *params), /* function: f(x, params)       */
              double ax,                             /* lower bound                  */
              double bx,                             /* upper bound                  */
              double *params,                        /* additional parameters for f  */
              double tol,                            /* tolerance on convergence     */
              int maxit) {                           /* maximum number of iterations */
  double a, b, c, fa, fb, fc;
  a = ax;
  b = bx;
  fa = (*f)(a, params);
  fb = (*f)(b, params);
  c = a;
  fc = fa;
  maxit++;
  while(maxit--) {
    double prev_step = b - a;
    double tol_act;
    double p;
    double q;
    double new_step;
    if(fabs(fc) < fabs(fb)) {
      a = b;
      b = c;
      c = a;
      fa = fb;
      fb = fc;
      fc = fa;
    }
    tol_act = 2 * DBL_EPSILON * fabs(b) + tol/2;
    new_step = (c - b)/2;
    if(fabs(new_step) <= tol_act || fb == (double)0)
      return b;
    if(fabs(prev_step) >= tol_act && fabs(fa) > fabs(fb) ) {
      register double t1, cb, t2;
      cb = c - b;
      if( a==c ) {
        t1 = fb/fa;
        p = cb * t1;
        q = 1.0 - t1;
      } else {
        q = fa/fc;
        t1 = fb/fc;
        t2 = fb/fa;
        p = t2 * ( cb * q * (q - t1) - (b - a) * (t1 - 1.0) );
        q = (q - 1.0) * (t1 - 1.0) * (t2 - 1.0);
      }
      if(p > 0) q = -q;
      else      p = -p;
      if(p < (0.75 * cb * q - fabs(tol_act * q)/2) && p < fabs(prev_step * q/2))
        new_step = p/q;
    }
    if(fabs(new_step) < tol_act)
      new_step = (new_step > 0) ? tol_act : -tol_act;
    a = b;
    fa = fb;
    b += new_step;
    fb = (*f)(b, params);
    if((fb > 0 && fc > 0) || (fb < 0 && fc < 0)) {
      c = a;
      fc = fa;
    }
  }
  /* failed! */
  return b;
}

void binom_bayes(int *x,
                 int *n,
                 double *a,
                 double *b,
                 double *alpha,
                 double *lcl,
                 double* ucl,
                 int *len,
                 int *maxit,
                 double *tol,
                 int *error) {
  int i, j, first, down;
  double lcl_x, ucl_x, lcl_y, ucl_y;
  double y1, y2, y3;
  double px1, px2, sig;
  double mode, xx;
  double x1, x2;
  double lx1, lx2, ux1, ux2;
  double p[3];
  for(j = 0; j < len[0]; j++) {
    lcl_x = lcl[j];
    ucl_x = ucl[j];
    lcl_y = dbeta(lcl_x, a[j], b[j], NCP);
    ucl_y = dbeta(ucl_x, a[j], b[j], NCP);
    y3 = fmax(lcl_y, ucl_y);
    y1 = 0;
    mode = (a[j] - 1)/(a[j] + b[j] - 2);
    first = (lcl_y < ucl_y ? 0 : 1);
    x1 = first ? mode : 0;
    x2 = first ? 1 : mode;
    p[0] = y3;
    p[1] = a[j];
    p[2] = b[j];
    xx = zeroin(dbeta_shift, x1, x2, p, tol[0], maxit[0]);
    if(first) {
      ucl_x = xx;
    } else {
      lcl_x = xx;
    }
    px1 = pbeta(lcl_x, a[j], b[j], LOWER_TAIL, LOG_P);
    px2 = pbeta(ucl_x, a[j], b[j], UPPER_TAIL, LOG_P);
    sig = px1 + px2;
    down = 0;
    i = 0;
    while(fabs(sig - 2 * alpha[j]) > tol[0] && i < maxit[0]) {
      y2 = (y1 + y3) * 0.5;
      if(down) {
        if(dbeta(lcl_x, a[j], b[j], 0) < y2)
	  lcl_x = mode;
        lx1 = 0;
	lx2 = lcl_x;
        if(dbeta(ucl_x, a[j], b[j], 0) < y2)
	  ucl_x = mode;
        ux1 = ucl_x;
	ux2 = 1;
      } else {
        if(dbeta(lcl_x, a[j], b[j], 0) > y2)
	  lcl_x = 0;
        lx1 = lcl_x;
	lx2 = mode;
        if(dbeta(ucl_x, a[j], b[j], 0) > y2)
	  ucl_x = 1;
        ux1 = mode;
	ux2 = ucl_x;
      }
      p[0] = y2;
      lcl_x = zeroin(dbeta_shift, lx1, lx2, p, tol[0], maxit[0]);
      ucl_x = zeroin(dbeta_shift, ux1, ux2, p, tol[0], maxit[0]);
      px1 = pbeta(lcl_x, a[j], b[j], LOWER_TAIL, LOG_P);
      px2 = pbeta(ucl_x, a[j], b[j], UPPER_TAIL, LOG_P);
      sig = px1 + px2;
      if(sig > 2 * alpha[j]) {
        down = 0;
        y3 = y2;
      } else {
        down = 1;
        y1 = y2;
      }
      i++;
    }
    error[j] = (i >= maxit[0] ? 1 : 0);
    lcl[j] = lcl_x;
    ucl[j] = ucl_x;
  }
}
