/*
 *  AUTHOR
 *    Catherine Loader, catherine@research.bell-labs.com.
 *    October 23, 2000.
 *
 *  Merge in to R:
 *	Copyright (C) 2000-2024, The R Core Team
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, a copy is available at
 *  https://www.R-project.org/Licenses/
 *
 *
 *  DESCRIPTION
 *
 *    Computes the log of the error term in Stirling's formula.
 *      For n > 15, uses the series 1/12n - 1/360n^3 + ...
 *      For n <=15, integers or half-integers, uses stored values.
 *      For other n < 15, uses lgamma directly (don't use this to
 *        write lgamma!)
 *
 * Merge in to R:
 * Copyright (C) 2000, The R Core Team
 * R has lgammafn, and lgamma is not part of ISO C
 */

#include "DPQpkg.h"

/* stirlerr(n) = log(n!) - log( sqrt(2*pi*n)*(n/e)^n )
 *             = log Gamma(n+1) - 1/2 * [log(2*pi) + log(n)] - n*[log(n) - 1]
 *             = log Gamma(n+1) - (n + 1/2) * log(n) + n - log(2*pi)/2
 *
 * see also lgammacor() in ./lgammacor.c  which computes almost the same!
 *
 * NB: stirlerr(n/2) & stirlerr((n+1)/2) are called from dt(x,n) for all real n > 0 ;
 *     stirlerr(x) from gammafn(x) when |x| > 10,  2|x| is integer, but |x| is *not* in {11:50}
 *     stirlerr(x)   from dpois_raw(x, lam) for any  x > 0; --> ./dpois.c [dpois_raw() called by many]
 *     stirlerr(n), stirlerr(x), stirlerr(n-x) from binom_raw(x, n, ..) for all possible 0 < x < n
 */

/* NB 2 --- in DPQ's R code, we have much improved and versatile stirlerr() ==> ../R/dgamma.R
 * --       This is *really* only for reproducibility of R's C level (non-API!) stirlerr()
 * _________ Now callable from R via R_dpq_stirlerr() below
 * also only used in gammafn_ver() ==> ./gamma-variants.c

_________ TODO:  add  'int version'  or 'order' or   double* cutoffs ________________
 */
double dpq_stirlerr(double n, stirlerr_version_t version)
{

#define S0 0.083333333333333333333       /* 1/12 */
#define S1 0.00277777777777777777778     /* 1/360 */
#define S2 0.00079365079365079365079365  /* 1/1260 */
#define S3 0.000595238095238095238095238 /* 1/1680 */
#define S4 0.0008417508417508417508417508/* 1/1188 */
/* From S5 on: ---> ../Misc/stirlerr-trms.R  ==> gmp::BernoulliQ() etc */
#define  S5   0.0019175269175269175269175262 // 691/360360
#define  S6   0.0064102564102564102564102561 // 1/156
#define  S7   0.029550653594771241830065352  // 3617/122400
#define  S8   0.17964437236883057316493850   // 43867/244188
#define  S9   1.3924322169059011164274315    // 174611/125400
#define S10  13.402864044168391994478957     // 77683/5796
#define S11 156.84828462600201730636509     // 236364091/1506960
#define S12 2193.1033333333333333333333     // 657931/300
#define S13 36108.771253724989357173269     // 3392780147/93960
#define S14 691472.26885131306710839498     // 1723168255201/2492028
#define S15 15238221.539407416192283370     // 7709321041217/505920
#define S16 382900751.39141414141414141     // 151628697551/396
#define S17 10882266035.784391089015145     // 26315271553053477373/2418179400
#define S18 347320283765.00225225225225     // 154210205991661/444
#define S19 12369602142269.274454251718     // 261082718496449122051/21106800

/*
  exact values for 0, 0.5, 1.0, 1.5, ..., 14.5, 15.0.
*/
    const static double sferr_halves[31] = {
	0.0, /* n=0 - wrong, place holder only */
	0.1534264097200273452913848,  /* 0.5 */
	0.0810614667953272582196702,  /* 1.0 */
	0.0548141210519176538961390,  /* 1.5 */
	0.0413406959554092940938221,  /* 2.0 */
	0.03316287351993628748511048, /* 2.5 */
	0.02767792568499833914878929, /* 3.0 */
	0.02374616365629749597132920, /* 3.5 */
	0.02079067210376509311152277, /* 4.0 */
	0.01848845053267318523077934, /* 4.5 */
	0.01664469118982119216319487, /* 5.0 */
	0.01513497322191737887351255, /* 5.5 */
	0.01387612882307074799874573, /* 6.0 */
	0.01281046524292022692424986, /* 6.5 */
	0.01189670994589177009505572, /* 7.0 */
	0.01110455975820691732662991, /* 7.5 */
	0.010411265261972096497478567, /* 8.0 */
	0.009799416126158803298389475, /* 8.5 */
	0.009255462182712732917728637, /* 9.0 */
	0.008768700134139385462952823, /* 9.5 */
	0.008330563433362871256469318, /* 10.0 */
	0.007934114564314020547248100, /* 10.5 */
	0.007573675487951840794972024, /* 11.0 */
	0.007244554301320383179543912, /* 11.5 */
	0.006942840107209529865664152, /* 12.0 */
	0.006665247032707682442354394, /* 12.5 */
	0.006408994188004207068439631, /* 13.0 */
	0.006171712263039457647532867, /* 13.5 */
	0.005951370112758847735624416, /* 14.0 */
	0.005746216513010115682023589, /* 14.5 */
	0.005554733551962801371038690  /* 15.0 */
    };

// for version = R_4_4
#define LARGE 1e99
#define nC 19
    const static double cutoffs[nC] = {
	5.25, //  (k=19 skipped, 4.9) but 5.25 == cutoffs[0] is crucially needed
	LARGE, //   k=18 skipped, 5.0,
	5.25, // k=17
	LARGE, //   k=16 skipped, was 5.25
	6.1,  // k=15
	LARGE, //   k=14 skipped, was 5.7
	6.6,
	LARGE, //   k=12 skipped, was 6.5,
	7.3,  // nC-11
	LARGE, //   k=10 skipped was 7.9,
	8.9,  // k = 9
	12.3, // nC-8
	12.8, // (was 13)  nC-7 <--> k=7
	//-- the first cutoff to  "high" vs "low" branch
	23.5, // nC-6;  was 20
	////
	27,  // nC-5;  15 (then   26)
	86,  // nC-4;  35 (then   60)
	205, // nC-3;  80 (then  200)
	6180,// nC-2; 500 (then 3300)
	15.7e6 // 8.37e6 for M1 mac; was 17.4e6 // nC-1;
    };

    double nn;
    int k;

  switch(version) { // stirlerr_version_t  <--> ./DPQpkg.h
  case R_3:   // "R<=3"
  case R_3_1: // "R4..1" : using lgamma1p()
    if (n <= 15.0) {
	nn = n + n;
	if (nn == (int)nn) return sferr_halves[(int)nn];
	return ((version == R_3) ? lgamma(n + 1.) : lgamma1p(n))
	    - (n + 0.5)*log(n) + n - M_LN_SQRT_2PI;
    }
    nn = n*n;
    if (n>500) return((S0-S1/nn)/n);
    if (n> 80) return((S0-(S1-S2/nn)/nn)/n);
    if (n> 35) return((S0-(S1-(S2-S3/nn)/nn)/nn)/n);
    /* 15 < n <= 35 : */
    return((S0-(S1-(S2-(S3-S4/nn)/nn)/nn)/nn)/n);

  case R_4_4: // "R4.4_0"  --- see R code in ../R/dgamma.R  <<<<<<
    k = nC-6;
    if (n <= cutoffs[k]) {// == 23.5 (was 20)
      nn = n + n;
      if (n <= 15. && (nn == (int)nn)) return sferr_halves[(int)nn];
      // else:
      if (n <= cutoffs[0]) {
	  if(n >= 1.) { // "MM2"; slightly more accurate than direct form
	      double l_n = log(n);	      // ldexp(u, -1) == u/2
	      return lgamma(n) + n*(1 - l_n) + ldexp(l_n - M_LN_2PI, -1);
	  }
	  else // n < 1
	      return lgamma1p(n) - (n + 0.5)*log(n) + n - M_LN_SQRT_2PI;
      }
      // else cutoffs[0] < n <= cutoffs[nC-6] == 23.5 :
      nn = n*n;
      if (n > cutoffs[--k]) return (S0-(S1-(S2-(S3-(S4-(S5 -S6/nn)/nn)/nn)/nn)/nn)/nn)/n; // k=7
      if (n > cutoffs[--k]) return (S0-(S1-(S2-(S3-(S4-(S5-(S6 -S7/nn)/nn)/nn)/nn)/nn)/nn)/nn)/n;
      if (n > cutoffs[--k]) return (S0-(S1-(S2-(S3-(S4-(S5-(S6-(S7 -S8/nn)/nn)/nn)/nn)/nn)/nn)/nn)/nn)/n; // k=9
      if (n > cutoffs[--k]) return (S0-(S1-(S2-(S3-(S4-(S5-(S6-(S7-(S8 -S9/nn)/nn)/nn)/nn)/nn)/nn)/nn)/nn)/nn)/n;
      if (n > cutoffs[--k]) return (S0-(S1-(S2-(S3-(S4-(S5-(S6-(S7-(S8-(S9-S10/nn)/nn)/nn)/nn)/nn)/nn)/nn)/nn)/nn)/nn)/n;
      if (n > cutoffs[--k]) return (S0-(S1-(S2-(S3-(S4-(S5-(S6-(S7-(S8-(S9-(S10-S11/nn)/nn)/nn)/nn)/nn)/nn)/nn)/nn)/nn)/nn)/nn)/n;
      if (n > cutoffs[--k]) return (S0-(S1-(S2-(S3-(S4-(S5-(S6-(S7-(S8-(S9-(S10-(S11-S12/nn)/nn)/nn)/nn)/nn)/nn)/nn)/nn)/nn)/nn)/nn)/nn)/n;
      if (n > cutoffs[--k]) return (S0-(S1-(S2-(S3-(S4-(S5-(S6-(S7-(S8-(S9-(S10-(S11-(S12-S13/nn)/nn)/nn)/nn)/nn)/nn)/nn)/nn)/nn)/nn)/nn)/nn)/nn)/n;
      if (n > cutoffs[--k]) return (S0-(S1-(S2-(S3-(S4-(S5-(S6-(S7-(S8-(S9-(S10-(S11-(S12-(S13-S14/nn)/nn)/nn)/nn)/nn)/nn)/nn)/nn)/nn)/nn)/nn)/nn)/nn)/nn)/n;
      if (n > cutoffs[--k]) return (S0-(S1-(S2-(S3-(S4-(S5-(S6-(S7-(S8-(S9-(S10-(S11-(S12-(S13-(S14-S15/nn)/nn)/nn)/nn)/nn)/nn)/nn)/nn)/nn)/nn)/nn)/nn)/nn)/nn)/nn)/n;
      if (n > cutoffs[--k]) return (S0-(S1-(S2-(S3-(S4-(S5-(S6-(S7-(S8-(S9-(S10-(S11-(S12-(S13-(S14-(S15-S16/nn)/nn)/nn)/nn)/nn)/nn)/nn)/nn)/nn)/nn)/nn)/nn)/nn)/nn)/nn)/nn)/n;
      if (n > cutoffs[--k]) return (S0-(S1-(S2-(S3-(S4-(S5-(S6-(S7-(S8-(S9-(S10-(S11-(S12-(S13-(S14-(S15-(S16-S17/nn)/nn)/nn)/nn)/nn)/nn)/nn)/nn)/nn)/nn)/nn)/nn)/nn)/nn)/nn)/nn)/nn)/n;
      if (n > cutoffs[--k]) return (S0-(S1-(S2-(S3-(S4-(S5-(S6-(S7-(S8-(S9-(S10-(S11-(S12-(S13-(S14-(S15-(S16-(S17-S18/nn)/nn)/nn)/nn)/nn)/nn)/nn)/nn)/nn)/nn)/nn)/nn)/nn)/nn)/nn)/nn)/nn)/nn)/n;
 /* if (n > cutoffs[--k])*/ return (S0-(S1-(S2-(S3-(S4-(S5-(S6-(S7-(S8-(S9-(S10-(S11-(S12-(S13-(S14-(S15-(S16-(S17-(S18-S19/nn)/nn)/nn)/nn)/nn)/nn)/nn)/nn)/nn)/nn)/nn)/nn)/nn)/nn)/nn)/nn)/nn)/nn)/nn)/n;

    } else { // n > cutoffs[nC-6] == 23.5
	nn = n*n;
	k = nC;
	if (n > cutoffs[--k])    return S0/n;
	if (n > cutoffs[--k])    return (S0 -S1/nn)/n;
	if (n > cutoffs[--k])    return (S0-(S1 -S2/nn)/nn)/n;
	if (n > cutoffs[--k])    return (S0-(S1-(S2 -S3/nn)/nn)/nn)/n;
	if (n > cutoffs[--k])    return (S0-(S1-(S2-(S3 -S4/nn)/nn)/nn)/nn)/n;
     /* if (n > cutoffs[--k]) */ return (S0-(S1-(S2-(S3-(S4 -S5/nn)/nn)/nn)/nn)/nn)/n;
    }

  default: error(_("switch() logic programming errors in '%s()'"), "dpq_stirlerr");
  } // switch()
}

SEXP R_dpq_stirlerr(SEXP x_, SEXP stirlerr_v) {
    PROTECT(x_ = isReal(x_) ? x_ : coerceVector(x_, REALSXP));
    R_xlen_t i, n = XLENGTH(x_);
    SEXP ans = PROTECT(allocVector(REALSXP, n));
    double *x = REAL(x_), *r = REAL(ans);
    stirlerr_version_t version = (stirlerr_version_t) asInteger(stirlerr_v);

    for(i=0; i < n; i++) {
	r[i] = dpq_stirlerr(x[i], version);
    }
    UNPROTECT(2);
    return ans;
}

