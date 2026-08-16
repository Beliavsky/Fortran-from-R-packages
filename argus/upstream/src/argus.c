
/*---------------------------------------------------------------------------*/
/* header files */

#include <R.h>
#include <Rmath.h>
#include <Rdefines.h>

#include "argus.h"




/*****************************************************************************/
/* API                                                                       */
/*****************************************************************************/

/*---------------------------------------------------------------------------*/




SEXP rargusRoU(SEXP sexp_n, SEXP sexp_chi)
/*---------------------------------------------------------------------------*/
/* Draw sample from argus distribution.                                      */
/* Wrapper for do_rargusRoU() with GetRNGstate() ... PutRNGstate()          */
/*                                                                           */
/* Parameters:                                                               */
/*   n ....... sample size (positive integer)                                */
/*   chi   ... parameter for distribution                                    */
/*                                                                           */
/* Return:                                                                   */
/*   random sample of size 'n*length(chi)'                                   */
/*---------------------------------------------------------------------------*/
{
  int n;          /* sample size */
  double *chi;                /* parameter */
  int nchi;                    /* length of chi */
  SEXP sexp_res;  /* results */


  /* extract sample size */
  n = *(INTEGER (AS_INTEGER (sexp_n)));

  /* extract chi vector */
  nchi = length(sexp_chi);
  PROTECT(sexp_chi = AS_NUMERIC(sexp_chi));
  chi = REAL(sexp_chi);

  /* Get state of R uniform PRNG */
  GetRNGstate();

  /* run generator */
  PROTECT(sexp_res = do_rargusRoU(n, chi, nchi));

  /* Return state of PRNG to R */
  PutRNGstate();

  /* return result to R */
  UNPROTECT(2);
  return (sexp_res);

} /* end of rargusRoU() */

/*---------------------------------------------------------------------------*/

SEXP do_rargusRoU(int n, double *chi, int nchi)
/*---------------------------------------------------------------------------*/
/* Draw sample from GIG distribution.                                        */
/* without calling GetRNGstate() ... PutRNGstate()                           */
/*                                                                           */
/* Parameters:                                                               */
/*   n ....... sample size (positive integer)                                */
/*   chi  ... parameter vector                                               */
/*   nchi ... length ofparameter vector                                      */
/*                                                                           */
/* Return:                                                                   */
/*   random sample of size 'n*nchi                                           */
/*---------------------------------------------------------------------------*/
{
  SEXP sexp_res;           /* results */
  double *res;
  int i,j;

  /* check sample size */
  if (n<0) {
    error("sample size 'n' must be non-negative integer.");
  }

  /* allocate array for random sample */
  PROTECT(sexp_res = NEW_NUMERIC(n*nchi));
  res = REAL(sexp_res);

/* setup forRoU */    
  double m,b,ap,am,xp,xm,U,V,Y;
  for (j=0; j<nchi; j++) {
    if(chi[j]<=1){
      m = chi[j]*chi[j]*0.5;
      b = sqrt(sqrt(0.5)*chi[j])*exp(-chi[j]*chi[j]*0.25);
/*    b = 2^-0.25*sqrt(chi[j])*exp(-chi[j]*chi[j]/4) */
      ap = 0.;
      xm = (chi[j]*chi[j]+5.-sqrt(chi[j]*chi[j]*(chi[j]*chi[j]+6.)+25.))*0.25;
      am = (xm-chi[j]*chi[j]*0.5)*sqrt(sqrt(xm)*exp(-xm));
    }else{
      m = 0.5;
      b = exp(-(1.+log(2.))/4.);
/*  b <- 2^-0.25*exp(-1/4) */
      xp = chi[j]*chi[j]*0.5;
      if(xp > 1.5+sqrt(2.))xp = 1.5+sqrt(2.); 
/*  xp <- min(1.5 + sqrt(2),chi*chi*0.5) */
      xm = 1.5 - sqrt(2.);
      ap = (xp - 0.5)*sqrt(sqrt(xp)*exp(-xp));
      am = (xm - 0.5)*sqrt(sqrt(xm)*exp(-xm));
    }
    /* run generator */
    for (i=0; i<n; i++) {
      do {
/*      ++count; */
        U = b * unif_rand();      
        V = am+unif_rand()*(ap-am);
        Y = V/U + m ;
      }                              /* Acceptance/Rejection */
      while ( Y<0. || Y > chi[j]*chi[j]*0.5 || U*U > sqrt(Y)*exp(-Y)  );
      /* store random point */
      res[n*j+i] = sqrt(1.-2.*Y/(chi[j]*chi[j]));
    }
  }


  /* return result */
  UNPROTECT(1);
  return sexp_res;

} /* end of do_rargusRoU() */

