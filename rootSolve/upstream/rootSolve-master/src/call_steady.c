/* Patterned on code call_lsoda.c from package odesolve */

#include <time.h>
#include <string.h>
#include "steady.h"   
#include "externalptr.h"

void F77_NAME(dsteady)(void (*)(int *, double *, double *, double *, double*, int*),
		     int *, int *, double *, double *, double *, double *,
		     int *, int *, int *, int*, double *, double *, double *, int*,
		     void (*)(int *, double *, double *, int *,
			            int *, double *, int *, double*, int*),
		     int *, int *, int *, int *, double *, double *, double *, int *,
         double *, int *, double *, int *);

static C_deriv_func_type *derivb;

                                        
static void C_steady_derivs (int *neq, double *t, double *y, double *ydot,
  double *yout, int *iout)
{

  int i;
  SEXP R_fcall, ans;     

  REAL(Time)[0] = *t;
  for (i = 0; i < *neq; i++)  REAL(Y)[i] = y[i];

  PROTECT(R_fcall = lang3(Rst_deriv_func,Time,Y)) ;
  PROTECT(ans = eval(R_fcall, Rst_envir))         ;

  for (i = 0; i < *neq; i++)	ydot[i] = REAL(VECTOR_ELT(ans,0))[i];
  UNPROTECT(2);      

}

/* assumes 1-D multi-species model; rearrange state vars and rates of change */

static void C_steady_derivs2(int *neq, double *t, double *y, double *ydot,
  double *yout, int *iout)
{
  int i, j, ii;
     ii = 0;
     for (i = 0; i < ndim; i++)         /* rearrange states */
      { for (j = 0; j < nspec; j++)
        /*y2[j* ndim+i] = y[i* nspec+j];   */
        y2[j* ndim+i] = y[ii++];   
        }
        
     derivb (neq, t, y2, dy2, yout, iout) ;
     ii = 0;
     for (i = 0; i < ndim; i++)         /* rearrange rates of change */
      { for (j = 0; j < nspec; j++)
       /*  ydot[i* nspec+j]=dy2[j* ndim+i];   */
         ydot[ii++]=dy2[j* ndim+i];            
         }
        
     }


static void C_steady_jac (int *neq, double *t, double *y, int *ml,
		    int *mu, double *pd, int *nrowpd, double *RPAR, int *IPAR)
{
  int i, j;
  SEXP R_fcall, ans;

  REAL(Time)[0] = *t;
  for (i = 0; i < *neq; i++)  REAL(Y)[i] = y[i];

  PROTECT(R_fcall = lang3(Rst_jac_func,Time,Y));  
  PROTECT(ans = eval(R_fcall, Rst_envir));        

  for (i = 0; i < *neq; i++)
    for (j = 0; j < *nrowpd; j++)
    {
      pd[i * (*nrowpd) + j] = REAL(ans)[i * (*neq) + j];
    }
  UNPROTECT(2);
}

typedef void C_jac_func_type(int *, double *, double *, int *,
		                  int *, double *, int *, double *, int *);

SEXP call_dsteady(SEXP y, SEXP time, SEXP func, SEXP parms, SEXP forcs, 
		SEXP chtol, SEXP atol, SEXP rtol, SEXP itol, SEXP rho, SEXP jacfunc, 
    SEXP initfunc, SEXP initforc,  
		SEXP verbose, SEXP mf, SEXP BU, SEXP BD, SEXP nIter, SEXP Posit, SEXP Pos,
    SEXP nOut,SEXP nAbd, SEXP nSpec, SEXP nDim, SEXP Rpar, SEXP Ipar)
{
  SEXP   yout, RWORK, IWORK;
  int    j, k, ny, maxit, isSteady;
  double *svar, *dy, *beta, *alpha, tin, *Atol, *Rtol, Chtol;
  double *copyvar, *delt, *precis, *ewt ;
  int    neq, bu, bd, jt, niter, mflag, nabd, posit, *pos, ipos, *indx, Itol;
  int    len, isDll, rearrange;
    
  C_deriv_func_type *derivs;
  C_jac_func_type   *jac=NULL;

  int nprot = 0;
  jt = INTEGER(mf)[0];        

  bu = INTEGER(BU)[0];        
  bd = INTEGER(BD)[0];  
  nabd = INTEGER(nAbd)[0];      
  ny = LENGTH(y);
  Itol = INTEGER(itol)[0];
  maxit = INTEGER(nIter)[0];  
  posit = INTEGER(Posit)[0];

  ipos = LENGTH(Pos);
  pos = (int *) R_alloc(ipos, sizeof(int));
    for (j = 0; j < ipos; j++) pos[j] = INTEGER(Pos)[j];

  neq = ny; 
  mflag = INTEGER(verbose)[0];

  rearrange = 0;
  if (jt == 0)   /* state variables and rate of changes need rearranging*/
  {
   jt =25;
   rearrange = 1;
  } 

  if (inherits(func, "NativeSymbol"))  /* function is a dll */
   isDll = 1;
  else                        /* function is not a dll */
   isDll = 0;

  /* initialise output ... */
  initOut(isDll, neq, nOut, Rpar, Ipar);

  /* initialise global variables... */

  PROTECT(Time = NEW_NUMERIC(1))                   ;    nprot++;
  PROTECT(Y = allocVector(REALSXP, neq))           ;    nprot++;

  /* copies of all variables that will be changed in the FORTRAN subroutine */

  indx = (int *) R_alloc(neq, sizeof(int));
    for (j = 0; j < ny; j++) indx[j] = 0;
 
  svar = (double *) R_alloc(neq, sizeof(double));
    for (j = 0; j < ny; j++) svar[j] = REAL(y)[j];

  if (jt > 23)
  {
   delt = (double *) R_alloc(neq, sizeof(double));
    for (j = 0; j < ny; j++) delt[j] = 0; 
  } else delt = (double *) R_alloc(1, sizeof(double));
  
  dy = (double *) R_alloc(neq, sizeof(double));
    for (j = 0; j < ny; j++) dy[j] = 0; 

  ewt = (double *) R_alloc(neq, sizeof(double));
    for (j = 0; j < ny; j++) ewt[j] = 0; 

  beta = (double *) R_alloc(neq, sizeof(double));
    for (j = 0; j < ny; j++) beta[j] = 0; 

  copyvar = (double *) R_alloc(neq, sizeof(double));
    for (j = 0; j < ny; j++) copyvar[j] = 0; 

  alpha = (double *) R_alloc(neq*nabd, sizeof(double));
    for (j = 0; j < neq*nabd; j++) alpha[j] = 0; 
  
  len = LENGTH(atol);  
  Atol = (double *) R_alloc(len, sizeof(double));
    for (j = 0; j < len; j++) Atol[j] = REAL(atol)[j];

  len = LENGTH(rtol);  
  Rtol = (double *) R_alloc(len, sizeof(double));
    for (j = 0; j < len; j++) Rtol[j] = REAL(rtol)[j];

  Chtol = REAL(chtol)[0];

  precis =(double *) R_alloc(maxit, sizeof(double));
    for (j = 0; j < maxit; j++) precis[j] = 0;
  
  PROTECT(yout = allocVector(REALSXP,ntot))    ;     nprot++;


 /* The initialisation routine */
  //initParms(initfunc, parms);
  if (initfunc != NA_STRING) {
    
    if (inherits(initfunc, "NativeSymbol"))  {
      C_init_func_type *initializer;
      
      PROTECT(st_gparms = parms);     nprot++;
      initializer = (C_init_func_type *) R_ExternalPtrAddrFn_(initfunc);
      initializer(Initstparms);
    }
  }
  
//  initForcs(initforc, forcs);
  if (initforc != NA_STRING) {
    if (inherits(initforc, "NativeSymbol"))  {
      C_init_func_type *initializer;
    
      PROTECT(st_gforcs = forcs);     nprot++;
      initializer = (C_init_func_type *) R_ExternalPtrAddrFn_(initforc);
      initializer(Initstforcs);
    }
  }

 /* pointers to functions derivs and jac, passed to the FORTRAN subroutine */

  if (isDll==1) 
    { 
    if (rearrange == 0)
      {
      derivs = (C_deriv_func_type *) R_ExternalPtrAddrFn_(func);
      } else {
       nspec=INTEGER(nSpec)[0];
       ndim =INTEGER(nDim)[0];
       derivs = (C_deriv_func_type *) C_steady_derivs2; 
       derivb = (C_deriv_func_type *) R_ExternalPtrAddrFn_(func);

       y2 = (double *) R_alloc(neq, sizeof(double));
       dy2 = (double *) R_alloc(neq, sizeof(double));   
      }
    } else /* not a DLL */
    {  derivs = (C_deriv_func_type *) C_steady_derivs;  
      PROTECT(Rst_deriv_func = func);     nprot++;
      PROTECT(Rst_envir = rho);    nprot++;
    } 
    
   if (!isNull(jacfunc))
    {
      if (inherits(jacfunc,"NativeSymbol"))
     	{
	    jac = (C_jac_func_type *) R_ExternalPtrAddrFn_(jacfunc);
	    } else {
	    Rst_jac_func = jacfunc;
	    jac = C_steady_jac;
	    }
    }

    tin = REAL(time)[0];

      
	  F77_CALL(dsteady) (derivs, &neq, &nabd, &tin, svar, beta, alpha,
			   &jt, &bu, &bd, &maxit, &Chtol, Atol, Rtol, &Itol, jac, &posit, pos, &ipos,
         &isSteady,delt, copyvar, ewt, indx, precis, &niter, out, ipar);

	  for (j = 0; j < ny; j++)
	    REAL(yout)[j] = svar[j];
   
	  if (isOut == 1) 
    {
        derivs (&neq, &tin, svar, dy, out, ipar) ;
	      for (j = 0; j < nout; j++)
	       REAL(yout)[j + ny] = out[j]; 
    }
 

  PROTECT(RWORK = allocVector(REALSXP, niter));    nprot++;
  for (k = 0;k<niter;k++) REAL(RWORK)[k] = precis[k];

  if (mflag == 1) Rprintf("mean residual derivative %g\n",precis[niter-1]);

  setAttrib(yout, install("precis"), RWORK);    

  PROTECT(IWORK = allocVector(INTSXP, 1));    nprot++;
                          INTEGER(IWORK)[0] = isSteady;
  
  setAttrib(yout, install("steady"), IWORK);    
       
  UNPROTECT(nprot);
  return(yout);
}
