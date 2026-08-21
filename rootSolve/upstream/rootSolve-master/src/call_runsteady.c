#include <time.h>
#include <string.h>
#include "steady.h"
#include "externalptr.h"

/* definition of the call  to the fortran functions - in file lsode.f*/
/* bug fix as suggested by Alejandro Morales 4/04/2014*/

void F77_NAME(dlsode)(void (*)(int *, double *, double *, double *, double *, int *),
		     int *, double *, double *, double *,
		     int *, double *, double *, int *, int *,
		     int *, double *,int *,int *, int *,
		     void (*)(int *, double *, double *, int *,
			      int *, double *, int *, double *, int *),
		     int *, double *, int *);
static C_deriv_func_type *derivb;

void F77_NAME(dlsodes)(void (*)(int *, double *, double *, double *, double *, int *),
             int *, double *, double *, double *,
             int *, double *, double *, int *, int *,
             int *, double *, int *, int *, int *,  int *,  /* extra 'double'; is integer in fortran, thpe: now u_work.iwk */
             void (*)(int *, double *, double *, int *,
                      int *, int *, double *, double *, int *),     /* jacvec */
             int *, double *, int *);

/* interface between fortran function call and R function 
   Fortran code calls C_ode_derivs(N, t, y, ydot, yout, iout) 
   R code called as lsode_deriv_func(time, y) and returns ydot 
   Note: passing of parameter values and "..." is done in R-function runsteady*/

static void C_ode_derivs (int *neq, double *t, double *y, 
                          double *ydot, double *yout, int *iout)
{
  int i;
  SEXP R_fcall, ans;

                              REAL(Time)[0] = *t;
  for (i = 0; i < *neq; i++)  REAL(Y)[i] = y[i];

  PROTECT(R_fcall = lang3(lsode_deriv_func,Time,Y));   
  PROTECT(ans = eval(R_fcall, lsode_envir));           

  for (i = 0; i < *neq; i++)   ydot[i] = REAL(VECTOR_ELT(ans,0))[i];

  UNPROTECT(2);
}

/* interface between fortran call to jacobian and R function */

/* assumes 1-D multi-species model; rearrange state vars and rates of change */

static void C_ode_derivs2(int *neq, double *t, double *y, double *ydot,
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


static void C_ode_jac (int *neq, double *t, double *y, int *ml,
		    int *mu, double *pd, int *nrowpd, double *yout, int *iout)
{
  int i;
  SEXP R_fcall, ans;

                             REAL(Time)[0] = *t;
  for (i = 0; i < *neq; i++) REAL(Y)[i] = y[i];

  PROTECT(R_fcall = lang3(lsode_jac_func,Time,Y));    
  PROTECT(ans = eval(R_fcall, lsode_envir));          

  for (i = 0; i < *neq * *nrowpd; i++)  pd[i] = REAL(ans)[i];

  UNPROTECT(2);
}


/* give name to data types */
typedef void C_jac_func_type  (int *, double *, double *, int *,
		                    int *, double *, int *, double *, int *);
typedef void C_jac_vec_type   (int *, double *, double *, int *,
                            int *, int *, double *, double *, int *);

/* MAIN C-FUNCTION, CALLED FROM R-code */

SEXP call_lsode(SEXP y, SEXP times, SEXP func, SEXP parms, SEXP forcs, 
    SEXP stol, SEXP rtol, SEXP atol, SEXP rho, SEXP tcrit, 
    SEXP jacfunc, SEXP initfunc, SEXP initforc,
		SEXP verbose, SEXP iTask, SEXP rWork, SEXP iWork, SEXP jT, 
    SEXP nOut, SEXP lRw, SEXP lIw, SEXP nSpec, SEXP nDim, 
    SEXP Rpar, SEXP Ipar, SEXP Solver )

{
/******************************************************************************/
/******                   DECLARATION SECTION                            ******/
/******************************************************************************/

/* These R-structures will be allocated and returned to R*/
  SEXP yout, RWORK, ISTATE;    

  int  i, j, k, latol, lrtol, lrw, liw, maxit, rearrange;
  double *xytmp, *rwork, tin, tout, *Atol, *Rtol, Stol, *dy, ss, sumder=0.;
  int neq, itol, itask, istate, iopt, *iwork, jt, mflag, is;
  int isDll, Steady, solver;
  
  C_deriv_func_type *derivs;
  C_jac_func_type   *jac=NULL;
  C_jac_vec_type    *jac_vec=NULL;

  /* memory overlay, KS, TP 2019-07-03 */
  union t_work {
	double * rwk;
    int * iwk;	
  } u_work;

/******************************************************************************/
/******                         STATEMENTS                               ******/
/******************************************************************************/

/*                      #### initialisation ####                              */    
  int nprot = 0;

  solver   = INTEGER(Solver)[0];

  jt  = INTEGER(jT)[0];         /* method flag */
  rearrange = 0;
  if (jt == 0)   /* state variables and rate of changes need rearranging*/
  {
   jt =25;
   rearrange = 1;
  } 

  neq = LENGTH(y);              /* number of equations */ 
  
  mflag = INTEGER(verbose)[0];
  tin  = REAL(times)[0];        /* start and end time */
  tout = REAL(times)[1];
  Stol = REAL(stol)[0];         /* steady-state tolerance */ 

  if (inherits(func, "NativeSymbol"))  /* function is a dll */
    isDll = 1;
  else                              /* function is not a dll */
   isDll = 0;

/* initialise output ... */
  initOut(isDll, neq, nOut, Rpar, Ipar);

/* copies of all variables that will be changed in the FORTRAN subroutine */

  xytmp = (double *) R_alloc(neq, sizeof(double));
  for (j = 0; j < neq; j++) xytmp[j] = REAL(y)[j];
 
  latol = LENGTH(atol);
  Atol = (double *) R_alloc((int) latol, sizeof(double));

  lrtol = LENGTH(rtol);
  Rtol = (double *) R_alloc((int) lrtol, sizeof(double));

  liw = INTEGER (lIw)[0];
  iwork = (int *) R_alloc(liw, sizeof(int));
     for (j=0; j<LENGTH(iWork); j++) iwork[j] = INTEGER(iWork)[j];
  maxit = iwork[5];
  
  lrw = INTEGER(lRw)[0];

  // ks, tp 2019-07-03
  //rwork = (double *) R_alloc(lrw, sizeof(double));
  //for (j=0; j<length(rWork); j++) rwork[j] = REAL(rWork)[j];
  u_work.rwk = (double *) R_alloc(lrw, sizeof(double));
  for (j=0; j<length(rWork); j++) u_work.rwk [j] = REAL(rWork)[j];
  rwork = u_work.rwk;

/* initialise global R-variables... */

  PROTECT(Time = NEW_NUMERIC(1));                  nprot++;
  PROTECT(Y = allocVector(REALSXP,(neq)));         nprot++;
  PROTECT(yout = allocVector(REALSXP,ntot));       nprot++;

 /* The initialisation routine */
  // initParms(initfunc, parms);
  if (initfunc != NA_STRING) {
    
    if (inherits(initfunc, "NativeSymbol"))  {
      C_init_func_type *initializer;
      
      PROTECT(st_gparms = parms);     nprot++;
      initializer = (C_init_func_type *) R_ExternalPtrAddrFn_(initfunc);
      initializer(Initstparms);
    }
   }
  
  // initForcs(initforc, forcs);
  if (initforc != NA_STRING) {
    if (inherits(initforc, "NativeSymbol"))  {
      C_init_func_type *initializer;
    
      PROTECT(st_gforcs = forcs);     nprot++;
      initializer = (C_init_func_type *) R_ExternalPtrAddrFn_(initforc);
      initializer(Initstforcs);
    }
  }

/* pointers to functions derivs, jac, passed to FORTRAN */
  dy = (double *) R_alloc(neq, sizeof(double));
  for (j = 0; j < neq; j++) dy[j] = 0.; 

  if (isDll ==1) 
    { /* DLL address passed to fortran */
    if (rearrange == 0)
      {
      derivs = (C_deriv_func_type *) R_ExternalPtrAddrFn_(func);  
      /* no need to communicate with R - but output variables set here */
      } else {
       nspec=INTEGER(nSpec)[0];
       ndim =INTEGER(nDim)[0];
       derivs = (C_deriv_func_type *) C_ode_derivs2; 
       derivb = (C_deriv_func_type *) R_ExternalPtrAddrFn_(func);
       y2 = (double *) R_alloc(neq, sizeof(double));
       dy2 = (double *) R_alloc(neq, sizeof(double));   
      }	  
    } else {
      /* interface function between fortran and R passed to Fortran*/ 
      derivs = (C_deriv_func_type *) C_ode_derivs; 
      /* needed to communicate with R */
      lsode_deriv_func = func;

    }
      lsode_envir = rho;

  if (!isNull(jacfunc)) 
    {
      if (isDll ==1)
	    {
	     jac = (C_jac_func_type *) R_ExternalPtrAddrFn_(jacfunc);
	    } else  {
	     lsode_jac_func = jacfunc;
	     jac = C_ode_jac;
	    }
    }


/* tolerance specifications */
  if (latol == 1 && lrtol == 1 ) itol = 1;
  if (latol  > 1 && lrtol == 1 ) itol = 2;
  if (latol == 1 && lrtol  > 1 ) itol = 3;
  if (latol  > 1 && lrtol  > 1 ) itol = 4;
  for (j = 0; j < lrtol; j++) Rtol[j] = REAL(rtol)[j];
  for (j = 0; j < latol; j++) Atol[j] = REAL(atol)[j];

  itask = INTEGER(iTask)[0];   
  istate = 1;

  iopt = 0;
  ss = 0.;
  is = 0 ;
  for (i = 5; i < 8 ; i++) ss = ss+rwork[i];
  for (i = 5; i < 10; i++) is = is+iwork[i];
  if (ss >0 || is > 0) iopt = 1; /* non-standard input */

  Steady = 0;
/*                     ####   main time loop   ####                           */    
  for (i = 0; i < maxit; i++)
	{  /* one step */
    if (solver == 0)
     F77_CALL(dlsode) (derivs, &neq, xytmp, &tin, &tout,
			   &itol, Rtol, Atol, &itask, &istate, &iopt, rwork,
			   &lrw, iwork, &liw, jac, &jt, out, ipar); 
    else
     F77_CALL(dlsodes) (derivs, &neq, xytmp, &tin, &tout,
         &itol, Rtol, Atol, &itask, &istate, &iopt, u_work.rwk,
         &lrw, iwork, &liw, u_work.iwk, jac_vec, &jt, out, ipar);    
     /* check steady-state */
    sumder = 0. ;
    derivs (&neq, &tin, xytmp, dy, out, ipar) ;
    for (j = 0; j < neq; j++) sumder = sumder+fabs(dy[j]); 

    if (sumder/neq < Stol) {
     Steady = 1;
     break;
    }
    if (tin >= tout) break;
    
    /* errors? */
	  if (istate == -2)
	    {
	      for (j = 0; j < lrtol; j++) Rtol[j] *= 10.0;
	      for (j = 0; j < latol; j++) Atol[j] *= 10.0;
	      warning("Excessive precision requested.  `rtol' and `atol' have been scaled upwards by the factor %g\n",10.0);
	      istate = 3;
	    }
	  if (istate == -1) 
     {
      warning("an excessive amount of work (> maxsteps ) was done, but integration was successful - increase maxsteps");
     }
    if (istate == -3)
	   {
	    UNPROTECT(nprot);
	    error("Illegal input to lsode\n");
	   }
	} 

  /* If here: hopefully steady-state is reached */  
 
	  for (j = 0; j < neq; j++)
	    REAL(yout)[j] = xytmp[j];
    if (isOut == 1) 
      {
      for (j = 0; j < nout; j++)  REAL(yout)[j + neq] = out[j];   
        /* above = bug fix as suggested by Alejandro Morales; was yout[j+neq+1]*/
	    }
/*                    ####  an error occurred   ####                          */    
  if (istate < 0 ) warning("Returning early.  Results are accurate, as far as they go\n");

  PROTECT(ISTATE = allocVector(INTSXP, 24));    nprot++;
  for (k = 0;k<22;k++) INTEGER(ISTATE)[k+1] = iwork[k];
  INTEGER(ISTATE)[0] = istate;
  INTEGER(ISTATE)[23] = Steady;
       
  PROTECT(RWORK = allocVector(REALSXP, 7));    nprot++;
  for (k = 0;k<5;k++) REAL(RWORK)[k] = rwork[k+10];

  REAL(RWORK)[5] = sumder/neq;
  REAL(RWORK)[6] = tin;
  if (mflag == 1) Rprintf("mean residual derivative %g\n",sumder/neq);
  setAttrib(yout, install("rstate"), RWORK);    
  setAttrib(yout, install("istate"), ISTATE);    

/*                       ####   termination   ####                            */    
  UNPROTECT(nprot);
  return(yout);
}

