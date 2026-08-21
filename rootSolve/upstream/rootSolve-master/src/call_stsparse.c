/* SPARSE steady-state solver */

#include <time.h>
#include <string.h>
#include "steady.h"   
#include "externalptr.h"
                           
void F77_NAME(dsparse)(void (*)(int *, double *, double *, double *, double*, int*),
             int *, int *, int *, double *, double *, double *, double *, double *,
             double *, double *, double *, int*, int*, int*, int*, int*, 
         int*, int*, int*, int*, 
         int *, double *, double *, double *, int *,
             int *, int *, int *, int *, double *, int*,
         int *, double *, int *, int*, int*);

void F77_NAME(dsparsekit)(void (*)(int *, double *, double *, double *, double*, int*),
             int *, int *, int *, double *, double *, double *, double *, double *,
             double *, double *, int*, int*, int*, int*, int*,  int*,  int*,  int*,  int*, 
         int *, double *, double *, double *, int *, int *,
             int *, int *, int *, double *, int*,
         int *, double *, int *, int*, double *, double *, int *, 
         int *, int*, double *, double *, int *);

/*==================================================
 extracting elements from a list
===================================================*/

SEXP getListElement(SEXP list, const char *str) {
  SEXP elmt = R_NilValue, names = getAttrib(list, R_NamesSymbol);
  int i;
  for (i = 0; i < length(list); i++)
    if (strcmp(CHAR(STRING_ELT(names, i)), str) == 0) {
      elmt = VECTOR_ELT(list, i);
      break;
    }
  return elmt;
}

/*==================================================
 interface FORTRAN and R-code
===================================================*/

static void C_stsparse_derivs (int *neq, double *t, double *y, double *ydot, 
                            double *yout, int *iout)
{
  int i;
  SEXP R_fcall, ans;     

  REAL(Time)[0] = *t;
  for (i = 0; i < *neq; i++)  REAL(Y)[i] = y[i];

  PROTECT(R_fcall = lang3(stsparse_deriv_func,Time,Y)) ;
  PROTECT(ans = eval(R_fcall, stsparse_envir))         ;

  for (i = 0; i < *neq; i++)    ydot[i] = REAL(VECTOR_ELT(ans,0))[i];
  UNPROTECT(2);      

}

SEXP call_stsparse(SEXP y, SEXP time, SEXP func, SEXP parms, SEXP forcs, SEXP chtol, 
        SEXP atol, SEXP rtol, SEXP itol, SEXP rho, SEXP initfunc, SEXP initforc,
        SEXP verbose, SEXP NNZ, SEXP NSP, SEXP MAXG, SEXP nIter, SEXP Posit,
    SEXP Pos, SEXP nOut, SEXP Rpar, SEXP Ipar, SEXP Type, SEXP Ian, SEXP Jan,
    SEXP Met, SEXP Option)
{
  SEXP   yout, RWORK, IWORK, IANX, JANX;
  int    j, k, ny, maxit, isSteady, method, lenplufac, lenplumx, lfill;
  double *svar, *dsvar, *beta, *alpha, tin, *Atol, *Rtol, Chtol;
  double *x, *precis, *ewt, droptol, permtol;
  int    neq, nnz, nsp, maxg, niter, mflag, posit, TotN, ipos, Itol, type;
  int    *ian, *jan, *igp, *jgp, *dims, *pos;
  int    len, isDll, ilumethod;

  double *rsp= NULL, *plu= NULL,  *rwork= NULL;
  int    *R= NULL, *C= NULL, *IC= NULL, *indDIM = NULL;
  int    *isp= NULL, *iwork= NULL, *iperm= NULL, *jlu= NULL, *ju= NULL; 

  C_deriv_func_type *derivs;
  int nprot = 0;

  nnz   = INTEGER(NNZ)[0];
  nsp   = INTEGER(NSP)[0];  
  maxg  = INTEGER(MAXG)[0];  
  ny    = LENGTH(y);
  Itol  = INTEGER(itol)[0];
  maxit = INTEGER(nIter)[0];  
  type  = INTEGER(Type)[0];
  method = INTEGER(Met)[0];

  posit = INTEGER(Posit)[0];   /* positivity of state variables: either specified at once, or via a vector..*/
  ipos = LENGTH(Pos);
  pos = (int *) R_alloc(ipos, sizeof(int));
    for (j = 0; j < ipos; j++) pos[j] = INTEGER(Pos)[j];

  neq   = ny; 
  mflag = INTEGER(verbose)[0];

  if (inherits(func, "NativeSymbol"))  /* function is a dll */
     isDll = 1;
  else
     isDll = 0;
  if (nout > 0) isOut = 1; 

  /* initialise output ... */
  initOut(isDll, neq, nOut, Rpar, Ipar);

  /* initialise global variables... */
  PROTECT(Time = NEW_NUMERIC(1))                   ;    nprot++;
  PROTECT(Y = allocVector(REALSXP, neq))           ;    nprot++;

  /* copies of all variables that will be changed in the FORTRAN subroutine */
  if ( (method == 1) | (method == 10)) {       /* yale sparse matrix solver */
    R = (int *) R_alloc(neq, sizeof(int));
     for (j = 0; j < ny; j++) R[j] = 0;
 
    C = (int *) R_alloc(neq, sizeof(int));
     for (j = 0; j < ny; j++) C[j] = 0;

    IC = (int *) R_alloc(neq, sizeof(int));
     for (j = 0; j < ny; j++) IC[j] = 0;

    rsp = (double *) R_alloc(nsp, sizeof(double));
     for (j = 0; j < nsp; j++) rsp[j] = 0.;

    isp = (int *) R_alloc(2*nsp, sizeof(int));
     for (j = 0; j < 2*nsp; j++) isp[j] = 0;
  } else {                                      /* sparskit matrix solver */
    /* get options */
    lenplufac = INTEGER(getListElement(Option, "lenplufac"))[0];
    lfill     = INTEGER(getListElement(Option, "fillin")   )[0];
    droptol   = REAL   (getListElement(Option, "droptol")  )[0];
    permtol   = REAL   (getListElement(Option, "permtol")  )[0];

    ilumethod = method - 1; /* 1 = ilut, 2 = ilutp */
    lenplumx = nnz + lenplufac*neq;
    
    jlu = (int *) R_alloc(lenplumx, sizeof(int));
     for (j = 0; j < lenplumx; j++) jlu[j] = 0;

    ju = (int *) R_alloc(neq, sizeof(int));
     for (j = 0; j < neq; j++) ju[j] = 0;

    iwork = (int *) R_alloc(2*neq, sizeof(int));
     for (j = 0; j < 2*neq; j++) iwork[j] = 0;

    iperm = (int *) R_alloc(2*neq, sizeof(int));
     for (j = 0; j < 2*neq; j++) iperm[j] = 0;
 
    plu = (double *) R_alloc(lenplumx, sizeof(double));
     for (j = 0; j < lenplumx; j++) plu[j] = 0.;

    rwork = (double *) R_alloc(neq, sizeof(double));
     for (j = 0; j < neq; j++) rwork[j] = 0.;
  }

  dims = (int *) R_alloc(7, sizeof(int));   /* 7 is maximal amount */
    for (j = 0; j < 7; j++) dims[j] = 0;

  svar = (double *) R_alloc(neq, sizeof(double));
    for (j = 0; j < ny; j++) svar[j] = REAL(y)[j];

  dsvar = (double *) R_alloc(neq, sizeof(double));
    for (j = 0; j < ny; j++) dsvar[j] = 0; 

  beta = (double *) R_alloc(neq, sizeof(double));
    for (j = 0; j < ny; j++) beta[j] = 0; 

  x = (double *) R_alloc(neq, sizeof(double));
    for (j = 0; j < ny; j++) x[j] = 0; 

  alpha = (double *) R_alloc(nnz, sizeof(double));
    for (j = 0; j < nnz; j++) alpha[j] = 0; 

  ewt = (double *) R_alloc(neq, sizeof(double));
    for (j = 0; j < ny; j++) ewt[j] = 0; 

  ian = (int *) R_alloc(neq+1, sizeof(int));
   if (type == 0) {
     for (j = 0; j < neq+1; j++) ian[j] = INTEGER(Ian)[j];} 
   else {
     for (j = 0; j < neq+1; j++) ian[j] = 0;}

  jan = (int *) R_alloc(nnz, sizeof(int));
   if (type == 0){ 
     for (j = 0; j < nnz; j++) jan[j] = INTEGER(Jan)[j];} 
   else {
     for (j = 0; j < nnz; j++) jan[j] = 0;}
   
  /* 1-D, 2-D, 3-D problem:  */
  if (type == 2)        /* 1=ncomp,2:dim(x), 3: cyclic(x)*/
    for (j = 0; j<3 ; j++) dims[j] = INTEGER(NNZ)[j+1];
  else if (type == 3)   /* 1=ncomp,2-3:dim(x,y), 4-5: cyclic(x,y)*/
    for (j = 0; j<5 ; j++) dims[j] = INTEGER(NNZ)[j+1];
  else if (type == 4)   /* 1=ncomp,2-4:dim(x,y,z), 5-7: cyclic(x,y,z)*/
    for (j = 0; j<7 ; j++) dims[j] = INTEGER(NNZ)[j+1];
  else if (type == 30)  { /* same as type 3 (2-D) but with mapping */
    for (j = 0; j<5 ; j++) dims[j] = INTEGER(NNZ)[j+1];
    TotN = INTEGER(NNZ)[6];
	  indDIM = (int *) R_alloc(TotN, sizeof(int));
    for (j = 0; j < TotN ; j++) indDIM[j] = INTEGER(NNZ)[j+7];
  }  else if (type == 40)  { /* same as type 4 (3-D) but with mapping */
    for (j = 0; j<7 ; j++) dims[j] = INTEGER(NNZ)[j+1];
    TotN = INTEGER(NNZ)[8];
	  indDIM = (int *) R_alloc(TotN, sizeof(int));
    for (j = 0; j < TotN ; j++) indDIM[j] = INTEGER(NNZ)[j+9];
  }

  igp = (int *) R_alloc(maxg+1, sizeof(int));
    for (j = 0; j < maxg+1; j++) igp[j] = 0;
  
  jgp = (int *) R_alloc(neq, sizeof(int));
    for (j = 0; j < neq; j++) jgp[j] = 0;

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

  if (isDll)
    {
      derivs = (C_deriv_func_type *) R_ExternalPtrAddrFn_(func);

    } else {  derivs = (C_deriv_func_type *) C_stsparse_derivs;  
      PROTECT(stsparse_deriv_func = func);     nprot++;
      PROTECT(stsparse_envir = rho);       nprot++;
    }
    
    tin = REAL(time)[0];
      
  if ((method == 1) |(method == 10)) {
      F77_CALL(dsparse)(derivs, &neq, &nnz, &nsp, &tin, svar, dsvar, beta, x,
         alpha, ewt, rsp, ian, jan, igp, jgp, &maxg, R, C, IC, isp,
         &maxit,  &Chtol, Atol, Rtol, &Itol, &posit, pos, &ipos, &isSteady,
         precis, &niter, dims, out, ipar, &type, indDIM);
  } else {
      F77_CALL(dsparsekit) (derivs, &neq, &nnz, &nsp, &tin, svar, dsvar, beta, x,
         alpha, ewt, ian, jan, igp, jgp, &maxg, jlu, ju, iwork, iperm,
         &maxit,  &Chtol, Atol, Rtol, &Itol, &posit, pos, &ipos, &isSteady,
         precis, &niter, dims, out, ipar, &type, &droptol, &permtol, &ilumethod,
         &lfill, &lenplumx, plu, rwork, indDIM);
  }
      for (j = 0; j < ny; j++)
        REAL(yout)[j] = svar[j];
   
      if (isOut == 1) 
    {
        derivs (&neq, &tin, svar, dsvar, out, ipar) ;
          for (j = 0; j < nout; j++)
           REAL(yout)[j + ny] = out[j]; 
    }
  
     if (method != 10) {      

       PROTECT(RWORK = allocVector(REALSXP, niter));    nprot++;
       for (k = 0;k<niter;k++) REAL(RWORK)[k] = precis[k];
       if (mflag == 1) Rprintf("mean residual derivative %g\n",precis[niter-1]);
       setAttrib(yout, install("precis"), RWORK);    

       PROTECT(IWORK = allocVector(INTSXP, 4));    nprot++;
                             INTEGER(IWORK)[0]   = isSteady;
       for (k = 0; k<3; k++) INTEGER(IWORK)[k+1] = dims[k];
       setAttrib(yout, install("steady"), IWORK);    
     } else {                  /* returns the sparsity vectors ian and jan */
       
       PROTECT(RWORK = allocVector(REALSXP, 1));     nprot++;
       REAL(RWORK)[0] = 0.;
       setAttrib(yout, install("precis"), RWORK);    

       PROTECT(IWORK = allocVector(INTSXP, 4));    nprot++;
                             INTEGER(IWORK)[0]   = 0;
       for (k = 0; k<3; k++) INTEGER(IWORK)[k+1] = dims[k];
       setAttrib(yout, install("steady"), IWORK);    
      
       PROTECT(IANX = allocVector(INTSXP, neq+1));    nprot++;
       for (k = 0; k<neq+1; k++) INTEGER(IANX)[k] = ian[k];
       setAttrib(yout, install("ian"), IANX);              
     
       PROTECT(JANX = allocVector(INTSXP, dims[0]));    nprot++;
       for (k = 0; k< dims[0]; k++) INTEGER(JANX)[k] = jan[k];
       setAttrib(yout, install("jan"), JANX);                
    } 
  UNPROTECT(nprot);
  return(yout);
}
