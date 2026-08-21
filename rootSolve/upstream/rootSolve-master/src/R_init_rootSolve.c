#ifndef R_R_H
# include <R.h>
#endif

#include <Rinternals.h>
#include <stdlib.h> // for NULL

#ifndef R_EXT_DYNLOAD_H_
# include <R_ext/Rdynload.h>
#endif

// define all the globals in this file (only)
#define extern
#include "steady.h"
#undef extern

/* register native routines ------------------------------------------------ */

/* 
   ToDo: 
   - consider replacing SEXP with REALSXP, INTSXP, STRSXP (character), VEXSXP (lists) etc.
   - unlock
*/

/* .Call calls */
extern SEXP call_dsteady(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, 
                       SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, 
                       SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP call_stsparse(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, 
                       SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, 
                       SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP call_lsode(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, 
                       SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, 
                       SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"call_dsteady",    (DL_FUNC) &call_dsteady,    26},
    {"call_stsparse",   (DL_FUNC) &call_stsparse,   27},
    {"call_lsode",      (DL_FUNC) &call_lsode,      26},
    {NULL, NULL, 0}
};  

void R_init_rootSolve(DllInfo *dll) {

  // thpe 2017-03-22, register entry points
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  
  // the following two lines protect against accidentially finding entry points
  R_useDynamicSymbols(dll, FALSE); // disable dynamic searching
}  
