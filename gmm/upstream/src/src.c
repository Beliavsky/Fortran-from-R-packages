#include <stdlib.h>
#include <R_ext/Rdynload.h>
#include "gmm.h"

static const R_FortranMethodDef fortranMethods[] = {
  {"wu", (DL_FUNC) &F77_SUB(wu), 9},
  {"lamcuep", (DL_FUNC) &F77_SUB(lamcuep), 9},
  {NULL, NULL, 0}
};

void R_init_gmm(DllInfo *dll)
   {
     R_registerRoutines(dll,
			NULL, NULL, 
			fortranMethods, NULL);
     R_useDynamicSymbols(dll, FALSE);
   }

