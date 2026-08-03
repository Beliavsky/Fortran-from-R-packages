#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME:
  Check these declarations against the C/Fortran source code.
*/

  /* .C calls */
extern void cvm2d(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void crosscor2d(void *, void *, void *, void *, void *, void *);
extern void cvm3d(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *,void *, void *, void *, void *, void *, void *, void *, void *, void *, void *,void *, void *, void *, void *, void *, void *, void *, void *, void *, void *,void *, void *, void *, void *);
extern void crosscor3d(void *, void *, void *, void *, void *, void *,void *, void *, void *, void *, void *, void *, void *, void *, void *);
static const R_CMethodDef CEntries[] = {
  {"cvm2d",           (DL_FUNC) &cvm2d,      10},
  {"crosscor2d",    (DL_FUNC) &crosscor2d,    6},
  {"cvm3d",           (DL_FUNC) &cvm3d,      34},
  {"crosscor3d",    (DL_FUNC) &crosscor3d,    15},
  {NULL, NULL, 0}
};

void R_init_IndGenErrors(DllInfo *dll)
{
  R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
