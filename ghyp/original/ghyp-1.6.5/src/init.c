#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .C calls */
extern void rgig_R(void *, void *, void *, void *, void *);

static const R_CMethodDef CEntries[] = {
    {"rgig_R", (DL_FUNC) &rgig_R, 5},
    {NULL, NULL, 0}
};

void R_init_ghyp(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
