#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .C calls */
extern void bdstest_main(void *, void *, void *, void *, void *, void *, void *);
extern void C2(void *, void *, void *, void *, void *, void *, void *);
extern void d2(void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void falseNearest(void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void find_nearest(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void follow_points(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void mutual(void *, void *, void *, void *, void *);
extern void stplot(void *, void *, void *, void *, void *, void *, void *, void *);

static const R_CMethodDef CEntries[] = {
    {"bdstest_main",  (DL_FUNC) &bdstest_main,   7},
    {"C2",            (DL_FUNC) &C2,             7},
    {"d2",            (DL_FUNC) &d2,             9},
    {"falseNearest",  (DL_FUNC) &falseNearest,   9},
    {"find_nearest",  (DL_FUNC) &find_nearest,  10},
    {"follow_points", (DL_FUNC) &follow_points, 11},
    {"mutual",        (DL_FUNC) &mutual,         5},
    {"stplot",        (DL_FUNC) &stplot,         8},
    {NULL, NULL, 0}
};

void R_init_fNonlinear(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
