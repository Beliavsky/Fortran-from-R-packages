#include <R_ext/RS.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .C calls */
extern void dwt(void *, void *, void *, void *, void *, void *, void *);
extern void hosking(void *, void *, void *);
extern void idwt(void *, void *, void *, void *, void *, void *, void *);
extern void imodwt(void *, void *, void *, void *, void *, void *, void *, void *);
extern void modwt(void *, void *, void *, void *, void *, void *, void *, void *);
extern void three_D_dwt(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void three_D_idwt(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void three_D_imodwt(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void three_D_modwt(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void two_D_dwt(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void two_D_idwt(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void two_D_imodwt(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void two_D_modwt(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);

static const R_CMethodDef CEntries[] = {
    {"dwt",            (DL_FUNC) &dwt,             7},
    {"hosking",        (DL_FUNC) &hosking,         3},
    {"idwt",           (DL_FUNC) &idwt,            7},
    {"imodwt",         (DL_FUNC) &imodwt,          8},
    {"modwt",          (DL_FUNC) &modwt,           8},
    {"three_D_dwt",    (DL_FUNC) &three_D_dwt,    15},
    {"three_D_idwt",   (DL_FUNC) &three_D_idwt,   15},
    {"three_D_imodwt", (DL_FUNC) &three_D_imodwt, 16},
    {"three_D_modwt",  (DL_FUNC) &three_D_modwt,  16},
    {"two_D_dwt",      (DL_FUNC) &two_D_dwt,      10},
    {"two_D_idwt",     (DL_FUNC) &two_D_idwt,     10},
    {"two_D_imodwt",   (DL_FUNC) &two_D_imodwt,   11},
    {"two_D_modwt",    (DL_FUNC) &two_D_modwt,    11},
    {NULL, NULL, 0}
};

void R_init_waveslim(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
