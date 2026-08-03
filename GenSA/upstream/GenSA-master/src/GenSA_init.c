#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* .Call calls */
extern SEXP createInstance(void);
extern SEXP execute(void *, void *, void *, void *, void *, void *, void *, void *);
extern SEXP getREnergy(void *);
extern SEXP getRNbFuncCall(void *);
extern SEXP getRTraceMat(void *, void *);
extern SEXP getRTraceMatSize(void *);
extern SEXP getRXMiniVector(void *);
extern SEXP releaseInstance(void *);

static const R_CallMethodDef CallEntries[] = {
    {"createInstance",   (DL_FUNC) &createInstance,   0},
    {"execute",          (DL_FUNC) &execute,          8},
    {"getREnergy",       (DL_FUNC) &getREnergy,       1},
    {"getRNbFuncCall",   (DL_FUNC) &getRNbFuncCall,   1},
    {"getRTraceMat",     (DL_FUNC) &getRTraceMat,     2},
    {"getRTraceMatSize", (DL_FUNC) &getRTraceMatSize, 1},
    {"getRXMiniVector",  (DL_FUNC) &getRXMiniVector,  1},
    {"releaseInstance",  (DL_FUNC) &releaseInstance,  1},
    {NULL, NULL, 0}
};

void R_init_GenSA(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
