#include <R.h>
#include <Rinternals.h>

#include "cluster.h"

#include <R_ext/Rdynload.h>

#define CDEF(name)  {#name, (DL_FUNC) &name, sizeof(name ## _typ)/sizeof(name ## _typ[0]), name ##_typ}
#define CALLDEF(name, n)  {#name, (DL_FUNC) &name, n}


static R_NativePrimitiveArgType R_bncoef_typ[3] = {
    INTSXP, REALSXP, REALSXP
};

static R_NativePrimitiveArgType cl_clara_typ[34] = {
    /*n:*/ INTSXP, INTSXP, INTSXP, REALSXP, INTSXP, INTSXP, REALSXP, INTSXP,
    /*valmd:*/ REALSXP, INTSXP, INTSXP, /* rng_R: */ LGLSXP, /* pam_like:*/ LGLSXP,
    /*d_flag:*/ INTSXP,
    /*nrepr: */ INTSXP, INTSXP, INTSXP, INTSXP, INTSXP,
    /*radus:*/ REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, INTSXP,
    /*obj: */ REALSXP, REALSXP, REALSXP, REALSXP,  INTSXP, INTSXP,
    /*tmp: */ REALSXP,INTSXP
};

static R_NativePrimitiveArgType cl_fanny_typ[27] = {
    INTSXP, INTSXP, INTSXP, REALSXP, REALSXP,
    /*jdyss: */ INTSXP, REALSXP, INTSXP,  INTSXP, INTSXP, INTSXP,
    /*negbr: */ INTSXP, /*syl: */ REALSXP, REALSXP, REALSXP, REALSXP,
    /*nfuzz: */ INTSXP, REALSXP, REALSXP, REALSXP, REALSXP,
    /*obj: */ REALSXP, INTSXP, REALSXP, REALSXP, REALSXP, INTSXP
};

#ifdef _UNUSED_C_pam
static R_NativePrimitiveArgType cl_pam_typ[24] = {
    INTSXP, INTSXP, INTSXP, REALSXP, REALSXP,
    /*jdyss: */ INTSXP, REALSXP, INTSXP,  INTSXP, INTSXP,
    /*nrepr: */ LGLSXP, INTSXP, /*radus: */ REALSXP, REALSXP, REALSXP, REALSXP,
    /*ttsyl: */ REALSXP, REALSXP, /*medoids*/ INTSXP, INTSXP,  REALSXP, REALSXP, INTSXP,
    /*optim: */ INTSXP
};
#endif

static R_NativePrimitiveArgType spannel_typ[12] = { // ./spannel.c :
    INTSXP, INTSXP, REALSXP, REALSXP, REALSXP, REALSXP,
    /*varss: */ REALSXP, REALSXP, REALSXP, REALSXP, INTSXP, INTSXP
};

static R_NativePrimitiveArgType sildist_typ[] = {
    REALSXP, INTSXP, INTSXP, INTSXP, REALSXP, INTSXP,
    /* si: */ REALSXP, INTSXP, LGLSXP
};

// .C(twins, ..)  called from R's  agnes() and diana():
static R_NativePrimitiveArgType twins_typ[18] = {
    INTSXP, INTSXP, REALSXP, REALSXP, REALSXP,
    /* jdiss: */ INTSXP, REALSXP,
    INTSXP, INTSXP, INTSXP, INTSXP,
    /* kwan: */ INTSXP, INTSXP, REALSXP, REALSXP,
    REALSXP, INTSXP, INTSXP
};

/* is only .C()-called from ../tests/sweep-ex.R : */
static R_NativePrimitiveArgType cl_sweep_typ[5] = {
    REALSXP, INTSXP, INTSXP, INTSXP, REALSXP
};

static const R_CMethodDef CEntries[]  = {
    CDEF(R_bncoef),
    CDEF(cl_clara),
    {"cl_daisy", (DL_FUNC) &cldaisy, 11},/* ./daisy.c */
    CDEF(cl_fanny),
    {"cl_mona", (DL_FUNC) &clmona, 9},/* ./mona.c */
    CDEF(spannel),
    CDEF(cl_sweep),
    CDEF(sildist),
    CDEF(twins),
    {NULL, NULL, 0}
};

static const R_CallMethodDef CallEntries[] = {
    CALLDEF(cl_Pam, 13),
    {NULL, NULL, 0}
};

void R_init_cluster(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, CallEntries, NULL/*FortEntries*/, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
