#include <R.h>
#include <R_ext/Lapack.h>
#include <Rinternals.h>
#include <R_ext/RS.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>


/* to get all functions:

   nm -g ./lib/spam/libs/spam.so | grep " T "

*/

/* .Fortran calls */
extern void F77_NAME(cholstepwise     )(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
// 16
extern void F77_NAME(updatefactor     )( void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
// 15
extern void F77_NAME(amub             )( void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
// 14
extern void F77_NAME(aplsb1           )( void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(closestdist      )( void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
// 13
extern void F77_NAME(aemub            )( void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(submat           )( void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(pivotforwardsolve)( void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(pivotbacksolve   )( void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(backsolves       )( void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(kroneckerf       )( void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(dn_eigen_f       )( void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(rbindf           )( void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
// 12
extern void F77_NAME(amask            )( void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(subass           )( void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(kroneckermult    )( void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(ds_eigen_f       )( void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
// 11
extern void F77_NAME(getblock         )( void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(cbindf           )( void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
// 10
extern void F77_NAME(getdia           )( void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(amubdg           )( void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(reducedim        )( void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(triplet3csr      )( void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(triplet2csr      )( void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
// 9
extern void F77_NAME(aplbdg           )( void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(getlines         )( void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(dperm            )( void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(backsolvef       )( void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(forwardsolvef    )( void *, void *, void *, void *, void *, void *, void *, void *, void *);
// 8
extern void F77_NAME(gmult_f          )( void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(spamdnscsr       )( void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(transpose        )( void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(colmeans         )( void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(amuxmat          )( void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(toeplitz         )( void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(notzero          )( void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(getallelem       )( void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(cperm            )( void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(rperm            )( void *, void *, void *, void *, void *, void *, void *, void *);
// 7
extern void F77_NAME(calcja           )( void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(spamback         )( void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(spamforward      )( void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(circulant        )( void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(setdiagmat       )( void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(diagaddmat       )( void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(getelem          )( void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(getl             )( void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(getu             )( void *, void *, void *, void *, void *, void *, void *);
// 6
extern void F77_NAME(rowmeans         )( void *, void *, void *, void *, void *, void *);
extern void F77_NAME(amux             )( void *, void *, void *, void *, void *, void *);
extern void F77_NAME(disttospam       )( void *, void *, void *, void *, void *, void *);
extern void F77_NAME(subfullsparse    )( void *, void *, void *, void *, void *, void *);
// 5
extern void F77_NAME(colsums          )( void *, void *, void *, void *, void *);
extern void F77_NAME(getdiag          )( void *, void *, void *, void *, void *);
extern void F77_NAME(spamcsrdns       )( void *, void *, void *, void *, void *);
extern void F77_NAME(addsparsefull    )( void *, void *, void *, void *, void *);
extern void F77_NAME(subsparsefull    )( void *, void *, void *, void *, void *);
extern void F77_NAME(cleanspam        )( void *, void *, void *, void *, void *);
extern void F77_NAME(getbwd           )( void *, void *, void *, void *, void *);
// 4
extern void F77_NAME(rowsums          )( void *, void *, void *, void *);
extern void F77_NAME(sortrows         )( void *, void *, void *, void *);
extern void F77_NAME(constructia      )( void *, void *, void *, void *);
extern void F77_NAME(diagmua          )( void *, void *, void *, void *);


static const R_FortranMethodDef FortranEntries[] = {
    {"cholstepwise",      (DL_FUNC) &F77_NAME(cholstepwise      ),20},
    {"updatefactor",      (DL_FUNC) &F77_NAME(updatefactor      ),16},
    {"amub",              (DL_FUNC) &F77_NAME(amub              ),15},

    {"aplsb1",            (DL_FUNC) &F77_NAME(aplsb1            ),14},
    {"closestdist",       (DL_FUNC) &F77_NAME(closestdist       ),14},

    {"pivotbacksolve",    (DL_FUNC) &F77_NAME(pivotbacksolve    ),13},
    {"pivotforwardsolve", (DL_FUNC) &F77_NAME(pivotforwardsolve ),13},
    {"backsolves",        (DL_FUNC) &F77_NAME(backsolves        ),13},
    {"submat",            (DL_FUNC) &F77_NAME(submat            ),13},
    {"aemub",             (DL_FUNC) &F77_NAME(aemub             ),13},
    {"kroneckerf",        (DL_FUNC) &F77_NAME(kroneckerf        ),13},
    {"dn_eigen_f",        (DL_FUNC) &F77_NAME(dn_eigen_f        ),13},
    {"rbindf",            (DL_FUNC) &F77_NAME(rbindf            ),13},

    {"subass",            (DL_FUNC) &F77_NAME(subass            ),12},
    {"amask",             (DL_FUNC) &F77_NAME(amask             ),12},
    {"kroneckermult",     (DL_FUNC) &F77_NAME(kroneckermult     ),12},
    {"ds_eigen_f",        (DL_FUNC) &F77_NAME(ds_eigen_f        ),12},

    {"getblock",          (DL_FUNC) &F77_NAME(getblock          ),11},
    {"cbindf",            (DL_FUNC) &F77_NAME(cbindf            ),11},

    {"getdia",            (DL_FUNC) &F77_NAME(getdia            ),10},
    {"amubdg",            (DL_FUNC) &F77_NAME(amubdg            ),10},
    {"reducedim",         (DL_FUNC) &F77_NAME(reducedim         ),10},
    {"triplet3csr",       (DL_FUNC) &F77_NAME(triplet3csr       ),10},
    {"triplet2csr",       (DL_FUNC) &F77_NAME(triplet2csr       ),10},

    {"aplbdg",            (DL_FUNC) &F77_NAME(aplbdg            ), 9},
    {"getlines",          (DL_FUNC) &F77_NAME(getlines          ), 9},
    {"dperm",             (DL_FUNC) &F77_NAME(dperm             ), 9},
    {"backsolvef",        (DL_FUNC) &F77_NAME(backsolvef        ), 9},
    {"forwardsolvef",     (DL_FUNC) &F77_NAME(forwardsolvef     ), 9},

    {"gmult_f",           (DL_FUNC) &F77_NAME(gmult_f           ), 8},
    {"spamdnscsr",        (DL_FUNC) &F77_NAME(spamdnscsr        ), 8},
    {"transpose",         (DL_FUNC) &F77_NAME(transpose         ), 8},
    {"colmeans",          (DL_FUNC) &F77_NAME(colmeans          ), 8},
    {"amuxmat",           (DL_FUNC) &F77_NAME(amuxmat           ), 8},
    {"toeplitz",          (DL_FUNC) &F77_NAME(toeplitz          ), 8},
    {"notzero",           (DL_FUNC) &F77_NAME(notzero           ), 8},
    {"getallelem",        (DL_FUNC) &F77_NAME(getallelem        ), 8},
    {"cperm",             (DL_FUNC) &F77_NAME(cperm             ), 8},
    {"rperm",             (DL_FUNC) &F77_NAME(rperm             ), 8},

    {"spamback",          (DL_FUNC) &F77_NAME(spamback          ), 7},
    {"spamforward",       (DL_FUNC) &F77_NAME(spamforward       ), 7},
    {"calcja",            (DL_FUNC) &F77_NAME(calcja            ), 7},
    {"circulant",         (DL_FUNC) &F77_NAME(circulant         ), 7},
    {"setdiagmat",        (DL_FUNC) &F77_NAME(setdiagmat        ), 7},
    {"diagaddmat",        (DL_FUNC) &F77_NAME(diagaddmat        ), 7},
    {"getelem",           (DL_FUNC) &F77_NAME(getelem           ), 7},
    {"getl",              (DL_FUNC) &F77_NAME(getl              ), 7},
    {"getu",              (DL_FUNC) &F77_NAME(getu              ), 7},

    {"rowmeans",          (DL_FUNC) &F77_NAME(rowmeans          ), 6},
    {"amux",              (DL_FUNC) &F77_NAME(amux              ), 6},
    {"disttospam",        (DL_FUNC) &F77_NAME(disttospam        ), 6},
    {"subfullsparse",     (DL_FUNC) &F77_NAME(subfullsparse     ), 6},

    {"colsums",           (DL_FUNC) &F77_NAME(colsums           ), 5},
    {"getdiag",           (DL_FUNC) &F77_NAME(getdiag           ), 5},
    {"spamcsrdns",        (DL_FUNC) &F77_NAME(spamcsrdns        ), 5},
    {"addsparsefull",     (DL_FUNC) &F77_NAME(addsparsefull     ), 5},
    {"subsparsefull",     (DL_FUNC) &F77_NAME(subsparsefull     ), 5},
    {"cleanspam",         (DL_FUNC) &F77_NAME(cleanspam         ), 5},
    {"getbwd",            (DL_FUNC) &F77_NAME(getbwd            ), 5},

    {"rowsums",           (DL_FUNC) &F77_NAME(rowsums           ), 4},
    {"sortrows",          (DL_FUNC) &F77_NAME(sortrows          ), 4},
    {"constructia",       (DL_FUNC) &F77_NAME(constructia       ), 4},
    {"diagmua",           (DL_FUNC) &F77_NAME(diagmua           ), 4},
    {NULL, NULL, 0}
};





/* .Call calls */
extern SEXP _spam_closestdistCpp(SEXP, SEXP, SEXP, SEXP , SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
  {"_spam_closestdistCpp", (DL_FUNC) &_spam_closestdistCpp, 6},
  {NULL, NULL, 0}
};



void R_init_spam(DllInfo *dll)
{
  R_registerRoutines(dll, NULL, CallEntries, FortranEntries, NULL);
  R_useDynamicSymbols(dll, FALSE);
}

