#include <R.h>
#include <Rdefines.h>
/* global variables */
extern SEXP Time, Y ;
extern SEXP st_gparms;
extern SEXP st_gforcs;

extern int ndim, nspec;
extern double *y2,*dy2;

/* steady_utils.c globals */
extern SEXP Rst_deriv_func;
extern SEXP Rst_jac_func;
extern SEXP Rst_envir;

extern SEXP stsparse_deriv_func;
extern SEXP stsparse_jac_func;
extern SEXP stsparse_envir;

/* runsteady globals*/
extern SEXP lsode_deriv_func;
extern SEXP lsode_jac_func;
extern SEXP lsode_envir;


extern void Initstparms(int *, double *);
extern void Initstforcs(int *, double *);
extern void initOut(int isDll, int neq, SEXP nOut, SEXP Rpar, SEXP Ipar);

/* output in DLL globals */
extern int nout, ntot, isOut, lrpar, lipar, *ipar;
extern double *out;

typedef void C_deriv_func_type(int *, double *, double *,double *,double *, int *);
typedef void C_init_func_type(void (*)(int *, double *));
