/*
** $Id $
**
** Prototypes of all the survival functions
**  Including this in each routine helps prevent mismatched argument errors
*/
#include "R.h"
#include "Rinternals.h"
#define ALLOC(a,b)  R_alloc(a,b)

double **dmatrix(double *array, int ncol, int nrow);


SEXP expc(SEXP   efac2,   SEXP edims2,
		  	      SEXP   ecut2,     SEXP   expect2,
		  	      SEXP   x2, 	SEXP   y2);


SEXP netwei(   SEXP   efac2,   SEXP edims2,
	      SEXP   ecut2,     SEXP   expect2,
	      SEXP   x2, 	SEXP   y2,  SEXP status2,    SEXP times2) ;

SEXP netfastpinter(   SEXP   efac2,   SEXP edims2,
	      SEXP   ecut2,     SEXP   expect2,
	      SEXP   x2, 	SEXP   y2, SEXP ys2,  SEXP status2,    SEXP times2) ;

SEXP cmpfast(   SEXP   efac2,   SEXP edims2,
	      SEXP   ecut2,     SEXP   expect2,
	      SEXP   x2, 	SEXP   y2, SEXP ys2,  SEXP status2,    SEXP times2) ;

SEXP netfastpinter2(   SEXP   efac2,   SEXP edims2,
	      SEXP   ecut2,     SEXP   expect2,
	      SEXP   x2, 	SEXP   y2,SEXP ys2, SEXP status2,     SEXP times2, SEXP myprec2) ;


double pystep(int nc,        int  *index,  int  *index2,   double *wt,
	      double *data,  int *fac,    int *dims,     double **cuts,
	      double step,   int  edge);

double pystep2(int nc,        int  *index,  int  *index2,   double *wt,
	      double *data,  int *fac,    int *dims,     double **cuts,
	      double step,   int  edge);
