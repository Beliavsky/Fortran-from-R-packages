/* Register routines with R ***************************************************/

#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

#include "VaR_bounds.h"


static const R_CallMethodDef callMethods[] = {
	/* {"indices_opp_ordered_to", (DL_FUNC) &indices_opp_ordered_to, 1}, */
	{"col_split", (DL_FUNC) &col_split, 1},
	{NULL, NULL, 0}
};

void R_init_qrmtools(DllInfo *dll)
{
	R_registerRoutines(dll, NULL, callMethods, NULL, NULL); /* s. WRE (2015, Section 5.4) */
	R_useDynamicSymbols(dll, FALSE);
}
