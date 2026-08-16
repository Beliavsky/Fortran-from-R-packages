
/*---------------------------------------------------------------------------*/
/* header files */

#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

#include "argus.h"
     
/*---------------------------------------------------------------------------*/

static const R_CallMethodDef CallEntries[] = {
    {"rargusRoU", (DL_FUNC) &rargusRoU, 2},
    {NULL, NULL, 0}
};

/*-------------------------------------------------------------------------*/

void 
R_init_argus (DllInfo *dll  ATTRIBUTE__UNUSED) 
/*---------------------------------------------------------------------------*/
/* Initialization routine when loading the DLL                               */
/*                                                                           */
/* Parameters:                                                               */
/*   dll ... passed by R and describes the DLL                               */
/*                                                                           */
/* Return:                                                                   */
/*   (void)                                                                  */
/*---------------------------------------------------------------------------*/
{
  /* Declare some C routines to be callable from other packages */ 
  R_RegisterCCallable("argus", "rargusRoU", (DL_FUNC) rargusRoU);
/*  R_RegisterCCallable("argus", "do_rgig", (DL_FUNC) do_rargusRoU); */

  /* Register native routines */
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}

/*---------------------------------------------------------------------------*/
