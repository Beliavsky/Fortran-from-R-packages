#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME:
  Check these declarations against the C/Fortran source code.
*/

  /* .C calls */
extern void estdep(void *, void *, void *, void *, void *, void *);
extern void estdep_serial(void *, void *, void *, void *, void *, void *);
extern void Sn_serial0(void *, void *, void *, void *);
extern void stats_nonserial(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void stats_serial(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void stats_serialVectors(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void stats_serial_bin(void *, void *, void *, void *, void *, void *, void *, void *);
extern void statsim(void *, void *, void *, void *, void *, void *, void *);
extern void prepare_data(void *, void *, void *, void *, void *, void *);
extern void Stat_A_serial(void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void Stat_A(void *, void *, void *, void *, void *, void *, void *, void *, void *);

static const R_CMethodDef CEntries[] = {
  {"estdep",           (DL_FUNC) &estdep,            6},
  {"estdep_serial",    (DL_FUNC) &estdep_serial,     6},
  {"Sn_serial0",       (DL_FUNC) &Sn_serial0,        4},
  {"stats_nonserial",  (DL_FUNC) &stats_nonserial,  10},
  {"stats_serial",     (DL_FUNC) &stats_serial,     10},
  {"stats_serialVectors",(DL_FUNC) &stats_serialVectors,     11},
  {"stats_serial_bin", (DL_FUNC) &stats_serial_bin,  8},
  {"statsim",          (DL_FUNC) &statsim,           7},
  {"prepare_data",     (DL_FUNC) &prepare_data,      6},
  {"Stat_A_serial",    (DL_FUNC) &Stat_A_serial,    9},
  {"Stat_A",           (DL_FUNC) &Stat_A,    9},
  {NULL, NULL, 0}
};

void R_init_MixedIndTests(DllInfo *dll)
{
  R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
