// This ensures registration -- see also useDynLib(...) in  ../NAMESPACE

#include "bigintegerR.h"
#include <R.h>
#include <Rinternals.h>

// include those that have an  extern "C" { .... } :
#include "apply.h"

#include "bigrationalR.h"
//
#include "extract_matrix.h"
#include "factor.h"
#include "matrix.h"
#include "matrixq.h"
#include "solve.h"

#include <R_ext/Rdynload.h>

#define CALLDEF(name, n)  {#name, (DL_FUNC) &name, n}

static R_CallMethodDef CallEntries[] = {
// apply.h :
  CALLDEF(gmpMatToListZ, 2),
  CALLDEF(gmpMatToListQ, 2),
// bigintegerR.h :
  CALLDEF(R_gmp_get_version, 0),
  CALLDEF(biginteger_add, 2),
  CALLDEF(biginteger_sub, 2),
  CALLDEF(biginteger_mul, 2),
  CALLDEF(biginteger_div, 2),
  CALLDEF(biginteger_divq, 2),
  CALLDEF(biginteger_mod, 2),
  CALLDEF(biginteger_pow, 2),
  CALLDEF(biginteger_inv, 2),
  CALLDEF(biginteger_gcd, 2),
  CALLDEF(biginteger_lcm, 2),
  // CALLDEF(biginteger_setmod, 2),
  CALLDEF(biginteger_get_at, 2),
  CALLDEF(biginteger_set_at, 3),
  CALLDEF(biginteger_as, 2),
  CALLDEF(biginteger_as_character, 2),
  CALLDEF(biginteger_as_numeric, 1),
  CALLDEF(biginteger_as_integer, 1),
  CALLDEF(biginteger_length, 1),
  CALLDEF(biginteger_setlength, 2),
  CALLDEF(biginteger_is_na, 1),
  CALLDEF(biginteger_sgn, 1),
  CALLDEF(biginteger_lt, 2),
  CALLDEF(biginteger_gt, 2),
  CALLDEF(biginteger_lte, 2),
  CALLDEF(biginteger_gte, 2),
  CALLDEF(biginteger_eq, 2),
  CALLDEF(biginteger_neq, 2),
  CALLDEF(biginteger_c, 1),
  CALLDEF(biginteger_cbind, 1),
  CALLDEF(biginteger_rep, 2),
  CALLDEF(biginteger_is_prime, 2),
  CALLDEF(biginteger_nextprime, 1),
  CALLDEF(biginteger_abs, 1),
  CALLDEF(biginteger_gcdex, 2),
  CALLDEF(biginteger_rand_u, 4),
  CALLDEF(biginteger_sizeinbase, 2),
  CALLDEF(bigI_frexp, 1),
  CALLDEF(bigI_choose, 2),
  CALLDEF(bigI_factorial, 1),
  CALLDEF(bigI_fibnum, 1),
  CALLDEF(bigI_fibnum2, 1),
  CALLDEF(bigI_lucnum, 1),
  CALLDEF(bigI_lucnum2, 1),
  CALLDEF(biginteger_max, 2),
  CALLDEF(biginteger_min, 2),
  CALLDEF(biginteger_cumsum, 1),
  CALLDEF(biginteger_sum, 1),
  CALLDEF(biginteger_prod, 1),
  CALLDEF(biginteger_powm, 3),
  CALLDEF(biginteger_log2, 1),
  CALLDEF(biginteger_log, 1),

// bigrationalR.h:
  CALLDEF(bigrational_add, 2),
  CALLDEF(bigrational_sub, 2),
  CALLDEF(bigrational_mul, 2),
  CALLDEF(bigrational_div, 2),
  CALLDEF(bigrational_pow, 2),
  CALLDEF(bigrational_num, 1),
  CALLDEF(bigrational_den, 1),
  CALLDEF(bigrational_get_at, 2),
  CALLDEF(bigrational_set_at, 3),
  CALLDEF(bigrational_as, 2),
  CALLDEF(bigrational_as_character, 2),
  CALLDEF(bigrational_as_numeric, 1),
  CALLDEF(bigrational_length, 1),
  CALLDEF(bigrational_setlength, 2),
  CALLDEF(bigrational_is_na, 1),
  CALLDEF(bigrational_is_int, 1),
  CALLDEF(bigrational_lt, 2),
  CALLDEF(bigrational_gt, 2),
  CALLDEF(bigrational_lte, 2),
  CALLDEF(bigrational_gte, 2),
  CALLDEF(bigrational_eq, 2),
  CALLDEF(bigrational_neq, 2),
  CALLDEF(bigrational_c, 1),
  CALLDEF(bigrational_cbind, 1),
  CALLDEF(bigrational_rep, 2),
  CALLDEF(bigrational_max, 2),
  CALLDEF(bigrational_min, 2),
  CALLDEF(bigrational_cumsum, 1),
  CALLDEF(bigrational_sum, 1),
  CALLDEF(bigrational_prod, 1),

// extract_matrix.h :
  CALLDEF(matrix_get_at_z, 3),
  CALLDEF(matrix_set_at_z, 4),
  CALLDEF(matrix_get_at_q, 3),
  CALLDEF(matrix_set_at_q, 4),

// factor.h :
  CALLDEF(factorR, 1),

// matrix.h :
  CALLDEF(is_matrix_zq, 1),
  CALLDEF(as_matrixz, 5),
  CALLDEF(bigint_transposeR, 1),
  CALLDEF(matrix_crossp_z, 2),
  CALLDEF(matrix_mul_z, 3),
  CALLDEF(biginteger_rbind, 1),

// matrixq.h :
  CALLDEF(as_matrixq, 5),
  CALLDEF(bigq_transposeR, 1),
  CALLDEF(matrix_crossp_q, 2),
  CALLDEF(matrix_mul_q, 3),
  CALLDEF(bigrational_rbind, 1),

// solve.h :
  CALLDEF(inverse_q, 1),
  CALLDEF(solve_q, 2),
  CALLDEF(inverse_z, 1),
  CALLDEF(solve_z, 2),

    {NULL, NULL, 0}
};

extern "C"
void R_init_gmp(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, (Rboolean)FALSE);
}

