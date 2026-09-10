#ifndef FACTORIZE_GMP_R
#define FACTORIZE_GMP_R 1
#include "bigvec.h"


/**
 * Factorize a prime number
 * t: number to factorize
 * factors [out]: the list of factors
 *
 * Note: this is adapted from demo "factorize.c" file from gmplib
 */
void factor (mpz_t t, bigvec & factors);


#endif
