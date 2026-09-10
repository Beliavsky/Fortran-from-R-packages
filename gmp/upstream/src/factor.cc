/*! \file factor.cc
 *  \brief C function used for factorization
 *
 *  \version 1
 *
 *  \date Created: 04/12/04
 *  \date Last modified: Time-stamp: <2023-02-03 12:04:13 (antoine)>
 *
 *  \author Antoine Lucas (help from Immanuel Scholz) (R adaptation)
 *          Original C code from libgmp.
 *
 *  \note Licence: GPL (>= 2)
 */


#include <stdexcept>
#include "factorize.h"
#include "Rgmp.h"

using namespace std;


#include "factor.h"

//
// \brief function that gets values from R and send to functions
// factor
//
SEXP factorR (SEXP n)
{
  try{
    bigvec v = bigintegerR::create_bignum(n), result;
    if(v.size() > 0) {
    mpz_t val;
    mpz_init(val);
    mpz_t_sentry val_s(val);
    mpz_set(val,v[0].getValue().getValueTemp());

    int sgn = mpz_sgn(val);
    if(sgn == 0){
      v.clear();
      throw invalid_argument(_("Cannot factorize 0"));
    }
    if(sgn<0)
      {
	mpz_abs(val,val);
	result.push_back(bigmod(biginteger(-1)));
      }
    //
    // function from gmplib, in demo/factorize.c
    //
    factor(val,result);
  }
  return bigintegerR::create_SEXP(result);
  }    catch(std::invalid_argument & e){
    Rf_error("%s",e.what());
    return Rf_mkString(0);
  }
}
