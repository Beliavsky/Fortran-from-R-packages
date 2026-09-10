/*! \file solve.cc
 *  \brief functions to solve matrix
 *
 *  \version 1
 *
 *  \date Created: 25/05/06
 *  \date Last modified: Time-stamp: <2023-01-28 15:47:35 (antoine)>
 *
 *  \author A. Lucas
 *
 *  \note Licence: GPL (>= 2)
 */

#include "bigintegerR.h"
#include "solve.h"

#include "bigrationalR.h"

using namespace std;

// inverse a rational matrix
SEXP inverse_q(SEXP A)
{
 try{
    bigvec_q a = bigrationalR::create_bignum(A);

    return(solve_gmp_R::inverse_q(a));
  } catch(std::invalid_argument & e){
    Rf_error("%s",e.what());
    return Rf_mkString(0);
  }
}

SEXP solve_gmp_R::inverse_q(bigvec_q a)
{
  if(a.nrow * a.nrow != (int) a.size()){
    a.clear();
    throw invalid_argument(_("Argument 1 must be a square matrix"));
  }
  bigvec_q b (a.size());
  b.nrow = a.nrow;

  // initialize b to identity
  for(int i=0; i<b.nrow ; ++i)
    for(int j=0; j<b.nrow ; ++j)
      b[i+j*b.nrow].setValue((i == j) ? 1 : 0);

  solve_gmp_R::solve(a,b);

  return(bigrationalR::create_SEXP(b));
}


SEXP inverse_z (SEXP A)
{
 try{
    bigvec a = bigintegerR::create_bignum(A);

    if(a.nrow * a.nrow != (int) a.size()){
      throw invalid_argument(_("Argument 1 must be a square matrix"));
    }

    if(a.getType() == TYPE_MODULUS::MODULUS_GLOBAL) {
      bigvec b (a.size() );
      b.nrow = a.nrow;
  
      // initialize b to identity
      for(int i=0; i<b.nrow ; ++i)
	for(int j=0; j<b.nrow ; ++j)
	  b[i+j*b.nrow].setValue((i == j) ? 1 : 0);

      b.setGlobalModulus(a.getGlobalModulus());
      solve_gmp_R::solve(a,b);

      return(bigintegerR::create_SEXP(b));
    }
    else {
      bigvec_q aq (a);
      return(solve_gmp_R::inverse_q(aq));
    }
  } catch(std::invalid_argument & e){
    Rf_error("%s",e.what());
    return Rf_mkString(0);
 }
}

// solve AX=B
SEXP solve_z(SEXP A,SEXP B)
{
  try
    {
      bigvec a = bigintegerR::create_bignum(A);
      bigvec b = bigintegerR::create_bignum(B);

      // case: b a vector
      if(b.nrow<1)
	b.nrow = b.size();
  
      if(a.nrow * a.nrow != (int) a.size()){
	throw invalid_argument(_("Argument 1 must be a square matrix"));
      }
	
      if(a.nrow != b.nrow){
	throw invalid_argument(_("Dimensions do not match"));
      }
  
      if(a.getType() == TYPE_MODULUS::MODULUS_GLOBAL && b.getType() != TYPE_MODULUS::MODULUS_BY_CELL ) {

	if(b.getType() == TYPE_MODULUS::NO_MODULUS){
	  b.setGlobalModulus(a.getGlobalModulus());
	}
    
	if(*a.getGlobalModulus() == *b.getGlobalModulus())
	  {
	    // solve in Z/nZ
	
	
	    solve_gmp_R::solve(a,b);
	
	    return(bigintegerR::create_SEXP(b));
	
	  }
      }
      // Solve as rational numbers.
      bigvec_q aq (a);
      bigvec_q bq (b);
      return(solve_gmp_R::solve_q(aq,bq));
    } catch(std::invalid_argument & e){
    Rf_error("%s",e.what());
    return Rf_mkString(0);
  }
}


// solve AX=B
SEXP solve_q(SEXP A,SEXP B)
{
 try{
    bigvec_q a = bigrationalR::create_bignum(A);
    bigvec_q b = bigrationalR::create_bignum(B);

    return(solve_gmp_R::solve_q(a,b));
  } catch(std::invalid_argument & e){
    Rf_error("%s",e.what());
    return Rf_mkString(0);
  }
}

// solve AX = B
SEXP solve_gmp_R::solve_q(bigvec_q a, bigvec_q b)
{
  if(a.nrow * a.nrow != (int) a.size()){
    a.clear();
    b.clear();
    throw invalid_argument(_("Argument 1 must be a square matrix"));
  }
  // case: b a vector
  if(b.nrow<0)
    b.nrow = b.size();

  if(a.nrow != b.nrow){
    a.clear();
    b.clear();
    throw invalid_argument(_("Dimensions do not match"));
  }
  
  solve_gmp_R::solve(a,b);

  return(bigrationalR::create_SEXP(b));
}
