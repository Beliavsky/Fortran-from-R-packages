#include <stdexcept>
#include "bigintegerR.h"
#include <R.h>
#include <Rinternals.h>
// and one thing from Rdefines.h :
#define NEW_LIST(n) Rf_allocVector(VECSXP,n)

#include "apply.h"

#include "bigrationalR.h"


// X a matrix, or a bigz integer
// line: true = we return a list of all lines
//
SEXP gmpMatToListZ(SEXP X, SEXP line)
{
  SEXP ans;
  // dangerous... no check => use with care in R
  int lines =  INTEGER(line)[0];
  
  bigvec matrix = bigintegerR::create_bignum(X);
  
  unsigned int ncol = matrix.size() / matrix.nrow;
  unsigned int nrow =  matrix.nrow;
  
  if(lines == 1)
    {
      // RETURN a list of all lines
      PROTECT (ans = NEW_LIST(matrix.nrow) );
      try{
	for(unsigned int i = 0; i < nrow; ++i)
	  {
	    bigvec oneLine ;
	    for(unsigned int j = 0; j < ncol; ++j)
	      {
		oneLine.push_back(matrix[i+j*nrow]);
	      }
	    
	    SET_VECTOR_ELT(ans, i,bigintegerR::create_SEXP(oneLine));
	    
	  }
      } catch(std::invalid_argument & e){
	perror(e.what());
      }
      UNPROTECT(1);
    }
  else
    {
      // RETURN a list of all rows !
      PROTECT (ans = NEW_LIST(ncol) );
      try{
	for(unsigned int j = 0; j < ncol; ++j)
	  {
	    bigvec oneLine ;
	    for(unsigned int i = 0; i < nrow; ++i)
	      {
		oneLine.push_back(matrix[i+j*nrow]);
		
	      }
	    
	    SET_VECTOR_ELT(ans, j,bigintegerR::create_SEXP(oneLine));
	    
	  }
      } catch(std::invalid_argument & e){
	perror(e.what());
      }
      
      UNPROTECT(1);
      
    }
   
  return(ans);  
}


// X a matrix, or a bigq rational
// line: true = we return a list of all lines
//
SEXP gmpMatToListQ(SEXP X, SEXP line)
{
  SEXP ans;
 
  // dangerous... no check => use with care in R
  int lines =  INTEGER(line)[0];

  bigvec_q matrix = bigrationalR::create_bignum(X);
  
  unsigned int ncol = matrix.size() / matrix.nrow;
  unsigned int nrow =  matrix.nrow;
  
  if(lines == 1)
    {
      // RETURN a list of all lines
      PROTECT (ans = NEW_LIST(matrix.nrow) );
      try
	{
	  for(unsigned int i = 0; i < nrow; ++i)
	    {
	      bigvec_q oneLine ;
	      for(unsigned int j = 0; j < ncol; ++j)
		{
		  oneLine.value.push_back(matrix.value[i+j*nrow]);
		}
	      SET_VECTOR_ELT(ans, i,bigrationalR::create_SEXP(oneLine));
	    }
	}
      catch(std::invalid_argument & e){
	perror(e.what());
      }
      UNPROTECT(1);
    }
  else
    {
      // RETURN a list of all rows !
      PROTECT (ans = NEW_LIST(ncol) );
      try
	{
	  for(unsigned int j = 0; j < ncol; ++j)
	    {
	      bigvec_q oneLine ;
	      for(unsigned int i = 0; i < nrow; ++i)
		{
		  oneLine.value.push_back(matrix.value[i+j*nrow]);
		}

	      SET_VECTOR_ELT(ans, j,bigrationalR::create_SEXP(oneLine));
	    }
	}
      catch(std::invalid_argument & e){
	perror(e.what());
      }
      UNPROTECT(1);

    }

  
     
  return(ans);
 
}


