#include "extract_matrix.h"

#include "bigrationalR.h"
#include <stdexcept>
using namespace std;

// for something like x = A[indi, indj], but also simply  A[ind]
SEXP matrix_get_at_q(SEXP A,SEXP INDI, SEXP INDJ)
{
  try{
    bigvec_q mat = bigrationalR::create_bignum(A);
    
    return(bigrationalR::create_SEXP(extract_gmp_R::get_at( mat,INDI,INDJ)));
  }
  catch(std::invalid_argument & e){
    Rf_error("%s",e.what());
    return Rf_mkString(0);
  }
}

// for something like x = A[indi, indj], but also simply  A[ind]
SEXP matrix_get_at_z(SEXP A,SEXP INDI, SEXP INDJ)
{
 //  printf("ici\n");
  try
    {
      bigvec mat = bigintegerR::create_bignum(A);
      bigvec mat2 = extract_gmp_R::get_at( mat,INDI,INDJ);
  /*
  // now modulus !
  // cell based modulus
  if(mat.getType() == TYPE_MODULUS:: MODULUS_BY_CELL)
    {
      for(unsigned int i = 0; i< mat.size(); ++i)
	mat.value[i] = mat.modulus[i];

      mat = extract_gmp_R::get_at( mat,INDI,INDJ);

      mat2.modulus.resize(mat.size());
      for(unsigned int i = 0; i< mat.size(); ++i)
	mat2.modulus[i] = mat.value[i];
    }
  // row base modulus
  else if((int)mat.modulus.size() == mat.nrow)
    {
      for(unsigned int i = 0; i< mat.size(); ++i)
	mat.value[i] = mat.modulus[i];

      mat.modulus.clear();

      mat = bigintegerR::biginteger_get_at_C(mat,INDI);

      mat2.modulus.resize(mat.size());
      for(unsigned int i = 0; i< mat.size(); ++i)
	mat2.modulus[i] = mat.value[i];

    }
  //global modulus
  else if(mat.modulus.size() == 1)
    {
      mat2.modulus.resize(1);
      mat2.modulus[0] = mat.modulus[0];
      }*/
      return(bigintegerR::create_SEXP(mat2) );
    }
  catch(std::invalid_argument & e){
    Rf_error("%s",e.what());
    return Rf_mkString(0);
 }
}



// for something like A[indi, indj] <- val
SEXP matrix_set_at_z(SEXP A, SEXP VAL, SEXP INDI, SEXP INDJ)
{
  try{
    bigvec mat = bigintegerR::create_bignum(A);
    bigvec val = bigintegerR::create_bignum(VAL);
    
    extract_gmp_R::set_at( mat,val,INDI,INDJ);

    return(bigintegerR::create_SEXP(mat));
  }
  catch(std::invalid_argument & e){
    Rf_error("%s",e.what());
    return Rf_mkString(0);
 }



}

// for something like A[indi, indj] <- val
SEXP matrix_set_at_q(SEXP A,SEXP VAL ,SEXP INDI, SEXP INDJ)
{
  try{
    bigvec_q mat = bigrationalR::create_bignum(A);
    bigvec_q val = bigrationalR::create_bignum(VAL);
    
    
    extract_gmp_R::set_at( mat,val,INDI,INDJ);
    return(bigrationalR::create_SEXP(mat));
  }
  catch(std::invalid_argument & e){
    Rf_error("%s",e.what());
    return Rf_mkString(0);
 }



}



//
// return a vector of integers corresponding to values that must be affected.
//
std::vector<int> extract_gmp_R::indice_get_at (unsigned int n , SEXP & IND)
{


  std::vector<int> vidx = bigintegerR::create_int(IND);
  std::vector<int> result;

  
  if(TYPEOF(IND) == NILSXP){
    //LOCICAL: return true
    for (unsigned int i = 0;  i< n ; i++){
      result.push_back(i);
    }
  }
  else if (TYPEOF(IND) == LGLSXP)
    {
      // boolean
      for(unsigned int i = 0; i< n; ++i)
	if (vidx[i % vidx.size() ] ) result.push_back(i);
    }
  else
    //INTEGERS
    {
      vidx.erase(std::remove(vidx.begin(), vidx.end(), 0L), vidx.end()); // remove all zeroes
      if(vidx.size() == 0) return result;
      //negatives integers: all except indices will be modified
      if (vidx[0] < 0)
	{
	  std::vector<bool> tempo(n,true);

	  for (std::vector<int>::const_iterator jt = vidx.begin(); jt != vidx.end(); ++jt)
	    {
	      if(*jt > 0)
		throw invalid_argument("only 0's may mix with negative subscripts");
		  //error(_("only 0's may mix with negative subscripts"));
	      if( (*jt != 0) && (*jt >= - static_cast<int>(n)))
		tempo[-(*jt)-1] = false;
	    }
	  for (unsigned int i =0 ; i < n ; i++){
	    if(tempo[i]) result.push_back(i);
	  }
	}
      else
	{
	  //INTEGERS (and positive)
	  for (std::vector<int>::const_iterator jt = vidx.begin(); jt != vidx.end(); ++jt) {
	    int i = *jt;
	    if(i < 0)
	      throw invalid_argument("only 0's may mix with negative subscripts");
	    //error(_("only 0's may mix with negative subscripts"));
	    result.push_back(i-1);
	  }

	}
    }

  return(result);

}//end of indice_set_at

