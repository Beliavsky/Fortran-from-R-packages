/************************************************************/
/*! \file bigintegerR.cc
 *  \brief C function to interface R and libgmp with big integer values
 *
 *  \version 1
 *
 *  \date Created: 27/10/04
 *  \date Last modified: $Id: bigintegerR.cc,v 1.37 2014-09-16 07:33:43 mmaechler Exp $
 *
 *  \author Immanuel Scholz (help from A. Lucas)
 *
 *  \note Licence: GPL (>= 2)
 */

#include <stdexcept>
#include <vector>
#include <algorithm>
#include "bigintegerR.h"
#include "Rgmp.h"


// only for  <bigz> ^ -|n| :
#include "bigrationalR.h"


#include "matrix.h"




using namespace std;

#include <Rmath.h>
#include "extract_matrix.h"
/* Globals variables */

static gmp_randstate_t seed_state;
static int seed_init=0;
namespace bigintegerR
{
  // \brief create a vector of bigvecs, all without a modulus.
  bigvec create_vector(const SEXP & param) {
    lockSexp lock (param);
    switch (TYPEOF(param)) {
    case NILSXP:
      return bigvec(); // = bigz(0)
    case RAWSXP:
      {
	// deserialise the vector. first int is the size.
	bigvec v;
	const char* raw = (char*)RAW(param);
	int pos = sizeof(int); // position in raw[]. Starting after header.
	int sizevec = ((int*)raw)[0];
	//std::cout << "nb element a lire " << sizevec << std::endl;
	v.resize(sizevec);
	for (int i = 0; i < sizevec; ++i) {
	  v[i].setValue(make_shared<biginteger>(&raw[pos]));
	  pos += v[i].getValue().raw_size(); // increment number of bytes read.
	}
	return v;
      }
    case REALSXP:
      {
	double* d = REAL(param);
	//bigvec v(d,d+LENGTH(param));
	bigvec v;
	v.resize(LENGTH(param));
	for (int j = 0; j < LENGTH(param); ++j) {
	  /// New:   numeric '+- Inf'  give  +- "Large" instead of NA
	  double dj = d[j];
	  if(R_FINITE(dj) || ISNAN(dj))
	    v[j].setValue(make_shared<biginteger>( dj ));
	  else { // dj is +- Inf : use LARGE ( =   +- 2 ^ 80000 -- arbitrarily )
	    mpz_t LARGE;
	    mpz_init(LARGE);
	    // FIXME: Keep 'LARGE' a static const; initialized only once

	    if(dj == R_PosInf){
	      mpz_ui_pow_ui (LARGE, (unsigned long int) 2, (unsigned long int) 8000);
	      v[j].setValue(make_shared<biginteger>(LARGE));
	    }
	    else if(dj == R_NegInf) {
	      mpz_t neg_L;
	      mpz_init(neg_L);
	      mpz_neg(neg_L, LARGE);
	      v[j].setValue( make_shared<biginteger>(neg_L));
	      mpz_clear(neg_L);
	    }
	    else// should never happen
	      v[j].setValue( make_shared<biginteger>(dj));
	      
	    mpz_clear(LARGE);

	  }
	}
	return v;
      }
    case INTSXP:
    case LGLSXP:
      {
	int* i = INTEGER(param);
	//bigvec v(i,i+LENGTH(param));
	bigvec v;
	v.resize(LENGTH(param));

	for (int j = 0; j < LENGTH(param); ++j)
	  v[j].setValue( make_shared<biginteger>(i[j]));

	return v;
      }
    case STRSXP:
      {
	bigvec v;
	v.resize(LENGTH(param));
	for (int i = 0; i < LENGTH(param); ++i) {
	  if (STRING_ELT(param,i) != NA_STRING)
	    v[i].setValue(make_shared<biginteger>(std::string(CHAR(STRING_ELT(param,i)))));
	}
	return v;
      }
    default:
      // no longer: can be fatal later! return bigvec();
      throw invalid_argument(_("only logical, numeric or character (atomic) vectors can be coerced to 'bigz'"));
    }
  }

  bigvec create_bignum(const SEXP & param) {
    lockSexp lock (param);
    SEXP
      modAttr = Rf_getAttrib(param, Rf_mkString("mod")),
      dimAttr = Rf_getAttrib(param, Rf_mkString("nrow"));

    // try to catch biz-nrow dimension value
    //std::cout << "import value" << std::endl;
    bigvec v = bigintegerR::create_vector(param);

    if (TYPEOF(dimAttr) == INTSXP){
      v.nrow = INTEGER(dimAttr)[0];
    }
    else {
      // catch to get std matrix dimensions value
      dimAttr = Rf_getAttrib(param, Rf_mkString("dim"));
      v.nrow = (TYPEOF(dimAttr) == INTSXP) ? INTEGER(dimAttr)[0] : -1;// -1: want support 0-row
    }

    if (TYPEOF(modAttr) != NILSXP)
      {
	//std::cout << "import value" << std::endl;
	bigvec vMod = bigintegerR::create_vector(modAttr);
	for (unsigned int i = 0 ; i < v.size(); i++){
	  v[i].setModulus(vMod[i % vMod.size()].getValuePtr());
	}
	v.setType(vMod.size() == 1 ? TYPE_MODULUS::MODULUS_GLOBAL :  TYPE_MODULUS::MODULUS_BY_CELL);
      }

    return v;

  }

  std::vector<int> create_int(const SEXP & param) {
    lockSexp lock (param);
    switch (TYPEOF(param)) {
    case REALSXP:
      {
	double* d = REAL(param);
	// copy vector manually to avoid stupid conversion warning in STL-code :-/
	vector<int> v;
	for (int i = 0; i < LENGTH(param); ++i)
	  v.push_back(static_cast<int>(d[i]));
	return v;
	//return vector<int>(d, d+LENGTH(param));
      }
    case INTSXP:
    case LGLSXP:
      {
	int* i = INTEGER(param);

	return std::vector<int>(i, i+LENGTH(param));
      }
    default:
      return std::vector<int>();
    }
  }


  SEXP create_SEXP(const bigvec & v , bigmod_accessor_fn fct, unsigned int size)
  {
    unsigned int i;
    int sizeRaw = sizeof(int); // starting with vector-size-header
    for (i = 0; i < size; ++i)
      sizeRaw += fct(v[i]).raw_size(); // adding each bigint's needed size
    SEXP ans = PROTECT(Rf_allocVector(RAWSXP, sizeRaw));
    // Rprintf("    o create_SEXP(vect<biginteger>): size=%d, v.size()=%d\n", size, v.size());
    char* r = (char*)(RAW(ans));
    ((int*)(r))[0] = size; // first int is vector-size-header
    int pos = sizeof(int); // current position in r[] (starting after vector-size-header)
    for (i = 0; i < size; ++i)
      pos += fct(v[i]).as_raw(&r[pos]);
    UNPROTECT(1);
    return(ans);
  }


  const biginteger & bigModToValue(const bigmod& v){
    return v.getValue();
  }

  const biginteger & bigModToModulus(const bigmod& v){
    return v.getModulus();
  }

  SEXP create_SEXP(const bigvec& v) {
    int size = v.size();
    SEXP ans = PROTECT(create_SEXP(v, bigModToValue, size));
    // set the class attribute to "bigz"
    Rf_setAttrib(ans, R_ClassSymbol, Rf_mkString("bigz"));

    // Rprintf("   o create_SEXP(<bigvec>): v.nrow=%d", v.nrow);
    // set the dim attribute
    if(v.nrow >= 0)  {
      SEXP nrowAttr  = Rf_mkString("nrow");
      PROTECT(nrowAttr);
      SEXP nrowValue = Rf_ScalarInteger((int) v.nrow);
      PROTECT(nrowValue);
      Rf_setAttrib(ans, nrowAttr,nrowValue);
      UNPROTECT(2);
    }
    // set the mod attribute
    if(v.getType() !=  TYPE_MODULUS::NO_MODULUS && v.size() > 0 ) {
      int sizeM = v.getType() ==  TYPE_MODULUS::MODULUS_GLOBAL ? 1 : size;
      SEXP mod = PROTECT(create_SEXP(v, bigModToModulus, sizeM )); // and set *its* class
      Rf_setAttrib(mod, R_ClassSymbol, Rf_mkString("bigz"));
      Rf_setAttrib(ans, Rf_mkString("mod"), mod);
      UNPROTECT(1);
    }
    UNPROTECT(1);
    return ans;
  }

  SEXP biginteger_binary_operation(const bigvec& va,const bigvec& vb, biginteger_binary_fn f)
  {
    
    int size = (va.size() ==0 || vb.size()==0 ) ? 0 : max(va.size(), vb.size());
 
    int nrow = matrixz::checkDims(va.nrow,vb.nrow);
    if(nrow == -2){
      throw invalid_argument(_("Matrix dimensions do not match"));
    }
    bigvec result;
    
    for (int i = 0; i < size; ++i){
      result.push_back(f(va[i%va.size()], vb[i%vb.size()]));
    }
    result.nrow = nrow;
	
    // Rprintf(" o bigI_b_op(.); size=%d -> nrow=%d\n", size, result.nrow);
    return bigintegerR::create_SEXP(result);
    
  }



  /**
   * \brief Main function of doing a binary operation on bigintegers.
   * It calls a function argument for doing the correct thing.
   * This could also be written as a class functor (template)
   * to save one function call, but then code bloat will happen.
   */
  SEXP biginteger_binary_operation(const SEXP& a,const SEXP& b, biginteger_binary_fn f)
  {
    try
      {
	bigvec va = bigintegerR::create_bignum(a);
	bigvec vb = bigintegerR::create_bignum(b);
	return biginteger_binary_operation(va, vb, f);
      } catch(std::invalid_argument & e){
      Rf_error("%s",e.what());
      return Rf_mkString(0);
    }
  }


  SEXP biginteger_logical_binary_operation(const SEXP & a,const SEXP & b, biginteger_logical_binary_fn f)
  {
    try{
      bigvec va = bigintegerR::create_bignum(a);
      bigvec vb = bigintegerR::create_bignum(b), result;

      int nrow = matrixz::checkDims(va.nrow,vb.nrow) ;
      if(nrow == -2){
	va.clear();
	vb.clear();
	throw invalid_argument (_("Matrix dimensions do not match"));
      }
    
      int size = (va.size()==0 || vb.size() == 0) ? 0 : max(va.size(), vb.size());
      //	int sizemod = max(va.modulus.size(), vb.modulus.size());
      SEXP ans = PROTECT(Rf_allocVector(LGLSXP, size));
      int *r = LOGICAL(ans);
      /* TODO: this kind of situation 5 == (5 %% 17)*/
      for (int i = 0; i < size; ++i) {
	biginteger & am = va[i % va.size()].getValue();
	biginteger & bm = vb[i % vb.size()].getValue();
	if (am.isNA() || bm.isNA())
	  r[i] = NA_LOGICAL;
	else
	  r[i] = f(am, bm) ? 1 : 0;
      }

      // Add dimension parameter when available
      if(nrow >= 0)
	{
	  SEXP dimVal;
	  PROTECT(dimVal = Rf_allocVector(INTSXP, 2));
	  INTEGER(dimVal)[0] = (int) nrow;
	  INTEGER(dimVal)[1] = (int) size / nrow;
	  Rf_setAttrib(ans, Rf_mkString("dim"), dimVal);
	  UNPROTECT(1);
	}

      UNPROTECT(1);
      return ans;
    }catch(std::invalid_argument & e){
      Rf_error("%s",e.what());
      return Rf_mkString(0);
    }
    
  }

  bool lt(const biginteger& lhs, const biginteger& rhs)
  {return mpz_cmp(lhs.getValueTemp(), rhs.getValueTemp()) < 0;}
  bool gt(const biginteger& lhs, const biginteger& rhs)
  {return mpz_cmp(lhs.getValueTemp(), rhs.getValueTemp()) > 0;}
  bool lte(const biginteger& lhs, const biginteger& rhs)
  {return mpz_cmp(lhs.getValueTemp(), rhs.getValueTemp()) <= 0;}
  bool gte(const biginteger& lhs, const biginteger& rhs)
  {return mpz_cmp(lhs.getValueTemp(), rhs.getValueTemp()) >= 0;}
  bool eq(const biginteger& lhs, const biginteger& rhs)
  {return mpz_cmp(lhs.getValueTemp(), rhs.getValueTemp()) == 0;}
  bool neq(const biginteger& lhs, const biginteger& rhs)
  {return mpz_cmp(lhs.getValueTemp(), rhs.getValueTemp()) != 0;}
  
}

/* End of namespace bigintegerR*/


SEXP R_gmp_get_version() {
  return Rf_mkString(gmp_version);
}




SEXP biginteger_add (SEXP a, SEXP b) {return bigintegerR::biginteger_binary_operation(a,b,operator+);}
SEXP biginteger_sub (SEXP a, SEXP b) {return bigintegerR::biginteger_binary_operation(a,b,operator-);}
SEXP biginteger_mul (SEXP a, SEXP b) {return bigintegerR::biginteger_binary_operation(a,b,operator*);}
SEXP biginteger_divq(SEXP a, SEXP b) {return bigintegerR::biginteger_binary_operation(a,b,operator/);}
SEXP biginteger_mod (SEXP a, SEXP b) {return bigintegerR::biginteger_binary_operation(a,b,operator%);}

SEXP biginteger_div (SEXP a, SEXP b) { // called from  "/.bigz" == div.bigz
  try{
    bigvec A = bigintegerR::create_bignum(a),
      B = bigintegerR::create_bignum(b);
    // Note: a or b may be simple numbers (e.g. '1') => work with (A,B)
    TYPE_MODULUS type_a = A.getType();
    TYPE_MODULUS type_b = B.getType();
    
    if(type_a == TYPE_MODULUS::NO_MODULUS && type_b == TYPE_MODULUS::NO_MODULUS) // deal with important case quickly:
      {
	return bigrationalR::bigrational_binary_operation(bigrationalR::create_bignum(a),
							  bigrationalR::create_bignum(b),
							  operator/);
      }
    else if(type_a == TYPE_MODULUS::NO_MODULUS) { // and  len_m_b > 0:
      // should work directly using b's "mod" --> compute  a * b^(-1)
    }
    else if(type_b == TYPE_MODULUS::NO_MODULUS) { // and  len_m_a > 0:
      // should use a's "mod" for b: div_via_inv() need's  b's modulus
      if(type_a == TYPE_MODULUS::MODULUS_GLOBAL){
	B.setGlobalModulus(A.getGlobalModulus());
      } else{
	for (unsigned int i = 0 ; i < B.size(); i++){
	  B[i].setModulus(A[i % A.size()].getModulusPtr());
	}
      }

      return bigintegerR::biginteger_binary_operation(A,
						      B,
						      div_via_inv);
    }
    else { // len_m_a > 0 and  len_m_b > 0:
      bool same_mod = true;// are the two mods the "same" (after recycling)?
      unsigned int len_m_a = A.getModulusSize();
      unsigned int len_m_b = B.getModulusSize();
      unsigned int m = (len_m_a < len_m_b) ? len_m_b : len_m_a; // = max(l..a, l..b)
      for(unsigned int i = 0; i < m; i++)
	if(A[i % len_m_a].getModulus() != B[i % len_m_b].getModulus())  {
	  same_mod = false; break;
	}
      if(same_mod) {
	// compute   a * b^(-1) ... should work w/o more
      } else {
	// use *rational* a/b  (not considering 'mod' anymore):
	A.clear();
	B.clear();
	return  bigrationalR::bigrational_binary_operation(bigrationalR::create_bignum(a),
							   bigrationalR::create_bignum(b),
							   operator/);
      }
    }

    return bigintegerR::biginteger_binary_operation(A,B, div_via_inv);
  }catch(std::invalid_argument & e){
    Rf_error("%s",e.what());
    return Rf_mkString(0);
  }
}

SEXP biginteger_pow (SEXP a, SEXP b) {
  try{
    bigvec v = bigintegerR::create_bignum(a),
      exp = bigintegerR::create_bignum(b);
    if(v.getType() == TYPE_MODULUS::NO_MODULUS) { /* has no modulus: now, if any b < 0, the
						     result must be (non-integer) bigrational */
      bool use_rat = FALSE;
      for (unsigned int i = 0; i < exp.size(); ++i) {
	if(mpz_sgn(exp[i].getValue().getValueTemp()) < 0) {
	  use_rat = TRUE;
	  break;
	}
      }
      if (use_rat) { // a ^ b  with some b negative --> rational result
	// 1)  a := as.bigq(a, 1)
	SEXP one = Rf_ScalarInteger(1);
	PROTECT(one);
	SEXP aq = bigrational_as(a, one);
	PROTECT(aq);
	// 2)  result =  <bigq a> ^ b:
	SEXP ans = bigrational_pow(aq, b);
	UNPROTECT(2);
	return ans;
      }
    }
    // else, either, a has a modulus, or (no modulus *and*  exp >= 0) :
    return bigintegerR::biginteger_binary_operation(a,b, pow); // -> pow() in ./bigmod.cc
  }catch(std::invalid_argument & e){
    Rf_error("%s",e.what());
    return Rf_mkString(0);
 }
}
SEXP biginteger_inv (SEXP a, SEXP b) {return bigintegerR::biginteger_binary_operation(a,b,inv);}
SEXP biginteger_gcd (SEXP a, SEXP b) {return bigintegerR::biginteger_binary_operation(a,b,gcd);}
SEXP biginteger_lcm (SEXP a, SEXP b) {return bigintegerR::biginteger_binary_operation(a,b,lcm);}
SEXP biginteger_as (SEXP a, SEXP mod){return bigintegerR::biginteger_binary_operation(a,mod,set_modulus);}
//								set_modulus :  -> ./bigmod.cc

SEXP biginteger_lt  (SEXP a, SEXP b) {return bigintegerR::biginteger_logical_binary_operation(a,b,bigintegerR::lt);}
SEXP biginteger_gt  (SEXP a, SEXP b) {return bigintegerR::biginteger_logical_binary_operation(a,b,bigintegerR::gt);}
SEXP biginteger_lte (SEXP a, SEXP b) {return bigintegerR::biginteger_logical_binary_operation(a,b,bigintegerR::lte);}
SEXP biginteger_gte (SEXP a, SEXP b) {return bigintegerR::biginteger_logical_binary_operation(a,b,bigintegerR::gte);}
SEXP biginteger_eq  (SEXP a, SEXP b) {return bigintegerR::biginteger_logical_binary_operation(a,b,bigintegerR::eq);}
SEXP biginteger_neq (SEXP a, SEXP b) {return bigintegerR::biginteger_logical_binary_operation(a,b,bigintegerR::neq);}



SEXP biginteger_as_character(SEXP a, SEXP b)
{
  try{
    bigvec v = bigintegerR::create_bignum(a);
    SEXP ans;
    int base =  Rf_asInteger(b);
    if (base < 2 || base > 36){
      v.clear();
      throw invalid_argument(_("select a base between 2 and 36"));
    }
  
    PROTECT(ans = Rf_allocVector(STRSXP, v.size()));
    for (unsigned int i = 0; i < v.size(); ++i)
      SET_STRING_ELT(ans, i, Rf_mkChar(v.str(i,base).c_str()));
    // matrix part
    if(v.nrow >= 0)
      {
	SEXP dim = PROTECT(Rf_allocVector(INTSXP, 2));
	INTEGER(dim)[0] = v.nrow;
	INTEGER(dim)[1] = v.size() / v.nrow;
	Rf_setAttrib(ans, Rf_mkString("dim"), dim);
	UNPROTECT(1);
      }

    UNPROTECT(1);
    return ans;
  }catch(std::invalid_argument & e){
    Rf_error("%s",e.what());
    return Rf_mkString(0);
 }
}

SEXP biginteger_as_numeric(SEXP a)
{
  try{
    bigvec v = bigintegerR::create_bignum(a);
    SEXP ans = PROTECT(Rf_allocVector(REALSXP,v.size()));
    double *r = REAL(ans);
    for (unsigned int i = 0; i < v.size(); ++i)
      r[i] = v[i].isNA() ? NA_REAL : v[i].getValue().as_double();
    UNPROTECT(1);
    return ans;
  }catch(std::invalid_argument & e){
    Rf_error("%s",e.what());
    return Rf_mkString(0);
 }
}

SEXP biginteger_as_integer(SEXP a)
{
  try{
    bigvec v = bigintegerR::create_bignum(a);
    SEXP ans = PROTECT(Rf_allocVector(INTSXP,v.size()));
    int *r = INTEGER(ans);
    for (unsigned int i = 0; i < v.size(); ++i) {
      if(v[i].isNA()) {
	r[i] = NA_INTEGER;
      }
      else if(!mpz_fits_slong_p(v[i].getValue().getValueTemp())) {
	Rf_warning("NAs introduced by coercion from big integer");
	r[i] = NA_INTEGER;
      } else {
	r[i] = mpz_get_si(v[i].getValue().getValueTemp());
      }
    }
    UNPROTECT(1);
    return ans;
  }catch(std::invalid_argument & e){
    Rf_error("%s",e.what());
    return Rf_mkString(0);
  }
}

SEXP biginteger_get_at(SEXP a, SEXP i)
{
  //result = a [i]
  try{
    bigvec va = bigintegerR::create_bignum(a);
    return(bigintegerR::create_SEXP(bigintegerR::biginteger_get_at_C(va,i)));
  }catch(std::invalid_argument & e){
    Rf_error("%s",e.what());
    return Rf_mkString(0);
 }

}

// also called from  matrix_get_at_z(.)  in ./extract_matrix.cc :
bigvec bigintegerR::biginteger_get_at_C(bigvec va, SEXP ind)
{
  bigvec result;
  vector<int> v_ind = extract_gmp_R::indice_get_at(va.size(),ind);
  
  for(unsigned int i = 0 ; i < v_ind.size(); i++){
    int indice = v_ind[i];
    if(indice < (int) va.size()){
      result.push_back(va[indice]);
    } else {
      result.push_back(bigmod());
    }
  }

  return (result);
}

SEXP biginteger_set_at(SEXP src, SEXP idx, SEXP value)
{
  // return = ( src[idx] <- value )
  try{

  
    bigvec vvalue = bigintegerR::create_bignum(value);
    bigvec result = bigintegerR::create_bignum(src);
    vector<int>     vidx = extract_gmp_R::indice_get_at(result.size(),idx);
 
    if(vidx.size() == 0) {
      return bigintegerR::create_SEXP(result);
    }

    if(vvalue.size() == 0) {
      vvalue.clear();
      result.clear();
      throw invalid_argument(_("replacement has length zero"));
    }
 
    int pos = 0;

    for(unsigned int i = 0 ; i < vidx.size(); i++){
      while(result.size() <=(unsigned int) (vidx[i] )) {
	result.push_back(bigmod());
      }
      result.set(vidx[i],vvalue[pos++ % vvalue.size()]);
    }
  
    return bigintegerR::create_SEXP(result);
  }catch(std::invalid_argument & e){
    Rf_error("%s",e.what());
    return Rf_mkString(0);
 }

}

SEXP biginteger_length(SEXP a)
{
  try{
    return Rf_ScalarInteger(bigintegerR::create_bignum(a).size());
  }catch(std::invalid_argument & e){
    Rf_error("%s",e.what());
    return Rf_mkString(0);
 }
}

SEXP biginteger_setlength(SEXP vec, SEXP value)
{
  try{
    int len = 0;
    switch (TYPEOF(value)) {
    case INTSXP:
    case LGLSXP:
      if (LENGTH(value) != 1)
	Rf_error("%s",_("invalid second argument"));
      len =  Rf_asInteger(value);
      if (len < 0)
	Rf_error("%s",_("vector size cannot be negative"));
      else if (len == NA_INTEGER)
	Rf_error("%s",_("vector size cannot be NA"));
      break;
    case REALSXP:
      if (LENGTH(value) != 1)
	Rf_error("%s",_("invalid second argument"));
      len = (int)*REAL(value);
      if (len < 0)
	Rf_error("%s",_("vector size cannot be negative"));
      else if (! (R_FINITE (len ) ))
	Rf_error("%s",_("vector size cannot be NA, NaN of Inf"));
      break;
    case STRSXP:
      // dunno why R spits out this strange error on "length(foo) <- -1"
      // but I always follow the holy standard ;-)
      Rf_error("%s",_("negative length vectors are not allowed"));
    default:
      Rf_error("%s",_("invalid second argument"));
    }
    bigvec v =bigintegerR::create_bignum(vec);
    v.resize(len);
    return bigintegerR::create_SEXP(v);
  }catch(std::invalid_argument & e){
    Rf_error("%s",e.what());
    return Rf_mkString(0);
  }
}

SEXP biginteger_is_na(SEXP a)
{
  try{
    bigvec v = bigintegerR::create_bignum(a);
    SEXP ans = PROTECT(Rf_allocVector(LGLSXP, v.size()));
    for (unsigned int i = 0; i < v.size(); ++i)
      LOGICAL(ans)[i] = v[i].getValue().isNA();
    UNPROTECT(1);
    return ans;
  }catch(std::invalid_argument & e){
    Rf_error("%s",e.what());
    return Rf_mkString(0);
 }
}


SEXP biginteger_sgn(SEXP a)
{
  try{
    bigvec v = bigintegerR::create_bignum(a);
    SEXP ans = PROTECT(Rf_allocVector(INTSXP, v.size()));
    int *r = INTEGER(ans);
    for (unsigned int i = 0; i < v.size(); ++i)
      r[i] = mpz_sgn(v[i].getValue().getValueTemp());
    UNPROTECT(1);
    return ans;
  }catch(std::invalid_argument & e){
    Rf_error("%s",e.what());
    return Rf_mkString(0);
  }
}

SEXP biginteger_c(SEXP args)
{
  try{
    // if(TYPEOF(args) != VECSXP) Rf_error("%s","%s",_("should be a list"));
    bigvec result;
    for(int i=0; i < LENGTH(args); i++) {
      bigvec v = bigintegerR::create_bignum(VECTOR_ELT(args,i));
      for(unsigned int j=0; j < v.size(); j++){
	result.push_back(v[j]);
      }
      v.clear();
    }
    return bigintegerR::create_SEXP(result);
  }catch(std::invalid_argument & e){
    Rf_error("%s",e.what());
    return Rf_mkString(0);
  }
}


SEXP biginteger_rep(SEXP x, SEXP times)
{
  try{
    bigvec v = bigintegerR::create_bignum(x),
      result;
    int rep =  Rf_asInteger(times);


    for(int i = 0 ; i < rep ; i++)
      for(unsigned int j = 0 ; j < v.size() ; j++)
	result.push_back(v[j]);

    return bigintegerR::create_SEXP(result);
  }catch(std::invalid_argument & e){
    Rf_error("%s",e.what());
    return Rf_mkString(0);
 }
}


SEXP biginteger_is_prime(SEXP a, SEXP reps)
{
  try{
    bigvec v = bigintegerR::create_bignum(a);
    vector<int> vb = bigintegerR::create_int(reps);
    unsigned int i;
    SEXP ans = PROTECT(Rf_allocVector(INTSXP, v.size()));
    int *r = INTEGER(ans);
    if(v.size() == vb.size())
      for (i = 0; i < v.size(); ++i)
	r[i] = v[i].getValue().isprime(vb[i]);
    else
      for (i = 0; i < v.size(); ++i)
	r[i] = v[i].getValue().isprime(vb[0]);
    UNPROTECT(1);
    return ans;
  }catch(std::invalid_argument & e){
    Rf_error("%s",e.what());
    return Rf_mkString(0);
  }
}


SEXP biginteger_nextprime(SEXP a)
{
  try{
    bigvec v = bigintegerR::create_bignum(a),
      result;


    mpz_t val;
    mpz_init(val);
    mpz_t_sentry val_s(val);
    //    result.reserve(v.size());
      for (unsigned int i = 0; i < v.size(); ++i) {
	mpz_nextprime(val,v[i].getValue().getValueTemp());
	result.push_back(bigmod(val));
      }
    return bigintegerR::create_SEXP(result);
  }catch(std::invalid_argument & e){
    Rf_error("%s",e.what());
    return Rf_mkString(0);
 }
}

SEXP biginteger_abs(SEXP a)
{
  try{
    bigvec v = bigintegerR::create_bignum(a),
      result;


    mpz_t val;
    mpz_init(val);
    mpz_t_sentry val_s(val);

    for (unsigned int i = 0; i < v.size(); ++i)
      {
	mpz_abs(val,v[i].getValue().getValueTemp());
	result.push_back(bigmod(make_shared<biginteger>(val),v[i].getModulusPtr()));
      }

    return bigintegerR::create_SEXP(result);
  }catch(std::invalid_argument & e){
    Rf_error("%s",e.what());
    return Rf_mkString(0);
 }
}


/** @brief Bezoult coefficients: compute g,s and t as as + bt = g
 *  @param a BigInteger
 *  @param b BigInteger
 */
SEXP biginteger_gcdex(SEXP a, SEXP b)
{
  try{
    bigvec
      va = bigintegerR::create_bignum(a),
      vb = bigintegerR::create_bignum(b),
      result;

    if(va.size() != vb.size())
      return bigintegerR::create_SEXP(bigvec());


    mpz_t g;
    mpz_t s;
    mpz_t t;
    mpz_init(g);
    mpz_init(s);
    mpz_init(t);
    mpz_t_sentry val_g(g);
    mpz_t_sentry val_s(s);
    mpz_t_sentry val_t(t);
    shared_ptr<biginteger> mod = make_shared<biginteger>(); // NA

    for(unsigned int i=0; i < va.size(); i++)
      {
	mpz_gcdext (g,s,t,va[i].getValue().getValueTemp(),vb[i].getValue().getValueTemp());
	result.push_back(bigmod(make_shared<biginteger>(g))); 
	result.push_back(bigmod(make_shared<biginteger>(s)));
	result.push_back(bigmod(make_shared<biginteger>(t)));
      
      }
    return bigintegerR::create_SEXP(result);
  }catch(std::invalid_argument & e){
    Rf_error("%s",e.what());
    return Rf_mkString(0);
 }
}

/** @brief Random number generation
    \note If seed is not initialised: generation of a new seed
    @param nb  Integer: number of number to generate
    @param length Integer number will be of length 2^length
    @param newseed Integer, seed initialisation (if exists)
    @param ok Integer 1: seed generation 0 not
*/
SEXP biginteger_rand_u (SEXP nb, SEXP length, SEXP newseed, SEXP ok)
{
  try{
    mpz_t    bz;
    int i,flag,len,size;
    bigvec result;

    //extern int seed_init;
    //extern gmp_randstate_t seed_state;


    /* store input data into appropriate mode */
    bigvec va = bigintegerR::create_bignum(newseed);
    PROTECT (ok = AS_INTEGER(ok));
    PROTECT (length = AS_INTEGER(length));
    PROTECT (nb = AS_INTEGER(nb));
    flag =  Rf_asInteger(ok);
    len =  Rf_asInteger(length);
    size =  Rf_asInteger(nb);
    UNPROTECT(3);


    /* Random seed initialisation */

    if(seed_init==0)
      {
	gmp_randinit_default(seed_state);
	Rprintf("Seed default initialisation\n");
      }
    if(flag == 1)
      {
	gmp_randseed(seed_state,va[0].getValue().getValueTemp());
	Rprintf("Seed initialisation\n");
      }

    seed_init = 1;

    mpz_init (bz);
    mpz_t_sentry val_s(bz);

    for(i= 0; i<size; i++)
      {
	/*  Random number generation  */
	mpz_urandomb(bz,seed_state,len);
	result.push_back(bigmod(bz));
      }
    return bigintegerR::create_SEXP(result);
  }catch(std::invalid_argument & e){
    Rf_error("%s",e.what());
    return Rf_mkString(0);
 }
}


/** @brief biginteger_sizeinbase return
 *  @param x BigInteger
 *  @param base BigInteger
 */
SEXP biginteger_sizeinbase(SEXP x, SEXP base)
{
  try{
    bigvec vx = bigintegerR::create_bignum(x);
    int basesize= Rf_asInteger(base);
    SEXP ans = PROTECT(Rf_allocVector(INTSXP,vx.size()));
    int *r = INTEGER(ans);
    for(unsigned int i=0; i < vx.size(); i++)
      r[i] = mpz_sizeinbase(vx[i].getValue().getValueTemp(), basesize);
    UNPROTECT(1);
    return ans;
  }catch(std::invalid_argument & e){
    Rf_error("%s",e.what());
    return Rf_mkString(0);
 }
}


/** @brief bigI_factorial returns n!
 *  @param n non-negative integer vector
 */
SEXP bigI_factorial(SEXP n)
{
  try{
    bigvec result;
    int *nn = INTEGER(AS_INTEGER(n)), size = Rf_length(n);
    result.resize(size);
    for (int i = 0; i < size; ++i) {
      result[i].getValue().NA(false);
      if(nn[i] != NA_INTEGER && nn[i] >= 0) {
	mpz_fac_ui(result[i].getValue().getValue(), (unsigned long int)nn[i]);
      }
    }
    return bigintegerR::create_SEXP(result);
  }catch(std::invalid_argument & e){
    Rf_error("%s",e.what());
    return Rf_mkString(0);
}
} // bigI_factorial

/** @brief bigI_choose(n, k) returns binomial coefficient (n \choose k)
 *  @param n integer, either R "integer" (non-negative), or a "bigz"
 *  @param k non-negative integer
 */
SEXP bigI_choose(SEXP n, SEXP k)
{
  try{
    bigvec result, n_ = bigintegerR::create_bignum(n);
    int *kk = INTEGER(AS_INTEGER(k)), n_k = Rf_length(k);
    int size = (n_.size() == 0 || n_k == 0) ? 0 :
      // else:  max(n_.value.size(), n_k)
      (((int)n_.size() <= n_k) ? n_k : n_.size());

    result.resize(size);
    // if(n_.getType() ==  MODULUS_GLOBAL){
    //   result.setGlobalModulus(n_.getGlobalModulus());
    //}
    for (int i = 0; i < size; ++i) {
      result[i].getValue().NA(false);
      int ik_i = kk[i % n_k];
      // check if k in range:
      if(ik_i != NA_INTEGER && ik_i >= 0) {
	unsigned long int k_i = (unsigned long int)ik_i;
	/* void mpz_bin_ui (mpz_t ROP, mpz_t N, unsigned long int K) */
	mpz_bin_ui(result[i].getMpValue(),
		   n_[i % n_.size()].getValueTemp(), k_i);
      }
    }
    return bigintegerR::create_SEXP(result);
  }catch(std::invalid_argument & e){
    Rf_error("%s",e.what());
    return Rf_mkString(0);
 }
}

/** @brief fibnum return nth Fibonacci number
 *  @param n integer
 */
SEXP bigI_fibnum(SEXP n)
{
  try{
    bigvec result;
    if(Rf_length(n) > 0)
      {
	int nn = Rf_asInteger(n);
	unsigned long int num = nn;
	if(nn < 0 || nn == NA_INTEGER)
	  throw invalid_argument (_("argument must be non-negative"));
	mpz_t val;
	mpz_init(val);
	mpz_t_sentry val_s(val);

	mpz_fib_ui(val,num);
	result.push_back(bigmod(val));
	//      result[0].value.setValue(val);
      }
    // else
    //   Rf_error("%s","%s",_("argument must not be an empty list"));

    return bigintegerR::create_SEXP(result);
  }catch(std::invalid_argument & e){
    Rf_error("%s",e.what());
  }
  return Rf_mkString(0);

}

/** @brief fibnum2 return nth and n-1th Fibonacci number
 *  @param n integer
 */
SEXP bigI_fibnum2(SEXP n)
{
  try{
    bigvec result;
    if(Rf_length(n) > 0)
      {
	int nn = Rf_asInteger(n);
	unsigned long int num = nn;
	if(nn < 0 || nn == NA_INTEGER)
	  throw invalid_argument(_("argument must be non-negative"));

	mpz_t val;
	mpz_init(val);
	mpz_t_sentry val_s(val);
	mpz_t val2;
	mpz_init(val2);
	mpz_t_sentry val_s2(val2);

	mpz_fib2_ui(val,val2, num);
	result.push_back(bigmod(val2));
	result.push_back(bigmod(val));
      }
    else
      throw invalid_argument(_("argument must not be an empty list"));

    return bigintegerR::create_SEXP(result);
  }catch(std::invalid_argument & e){
    Rf_error("%s",e.what());
    return Rf_mkString(0);

  }
}

/** @brief lucnum return nth lucas number
 *  @param n integer
 */
SEXP bigI_lucnum(SEXP n)
{
  try{
    bigvec result;
    if(Rf_length(n) > 0)
      {
	int nn = Rf_asInteger(n);
	unsigned long int num = nn;
	if(nn < 0 || nn == NA_INTEGER)
	  throw invalid_argument(_("argument must be non-negative"));

	mpz_t val;
	mpz_init(val);
	mpz_t_sentry val_s(val);

	mpz_lucnum_ui(val,num);
	result.push_back(bigmod(val));
      }
    // else
    //   Rf_error("%s","%s",_("argument must not be an empty list"));

    return bigintegerR::create_SEXP(result);
  }catch(std::invalid_argument & e){
    Rf_error("%s",e.what());
    return Rf_mkString(0);
 }
}

/** @brief lucnum2 return nth and n-1th lucas number
 *  @param n integer
 */
SEXP bigI_lucnum2(SEXP n)
{
 try{
    bigvec result;

    if(Rf_length(n) > 0) {
      int nn = Rf_asInteger(n);
      unsigned long int num = nn;
      if(nn < 0 || nn == NA_INTEGER)
	throw invalid_argument(_("argument must be non-negative"));
      mpz_t val;
      mpz_init(val);
      mpz_t_sentry val_s(val);
      mpz_t val2;
      mpz_init(val2);
      mpz_t_sentry val_s2(val2);

      mpz_lucnum2_ui(val,val2,num);
      result.push_back(bigmod(val2));
      result.push_back(bigmod(val));
    }
    else
      throw invalid_argument(_("argument must not be an empty list"));

    return bigintegerR::create_SEXP(result);
  }catch(std::invalid_argument & e){
    Rf_error("%s",e.what());
    return Rf_mkString(0);
 }
}



//Return max

SEXP biginteger_max(SEXP a, SEXP narm)
{
  try{
    bigvec result;
    bigvec va = bigintegerR::create_bignum(a);

    if( ! va.size())
      return bigintegerR::create_SEXP(result);

    unsigned int maximum = 0;
    int na_remove = Rf_asInteger(narm);

    for(unsigned int i = 1 ; i < va.size(); ++i)
      {
	if(va[i].isNA() && !na_remove)
	  return(bigintegerR::create_SEXP(result));
	else
	  if(va[i].getValue() > va[maximum].getValue() )
	    maximum = i; // if va.value[maximum = 0] is NA => false for the "<" => maximum changed = good
      }

    result.push_back(va[maximum]);

    if(va.getType() == TYPE_MODULUS::MODULUS_BY_CELL) result[0].getModulus().NA(TRUE);

    return(bigintegerR::create_SEXP(result));
  }catch(std::invalid_argument & e){
    Rf_error("%s",e.what());
    return Rf_mkString(0);
 }
}


// Return min
SEXP biginteger_min(SEXP a, SEXP narm)
{
  try{
    bigvec result;
    bigvec va = bigintegerR::create_bignum(a);

    if( ! va.size())
      return bigintegerR::create_SEXP(result);

    unsigned int minimum = 0;
    int na_remove = Rf_asInteger(narm);

    for(unsigned int i = 1 ; i < va.size(); ++i)
      {
	if(va[i].isNA() && !na_remove)
	  return bigintegerR::create_SEXP(result);
	else
	  if(va[i].getValue() <  va[minimum].getValue() )
	    minimum = i;
      }

    result.push_back(va[minimum]);

    // now the modulus !
    if(va.getType() == TYPE_MODULUS::MODULUS_BY_CELL) result[0].getModulus().NA(TRUE);

    return bigintegerR::create_SEXP(result);
  }catch(std::invalid_argument & e){
    Rf_error("%s",e.what());
    return Rf_mkString(0);
 }
}


SEXP biginteger_cumsum(SEXP a)
{
  try{
    bigvec result, va = bigintegerR::create_bignum(a);

    result.resize(va.size());

    mpz_t val;
    mpz_init(val);
    mpz_t_sentry val_s(val);

    bool hasmodulus = va.getType() == TYPE_MODULUS::MODULUS_GLOBAL;

 
    for(unsigned int i = 0 ; i < va.size(); ++i)
      {
	{
	  if(va[i].isNA() )
	    {
	      break; // all last values are NA.
	    }

	  mpz_add(val,val,va[i].getValue().getValueTemp());

	  if(hasmodulus){
	    mpz_mod(val,val,va.getGlobalModulus()->getValueTemp() );
	    result[i].setModulus(va.getGlobalModulus());
	  }
	  result[i].getValue().setValue(val);
	
	}
      }

    return(bigintegerR::create_SEXP(result));
  }catch(std::invalid_argument & e){
    Rf_error("%s",e.what());
    return Rf_mkString(0);
 }
}



SEXP biginteger_sum(SEXP a)
{
  try{
    bigvec result, va = bigintegerR::create_bignum(a);

    result.resize(1);

    mpz_t val;
    mpz_init(val);
    mpz_t_sentry val_s(val);

    bool hasmodulus = va.getType() == TYPE_MODULUS::MODULUS_GLOBAL;



    for(unsigned int i = 0 ; i < va.size(); ++i)
      {
	{
	  if(va[i].isNA() )
	    {
	      break; // all last values are NA.
	    }

	  mpz_add(val,val,va[i].getValueTemp());

	  if(hasmodulus)
	    mpz_mod(val,val,va.getGlobalModulus()->getValueTemp() );
	}
      }

    result[0].setValue(val);
    if(hasmodulus) result[0].setModulus(va.getGlobalModulus());

    return(bigintegerR::create_SEXP(result));
  }catch(std::invalid_argument & e){
    Rf_error("%s",e.what());
    return Rf_mkString(0);
 }
}


SEXP biginteger_prod(SEXP a)
{
  try{
    bigvec result;
    bigvec va = bigintegerR::create_bignum(a);

    result.resize(1);

    mpz_t val;
    mpz_init(val);
    mpz_set_ui(val,1);
    mpz_t_sentry val_s(val);

    bool hasmodulus = va.getType() == TYPE_MODULUS::MODULUS_GLOBAL;

  
    for(unsigned int i = 0 ; i < va.size(); ++i)
      {
	if(va[i].isNA() )
	  {
	    return (bigintegerR::create_SEXP(result));
	  }

	mpz_mul(val,val,va[i].getValue().getValueTemp());

	if(hasmodulus)
	  mpz_mod(val,val,va.getGlobalModulus()->getValueTemp() );
      }

    result[0].setValue(val);
    if(hasmodulus) result[0].setModulus(va.getGlobalModulus());

    return(bigintegerR::create_SEXP(result));
  }catch(std::invalid_argument & e){
    Rf_error("%s",e.what());
    return Rf_mkString(0);
  }
}

// return x ^ y [n]
SEXP biginteger_powm(SEXP x, SEXP y, SEXP n)
{
  try{
    bigvec result;

    bigvec vx = bigintegerR::create_bignum(x);
    bigvec vy = bigintegerR::create_bignum(y);
    bigvec vn = bigintegerR::create_bignum(n);

    result.resize(vx.size());
    // this cause pb to package
    // homeomopheR result.setGlobalModulus(vn[0].getValuePtr());
    // ...
    for (unsigned int i = 0 ; i < vx.size(); i++)
      {

	result[i].getValue().NA(false);
	// check if n != 0
	if(mpz_sgn(vn[i % vn.size()].getValueTemp()) != 0)
	  mpz_powm(result[i].getMpValue(),
		   vx[i].getValue().getValueTemp(),
		   vy[i % vy.size()].getValueTemp(),
		   vn[i % vn.size()].getValueTemp());

      }

    return bigintegerR::create_SEXP(result);
  }catch(std::invalid_argument & e){
    Rf_error("%s",e.what());
    return Rf_mkString(0);
 }
} // ..._powm()


// TODO:  A version that only returns  'ex'  {i.e. the binary precision}
// ----   GMP manual suggests that     size_t mpz_sizeinbase (mpz_t OP, int BASE)
//        i.e.,  mpz_sizeinbase (OP, 2)  would give that

// (d, ex) where  x = d * 2^ex, and  0.5 <= |d| < 1
SEXP bigI_frexp(SEXP x)
{
  // double mpz_get_d_2exp (signed long int *exp, mpz t op )

  // Convert op to a double, truncating if necessary (ie. rounding towards zero), and returning
  // the exponent separately.
  // The return value is in the range 0.5 ≤ |d| < 1 and the exponent is stored to *exp . d ∗ 2exp is
  // the (truncated) op value. If op is zero, the return is 0.0 and 0 is stored to *exp .
  // This is similar to the standard C frexp function.
  
  bigvec vx;
  try{
    vx = bigintegerR::create_bignum(x);
  }catch(std::invalid_argument & e){
    Rf_error("%s",e.what());
  }
  
  const char *nms[] = {"d", "exp", ""};
  SEXP ans, d_R, exp_R;
  int n = vx.size();

  PROTECT(ans = Rf_mkNamed(VECSXP, nms)); // =~= R  list(d = . , exp= .)
  d_R   = Rf_allocVector(REALSXP, n); SET_VECTOR_ELT(ans, 0, d_R);
  exp_R = Rf_allocVector(INTSXP,  n); SET_VECTOR_ELT(ans, 1, exp_R);
  double *d_ = REAL(d_R);
  int  *exp_ = INTEGER(exp_R);
  try{
    for (int i = 0 ; i < n; i++) {
      signed long int ex;
      d_[i] = mpz_get_d_2exp(&ex, vx[i].getValueTemp());
      if(abs(ex) < INT_MAX){
	exp_[i] = (int) ex;
      }
      else{
	vx.clear();
	throw invalid_argument(_("exponent too large to fit into an integer"));
      }
    }
  }catch(std::invalid_argument & e){
    Rf_error("%s",e.what());
 }
  UNPROTECT(1);
  return ans;
  
} // bigI_frexp()


SEXP biginteger_log2(SEXP x)
{
  bigvec v;
  try{
    v = bigintegerR::create_bignum(x);
  }
  catch(std::invalid_argument & e){
    Rf_error("%s",e.what());
  }
  SEXP ans = PROTECT(Rf_allocVector(REALSXP,v.size()));
  double *r = REAL(ans);
  
  for (unsigned int i = 0; i < v.size(); ++i) {
    signed long int ex;
    double di = mpz_get_d_2exp(&ex, v[i].getValueTemp());
    // xi = di * 2 ^ ex  ==> log2(xi) = log2(di) + ex :
    r[i] = log(di) / M_LN2 + (double)ex;
  }
  UNPROTECT(1);
  return ans;
  
}

SEXP biginteger_log(SEXP x)
{
  bigvec v ;
  try{
    v = bigintegerR::create_bignum(x);
  }catch(std::invalid_argument & e){
    Rf_error("%s",e.what());
  }
  SEXP ans = PROTECT(Rf_allocVector(REALSXP,v.size()));
  double *r = REAL(ans);
  
  for (unsigned int i = 0; i < v.size(); ++i) {
    signed long int ex;
    double di = mpz_get_d_2exp(&ex, v[i].getValueTemp());
    // xi = di * 2 ^ ex  ==> log(xi) = log(di) + ex*log(2) :
    r[i] = log(di) + M_LN2*(double)ex;
  }
  UNPROTECT(1);
  return ans;
 
}
