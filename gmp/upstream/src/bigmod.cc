#include "bigmod.h"

#define DEBUG_bigmod
#undef DEBUG_bigmod

#ifdef DEBUG_bigmod
#include <R_ext/Print.h>
#endif

// for GetOption etc:
#include <R.h>
#include <Rinternals.h>
using namespace std;
#include <stdexcept>

// string representation of (.) wrt base 'b' :
std::string bigmod::str(int b) const
{
  if (value->isNA())
    return "NA";

  std::string s; // sstream seems to collide with libgmp :-(
  if (!modulus->isNA())
    s = "(";
  s += value->str(b);
  if (!modulus->isNA()) {
    s += " %% ";
    s += modulus->str(b);
    s += ")";
  }
  return s;
}

bigmod & bigmod::operator= (const bigmod& rhs)
{
  if(this != &rhs)
    {
      modulus = make_shared<biginteger>(rhs.getModulus());
      value = make_shared<biginteger>(rhs.getValue() );
    }
  return(*this);
}

bigmod bigmod::inv () const
{
  if(value->isNA() || modulus->isNA()) {
    return bigmod();
  }
    
  mpz_t val;
  mpz_init(val);
  mpz_t_sentry val_s(val);
  if (mpz_invert(val, getValue().getValueTemp(), getModulus().getValueTemp()) == 0) {
    SEXP wOpt = Rf_GetOption1(Rf_install("gmp:warnNoInv"));
    if(wOpt != R_NilValue && Rf_asInteger(wOpt))
      Rf_warning("inv(x) returning NA as x has no inverse");
    
    return bigmod(); // return NA; was
  }

  return bigmod(std::make_shared<biginteger>(val), std::make_shared<biginteger>(modulus->getValueTemp()));
}


bool operator!=(const bigmod& rhs, const bigmod& lhs)
{
  if(rhs.getValue() != lhs.getValue())
    return(true);
  return(rhs.getModulus() != lhs.getModulus());
}

bool operator==(const bigmod& rhs, const bigmod& lhs)
{
  
  if(rhs.getValue() != lhs.getValue())
    return(false);
  // if(rhs.getModulusPtr().get() == lhs.getModulusPtr().get()) return true; // same pointer
  return(!(rhs.getModulus() != lhs.getModulus()));
}


bigmod operator+(const bigmod& lhs, const bigmod& rhs)
{
  return create_bigmod(lhs, rhs, mpz_add);
}

bigmod operator-(const bigmod& lhs, const bigmod& rhs)
{
  return create_bigmod(lhs, rhs, mpz_sub);
}

bigmod operator*(const bigmod& lhs, const bigmod& rhs)
{
  return create_bigmod(lhs, rhs, mpz_mul);
}

/* called via biginteger_binary_operation(.) from biginteger_div() in
 * ./bigintegerR.cc  from R's  .Call(biginteger_div,  a, b)
 *   ~~~~~~~~~~~~~~
 * itself called from  "/.bigz" = div.bigz()
 */
bigmod div_via_inv(const bigmod& a, const bigmod& b) {
    // compute  a/b  as  a * b^(-1)
    return operator*(a, pow(b, bigmod(-1)));
}


void integer_div(mpz_t result,const mpz_t a, const mpz_t b) {
  mpz_tdiv_q(result,a,b);
  //
  // si resulat < 0 et module != 0: on enleve 1.
  // i.e. 3 / -4 = -1
  if (mpz_sgn(a) * mpz_sgn(b) == -1) {
    mpz_t val;
    mpz_init(val);
    mpz_t_sentry val_s(val);
    mpz_mod(val, a, b);
    if (mpz_cmp_ui(val, 0) != 0) 
      {	
	mpz_sub_ui(result, result, 1);
      }
  }
}


/* called via biginteger_binary_operation(.) from R's
 * .Call(biginteger_divq, a, b) , itself called from '%/%.bigz' = divq.bigz()
 */
bigmod operator/(const bigmod& lhs, const bigmod& rhs) {
  return create_bigmod(lhs, rhs, integer_div, false);
}

bigmod operator%(const bigmod& lhs, const bigmod& rhs)
{
  if (lhs.getValue().isNA() || rhs.getValue().isNA())
    return bigmod();
  if (mpz_sgn(rhs.getValue().getValueTemp()) == 0) {
    Rf_warning("biginteger division by zero: returning NA");
    return bigmod();
  }
  biginteger mod;
  if (!lhs.getModulus().isNA() || !rhs.getModulus().isNA())
    mod = rhs.getValue();

  mpz_t val;
  mpz_init(val);
  mpz_t_sentry val_s(val);
  mpz_mod(val, lhs.getValue().getValueTemp(), rhs.getValue().getValueTemp());
  return bigmod(val, mod);
}


// Either 'base' has a modulus, or it has not *and* exp >= 0 :

bigmod pow(const bigmod& base, const bigmod& exp)
{
  biginteger mod = get_modulus(base, exp);
#ifdef DEBUG_bigmod
  if(mod.isNA() && !mpz_cmp_si(base.getValue().getValueTemp(), 1))
    Rprintf("bigmod pow(1, exp=%d)\n", mpz_get_si(exp.getValue().getValueTemp()));
  else if(mod.isNA() && !mpz_cmp_si(exp.getValue().getValueTemp(), 0))
    Rprintf("bigmod pow(base=%d, 0)\n", mpz_get_si(base.getValue().getValueTemp()));
#endif

  // if (base == 1  or  exp == 0)  return 1
  if(mod.isNA() &&
     ((!base.getValue().isNA() && !mpz_cmp_si(base.getValue().getValueTemp(), 1)) ||
      (! exp.getValue().isNA() && !mpz_cmp_si( exp.getValue().getValueTemp(), 0))))
    return bigmod(biginteger(1));
  if (base.getValue().isNA() || exp.getValue().isNA())
    return bigmod();
  int sgn_exp = mpz_sgn(exp.getValue().getValueTemp());
  bool neg_exp = (sgn_exp < 0); // b ^ -|e| =  1 / b^|e|
  mpz_t val; mpz_init(val); mpz_t_sentry val_s(val);
#ifdef DEBUG_bigmod
  Rprintf("bigmod pow(base=%3s, exp=%3s [mod=%3s]) ..\n",
          base.getValue().str(10).c_str(), exp.getValue().str(10).c_str(),
	  mod.str(10).c_str());
#endif
  if (mod.isNA()) { // <==> (both have no mod || both have mod. but differing)
    if(neg_exp){
      throw invalid_argument(_("** internal error (negative powers for Z/nZ), please report!"));
    }
    if (!mpz_fits_ulong_p(exp.getValue().getValueTemp()))
      throw invalid_argument(_("exponent e too large for pow(z,e) = z^e"));// FIXME? return( "Inf" )
    // else :
    mpz_pow_ui(val, base.getValue().getValueTemp(),
 	       mpz_get_ui(exp.getValue().getValueTemp()));
  }
  else if( mpz_sgn(mod.getValueTemp()) != 0) { // check modulus non-zero
    if(neg_exp) { // negative exponent -- only ok if inverse exists
      if (mpz_invert(val, base.getValue().getValueTemp(), mod.getValueTemp()) == 0) {
	SEXP wOpt = Rf_GetOption1(Rf_install("gmp:warnNoInv"));
	if(wOpt != R_NilValue && Rf_asInteger(wOpt))
	  Rf_warning("pow(x, -|n|) returning NA as x has no inverse wrt modulus");
	return(bigmod()); // return NA; was
      } // else: val = x^(-1) already: ==> result =  val ^ |exp| =  val ^ (-exp) :
      // nExp := - exp
      mpz_t nExp;
      mpz_init(nExp);
      mpz_t_sentry val_ex(nExp);
      mpz_neg(nExp, exp.getValue().getValueTemp());
      mpz_powm(val, val, nExp, mod.getValueTemp());
    } else { // non-negative exponent
      mpz_powm(val, base.getValue().getValueTemp(), exp.getValue().getValueTemp(), mod.getValueTemp());
    }
  }
  return bigmod(val, mod);
}

bigmod inv(const bigmod& x, const bigmod& m)
{
  if (x.getValue().isNA() || m.getValue().isNA())
    return bigmod();
  SEXP wOpt = Rf_GetOption1(Rf_install("gmp:warnNoInv"));
  bool warnI = (wOpt != R_NilValue && Rf_asInteger(wOpt));
  if (mpz_sgn(m.getValue().getValueTemp()) == 0) {
    if(warnI) Rf_warning("inv(0) returning NA");
    return bigmod();
  }
  biginteger mod = get_modulus(x, m);
  mpz_t val;
  mpz_init(val);
  mpz_t_sentry val_s(val);
  if (mpz_invert(val, x.getValue().getValueTemp(), m.getValue().getValueTemp()) == 0) {
    if(warnI) Rf_warning("inv(x,m) returning NA as x has no inverse modulo m");
    return(bigmod()); // return NA; was
  }
  return bigmod(val, mod);
}

// R  as.bigz() :
bigmod set_modulus(const bigmod& x, const bigmod& m)
{
  if (!m.getValue().isNA() && mpz_sgn(m.getValue().getValueTemp()) == 0){
    throw invalid_argument("modulus 0 is invalid");
  }
  //    if (!m.getValue().isNA() && mpz_cmp(x.getValue().getValueTemp(),m.getValue().getValueTemp())>=0) {
  if (!m.getValue().isNA() ) {
    bigmod t(x%m);
    return bigmod(t.getValue(), m.getValue());
  } else
    return bigmod(x.getValue(), m.getValue());
}

bigmod gcd(const bigmod& lhs, const bigmod& rhs)
{
  return create_bigmod(lhs, rhs, mpz_gcd);
}

bigmod lcm(const bigmod& lhs, const bigmod& rhs)
{
  return create_bigmod(lhs, rhs, mpz_lcm);
}


// return the modulus to use for the two bigmods.
// NA if incompatible.
biginteger get_modulus(const bigmod& b1, const bigmod& b2)
{
  if (b1.getModulus().isNA()) // NA: means "no modulus" <==> R's is.null(modulus(.))
    return b2.getModulus(); // if b2 is NA too, the return is correct: NA
  else if (b2.getModulus().isNA())
    return b1.getModulus();
  else if (mpz_cmp(b1.getModulus().getValueTemp(), b2.getModulus().getValueTemp())) {
    SEXP wOpt = Rf_GetOption1(Rf_install("gmp:warnModMismatch"));
    if(wOpt != R_NilValue && Rf_asInteger(wOpt))
      Rf_warning("modulus mismatch in bigz.* arithmetic");
    return biginteger(); // i.e. NA
  } else // equal
    return b1.getModulus();
}


//    typedef void (*gmp_binary)(mpz_t, const mpz_t, const mpz_t);




// Create a bigmod from a binary combination of two other bigmods
bigmod create_bigmod(const bigmod& lhs, const bigmod& rhs, gmp_binary f,
		     bool zeroRhsAllowed) {
  if (lhs.getValue().isNA() || rhs.getValue().isNA())
    return bigmod();
  if (!zeroRhsAllowed && mpz_sgn(rhs.getValue().getValueTemp()) == 0) {
    Rf_warning("returning NA  for (modulus) 0 in RHS");
    return bigmod();
  }
  biginteger mod = get_modulus(lhs, rhs);
  mpz_t val;
  mpz_init(val);
  mpz_t_sentry val_s(val);
  f(val, lhs.getValue().getValueTemp(), rhs.getValue().getValueTemp());
  //--- val := f(lhs, rhs)
#ifdef DEBUG_bigmod
  bool iNA = biginteger(val).isNA();
  char* buf;
  if(iNA)
    buf = NULL;
  else {
    buf = new char[mpz_sizeinbase(val, 10)+2];
    // possible minus sign, size of number + '\0'
    mpz_get_str(buf, 10, val);
  }
  Rprintf("create_bigmod(lhs=%3s, rhs=%3s [mod=%3s]) = %s%s",
	  lhs.getValue().str(10).c_str(),
	  rhs.getValue().str(10).c_str(),
	  mod.str(10).c_str(),
	  (iNA)? "NA" : buf,
	  (mod.isNA())? "\n" : " {before 'mod'}");
#endif
  if (!mod.isNA()) {
    mpz_mod(val, val, mod.getValueTemp());
#ifdef DEBUG_bigmod
    if(biginteger(val).isNA())
      Rprintf(" -> val = NA\n");
    else {
      mpz_get_str(buf, 10, val);
      Rprintf(" -> val = %s\n", buf);
    }
#endif
  }
  return bigmod(val, mod);
}
