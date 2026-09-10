/*! \file bigrational.cc
 *  \brief Implementation of class bigrational
 *
 *  \version 1
 *
 *  \date Created: 12/12/04
 *  \date Last modified: Time-stamp: <2023-01-24 19:01:08 (antoine)>
 *
 *  \author Antoine Lucas (adapted from biginteger class made by
 *                         Immanuel Scholz)
 *
 *  \note Licence: GPL (>= 2)
 */

#include <cstdio>

#include "bigrational.h"
#include "bigrationalR.h"

using std::string;

static int count=0;

/**
 * Construct a "NA" bigrational.
 */
bigrational::bigrational() :
  value(),
  na(true) {
  count++;
  mpq_init(value);
}


// set numerator
bigrational::bigrational(void* raw):
  value(),
  na(true)
{
  count++;
  mpz_t tmpVal;
  mpz_init(tmpVal);
  mpz_t_sentry val_s(tmpVal);

  mpq_init(value);

  int* r = (int*)raw;
  if (r[0]>0)
    {
      mpz_import(tmpVal, r[0], 1, sizeof(int), 0, 0, &r[2]);
      if(r[1]==-1)
	mpz_neg(tmpVal,tmpVal);
      na = false;
      mpq_set_z(value,tmpVal);
    }

}

/**
 * Create a bigrational from a value. Remember to free the
 * parameter's mpz_t if you allocated them by yourself -
 * biginteger will copy the value.
 */
bigrational::bigrational(const mpq_t& value_) :
  value(),
  na(false)
{
  count++;
  mpq_init(value);
  mpq_set(value, value_);
}

/**
 * \brief create a rational from an [big] integer
 */
bigrational::bigrational(const mpz_t& value_) :
  value(),
  na(false)
{
  count++;
  mpq_init(value);
  mpq_set_z(value, value_);
}

/**
 * Construct a bigrational from a long value.
 */
bigrational::bigrational(int value_) :
  value(),
  na(false) {
  count++;
  mpq_init(value);
  if(value_ ==  NA_INTEGER)
    na = true  ;
  else
    mpq_set_si(value, value_,1);
}

/**
 * Construct a bigrational from a long value.
 */
bigrational::bigrational(int num_, int den_) :
  value(),
  na(false) {
  count++;
  mpq_init(value);
  if((num_ ==  NA_INTEGER) || (den_ == NA_INTEGER) )
    na = true  ;
  else
    mpq_set_si(value, num_,den_);}

/**
 * Construct a bigrational from a double value.
 */
bigrational::bigrational(double value_) :
  value(),
  na(false)
{
  count++;
  mpq_init(value);
  if(R_FINITE( value_ ) )
    mpq_set_d(value, value_);
  else // FIXME: consider  "1/0" and "(-1)/0" for  +- Inf
    na = true  ;
}

/**
 * Construct a bigrational from a string value. it can be "4343" or "2322/4343"
 */
bigrational:: bigrational(const std::string& value_) :
  value(),
  na(false)
{
  count++;
  mpq_init(value);
  /* mpz_init.. return -1 when error, 0: ok */
  if(mpq_set_str(value, value_.c_str(), 0))
    na=true;
  /*	if(mpz_init_set_str(value, value_.c_str(), 0) == -1)
	Rf_error("Not a valid number");    */
}

/**
 *  Copy constructor (mpz_t aren't standard-copyable)
 */
bigrational::bigrational(const biginteger & rhs) :
  value(),
  na(rhs.isNA())
{
  count++;
  mpq_init(value);
  mpq_set_z(value, rhs.getValueTemp());
}

/**
 *  Copy constructor (mpz_t aren't standard-copyable)
 */
bigrational::bigrational(const bigrational & rhs) :
  value(),
  na(rhs.na)
{
  count++;
  mpq_init(value);
  mpq_set(value, rhs.value);
}

/**
 * Free the owned mpz_t structs
 */
bigrational::~bigrational() {
  count--;
  //printf("bigq count %d \n",count);
  mpq_clear(value);}

bigrational & bigrational::operator= (const bigrational& rhs)
{
  if(this != &rhs)
    {
      mpq_set(value, rhs.getValueTemp());
      na= rhs.na;
    }
  return(*this);

}


/**
 * \brief Print value
 */
std::string bigrational::str(int b) const
{
  if (isNA())
    return "NA";

  unsigned int totSize = mpz_sizeinbase(mpq_numref(value),b) + \
    mpz_sizeinbase(mpq_denref(value),b) + 3 ;
  char* buf = new char[totSize];
  mpq_get_str(buf, b, value);
  std::string s = buf;
  delete [] buf;
  return s;

}


/* \brief  simplify n/d (use of mpq_canonical)
 *
 */
void  bigrational::simplify ()
{
  mpq_canonicalize(value);
}


// size of numerator !!
size_t bigrational::raw_size() const
{
  if (isNA())
    return sizeof(int);

  int numb = 8*sizeof(int);

  return sizeof(int) * (2 + (mpz_sizeinbase(mpq_numref(value),2)+numb-1) / numb);
}

bigrational operator+(const bigrational& lhs, const bigrational& rhs)
{
  return bigrationalR::create_bigrational(lhs, rhs, mpq_add);
}

/**
 * \brief Return  a - b
 */
bigrational operator-(const bigrational& lhs, const bigrational& rhs)
{
  return bigrationalR::create_bigrational(lhs, rhs, mpq_sub);
}

/**
 * \brief Return  a * b
 */
bigrational operator*(const bigrational& lhs, const bigrational& rhs)
{
  return bigrationalR::create_bigrational(lhs, rhs, mpq_mul);
}

/**
 * \brief Return  a / b
 */
bigrational operator/(const bigrational& lhs, const bigrational& rhs)
{

  return bigrationalR::create_bigrational(lhs, rhs, mpq_div, false);
}

/**
 * \brief Return  a ^ b
 */
bigrational operator^(const bigrational& lhs, const biginteger& rhs)
{
  // if (base == 1  or  exp == 0)  return 1
  if((!lhs.isNA() && !mpq_cmp_si(lhs.getValueTemp(), 1,1)) ||
     (!rhs.isNA() && !mpz_cmp_si(rhs.getValueTemp(), 0)))
    return bigrational(1);
  if (lhs.isNA() || rhs.isNA())
    return bigrational();

  return bigrationalR::create_bigrational_z(lhs, rhs, bigrationalR::mpqz_pow);
}

//
bool operator!=(const bigrational& lhs, const bigrational& rhs)
{
  if(rhs.isNA() || lhs.isNA())
    return(false); // SHOULD RETURN NA

  return(mpq_cmp(lhs.getValueTemp(),rhs.getValueTemp()) != 0);
}

bool operator>(const bigrational& lhs, const bigrational& rhs)
{
  if(rhs.isNA() || lhs.isNA())
    return(false); // SHOULD RETURN NA

  return(mpq_cmp(lhs.getValueTemp(),rhs.getValueTemp()) > 0);
}

bool operator<(const bigrational& lhs, const bigrational& rhs)
{
  if(rhs.isNA() || lhs.isNA())
    return(false); // SHOULD RETURN NA

  return(mpq_cmp(lhs.getValueTemp(),rhs.getValueTemp()) < 0);
}

/**
 * \brief Well... an heritage from biginteger class, this should be
 * integrated earlier... put denominator & simplify if there is not.
 */
// R  as.bigq() :
bigrational set_denominator(const bigrational& lhs, const bigrational& rhs)
{
  return bigrationalR::create_bigrational(lhs, rhs, mpq_div, false);
}

// return 1/x
bigrational bigrational::inv()
{
  if(isNA())
    return(bigrational());

  mpq_t tmpVal;
  mpq_init(tmpVal);
  mpq_t_sentry val_s(tmpVal);

  mpq_inv(tmpVal,value);

  return(bigrational(tmpVal));
}
