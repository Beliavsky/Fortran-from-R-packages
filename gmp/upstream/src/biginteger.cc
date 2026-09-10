/*! \file biginteger.cc
 *  \brief C function for class biginteger & bigmod
 *
 *  \version 1
 *
 *  \date Created: 27/10/04
 *  \date Last modified: Time-stamp: <2023-01-28 15:51:50 (antoine)>
 *
 *  \author Immanuel Scholz
 *
 *  \note Licence: GPL (>= 2)
 */

#include <cstdio>
#include <iostream>

#include "biginteger.h"

using std::string;
static int count=0;
static int countALL=0;

biginteger::biginteger() 
  : na(true) {
  count++;
  countALL++;
  mpz_init(value);}



/**
 * Construct a biginteger from a int value.
 */
biginteger::biginteger(const int value_) : na(false) {
  count++;
  countALL++;
  if(value_ ==  NA_INTEGER)
    {mpz_init(value); na = true  ;}
  else
    mpz_init_set_si(value, value_);}

/**
 * Construct a biginteger from a long value.
 */
biginteger::biginteger(const long int value_) : na(false) {
  count++;
  countALL++;
  if(value_ ==  NA_INTEGER)
    {mpz_init(value); na = true  ;}
  else
    mpz_init_set_si(value, value_);}

/**
 * Construct a biginteger from a unsigned long value.
 */
biginteger:: biginteger(const unsigned long int value_) : na(false) {
  count++;
  countALL++;
  if(value_ == (unsigned long int) NA_INTEGER)
    {mpz_init(value); na = true  ;}
  else
    mpz_init_set_ui(value, value_);}

/**
 * Construct a biginteger from a double value.
 */
biginteger::biginteger(const double value_) : na(false) {
  count++;
  countALL++;
  if( R_FINITE(value_) )
    mpz_init_set_d(value, value_);
  else
    {mpz_init(value); na = true  ;}
}

/**
 * Construct a biginteger from a string value.
 */
biginteger::biginteger(const std::string& value_) : na(false)
{
  count++;
  countALL++;
  /* mpz_init.. return -1 when error, 0: ok */
  if(mpz_init_set_str(value, value_.c_str(), 0))
    {
      mpz_set_si(value, 0);
      na=true;
    }
  /*	if(mpz_init_set_str(value, value_.c_str(), 0) == -1)
	Rf_error("Not a valid number");    */
}

/**
 *  Copy constructor (mpz_t aren't standard-copyable)
 */
biginteger::biginteger(const biginteger& rhs) : na(rhs.na)
{
  count++;
  countALL++;
  mpz_init_set(value, rhs.getValueTemp());
}


  /**
   * Free the owned mpz_t structs
   */
biginteger::~biginteger() {
  count--;
  //    printf("nb : %d %d \n",count,countALL);
  mpz_clear(value);
}


int biginteger::getCount(){
  return count;
}

biginteger::biginteger(const char* raw)
{
  count++;
  countALL++;
  mpz_init(value);
  na = true;
  const int* r = (int*)(raw);
  if (r[0]>0)
    {
      mpz_import(value, r[0], 1, sizeof(int) , 0, 0, (void*)&(r[2] ));
      if(r[1]==-1)
	mpz_neg(value,value);
      na = false;
    }
  else
    mpz_set_si(value, 0);

  //std::cout << "lu: "<<  this->str(10) << std::endl;
}

  /**
   * Create a biginteger from a value. Remember to free the
   * parameter's mpz_t if you allocated them by yourself -
   * biginteger will copy the value.
   */
biginteger::biginteger(const mpz_t value_)
    : na(false)
{
  count++;
  countALL++;
  mpz_init_set(value, value_);
}

/*
 * Convert to string in base b; b from 2 to 36
 */

string biginteger::str(int b) const
{
  if (isNA())
    return "NA";

  // possible minus sign, size of number + '\0'
  char* buf = new char[mpz_sizeinbase(value, b)+2];
  mpz_get_str(buf, b, value);
  string s = buf;
  delete [] buf;
  return s;
}


/**
 * \brief export mpz to R raw value
 */
int biginteger::as_raw(char* raw) const
{
  int totals = raw_size() ;
  memset(raw, 0, totals );

  int* r = (int*)raw;
  r[0] = totals/sizeof(int) - 2;

  if (!isNA())
    {
      r[1] = (int) mpz_sgn(value);
      mpz_export(&r[2], 0, 1, sizeof(int), 0, 0, value);
    }


  return totals;
}

/**
 * \brief export mpz to R raw value
 *
 * return number of byte used (if na => 2*sizeofint, else
 * the two int + size of the big integer)
 *
 * \note IF size of big integer exceed value of one int => cannot store value
 */
int as_raw(char* raw,mpz_t value, bool na)
{
  int numb = 8*sizeof(int);
  int totals = sizeof(int);
  if(!na)
    totals =  sizeof(int) * (2 + (mpz_sizeinbase(value,2)+numb-1) / numb);
  memset(raw, 0, totals );

  int* r = (int*)raw;
  r[0] = totals/sizeof(int) - 2;

  if (!na)
    {
      r[1] = (int) mpz_sgn(value);
      mpz_export(&r[2], 0, 1, sizeof(int), 0, 0, value);
    }

  return totals;

}


// number of int used in R by on biginteger
size_t biginteger::raw_size() const
{
  if (isNA())
    return sizeof(int);

  int numb = 8*sizeof(int);
  return sizeof(int) * (2 + (mpz_sizeinbase(value,2)+numb-1) / numb);

  //  return (sizeof(int) * ( 2 + (mpz_sizeinbase(value,2)/(8*sizeof(int)))  ) );
}

void biginteger::swap(biginteger& other)
{
  mpz_swap(value, other.value);
  bool n = na;
  na = other.na;
  other.na = n;
}


biginteger & biginteger::operator= (const biginteger& rhs)
{
  if(this != &rhs)
    {
      mpz_set(value,rhs.getValueTemp());
      na = rhs.na;
    }
  return(*this);
}


// comparison
bool operator!=(const biginteger& rhs, const biginteger& lhs)
{
  if(rhs.isNA() || lhs.isNA())
    return(false); // SHOULD RETURN NA

  return(mpz_cmp(rhs.getValueTemp(),lhs.getValueTemp())!=0);
}



// comparison
bool operator<(const biginteger& rhs, const biginteger& lhs)
{
  if(rhs.isNA() || lhs.isNA())
    return(false); // SHOULD RETURN NA

  return(mpz_cmp(rhs.getValueTemp(),lhs.getValueTemp())<0);
}



// comparison
bool operator>(const biginteger& rhs, const biginteger& lhs)
{
  if(rhs.isNA() || lhs.isNA())
    return(false); // SHOULD RETURN NA

  return(mpz_cmp(rhs.getValueTemp(),lhs.getValueTemp())>0);
}



// addition
biginteger operator* (const biginteger& rhs, const biginteger& lhs)
{
  // one of them NA: return NA
  if(rhs.isNA() || lhs.isNA())
    return(biginteger());

  mpz_t result;
  mpz_init(result);
  mpz_t_sentry val_s(result);
  mpz_mul(result,rhs.getValueTemp(),lhs.getValueTemp());
  return(biginteger(result));

}

// substraction
biginteger operator- (const biginteger& rhs, const biginteger& lhs)
{
  // one of them NA: return NA
  if(rhs.isNA() || lhs.isNA())
    return(biginteger());

  mpz_t result;
  mpz_init(result);
  mpz_t_sentry val_s(result);
  mpz_sub(result,rhs.getValueTemp(),lhs.getValueTemp());
  return(biginteger(result));
}

// modulus
biginteger operator% (const biginteger& rhs, const biginteger& lhs)
{
  // one of them NA: return NA
  if(rhs.isNA() || lhs.isNA())
    return(biginteger());

  mpz_t result;
  mpz_init(result);
  mpz_t_sentry val_s(result);
  mpz_mod(result,rhs.getValueTemp(),lhs.getValueTemp());
  return(biginteger(result));
}
