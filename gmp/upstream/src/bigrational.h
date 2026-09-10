/*! \file bigrational.h
 *  \brief Description of class bigrational
 *
 *  \date Created: 22/05/06
 *  \date Last modified: Time-stamp: <2023-01-24 18:42:58 (antoine)>
 *
 *  \author Antoine Lucas (adapted from biginteger class made by
 *                         Immanuel Scholz)
 *
 *  \note Licence: GPL (>= 2)
 */

#ifndef BIGRATIONAL_HEADER_
#define BIGRATIONAL_HEADER_ 1

#include <string>

#include "biginteger.h"

/**
 * \brief Class for rational value.
 *
 * A big rational. Actually a wrapper for mpq_t to work with plus
 * some special stuff.
 *
 * The bigrational special state "NA" means, no value is assigned.
 * This does not mean, the internal state is not constructed, but
 * the value explicit is "not available".
 */
class bigrational
{

 private:
  /**
   * The rational = n/d .
   */
  mpq_t value ;

  /**
   * \brief a flag: true => NA value.
   */
  bool na;

 public:
  /**
   * Construct a "NA" bigrational.
   */
  bigrational();

  /**
   * Construct a bigrational from a raw expression.
   * in fact in raw: we will store an mpz_export of
   * either numerator, or denominator
   */
  bigrational(void* raw);

  /**
   * Create a bigrational from a value. Remember to free the
   * parameter's mpz_t if you allocated them by yourself -
   * biginteger will copy the value.
   */
  bigrational(const mpq_t& value_);
  
  /**
   * \brief create a rational from an [big] integer
   */
  bigrational(const mpz_t& value_);

  /**
   * Construct a bigrational from a long value.
   */
  bigrational(int value_);

  /**
   * Construct a bigrational from a long value.
   */
  bigrational(int num_, int den_);

  /**
   * Construct a bigrational from a double value.
   */
  bigrational(double value_);

  /**
   * Construct a bigrational from a string value. it can be "4343" or "2322/4343"
   */
  bigrational(const std::string& value_);

  /**
   *  Copy constructor (mpz_t aren't standard-copyable)
   */
  bigrational(const biginteger & rhs) ;

  /**
   *  Copy constructor (mpz_t aren't standard-copyable)
   */
  bigrational(const bigrational & rhs);

  /**
   * Free the owned mpz_t structs
   */
  virtual ~bigrational();

  /** \brief overload affectation operator
   *
   */
  bigrational & operator= (const bigrational& rhs);

  /**
   * Set the bigrational to state "NA".
   */
  void setValue() {mpq_set_si(value, 0,1); na = true;}

  /**
   * Set the bigrational to a specific value.
   */
  void setValue(const mpq_t value_ ) {
    mpq_set(value, value_); na = false;
  }

  /** \brief Set Denominator: return value = value / value_ */
  void setDenValue(const mpq_t value_ ) {
    if(!na)
      mpq_div(value,value, value_);
  }

  /**
   * \brief set value from integer
   */
  void setValue(int value_) {
    if(value_ == NA_INTEGER)
      {mpq_set_ui(value, 0,1); na = true  ;}
    else
      {
	mpq_set_si(value, value_,1);
	na = false;
      }
  }

  /** \brief set value from unsigned int
   */
  void setValue(unsigned long int value_) {
    if((int)value_ == NA_INTEGER)
      {mpq_set_ui(value, 0,1); na = true  ;}
    else
      {mpq_set_ui(value, value_,1); na = false;}
  }

  /** \brief set value from float (double)
   */
  void setValue(double value_) {
    if( R_FINITE (value_ ) )
      {mpq_set_d(value, value_); na = false;}
    else
      {mpq_set_ui(value, 0,1); na = true  ;}
  }


  /**
   * \brief set value from big integer
   */
  void setValue(const mpz_t & value_) {
    mpq_set_z(value,value_);
    na = false;
  }

  /**
   * \brief set value from biginteger value
   */
  void setValue(const biginteger & value_)
    {
      mpq_set_z(value,value_.getValueTemp());
      na = value_.isNA();
    }

  /**
   * \brief set value from biginteger value
   */
  void setValue(const bigrational & value_)
    {
      mpq_set(value,value_.getValueTemp());
      na = value_.isNA();
    }

  /**
   * For const-purposes, return the value. Remember, that the return value
   * only lives as long as this class live, so do not call getValueTemp on
   * temporary objects.
   */
  const mpq_t& getValueTemp() const {return value;}

  /**
   * Return true, if the value is NA.
   */
  bool isNA() const {return na;}

  /**
   * set NA value
   */
  void NA(bool value_p)  {na = value_p;}

  /**
   * Return 1, if the value is > 0;  -1 when negative, 0 when 0.
   */
  int sgn() const {return mpq_sgn(value);}

  /**
   *  Convert the bigrational into a standard string.
   */
  std::string str(int b) const;

  /**
   * \brief Convert the bigrational into a double value
   * (you may loose precision)
   */
  double as_double() const {return mpq_get_d(value);}

  /**
   * Convert the bigrational to a raw memory block. Obtain the size needed
   * from biginteger_raw_size() first and make sure, the buffer provided is
   * large enough to hold the data.
   *
   * Also remember, that the modulus is not saved this way. To obtain a
   * modulus raw byte use get_modulus().as_raw(void*).
   *
   * @return number of bytes used (same as raw_size())
   */
  //  int as_raw(void* raw) const;

  /**
   * Return the number of bytes needed to store the numerator only in a
   * continous memory block.
   */
  size_t raw_size() const;

  /** \brief  simplify n/d (use of mpq_canonical)
   *
   */
  void  simplify ();

  /**
   * \brief  return 1/x
   */
  bigrational inv();


};


/** \brief mpq struct.
 *
 * Use this to clear mpq_t structs at end-of-function automatically
 */
struct mpq_t_sentry {
  /** the bigrational */
  mpq_t& value;
  /** initialisation */
  mpq_t_sentry(mpq_t& v): value(v) {}
  /** free all */
  ~mpq_t_sentry()
  {
    mpq_clear(value);
  }
};


/**
 * \brief set denominator to a bigrational (in fact divide
 * bigrational by value m.
 */
bigrational set_denominator(const bigrational& x, const bigrational& m);


/**
 * \brief Add two bigrationals together.
 *
 */
bigrational operator+(const bigrational& rhs, const bigrational& lhs);

/**
 * \brief Subtract two bigrationals.
 */
bigrational operator-(const bigrational& rhs, const bigrational& lhs);

/**
 * \brief Multiply two bigrationals.
 */
bigrational operator*(const bigrational& rhs, const bigrational& lhs);

/**
 *\brief  Divide two bigrationals.
 */
bigrational operator/(const bigrational& rhs, const bigrational& lhs);

/**
 *\brief  x^y = pow(x,y) :  bigrational ^ biginteger
 */
bigrational operator^(const bigrational& rhs, const biginteger& lhs);


/** \brief comparison operator
 */
bool operator> (const bigrational& rhs, const bigrational& lhs);

/** \brief comparison operator
 */
bool operator< (const bigrational& rhs, const bigrational& lhs);


/** \brief comparison operator
 */
bool operator!= (const bigrational& rhs, const bigrational& lhs);

#endif
