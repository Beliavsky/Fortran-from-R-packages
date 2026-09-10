/*! \file biginteger.h
 *  \brief Description of class biginteger
 *
 *  \date Created: 2004
 *  \date Last modified: Time-stamp: <2023-01-28 15:50:33 (antoine)>
 *
 *  \author Immanuel Scholz
 *
 *  \note Licence: GPL (>= 2)
 */

#ifndef biginteger_HEADER
#define biginteger_HEADER 1

#include <string>
#include <vector>

#include "Rgmp.h"

/** \mainpage Gmp package for R language.
 *
 * Theses pages made to help developpers to enter into the C++ code related to
 * the gmp package.
 *
 * \section sec_intro Introduction
 *
 * The gmp R package uses directly the gmp C/C++ library (http://gmplib.org)
 * to provide powerful computation of big integers and rational within R.
 * This package is more than a simple wrapper as it allows to handle
 *    - easy computation in Z/nZ (it is an option of bigz class)
 *    - matrix computation
 *
 */


/** \brief mpz struct
 *
 * Use this to clear mpz_t structs at end-of-function automatically
 */
struct mpz_t_sentry {
  /** \brief value */
  mpz_t& value;
  /** \brief constructor */
  mpz_t_sentry(mpz_t& v): value(v) {}
  /** \brief detructor (and clear mpz) */
  ~mpz_t_sentry() {mpz_clear(value);}
};



/** \brief Class biginteger
 *
 * A big integer. Actually a wrapper for mpz_t to work with plus
 * some special stuff.
 *
 * The biginteger special state "NA" means, no value is assigned.
 * This does not mean, the internal state is not constructed, but
 * the value explicit is "not available".
 */
class biginteger
{
 private:
  /**
   * The actual integer value.
   */
  mpz_t value;

  /**
   * True, if the value is "NA".
   */
  bool na;

 public:
  /**
   * Construct a "NA" biginteger.
   */
  biginteger() ;

  /**
   * Construct a biginteger from a raw expression.
   */
  biginteger(const char* raw);

  /**
   * Create a biginteger from a value. Remember to free the
   * parameter's mpz_t if you allocated them by yourself -
   * biginteger will copy the value.
   */
  biginteger(const mpz_t value_);

  /**
   * Construct a biginteger from a int value.
   */
  biginteger(const int value_);

  /**
   * Construct a biginteger from a long value.
   */
  biginteger(const long int value_);

  /**
   * Construct a biginteger from a unsigned long value.
   */
  biginteger(const unsigned long int value_);

  /**
   * Construct a biginteger from a double value.
   */
  biginteger(const double value_);

  /**
   * Construct a biginteger from a string value.
   */
  biginteger(const std::string& value_);

  /**
   *  Copy constructor (mpz_t aren't standard-copyable)
   */
  biginteger(const biginteger& rhs) ;


  /**
   * Free the owned mpz_t structs
   */
  virtual ~biginteger();


  /**
   * Set the biginteger to state "NA".
   */
  inline void setValue() {mpz_set_si(value, 0); na = true;}

  /**
   * Set the biginteger to a specific value.
   */
  inline void setValue(const mpz_t & value_ ) {
    mpz_set(value, value_); na = false;
  }

  /** \brief set value from an integer
   */
  inline void setValue(int value_) {
    if(value_ == NA_INTEGER)
      {mpz_set_ui(value, 0); na = true  ;}
    else
      {
	mpz_set_si(value, value_);
	na = false;
      }
  }
  /** \brief set value from a long integer
   */
  inline void setValue(long int value_) {
    if(value_ == NA_INTEGER)
      {mpz_set_ui(value, 0); na = true  ;}
    else
      {
	mpz_set_si(value, value_);
	na = false;
      }
  }

  /** \brief set value from an unsigned int
   */
  inline void setValue(unsigned long int value_) {
    if((int)value_ == NA_INTEGER)
      {mpz_set_ui(value, 0); na = true  ;}
    else
      {
	mpz_set_ui(value, value_);
	na = false;
      }
  }

  /** \brief set value from a float
   */
  inline void setValue(double value_) {
    if(R_FINITE (value_) )
      {mpz_set_d(value, value_); na = false;}
    else
      {mpz_set_ui(value, 0); na = true  ;}
  }


  /** \brief set value from a biginteger
   */
  void setValue(const biginteger  & value_) {
    setValue(value_.getValueTemp());
    na = value_.isNA();

  }
  /**
   * For const-purposes, return the value. Remember, that the return value
   * only lives as long as this class live, so do not call getValueTemp on
   * temporary objects.
   */
  inline const mpz_t& getValueTemp() const {return value;}


  /** \brief accessor on value
   */
  inline mpz_t & getValue()
  {
    return value;
  }

  /**
   * Return true, if the value is NA.
   */
  inline bool isNA() const {return na;}

  /**
   * set NA value
   */
  inline void NA(bool value_p)  {na = value_p;}

  /**
   * Return true, if the value is 0.
   */
  inline int sgn() const {return mpz_sgn(value);}

  /**
   *  Convert the biginteger into a standard string.
   */
  std::string str(int b) const;

  /**
   * Convert the biginteger into a long value (cut off the msb's if it don't
   * fit).
   */
  inline long as_long() const {return mpz_get_si(value);}

  /**
   * \brief Convert the biginteger into a double value
   * (you may loose precision)
   */
  inline double as_double() const {return mpz_get_d(value);}

  /**
   * Convert the biginteger to a raw memory block. Obtain the size needed
   * from biginteger_raw_size() first and make sure, the buffer provided is
   * large enough to hold the data.
   *
   * Also remember, that the modulus is not saved this way. To obtain a
   * modulus raw byte use get_modulus().as_raw(void*).
   *
   * @return number of bytes used (same as raw_size())
   */
  int as_raw(char* raw) const;

  /**
   * Return the number of bytes needed to store this biginteger in a
   * continous memory block.
   */
  size_t raw_size() const;

  /**
   * Swap values with the argument
   */
  void swap(biginteger& other);

  /**
   * Test prime numbers
   */
  int isprime(int reps) const {return  mpz_probab_prime_p(value,reps);} 


  /** \brief overload affectation operator
   */
  biginteger & operator= (const biginteger& rhs);


  static int getCount();
};



/** \brief comparison operator
 */
bool operator!= (const biginteger& rhs, const biginteger& lhs);
inline bool operator== (const biginteger& rhs, const biginteger& lhs){
  return !(rhs != lhs);
}


/** \brief comparison operator
 */
bool operator> (const biginteger& rhs, const biginteger& lhs);

/** \brief comparison operator
 */
bool operator< (const biginteger& rhs, const biginteger& lhs);

/** \brief standard operator */
biginteger operator* (const biginteger& rhs, const biginteger& lhs);

/** \brief standard operator */
biginteger operator- (const biginteger& rhs, const biginteger& lhs);

/** \brief standard operator */
biginteger operator% (const biginteger& rhs, const biginteger& lhs);

/** \brief function used by rational that should export
 *   mpz value
 */
int as_raw(char* raw,mpz_t value,bool na);



  





#endif
