/*! \file bigmod.h
 *  \brief Description of class bigmod
 *
 *  \date Created: 22/05/06
 *  \date Last modified: Time-stamp: <2023-01-24 19:36:06 (antoine)>
 *
 *  \author Immanuel Scholz
 *
 *  \note Licence: GPL (>= 2)
 */


#ifndef BIGMOD_HEADER_
#define BIGMOD_HEADER_ 1

#include <memory>
#include "biginteger.h"

typedef void (*gmp_binary)(mpz_t, const mpz_t, const mpz_t);

extern "C" {
  /**
   * division
   * result = a / b
   */
  void integer_div(mpz_t result,const mpz_t a, const mpz_t b);
}



/**
 * \brief class for bigmod values. Represent any integer in Z/nZ
 *
 * Represents two biginteger: a value and a modulus. These both are used
 * to operate arithmetic functions on it. If the modulus is NA, no modulus
 * to the operation result is applied. If the value is NA, the result is always NA.
 */
class bigmod {
 private:
  /** \brief  Value of our bigmod */
  std::shared_ptr<biginteger> value;
  /** \brief  modulus of our bigmod representation: value %% modulus */
  std::shared_ptr<biginteger> modulus;
 
 public:

  /** keep both pointers value / modulus
   */
  bigmod(biginteger * value_,
	 biginteger * modulus_)  :
    value(value_),
    modulus(modulus_) {};

 /** keep references value / modulus is new object.
   */
 bigmod(biginteger* value_)  :
   value(value_),
   modulus(std::make_shared<biginteger>()) {};


  /**
   * create 2 new objects valus / modulus.
   */
 bigmod(const biginteger& value_,
	 const biginteger& modulus_)  :
   value(std::make_shared<biginteger>(value_)),
   modulus(std::make_shared<biginteger>(modulus_)) {};


 bigmod(const biginteger& value_)  :
   value(std::make_shared<biginteger>(value_)),
   modulus(std::make_shared<biginteger>()) {};


 bigmod()  :
   value(std::make_shared<biginteger>()),
   modulus(std::make_shared<biginteger>()) {};

 bigmod(std::shared_ptr<biginteger> value_, std::shared_ptr<biginteger> modulus_ )  :
    value(),
    modulus()
  {
    value = value_;
    modulus = modulus_;
  };
 
  bigmod(std::shared_ptr<biginteger> value_ )  :
    value(),
    modulus(std::make_shared<biginteger>())
  {
    value = value_;
  };
 
 
  /** \brief copy operator  */
   bigmod(const bigmod & rhs) : 
    value(),
    modulus()
  {
    value = rhs.value;
    modulus = rhs.modulus;
  };
 


  virtual ~bigmod(){
  
  };

  /**
   * \brief  Return as a human readible string
   */
  std::string str(int b) const;

  /** \brief assignement operator */
  bigmod & operator= (const bigmod& rhs);

  /** \brief return sign (-1 if negative, 0 if 0; +1 if positive)
   */
  inline int sgn() const
    {
      return(mpz_sgn(getValue().getValueTemp()));
    }

  bigmod  inv () const;

 
  inline biginteger & getValue() {
    return *value;
  }
  
  inline mpz_t & getMpValue()
  {
    return value->getValue();
  }
  
  inline const mpz_t& getValueTemp() const
  {
    return value->getValueTemp();
  }

  biginteger & getModulus() {
    return *modulus;
  }


  inline bool isNA(){
    return value->isNA();
  }
  
  std::shared_ptr<biginteger> & getValuePtr()  {
    return value;
  }
  
  const std::shared_ptr<biginteger> & getModulusPtr() const {
    return modulus;
  }

  inline void setValue(const std::shared_ptr<biginteger>  & v){
    value = v;
  }
  
  inline void setValue(const bigmod  & v){
    value = std::make_shared<biginteger>(v.getValue());
    modulus = std::make_shared<biginteger>(v.getModulus());
  }

  inline void setValue( bigmod  & v){
    value = v.getValuePtr();
    modulus = v.getModulusPtr();
  }

  inline void setValue(const mpz_t  & v){
    value->setValue(v);
  }

  inline void setValue(int v){
    value->setValue(v);
  }

  inline void setValue(double v){
    value->setValue(v);
  }

  inline void setValue(){
    value->setValue();
    modulus->setValue();
  }

  inline void setModulus(const std::shared_ptr<biginteger> & v){
    modulus = v;
  }
  
  const biginteger & getValue() const{
    return *value;
  }

  const biginteger & getModulus() const {
    return *modulus;
  }

};





/** \brief comparison operator
 */
bool operator!= (const bigmod& rhs, const bigmod& lhs);

/** \brief comparison operator
 */
bool operator== (const bigmod& rhs, const bigmod& lhs);



/**
 * \brief Add two bigmods together.
 *
 * If only one has a modulus set, the result will have this
 * modulus. If both bigmods disagree with the modulus, the result will not have
 * a modulus set. If none modulus for either bigmod is set, the result will not
 * have a modulus as well.
 */
bigmod operator+(const bigmod& rhs, const bigmod& lhs);

/**
 * \brief Subtract two bigmods.
 *
 * For modulus description, see operator+(bigmod, bigmod)
 */
bigmod operator-(const bigmod& rhs, const bigmod& lhs);

/**
 * \brief Multiply two bigmods.
 *
 * For modulus description, see operator+(bigmod, bigmod)
 */
bigmod operator*(const bigmod& rhs, const bigmod& lhs);

/**
 * \brief Divide two bigmods   a / b  :=  a * b^(-1)
 */
bigmod div_via_inv(const bigmod& a, const bigmod& b);

/**
 * \brief Divide two bigmods.
 *
 * For modulus description, see operator+(bigmod, bigmod)
 */
bigmod operator/(const bigmod& rhs, const bigmod& lhs);

/**
 * \brief Calculate the modulus (remainder) of two bigmods.
 *
 * The resulting bigmod will have set the intern modulus to
 * the value of lhs, no matter what rhs.modulus or lhs.modulus
 * was before, except if rhs and lhs has both no modulus set,
 * in which case the resulting modulus will be unset too.
 */
bigmod operator%(const bigmod& rhs, const bigmod& lhs);

/**
 * \brief Return the power of "exp" to the base of "base" (return = base^exp).
 *
 * If both moduli are unset or unequal, this may EAT your memory alive,
 * since then the infinite "pow" is used instead of the modulus "powm".
 * You  may not try to pow a value this way with an exponent that does
 * not fit into a long value.
 *
 * For other modulus description, see operator+(bigmod, bigmod)
 */
bigmod pow(const bigmod& base, const bigmod& exp);

/**
 * \brief Return the modulo inverse to x mod m. (return = x^-1 % m)
 *
 * For modulus description, see operator+(bigmod, bigmod)
 */
bigmod inv(const bigmod& x, const bigmod& m);

/**
 * \brief Return a bigmod with value (x % m) and the intern modulus set to m.
 * Intern modulus settings of x and m are ignored.
 *
 * Do not confuse this with operator%(bigmod, bigmod).
 */
bigmod set_modulus(const bigmod& x, const bigmod& m);


biginteger get_modulus(const bigmod& b1, const bigmod& b2);
/**
 * \brief Return the greatest common divisor of both parameters
 *
 * For modulus description, see operator+(bigmod, bigmod)
 */
bigmod gcd(const bigmod& rhs, const bigmod& lhs);

/**
 * \brief  Return the least common multiply of both parameter.
 *
 * For modulus description, see operator+(bigmod, bigmod)
 */
bigmod lcm(const bigmod& rhs, const bigmod& lhs);

/**
 * \brief function used to make any binary operation between
 * two bigmod that return a bigmod (addition substraction... )
 */
bigmod create_bigmod(const bigmod& lhs, const bigmod& rhs, gmp_binary f,
		     bool zeroRhsAllowed = true) ;

#endif
