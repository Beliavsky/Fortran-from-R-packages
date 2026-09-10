/*! \file bigvec.h
 *  \brief bigvec class definition
 *
 *  \version 1
 *
 *  \date Created: 2005
 *  \date Last modified: Time-stamp: <2023-02-03 08:59:56 (antoine)>
 *
 *
 *  \note Licence: GPL (>= 2)
 */


#ifndef BIGVEC_HEADER_
#define BIGVEC_HEADER_ 1

#include <memory>

#include "bigmod.h"
#include "templateMatrix.h"

enum TYPE_MODULUS{
		  NO_MODULUS,
		  MODULUS_GLOBAL,
		  MODULUS_BY_CELL
};

/** \brief class bigvec
 *
 * It a class composed of 2 vectors, (value & modulus) that
 * can be of different size and a nrow
 * parameter (for matrix support)
 */
class bigvec : public math::Matrix<bigmod> {
  
protected:
  /** array with all bigmod, that are references to values in vector. */
  std::vector<bigmod > values;
  TYPE_MODULUS type;
  std::shared_ptr<biginteger> modulus;
  
public:
  /** \brief optional parameter used with matrix -- set to -1 for non-matrix */
  int nrow ;

  /** \brief initialize value to size i
   */
  bigvec(unsigned int i = 0);

  /**
   * \brief copy constructor
   */
  bigvec(const bigvec & vecteur);

  virtual ~bigvec();

  inline bool isVector() const{
    return nrow < 0 ;
  }

  /**
   * \brief construct a bigmod at indice i
   *
   * It gets values from value[i] & modulus[i % modulus.size()]
   *
   * \note should  not used for assignement
   */
  const bigmod & operator[] (unsigned int i) const;

  bigmod & operator[] (unsigned int i) ;

  /**
   * \brief assign a value at indice i
   */
  void set(unsigned int i,const bigmod & val);

  void set(unsigned int row, unsigned int col, const bigmod & val) ;

  bigmod & get(unsigned int row, unsigned int col) ;

  /**
   * \brief extend our vectors.
   *
   * This function will check if modulus should be added.
   * Modulus can be set "globaly" i.e. one modulus for the
   * whole vector/matrix.
   * Or set by row (a constant modulus for each row)
   * Or set by cell (as many modulus as value)
   */
  void push_back(const bigmod &i);

  /** 
   * insert int value
   */
  void push_back(int value_p);

  inline void erase(unsigned int index){
    values.erase(values.begin() + index);
  }

  /**
   * Insert Big Integer value
   */
  void push_back(biginteger & value_p);
  void push_back(const __mpz_struct* value_p);

  /**
   * \brief return size of vector value
   */
  unsigned int size() const ;

  unsigned int nRows() const;

  /**
   * \brief extend vector value.
   */
  void resize(unsigned int i);

  /**
   * \brief clear all
   */
  void clear();

  /**
   * \brief Return as a human readible string
   */
  std::string str(int i, int b) const;

  /**
   * \brief assignement operator
   */
  bigvec & operator= (const bigvec & rhs);

  /** \brief print matrix to stdout
   *
   * use for debug purpose
   */
  void print();

  inline unsigned int getModulusSize(){
    if(type == NO_MODULUS) return 0;
    if(type == MODULUS_GLOBAL) return 1;
    return size();
  }

  void setGlobalModulus(std::shared_ptr<biginteger> & modulus);

  inline std::shared_ptr<biginteger> & getGlobalModulus(){
    return modulus;
  }
  
  const  TYPE_MODULUS getType() const {
    return type;
  }

  inline void setType(TYPE_MODULUS val){
    type = val;
    if(val == MODULUS_GLOBAL && size() > 0){
      modulus = (*this)[0].getModulusPtr();
    }
  }

  inline void reserve(unsigned int size){
    values.reserve(size);
  }
  /** return a global modulus if possible, null if not */
  static std::shared_ptr<biginteger> getGlobalModulus(bigvec & first,  bigvec & second);
  
 private:
  void checkValuesMod() ;
  void clearValuesMod() ;


};
//


/** \brief comparison operator
 */
bool operator!= (const bigvec& rhs, const bigvec& lhs);



#endif
