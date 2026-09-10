
#include "bigvec.h"

static int count = 0;
static int countAll = 0;


/** \brief constructor
 *
 */
bigvec::bigvec(unsigned int size) :
  math::Matrix<bigmod>(),
  values(),
  type(NO_MODULUS),
  modulus(),
  nrow(-1)
{
  count++;
  countAll++;
 
  for (unsigned int i=0 ; i < size; i++){
    values.push_back(bigmod());
  }
}


bigvec::bigvec(const bigvec & vecteur) :
  math::Matrix<bigmod>(),
  values(),
  type(vecteur.type),
  modulus(vecteur.modulus),
  nrow(vecteur.nrow)
{
  count++;
  countAll++;
   //  *this = vecteur;
  values.reserve(vecteur.size());
  for(std::vector<bigmod>::const_iterator it= vecteur.values.begin(); it !=  vecteur.values.end(); ++it)
    {
      values.push_back(*it);
    }
}

bigvec::~bigvec(){
  count--;
  //  printf("bv %d %d \n",count, countAll);
  clear();
    
}

//
std::string bigvec::str(int i,int b) const
{
   return values[i].str(b);
}

bigmod & bigvec::get(unsigned int row, unsigned int col) {
  return (*this)[row + col*nRows() % size()];
}


 bigmod & bigvec::operator[] (unsigned int i)
{
  return values[i];
}


const bigmod & bigvec::operator[] (unsigned int i) const
{
  return values[i];
}

void bigvec::set(unsigned int row, unsigned int col, const  bigmod & val) {
  set( row + col*nRows(),val);
}


void bigvec::set(unsigned int i,const bigmod & val)
{
  values[i] = val;


  if(type == NO_MODULUS){
    if(val.getModulus().isNA()) return;
    if(i == 0 && values.size() == 1) {
      type = MODULUS_GLOBAL;
      modulus = val.getModulusPtr();
    } else {
      type = MODULUS_BY_CELL;
    }
  }
  if(type == MODULUS_GLOBAL){
    if(values.size() == 1){
      modulus = val.getModulusPtr();
    } else 	if(val.getModulus() != *modulus ){
      type = MODULUS_BY_CELL;
    }
  }
}

void bigvec::push_back(const bigmod & number)
{
  values.push_back(bigmod());
  set(values.size()-1, number);
}

/**
 * insert int value
 */
void bigvec::push_back(int value_p)
{
  push_back(biginteger(value_p));
}

/**
 * insert int value
 */
void bigvec::push_back(biginteger & value_p)
{
  push_back(bigmod(value_p));
}

/**
 * Insert Big Integer value
 */
void bigvec::push_back(const __mpz_struct * value_p)
{
  push_back(biginteger(value_p));
}


// return size of value
unsigned int bigvec::size() const
{
  return(values.size());
}


unsigned int  bigvec::nRows() const {
   return abs(nrow);
}


// hummm. to avoid !
void bigvec::resize(unsigned int i)
{

  values.resize(i);
}

// clear all
void bigvec::clear()
{
  values.clear();
  type=NO_MODULUS;
  modulus=nullptr;
  nrow = -1;
}


// assignment operator
bigvec & bigvec::operator= (const bigvec & rhs)
{
  if(this != &rhs)
    {
      values.resize(0);
      modulus = rhs.modulus;
      type=rhs.type;
      
      for (unsigned int i = 0 ;  i < rhs.size(); i++){
	values.push_back(rhs[i]);
      }
      
      nrow = rhs.nrow;
      
    }
  return(*this);
}



// Comparison operator
bool operator!=(const bigvec & rhs, const bigvec& lhs)
{


  if( (rhs.size() != lhs.size()) || \
      (rhs.nrow != lhs.nrow )  )
    return(false);


  
  // check value
  for (unsigned int i = 0 ; i < rhs.size(); i++)
    {
      if(rhs[i] != lhs[i])
	return(false);
    }

  return(true);
}



// never used
void bigvec::print()
{
  if(nrow > 0) {
    for(int i=0; i < nrow; ++i)
      {
	for(unsigned int j=0; j < (values.size() / nrow); ++j)
	  Rprintf("%s\t", values[i+j* nrow].str(10).c_str() );
	Rprintf("\n");
      }
  }
  else {
    for(unsigned int i=0; i < values.size(); ++i)
      Rprintf("%s\t", values[i].str(10).c_str() );
    Rprintf("\n");
  }
}


void bigvec::setGlobalModulus(std::shared_ptr<biginteger> & val){
  modulus = val;
  type = MODULUS_GLOBAL;
  for (unsigned int i = 0 ; i < values.size(); i++){
    values[i].setModulus(val);
  }
}


std::shared_ptr<biginteger> bigvec::getGlobalModulus(bigvec & first,  bigvec & second){
  std::shared_ptr<biginteger> empty(nullptr);
  if(first.getType() == MODULUS_GLOBAL && second.getType() ==  NO_MODULUS) return first.getGlobalModulus();
  if(first.getType() == NO_MODULUS && second.getType() ==  MODULUS_GLOBAL) return second.getGlobalModulus();
  if(first.getType() == MODULUS_GLOBAL && second.getType() ==  MODULUS_GLOBAL) {
    return (*first.getGlobalModulus()) == (*second.getGlobalModulus()) ? first.getGlobalModulus() : empty;
  }
  return empty;
}
