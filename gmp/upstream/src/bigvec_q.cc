#include "bigvec_q.h"



bigvec_q::bigvec_q(const bigvec_q & rhs):
  value(rhs.value.size()),
  nrow(0)
{
  *this = rhs;
}


bigvec_q::bigvec_q(const bigvec & rhs):
  value(rhs.size()),
  nrow(rhs.nrow)
{
  for(unsigned int i=0; i< rhs.size(); ++i)
    {
      value[i].setValue(rhs[i].getValue());
    }
}

bigvec_q::~bigvec_q(){
  value.clear();
}


bigvec_q & bigvec_q::operator= (const bigvec_q & rhs)

{
  if(this != &rhs)
    {
      nrow = rhs.nrow;
      value.resize(rhs.value.size());
      std::vector<bigrational>::iterator it = value.begin();
      std::vector<bigrational>::const_iterator jt = rhs.value.begin();
      while(it != value.end())
	{
	  *it = *jt;
	  ++it;
	  ++jt;
	}
    }
  return(*this);

}


const bigrational & bigvec_q::operator[] (unsigned int i) const
{
  return(value[i]);
}

bigrational & bigvec_q::operator[] (unsigned int i) 
{
  return(value[i]);
}

bigrational & bigvec_q::get(unsigned int row, unsigned int col) {
  bigrational & val = (*this)[row + col*nRows()];
  return val;
}


void bigvec_q::set(unsigned int row, unsigned int col, const bigrational & val) {
  set( row + col*nRows(),val);
}



void bigvec_q::set(unsigned int i,const bigrational & val)
{
  value[i] = val;
}
void bigvec_q::set(unsigned int i,const mpq_t & val)
{

  if(i>=value.size())
    {
      Rprintf("ERROR at bigvec_q_set_mpq __LINE__ \n");
      return;
    }

  value[i].setValue(val);

}

void bigvec_q::push_back(const bigrational &i)
{
  value.push_back(i);
}

unsigned int bigvec_q::size() const
{
  return(value.size());
}

unsigned int bigvec_q::nRows() const {
  return abs(nrow);
}

void bigvec_q::resize(unsigned int n)
{
  value.resize(n);
}

void bigvec_q::clear()
{
  value.clear();
  nrow=0;
}




// never used
void bigvec_q::print()
{
  if(nrow>0) {
    for(int i=0; i < nrow; ++i)
      {
	for(unsigned int j=0; j < (value.size() / nrow); ++j)
	  Rprintf("%s\t", value[i+j* nrow].str(10).c_str() );
	Rprintf("\n");
      }
  }
  else {
    for(unsigned int i=0; i < value.size(); ++i)
      Rprintf("%s\t", value[i].str(10).c_str() );
    Rprintf("\n");
  }
}
