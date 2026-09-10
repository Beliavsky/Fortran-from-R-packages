/*! \file extract_matrix.h
 *  \brief functions to extract parts of matrix
 *
 *  \version 1
 *
 *  \date Created: 25/06/11
 *  \date Modifications: see cvs log
 *
 *  \author A. Lucas
 *
 * \note
 *  as usually, matrix x[i,j] (n x p) is represented by a vector
 *              x[i + j x n]  (i=0..n-1 ; j=0..p-1)
 *
 *  \note Licence: GPL (>= 2)
 */


#ifndef EXTRACT_MATRIX_HEADER_GMP_R_
#define EXTRACT_MATRIX_HEADER_GMP_R_ 1

#include <functional>
#include <algorithm>

#include "bigvec_q.h"
#include "bigintegerR.h"
#include <stdexcept>

#include <stdio.h>

extern "C"
{

  /** \brief get subsets of a matrix */
  SEXP matrix_get_at_z(SEXP A,SEXP INDI, SEXP INDJ);

  /** \brief set subsets of a matrix */
  SEXP matrix_set_at_z(SEXP A,SEXP VAL ,SEXP INDI, SEXP INDJ);

  /** \brief get subsets of a matrix */
  SEXP matrix_get_at_q(SEXP A, SEXP INDI, SEXP INDJ);

  /** \brief set subsets of a matrix */
  SEXP matrix_set_at_q(SEXP A,SEXP VAL, SEXP INDI, SEXP INDJ);

}


namespace extract_gmp_R
{

  /** \brief Change R indices (in function like A[IND]) from
   * R (it can be logical/positive/negative) to a vector of size n
   * containing boolean: true: we update line/row, flase : we do not
   * change anything
   *
   * \return a vector of n boolean corresponding to values that must be affected.
   *
   */
  std::vector<int> indice_get_at (unsigned int n , SEXP & IND);


  /**
   * \brief tranform a matrix from bigvec or bigvec_q format to
   * vector<vector< big > > (a vector of columns)
   *
   * \note for bigvec: it does not take into account modulus.
   *
   */
  template< class T> void  toVecVec(T& A, std::vector<T * > & retour )
  {
    //typename std::vector<T> retour;
    unsigned int i;

    // case: vector
    if(A.nrow < 0)
      A.nrow = A.size();
    else // check that size is a multiple of row
      {
	if((A.size() / A.nrow) != static_cast<float>(A.size()) / static_cast<float>(A.nrow)){
	  A.clear();
	  Rf_error("malformed matrix");
	}
      }
    
    retour.resize(A.size() / A.nrow);
    for(i = 0; i < retour.size();  ++i)
      {
	retour[i] = new T();
	retour[i]->resize(A.nrow);
      }
    // go !
    for(i= 0 ; i < A.size(); ++i)
      // retour[col       ]  [row        ]
      (*retour[i / A.nrow ])[ i % A.nrow].setValue(A[i]);

    //return(retour);

  }




  template< class T> void clearVec(typename std::vector<T*> & vec )
    {
      for (typename std::vector<T*>::const_iterator it = vec.begin();
	   it != vec.end();
	   ++it)
	delete *it;
    }


  /**
   * \brief extract parts of a matrix
   *
   * \param A matrix (class bigvec or bigvec_q)
   * \param INDI,INDJ indices: either "LOGICAL" (true/false) or
   *  numbers:
   *      - positives: we return row/col in INDI/INDJ
   *      - negatives: we retun all except row/col in INDI/INJ
   */
  template <class T> T get_at (T & A, SEXP& INDI, SEXP& INDJ)
  {
    //DEBUG
    //return A;
    // result = A[indi,indj]
    int oldnrow = A.nrow;
    std::vector<int> indJ;

    typename std::vector<T*> matricetmp ;
    typename std::vector<T*> matricetmp2;
    
    toVecVec(A,matricetmp);
    typename std::vector<T*> copyAdress(matricetmp);

    // only pointers
    typename std::vector<T*> * matrice = &matricetmp;
    typename std::vector<T*> * matricework = &matricetmp2;

    T retour;

    unsigned int i,j, newnrow;
    std::vector<int>::iterator it;

    // --------------------------
    // PART 1:  COLUMNS
    // --------------------------

    if(A.size()==0)
      {
	clearVec<T>(copyAdress);
	return(retour);
      }

    if(TYPEOF(INDJ) != NILSXP) {
      indJ = bigintegerR::create_int(INDJ);

      if (TYPEOF(INDJ) == LGLSXP) // LOGICAL
	{

	  // for all columns
	  unsigned int delta=0;
	  for (i = 0; i < (*matrice)[0]->size(); ++i)
	    {
	      if (! indJ[i+delta% indJ.size()])
		{
		  // remove columns i
		  matrice->erase(i+ matrice->begin());
		  --i; // indice in modified matrix
		  ++delta; // i+delta: old indice
		}
	    }

	}
      else //INDJ: numbers
	{
	  indJ.erase(std::remove(indJ.begin(), indJ.end(), 0L), indJ.end()); // remove all zeroes

	  if (indJ.empty())
	    {
	      clearVec<T>(copyAdress);
	      return retour;
	    }

	  // case: a[-b]
	  // negatives...
	  if(indJ[0] < 0)
	    {
	      // sort in reverse order
	      std::sort(indJ.begin(),indJ.end(),std::greater<int>() );

	      // we should have indJ like -1 -3 -7 -7 -12 ...

	      // remove consecutive  duplicates
	      it = std::unique(indJ.begin(),indJ.end());
	      indJ.erase(it,indJ.end());

	      if ( indJ.back() > 0){
		clearVec<T>(copyAdress);
		A.clear();
		Rf_error("only 0's may mix with negative subscripts");
	      }


	      it=indJ.begin();
	      unsigned int delta=0;
	      // for all columns
	      for (j = 0; j < matrice->size(); ++j)
		{
		  if(it == indJ.end())
		    break;

		  if (*it == - static_cast<int>(j+1+delta) )
		    {
		      matrice->erase(j+ matrice->begin());
		      ++it;
		      ++delta;
		      --j;
		    }

		}

	    }
	  else
	    // case: positive values: 1;5;7;7;9;10...
	    {
	      // note : can have duplicates and indices are not sorted

	      // allocate new matrix (all will be copied)
	      // number of columns
	      matricework->reserve(indJ.size());

	      // for all [new] rows
	      for( it=indJ.begin(); it != indJ.end(); it++)
		{
		  if (*it < 0){
		    clearVec<T>(copyAdress);
		    A.clear();		    
		    Rf_error("only 0's may mix with negative subscripts");
		  }
		  if( static_cast<unsigned int>(*it-1) < matrice->size() ) {
		      //printf("on sort %s",(*matrice)[(*it)-1][0].str(10).c_str());
		      matricework->push_back( (*matrice)[*it-1] );
		  } else {
		    clearVec<T>(copyAdress);
		    A.clear();
		    Rf_error("column subscript out of bounds");
		  }
		}

	      // change addresses
	      matrice = &matricetmp2;
	      matricework = &matricetmp;

	    }//end of case: int+positive values

	}
      } // INDJ "exists"

    if(matrice->size()==0)
      {
	clearVec<T>(copyAdress);
	return(retour);
      }

    // --------------------------
    // PART 2:  ROWS
    // --------------------------

    std::vector<int> indI = bigintegerR::create_int(INDI);
    if(TYPEOF(INDI) != NILSXP) {
	if (TYPEOF(INDI) == LGLSXP) { // LOGICAL
	  // for all rows
	  unsigned int delta = 0;
	  for (i = 0; i < (*matrice)[0]->size(); ++i)
	    {
	      if (! indI[(i+delta)% indI.size()])
		{
		  // for all columns j delete row i
		  for (j = 0; j < matrice->size(); ++j)
		    (*matrice)[j]->erase(i);

		  //++newnrow;
		  --i; // i: new indice in modified matrix
		  ++delta; // i+delta = old indices
		}
	    }

	}
	else { // INDI : numbers
	  // remove zeros:
	  indI.erase(std::remove(indI.begin(), indI.end(), 0L), indI.end());

	  if (indI.empty())
	    {
	      clearVec<T>(copyAdress);
	      return retour;
	    }

	  // case: a[-b]
	  // negatives...
	  if(indI[0] < 0)
	    {
	      std::sort(indI.begin(),indI.end(),std::greater<int>() );
	      // we should have indI like -1 -3 -7 -7 -12 ...

	      // remove duplicates
	      it= std::unique(indI.begin(),indI.end());
	      indI.erase(it,indI.end());

	      if ( indI.back() > 0){
		clearVec<T>(copyAdress);
		A.clear();	
		Rf_error("only 0's may mix with negative subscripts");
	      }


	      //newnrow = A.nrow;
	      it=indI.begin();
	      // for all rows
	      unsigned int delta = 0;
	      for (i = 0; i < (*matrice)[0]->size(); ++i)
		{

		  if(it != indI.end() )
		    if (*it == - static_cast<int>(i+1+delta) )
		      {
			// for all columns j remove row i
			for (j = 0; j < matrice->size(); ++j)
			  {
			    (*matrice)[j]->erase(i);
			  }
			//--newnrow;
			--i; // i: new indice in modified matrix
			++delta; // i+delta = old indices
			++it;
		      }

		}

	    }
	  else
	    {
	    // case: positive values: 1;5;7;7;9;10...
	    
	      // delete too large values or give error iff  INDJ is non-NULL
	      for(it = indI.begin(); it != indI.end(); ++it) 
		{
		  if(*it > static_cast<int>((*matrice)[0]->size())) 
		    {
		      if(oldnrow < 0) { // vector case: out-of-bound --> NA
			/* it = indI.erase(it); */
			/* --it; */
		      } else { // matrix case:
			clearVec<T>(copyAdress);
			A.clear();
			Rf_error("subscript out of bounds");
		      }
		    }
		}
	      // note : can have duplicates and indices are not sorted

	      //newnrow = indI.size();

	      // allocate new matrix (all will be copied)
	      // number of columns

	      matricework->resize( matrice->size());
	      for (typename std::vector<T*>::iterator it = matricework->begin();
		   it != matricework->end();
		   ++it)
		{
		  *it = new T();
		  copyAdress.push_back(*it);
		}

	      // number of row
	      for (j = 0; j < matricework->size(); ++j)
		(*matricework)[j]->resize( indI.size() );


	      // for all [new] rows
	      for( i=0; i < indI.size(); ++i)
		{
		  if (indI[i] < 0){
		    clearVec<T>(copyAdress);
		    A.clear();
		    Rf_error("only 0's may mix with negative subscripts");
		  }
		  if( static_cast<unsigned int>(indI[i]-1) < (*matrice)[0]->size() )
		    {
		      // for all columns
		      for (j = 0; j < matricework->size(); ++j)
			//newmat[i,j] = oldmat[ind[i],j]
			( *(*matricework)[j])[i] = (*(*matrice)[j])[indI[i]-1];
		    }
		  else
		    for (j = 0; j < matricework->size(); ++j)
		      ( *(*matricework)[j])[i].setValue();
		}

	      matrice = matricework; // change addresses
	    }//end of case: int+positive values

	}
      }

    // --------------------------
    // PART 3:  COPY RESULT
    // --------------------------

    newnrow = (*matrice)[0]->size();
    retour.resize(matrice->size() * newnrow);
    for(i=0; i< newnrow ; ++i)
      for(j=0; j <  matrice->size() ; ++j)
	retour[i + j* newnrow ]  =  (*(*matrice)[j])[i] ;

    retour.nrow = (oldnrow < 0) ? -1 : newnrow;

    clearVec<T>(copyAdress);
    return(retour);

  }  // end get_at


  /** \brief set a matrix: for R function src[idx,jdx] <- value
   *
   */
  template<class T> void set_at(T & src ,const T & value, SEXP & IDX, SEXP & JDX)
  {
    // case: vector
    if(src.nrow < 0)
      src.nrow = src.size();

    // check that size is a multiple of row
    if((src.size() / src.nrow) != static_cast<float>(src.size()) / static_cast<float>(src.nrow)){
      src.clear();
      T & vval = const_cast<T&>(value);
      vval.clear();
      throw std::invalid_argument("malformed matrix");
      //Rf_error("malformed matrix");
    }
    unsigned int ncol = src.size() / src.nrow; // number of col
    std::vector<int> vidx =  indice_get_at ( src.nrow, IDX);
    std::vector<int> vjdx =  indice_get_at ( ncol, JDX);
   
    unsigned int k=0;

    for(unsigned int j = 0 ; j < vjdx.size() ; ++j)
      {
	
	for(unsigned int i = 0 ; i < vidx.size(); ++i)
	  {
	    unsigned int index = vidx[i] + vjdx[j] * src.nrow;
	    if(index >= src.size()) {
	      src.clear();
	      T & vval = const_cast<T&>(value);
	      vval.clear();
	      throw std::invalid_argument("indice out of bounds");

	      //	      Rf_error("indice out of bounds");
	    }

	    src.set(index, value[k % value.size()] );
	    ++k;
	  }
      }

    return;

  }//end set_at

}// end namespace






#endif
