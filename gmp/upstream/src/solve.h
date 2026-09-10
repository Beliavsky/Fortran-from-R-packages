/*! \file solve.h
 *  \brief functions to solve matrix
 *
 *  \version 1
 *
 *  \date Created: 25/05/06
 *  \date Last modified: $Date: 2012-01-12 17:46:59 $
 *
 *  \author A. Lucas
 *
 * \note
 *  as usually, matrix x[i,j] (n x p) is represented by a vector
 *              x[i + j x n]  (i=0..n-1 ; j=0..p-1)
 *
 *  \note Licence: GPL (>= 2)
 */


#ifndef SOLVE_HEADER_GMP_R_
#define SOLVE_HEADER_GMP_R_ 1
#include <stdexcept>
#include <R.h>
#include <Rinternals.h>


#include "bigvec_q.h"

extern "C"
{

  /** \brief inverse a rational matrix */
  SEXP inverse_q(SEXP A);

  /** \brief solve matrix system AX=B */
  SEXP solve_q(SEXP A,SEXP B);

  /** \brief inverse a matrix in Z/nZ */
  SEXP inverse_z(SEXP A);

  /** \brief solve matrix system AX=B in Z/nZ*/
  SEXP solve_z(SEXP A,SEXP B);

}


namespace solve_gmp_R
{

  /** \brief solve matrix
   *
   * solve A X = B  (return X) with A & B matrix (of rational or biginteger)
   *
   * A is of dimension nxn X nxm and B nxm (X will be return a B address)
   * We use the Gauss algorithm
   */
  template< class T> void solve (math::Matrix<T> & A , math::Matrix<T> & B)
    {

      // A [ i ,j] = A[ i + j * A.nrow]
      for(unsigned int k = 0 ; k < A.nRows(); ++k)
	{
	  if(A.get(k, k).sgn() == 0 ){
	    A.clear();
	    B.clear();
	    throw std::invalid_argument("System is singular");
	  }
	  
	  // l_k <- (1/akk) l_k
	  T tmpValeur =A.get(k , k).inv() ;
	  A.mulLine(k,tmpValeur);
	  B.mulLine(k,tmpValeur);

	  for(unsigned int i = 0; i < A.nRows(); ++i)
	    {
	      if(i == k)
		continue;
	      // l_i <- l_i - a_ik l_k
	      tmpValeur= A.get(i, k) ;
	      A.subLine(i,k, tmpValeur) ;
	      B.subLine(i,k, tmpValeur) ;
	    }
	}

    }


  /** \brief inverse a rational matrix
   *
   */
  SEXP inverse_q(bigvec_q A);

  /** \brief solve matrix equation AX = B
   */
  SEXP solve_q(bigvec_q A,bigvec_q B);

}






#endif
