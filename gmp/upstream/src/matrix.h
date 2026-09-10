/*! \file matrix.h
 *  \brief C++ function to add matrix support
 *
 *  \version 1
 *
 *  \date Created: 25/05/06
 *  \date Last modified: Time-stamp: <2022-12-09 14:56:06 (antoine)>
 *
 *  \author A. Lucas
 *
 * \note
 *  as usually, matrix x[i,j] (n x p) is represented by a vector
 *              x[i + j x n]  (i=0..n-1 ; j=0..p-1)
 *
 *  \note Licence: GPL (>= 2)
 */

#ifndef MATRIXZ_HEADER_
#define MATRIXZ_HEADER_ 1



extern "C"
{
  /** \brief is x a bigz or bigz matrix ? */
  SEXP is_matrix_zq(SEXP x);

  /**
   * \brief build a matrix x with dimensions p&q byrow: 0 or 1
   */
  SEXP as_matrixz(SEXP x, SEXP p, SEXP q, SEXP byrow, SEXP mod);

  /**
   * \brief transpose a "matrix" n x p into the transpose p x n
   */
  SEXP bigint_transposeR(SEXP x);

  /** \brief  matrix cross product */
  SEXP matrix_crossp_z (SEXP a, SEXP trans);
  /** \brief  matrix multiplication */
  SEXP matrix_mul_z (SEXP a, SEXP b, SEXP op);


  /** \brief for function rbind
   */
  SEXP biginteger_rbind(SEXP args) ;

  /** \brief for function rbind
   */
  SEXP biginteger_cbind(SEXP args) ;

}


/** \brief  C functions for matrix */
namespace matrixz{

  /** \brief C function use to transpose a matrix */
  bigvec bigint_transpose ( bigvec & mat);

  /** \brief Check dimension compatibility */
  int checkDims(int  dima,int  dimb);
}


#endif
