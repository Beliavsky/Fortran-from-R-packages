/************************************************************/
/*! \file matrix.cc
 *  \brief C++ function to add matrix support
 *
 *  \version 1
 *
 *  \date Created: 19/02/06
 *  \date Last modified: Time-stamp: <2023-02-03 12:03:33 (antoine)>
 *
 *  \author A. Lucas
 *
 * \note
 *  as usually, matrix x[i,j] (n x p) is represented by a vector
 *              x[i + j x n]  (i=0..n-1 ; j=0..p-1)
 *
 *  \note Licence: GPL (>= 2)
 */

#include <vector>
#include <stdexcept>

using namespace std;

#include "bigintegerR.h"
#include "matrix.h"
// need to call matrix_mul_q()
#include "matrixq.h"

// given that x is "bigz" or "bigq",
// return TRUE if x is a bigz/q *matrix*: R's is.matrixZQ(.)
SEXP is_matrix_zq(SEXP x) {
  SEXP nrowSexp = Rf_mkString("nrow");
  PROTECT(nrowSexp);
  SEXP attributeRow = Rf_getAttrib(x,nrowSexp );
  PROTECT(attributeRow);
  SEXP ans = Rf_ScalarLogical(attributeRow != R_NilValue);
  UNPROTECT(2);
  return ans;
}

// C++ side of R function matrix.bigz()
SEXP as_matrixz (SEXP x, SEXP nrR, SEXP ncR, SEXP byrowR, SEXP mod)
{
  try{
    int
      nc=INTEGER(ncR)[0],
      nr=INTEGER(nrR)[0],
      byrow=INTEGER(byrowR)[0];

    // get "bigz" vector, this makes all conversion int to bigz etc...
    bigvec mat = bigintegerR::create_bignum(x);
    int lendat = mat.size();
    // int sizemod = mat.modulus.size();
    // when modulus specified
    bigvec modulus = bigintegerR::create_bignum(mod);

    // A copy from R project to get correct dimension
    // all warnings...
    //
    if (nr == NA_INTEGER){ // This is < 0
      throw invalid_argument(_("matrix: invalid 'nrow' value (too large or NA)"));
    }
    if (nr < 0){
      throw invalid_argument(_("matrix: invalid 'nrow' value (< 0)"));
    }
    if (nc < 0){
      throw invalid_argument(_("matrix: invalid 'ncol' value (< 0)"));
    }
    if (nc == NA_INTEGER){
      throw invalid_argument(_("matrix: invalid 'ncol' value (too large or NA)"));
    }

    if(lendat > 0 ) {
      if (lendat > 1 && (nr * nc) % lendat != 0) {
	if (((lendat > nr) && (lendat / nr) * nr != lendat) ||
	    ((lendat < nr) && (nr / lendat) * lendat != nr))
	  Rf_warning("data length [%d] is not a sub-multiple or multiple of the number of rows [%d] in matrix", lendat, nr);
	else if (((lendat > nc) && (lendat / nc) * nc != lendat) ||
		 ((lendat < nc) && (nc / lendat) * lendat != nc))
	  Rf_warning("data length [%d] is not a sub-multiple or multiple of the number of columns [%d] in matrix", lendat, nc);
      }
      else if ((lendat > 1) && (nr * nc == 0)){
	Rf_warning("data length exceeds size of matrix");
      }
    }

    // update dimension parameters
    if(nr == 1)
      nr = (int)ceil(lendat / (double) nc);
    if(nc == 1)
      nc = (int)ceil(lendat / (double)nr);

    // when we extend "x"
    if(nc*nr > lendat)
      {
	mat.resize(nr*nc);
	for(int i = lendat; i < nr*nc; i++)
	  mat[i] = mat[i % lendat];
      }
    mat.nrow = nr;


  
    if(modulus.size()>0) {// should be allways the case
      if(modulus[0].isNA()) {
	// nothing but mat could have already a modulus
      } else {
	for (unsigned int i = 0; i < mat.size(); i++){
	  mat[i].setModulus(modulus[i % modulus.size()].getValuePtr());
	}
	if(modulus.size() == 1){
	  mat.setGlobalModulus(modulus[0].getValuePtr());
	}
	mat.setType( modulus.size() == 1 ? TYPE_MODULUS::MODULUS_GLOBAL : TYPE_MODULUS::MODULUS_BY_CELL);
	// sizemod = modulus.size();
      }
    }
  
    if(byrow)
      {
	bigvec mat2 = matrixz::bigint_transpose (mat);
	return( bigintegerR::create_SEXP (mat2));
      }

    return( bigintegerR::create_SEXP (mat));
  }
  catch(std::invalid_argument & e){
    Rf_error("%s",e.what());
    return Rf_mkString(0);
  }
}


/*
 * Transposition
 */
SEXP bigint_transposeR(SEXP x)
{
  try{
    SEXP dimKey =Rf_mkString("nrow");
    PROTECT(dimKey);
    SEXP dimAttr = Rf_getAttrib(x,dimKey );
    PROTECT(dimAttr);
    bigvec mat = bigintegerR::create_bignum(x);
    int nr, n = mat.size();

    if (dimAttr == R_NilValue) { // vector
      nr = n;
    } else if (TYPEOF(dimAttr) == INTSXP) {
      nr = INTEGER(dimAttr)[0];
    } else { nr = -1;// -Wall
      mat.clear();
      throw invalid_argument(_("argument must be a matrix of class \"bigz\""));
    }
    UNPROTECT(2);

    mat.nrow = nr;
    // Rprintf(" o bigI_tr(<%d x %d>) ..\n", nr,nc);
    return( bigintegerR::create_SEXP(matrixz::bigint_transpose(mat)));
  }
  catch(std::invalid_argument & e){
    Rf_error("%s",e.what());
    return Rf_mkString(0);
 }
}


/* \brief  matrix cross product
 *
 * \param a matrix (n x p)
 * \param trans if(trans), compute tcrossprod(), else crossprod()
 * \return  crossprod(a) := t(a) %*% a  [p x p]  or
 *         tcrossprod(a) := a %*% t(a)  [n x n]
 */
SEXP matrix_crossp_z (SEXP a, SEXP trans)
{
  try{
    bool useMod = FALSE,
      tr = (bool)Rf_asLogical(trans);
    bigvec mat_a = bigintegerR::create_bignum(a);
    int sizemod = mat_a.getModulusSize(),
      a_nrow = mat_a.nrow,
      a_len = mat_a.size();

    // in case of a vector; crossprod() returns scalar product,
    // whereas             tcrossprod() gives  n x n matrix.
    if(a_nrow < 0)
      a_nrow = a_len;
    int a_ncol = a_len / a_nrow;

    // Result R is R[1..m, 1..m] -- and R_{ij} = sum_{k=1}^p  A.. B..
    int m, p;
    if(tr) { // tcrossprod()
      m= a_nrow; p= a_ncol;
    } else { //  crossprod()
      m= a_ncol; p= a_nrow;
    }
    bigvec res(m*m);
    res.nrow= m;

    mpz_t R_ij, tt;
    mpz_init(R_ij);
    mpz_init(tt);
    mpz_t common_modulus; mpz_init(common_modulus);

    if(sizemod == 1) {
      mpz_set(common_modulus, mat_a.getGlobalModulus()->getValueTemp());
      useMod = TRUE;
    }
 

    // here the computation:
    for(int i=0; i < m; i++)
      for(int j=0; j < m; j++) {
	mpz_set_ui(R_ij, 0);
	bool isna = false;
#define K_LOOP								\
	for(int k=0; k < p; k++) {					\
	  /* R_ij = \sum_{k=1}^p  a_{ik} b_{kj} */			\
	  if( !(A_I_K.isNA() || B_K_J.isNA())) {			\
	    mpz_mul(tt, A_I_K.getValueTemp(), B_K_J.getValueTemp());	\
	    mpz_add(R_ij, tt,R_ij);					\
	  }								\
	  else {							\
	    isna = true; break;						\
	  }								\
	}

	if(tr) {//------------- tcrossprod ---------------------------

#define A_I_K  mat_a [ i + k *a_nrow]
#define B_K_J  mat_a [ j + k *a_nrow]
	  K_LOOP
#undef A_I_K
#undef B_K_J

	    } else {//------------- crossprod ---------------------------

#define A_I_K  mat_a [ k + i *a_nrow]
#define B_K_J  mat_a [ k + j *a_nrow]
	  K_LOOP
#undef A_I_K
#undef B_K_J

	    }

	if(isna) {
	  res[i + j*m].setValue(0);
	  res[i + j*m].getValue().NA(true);
	}
	else
	  res[i + j*m].setValue(R_ij);

      } // for(i ..)  for(j ..)

    if(useMod) {
      std::shared_ptr<biginteger> ptr = std::make_shared<biginteger>(common_modulus);
      res.setGlobalModulus(ptr);
    }
    mpz_clear(R_ij);
    mpz_clear(tt);
    mpz_clear(common_modulus);

    return( bigintegerR::create_SEXP (res));
  }
  catch(std::invalid_argument & e){
    Rf_error("%s",e.what());
    return Rf_mkString(0);
  }
  
} // matrix_crossp_z()
#undef K_LOOP

/* \brief  matrix multiplication
 *
 * returns matrix multiplication  T(a) %*% b  or  b %*% T(a)
 * \param a matrix
 * \param b matrix
 * \param op operation code: 0: %*%,  1: crossprod,  2: tcrossprod
 *     (same codes as in R's do_matprod() in src/main/array.c )
 */
SEXP matrix_mul_z (SEXP a, SEXP b, SEXP op)
{
  try{
    if(!strcmp(class_P(b), "bigq")) { // b  "bigq",  --> use q arithm:
      return(matrix_mul_q(a, b, op));
    }
    // FIXME: we may know that a is 'bigz' - but we don't know at all about b !!
    // -----  create_bignum(.) should be much more careful (better: have a careful option!)
    bool useMod = FALSE;// if(useMod)  use a *common* modulus
    int o_ = Rf_asInteger(op); // INTEGER(op)[0]
    bigvec mat_a = bigintegerR::create_bignum(a),
      mat_b = bigintegerR::create_bignum(b);

  


    int
      a_nrow = mat_a.nrow, a_len = mat_a.size(),
      b_nrow = mat_b.nrow, b_len = mat_b.size(),
      a_ncol = -1, b_ncol = -1;// -Wall

    // distinguish cases of vectors / matrices ---------------------
    if(a_nrow < 0) {
      if(b_nrow < 0) { // *both* are vectors
	if(o_ == 0) {
	  a_nrow = 1;
	  a_ncol = a_len;
	} else {
	  a_nrow = a_len;
	  a_ncol = 1;
	}
	b_nrow = b_len;
	b_ncol = 1;

      } else { // a : vector,   b : matrix
	b_ncol = b_len / b_nrow;
	if(o_ == 0) {
	  if (a_len == b_nrow) {	/* x as row vector */
	    a_nrow = 1;
	    a_ncol = b_nrow; /* == a_len */
	  }
	  else if (b_nrow == 1) {	/* x as col vector */
	    a_nrow = a_len;
	    a_ncol = 1;
	  }
	} else if(o_ == 1) { /* crossprod() */
	  if (a_len == b_nrow) {	/* x is a col vector */
	    a_nrow = b_nrow; /* == a_len */
	    a_ncol = 1;
	  }
	  /* else if (b_nrow == 1) ... not being too tolerant
	     to treat x as row vector, as t(x) *is* row vector */
	} else { // o_ == 2 -- tcrossprod()
	  if (a_len == b_ncol) {	/* x as row vector */
	    a_nrow = 1;
	    a_ncol = b_ncol; /* == a_len */
	  }
	  else if (b_ncol == 1) {	/* x as col vector */
	    a_nrow = a_len;
	    a_ncol = 1;
	  }
	}
      }
    }
    else if (b_nrow < 0) { // a : matrix,   b : vector
      a_ncol = a_len / a_nrow;
      if (o_ == 0) {
	if (b_len == a_ncol) {	/* y as col vector */
	  b_nrow = a_ncol;
	  b_ncol = 1;
	}
	else if (a_ncol == 1) {	/* y as row vector */
	  b_nrow = 1;
	  b_ncol = b_len;
	}
      }
      else if (o_ == 1) { /* crossprod() */
	if (b_len == a_nrow) {	/* y is a col vector */
	  b_nrow = a_nrow;
	  b_ncol = 1;
	}
      }
      else { /* tcrossprod --	   y is a col vector */
	b_nrow = b_len;
	b_ncol = 1;
      }

    } else { // a, b  *both* matrices
      a_ncol = a_len / a_nrow;
      b_ncol = b_len / b_nrow;
    }

    if(((o_ == 0) && (a_ncol != b_nrow)) ||
       ((o_ == 1) && (a_nrow != b_nrow)) || // crossprod()
       ((o_ == 2) && (a_ncol != b_ncol))    // tcrossprod()
       ){
      mat_a.clear();
      mat_b.clear();
      throw invalid_argument(_("Matrix dimensions do not match"));
    }
    // Result R is R[1..n, 1..m] -- and R_{ij} = sum_{k=1} ^ p  A.. B..
    int n,m, p;
    if(o_ == 0) {
      n= a_nrow; m= b_ncol; p= a_ncol;// = b_nrow
    }else if (o_ == 1) {
      n= a_ncol; m= b_ncol; p= a_nrow;// = b_nrow
    }else if (o_ == 2) {
      n= a_nrow; m= b_nrow; p= a_ncol;// = b_ncol
    }else {
      mat_a.clear();
      mat_b.clear();
      throw invalid_argument(_("invalid 'op' code in matrix_mul_z()"));
      n = m = p = -1;// -Wall
    }

    bigvec res(n*m);
    res.nrow=n;

    mpz_t common_modulus, tt;
    mpz_init(tt);
    mpz_init(common_modulus);

    /* modulus when modulus are "global" (i.e. of size 1) and
     * either are the same, or only one of a or b is specified
     */
    std::shared_ptr<biginteger> globalMod = bigvec::getGlobalModulus(mat_a,mat_b);
    useMod = globalMod.get() != nullptr;
    if(useMod){
      mpz_init_set(common_modulus,globalMod->getValueTemp());
    }
    // bigmod tmp;

    // here the computation:
    for(int i=0; i < n; i++)
      for(int j=0; j < m; j++)
	{
#define	R_IJ res[ i + j*n]
#define K_LOOP								\
	  for(int k=0; k < p; k++)					\
	    {								\
	      if(A_I_K.isNA() || B_K_J.isNA()) {			\
		R_IJ.setValue(0); R_IJ.getValue().NA(true);		\
		break;							\
	      }								\
	      /* Z = A_I_K * B_K_J */					\
	      mpz_mul(tt, A_I_K.getValueTemp(), B_K_J.getValueTemp());	\
	      /* R_IJ = R_IJ + A_I_K * B_K_J  */			\
	      mpz_add(tt, tt, R_IJ.getValueTemp());			\
	      if(useMod)						\
		mpz_mod(tt,tt,common_modulus);				\
	      R_IJ.setValue(tt);					\
	    }

	  R_IJ.setValue(0);

	  if(o_ == 0) { //------------- %*% ---------------------------

#define A_I_K  mat_a [ i + k *a_nrow]
#define B_K_J  mat_b [ k + j *b_nrow]
	    K_LOOP
#undef A_I_K
#undef B_K_J

	      } else if(o_ == 1){//------------- crossprod ---------------------------

#define A_I_K  mat_a [ k + i *a_nrow]
#define B_K_J  mat_b [ k + j *b_nrow]
	    K_LOOP
#undef A_I_K
#undef B_K_J

	      } else {//(o_ == 2) ------------- tcrossprod ---------------------------

#define A_I_K  mat_a [ i + k *a_nrow]
#define B_K_J  mat_b [ j + k *b_nrow]
	    K_LOOP
#undef A_I_K
#undef B_K_J

	      }
	}
    if(useMod){
      std::shared_ptr<biginteger> ptr = std::make_shared<biginteger>(common_modulus);
      res.setGlobalModulus(ptr);
    }
    mpz_clear(tt);
    mpz_clear(common_modulus);

    return( bigintegerR::create_SEXP (res));
  }
  catch(std::invalid_argument & e){
    Rf_error("%s",e.what());
    return Rf_mkString(0);
 }
} // matrix_mul_z()
#undef R_IJ
#undef K_LOOP



/**
 * arg = v1, v2, v3.... lines
 * return matrix
 * out =[ v1 
 v2
 v3 ]
*/
SEXP biginteger_rbind(SEXP args)
{
  try{
    bigvec result;
    vector<bigvec*> source;
    unsigned int maxSize=0;
    for(int i = 0 ; i < LENGTH(args) ; i++){
      bigvec v = bigintegerR::create_bignum(VECTOR_ELT(args,i));

      if(v.size() == 0) continue;
      for (unsigned int row = 0 ; row < v.nRows(); row++){
	bigvec * line = new bigvec();
	for(unsigned int col = 0 ; col < v.nCols(); col++){
	  line->push_back(v.get(row,col));
	}
	source.push_back(line);	
	maxSize = std::max(maxSize,line->size());
      }
    }


  
    for (unsigned int j = 0 ; j < maxSize; j++){
      for(unsigned int i = 0 ; i < source.size() ; i++){
	bigvec *  u = source[i];
	if(u->size() == 0) result.push_back(bigmod());
	else    result.push_back((*u)[j % u->size()]); 
      }
    }
    result.nrow =  source.size();
    for (unsigned int i = 0 ; i < source.size() ; i++){
      if(source[i] != nullptr) delete source[i];
      source[i] = nullptr;
    }
    return bigintegerR::create_SEXP(result);
  }
  catch(std::invalid_argument & e){
    Rf_error("%s",e.what());
    return Rf_mkString(0);
 }
}

/**
 * arg = v1, v2, v3.... row
 * return matrix
 * out =[ v1 v2 v3 ]
 */
SEXP biginteger_cbind(SEXP args)
{
  try{
    bigvec result;
    vector<bigvec*> source;
    unsigned int maxSize=0;
    for(int i = 0 ; i < LENGTH(args) ; i++){
      bigvec v = bigintegerR::create_bignum(VECTOR_ELT(args,i));
 
      if(v.size() == 0) continue;
      if(v.nrow <0) v.nrow = v.size();
      for(unsigned int col = 0 ; col < v.nCols(); col++){
	bigvec * column = new bigvec();
	for (unsigned int row = 0 ; row < v.nRows(); row++){
	  column->push_back(v.get(row,col));
	}
	source.push_back(column);	
	maxSize = std::max(maxSize,column->size());
      }
    }

    for(unsigned int i = 0 ; i < source.size() ; i++){
      bigvec *  u = source[i];
      for (unsigned int j = 0 ; j < maxSize; j++){
	if(u->size() == 0) result.push_back(bigmod());
	else    result.push_back((*u)[j % u->size()]); 
      }
    }
    result.nrow = result.size() /  source.size();
    for (unsigned int i = 0 ; i < source.size() ; i++){
      if(source[i] != nullptr) delete source[i];
      source[i] = nullptr;
    }
  
    return bigintegerR::create_SEXP(result);
  } catch(std::invalid_argument & e){
    Rf_error("%s",e.what());
    return Rf_mkString(0);
  }
}



namespace matrixz
{
  bigvec bigint_transpose ( bigvec & mat)
  {
    bigvec matbis (mat.size());
    matbis.nrow = mat.nCols();

    if(mat.getType() == MODULUS_GLOBAL){
      matbis.setGlobalModulus(mat.getGlobalModulus());
    }
    
    /* we compute transpose */
    for(unsigned int i =0; i<mat.nRows(); i++)
      for(unsigned int j =0; j<mat.nCols(); j++)
	matbis.set(j+i*mat.nCols(),mat[i+j*mat.nRows()]);

    return matbis;
  }

  /* return dimension in dima return -2 if dimension does not match */
  int checkDims(int dima, int dimb)
  {
    if(dima > 0 && dimb > 0) {
      if (dimb != dima) return -2;
      //error(_("Matrix dimensions do not match"));
    }
    else { /* either a or b is a matrix */
      if(dima == -1)
	return(dimb);
    }
    return(dima);
  }
}
