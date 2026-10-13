/*

  $Id: permdist.c,v 1.27.4.1 2003/08/30 14:16:53 hothorn Exp $
  
  permdist : Distribution of Permutation Tests by Streitberg & Roehmel
  Copyright (C) 2000-2003  Torsten Hothorn 
                           <Torsten.Hothorn@rzmail.uni-erlangen.de>
    
  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or (at
  your option) any later version.
            
  This program is distributed in the hope that it will be useful, but
  WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  General Public License for more details.
                   
  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA.
                          
  SYNOPSIS

    SEXP cpermdist1(SEXP scores)                    
    SEXP cpermdist2(SEXP m_a,  SEXP m_b, 
                    SEXP score_a, SEXP score_b, SEXP retProb)
                                             
  DESCRIPTION
                             
    cpermdist1	The density of the permutation distribution for 
                the symmetry problem.
    cpermdist2	The density of the permutation distribution for 
  		the independent two sample problem.
                                                                  
  REFERENCES

    Bernd Streitberg & Joachim R\"ohmel (1986),
    Exact distributions for permutations and rank tests:
    An introduction to some recently published algorithms. 
    Statistical Software Newsletter 12(1), 10-17.
                         
    Bernd Streitberg & Joachim R\"ohmel (1987),
    Exakte Verteilungen f\"ur Rang- und Randomisierungstests
    im allgemeinen $c$-Stichprobenfall.
    EDV in Medizin und Biologie 18(1), 12-19 (in german).

*/

                                                                   
#include <R.h>
#include <Rmath.h>
#include <Rdefines.h>

/*
	length(scores) <= 1.000.000 observations only.
*/

#define PERM_MAX_N 1000000


int aindx(int i, int j, int n) {
  /* 
    array indexing for vectors: for a (m x n) matrix get position of element
    (j,i)
  */
  return(i*(n + 1) + j);
}



SEXP cpermdist1(SEXP scores) {

  /*
    compute the permutation distribution of the 
    sum of the absolute values of the positive elements of `scores'
  */ 

  int N; 	/* number of observations */ 
  SEXP H;	/* vector giving the density of statistics 0:sum(scores) */
  
  int i, k, sum_a = 0, s_a = 0; /* little helpers */
  double msum = 0.0;
	
  if (!isVector(scores)) {
    error("scores is not a vector");
  }
	         
  N = LENGTH(scores);
	               
  if (N > PERM_MAX_N)
    error("N > %d in cpermdistr1", PERM_MAX_N); 
	
  for (i = 0; i < N; i++) {
    if (INTEGER(scores)[i] < 0)
      error("score for observation number %d is negative", i);
    sum_a += INTEGER(scores)[i];
  }

  /*
    Initialize H
  */

  PROTECT(H = allocVector(REALSXP, sum_a + 1));
  for (i = 0; i <= sum_a; i++) REAL(H)[i] = 0.0;

  /*
    start the Shift-Algorithm with H[0] = 1.0
  */
		
  REAL(H)[0] = 1.0;
	
  for (k = 0; k < N; k++) {
    s_a = s_a + INTEGER(scores)[k];
    for (i = s_a; i >= INTEGER(scores)[k]; i--)
      REAL(H)[i] = REAL(H)[i] + REAL(H)[i - INTEGER(scores)[k]];
  }


  /* 
    get the number of permutations
  */

  for (i = 0; i <= sum_a; i++)
    msum += REAL(H)[i];
	
  /*
    compute probabilities and return the density H to R
    [dpq] stuff is done in R
  */ 
	
  for (i = 0; i <= sum_a; i++)
    REAL(H)[i] = REAL(H)[i]/msum;	/* 0 is a possible realization */

  UNPROTECT(1);	
  return(H);
}

SEXP cpermdist2(SEXP m_a,  SEXP m_b, 
                SEXP score_a, SEXP score_b, SEXP retProb) {
  /*
    compute the joint permutation distribution of the 
    sum of the first m_a elements of score_a and score_b
    (usualy score_a = rep(1, length(score_a)) and 
            score_b = Data scores, Wilcoxon, Ansari ...).
    In this case the exact conditional distribution 
    in the simple independent two-sample problem is computed.
    The logical retProb indicates wheater the density for sample 
    size m_a should be returned or the whole matrix of all 
    possible permutations for sample sizes 0:m_a.
  */ 

  int N, m, c;		/* number of observations */

  SEXP H, x;		/* matrix of permutations and vector 
                           of probabilities */ 
  
  int i, j, k, sum_a = 0, sum_b = 0, s_a = 0, s_b = 0;
  double msum = 0.0; 	/* little helpers */

  if (!isVector(score_a)) {
    error("score_a is not a vector");
  }

  N = LENGTH(score_a);

  if (!isVector(score_b)) {
    error("score_b is not a vector");
  }
        
  if (LENGTH(score_b) != N) {
    error("length of score_a and score_b differ");
  }
        
  if (TYPEOF(retProb) != LGLSXP)
    error("retProb is not a logical");                                  

  m = INTEGER(m_a)[0];  /* cosmetics only */
  c = INTEGER(m_b)[0];
        

  if (N > PERM_MAX_N)
    error("N > %d in cpermdistr2", PERM_MAX_N); 

  /* compute the total sum of the scores and check if they are >= 0 */
	
  for (i = 0; i < N; i++) {
    if (INTEGER(score_a)[i] < 0) 
      error("score_a for observation number %d is negative", i);
    if (INTEGER(score_b)[i] < 0) 
      error("score_b for observation number %d is negative", i);
    sum_a += INTEGER(score_a)[i];
    sum_b += INTEGER(score_b)[i];
  }

  /*
    optimization according to Streitberg & Roehmel
  */
	
  sum_a = imin2(sum_a, m);
  sum_b = imin2(sum_b, c);

  /*
    initialize H
    note: there is no advantage of using a matrix since we need to use vector
          indexing anyway (see aindx).
  */

  PROTECT(H = allocVector(REALSXP, (sum_a + 1) * (sum_b + 1)));

  for (i = 0; i <= sum_a; i++) {
    for (j = 0; j <= sum_b; j++)
      REAL(H)[aindx(i,j,sum_b)] = 0.0;
  }
		
  /*
    start the Shift-Algorithm with H[0][0] = 1
  */
		
  REAL(H)[aindx(0,0,sum_b)] = 1.0;
	
  for (k = 0; k < N; k++) {
    s_a = s_a + INTEGER(score_a)[k];
    s_b = s_b + INTEGER(score_b)[k];

    /*
      compute H up to row m and column c
      note: 
        sum_a = min(sum_a, m)
        sum_b = min(sum_b, c)
    */
		
    for (i = imin2(m, s_a); i >= INTEGER(score_a)[k]; i--) {
      for (j = imin2(c,s_b); j >= INTEGER(score_b)[k]; j--) {
        REAL(H)[aindx(i,j,sum_b)] = REAL(H)[aindx(i,j,sum_b)] + 
          REAL(H)[aindx(i - INTEGER(score_a)[k],
                        j - INTEGER(score_b)[k], sum_b)];
      }
    }
  }

  /*
    return the whole matrix H 
  */ 

  if (!LOGICAL(retProb)[0]) {
    UNPROTECT(1);
    return(H);
    /* note: use matrix(H, nrow=m_a+1, byrow=TRUE) in R */
  } else {	
    PROTECT(x = allocVector(REALSXP, sum_b));

    /* 
      get the values for sample size m_a (in row m) and sum it up
    */

    for (j = 0; j < sum_b; j++) {
      REAL(x)[j] = REAL(H)[aindx(m, j+1, sum_b)];
      msum += REAL(x)[j];
    }
	
    /*
      compute probabilities and return the density H to R
      [dpq] stuff is done in R
    */
            	
    for (j = 0; j < sum_b; j++)
      REAL(x)[j] = REAL(x)[j]/msum;
		
    UNPROTECT(2);
    return(x);
  }
}
