/*

  $Id: spermdist.c,v 1.5.4.1 2003/08/30 14:16:53 hothorn Exp $
  
  spermdist : Simulate Distribution of Permutation Test Statistics
  Copyright (C) 2003  Torsten Hothorn 
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
               
    SEXP sim2is(SEXP scores, SEXP mfirst, SEXP Nsim)
                                             
  DESCRIPTION
                                                  
    sim2is	Simulate permutation distribution for the symmetry and 
                independent two-sample problems.
                                                                  
*/

                                                                   
#include <R.h>
#include <Rmath.h>
#include <Rdefines.h>

SEXP sim2is(SEXP scores, SEXP mfirst, SEXP Nsim) {

  int N;		/* sample size: length(scores) */
  int m;		/* number of observations in the first group:
  			   if m == N, paired observations are assumed: mfirst 
  			*/
  int ns;		/* number of simulation runs: Nsim */
  
  SEXP distr;		/* vector of length ns for storing the 
  			   simulated statistics */
  SEXP counts;		/* vector of length ns for counting duplicates */
  SEXP dlist;		/* a list of two elements: */
  SEXP T, prob;		/* the vector of possible realizations T and
  			   the vector of corresponding probabilities */
   

  /*
  	little helpers
  */
  
  int i, j, lasti, countus = 0, this, k;
  double stat = 0.0, cut;
  double *urand, *help;
              
  /*
  	check the arguments 
  */
              
  if (!isVector(scores)) {
    error("scores is not a vector");
  }

  /* for readability reasons only */

  m = INTEGER(mfirst)[0];
  N = LENGTH(scores);
  ns = INTEGER(Nsim)[0];
  
  /* needed for sampling without replacement */
  
  urand = (double *) R_alloc(N, sizeof(double));
  help = (double *) R_alloc(N, sizeof(double));
         
  PROTECT(distr = allocVector(REALSXP, ns));
  PROTECT(counts = allocVector(INTSXP, ns));

  /* initialize R's random number generator */

  GetRNGstate();

  /* ok, lets walk */
  
  for (i=0; i < ns; i++) {

    INTEGER(counts)[i] = 0;

    /* get N uniform random numbers */

    for (k = 0; k < N; k++) {
      urand[k] = unif_rand();
      help[k] = urand[k];
    }

    /* for the independent sample case only */

    if (m < N) {

      /* sort them and get the m smallest element */

      R_rsort(urand, N);
      cut = urand[m];
        
    } else {
      cut = 0.5;
    }                   
      
    stat = 0.0;
    
    /* 
      the statistic is the sum over m randomly choosen
      observations 
    */
    
    for (j = 0; j < N; j++)
      if (help[j] < cut) 
        stat += REAL(scores)[j];
    REAL(distr)[i] = stat;
  }
  
  PutRNGstate();

  /* we need to compute the distribution based on the simulation */
  
  /* 
     <FIXME>
     all of this is the overkill for simulating p-values only
     but currently I'll have the whole distribution at hand 
     (for conf.int's for example)
     </FIXME>
  */
    
  R_rsort(REAL(distr), ns);

  /* determine duplicates and its frequency */
  
  stat = REAL(distr)[0];
  lasti = 0;
  for (i = 0; i < ns; i++) {
    if (stat == REAL(distr)[i]) {
      INTEGER(counts)[lasti]++;
    } else {
      INTEGER(counts)[i]++;
      lasti = i;
    }
    if (INTEGER(counts)[i] == 0) countus++;
    stat = REAL(distr)[i];
  }
  
  /* 
     allocate memory for the return values: we do have ns-countus 
     unique statistics
  */
  
  countus = ns - countus;
  PROTECT(dlist = allocVector(VECSXP, 2));
  PROTECT(T = allocVector(REALSXP, countus));  
  PROTECT(prob = allocVector(REALSXP, countus));  

  /* we are ready to compute the (estimated) density now */
  
  this = 0;
  for (i = 0; i < ns; i++) {
    if (INTEGER(counts)[i] != 0) {
      REAL(T)[this] = REAL(distr)[i];
      REAL(prob)[this] = (double) INTEGER(counts)[i]/ns;
      this++;
    }
  }
  
  /* store the distribution in a lists with two elements */
  
  SET_VECTOR_ELT(dlist, 0, T);
  SET_VECTOR_ELT(dlist, 1, prob);
  
  UNPROTECT(5);
  return(dlist);
}
