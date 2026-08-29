# include <RcppArmadillo.h>
// [[ Rcpp :: depends ( RcppArmadillo )]]


//' Computes the initial distribution and sub-intensity of the sum of two 
//'  phase-type distributed random variables. 
//' 
//' @param alpha1 Initial distribution.
//' @param S1 Sub-intensity matrix.
//' @param alpha2 Initial distribution.
//' @param S2 Sub-intensity matrix.
//' 
// [[Rcpp::export]]
Rcpp::List sum_ph(arma::rowvec alpha1, arma::mat S1, arma::rowvec alpha2, arma::mat S2) {
  unsigned p1{S1.n_cols};
  unsigned p2{S2.n_cols};
  
  arma::rowvec alpha_sum(p1 + p2);
  arma::mat S_sum(p1 + p2, p1 + p2);
  
  arma::mat e;
  e.ones(S1.n_cols, 1);
  arma::mat exit_vect = (S1 * (-1)) * e;
  
  arma::mat aux_mat = exit_vect * alpha2;
  
  for (int i{0}; i < p1 + p2; ++i) {
    if (i < p1) {
      alpha_sum[i] = alpha1[i];
    }
    
    for (int j{0}; j < p1 + p2; ++j) {
      if (i < p1) {
        if (j < p1) {
          S_sum(i,j) = S1(i,j);
        }
        else {
          S_sum(i,j) = aux_mat(i,j - p1);
        }
      }
      else if (i >= p1 && j>= p1) {
        S_sum(i,j) = S2(i - p1,j - p1);
      }
    }
  }
  
  Rcpp::List L = Rcpp::List::create(Rcpp::Named("alpha") = alpha_sum, Rcpp::Named("S") = S_sum);
  return L;
}


//' Computes the initial distribution and sub-intensity of the sum of two 
//'  discrete phase-type distributed random variables
//' 
//' @param alpha1 Initial distribution.
//' @param S1 Sub-transition matrix.
//' @param alpha2 Initial distribution.
//' @param S2 Sub-transition matrix.
//' 
// [[Rcpp::export]]
Rcpp::List sum_dph(arma::rowvec alpha1, arma::mat S1, arma::rowvec alpha2, arma::mat S2) {
  unsigned p1{S1.n_cols};
  unsigned p2{S2.n_cols};
  
  arma::rowvec alpha_sum(p1 + p2);
  arma::mat S_sum(p1 + p2, p1 + p2);
  
  arma::mat e;
  e.ones(S1.n_cols, 1);
  arma::mat exit_vect = e - (S1 * e);
  
  arma::mat aux_mat = exit_vect * alpha2;
  
  for (int i{0}; i < p1 + p2; ++i) {
    if (i < p1) {
      alpha_sum[i] = alpha1[i];
    }
    
    for (int j{0}; j < p1 + p2; ++j) {
      if (i < p1) {
        if (j < p1) {
          S_sum(i,j) = S1(i,j);
        }
        else {
          S_sum(i,j) = aux_mat(i,j - p1);
        }
      }
      else if (i >= p1 && j>= p1) {
        S_sum(i,j) = S2(i - p1,j - p1);
      }
    }
  }
  
  Rcpp::List L = Rcpp::List::create(Rcpp::Named("alpha") = alpha_sum, Rcpp::Named("S") = S_sum);
  return L;
}

