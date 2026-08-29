#include <RcppArmadillo.h>
#include "transformations.h"
#include "Simulation.h"

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins("cpp11")]]

//' Find how many states have positive reward
//'
//' @param R reward vector
//'
//' @return The number of states with positive rewards
//'
// [[Rcpp::export]]
int n_pos(arma::vec R) {
  int d{0};
  
  for (int i{0}; i < R.size(); ++i) {
    if (R(i) > 0){
      ++d;
    }
  }
  return d;
}


//' Find which states have positive reward
//'
//' @param R reward vector
//'
//' @return A vector with the states (number) that are associated with positive rewards
//'
// [[Rcpp::export]]
arma::vec plus_states(arma::vec R) {
  int d = n_pos(R);
  int j{0};
  
  arma::vec plusState(d);
  
  for (int i{0}; i < R.size(); ++i) {
    if (R(i) > 0) {
      plusState(j) = i + 1;
      ++j;
    }
  }
  return plusState;
}


//' Performs TVR for phase-type distributions
//'
//' @param alpha Initial distribution vector.
//' @param S Sub-intensity matrix.
//' @param R Reward vector.
//'
//' @return A list of phase-type parameters.
//'
// [[Rcpp::export]]
Rcpp::List tvr_ph(arma::vec alpha, arma::mat S, arma::vec R) {
  unsigned p{S.n_rows};
  
  unsigned n0{0};
  std::vector<int> delete_rows; 
  std::vector<int> keep_rows;
  
  for (int j{0}; j < p; ++j) {
    if (R[j] == 0) {
      delete_rows.push_back(j);
      ++n0;
    }
    else {
      keep_rows.push_back(j);
    }
  }
  
  arma::mat Q_aux = embedded_mc(S);
  
  unsigned np{p - n0};
  
  arma::mat Qpp(np,np);
  arma::mat Qp0(np,n0);
  arma::mat Q0p(n0,np);
  arma::mat Q00(n0,n0);
  
  arma::rowvec alpha0(n0);
  arma::rowvec alphap(np);
  
  for (int i{0}; i < np; i++) {
    for (int j = 0; j < np; j++) {
      Qpp(i,j) = Q_aux(keep_rows[i],keep_rows[j]);
    }
    for (int j{0}; j < n0; j++) {
      Qp0(i,j) = Q_aux(keep_rows[i],delete_rows[j]);
    }
    alphap[i] = alpha[keep_rows[i]];
  }
  for (int i{0}; i < n0; i++) {
    for (int j{0}; j < np; j++) {
      Q0p(i,j) = Q_aux(delete_rows[i],keep_rows[j]);
    }
    for (int j{0}; j < n0; j++){
      Q00(i,j) = Q_aux(delete_rows[i],delete_rows[j]);
    }
    alpha0[i] = alpha[delete_rows[i]];
  }
  
  arma::rowvec alpha_trans = alphap + alpha0 * inv(arma::eye(n0,n0) - Q00) * Q0p;
  arma::mat P_trans = Qpp + Qp0 * inv(arma::eye(n0,n0) - Q00) * Q0p;
  
  arma::mat e;
  e.ones(np, 1);
  arma::mat exit_vect = e - (P_trans * e) ;
  
  arma::vec row_sum(np);
  arma::mat S_trans(np,np);
  
  for (int i{0}; i < np; i++) {
    for (int j = 0; j < np; j++) {
      if (i != j) {
        S_trans(i,j) = - S(keep_rows[i],keep_rows[i]) * P_trans(i,j) / R[keep_rows[i]];
        row_sum[i] += S_trans(i,j);
      }
    }
    S_trans(i,i) = - row_sum[i] + S(keep_rows[i],keep_rows[i]) * exit_vect[i] / R[keep_rows[i]];
  }
  
  Rcpp::List x = Rcpp::List::create(Rcpp::Named("alpha") = alpha_trans, Rcpp::_["S"] = S_trans);
  return (x);
}


//' Performs TVR for discrete phase-type distributions
//'
//' @param alpha Initial distribution vector.
//' @param S Sub-intensity matrix.
//' @param R Reward vector.
//'
//' @return A list of PH parameters.
//' 
// [[Rcpp::export]]
Rcpp::List tvr_dph(arma::vec alpha, arma::mat S, arma::vec R) {
  unsigned p{S.n_rows};
  
  unsigned n0{0};
  std::vector<int> delete_rows; 
  std::vector<int> keep_rows;
  
  for (int j{0}; j < p; ++j) {
    if (R[j] == 0) {
      delete_rows.push_back(j);
      ++n0;
    }
    else {
      keep_rows.push_back(j);
    }
  }
  
  unsigned np{p - n0};
  
  arma::rowvec alpha_trans(np);
  arma::mat S_trans(np,np);
  
  arma::mat Spp(np,np);
  arma::mat Sp0(np,n0);
  arma::mat S0p(n0,np);
  arma::mat S00(n0,n0);
    
  arma::rowvec alpha0(n0);
  arma::rowvec alphap(np);
    
  for (int i{0}; i < np; i++) {
    for (int j = 0; j < np; j++) {
      Spp(i,j) = S(keep_rows[i],keep_rows[j]);
    }
    for (int j{0}; j < n0; j++) {
      Sp0(i,j) = S(keep_rows[i],delete_rows[j]);
    }
    alphap[i] = alpha[keep_rows[i]];
  }
  for (int i{0}; i < n0; i++) {
    for (int j{0}; j < np; j++) {
      S0p(i,j) = S(delete_rows[i],keep_rows[j]);
    }
    for (int j{0}; j < n0; j++){
      S00(i,j) = S(delete_rows[i],delete_rows[j]);
    }
    alpha0[i] = alpha[delete_rows[i]];
  }
  
  alpha_trans = alphap + alpha0 * inv(arma::eye(n0,n0) - S00) * S0p;
  S_trans = Spp + Sp0 * inv(arma::eye(n0,n0) - S00) * S0p;
  
  Rcpp::List x = Rcpp::List::create(Rcpp::Named("alpha") = alpha_trans, Rcpp::_["S"] = S_trans);
  return (x);
}


//' Computes PH parameters of a linear combination of vector from MPHstar
//'
//' @param w Vector with weights.
//' @param alpha Initial distribution vector.
//' @param S Sub-intensity matrix.
//' @param R Reward matrix.
//'
//' @return A list of PH parameters.
//' 
// [[Rcpp::export]]
Rcpp::List linear_combination(arma::vec w, arma::vec alpha, arma::mat S, arma::mat R) {
  unsigned p{S.n_rows};
  
  unsigned n0{0};
  std::vector<int> delete_rows; 
  std::vector<int> keep_rows;
  
  arma::mat Rw_m = R * w;
  arma::vec Rw = Rw_m.col(0);
  
  for (int j{0}; j < p; ++j) {
    if (Rw[j] == 0) {
      delete_rows.push_back(j);
      ++n0;
    }
    else {
      keep_rows.push_back(j);
    }
  }
  
  unsigned np{p - n0};
  
  arma::rowvec alpha_trans(np);
  arma::mat S_trans(np,np);
  
  if (n0 == 0) {
    S_trans = diagmat(Rw) * S;
    alpha_trans = alpha.t();
  } else {
    arma::mat Spp(np,np);
    arma::mat Sp0(np,n0);
    arma::mat S0p(n0,np);
    arma::mat S00(n0,n0);
    
    arma::rowvec alpha0(n0);
    arma::rowvec alphap(np);
    
    for (int i{0}; i < np; i++) {
      for (int j = 0; j < np; j++) {
        Spp(i,j) = S(keep_rows[i],keep_rows[j]);
      }
      for (int j{0}; j < n0; j++) {
        Sp0(i,j) = S(keep_rows[i],delete_rows[j]);
      }
      alphap[i] = alpha[keep_rows[i]];
    }
    for (int i{0}; i < n0; i++) {
      for (int j{0}; j < np; j++) {
        S0p(i,j) = S(delete_rows[i],keep_rows[j]);
      }
      for (int j{0}; j < n0; j++){
        S00(i,j) = S(delete_rows[i],delete_rows[j]);
      }
      alpha0[i] = alpha[delete_rows[i]];
    }
    
    alpha_trans = alphap + alpha0 * inv(S00 * (-1)) * S0p;
    arma::mat S_aux = Spp + Sp0 * inv(S00 * (-1)) * S0p;
    
    arma::mat diagonal(np, np);
    for (int i{0}; i < np; ++i) {
      diagonal(i,i) = 1.0 / Rw[keep_rows[i]];
    }
    S_trans = diagonal * S_aux;
  }
  
  Rcpp::List x = Rcpp::List::create(Rcpp::Named("alpha") = alpha_trans, Rcpp::_["S"] = S_trans);
  return (x);
}
