#include "m_exp.h"
# include <RcppArmadillo.h>
// [[ Rcpp :: depends ( RcppArmadillo )]]

////////////////////////////////////////////
//    EM using uniformization            ///
////////////////////////////////////////////

//' Computes the elements S^n / n! until the a given size
//' 
//' @param vect A vector.
//' @param S Sub-intensity matrix.
//' @param a A number.
//' @param vect_size Size of vector.
//' 
// [[Rcpp::export]]
void vector_of_matrices(std::vector<arma::mat> & vect, const arma::mat & S, double a, int vect_size) {
  arma::mat I;
  I.eye(size(S));
  
  arma::mat P = I + S * (1 / a);
  
  vect.push_back(I);
  
  for (int k{1}; k <= vect_size; ++k) {
    vect.push_back((P * (1.0 / k) ) * vect[k - 1]);
  }
}


//' Computes exp(Sx) via series representation
//' 
//' @param x A number.
//' @param n An integer.
//' @param pow_vector A vector.
//' @param a A number.
//' 
// [[Rcpp::export]]
arma::mat m_exp_sum(double x, int n, const std::vector<arma::mat> & pow_vector, double a) {
  arma::mat res_mat = pow_vector[0];
  
  for (int i{1}; i <= n; ++i) {
    res_mat = res_mat + pow_vector[i] * exp(i * std::log(a * x));
  }
  res_mat = res_mat * exp(-a * x);
  
  return res_mat;
}


//' Computes A^(2^n)
//' 
//' @param n An integer.
//' @param A A matrix.
//' @return A^(2^n).
//' 
// [[Rcpp::export]]
void pow2_matrix(int n , arma::mat & A) {
  arma::mat aux_mat(size(A));
  
  for (int i{1}; i <= n; ++i) {
    aux_mat = A * A;
    A = aux_mat;
  }
}


//' Find n such that P(N > n) = h with N Poisson distributed
//'
//' @param h Probability.
//' @param lambda Mean of Poisson random variable.
//' @return Integer satisfying condition.
//'
// [[Rcpp::export]] 
int find_n(double h, double lambda) {
  int n{0};
  double cum_prob{0.0};
  
  do {
    cum_prob += R::dpois(n, lambda, false);
    ++n;
  } while (cum_prob < 1.0 - h);
  
  return (n - 1);
}


//' EM for phase-type using uniformization for matrix exponential
//' 
//' @param h Positive parameter.
//' @param alpha Initial probabilities.
//' @param S Sub-intensity matrix.
//' @param obs The observations.
//' @param weight The weights for the observations.
//' @param rcens Censored observations.
//' @param rcweight The weights for the censored observations.
//' 
// [[Rcpp::export]]
void EMstep_UNI(double h, arma::vec & alpha, arma::mat & S, const Rcpp::NumericVector & obs, const Rcpp::NumericVector & weight, const Rcpp::NumericVector & rcens, const Rcpp::NumericVector & rcweight) {
  unsigned p{S.n_rows};
  
  arma::mat e;
  e.ones(S.n_cols, 1);
  arma::mat exit_vect = (S * (-1)) * e;
  
  arma::mat Bmean = arma::zeros(p,1);
  arma::mat Zmean = arma::zeros(p,1);
  arma::mat Nmean = arma::zeros(p,p + 1);
  
  arma::mat avector(1,p);
  arma::mat bvector(p,1);
  arma::mat cmatrix(p,p);
  arma::mat aux_exp(p,p);
  
  arma::mat aux_mat(1,1);
  
  arma::mat J(2 * p,2 * p);
  arma::mat s_prod_alpha(p,p);
  s_prod_alpha = exit_vect * alpha.t();
  
  J = matrix_vanloan(S, S, s_prod_alpha);
  
  double a = max_diagonal(J * (-1));
  
  int m{find_n(h, 1)};
  
  std::vector<arma::mat> aux_vect;
  
  vector_of_matrices(aux_vect, J, a, m);
  
  double sum_weights{0.0};
  double density{0.0};
  
  // E-step
  // Uncensored data
  for (int k{0}; k < obs.size(); ++k) {
    sum_weights += weight[k];
    
    double x{obs[k]};
    
    if (x * a <= 1.0) {
      J = m_exp_sum(x, m, aux_vect, a);
    }
    else {
      int n{};
      n = std::log(a * x) / std::log(2.0);
      ++n;
      
      J = m_exp_sum(x / pow(2.0, n), m, aux_vect, a);
      
      pow2_matrix(n, J);
    }
    
    for (int i{0}; i < p; ++i) {
      for (int j{0}; j < p; ++j) {
        aux_exp(i,j) = J(i,j);
        cmatrix(i,j) = J(i,j + p);
      }
    }
    
    avector = alpha.t() * aux_exp;
    bvector = aux_exp * exit_vect;
    aux_mat = alpha.t() * bvector;
    density = aux_mat(0,0);
    
    // E-step
    for (int i{0}; i < p; ++i) {
      Bmean(i,0) += alpha[i] * bvector(i,0) * weight[k] / density;
      Nmean(i,p) += avector(0,i) * exit_vect(i,0) * weight[k] / density;
      Zmean(i,0) += cmatrix(i,i) * weight[k] / density;
      for (int j{0}; j < p; ++j) {
        Nmean(i,j) += S(i,j) * cmatrix(j,i) * weight[k] / density;
      }
    }
  }
  // Right-Censored Data
  double sum_censored{0.0};
  if (rcens.size() > 0) {
    s_prod_alpha = e * alpha.t();
    J = matrix_vanloan(S, S, s_prod_alpha);
    aux_vect.clear();
    vector_of_matrices(aux_vect, J, a, m);
  }
  for (int k{0}; k < rcens.size(); ++k) {
    sum_censored += rcweight[k];
    
    double x{rcens[k]};
    
    if (x * a <= 1.0) {
      J = m_exp_sum(x, m, aux_vect, a);
    }
    else {
      int n{};
      n = std::log(a * x) / std::log(2.0);
      ++n;
      
      J = m_exp_sum(x / pow(2.0, n), m, aux_vect, a);
      
      pow2_matrix(n, J);
    }
    
    for (int i{0}; i < p; ++i) {
      for (int j{0}; j < p; ++j) {
        aux_exp(i,j) = J(i,j);
        cmatrix(i,j) = J(i,j + p);
      }
    }
    
    bvector = aux_exp * e;
    aux_mat = alpha.t() * bvector;
    density = aux_mat(0,0);
    
    // E-step
    for (int i{0}; i < p; ++i) {
      Bmean(i,0) += alpha[i] * bvector(i,0) * rcweight[k] / density;
      Zmean(i,0) += cmatrix(i,i) * rcweight[k] / density;
      for (int j{0}; j < p; ++j) {
        Nmean(i,j) += S(i,j) * cmatrix(j,i) * rcweight[k] / density;
      }
    }
  }
  
  // M step
  for (int i{0}; i < p; ++i) {
    alpha[i] = Bmean(i,0) / (sum_weights + sum_censored);
    if (alpha[i] < 0) {
      alpha[i] = 0;
    }
    exit_vect(i,0) = Nmean(i,p) / Zmean(i,0);
    if (exit_vect(i,0) < 0) {
      exit_vect(i,0) = 0;
    }
    S(i,i) = -exit_vect(i,0);
    for (int j{0}; j < p; ++j) {
      if (i != j) {
        S(i,j) = Nmean(i,j) / Zmean(i,0);
        if (S(i,j) < 0) {
          S(i,j) = 0;
        }
        S(i,i) -= S(i,j);
      }
    }
  }
}


////////////////////////////////////////////
// Log-likelihoods
////////////////////////////////////////////

//' Loglikelihood of phase-type using uniformization
//' 
//' Loglikelihood for a sample.
//' 
//' @param h Positive parameter.
//' @param alpha Initial probabilities.
//' @param S Sub-intensity matrix.
//' @param obs The observations.
//' @param weight Weights of the observations.
//' @param rcens Censored observations.
//' @param rcweight Weights of the censored observations.
//' 
// [[Rcpp::export]]
double logLikelihoodPH_UNI(double h, arma::vec & alpha, arma::mat & S, const Rcpp::NumericVector & obs, const Rcpp::NumericVector & weight, const Rcpp::NumericVector & rcens, const Rcpp::NumericVector & rcweight) {
  arma::mat e;
  e.ones(S.n_cols, 1);
  arma::mat exit_vect = (S * (-1)) * e;
  
  arma::mat expm(size(S));
  
  double a = max_diagonal(S * (-1));
  
  int m{find_n(h, 1)};
  
  std::vector<arma::mat> aux_vect;
  
  vector_of_matrices(aux_vect, S, a, m);
  
  arma::mat aux_mat(1,1);
  
  double density{0.0};
  
  double logLh{0.0};
  
  // Non censored data
  for (int k{0}; k < obs.size(); ++k) {
    double x{obs[k]};
    
    if (x * a <= 1.0) {
      expm = m_exp_sum(x, m, aux_vect, a);
    }
    else {
      int n{};
      n = std::log(a * x) / std::log(2.0);
      ++n;
      
      expm = m_exp_sum(x / pow(2.0, n), m, aux_vect, a);
      
      pow2_matrix(n, expm);
    }
    aux_mat = alpha.t() * expm * exit_vect;
    density = aux_mat(0,0);
    logLh += weight[k] * std::log(density);
  }
  //Right censored data
  for (int k{0}; k < rcens.size(); ++k) {
    double x{rcens[k]};
    
    if (x * a <= 1.0) {
      expm = m_exp_sum(x, m, aux_vect, a);
    }
    else {
      int n{};
      n = std::log(a * x) / std::log(2.0);
      ++n;
      
      expm = m_exp_sum(x / pow(2.0, n), m, aux_vect, a);
      
      pow2_matrix(n, expm);
    }
    aux_mat = alpha.t() * expm * e;
    density = aux_mat(0,0);
    logLh += rcweight[k] * std::log(density);
  }
  
  return logLh;
}


//' Loglikelihood of matrix-Weibull using uniformization
//' 
//' Loglikelihood for a sample.
//' 
//' @param h Positive parameter.
//' @param alpha Initial probabilities.
//' @param S Sub-intensity matrix.
//' @param beta Parameter of transformation.
//' @param obs The observations.
//' @param weight Weights of the observations.
//' @param rcens censored observations.
//' @param rcweight Weights of the censored observations.
//' 
// [[Rcpp::export]]
double logLikelihoodMweibull_UNI(double h, arma::vec & alpha, arma::mat & S, double beta , const Rcpp::NumericVector & obs, const Rcpp::NumericVector & weight, const Rcpp::NumericVector & rcens, const Rcpp::NumericVector & rcweight) {
  if(beta < 0) return NA_REAL;
  
  arma::mat e;
  e.ones(S.n_cols, 1);
  arma::mat exit_vect = (S * (-1)) * e;
  
  arma::mat expm(size(S));
  
  double a = max_diagonal(S * (-1));
  
  int m{find_n(h, 1)};
  
  std::vector<arma::mat> aux_vect;
  
  vector_of_matrices(aux_vect, S, a, m);
  
  arma::mat aux_mat(1,1);
  
  double density{0.0};
  
  double logLh{0.0};
  
  // Non censored data
  for (int k{0}; k < obs.size(); ++k) {
    double x{pow(obs[k], beta)};
    
    if (x * a <= 1.0) {
      expm = m_exp_sum(x, m, aux_vect, a);
    }
    else {
      int n{};
      n = std::log(a * x) / std::log(2.0);
      ++n;
      
      expm = m_exp_sum(x / pow(2.0, n), m, aux_vect, a);
      
      pow2_matrix(n, expm);
    }
    aux_mat = alpha.t() * expm * exit_vect;
    density = aux_mat(0,0);
    logLh += weight[k] * (std::log(density) + std::log(beta) + (beta - 1) * std::log(obs[k]));
  }
  // Right censored data
  for (int k{0}; k < rcens.size(); ++k) {
    double x{pow(rcens[k], beta)};
    
    if (x * a <= 1.0) {
      expm = m_exp_sum(x, m, aux_vect, a);
    }
    else {
      int n{};
      n = std::log(a * x) / std::log(2.0);
      ++n;
      
      expm = m_exp_sum(x / pow(2.0, n), m, aux_vect, a);
      
      pow2_matrix(n, expm);
    }
    aux_mat = alpha.t() * expm * e;
    density = aux_mat(0,0);
    logLh += rcweight[k] * std::log(density);
  }
  
  return logLh;
}


//' Loglikelihood of matrix-Pareto using uniformization
//' 
//' Loglikelihood for a sample.
//' 
//' @param h Positive parameter.
//' @param alpha Initial probabilities.
//' @param S Sub-intensity matrix.
//' @param beta Parameter of transformation.
//' @param obs The observations.
//' @param weight Weights of the observations.
//' @param rcens censored observations.
//' @param rcweight Weights of the censored observations.
//' 
// [[Rcpp::export]]
double logLikelihoodMpareto_UNI(double h, arma::vec & alpha, arma::mat & S, double beta , const Rcpp::NumericVector & obs, const Rcpp::NumericVector & weight, const Rcpp::NumericVector & rcens, const Rcpp::NumericVector & rcweight) {
  if(beta < 0) return NA_REAL;
  
  arma::mat e;
  e.ones(S.n_cols, 1);
  arma::mat exit_vect = (S * (-1)) * e;
  
  arma::mat expm(size(S));
  
  double a = max_diagonal(S * (-1));
  
  int m{find_n(h, 1)};
  
  std::vector<arma::mat> aux_vect;
  
  vector_of_matrices(aux_vect, S, a, m);
  
  arma::mat aux_mat(1,1);
  
  double density{0.0};
  
  double logLh{0.0};
  
  // Non censored data
  for (int k{0}; k < obs.size(); ++k) {
    double x{std::log(obs[k] / beta + 1)};
    
    if (x * a <= 1.0) {
      expm = m_exp_sum(x, m, aux_vect, a);
    }
    else {
      int n{};
      n = std::log(a * x) / std::log(2.0);
      ++n;
      
      expm = m_exp_sum(x / pow(2.0, n), m, aux_vect, a);
      
      pow2_matrix(n, expm);
    }
    aux_mat = alpha.t() * expm * exit_vect;
    density = aux_mat(0,0);
    logLh += weight[k] * (std::log(density) - std::log(obs[k] + beta));
  }
  // Right censored data
  for (int k{0}; k < rcens.size(); ++k) {
    double x{std::log(rcens[k] / beta + 1)};
    
    if (x * a <= 1.0) {
      expm = m_exp_sum(x, m, aux_vect, a);
    }
    else {
      int n{};
      n = std::log(a * x) / std::log(2.0);
      ++n;
      
      expm = m_exp_sum(x / pow(2.0, n), m, aux_vect, a);
      
      pow2_matrix(n, expm);
    }
    aux_mat = alpha.t() * expm * e;
    density = aux_mat(0,0);
    logLh += rcweight[k] * std::log(density);
  }
  
  return logLh;
}


//' Loglikelihood of matrix-lognormal using uniformization
//' 
//' Loglikelihood for a sample.
//' 
//' @param h Positive parameter.
//' @param alpha Initial probabilities.
//' @param S Sub-intensity matrix.
//' @param beta Parameter of transformation.
//' @param obs The observations.
//' @param weight Weights of the observations.
//' @param rcens censored observations.
//' @param rcweight Weights of the censored observations.
//' 
// [[Rcpp::export]]
double logLikelihoodMlognormal_UNI(double h, arma::vec & alpha, arma::mat & S, double beta , const Rcpp::NumericVector & obs, const Rcpp::NumericVector & weight, const Rcpp::NumericVector & rcens, const Rcpp::NumericVector & rcweight) {
  if(beta < 0) return NA_REAL;
  
  arma::mat e;
  e.ones(S.n_cols, 1);
  arma::mat exit_vect = (S * (-1)) * e;
  
  arma::mat expm(size(S));
  
  double a = max_diagonal(S * (-1));
  
  int m{find_n(h, 1)};
  
  std::vector<arma::mat> aux_vect;
  
  vector_of_matrices(aux_vect, S, a, m);
  
  arma::mat aux_mat(1,1);
  
  double density{0.0};
  
  double logLh{0.0};
  
  // Non censored data
  for (int k{0}; k < obs.size(); ++k) {
    double x{pow(std::log(obs[k] + 1), beta)};
    
    if (x * a <= 1.0) {
      expm = m_exp_sum(x, m, aux_vect, a);
    }
    else {
      int n{};
      n = std::log(a * x) / std::log(2.0);
      ++n;
      
      expm = m_exp_sum(x / pow(2.0, n), m, aux_vect, a);
      
      pow2_matrix(n, expm);
    }
    aux_mat = alpha.t() * expm * exit_vect;
    density = aux_mat(0,0);
    logLh += weight[k] * (std::log(density) + std::log(beta) + (beta -1) * std::log(std::log(obs[k] + 1)) - std::log(obs[k] + 1));
  }
  // Right censored data
  for (int k{0}; k < rcens.size(); ++k) {
    double x{pow(std::log(rcens[k] + 1), beta)};
    
    if (x * a <= 1.0) {
      expm = m_exp_sum(x, m, aux_vect, a);
    }
    else {
      int n{};
      n = std::log(a * x) / std::log(2.0);
      ++n;
      
      expm = m_exp_sum(x / pow(2.0, n), m, aux_vect, a);
      
      pow2_matrix(n, expm);
    }
    aux_mat = alpha.t() * expm * e;
    density = aux_mat(0,0);
    logLh += rcweight[k] * std::log(density);
  }
  
  return logLh;
}


//' Loglikelihood of matrix-loglogistic using uniformization
//' 
//' Loglikelihood for a sample.
//' 
//' @param h Positive parameter.
//' @param alpha Initial probabilities.
//' @param S Sub-intensity matrix.
//' @param beta Parameter of transformation.
//' @param obs The observations.
//' @param weight Weights of the observations.
//' @param rcens censored observations.
//' @param rcweight Weights of the censored observations.
//' 
// [[Rcpp::export]]
double logLikelihoodMloglogistic_UNI(double h, arma::vec & alpha, arma::mat & S, Rcpp::NumericVector beta , const Rcpp::NumericVector & obs, const Rcpp::NumericVector & weight, const Rcpp::NumericVector & rcens, const Rcpp::NumericVector & rcweight) {
  if(beta[0] < 0 || beta[1] < 0) return NA_REAL;
  
  arma::mat e;
  e.ones(S.n_cols, 1);
  arma::mat exit_vect = (S * (-1)) * e;
  
  arma::mat expm(size(S));
  
  double a = max_diagonal(S * (-1));
  
  int m{find_n(h, 1)};
  
  std::vector<arma::mat> aux_vect;
  
  vector_of_matrices(aux_vect, S, a, m);
  
  arma::mat aux_mat(1,1);
  
  double density{0.0};
  
  double logLh{0.0};
  
  // Non censored data
  for (int k{0}; k < obs.size(); ++k) {
    double x{std::log(pow(obs[k] / beta[0], beta[1]) + 1)};
    
    if (x * a <= 1.0) {
      expm = m_exp_sum(x, m, aux_vect, a);
    }
    else {
      int n{};
      n = std::log(a * x) / std::log(2.0);
      ++n;
      
      expm = m_exp_sum(x / pow(2.0, n), m, aux_vect, a);
      
      pow2_matrix(n, expm);
    }
    aux_mat = alpha.t() * expm * exit_vect;
    density = aux_mat(0,0);
    logLh += weight[k] * (std::log(density) + std::log(beta[1]) - std::log(beta[0]) + (beta[1] - 1) * (std::log(obs[k]) - std::log(beta[0])) - std::log(pow(obs[k] / beta[0], beta[1]) + 1));
  }
  // Right censored data
  for (int k{0}; k < rcens.size(); ++k) {
    double x{std::log(pow(rcens[k] / beta[0], beta[1]) + 1)};
    
    if (x * a <= 1.0) {
      expm = m_exp_sum(x, m, aux_vect, a);
    }
    else {
      int n{};
      n = std::log(a * x) / std::log(2.0);
      ++n;
      
      expm = m_exp_sum(x / pow(2.0, n), m, aux_vect, a);
      
      pow2_matrix(n, expm);
    }
    aux_mat = alpha.t() * expm * e;
    density = aux_mat(0,0);
    logLh += rcweight[k] * std::log(density);
  }
  
  return logLh;
}


//' Loglikelihood of matrix-Gompertz using uniformization
//' 
//' Loglikelihood for a sample.
//' 
//' @param h Positive parameter.
//' @param alpha Initial probabilities.
//' @param S Sub-intensity matrix.
//' @param beta Parameter of transformation.
//' @param obs The observations.
//' @param weight Weights of the observations.
//' @param rcens censored observations.
//' @param rcweight Weights of the censored observations.
//' 
// [[Rcpp::export]]
double logLikelihoodMgompertz_UNI(double h, arma::vec & alpha, arma::mat & S, double beta , const Rcpp::NumericVector & obs, const Rcpp::NumericVector & weight, const Rcpp::NumericVector & rcens, const Rcpp::NumericVector & rcweight) {
  if(beta < 0) return NA_REAL;
  
  arma::mat e;
  e.ones(S.n_cols, 1);
  arma::mat exit_vect = (S * (-1)) * e;
  
  arma::mat expm(size(S));
  
  double a = max_diagonal(S * (-1));
  
  int m{find_n(h, 1)};
  
  std::vector<arma::mat> aux_vect;
  
  vector_of_matrices(aux_vect, S, a, m);
  
  arma::mat aux_mat(1,1);
  
  double density{0.0};
  
  double logLh{0.0};
  
  // Non censored data
  for (int k{0}; k < obs.size(); ++k) {
    double x{(exp(obs[k] * beta) - 1) / beta};
    
    if (x * a <= 1.0) {
      expm = m_exp_sum(x, m, aux_vect, a);
    }
    else {
      int n{};
      n = std::log(a * x) / std::log(2.0);
      ++n;
      
      expm = m_exp_sum(x / pow(2.0, n), m, aux_vect, a);
      
      pow2_matrix(n, expm);
    }
    aux_mat = alpha.t() * expm * exit_vect;
    density = aux_mat(0,0);
    logLh += weight[k] * (std::log(density) + obs[k] * beta);
  }
  // Right censored data
  for (int k{0}; k < rcens.size(); ++k) {
    double x{(exp(rcens[k] * beta) - 1) / beta};
    
    if (x * a <= 1.0) {
      expm = m_exp_sum(x, m, aux_vect, a);
    }
    else {
      int n{};
      n = std::log(a * x) / std::log(2.0);
      ++n;
      
      expm = m_exp_sum(x / pow(2.0, n), m, aux_vect, a);
      
      pow2_matrix(n, expm);
    }
    aux_mat = alpha.t() * expm * e;
    density = aux_mat(0,0);
    logLh += rcweight[k] * std::log(density);
  }
  
  return logLh;
}


//' Loglikelihood of matrix-GEV using uniformization
//' 
//' Loglikelihood for a sample.
//' 
//' @param h Positive parameter.
//' @param alpha Initial probabilities.
//' @param S Sub-intensity matrix.
//' @param beta Parameter of transformation.
//' @param obs The observations.
//' @param weight Weights of the observations.
//' @param rcens censored observations.
//' @param rcweight Weights of the censored observations.
//' 
// [[Rcpp::export]]
double logLikelihoodMgev_UNI(double h, arma::vec & alpha, arma::mat & S, Rcpp::NumericVector beta , const Rcpp::NumericVector & obs, const Rcpp::NumericVector & weight, const Rcpp::NumericVector & rcens, const Rcpp::NumericVector & rcweight) {
  if(beta[1] < 0) return NA_REAL;
  
  arma::mat e;
  e.ones(S.n_cols, 1);
  arma::mat exit_vect = (S * (-1)) * e;
  
  arma::mat expm(size(S));
  
  double a = max_diagonal(S * (-1));
  
  int m{find_n(h, 1)};
  
  std::vector<arma::mat> aux_vect;
  
  vector_of_matrices(aux_vect, S, a, m);
  
  arma::mat aux_mat(1,1);
  
  double density{0.0};
  
  double logLh{0.0};
  
  if (beta[2] == 0) {
    // Non censored data
    for (int k{0}; k < obs.size(); ++k) {
      double x{exp(-(obs[k] - beta[0]) / beta[1])};
      
      if (x * a <= 1.0) {
        expm = m_exp_sum(x, m, aux_vect, a);
      }
      else {
        int n{};
        n = std::log(a * x) / std::log(2.0);
        ++n;
        
        expm = m_exp_sum(x / pow(2.0, n), m, aux_vect, a);
        
        pow2_matrix(n, expm);
      }
      aux_mat = alpha.t() * expm * exit_vect;
      density = aux_mat(0,0);
      logLh += weight[k] * (std::log(density) - std::log(beta[1]) - (obs[k] - beta[0]) / beta[1]);
    }
    // Right censored data
    for (int k{0}; k < rcens.size(); ++k) {
      double x{exp(-(rcens[k] - beta[0]) / beta[1])};
      
      if (x * a <= 1.0) {
        expm = m_exp_sum(x, m, aux_vect, a);
      }
      else {
        int n{};
        n = std::log(a * x) / std::log(2.0);
        ++n;
        
        expm = m_exp_sum(x / pow(2.0, n), m, aux_vect, a);
        
        pow2_matrix(n, expm);
      }
      aux_mat = alpha.t() * expm * e;
      density = aux_mat(0,0);
      logLh += rcweight[k] * std::log(density);
    }
  }
  else {
    // Non censored data
    for (int k{0}; k < obs.size(); ++k) {
      double x{pow(1 + (beta[2] / beta[1]) * (obs[k] - beta[0]) , - 1 / beta[2])};
      
      if (x * a <= 1.0) {
        expm = m_exp_sum(x, m, aux_vect, a);
      }
      else {
        int n{};
        n = std::log(a * x) / std::log(2.0);
        ++n;
        
        expm = m_exp_sum(x / pow(2.0, n), m, aux_vect, a);
        
        pow2_matrix(n, expm);
      }
      aux_mat = alpha.t() * expm * exit_vect;
      density = aux_mat(0,0);
      logLh += weight[k] * (std::log(density) - std::log(beta[1]) - (1 + 1 / beta[2]) * std::log(1 + (beta[2] / beta[1]) * (obs[k] - beta[0])));
    }
    // Right censored data
    for (int k{0}; k < rcens.size(); ++k) {
      double x{pow(1 + (beta[2] / beta[1]) * (rcens[k] - beta[0]) , - 1 / beta[2])};
      
      if (x * a <= 1.0) {
        expm = m_exp_sum(x, m, aux_vect, a);
      }
      else {
        int n{};
        n = std::log(a * x) / std::log(2.0);
        ++n;
        
        expm = m_exp_sum(x / pow(2.0, n), m, aux_vect, a);
        
        pow2_matrix(n, expm);
      }
      aux_mat = alpha.t() * expm * e;
      density = aux_mat(0,0);
      logLh += rcweight[k] * std::log(density);
    }
  }
  return logLh;
}


/////////////////////////////////////////////////////////
// Scaled versions of loglikelihoods (for regression)  //
/////////////////////////////////////////////////////////

//' Loglikelihood of PI with phase-type using uniformization
//' 
//' Loglikelihood for a sample.
//' 
//' @param h Positive parameter.
//' @param alpha Initial probabilities.
//' @param S Sub-intensity matrix.
//' @param obs The observations.
//' @param weight Weights of the observations.
//' @param rcens Censored observations.
//' @param rcweight Weights of the censored observations.
//' @param scale1 Scale for observations.
//' @param scale2 Scale for censored observations.
//' 
// [[Rcpp::export]]
double logLikelihoodPH_UNIs(double h, arma::vec & alpha, arma::mat & S, const Rcpp::NumericVector & obs, const Rcpp::NumericVector & weight, const Rcpp::NumericVector & rcens, const Rcpp::NumericVector & rcweight, const Rcpp::NumericVector & scale1, const Rcpp::NumericVector & scale2) {
  arma::mat e;
  e.ones(S.n_cols, 1);
  arma::mat exit_vect = (S * (-1)) * e;
  
  arma::mat expm(size(S));
  
  double a = max_diagonal(S * (-1));
  
  int m{find_n(h, 1)};
  
  std::vector<arma::mat> aux_vect;
  
  vector_of_matrices(aux_vect, S, a, m);
  
  arma::mat aux_mat(1,1);
  
  double density{0.0};
  
  double logLh{0.0};
  
  // Non censored data
  for (int k{0}; k < obs.size(); ++k) {
    double x{scale1[k] * obs[k]};
    
    if (x * a <= 1.0) {
      expm = m_exp_sum(x, m, aux_vect, a);
    }
    else {
      int n{};
      n = std::log(a * x) / std::log(2.0);
      ++n;
      
      expm = m_exp_sum(x / pow(2.0, n), m, aux_vect, a);
      
      pow2_matrix(n, expm);
    }
    aux_mat = alpha.t() * expm * exit_vect;
    density = aux_mat(0,0);
    logLh += weight[k] * (std::log(density) + std::log(scale1[k]));
  }
  // Right censored data
  for (int k{0}; k < rcens.size(); ++k) {
    double x{scale2[k] * rcens[k]};
    
    if (x * a <= 1.0) {
      expm = m_exp_sum(x, m, aux_vect, a);
    }
    else {
      int n{};
      n = std::log(a * x) / std::log(2.0);
      ++n;
      
      expm = m_exp_sum(x / pow(2.0, n), m, aux_vect, a);
      
      pow2_matrix(n, expm);
    }
    aux_mat = alpha.t() * expm * e;
    density = aux_mat(0,0);
    logLh += rcweight[k] * std::log(density);
  }
  
  return logLh;
}


//' Loglikelihood of PI with matrix-Weibull using uniformization
//' 
//' Loglikelihood for a sample.
//' 
//' @param h Positive parameter.
//' @param alpha Initial probabilities.
//' @param S Sub-intensity matrix.
//' @param beta Parameter of transformation.
//' @param obs The observations.
//' @param weight Weights of the observations.
//' @param rcens Censored observations.
//' @param rcweight Weights of the censored observations.
//' @param scale1 Scale for observations.
//' @param scale2 Scale for censored observations.
//' 
// [[Rcpp::export]]
double logLikelihoodMweibull_UNIs(double h, arma::vec & alpha, arma::mat & S, double beta , const Rcpp::NumericVector & obs, const Rcpp::NumericVector & weight, const Rcpp::NumericVector & rcens, const Rcpp::NumericVector & rcweight, const Rcpp::NumericVector & scale1, const Rcpp::NumericVector & scale2) {
  if(beta < 0) return NA_REAL;
  
  arma::mat e;
  e.ones(S.n_cols, 1);
  arma::mat exit_vect = (S * (-1)) * e;
  
  arma::mat expm(size(S));
  
  double a = max_diagonal(S * (-1));
  
  int m{find_n(h, 1)};
  
  std::vector<arma::mat> aux_vect;
  
  vector_of_matrices(aux_vect, S, a, m);
  
  arma::mat aux_mat(1,1);
  
  double density{0.0};
  
  double logLh{0.0};
  
  // Non censored data
  for (int k{0}; k < obs.size(); ++k) {
    double x{scale1[k] * pow(obs[k], beta)};
    
    if (x * a <= 1.0) {
      expm = m_exp_sum(x, m, aux_vect, a);
    }
    else {
      int n{};
      n = std::log(a * x) / std::log(2.0);
      ++n;
      
      expm = m_exp_sum(x / pow(2.0, n), m, aux_vect, a);
      
      pow2_matrix(n, expm);
    }
    aux_mat = alpha.t() * expm * exit_vect;
    density = aux_mat(0,0);
    logLh += weight[k] * (std::log(density) + std::log(scale1[k]) + std::log(beta) + (beta -1) * std::log(obs[k]));
  }
  // Right censored data
  for (int k{0}; k < rcens.size(); ++k) {
    double x{scale2[k] * pow(rcens[k], beta)};
    
    if (x * a <= 1.0) {
      expm = m_exp_sum(x, m, aux_vect, a);
    }
    else {
      int n{};
      n = std::log(a * x) / std::log(2.0);
      ++n;
      
      expm = m_exp_sum(x / pow(2.0, n), m, aux_vect, a);
      
      pow2_matrix(n, expm);
    }
    aux_mat = alpha.t() * expm * e;
    density = aux_mat(0,0);
    logLh += rcweight[k] * std::log(density);
  }
  
  return logLh;
}


//' Loglikelihood of PI with matrix-Pareto using uniformization
//' 
//' Loglikelihood for a sample.
//' 
//' @param h Positive parameter.
//' @param alpha Initial probabilities.
//' @param S Sub-intensity matrix.
//' @param beta Parameter of transformation.
//' @param obs The observations.
//' @param weight Weights of the observations.
//' @param rcens Censored observations.
//' @param rcweight Weights of the censored observations.
//' @param scale1 Scale for observations.
//' @param scale2 Scale for censored observations.
//' 
// [[Rcpp::export]]
double logLikelihoodMpareto_UNIs(double h, arma::vec & alpha, arma::mat & S, double beta , const Rcpp::NumericVector & obs, const Rcpp::NumericVector & weight, const Rcpp::NumericVector & rcens, const Rcpp::NumericVector & rcweight, const Rcpp::NumericVector & scale1, const Rcpp::NumericVector & scale2) {
  if(beta < 0) return NA_REAL;
  
  arma::mat e;
  e.ones(S.n_cols, 1);
  arma::mat exit_vect = (S * (-1)) * e;
  
  arma::mat expm(size(S));
  
  double a = max_diagonal(S * (-1));
  
  int m{find_n(h, 1)};
  
  std::vector<arma::mat> aux_vect;
  
  vector_of_matrices(aux_vect, S, a, m);
  
  arma::mat aux_mat(1,1);
  
  double density{0.0};
  
  double logLh{0.0};
  
  // Non censored data
  for (int k{0}; k < obs.size(); ++k) {
    double x{scale1[k] * std::log(obs[k] / beta + 1)};
    
    if (x * a <= 1.0) {
      expm = m_exp_sum(x, m, aux_vect, a);
    }
    else {
      int n{};
      n = std::log(a * x) / std::log(2.0);
      ++n;
      
      expm = m_exp_sum(x / pow(2.0, n), m, aux_vect, a);
      
      pow2_matrix(n, expm);
    }
    aux_mat = alpha.t() * expm * exit_vect;
    density = aux_mat(0,0);
    logLh += weight[k] * (std::log(density) + std::log(scale1[k]) - std::log(obs[k] + beta));
  }
  // Right censored data
  for (int k{0}; k < rcens.size(); ++k) {
    double x{scale2[k] * std::log(rcens[k] / beta + 1)};
    
    if (x * a <= 1.0) {
      expm = m_exp_sum(x, m, aux_vect, a);
    }
    else {
      int n{};
      n = std::log(a * x) / std::log(2.0);
      ++n;
      
      expm = m_exp_sum(x / pow(2.0, n), m, aux_vect, a);
      
      pow2_matrix(n, expm);
    }
    aux_mat = alpha.t() * expm * e;
    density = aux_mat(0,0);
    logLh += rcweight[k] * std::log(density);
  }
  
  return logLh;
}


//' Loglikelihood of PI with matrix-lognormal using uniformization
//' 
//' Loglikelihood for a sample.
//' 
//' @param h Positive parameter.
//' @param alpha Initial probabilities.
//' @param S Sub-intensity matrix.
//' @param beta Parameter of transformation.
//' @param obs The observations.
//' @param weight Weights of the observations.
//' @param rcens Censored observations.
//' @param rcweight Weights of the censored observations.
//' @param scale1 Scale for observations.
//' @param scale2 Scale for censored observations.
//' 
// [[Rcpp::export]]
double logLikelihoodMlognormal_UNIs(double h, arma::vec & alpha, arma::mat & S, double beta , const Rcpp::NumericVector & obs, const Rcpp::NumericVector & weight, const Rcpp::NumericVector & rcens, const Rcpp::NumericVector & rcweight, const Rcpp::NumericVector & scale1, const Rcpp::NumericVector & scale2) {
  if(beta < 0) return NA_REAL;
  
  arma::mat e;
  e.ones(S.n_cols, 1);
  arma::mat exit_vect = (S * (-1)) * e;
  
  arma::mat expm(size(S));
  
  double a = max_diagonal(S * (-1));
  
  int m{find_n(h, 1)};
  
  std::vector<arma::mat> aux_vect;
  
  vector_of_matrices(aux_vect, S, a, m);
  
  arma::mat aux_mat(1,1);
  
  double density{0.0};
  
  double logLh{0.0};
  
  // Non censored data
  for (int k{0}; k < obs.size(); ++k) {
    double x{scale1[k] * pow(std::log(obs[k] + 1), beta)};
    
    if (x * a <= 1.0) {
      expm = m_exp_sum(x, m, aux_vect, a);
    }
    else {
      int n{};
      n = std::log(a * x) / std::log(2.0);
      ++n;
      
      expm = m_exp_sum(x / pow(2.0, n), m, aux_vect, a);
      
      pow2_matrix(n, expm);
    }
    aux_mat = alpha.t() * expm * exit_vect;
    density = aux_mat(0,0);
    logLh += weight[k] * (std::log(density) + std::log(scale1[k]) + std::log(beta) + (beta - 1) * std::log(std::log(obs[k] + 1)) - std::log(obs[k] + 1));
  }
  // Right censored data
  for (int k{0}; k < rcens.size(); ++k) {
    double x{scale2[k] * pow(std::log(rcens[k] + 1), beta)};
    
    if (x * a <= 1.0) {
      expm = m_exp_sum(x, m, aux_vect, a);
    }
    else {
      int n{};
      n = std::log(a * x) / std::log(2.0);
      ++n;
      
      expm = m_exp_sum(x / pow(2.0, n), m, aux_vect, a);
      
      pow2_matrix(n, expm);
    }
    aux_mat = alpha.t() * expm * e;
    density = aux_mat(0,0);
    logLh += rcweight[k] * std::log(density);
  }
  
  return logLh;
}


//' Loglikelihood of PI with matrix-loglogistic using uniformization
//' 
//' Loglikelihood for a sample.
//' 
//' @param h Positive parameter.
//' @param alpha Initial probabilities.
//' @param S Sub-intensity matrix.
//' @param beta Parameter of transformation.
//' @param obs The observations.
//' @param weight Weights of the observations.
//' @param rcens Censored observations.
//' @param rcweight Weights of the censored observations.
//' @param scale1 Scale for observations.
//' @param scale2 Scale for censored observations.
//' 
// [[Rcpp::export]]
double logLikelihoodMloglogistic_UNIs(double h, arma::vec & alpha, arma::mat & S, Rcpp::NumericVector beta , const Rcpp::NumericVector & obs, const Rcpp::NumericVector & weight, const Rcpp::NumericVector & rcens, const Rcpp::NumericVector & rcweight, const Rcpp::NumericVector & scale1, const Rcpp::NumericVector & scale2) {
  if(beta[0] < 0 || beta[1] < 0) return NA_REAL;
  
  arma::mat e;
  e.ones(S.n_cols, 1);
  arma::mat exit_vect = (S * (-1)) * e;
  
  arma::mat expm(size(S));
  
  double a = max_diagonal(S * (-1));
  
  int m{find_n(h, 1)};
  
  std::vector<arma::mat> aux_vect;
  
  vector_of_matrices(aux_vect, S, a, m);
  
  arma::mat aux_mat(1,1);
  
  double density{0.0};
  
  double logLh{0.0};
  
  // Non censored data
  for (int k{0}; k < obs.size(); ++k) {
    double x{scale1[k] * std::log(pow(obs[k] / beta[0], beta[1]) + 1)};
    
    if (x * a <= 1.0) {
      expm = m_exp_sum(x, m, aux_vect, a);
    }
    else {
      int n{};
      n = std::log(a * x) / std::log(2.0);
      ++n;
      
      expm = m_exp_sum(x / pow(2.0, n), m, aux_vect, a);
      
      pow2_matrix(n, expm);
    }
    aux_mat = alpha.t() * expm * exit_vect;
    density = aux_mat(0,0);
    logLh += weight[k] * (std::log(density) + std::log(scale1[k]) + std::log(beta[1]) - std::log(beta[0]) + (beta[1] - 1) * (std::log(obs[k]) - std::log(beta[0])) - std::log(pow(obs[k] / beta[0], beta[1]) + 1));
  }
  // Right censored data
  for (int k{0}; k < rcens.size(); ++k) {
    double x{scale2[k] * std::log(pow(rcens[k] / beta[0], beta[1]) + 1)};
    
    if (x * a <= 1.0) {
      expm = m_exp_sum(x, m, aux_vect, a);
    }
    else {
      int n{};
      n = std::log(a * x) / std::log(2.0);
      ++n;
      
      expm = m_exp_sum(x / pow(2.0, n), m, aux_vect, a);
      
      pow2_matrix(n, expm);
    }
    aux_mat = alpha.t() * expm * e;
    density = aux_mat(0,0);
    logLh += rcweight[k] * std::log(density);
  }
  
  return logLh;
}


//' Loglikelihood of PI with matrix-Gompertz using Uniformization
//' 
//' Loglikelihood for a sample.
//' 
//' @param h Positive parameter.
//' @param alpha Initial probabilities.
//' @param S Sub-intensity matrix.
//' @param beta Parameter of transformation.
//' @param obs The observations.
//' @param weight Weights of the observations.
//' @param rcens Censored observations.
//' @param rcweight Weights of the censored observations.
//' @param scale1 Scale for observations.
//' @param scale2 Scale for censored observations.
//' 
// [[Rcpp::export]]
double logLikelihoodMgompertz_UNIs(double h, arma::vec & alpha, arma::mat & S, double beta , const Rcpp::NumericVector & obs, const Rcpp::NumericVector & weight, const Rcpp::NumericVector & rcens, const Rcpp::NumericVector & rcweight, const Rcpp::NumericVector & scale1, const Rcpp::NumericVector & scale2) {
  if(beta < 0) return NA_REAL;
  
  arma::mat e;
  e.ones(S.n_cols, 1);
  arma::mat exit_vect = (S * (-1)) * e;
  
  arma::mat expm(size(S));
  
  double a = max_diagonal(S * (-1));
  
  int m{find_n(h, 1)};
  
  std::vector<arma::mat> aux_vect;
  
  vector_of_matrices(aux_vect, S, a, m);
  
  arma::mat aux_mat(1,1);
  
  double density{0.0};
  
  double logLh{0.0};
  
  // Non censored data
  for (int k{0}; k < obs.size(); ++k) {
    double x{scale1[k] * (exp(obs[k] * beta) - 1) / beta};
    
    if (x * a <= 1.0) {
      expm = m_exp_sum(x, m, aux_vect, a);
    }
    else {
      int n{};
      n = std::log(a * x) / std::log(2.0);
      ++n;
      
      expm = m_exp_sum(x / pow(2.0, n), m, aux_vect, a);
      
      pow2_matrix(n, expm);
    }
    aux_mat = alpha.t() * expm * exit_vect;
    density = aux_mat(0,0);
    logLh += weight[k] * (std::log(density) + std::log(scale1[k]) + obs[k] * beta);
  }
  // Right censored data
  for (int k{0}; k < rcens.size(); ++k) {
    double x{scale2[k] * (exp(rcens[k] * beta) - 1) / beta};
    
    if (x * a <= 1.0) {
      expm = m_exp_sum(x, m, aux_vect, a);
    }
    else {
      int n{};
      n = std::log(a * x) / std::log(2.0);
      ++n;
      
      expm = m_exp_sum(x / pow(2.0, n), m, aux_vect, a);
      
      pow2_matrix(n, expm);
    }
    aux_mat = alpha.t() * expm * e;
    density = aux_mat(0,0);
    logLh += rcweight[k] * std::log(density);
  }
  
  return logLh;
}

////////////////// Additional function for speed:

//' expm terms of phase-type likelihood using uniformization
//' 
//' @param h Positive parameter.
//' @param S Sub-intensity matrix.
//' @param obs The observations.
//' 
// [[Rcpp::export]]
Rcpp::List expm_terms(double h, arma::mat & S, const Rcpp::NumericVector & obs) {
  Rcpp::List out(obs.size());
  arma::mat expm(size(S));
  
  double a = max_diagonal(S * (-1));
  int m{find_n(h, 1)};
  
  std::vector<arma::mat> aux_vect;
  vector_of_matrices(aux_vect, S, a, m);
  
  for (int k{0}; k < obs.size(); ++k) {
    double x{obs[k]};
    if (x * a <= 1.0) {
      expm = m_exp_sum(x, m, aux_vect, a);
    }
    else {
      int n{};
      n = std::log(a * x) / std::log(2.0);
      ++n;
      expm = m_exp_sum(x / pow(2.0, n), m, aux_vect, a);
      pow2_matrix(n, expm);
    }
    out[k] = expm;
  }
  return out;
}
