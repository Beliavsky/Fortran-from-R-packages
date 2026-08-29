#include "m_exp.h"
# include <RcppArmadillo.h>
// [[ Rcpp :: depends ( RcppArmadillo )]]

////////////////////////////////////////////
//    EM using Pade approximation        ///
////////////////////////////////////////////

//' Computes the elements S^n / n! until given value of n
//' 
//' @param vect A vector.
//' @param S Sub-intensity matrix.
//' @param vect_size Size of vector.
//' 
// [[Rcpp::export]]
void vector_of_matrices_2(std::vector<arma::mat> & vect, const arma::mat & S, int vect_size) {
  arma::mat I;
  I.eye(size(S));
  
  vect.push_back(I);
  
  for (int k{1}; k <= vect_size; ++k) {
    vect.push_back(S * vect[k - 1]);
  }
}


//' EM for phase-type distributions using Pade approximation for matrix exponential
//' 
//' @param h Nuisance parameter.
//' @param alpha Initial probabilities.
//' @param S Sub-intensity matrix.
//' @param obs The observations.
//' @param weight The weights for the observations.
//' @param rcens Censored observations.
//' @param rcweight The weights for the censored observations.
//' 
// [[Rcpp::export]]
void EMstep_PADE(double h, arma::vec & alpha,  arma::mat & S, const Rcpp::NumericVector & obs, const Rcpp::NumericVector & weight, const Rcpp::NumericVector & rcens, const Rcpp::NumericVector & rcweight) {
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
  
  double JNorm{inf_norm(J)};
  
  std::vector<arma::mat> aux_vect;
  
  vector_of_matrices_2(aux_vect, J, 6);
  
  arma::mat X(2 * p,2 * p);
  arma::mat D(2 * p,2 * p);
  
  const int q{6};
  int s{};
  double xmod{};
  double c{};
  
  double sum_weights{0.0};
  double density{0.0};
  
  //E-step
  //  Uncensored data
  for (int k{0}; k < obs.size(); ++k) {
    sum_weights += weight[k];
    
    // Matrix exponential
    int pind{1};
    int ee{static_cast<int>(log2(JNorm  * obs[k])) + 1};
    s = std::max(0, ee + 1);
    xmod = obs[k] / pow(2.0, s);
    c = 0.5;
    X = aux_vect[1] * (c * xmod);
    
    J = aux_vect[0] + X;
    D = aux_vect[0] - X;
    
    for (int l{2}; l <= q; ++l) {
      c = c * static_cast<double>(q - l + 1) / static_cast<double>(l * (2 * q - l + 1));
      X = aux_vect[l] * (c * pow(xmod,l));
      J = J + X;
      if (pind) {
        D =  D + X;
      }
      else {
        D = D - X;
      }
      pind = !pind;
    }
    J = inv(D) * J;
    for (int l = 1; l <= s; ++l) {
      J = J * J;
    }
    
    // Separate matrix
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
    
    //E-step
    for (int i{0}; i < p; ++i) {
      Bmean(i,0) += alpha[i] * bvector(i,0) * weight[k] / density;
      Nmean(i,p) += avector(0,i) * exit_vect(i,0) * weight[k] / density;
      Zmean(i,0) += cmatrix(i,i) * weight[k] / density;
      for (int j{0}; j < p; ++j) {
        Nmean(i,j) += S(i,j) * cmatrix(j,i) * weight[k] / density;
      }
    }
  }
  //  Right-Censored Data
  double sum_censored{0.0};
  if (rcens.size() > 0) {
    s_prod_alpha = e * alpha.t();
    J = matrix_vanloan(S, S, s_prod_alpha);
    JNorm = inf_norm(J);
    aux_vect.clear();
    vector_of_matrices_2(aux_vect, J, 6);
  }
  for (int k{0}; k < rcens.size(); ++k) {
    sum_censored += rcweight[k];
    
    // Matrix exponential
    int pind{1};
    int ee{static_cast<int>(log2(JNorm  * rcens[k])) + 1};
    s = std::max(0, ee + 1);
    xmod = rcens[k] / pow(2.0, s);
    c = 0.5;
    X = aux_vect[1] * (c * xmod);
    
    J = aux_vect[0] + X;
    D = aux_vect[0] - X;
    
    for (int l{2}; l <= q; ++l) {
      c = c * static_cast<double>(q - l + 1) / static_cast<double>(l * (2 * q - l + 1));
      X = aux_vect[l] * (c * pow(xmod,l));
      J = J + X;
      if (pind) {
        D =  D + X;
      }
      else {
        D = D - X;
      }
      pind = !pind;
    }
    J = inv(D) * J;
    for (int l = 1; l <= s; ++l) {
      J = J * J;
    }
    
    // Separate matrix
    for (int i{0}; i < p; ++i) {
      for (int j{0}; j < p; ++j) {
        aux_exp(i,j) = J(i,j);
        cmatrix(i,j) = J(i,j + p);
      }
    }
    
    bvector = aux_exp * e;
    aux_mat = alpha.t() * bvector;
    density = aux_mat(0,0);
    
    //E-step
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


//' EM for PH-MoE
//' 
//' No recycling of information
//' 
//' @param alpha Initial probabilities.
//' @param S Sub-intensity matrix.
//' @param obs The observations.
//' @param weight The weights for the observations.
//' @param rcens Censored observations.
//' @param rcweight The weights for the censored observations.
//' 
// [[Rcpp::export]]
Rcpp::List EMstep_MoE_PADE(arma::mat & alpha,  arma::mat & S, const Rcpp::NumericVector & obs, const Rcpp::NumericVector & weight, const Rcpp::NumericVector & rcens, const Rcpp::NumericVector & rcweight) {
  unsigned p{S.n_rows}; 
  
  arma::mat Bmatrix(obs.size(), p);
  arma::mat Bmatrix_cens;
  
  arma::mat e; 
  e.ones(S.n_cols, 1);
  arma::mat exit_vect = (S * (-1)) * e;
  
  arma::rowvec alpha_aux(alpha.row(0));
  
  arma::mat Zmean = arma::zeros(p,1);
  arma::mat Nmean = arma::zeros(p,p + 1);
  
  arma::mat avector(1,p);
  arma::mat bvector(p,1);
  arma::mat cmatrix(p,p);
  arma::mat aux_exp(p,p);
  
  arma::mat aux_mat(1,1);
  
  arma::mat J(2 * p,2 * p);
  arma::mat s_prod_alpha(p,p);
  
  double density{0.0};
  
  //E-step
  //  Uncensored data
  for (int k{0}; k < obs.size(); ++k) {
    alpha_aux = alpha.row(k);
    s_prod_alpha = exit_vect * alpha_aux;
    J = matrix_exponential(matrix_vanloan(S, S, s_prod_alpha) * obs[k]);
    
    for (int i{0}; i < p; ++i) {
      for (int j{0}; j < p; ++j) {
        aux_exp(i,j) = J(i,j);
        cmatrix(i,j) = J(i,j + p);
      }
    }
    
    avector = alpha_aux * aux_exp;
    bvector = aux_exp * exit_vect;
    aux_mat = alpha_aux * bvector;
    density = aux_mat(0,0);
    
    //E-step
    for (int i{0}; i < p; ++i) {
      Bmatrix(k,i) = alpha_aux[i] * bvector(i,0) * weight[k] / density;
      Nmean(i,p) += avector(0,i) * exit_vect(i,0) * weight[k] / density;
      Zmean(i,0) += cmatrix(i,i) * weight[k] / density;
      for (int j{0}; j < p; ++j) {
        Nmean(i,j) += S(i,j) * cmatrix(j,i) * weight[k] / density;
      }
    }
  }
  
  //  Right-Censored Data
  if (rcens.size() > 0) {
    Bmatrix_cens.set_size(rcens.size(), p);
  }
  for (int k{0}; k < rcens.size(); ++k) {
    alpha_aux = alpha.row(k + obs.size());
    s_prod_alpha = e * alpha_aux;
    J = matrix_exponential(matrix_vanloan(S, S, s_prod_alpha) * rcens[k]);
    
    for (int i{0}; i < p; ++i) {
      for (int j{0}; j < p; ++j) {
        aux_exp(i,j) = J(i,j);
        cmatrix(i,j) = J(i,j + p);
      }
    }
    
    bvector = aux_exp * e;
    aux_mat = alpha_aux * bvector;
    density = aux_mat(0,0);
    
    //E-step
    for (int i{0}; i < p; ++i) {
      Bmatrix_cens(k,i) = alpha_aux[i] * bvector(i,0) * rcweight[k] / density;
      Zmean(i,0) += cmatrix(i,i) * rcweight[k] / density;
      for (int j{0}; j < p; ++j) {
        Nmean(i,j) += S(i,j) * cmatrix(j,i) * rcweight[k] / density;
      }
    }
  }
  
  // M step
  for (int i{0}; i < p; ++i) {
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
  return Rcpp::List::create(
    Rcpp::Named("obs") = Bmatrix,
    Rcpp::Named("rcens") = Bmatrix_cens
  );
}


//' EM for bivariate phase-type distributions using Pade for matrix exponential
//'
//' @param alpha Initial probabilities.
//' @param S11 Sub-intensity.
//' @param S12 A matrix.
//' @param S22 Sub-intensity.
//' @param obs The observations.
//' @param weight The weights for the observations.
//' @return Fitted alpha, S11, S12 and S22 after one iteration.
//' 
// [[Rcpp::export]]
void EMstep_bivph(arma::vec & alpha, arma::mat & S11, arma::mat & S12, arma::mat & S22, const Rcpp::NumericMatrix & obs, const Rcpp::NumericVector & weight) {
  unsigned p1{S11.n_rows};
  unsigned p2{S22.n_rows};
  unsigned p{p1 + p2};
  
  double density{0};
  double aux;
  
  arma::mat Bmean = arma::zeros(p1,1);
  arma::mat Zmean = arma::zeros(p,1);
  arma::mat Nmean = arma::zeros(p,p + 1);
  
  arma::mat bmatrix1(p1,p1);
  arma::mat cmatrix1(p1,p1);
  arma::mat aux_exp1(p1,p1);
  
  arma::mat bmatrix2(p2,p2);
  arma::mat cmatrix2(p2,p2);
  arma::mat aux_exp2(p2,p2);
  
  arma::mat g1(2 * p1,2 * p1);
  arma::mat g2(2 * p2,2 * p2);
  
  arma::mat e;
  e.ones(p2, 1);
  arma::mat exit_vec = (S22 * (-1)) * e;
  
  arma::mat aux_matrix1(p1,1);
  arma::mat aux_matrix2(p2,p1);
  arma::mat aux_matrix3(1,p2);
  
  arma::mat aux_mat(1,1);
  
  double sum_weights{0.0};
  //E step
  for (int k{0}; k < obs.nrow(); ++k) {
    sum_weights += weight[k];
    
    aux_exp2 = matrix_exponential(S22 * obs(k,1));
    bmatrix1 = S12 * aux_exp2 * exit_vec * alpha.t();
    
    g1 = matrix_exponential(matrix_vanloan(S11, S11, bmatrix1) * obs(k,0));
    
    for (int i{0}; i < p1; ++i) {
      for (int j{0}; j < p1; ++j) {
        aux_exp1(i,j) = g1(i,j);
        cmatrix1(i,j) = g1(i,j + p1);
      }
    }
    
    bmatrix2 = exit_vec * alpha.t() * aux_exp1 * S12;
    
    g2 = matrix_exponential(matrix_vanloan(S22, S22, bmatrix2) * obs(k,1));
    
    for (int i{0}; i < p2; ++i) {
      for (int j{0}; j < p2; ++j) {
        cmatrix2(i,j) = g2(i,j + p2);
      }
    }
    
    aux_mat = alpha.t() * aux_exp1 * S12 * aux_exp2 * exit_vec;
    density = aux_mat(0,0);
    
    //E-step
    aux_matrix1 = aux_exp1 * S12 * aux_exp2 * exit_vec;
    aux_matrix2 =  aux_exp2 * exit_vec * alpha.t() * aux_exp1;
    
    for (int i{0}; i < p1; ++i) {
      aux = aux_matrix1(i,0);
      Bmean(i,0) += alpha[i] * aux * weight[k] / density;
      Zmean(i,0) += cmatrix1(i,i) * weight[k] / density;
      for (int j{0}; j < p1; ++j) {
        Nmean(i,j) += S11(i,j) * cmatrix1(j,i) * weight[k] / density;
      }
      for (int j{0}; j < p2; ++j) {
        aux = aux_matrix2(j,i);
        Nmean(i,j + p1) += S12(i,j) * aux * weight[k] / density;
      }
    }
    
    aux_matrix3 = alpha.t() * aux_exp1 * S12 * aux_exp2;
    
    for (int i{0}; i < p2; ++i) {
      Zmean(i + p1,0) += cmatrix2(i,i) * weight[k] / density;
      aux = aux_matrix3(0,i);
      Nmean(i + p1,p) += aux * exit_vec(i,0) * weight[k] / density;
      for (int j{0}; j < p2; ++j){
        Nmean(i + p1,j + p1) += S22(i,j) * cmatrix2(j,i) * weight[k] / density;
      }
    }
  }
  
  // M step
  for (int i{0}; i < p1; ++i) {
    alpha[i] = Bmean(i,0) / sum_weights;
    if (alpha[i] < 0) {
      alpha[i] = 0;
    }
    S11(i,i) = 0;
    for (int j{0}; j < p1; ++j) {
      if (i != j) {
        S11(i,j) = Nmean(i,j) / Zmean(i,0);
        if (S11(i,j) < 0) {
          S11(i,j) = 0;
        }
        S11(i,i) -= S11(i,j);
      }
    }
    for (int j{0}; j < p2; ++j) {
      S12(i,j) = Nmean(i,j + p1) / Zmean(i,0);
      if (S12(i,j) < 0){
        S12(i,j) = 0;
      }
      S11(i,i) -= S12(i,j);
    }
  }
  for (int i{0}; i < p2; ++i) {
    exit_vec(i,0) = Nmean(i + p1,p) / Zmean(i + p1,0);
    if (exit_vec(i,0) < 0) {
      exit_vec(i,0) = 0;
    }
    S22(i,i) = -exit_vec(i,0);
    for (int j{0}; j < p2; ++j) {
      if (i != j) {
        S22(i,j) = Nmean(i + p1,j + p1) / Zmean(i + p1,0);
        if (S22(i,j) < 0){
          S22(i,j) = 0;
        }
        S22(i,i) -= S22(i,j);
      }
    }
  }
}


////////////////////////////////////////////
// Loglikelihoods
////////////////////////////////////////////

//' Loglikelihood for Bivariate PH
//' 
//' @param alpha Vector of initial probabilities.
//' @param S11 Sub-intensity matrix.
//' @param S12 Matrix.
//' @param S22 Sub-intensity matrix.
//' @param obs The observations.
//' @param weight The weights of the observations.
//' 
// [[Rcpp::export]]
double logLikelihoodbivPH(arma::vec alpha, arma::mat S11, arma::mat S12, arma::mat S22, const Rcpp::NumericMatrix & obs, const Rcpp::NumericVector & weight) {
  arma::mat e;
  e.ones(S22.n_cols, 1);
  arma::mat exit_vect = (S22 * (-1)) * e;
  
  arma::mat aux_mat(1,1);
  
  double density{0.0};
  
  double logLh{0.0};
  
  // Non censored data
  for (int k{0}; k < obs.nrow(); ++k) {
    aux_mat = alpha.t() * matrix_exponential(S11 * obs(k,0)) * S12 * matrix_exponential(S22 * obs(k,1)) * exit_vect;
    density = aux_mat(0,0);
    logLh += weight[k] * std::log(density);
  }
  
  return logLh;
}


//' Loglikelihood for PH-MoE
//' 
//' @param alpha1 Initial probabilities for non-censored data.
//' @param alpha2 Initial probabilities for censored data.
//' @param S Sub-intensity matrix. 
//' @param obs The observations.
//' @param weight The weights of the observations.
//' @param rcens Censored observations.
//' @param rcweight The weights of the censored observations.
//' 
// [[Rcpp::export]]
double logLikelihoodPH_MoE(arma::mat & alpha1, arma::mat & alpha2, arma::mat & S, const Rcpp::NumericVector & obs, const Rcpp::NumericVector & weight, const Rcpp::NumericVector & rcens, const Rcpp::NumericVector & rcweight) {
  arma::mat e; 
  e.ones(S.n_cols, 1);
  arma::mat exit_vect = (S * (-1)) * e;
  
  arma::mat aux_mat(1,1);
  
  double density{0.0};
  
  double logLh{0.0};
  
  // Non censored data
  for (int k{0}; k < obs.size(); ++k) {
    arma::rowvec alpha_aux(alpha1.row(k));
    aux_mat = alpha_aux * matrix_exponential(S * obs[k]) * exit_vect;
    density = aux_mat(0,0);
    logLh += weight[k] * std::log(density);
  }
  //Right censored data
  for (int k{0}; k < rcens.size(); ++k) {
    arma::rowvec alpha_aux(alpha2.row(k));
    aux_mat = alpha_aux * matrix_exponential(S * rcens[k])  * e;
    density = aux_mat(0,0);
    logLh += rcweight[k] * std::log(density);
  }
  
  return logLh;
}


//' Loglikelihood of phase-type using Pade approximation
//' 
//' Loglikelihood for a sample.
//' 
//' @param h Nuisance parameter.
//' @param alpha Initial probabilities.
//' @param S Sub-intensity matrix.
//' @param obs The observations.
//' @param weight The weights of the observations.
//' @param rcens Censored observations.
//' @param rcweight The weights of the censored observations.
//' 
// [[Rcpp::export]]
double logLikelihoodPH_PADE(double h, arma::vec & alpha, arma::mat & S, const Rcpp::NumericVector & obs, const Rcpp::NumericVector & weight, const Rcpp::NumericVector & rcens, const Rcpp::NumericVector & rcweight) {
  unsigned p{S.n_rows};
  arma::mat e;
  e.ones(S.n_cols, 1);
  arma::mat exit_vect = (S * (-1)) * e;
  
  arma::mat expm(size(S));
  
  double s_norm{inf_norm(S)};
  
  std::vector<arma::mat> aux_vect;
  
  vector_of_matrices_2(aux_vect, S, 6);
  
  arma::mat X(p,p);
  arma::mat D(p,p);
  
  const int q{6};
  int s{};
  double xmod{};
  double c{};
  
  arma::mat aux_mat(1,1);
  
  double density{0.0};
  
  double logLh{0.0};
  
  // Non censored data
  for (int k{0}; k < obs.size(); ++k) {
    // Matrix exponential
    int pind{1};
    int ee{static_cast<int>(log2(s_norm  * obs[k])) + 1};
    s = std::max(0, ee + 1);
    xmod = obs[k] / pow(2.0, s);
    c = 0.5;
    X = aux_vect[1] * (c * xmod);
    
    expm = aux_vect[0] + X;
    D = aux_vect[0] - X;
    
    for (int l{2}; l <= q; ++l) {
      c = c * static_cast<double>(q - l + 1) / static_cast<double>(l * (2 * q - l + 1));
      X = aux_vect[l] * (c * pow(xmod,l));
      expm = expm + X;
      if (pind) {
        D =  D + X;
      }
      else {
        D = D - X;
      }
      pind = !pind;
    }
    expm = inv(D) * expm;
    for (int l = 1; l <= s; ++l) {
      expm = expm * expm;
    }
    
    aux_mat = alpha.t() * expm * exit_vect;
    density = aux_mat(0,0);
    logLh += weight[k] * std::log(density);
  }
  //Right censored data
  for (int k{0}; k < rcens.size(); ++k) {
    // Matrix exponential
    int pind{1};
    int ee{static_cast<int>(log2(s_norm  * rcens[k])) + 1};
    s = std::max(0, ee + 1);
    xmod = rcens[k] / pow(2.0, s);
    c = 0.5;
    X = aux_vect[1] * (c * xmod);
    
    expm = aux_vect[0] + X;
    D = aux_vect[0] - X;
    
    for (int l{2}; l <= q; ++l) {
      c = c * static_cast<double>(q - l + 1) / static_cast<double>(l * (2 * q - l + 1));
      X = aux_vect[l] * (c * pow(xmod,l));
      expm = expm + X;
      if (pind) {
        D =  D + X;
      }
      else {
        D = D - X;
      }
      pind = !pind;
    }
    expm = inv(D) * expm;
    for (int l = 1; l <= s; ++l) {
      expm = expm * expm;
    }
    
    aux_mat = alpha.t() * expm * e;
    density = aux_mat(0,0);
    logLh += rcweight[k] * std::log(density);
  }
  return logLh;
}


//' Loglikelihood of matrix-Weibull using Pade
//' 
//' Loglikelihood for a sample.
//' 
//' @param h Nuisance parameter.
//' @param alpha Initial probabilities.
//' @param S Sub-intensity matrix.
//' @param beta Inhomogeneity parameter.
//' @param obs The observations.
//' @param weight The weights of the observations.
//' @param rcens Censored observations.
//' @param rcweight The weights of the censored observations.
//' 
// [[Rcpp::export]]
double logLikelihoodMweibull_PADE(double h, arma::vec & alpha, arma::mat & S, double beta, const Rcpp::NumericVector & obs, const Rcpp::NumericVector & weight, const Rcpp::NumericVector & rcens, const Rcpp::NumericVector & rcweight) {
  if(beta < 0) return NA_REAL;
  
  unsigned p{S.n_rows};
  arma::mat e; 
  e.ones(S.n_cols, 1);
  arma::mat exit_vect = (S * (-1)) * e;
  
  arma::mat expm(size(S));
  
  double s_norm{inf_norm(S)};
  
  std::vector<arma::mat> aux_vect;
  
  vector_of_matrices_2(aux_vect, S, 6);
  
  arma::mat X(p,p);
  arma::mat D(p,p);
  
  const int q{6};
  int s{};
  double xmod{};
  double c{};
  
  arma::mat aux_mat(1,1);
  
  double density{0.0};
  
  double logLh{0.0};
  
  // Non censored data
  for (int k{0}; k < obs.size(); ++k) {
    // Matrix exponential
    int pind{1};
    int ee{static_cast<int>(log2(s_norm  * pow(obs[k], beta))) + 1};
    s = std::max(0, ee + 1);
    xmod = pow(obs[k], beta) / pow(2.0, s);
    c = 0.5;
    X = aux_vect[1] * (c * xmod);
    
    expm = aux_vect[0] + X;
    D = aux_vect[0] - X;
    
    for (int l{2}; l <= q; ++l) {
      c = c * static_cast<double>(q - l + 1) / static_cast<double>(l * (2 * q - l + 1));
      X = aux_vect[l] * (c * pow(xmod,l));
      expm = expm + X;
      if (pind) {
        D =  D + X;
      }
      else {
        D = D - X;
      }
      pind = !pind;
    }
    expm = inv(D) * expm;
    for (int l = 1; l <= s; ++l) {
      expm = expm * expm;
    }
    
    aux_mat = alpha.t() * expm * exit_vect;
    density = aux_mat(0,0);
    logLh += weight[k] * (std::log(density) + std::log(beta) + (beta - 1) * std::log(obs[k]));
  }
  //Right censored data
  for (int k{0}; k < rcens.size(); ++k) {
    // Matrix exponential
    int pind{1};
    int ee{static_cast<int>(log2(s_norm  * pow(rcens[k], beta))) + 1};
    s = std::max(0, ee + 1);
    xmod = pow(rcens[k], beta) / pow(2.0, s);
    c = 0.5;
    X = aux_vect[1] * (c * xmod);
    
    expm = aux_vect[0] + X;
    D = aux_vect[0] - X;
    
    for (int l{2}; l <= q; ++l) {
      c = c * static_cast<double>(q - l + 1) / static_cast<double>(l * (2 * q - l + 1));
      X = aux_vect[l] * (c * pow(xmod,l));
      expm = expm + X;
      if (pind) {
        D =  D + X;
      }
      else {
        D = D - X;
      }
      pind = !pind;
    }
    expm = inv(D) * expm;
    for (int l = 1; l <= s; ++l) {
      expm = expm * expm;
    }
    
    aux_mat = alpha.t() * expm * e;
    density = aux_mat(0,0);
    logLh += rcweight[k] * std::log(density);
  }
  return logLh;
}


//' Loglikelihood of matrix-Pareto using Pade
//' 
//' Loglikelihood for a sample.
//' 
//' @param h Nuisance parameter.
//' @param alpha Initial probabilities.
//' @param S Sub-intensity matrix.
//' @param beta Inhomogeneity parameter.
//' @param obs The observations.
//' @param weight The weights of the observations.
//' @param rcens Censored observations.
//' @param rcweight The weights of the censored observations.
//' 
// [[Rcpp::export]]
double logLikelihoodMpareto_PADE(double h, arma::vec & alpha, arma::mat & S, double beta, const Rcpp::NumericVector & obs, const Rcpp::NumericVector & weight, const Rcpp::NumericVector & rcens, const Rcpp::NumericVector & rcweight) {
  if(beta < 0) return NA_REAL;
  
  unsigned p{S.n_rows};
  arma::mat e;
  e.ones(S.n_cols, 1);
  arma::mat exit_vect = (S * (-1)) * e;
  
  arma::mat expm(size(S));
  
  double s_norm{inf_norm(S)};
  
  std::vector<arma::mat> aux_vect;
  
  vector_of_matrices_2(aux_vect, S, 6);
  
  arma::mat X(p,p);
  arma::mat D(p,p);
  
  const int q{6};
  int s{};
  double xmod{};
  double c{};
  
  arma::mat aux_mat(1,1);
  
  double density{0.0};
  
  double logLh{0.0};
  
  // Non censored data
  for (int k{0}; k < obs.size(); ++k) {
    // Matrix exponential
    int pind{1};
    int ee{static_cast<int>(log2(s_norm  * std::log(obs[k] / beta + 1))) + 1};
    s = std::max(0, ee + 1);
    xmod = std::log(obs[k] / beta + 1) / pow(2.0, s);
    c = 0.5;
    X = aux_vect[1] * (c * xmod);
    
    expm = aux_vect[0] + X;
    D = aux_vect[0] - X;
    
    for (int l{2}; l <= q; ++l) {
      c = c * static_cast<double>(q - l + 1) / static_cast<double>(l * (2 * q - l + 1));
      X = aux_vect[l] * (c * pow(xmod,l));
      expm = expm + X;
      if (pind) {
        D =  D + X;
      }
      else {
        D = D - X;
      }
      pind = !pind;
    }
    expm = inv(D) * expm;
    for (int l = 1; l <= s; ++l) {
      expm = expm * expm;
    }
    
    aux_mat = alpha.t() * expm * exit_vect;
    density = aux_mat(0,0);
    logLh += weight[k] * (std::log(density) - std::log(obs[k] + beta));
  }
  //Right censored data
  for (int k{0}; k < rcens.size(); ++k) {
    // Matrix exponential
    int pind{1};
    int ee{static_cast<int>(log2(s_norm  * std::log(rcens[k] / beta + 1))) + 1};
    s = std::max(0, ee + 1);
    xmod = std::log(rcens[k] / beta + 1) / pow(2.0, s);
    c = 0.5;
    X = aux_vect[1] * (c * xmod);
    
    expm = aux_vect[0] + X;
    D = aux_vect[0] - X;
    
    for (int l{2}; l <= q; ++l) {
      c = c * static_cast<double>(q - l + 1) / static_cast<double>(l * (2 * q - l + 1));
      X = aux_vect[l] * (c * pow(xmod,l));
      expm = expm + X;
      if (pind) {
        D =  D + X;
      }
      else {
        D = D - X;
      }
      pind = !pind;
    }
    expm = inv(D) * expm;
    for (int l = 1; l <= s; ++l) {
      expm = expm * expm;
    }
    
    aux_mat = alpha.t() * expm * e;
    density = aux_mat(0,0);
    logLh += rcweight[k] * std::log(density);
  }
  
  return logLh;
}


//' Loglikelihood of matrix-lognormal using Pade
//' 
//' Loglikelihood for a sample.
//' 
//' @param h Nuisance parameter.
//' @param alpha Initial probabilities.
//' @param S Sub-intensity matrix.
//' @param beta Inhomogeneity parameter.
//' @param obs The observations.
//' @param weight The weights of the observations.
//' @param rcens Censored observations.
//' @param rcweight The weights of the censored observations.
//' 
// [[Rcpp::export]]
double logLikelihoodMlognormal_PADE(double h, arma::vec & alpha, arma::mat & S, double beta, const Rcpp::NumericVector & obs, const Rcpp::NumericVector & weight, const Rcpp::NumericVector & rcens, const Rcpp::NumericVector & rcweight) {
  if(beta < 0) return NA_REAL;
  
  unsigned p{S.n_rows};
  arma::mat e;
  e.ones(S.n_cols, 1);
  arma::mat exit_vect = (S * (-1)) * e;
  
  arma::mat expm(size(S));
  
  double s_norm{inf_norm(S)};
  
  std::vector<arma::mat> aux_vect;
  
  vector_of_matrices_2(aux_vect, S, 6);
  
  arma::mat X(p,p);
  arma::mat D(p,p);
  
  const int q{6};
  int s{};
  double xmod{};
  double c{};
  
  arma::mat aux_mat(1,1);
  
  double density{0.0};
  
  double logLh{0.0};
  
  // Non censored data
  for (int k{0}; k < obs.size(); ++k) {
    // Matrix exponential
    int pind{1};
    int ee{static_cast<int>(log2(s_norm  * pow(std::log(obs[k] + 1), beta))) + 1};
    s = std::max(0, ee + 1);
    xmod = pow(std::log(obs[k] + 1), beta) / pow(2.0, s);
    c = 0.5;
    X = aux_vect[1] * (c * xmod);
    
    expm = aux_vect[0] + X;
    D = aux_vect[0] - X;
    
    for (int l{2}; l <= q; ++l) {
      c = c * static_cast<double>(q - l + 1) / static_cast<double>(l * (2 * q - l + 1));
      X = aux_vect[l] * (c * pow(xmod,l));
      expm = expm + X;
      if (pind) {
        D =  D + X;
      }
      else {
        D = D - X;
      }
      pind = !pind;
    }
    expm = inv(D) * expm;
    for (int l = 1; l <= s; ++l) {
      expm = expm * expm;
    }
    
    aux_mat = alpha.t() * expm * exit_vect;
    density = aux_mat(0,0);
    logLh += weight[k] * (std::log(density) + std::log(beta) + (beta -1) * std::log(std::log(obs[k] + 1)) - std::log(obs[k] + 1));
  }
  //Right censored data
  for (int k{0}; k < rcens.size(); ++k) {
    // Matrix exponential
    int pind{1};
    int ee{static_cast<int>(log2(s_norm  * pow(std::log(rcens[k] + 1), beta))) + 1};
    s = std::max(0, ee + 1);
    xmod = pow(std::log(rcens[k] + 1), beta) / pow(2.0, s);
    c = 0.5;
    X = aux_vect[1] * (c * xmod);
    
    expm = aux_vect[0] + X;
    D = aux_vect[0] - X;
    
    for (int l{2}; l <= q; ++l) {
      c = c * static_cast<double>(q - l + 1) / static_cast<double>(l * (2 * q - l + 1));
      X = aux_vect[l] * (c * pow(xmod,l));
      expm = expm + X;
      if (pind) {
        D =  D + X;
      }
      else {
        D = D - X;
      }
      pind = !pind;
    }
    expm = inv(D) * expm;
    for (int l = 1; l <= s; ++l) {
      expm = expm * expm;
    }
    
    aux_mat = alpha.t() * expm * e;
    density = aux_mat(0,0);
    logLh += rcweight[k] * std::log(density);
  }
  
  return logLh;
}


//' Loglikelihood of matrix-loglogistic using Pade
//' 
//' Loglikelihood for a sample.
//' 
//' @param h Nuisance parameter.
//' @param alpha Initial probabilities.
//' @param S Sub-intensity matrix.
//' @param beta Inhomogeneity parameter.
//' @param obs The observations.
//' @param weight The weights of the observations.
//' @param rcens Censored observations.
//' @param rcweight The weights of the censored observations.
//' 
// [[Rcpp::export]]
double logLikelihoodMloglogistic_PADE(double h, arma::vec & alpha, arma::mat & S, Rcpp::NumericVector beta, const Rcpp::NumericVector & obs, const Rcpp::NumericVector & weight, const Rcpp::NumericVector & rcens, const Rcpp::NumericVector & rcweight) {
  if(beta[0] < 0 || beta[1] < 0) return NA_REAL;
  
  unsigned p{S.n_rows};
  arma::mat e;
  e.ones(S.n_cols, 1);
  arma::mat exit_vect = (S * (-1)) * e;
  
  arma::mat expm(size(S));
  
  double s_norm{inf_norm(S)};
  
  std::vector<arma::mat> aux_vect;
  
  vector_of_matrices_2(aux_vect, S, 6);
  
  arma::mat X(p,p);
  arma::mat D(p,p);
  
  const int q{6};
  int s{};
  double xmod{};
  double c{};
  
  arma::mat aux_mat(1,1);
  
  double density{0.0};
  
  double logLh{0.0};
  
  // Non censored data
  for (int k{0}; k < obs.size(); ++k) {
    // Matrix exponential
    int pind{1};
    int ee{static_cast<int>(log2(s_norm  * std::log(pow(obs[k] / beta[0], beta[1]) + 1))) + 1};
    s = std::max(0, ee + 1);
    xmod = std::log(pow(obs[k] / beta[0], beta[1]) + 1) / pow(2.0, s);
    c = 0.5;
    X = aux_vect[1] * (c * xmod);
    
    expm = aux_vect[0] + X;
    D = aux_vect[0] - X;
    
    for (int l{2}; l <= q; ++l) {
      c = c * static_cast<double>(q - l + 1) / static_cast<double>(l * (2 * q - l + 1));
      X = aux_vect[l] * (c * pow(xmod,l));
      expm = expm + X;
      if (pind) {
        D =  D + X;
      }
      else {
        D = D - X;
      }
      pind = !pind;
    }
    expm = inv(D) * expm;
    for (int l = 1; l <= s; ++l) {
      expm = expm * expm;
    }
    
    aux_mat = alpha.t() * expm * exit_vect;
    density = aux_mat(0,0);
    logLh += weight[k] * (std::log(density) + std::log(beta[1]) - std::log(beta[0]) + (beta[1] - 1) * (std::log(obs[k]) - std::log(beta[0])) - std::log(pow(obs[k] / beta[0], beta[1]) + 1));
  }
  //Right censored data
  for (int k{0}; k < rcens.size(); ++k) {
    // Matrix exponential
    int pind{1};
    int ee{static_cast<int>(log2(s_norm  * std::log(pow(rcens[k] / beta[0], beta[1]) + 1))) + 1};
    s = std::max(0, ee + 1);
    xmod = std::log(pow(rcens[k] / beta[0], beta[1]) + 1) / pow(2.0, s);
    c = 0.5;
    X = aux_vect[1] * (c * xmod);
    
    expm = aux_vect[0] + X;
    D = aux_vect[0] - X;
    
    for (int l{2}; l <= q; ++l) {
      c = c * static_cast<double>(q - l + 1) / static_cast<double>(l * (2 * q - l + 1));
      X = aux_vect[l] * (c * pow(xmod,l));
      expm = expm + X;
      if (pind) {
        D =  D + X;
      }
      else {
        D = D - X;
      }
      pind = !pind;
    }
    expm = inv(D) * expm;
    for (int l = 1; l <= s; ++l) {
      expm = expm * expm;
    }
    
    aux_mat = alpha.t() * expm * e;
    density = aux_mat(0,0);
    logLh += rcweight[k] * std::log(density);
  }
  
  return logLh;
}


//' Loglikelihood of matrix-Gompertz using Pade
//' 
//' Loglikelihood for a sample.
//' 
//' @param h Nuisance parameter.
//' @param alpha Initial probabilities.
//' @param S Sub-intensity matrix.
//' @param beta Inhomogeneity parameter.
//' @param obs The observations.
//' @param weight The weights of the observations.
//' @param rcens Censored observations.
//' @param rcweight The weights of the censored observations.
//' 
// [[Rcpp::export]]
double logLikelihoodMgompertz_PADE(double h, arma::vec & alpha, arma::mat & S, double beta, const Rcpp::NumericVector & obs, const Rcpp::NumericVector & weight, const Rcpp::NumericVector & rcens, const Rcpp::NumericVector & rcweight) {
  if(beta < 0) return NA_REAL;
  
  unsigned p{S.n_rows};
  arma::mat e;
  e.ones(S.n_cols, 1);
  arma::mat exit_vect = (S * (-1)) * e;
  
  arma::mat expm(size(S));
  
  double s_norm{inf_norm(S)};
  
  std::vector<arma::mat> aux_vect;
  
  vector_of_matrices_2(aux_vect, S, 6);
  
  arma::mat X(p,p);
  arma::mat D(p,p);
  
  const int q{6};
  int s{};
  double xmod{};
  double c{};
  
  arma::mat aux_mat(1,1);
  
  double density{0.0};
  
  double logLh{0.0};
  
  // Non censored data
  for (int k{0}; k < obs.size(); ++k) {
    // Matrix exponential
    int pind{1};
    int ee{static_cast<int>(log2(s_norm  * (exp(obs[k] * beta) - 1) / beta)) + 1};
    s = std::max(0, ee + 1);
    xmod = (exp(obs[k] * beta) - 1) / beta / pow(2.0, s);
    c = 0.5;
    X = aux_vect[1] * (c * xmod);
    
    expm = aux_vect[0] + X;
    D = aux_vect[0] - X;
    
    for (int l{2}; l <= q; ++l) {
      c = c * static_cast<double>(q - l + 1) / static_cast<double>(l * (2 * q - l + 1));
      X = aux_vect[l] * (c * pow(xmod,l));
      expm = expm + X;
      if (pind) {
        D =  D + X;
      }
      else {
        D = D - X;
      }
      pind = !pind;
    }
    expm = inv(D) * expm;
    for (int l = 1; l <= s; ++l) {
      expm = expm * expm;
    }
    
    aux_mat = alpha.t() * expm * exit_vect;
    density = aux_mat(0,0);
    logLh += weight[k] * (std::log(density) + obs[k] * beta);
  }
  //Right censored data
  for (int k{0}; k < rcens.size(); ++k) {
    // Matrix exponential
    int pind{1};
    int ee{static_cast<int>(log2(s_norm  * (exp(rcens[k] * beta) - 1) / beta)) + 1};
    s = std::max(0, ee + 1);
    xmod = (exp(rcens[k] * beta) - 1) / beta / pow(2.0, s);
    c = 0.5;
    X = aux_vect[1] * (c * xmod);
    
    expm = aux_vect[0] + X;
    D = aux_vect[0] - X;
    
    for (int l{2}; l <= q; ++l) {
      c = c * static_cast<double>(q - l + 1) / static_cast<double>(l * (2 * q - l + 1));
      X = aux_vect[l] * (c * pow(xmod,l));
      expm = expm + X;
      if (pind) {
        D =  D + X;
      }
      else {
        D = D - X;
      }
      pind = !pind;
    }
    expm = inv(D) * expm;
    for (int l = 1; l <= s; ++l) {
      expm = expm * expm;
    }
    
    aux_mat = alpha.t() * expm * e;
    density = aux_mat(0,0);
    logLh += rcweight[k] * std::log(density);
  }
  
  return logLh;
}


//' Loglikelihood of matrix-GEV using Pade
//' 
//' Loglikelihood for a sample 
//' 
//' @param h Nuisance parameter.
//' @param alpha Initial probabilities.
//' @param S sub-intensity matrix.
//' @param beta Inhomogeneity parameter.
//' @param obs The observations.
//' @param weight The weights of the observations.
//' @param rcens Censored observations.
//' @param rcweight The weights of the censored observations.
//' 
// [[Rcpp::export]]
double logLikelihoodMgev_PADE(double h, arma::vec & alpha, arma::mat & S, Rcpp::NumericVector beta, const Rcpp::NumericVector & obs, const Rcpp::NumericVector & weight, const Rcpp::NumericVector & rcens, const Rcpp::NumericVector & rcweight) {
  if(beta[1] < 0) return NA_REAL;
  
  unsigned p{S.n_rows};
  arma::mat e;
  e.ones(S.n_cols, 1);
  arma::mat exit_vect = (S * (-1)) * e;
  
  arma::mat expm(size(S));
  
  double s_norm{inf_norm(S)};
  
  std::vector<arma::mat> aux_vect;
  
  vector_of_matrices_2(aux_vect, S, 6);
  
  arma::mat X(p,p);
  arma::mat D(p,p);
  
  const int q{6};
  int s{};
  double xmod{};
  double c{};
  
  arma::mat aux_mat(1,1);
  
  double density{0.0};
  
  double logLh{0.0};
  if (beta[2] == 0) {
    // Non censored data
    for (int k{0}; k < obs.size(); ++k) {
      // Matrix exponential
      int pind{1};
      int ee{static_cast<int>(log2(s_norm  * exp(-(obs[k] - beta[0]) / beta[1]))) + 1};
      s = std::max(0, ee + 1);
      xmod = exp(-(obs[k] - beta[0]) / beta[1])/ pow(2.0, s);
      c = 0.5;
      X = aux_vect[1] * (c * xmod);
      
      expm = aux_vect[0] + X;
      D = aux_vect[0] - X;
      
      for (int l{2}; l <= q; ++l) {
        c = c * static_cast<double>(q - l + 1) / static_cast<double>(l * (2 * q - l + 1));
        X = aux_vect[l] * (c * pow(xmod,l));
        expm = expm + X;
        if (pind) {
          D =  D + X;
        }
        else {
          D = D - X;
        }
        pind = !pind;
      }
      expm = inv(D) * expm;
      for (int l = 1; l <= s; ++l) {
        expm = expm * expm;
      }
      
      aux_mat = alpha.t() * expm * exit_vect;
      density = aux_mat(0,0);
      logLh += weight[k] * (std::log(density) - std::log(beta[1]) - (obs[k] - beta[0]) / beta[1]);
    }
    //Right censored data
    for (int k{0}; k < rcens.size(); ++k) {
      // Matrix exponential
      int pind{1};
      int ee{static_cast<int>(log2(s_norm  * exp(-(rcens[k] - beta[0]) / beta[1]))) + 1};
      s = std::max(0, ee + 1);
      xmod = exp(-(rcens[k] - beta[0]) / beta[1]) / pow(2.0, s);
      c = 0.5;
      X = aux_vect[1] * (c * xmod);
      
      expm = aux_vect[0] + X;
      D = aux_vect[0] - X;
      
      for (int l{2}; l <= q; ++l) {
        c = c * static_cast<double>(q - l + 1) / static_cast<double>(l * (2 * q - l + 1));
        X = aux_vect[l] * (c * pow(xmod,l));
        expm = expm + X;
        if (pind) {
          D =  D + X;
        }
        else {
          D = D - X;
        }
        pind = !pind;
      }
      expm = inv(D) * expm;
      for (int l = 1; l <= s; ++l) {
        expm = expm * expm;
      }
      
      aux_mat = alpha.t() * expm * e;
      density = aux_mat(0,0);
      logLh += rcweight[k] * std::log(density);
    }
  }
  else{
    // Non censored data
    for (int k{0}; k < obs.size(); ++k) {
      // Matrix exponential
      int pind{1};
      int ee{static_cast<int>(log2(s_norm  * pow(1 + (beta[2] / beta[1]) * (obs[k] - beta[0]) , - 1 / beta[2]))) + 1};
      s = std::max(0, ee + 1);
      xmod = pow(1 + (beta[2] / beta[1]) * (obs[k] - beta[0]) , - 1 / beta[2]) / pow(2.0, s);
      c = 0.5;
      X = aux_vect[1] * (c * xmod);
      
      expm = aux_vect[0] + X;
      D = aux_vect[0] - X;
      
      for (int l{2}; l <= q; ++l) {
        c = c * static_cast<double>(q - l + 1) / static_cast<double>(l * (2 * q - l + 1));
        X = aux_vect[l] * (c * pow(xmod,l));
        expm = expm + X;
        if (pind) {
          D =  D + X;
        }
        else {
          D = D - X;
        }
        pind = !pind;
      }
      expm = inv(D) * expm;
      for (int l = 1; l <= s; ++l) {
        expm = expm * expm;
      }
      
      aux_mat = alpha.t() * expm * exit_vect;
      density = aux_mat(0,0);
      logLh += weight[k] * (std::log(density) - std::log(beta[1]) - (1 + 1 / beta[2]) * std::log(1 + (beta[2] / beta[1]) * (obs[k] - beta[0])));
    }
    //Right censored data
    for (int k{0}; k < rcens.size(); ++k) {
      // Matrix exponential
      int pind{1};
      int ee{static_cast<int>(log2(s_norm  * pow(1 + (beta[2] / beta[1]) * (rcens[k] - beta[0]) , - 1 / beta[2]))) + 1};
      s = std::max(0, ee + 1);
      xmod = pow(1 + (beta[2] / beta[1]) * (rcens[k] - beta[0]) , - 1 / beta[2]) / pow(2.0, s);
      c = 0.5;
      X = aux_vect[1] * (c * xmod);
      
      expm = aux_vect[0] + X;
      D = aux_vect[0] - X;
      
      for (int l{2}; l <= q; ++l) {
        c = c * static_cast<double>(q - l + 1) / static_cast<double>(l * (2 * q - l + 1));
        X = aux_vect[l] * (c * pow(xmod,l));
        expm = expm + X;
        if (pind) {
          D =  D + X;
        }
        else {
          D = D - X;
        }
        pind = !pind;
      }
      expm = inv(D) * expm;
      for (int l = 1; l <= s; ++l) {
        expm = expm * expm;
      }
      
      aux_mat = alpha.t() * expm * e;
      density = aux_mat(0,0);
      logLh += rcweight[k] * std::log(density);
    }
  }
  return logLh;
}


////////////////////////////////////////////
// Scaled versions of loglikelihoods (for regression):
////////////////////////////////////////////

//' Loglikelihood of PI with phase-type using Pade
//' 
//' Loglikelihood for a sample.
//' 
//' @param h Nuisance parameter.
//' @param alpha Initial probabilities.
//' @param S Sub-intensity matrix.
//' @param obs The observations.
//' @param weight The weights of the observations.
//' @param rcens Censored observations.
//' @param rcweight The weights of the censored observations.
//' @param scale1 Scale for observations.
//' @param scale2 Scale for censored observations.
//' 
// [[Rcpp::export]]
double logLikelihoodPH_PADEs(double h, arma::vec & alpha, arma::mat & S, const Rcpp::NumericVector & obs, const Rcpp::NumericVector & weight, const Rcpp::NumericVector & rcens, const Rcpp::NumericVector & rcweight, const Rcpp::NumericVector & scale1, const Rcpp::NumericVector & scale2) {
  unsigned p{S.n_rows};
  arma::mat e;
  e.ones(S.n_cols, 1);
  arma::mat exit_vect = (S * (-1)) * e;
  
  arma::mat expm(size(S));
  
  double s_norm{inf_norm(S)};
  
  std::vector<arma::mat> aux_vect;
  
  vector_of_matrices_2(aux_vect, S, 6);
  
  arma::mat X(p,p);
  arma::mat D(p,p);
  
  const int q{6};
  int s{};
  double xmod{};
  double c{};
  
  arma::mat aux_mat(1,1);
  
  double density{0.0};
  
  double logLh{0.0};
  
  // Non censored data
  for (int k{0}; k < obs.size(); ++k) {
    // Matrix exponential
    int pind{1};
    int ee{static_cast<int>(log2(s_norm  * scale1[k] * obs[k])) + 1};
    s = std::max(0, ee + 1);
    xmod = scale1[k] * obs[k] / pow(2.0, s);
    c = 0.5;
    X = aux_vect[1] * (c * xmod);
    
    expm = aux_vect[0] + X;
    D = aux_vect[0] - X;
    
    for (int l{2}; l <= q; ++l) {
      c = c * static_cast<double>(q - l + 1) / static_cast<double>(l * (2 * q - l + 1));
      X = aux_vect[l] * (c * pow(xmod,l));
      expm = expm + X;
      if (pind) {
        D =  D + X;
      }
      else {
        D = D - X;
      }
      pind = !pind;
    }
    expm = inv(D) * expm;
    for (int l = 1; l <= s; ++l) {
      expm = expm * expm;
    }
    
    aux_mat = alpha.t() * expm * exit_vect;
    density = aux_mat(0,0);
    logLh += weight[k] * (std::log(density) + std::log(scale1[k]));
  }
  //Right censored data
  for (int k{0}; k < rcens.size(); ++k) {
    // Matrix exponential
    int pind{1};
    int ee{static_cast<int>(log2(s_norm  * scale2[k] * rcens[k])) + 1};
    s = std::max(0, ee + 1);
    xmod = scale2[k] * rcens[k] / pow(2.0, s);
    c = 0.5;
    X = aux_vect[1] * (c * xmod);
    
    expm = aux_vect[0] + X;
    D = aux_vect[0] - X;
    
    for (int l{2}; l <= q; ++l) {
      c = c * static_cast<double>(q - l + 1) / static_cast<double>(l * (2 * q - l + 1));
      X = aux_vect[l] * (c * pow(xmod,l));
      expm = expm + X;
      if (pind) {
        D =  D + X;
      }
      else {
        D = D - X;
      }
      pind = !pind;
    }
    expm = inv(D) * expm;
    for (int l = 1; l <= s; ++l) {
      expm = expm * expm;
    }
    
    aux_mat = alpha.t() * expm * e;
    density = aux_mat(0,0);
    logLh += rcweight[k] * std::log(density);
  }
  
  return logLh;
}


//' Loglikelihood of PI with matrix-Weibull using Pade
//' 
//' Loglikelihood for a sample.
//' 
//' @param h Nuisance parameter.
//' @param alpha Initial probabilities.
//' @param S Sub-intensity matrix.
//' @param beta Inhomogeneity parameter.
//' @param obs The observations.
//' @param weight The weights of the observations.
//' @param rcens Censored observations.
//' @param rcweight The weights of the censored observations.
//' @param scale1 Scale for observations.
//' @param scale2 Scale for censored observations.
//' 
// [[Rcpp::export]]
double logLikelihoodMweibull_PADEs(double h, arma::vec & alpha, arma::mat & S, double beta, const Rcpp::NumericVector & obs, const Rcpp::NumericVector & weight, const Rcpp::NumericVector & rcens, const Rcpp::NumericVector & rcweight, const Rcpp::NumericVector & scale1, const Rcpp::NumericVector & scale2) {
  if(beta < 0) return NA_REAL;
  
  unsigned p{S.n_rows};
  arma::mat e;
  e.ones(S.n_cols, 1);
  arma::mat exit_vect = (S * (-1)) * e;
  
  arma::mat expm(size(S));
  
  double s_norm{inf_norm(S)};
  
  std::vector<arma::mat> aux_vect;
  
  vector_of_matrices_2(aux_vect, S, 6);
  
  arma::mat X(p,p);
  arma::mat D(p,p);
  
  const int q{6};
  int s{};
  double xmod{};
  double c{};
  
  arma::mat aux_mat(1,1);
  
  double density{0.0};
  
  double logLh{0.0};
  
  // Non censored data
  for (int k{0}; k < obs.size(); ++k) {
    // Matrix exponential
    int pind{1};
    int ee{static_cast<int>(log2(s_norm  * scale1[k] * pow(obs[k], beta))) + 1};
    s = std::max(0, ee + 1);
    xmod = scale1[k] * pow(obs[k], beta) / pow(2.0, s);
    c = 0.5;
    X = aux_vect[1] * (c * xmod);
    
    expm = aux_vect[0] + X;
    D = aux_vect[0] - X;
    
    for (int l{2}; l <= q; ++l) {
      c = c * static_cast<double>(q - l + 1) / static_cast<double>(l * (2 * q - l + 1));
      X = aux_vect[l] * (c * pow(xmod,l));
      expm = expm + X;
      if (pind) {
        D =  D + X;
      }
      else {
        D = D - X;
      }
      pind = !pind;
    }
    expm = inv(D) * expm;
    for (int l = 1; l <= s; ++l) {
      expm = expm * expm;
    }
    
    aux_mat = alpha.t() * expm * exit_vect;
    density = aux_mat(0,0);
    logLh += weight[k] * (std::log(density) + std::log(scale1[k]) + std::log(beta) + (beta - 1) * std::log(obs[k]));
  }
  //Right censored data
  for (int k{0}; k < rcens.size(); ++k) {
    // Matrix exponential
    int pind{1};
    int ee{static_cast<int>(log2(s_norm  * scale2[k] * pow(rcens[k], beta))) + 1};
    s = std::max(0, ee + 1);
    xmod = scale2[k] * pow(rcens[k], beta) / pow(2.0, s);
    c = 0.5;
    X = aux_vect[1] * (c * xmod);
    
    expm = aux_vect[0] + X;
    D = aux_vect[0] - X;
    
    for (int l{2}; l <= q; ++l) {
      c = c * static_cast<double>(q - l + 1) / static_cast<double>(l * (2 * q - l + 1));
      X = aux_vect[l] * (c * pow(xmod,l));
      expm = expm + X;
      if (pind) {
        D =  D + X;
      }
      else {
        D = D - X;
      }
      pind = !pind;
    }
    expm = inv(D) * expm;
    for (int l = 1; l <= s; ++l) {
      expm = expm * expm;
    }
    
    aux_mat = alpha.t() * expm * e;
    density = aux_mat(0,0);
    logLh += rcweight[k] * std::log(density);
  }
  
  return logLh;
}


//' Loglikelihood of PI with matrix-Pareto using Pade
//' 
//' Loglikelihood for a sample.
//' 
//' @param h Nuisance parameter.
//' @param alpha Initial probabilities.
//' @param S Sub-intensity matrix.
//' @param beta Inhomogeneity parameter.
//' @param obs The observations.
//' @param weight Weights of the observations.
//' @param rcens Censored observations.
//' @param rcweight Weights of the censored observations.
//' @param scale1 Scale for observations.
//' @param scale2 Scale for censored observations.
//' 
// [[Rcpp::export]]
double logLikelihoodMpareto_PADEs(double h, arma::vec & alpha, arma::mat & S, double beta, const Rcpp::NumericVector & obs, const Rcpp::NumericVector & weight, const Rcpp::NumericVector & rcens, const Rcpp::NumericVector & rcweight, const Rcpp::NumericVector & scale1, const Rcpp::NumericVector & scale2) {
  if(beta < 0) return NA_REAL;
  
  unsigned p{S.n_rows};
  arma::mat e;
  e.ones(S.n_cols, 1);
  arma::mat exit_vect = (S * (-1)) * e;
  
  arma::mat expm(size(S));
  
  double s_norm{inf_norm(S)};
  
  std::vector<arma::mat> aux_vect;
  
  vector_of_matrices_2(aux_vect, S, 6);
  
  arma::mat X(p,p);
  arma::mat D(p,p);
  
  const int q{6};
  int s{};
  double xmod{};
  double c{};
  
  arma::mat aux_mat(1,1);
  
  double density{0.0};
  
  double logLh{0.0};
  
  // Non censored data
  for (int k{0}; k < obs.size(); ++k) {
    // Matrix exponential
    int pind{1};
    int ee{static_cast<int>(log2(s_norm  * scale1[k] * std::log(obs[k] / beta + 1))) + 1};
    s = std::max(0, ee + 1);
    xmod = scale1[k] * std::log(obs[k] / beta + 1) / pow(2.0, s);
    c = 0.5;
    X = aux_vect[1] * (c * xmod);
    
    expm = aux_vect[0] + X;
    D = aux_vect[0] - X;
    
    for (int l{2}; l <= q; ++l) {
      c = c * static_cast<double>(q - l + 1) / static_cast<double>(l * (2 * q - l + 1));
      X = aux_vect[l] * (c * pow(xmod,l));
      expm = expm + X;
      if (pind) {
        D =  D + X;
      }
      else {
        D = D - X;
      }
      pind = !pind;
    }
    expm = inv(D) * expm;
    for (int l = 1; l <= s; ++l) {
      expm = expm * expm;
    }
    
    aux_mat = alpha.t() * expm * exit_vect;
    density = aux_mat(0,0);
    logLh += weight[k] * (std::log(density) + std::log(scale1[k]) - std::log(obs[k] + beta));
  }
  //Right censored data
  for (int k{0}; k < rcens.size(); ++k) {
    // Matrix exponential
    int pind{1};
    int ee{static_cast<int>(log2(s_norm  * scale2[k] * std::log(rcens[k] / beta + 1))) + 1};
    s = std::max(0, ee + 1);
    xmod = scale2[k] * std::log(rcens[k] / beta + 1) / pow(2.0, s);
    c = 0.5;
    X = aux_vect[1] * (c * xmod);
    
    expm = aux_vect[0] + X;
    D = aux_vect[0] - X;
    
    for (int l{2}; l <= q; ++l) {
      c = c * static_cast<double>(q - l + 1) / static_cast<double>(l * (2 * q - l + 1));
      X = aux_vect[l] * (c * pow(xmod,l));
      expm = expm + X;
      if (pind) {
        D =  D + X;
      }
      else {
        D = D - X;
      }
      pind = !pind;
    }
    expm = inv(D) * expm;
    for (int l = 1; l <= s; ++l) {
      expm = expm * expm;
    }
    
    aux_mat = alpha.t() * expm * e;
    density = aux_mat(0,0);
    logLh += rcweight[k] * std::log(density);
  }
  
  return logLh;
}


//' Loglikelihood of PI with matrix-lognormal using Pade
//' 
//' Loglikelihood for a sample.
//' 
//' @param h Nuisance parameter.
//' @param alpha Initial probabilities.
//' @param S Sub-intensity matrix.
//' @param beta Inhomogeneity parameter.
//' @param obs The observations.
//' @param weight Weights of the observations.
//' @param rcens Censored observations.
//' @param rcweight Weights of the censored observations.
//' @param scale1 Scale for observations.
//' @param scale2 Scale for censored observations.
//' 
// [[Rcpp::export]]
double logLikelihoodMlognormal_PADEs(double h, arma::vec & alpha, arma::mat & S, double beta, const Rcpp::NumericVector & obs, const Rcpp::NumericVector & weight, const Rcpp::NumericVector & rcens, const Rcpp::NumericVector & rcweight, const Rcpp::NumericVector & scale1, const Rcpp::NumericVector & scale2) {
  if(beta < 0) return NA_REAL;
  
  unsigned p{S.n_rows};
  arma::mat e;
  e.ones(S.n_cols, 1);
  arma::mat exit_vect = (S * (-1)) * e;
  
  arma::mat expm(size(S));
  
  double s_norm{inf_norm(S)};
  
  std::vector<arma::mat> aux_vect;
  
  vector_of_matrices_2(aux_vect, S, 6);
  
  arma::mat X(p,p);
  arma::mat D(p,p);
  
  const int q{6};
  int s{};
  double xmod{};
  double c{};
  
  arma::mat aux_mat(1,1);
  
  double density{0.0};
  
  double logLh{0.0};
  
  // Non censored data
  for (int k{0}; k < obs.size(); ++k) {
    // Matrix exponential
    int pind{1};
    int ee{static_cast<int>(log2(s_norm  * scale1[k] * pow(std::log(obs[k] + 1), beta))) + 1};
    s = std::max(0, ee + 1);
    xmod = scale1[k] * pow(std::log(obs[k] + 1), beta) / pow(2.0, s);
    c = 0.5;
    X = aux_vect[1] * (c * xmod);
    
    expm = aux_vect[0] + X;
    D = aux_vect[0] - X;
    
    for (int l{2}; l <= q; ++l) {
      c = c * static_cast<double>(q - l + 1) / static_cast<double>(l * (2 * q - l + 1));
      X = aux_vect[l] * (c * pow(xmod,l));
      expm = expm + X;
      if (pind) {
        D =  D + X;
      }
      else {
        D = D - X;
      }
      pind = !pind;
    }
    expm = inv(D) * expm;
    for (int l = 1; l <= s; ++l) {
      expm = expm * expm;
    }
    
    aux_mat = alpha.t() * expm * exit_vect;
    density = aux_mat(0,0);
    logLh += weight[k] * (std::log(density) + std::log(scale1[k]) + std::log(beta) + (beta -1) * std::log(std::log(obs[k] + 1)) - std::log(obs[k] + 1));
  }
  //Right censored data
  for (int k{0}; k < rcens.size(); ++k) {
    // Matrix exponential
    int pind{1};
    int ee{static_cast<int>(log2(s_norm  * scale2[k] * pow(std::log(rcens[k] + 1), beta))) + 1};
    s = std::max(0, ee + 1);
    xmod = scale2[k] * pow(std::log(rcens[k] + 1), beta) / pow(2.0, s);
    c = 0.5;
    X = aux_vect[1] * (c * xmod);
    
    expm = aux_vect[0] + X;
    D = aux_vect[0] - X;
    
    for (int l{2}; l <= q; ++l) {
      c = c * static_cast<double>(q - l + 1) / static_cast<double>(l * (2 * q - l + 1));
      X = aux_vect[l] * (c * pow(xmod,l));
      expm = expm + X;
      if (pind) {
        D =  D + X;
      }
      else {
        D = D - X;
      }
      pind = !pind;
    }
    expm = inv(D) * expm;
    for (int l = 1; l <= s; ++l) {
      expm = expm * expm;
    }
    
    aux_mat = alpha.t() * expm * e;
    density = aux_mat(0,0);
    logLh += rcweight[k] * std::log(density);
  }
  
  return logLh;
}


//' Loglikelihood of PI with matrix-loglogistic using Pade
//' 
//' Loglikelihood for a sample.
//' 
//' @param h Nuisance parameter.
//' @param alpha Initial probabilities.
//' @param S Sub-intensity matrix.
//' @param beta Inhomogeneity parameter.
//' @param obs The observations.
//' @param weight Weights of the observations.
//' @param rcens Censored observations.
//' @param rcweight Weights of the censored observations.
//' @param scale1 Scale for observations.
//' @param scale2 Scale for censored observations.
//' 
// [[Rcpp::export]]
double logLikelihoodMloglogistic_PADEs(double h, arma::vec & alpha, arma::mat & S, Rcpp::NumericVector beta, const Rcpp::NumericVector & obs, const Rcpp::NumericVector & weight, const Rcpp::NumericVector & rcens, const Rcpp::NumericVector & rcweight, const Rcpp::NumericVector & scale1, const Rcpp::NumericVector & scale2) {
  if(beta[0] < 0 || beta[1] < 0) return NA_REAL;
  
  unsigned p{S.n_rows};
  arma::mat e;
  e.ones(S.n_cols, 1);
  arma::mat exit_vect = (S * (-1)) * e;
  
  arma::mat expm(size(S));
  
  double s_norm{inf_norm(S)};
  
  std::vector<arma::mat> aux_vect;
  
  vector_of_matrices_2(aux_vect, S, 6);
  
  arma::mat X(p,p);
  arma::mat D(p,p);
  
  const int q{6};
  int s{};
  double xmod{};
  double c{};
  
  arma::mat aux_mat(1,1);
  
  double density{0.0};
  
  double logLh{0.0};
  
  // Non censored data
  for (int k{0}; k < obs.size(); ++k) {
    // Matrix exponential
    int pind{1};
    int ee{static_cast<int>(log2(s_norm  * scale1[k] * std::log(pow(obs[k] / beta[0], beta[1]) + 1))) + 1};
    s = std::max(0, ee + 1);
    xmod = scale1[k] * std::log(pow(obs[k] / beta[0], beta[1]) + 1) / pow(2.0, s);
    c = 0.5;
    X = aux_vect[1] * (c * xmod);
    
    expm = aux_vect[0] + X;
    D = aux_vect[0] - X;
    
    for (int l{2}; l <= q; ++l) {
      c = c * static_cast<double>(q - l + 1) / static_cast<double>(l * (2 * q - l + 1));
      X = aux_vect[l] * (c * pow(xmod,l));
      expm = expm + X;
      if (pind) {
        D =  D + X;
      }
      else {
        D = D - X;
      }
      pind = !pind;
    }
    expm = inv(D) * expm;
    for (int l = 1; l <= s; ++l) {
      expm = expm * expm;
    }
    
    aux_mat = alpha.t() * expm * exit_vect;
    density = aux_mat(0,0);
    logLh += weight[k] * (std::log(density) + std::log(scale1[k]) + std::log(beta[1]) - std::log(beta[0]) + (beta[1] - 1) * (std::log(obs[k]) - std::log(beta[0])) - std::log(pow(obs[k] / beta[0], beta[1]) + 1));
  }
  //Right censored data
  for (int k{0}; k < rcens.size(); ++k) {
    // Matrix exponential
    int pind{1};
    int ee{static_cast<int>(log2(s_norm  * scale2[k] * std::log(pow(rcens[k] / beta[0], beta[1]) + 1))) + 1};
    s = std::max(0, ee + 1);
    xmod = scale2[k] * std::log(pow(rcens[k] / beta[0], beta[1]) + 1) / pow(2.0, s);
    c = 0.5;
    X = aux_vect[1] * (c * xmod);
    
    expm = aux_vect[0] + X;
    D = aux_vect[0] - X;
    
    for (int l{2}; l <= q; ++l) {
      c = c * static_cast<double>(q - l + 1) / static_cast<double>(l * (2 * q - l + 1));
      X = aux_vect[l] * (c * pow(xmod,l));
      expm = expm + X;
      if (pind) {
        D =  D + X;
      }
      else {
        D = D - X;
      }
      pind = !pind;
    }
    expm = inv(D) * expm;
    for (int l = 1; l <= s; ++l) {
      expm = expm * expm;
    }
    
    aux_mat = alpha.t() * expm * e;
    density = aux_mat(0,0);
    logLh += rcweight[k] * std::log(density);
  }
  
  return logLh;
}


//' Loglikelihood of PI with matrix-Gompertz using Pade
//' 
//' Loglikelihood for a sample.
//' 
//' @param h Nuisance parameter.
//' @param alpha Initial probabilities.
//' @param S Sub-intensity.
//' @param beta Inhomogeneity parameter.
//' @param obs The observations.
//' @param weight Weights of the observations.
//' @param rcens Censored observations.
//' @param rcweight Weights of the censored observations.
//' @param scale1 Scale for observations.
//' @param scale2 Scale for censored observations.
//' 
// [[Rcpp::export]]
double logLikelihoodMgompertz_PADEs(double h, arma::vec & alpha, arma::mat & S, double beta, const Rcpp::NumericVector & obs, const Rcpp::NumericVector & weight, const Rcpp::NumericVector & rcens, const Rcpp::NumericVector & rcweight, const Rcpp::NumericVector & scale1, const Rcpp::NumericVector & scale2) {
  if(beta < 0) return NA_REAL;
  
  unsigned p{S.n_rows};
  arma::mat e;
  e.ones(S.n_cols, 1);
  arma::mat exit_vect = (S * (-1)) * e;
  
  arma::mat expm(size(S));
  
  double s_norm{inf_norm(S)};
  
  std::vector<arma::mat> aux_vect;
  
  vector_of_matrices_2(aux_vect, S, 6);
  
  arma::mat X(p,p);
  arma::mat D(p,p);
  
  const int q{6};
  int s{};
  double xmod{};
  double c{};
  
  arma::mat aux_mat(1,1);
  
  double density{0.0};
  
  double logLh{0.0};
  
  // Non censored data
  for (int k{0}; k < obs.size(); ++k) {
    // Matrix exponential
    int pind{1};
    int ee{static_cast<int>(log2(s_norm  * scale1[k] * (exp(obs[k] * beta) - 1) / beta)) + 1};
    s = std::max(0, ee + 1);
    xmod = scale1[k] * (exp(obs[k] * beta) - 1) / beta / pow(2.0, s);
    c = 0.5;
    X = aux_vect[1] * (c * xmod);
    
    expm = aux_vect[0] + X;
    D = aux_vect[0] - X;
    
    for (int l{2}; l <= q; ++l) {
      c = c * static_cast<double>(q - l + 1) / static_cast<double>(l * (2 * q - l + 1));
      X = aux_vect[l] * (c * pow(xmod,l));
      expm = expm + X;
      if (pind) {
        D =  D + X;
      }
      else {
        D = D - X;
      }
      pind = !pind;
    }
    expm = inv(D) * expm;
    for (int l = 1; l <= s; ++l) {
      expm = expm * expm;
    }
    
    aux_mat = alpha.t() * expm * exit_vect;
    density = aux_mat(0,0);
    logLh += weight[k] * (std::log(density) + std::log(scale1[k]) + obs[k] * beta);
  }
  //Right censored data
  for (int k{0}; k < rcens.size(); ++k) {
    // Matrix exponential
    int pind{1};
    int ee{static_cast<int>(log2(s_norm  * scale2[k] * (exp(rcens[k] * beta) - 1) / beta)) + 1};
    s = std::max(0, ee + 1);
    xmod = scale2[k] * (exp(rcens[k] * beta) - 1) / beta / pow(2.0, s);
    c = 0.5;
    X = aux_vect[1] * (c * xmod);
    
    expm = aux_vect[0] + X;
    D = aux_vect[0] - X;
    
    for (int l{2}; l <= q; ++l) {
      c = c * static_cast<double>(q - l + 1) / static_cast<double>(l * (2 * q - l + 1));
      X = aux_vect[l] * (c * pow(xmod,l));
      expm = expm + X;
      if (pind) {
        D =  D + X;
      }
      else {
        D = D - X;
      }
      pind = !pind;
    }
    expm = inv(D) * expm;
    for (int l = 1; l <= s; ++l) {
      expm = expm * expm;
    }
    
    aux_mat = alpha.t() * expm * e;
    density = aux_mat(0,0);
    logLh += rcweight[k] * std::log(density);
  }
  
  return logLh;
}
