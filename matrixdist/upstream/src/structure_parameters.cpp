#include <Rcpp.h>
using namespace Rcpp;


//' Random structure of a phase-type
//' 
//' Generates random parameters \code{alpha} and \code{S} of a phase-type 
//'  distribution of dimension \code{p} with chosen structure.
//'  
//' @param p Dimension of the phase-type.
//' @param structure Type of structure: "general", "hyperexponential", "gerlang",
//'  "coxian" or "gcoxian".
//' @param scale_factor A factor that multiplies the sub-intensity matrix.
//' @return Random parameters \code{alpha} and \code{S} of a phase-type.
//' 
// [[Rcpp::export]]
List random_structure(int p, String structure = "general", double scale_factor = 1) {
  // Structure of alpha and S
  NumericVector alpha_legal(p);
  NumericMatrix S_legal(p, p);
  
  if (structure == "general") {
    for (int i{0}; i < p; ++i) {
      alpha_legal[i] = 1;
      for (int j = 0; j < p; ++j) {
        S_legal(i, j) = 1;
      }
    }
  }
  else if (structure == "hyperexponential") {
    for (int i{0}; i < p; i++) {  
      alpha_legal[i] = 1;
      S_legal(i, i) = 1;
    }
  }
  else if (structure == "gerlang") {
    alpha_legal[0] = 1;
    for (int i{0}; i < p - 1; ++i) {
      S_legal(i, i + 1) = 1;
    } 
    S_legal(p - 1, p - 1) = 1;
  }
  else if (structure == "coxian") {
    alpha_legal[0] = 1;
    for (int i{0}; i < p - 1; ++i) { 
      S_legal(i, i) = 1;
      S_legal(i, i + 1) = 1;
    }
    S_legal(p - 1, p - 1) = 1; 
  }
  else if (structure == "gcoxian") {
    for (int i{0}; i < p - 1; ++i) { 
      alpha_legal[i] = 1;
      S_legal(i, i) = 1;
      S_legal(i, i + 1) = 1;
    }
    S_legal(p - 1, p - 1) = 1;
    alpha_legal[p - 1] = 1;
  }
  else{
    Rcerr << "non-existent structure\n";
  }
  
  // Random initialization
  NumericVector alpha_rnd(p);
  NumericMatrix S_rnd(p, p);
  
  double sum{0.0};
  
  for (int i{0}; i < p; ++i) {
    if (alpha_legal[i] == 1) {
      alpha_rnd[i] = runif(1)[0];
      sum += alpha_rnd[i];
    }
  }
  
  for (int i{0}; i < p; ++i) {
    alpha_rnd[i] = alpha_rnd[i] / sum;
  }
  
  for (int i{0}; i < p; ++i) {
    for (int j{0}; j < p; ++j) {
      if ((i != j) && (S_legal(i,j) == 1)) {
        S_rnd(i,j) = runif(1)[0];
        S_rnd(i,i) -= S_rnd(i,j);
      }
    }
  }
  
  // exit rate
  for (int i{0}; i < p; ++i) {
    if (S_legal(i,i) == 1) {
      S_rnd(i,i) -= runif(1)[0];
    }
  }
  S_rnd = S_rnd * scale_factor;
  
  List L = List::create(Named("alpha") = alpha_rnd , _["S"] = S_rnd);
  
  return L;
}


//' Random structure of a bivariate phase-type
//'
//' Generates random parameters \code{alpha}, \code{S11}, \code{S12}, and \code{S22}
//' of a bivariate phase-type distribution of dimension \code{p  = p1 + p2}.
//'
//' @param p1 Dimension of the first block.
//' @param p2 Dimension of the second block.
//' @param scale_factor A factor that multiplies the sub-intensity matrix.
//' @return Random parameters  \code{alpha}, \code{S11}, \code{S12}, and \code{S22}
//'  of a bivariate phase-type.
//'
// [[Rcpp::export]]
List random_structure_bivph(int p1, int p2, double scale_factor = 1) {
  NumericVector alpha(p1);
  NumericMatrix S11(p1, p1);
  NumericMatrix S12(p1, p2);
  NumericMatrix S22(p2, p2);
  
  double sum{0.0};
  
  for (int i{0}; i < p1; ++i) {
    alpha[i] = runif(1)[0];
    sum += alpha[i];
  }
  alpha = alpha / sum;
  
  for (int i{0}; i < p1; ++i) {
    for (int j{0}; j < p1; ++j) {
      if (i != j) {
        S11(i,j) = runif(1)[0];
        S11(i,i) -= S11(i,j);
      }
    }
    for (int j{0}; j < p2; ++j) {
      S12(i,j) = runif(1)[0];
      S11(i,i) -= S12(i,j);
    }
  }
  
  for (int i{0}; i < p2; ++i) {
    for (int j{0}; j < p2; ++j) {
      if (i != j) {
        S22(i,j) = runif(1)[0];
        S22(i,i) -= S22(i,j);
      }
    }
  }
  
  for (int i{0}; i < p2; ++i) {
    S22(i,i) -= runif(1)[0];
  }
  
  S11 = S11 * scale_factor;
  S12 = S12 * scale_factor;
  S22 = S22 * scale_factor;
  
  List L = List::create(Named("alpha") = alpha, _["S11"] = S11, _["S12"] = S12, _["S22"] = S22);
  
  return L;
}


//' Merges the matrices S11, S12 and S22 into a sub-intensity matrix
//'
//' @param S11 A sub-intensity matrix.
//' @param S12 A matrix.
//' @param S22 A sub-intensity matrix.
//' @return A sub-intensity matrix.
//'
// [[Rcpp::export]]
NumericMatrix merge_matrices(NumericMatrix S11, NumericMatrix S12, NumericMatrix S22) {
  long p1{S11.nrow()};
  long p2{S22.nrow()};
  NumericMatrix S(p1 + p2, p1 + p2);
  
  for (int i{0}; i < p1; ++i) {
    for (int j{0}; j < p1; ++j) {
      S(i,j) = S11(i,j);
    }
    for (int j{0}; j < p2; ++j) {
      S(i,j + p1) = S12(i,j);
    }
  }
  for (int i{0}; i < p2; ++i) {
    for (int j{0}; j < p2; ++j) {
      S(i + p1,j + p1) = S22(i,j);
    }
  }
  return S;
}
