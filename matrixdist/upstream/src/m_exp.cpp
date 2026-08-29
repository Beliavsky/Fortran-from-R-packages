#include "m_exp.h"
# include <RcppArmadillo.h>
// [[ Rcpp :: depends ( RcppArmadillo )]]


//' L inf norm of a matrix
//' 
//' Computes the L inf norm of a matrix \code{A}, which is defined as:
//' L_inf(A) =  max(1 <= i <= M) sum(1 <= j <= N) abs(A(i,j)).
//' 
//' @param A A matrix.
//' @return The L inf norm.
//' 
// [[Rcpp::export]]
double inf_norm(arma::mat A) {
  double value{0.0};
  
  for (int i{0}; i < A.n_rows; ++i) {
    double row_sum{0.0};
    for (int j{0}; j < A.n_cols; ++j) {
      row_sum += fabs(A(i,j));
    }
    value = std::max(value, row_sum);
  }
  return value;
}


//' Creates the matrix  (A1, B1 ; 0, A2)
//' 
//' @param A1 Matrix.
//' @param A2 Matrix.
//' @param B1 Matrix.
//' @return Computes (A1, B1 ; 0, A2).
//' 
// [[Rcpp::export]]
arma::mat matrix_vanloan(arma::mat A1, arma::mat A2, arma::mat B1) {
  unsigned p1{A1.n_rows};
  unsigned p2{A2.n_rows};
  unsigned p{p1 + p2};
  
  arma::mat aux_mat(p, p);
  
  for (int i{0}; i < p; ++i) {
    for (int j{0}; j < p; ++j) {
      if (i < p1 && j < p1) {
        aux_mat(i,j) = A1(i,j);
      }
      else if (i >= p1 && j < p1) {
        aux_mat(i,j) = 0;
      }
      else if (i < p1 && j >= p1) {
        aux_mat(i,j) = B1(i,j - p1);
      }
      else {
        aux_mat(i,j) = A2(i - p1,j - p1);
      }
    }
  }
  return aux_mat;
}


//' Maximum diagonal element of a matrix
//' 
//' @param A Matrix.
//' @return The maximum value in the diagonal.
//' 
// [[Rcpp::export]]
double max_diagonal(const arma::mat & A) {
  double maximum{A(0,0)};
  for (int i{0}; i < A.n_rows; ++i) {
    if (A(i,i) > maximum) {
      maximum = A(i,i);
    }
  }
  return maximum;
}


//' Matrix exponential
//' 
//' MATLAB's built-in algorithm for matrix exponential - Pade approximation.
//' 
//' @param A A matrix.
//' @return exp(A).
//' 
// [[Rcpp::export]]
arma::mat matrix_exponential(arma::mat A) {
  const int q{6};
  
  arma::mat expm(A.n_rows,A.n_cols);
  
  double a_norm{inf_norm(A)};
  
  int ee{static_cast<int>(log2(a_norm)) + 1};
  
  int s{std::max(0, ee + 1)};
  
  double t{1.0 / pow(2.0, s)};
  
  arma::mat a2 = A * t;
  arma::mat x = a2;
  
  double c{0.5};
  
  expm.eye(size(A));
  expm = expm + (a2 * c);
  
  arma::mat d; 
  d.eye(size(A));
  d = (d + a2 * (-c));
  
  int p{1};
  
  for (int k{2}; k <= q; ++k) {
    c = c * static_cast<double>(q - k + 1) / static_cast<double>(k * (2 * q - k + 1));
    x = (a2 * x);
    expm = (x * c) + expm;
    if (p) {
      d = (x * c) + d;
    }
    else {
      d = (x * (-c)) + d;
    }
    p = !p;
  }
  //  E -> inverse(D) * E
  expm = inv(d) * expm;
  //  E -> E^(2*s)
  for (int k = 1; k <= s; ++k) {
    expm = expm * expm;
  }
  return(expm);
}


//' Matrix exponential
//' 
//' Armadillo matrix exponential implementation.
//' 
//' @param A A matrix.
//' @return exp(A).
//' 
// [[Rcpp::export]]
arma::mat expmat(arma::mat A) {
  arma::mat B = arma::expmat(A);  
  return(B);
}


//' Computes A^n
//'
//' @param A A matrix.
//' @param n An integer.
//' @return A^n.
//' @export
// [[Rcpp::export]]
arma::mat matrix_power(int n, arma::mat A) {
  if (n == 0) {
    arma::mat d;
    d.eye(size(A));
    return(d);
  }
  else if (n == 1) {
    return A;
  }
  else if (n == 2) {
    return A * A;
  }
  arma::mat previous_matrix = A * A;
  arma::mat new_matrix = A * previous_matrix;
  for (int i{4}; i <= n; ++i) {
    previous_matrix = new_matrix;
    new_matrix = A * previous_matrix;
  }
  return new_matrix;
}


//' Computes elements A^n until the given size
//'
//' @param A A matrix.
//' @param vect_size Size of vector.
//'
// [[Rcpp::export]]
std::vector<arma::mat> vector_of_powers(const arma::mat & A, int vect_size) {
  arma::mat Id;
  Id.eye(size(A));
  
  std::vector<arma::mat> vect;
  
  vect.push_back(Id);
  
  for (int k{1}; k <= vect_size; ++k) {
    vect.push_back(A * vect[k - 1]);
  }
  return vect;
}

