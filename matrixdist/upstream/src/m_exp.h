# include <RcppArmadillo.h>

double inf_norm(arma::mat A);

arma::mat matrix_vanloan(arma::mat A1, arma::mat A2, arma::mat B1);

double max_diagonal(const arma::mat & A);

arma::mat matrix_exponential(arma::mat A);

arma::mat matrix_power(int n, arma::mat A);

std::vector<arma::mat> vector_of_powers(const arma::mat & A, int vect_size);
