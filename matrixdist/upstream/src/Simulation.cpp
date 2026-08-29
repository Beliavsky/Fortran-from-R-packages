#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]


//' Embedded Markov chain of a sub-intensity matrix
//' 
//' Returns the transition probabilities of the embedded Markov chain determined
//'  the sub-intensity matrix.
//'  
//' @param S A sub-intensity matrix.
//' @return The embedded Markov chain.
//' 
// [[Rcpp::export]]
arma::mat embedded_mc(arma::mat S) {
  unsigned p{S.n_rows};
  arma::mat Q(p + 1, p + 1);
  
  arma::mat e;
  e.ones(S.n_cols, 1);
  arma::mat exit_vect = (S * (-1)) * e;
  
  for (int i{0}; i < p; ++i) {
    for (int j{0}; j < p + 1; ++j) {
      if (j != i && j < p) {
        Q(i,j) = -1.0 * S(i,j) / S(i,i);
      }
      else if(j == p) {
        Q(i,j) = -1.0 * exit_vect(i,0) / S(i,i);
      }
    }
  }
  Q(p,p) = 1;
  
  return (Q);
}


//' Cumulate matrix
//'
//' Creates a new matrix with entries the cumulated rows of \code{A}.
//' 
//' @param A A matrix.
//' @return The cumulated matrix.
//'
// [[Rcpp::export]]
arma::mat cumulate_matrix(arma::mat A) {
  unsigned p1{A.n_rows};
  unsigned p2{A.n_cols};
  
  arma::mat cumulated(p1, p2);
  
  for (int i{0}; i < p1; ++i) {
    for (int j{0}; j < p2; ++j) {
      if (j == 0) {
        cumulated(i,j) = A(i,j);
      }
      else {
        cumulated(i,j) = cumulated(i,j - 1) + A(i,j);
      }
    }
  }
  return cumulated;
}


//' Cumulate vector
//'
//' Creates a new vector with entries the cumulated entries of \code{A}.
//' 
//' @param A A vector.
//' @return The cumulated vector.
//'
// [[Rcpp::export]]
arma::vec cumulate_vector(arma::vec A) {
  unsigned p{A.size()};
  
  arma::vec cumulated(p);
  
  for (int i{0}; i < p; ++i) {
    if (i == 0) {
      cumulated[i] = A[i];
    }
    else {
      cumulated[i] = cumulated[i - 1] + A[i];
    }
  }
  return cumulated;
}


//' Initial state of Markov jump process
//'
//' Given the accumulated values of the initial probabilities \code{alpha} and a
//'  uniform value \code{u}, it returns the initial state of a Markov jump process.
//' This corresponds to the states satisfying cum_alpha_(k-1) < u < cum_alpha_(k).
//' 
//' @param cum_alpha A cummulated vector of initial probabilities.
//' @param u Random value in (0,1).
//' @return Initial state of the Markov jump process.
//'
// [[Rcpp::export]]
long initial_state(arma::vec cum_alpha, double u) {
  if (u <= cum_alpha[0]) {
    return 0;
  }
  
  for( int i{1}; i < cum_alpha.size(); ++i) {
    if (cum_alpha[i - 1] < u && u <= cum_alpha[i]) {
      return i;
    }
  }
  return 0;
}


//' New state in a Markov jump process
//'
//' Given a transition matrix \code{Q}, a uniform value \code{u}, and a previous
//'  state \code{k}, it returns the new state of a Markov jump process.
//'  
//' @param prev_state Previous state of the Markov jump process.
//' @param cum_embedded_mc Transition matrix.
//' @param u Random value in (0,1).
//' @return Next state of the Markov jump process.
//'
// [[Rcpp::export]]
long new_state(long prev_state, arma::mat cum_embedded_mc, double u) {
  if (u <= cum_embedded_mc(prev_state,0)) {
    return 0;
  }
  
  for (int i{1}; i < cum_embedded_mc.n_cols; ++i) {
    if (cum_embedded_mc(prev_state,i - 1) < u && u <= cum_embedded_mc(prev_state,i)) {
      return i;
    }
  }
  return 0;
}


//' Simulate phase-type
//'
//' Generates a sample of size \code{n} from a phase-type distribution with
//' parameters \code{alpha} and \code{S}.
//' 
//' @param n Sample size.
//' @param alpha Vector of initial probabilities.
//' @param S Sub-intensity matrix.
//' @return Simulated sample.
//'
// [[Rcpp::export]]
Rcpp::NumericVector rphasetype(int n, arma::vec alpha, arma::mat S) {
  Rcpp::NumericVector sample(n);
  
  arma::mat cum_embedded_mc = cumulate_matrix(embedded_mc(S));
  arma::vec cum_alpha = cumulate_vector(alpha);
  
  unsigned p{alpha.size()};
  long state{0};
  for (int i{0}; i < n; ++i) {
    double time{0.0};
    state = initial_state(cum_alpha, Rcpp::runif(1)[0]);
    while (state != p) {
      time += log(1.0 - Rcpp::runif(1)[0]) / S(state,state);
      state = new_state(state, cum_embedded_mc, Rcpp::runif(1)[0]);
    }
    sample[i] = time;
  }
  return sample;
}


//' Random inhomogeneous phase-type
//' 
//' Generates a sample of size \code{n} from an inhomogeneous phase-type 
//' distribution with parameters \code{alpha}, \code{S} and \code{beta}.
//' 
//' @param n Sample size.
//' @param dist_type Type of IPH.
//' @param alpha Initial probabilities.
//' @param S Sub-intensity matrix.
//' @param beta Parameter of the transformation.
//' @return The simulated sample.
//' 
// [[Rcpp::export]]
Rcpp::NumericVector riph(int n, Rcpp::String dist_type, arma::vec alpha, arma::mat S, Rcpp::NumericVector beta) {
  Rcpp::NumericVector sample(n);
  
  arma::mat cum_embedded_mc = cumulate_matrix(embedded_mc(S));
  arma::vec cum_alpha = cumulate_vector(alpha);
  
  int p = alpha.size();
  long state{0};
  for (int i{0}; i < n; ++i) {
    double time{0.0};
    state = initial_state(cum_alpha, Rcpp::runif(1)[0]);
    while (state != p) {
      time += log(1.0 - Rcpp::runif(1)[0]) / S(state,state);
      state = new_state(state, cum_embedded_mc, Rcpp::runif(1)[0]);
    }
    if (dist_type == "pareto") {
      time = beta[0] * (exp(time) - 1);
    }
    else if (dist_type == "weibull") {
      time = pow(time, 1.0 / beta[0]);
    }
    else if (dist_type == "lognormal") {
      time =  exp(pow(time, 1.0 / beta[0])) - 1;
    }
    else if (dist_type == "loglogistic") {
      time = beta[0] * pow(exp(time) - 1, 1 / beta[1]);
    }
    else if (dist_type == "gompertz") {
      time = log(beta[0] * time + 1) / beta[0];
    }
    sample[i] = time;
  }
  return (sample);
}


//' Random matrix GEV
//' 
//' Generates a sample of size \code{n} from an inhomogeneous phase-type 
//' distribution with parameters \code{alpha}, \code{S} and \code{beta}.
//' 
//' @param n Sample size.
//' @param alpha Initial probabilities.
//' @param S Sub-intensity matrix.
//' @param mu Location parameter.
//' @param sigma Scale parameter.
//' @param xi Shape parameter: Default 0 which corresponds to the Gumbel case.
//' @return The simulated sample.
//' 
// [[Rcpp::export]]
Rcpp::NumericVector rmatrixgev(int n, arma::vec alpha, arma::mat S, double mu, double sigma, double xi = 0) {
  Rcpp::NumericVector sample(n);
  
  arma::mat cum_embedded_mc = cumulate_matrix(embedded_mc(S));
  arma::vec cum_alpha = cumulate_vector(alpha);
  
  int p = alpha.size();
  long state{0};
  for (int i{0}; i < n; ++i) {
    double time{0.0};
    state = initial_state(cum_alpha, Rcpp::runif(1)[0]);
    while (state != p) {
      time += log(1.0 - Rcpp::runif(1)[0]) / S(state,state);
      state = new_state(state, cum_embedded_mc, Rcpp::runif(1)[0]);
    }
    if (xi == 0) {
      time = mu - sigma * log(time);
    }
    else {
      time = mu + sigma  * (pow(time, - xi) - 1) / xi;
    }
    sample[i] = time;
  }
  return (sample);
}


//' Simulate a MPH* random vector
//'
//' Generates a sample of size \code{n} from a MPH* distribution with parameters
//'  \code{alpha}, \code{S} and \code{R}.
//'
//' @param n Sample size.
//' @param alpha Initial probabilities.
//' @param S Sub-intensity matrix.
//' @param R Reward matrix.
//' @return The simulated sample.
//' 
// [[Rcpp::export]]
Rcpp::NumericMatrix rMPHstar(int n, arma::vec alpha, arma::mat S, arma::mat R) {
  unsigned dim{R.n_cols};
  
  Rcpp::NumericMatrix sample(n, dim);
  
  arma::mat cum_embedded_mc = cumulate_matrix(embedded_mc(S));
  arma::vec cum_alpha = cumulate_vector(alpha);
  
  unsigned p{alpha.size()};
  long state{0};
  double time{0.0};
  for (int i = 0; i < n; ++i) {
    state = initial_state(cum_alpha, Rcpp::runif(1)[0]);
    while (state != p) {
      time = log(1.0 - Rcpp::runif(1)[0]) / S(state,state);
      for (int j{0}; j < dim; ++j) {
        sample(i,j) += R(state, j) * time;
      }
      state = new_state(state, cum_embedded_mc, Rcpp::runif(1)[0]);
    }
  }
  return (sample);
}


//' Simulate a MIPH* random vector
//'
//' Generates a sample of size \code{n} from a MIPH* distribution with parameters
//'  \code{alpha}, \code{S} and \code{R}.
//'
//' @param n Sample size.
//' @param alpha Initial probabilities.
//' @param S Sub-intensity matrix.
//' @param R Reward matrix.
//' @param gfun Vector with transformations names.
//' @param gfun_par List with transformations parameters.
//' @return The simulated sample.
//' 
// [[Rcpp::export]]
Rcpp::NumericMatrix rMIPHstar(int n, arma::vec alpha, arma::mat S, arma::mat R, Rcpp::StringVector gfun, Rcpp::List gfun_par) {
  unsigned dim{R.n_cols};
  
  Rcpp::NumericMatrix sample(n, dim);
  
  arma::mat cum_embedded_mc = cumulate_matrix(embedded_mc(S));
  arma::vec cum_alpha = cumulate_vector(alpha);
  
  unsigned p{alpha.size()};
  long state{0};
  double time{0.0};
  for (int i = 0; i < n; ++i) {
    state = initial_state(cum_alpha, Rcpp::runif(1)[0]);
    while (state != p) {
      time = log(1.0 - Rcpp::runif(1)[0]) / S(state,state);
      for (int j{0}; j < dim; ++j) {
        sample(i,j) += R(state, j) * time;
      }
      state = new_state(state, cum_embedded_mc, Rcpp::runif(1)[0]);
    }
    for (int j{0}; j < dim; ++j) {
      arma::vec beta = gfun_par[j];
      if (gfun[j] == "pareto") {
        sample(i,j) = beta[0] * (exp(sample(i,j)) - 1);
      }
      else if (gfun[j] == "weibull") {
        sample(i,j) = pow(sample(i,j), 1.0 / beta[0]);
      }
      else if (gfun[j] == "lognormal") {
        sample(i,j) =  exp(pow(sample(i,j), 1.0 / beta[0])) - 1;
      }
      else if (gfun[j] == "loglogistic") {
        sample(i,j) = beta[0] * pow(exp(sample(i,j)) - 1, 1 / beta[1]);
      }
      else if (gfun[j] == "gompertz") {
        sample(i,j) = log(beta[0] * sample(i,j) + 1) / beta[0];
      }
      else if (gfun[j] == "gev") {
        if (beta[2] == 0) {
          sample(i,j) = beta[0] - beta[1] * log(sample(i,j));
        }
        else {
          sample(i,j) = beta[0] + beta[1]  * (pow(sample(i,j), - beta[2]) - 1) / beta[2];
        }
      }
    }
  }
  return (sample);
}


//' Simulate discrete phase-type
//'
//' Generates a sample of size \code{n} from a discrete phase-type distribution with
//' parameters \code{alpha} and \code{S}.
//' 
//' @param n Sample size.
//' @param alpha Vector of initial probabilities.
//' @param S Sub-transition matrix.
//' @return Simulated sample.
//'
// [[Rcpp::export]]
Rcpp::NumericVector rdphasetype(int n, arma::vec alpha, arma::mat S) {
  Rcpp::NumericVector sample(n);
  
  arma::mat cum_trans = cumulate_matrix(S);
  arma::vec cum_alpha = cumulate_vector(alpha);
  
  unsigned p{alpha.size()};
  long state{0};
  for (int i{0}; i < n; ++i) {
    double time{0.0};
    state = initial_state(cum_alpha, Rcpp::runif(1)[0]);
    while (state != p) {
      time += 1;
      state = new_state(state, cum_trans, Rcpp::runif(1)[0]);
    }
    sample[i] = time;
  }
  return sample;
}


//' Simulate MDPH*
//'
//' Generates a sample of size \code{n} from a MDPH* distribution with
//' parameters \code{alpha}, \code{S}, and \code{R}.
//' 
//' @param n Sample size.
//' @param alpha Vector of initial probabilities.
//' @param S Sub-transition matrix.
//' @param R Reward matrix.
//' @return Simulated sample.
//'
// [[Rcpp::export]]
Rcpp::NumericMatrix rMDPHstar(int n, arma::vec alpha, arma::mat S, arma::mat R) {
  unsigned dim{R.n_cols};
  
  Rcpp::NumericMatrix sample(n, dim);
  
  arma::mat cum_trans = cumulate_matrix(S);
  arma::vec cum_alpha = cumulate_vector(alpha);
  
  unsigned p{alpha.size()};
  long state{0};
  for (int i{0}; i < n; ++i) {
    state = initial_state(cum_alpha, Rcpp::runif(1)[0]);
    while (state != p) {
      for (int j{0}; j < dim; ++j) {
        sample(i,j) += R(state, j);
      }
      state = new_state(state, cum_trans, Rcpp::runif(1)[0]);
    }
  }
  return sample;
}
