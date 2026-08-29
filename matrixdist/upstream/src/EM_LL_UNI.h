#ifndef EM_LL_UNI   // if x.h hasn't been included yet...
#define EM_LL_UNI   //   #define this so the compiler knows it has been included

#include <RcppArmadillo.h>

void vector_of_matrices(std::vector<arma::mat> & vect, const arma::mat & S, double a, int vect_size);

arma::mat m_exp_sum(double x, int n, const std::vector<arma::mat> & pow_vector, double a);

void pow2_matrix(int n , arma::mat & A);

int find_n(double h, double lambda);
    
void EMstep_UNI(double h, arma::vec & alpha, arma::mat & S, const Rcpp::NumericVector & obs, const Rcpp::NumericVector & weight, const Rcpp::NumericVector & rcens, const Rcpp::NumericVector & rcweight);

////////////////////////////////////////////
// Log-likelihoods
////////////////////////////////////////////

double logLikelihoodPH_UNI(double h, arma::vec & alpha, arma::mat & S, const Rcpp::NumericVector & obs, const Rcpp::NumericVector & weight, const Rcpp::NumericVector & rcens, const Rcpp::NumericVector & rcweight);
  
double logLikelihoodMweibull_UNI(double h, arma::vec & alpha, arma::mat & S, double beta , const Rcpp::NumericVector & obs, const Rcpp::NumericVector & weight, const Rcpp::NumericVector & rcens, const Rcpp::NumericVector & rcweight);
    
double logLikelihoodMpareto_UNI(double h, arma::vec & alpha, arma::mat & S, double beta , const Rcpp::NumericVector & obs, const Rcpp::NumericVector & weight, const Rcpp::NumericVector & rcens, const Rcpp::NumericVector & rcweight);
      
double logLikelihoodMlognormal_UNI(double h, arma::vec & alpha, arma::mat & S, double beta , const Rcpp::NumericVector & obs, const Rcpp::NumericVector & weight, const Rcpp::NumericVector & rcens, const Rcpp::NumericVector & rcweight);
        
double logLikelihoodMloglogistic_UNI(double h, arma::vec & alpha, arma::mat & S, Rcpp::NumericVector beta , const Rcpp::NumericVector & obs, const Rcpp::NumericVector & weight, const Rcpp::NumericVector & rcens, const Rcpp::NumericVector & rcweight);

double logLikelihoodMgompertz_UNI(double h, arma::vec & alpha, arma::mat & S, double beta , const Rcpp::NumericVector & obs, const Rcpp::NumericVector & weight, const Rcpp::NumericVector & rcens, const Rcpp::NumericVector & rcweight);

double logLikelihoodMgev_UNI(double h, arma::vec & alpha, arma::mat & S, Rcpp::NumericVector beta , const Rcpp::NumericVector & obs, const Rcpp::NumericVector & weight, const Rcpp::NumericVector & rcens, const Rcpp::NumericVector & rcweight);

/////////////////////////////////////////////////////////
// Scaled versions of loglikelihoods (for regression)  //
/////////////////////////////////////////////////////////

double logLikelihoodPH_UNIs(double h, arma::vec & alpha, arma::mat & S, const Rcpp::NumericVector & obs, const Rcpp::NumericVector & weight, const Rcpp::NumericVector & rcens, const Rcpp::NumericVector & rcweight, const Rcpp::NumericVector & scale1, const Rcpp::NumericVector & scale2);

double logLikelihoodMweibull_UNIs(double h, arma::vec & alpha, arma::mat & S, double beta , const Rcpp::NumericVector & obs, const Rcpp::NumericVector & weight, const Rcpp::NumericVector & rcens, const Rcpp::NumericVector & rcweight, const Rcpp::NumericVector & scale1, const Rcpp::NumericVector & scale2);

double logLikelihoodMpareto_UNIs(double h, arma::vec & alpha, arma::mat & S, double beta , const Rcpp::NumericVector & obs, const Rcpp::NumericVector & weight, const Rcpp::NumericVector & rcens, const Rcpp::NumericVector & rcweight, const Rcpp::NumericVector & scale1, const Rcpp::NumericVector & scale2);
  
double logLikelihoodMlognormal_UNIs(double h, arma::vec & alpha, arma::mat & S, double beta , const Rcpp::NumericVector & obs, const Rcpp::NumericVector & weight, const Rcpp::NumericVector & rcens, const Rcpp::NumericVector & rcweight, const Rcpp::NumericVector & scale1, const Rcpp::NumericVector & scale2);
    
double logLikelihoodMloglogistic_UNIs(double h, arma::vec & alpha, arma::mat & S, Rcpp::NumericVector beta , const Rcpp::NumericVector & obs, const Rcpp::NumericVector & weight, const Rcpp::NumericVector & rcens, const Rcpp::NumericVector & rcweight, const Rcpp::NumericVector & scale1, const Rcpp::NumericVector & scale2);

double logLikelihoodMgompertz_UNIs(double h, arma::vec & alpha, arma::mat & S, double beta , const Rcpp::NumericVector & obs, const Rcpp::NumericVector & weight, const Rcpp::NumericVector & rcens, const Rcpp::NumericVector & rcweight, const Rcpp::NumericVector & scale1, const Rcpp::NumericVector & scale2);

  
#endif
