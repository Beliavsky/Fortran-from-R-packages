#include <RcppCommon.h>
#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

// Input parameters of the procedure
struct params {
  int n;              // Number of observations
  int p;              // Number of features
  double alpha;       // Trimming level
  int trimm;          // Number of observations to trim
  int no_trim;        // Number of observations to cluster
  int k;              // Number of clusters
  bool equal_weights; // Equal population proportion for all clusters
  double zero_tol;    // Tolerance that substitutes 0 (to avoid some numerical issues)
  int restrC;         // Restriction type (0 = eigenvalues or determinants)  
  bool deterC;        // Determinants restriction  
  double restr_fact;  // Level of constraint (of eigenvalues or determinants)
  double cshape;      // Level of constraint of the eigenvalues when deterC=true
  Rcpp::String opt;   // Model of estimate, can be "HARD" for hard assignment or "MIXTURE" for mixture assignment
};

// Cluster information
struct iteration {
  arma::mat centers;    // Cluster centers
  arma::cube cov;       // Cluster covariance matrices
  arma::uvec cluster;   // Cluster assignation indices
  arma::vec disttom;    // Distances to cluster center
  double obj;           // Value of the objective function
  double NlogL;         // Value of the the negative of the CLASSIFICATION LOG-LIKELIHOOD  of the untrimmed units
                        // NlogL = -sum(max(ll(untrimmed units,[],2));
  arma::vec size;       // Cluster sizes
  arma::vec weights;    // Cluster weights
  int code;             // A return code signaling particular situations (e.g. data are aligned)
  arma::mat posterior;  // Cluster assignment given by 0/1 columns
  arma::vec lmd;
  arma::cube OMG;
};

// Input parameters of the GPCM constraints
struct GPCMPars {

    int p;                      // number of variables
    int k;                      // number of clusters

    std::string pars;           // type of Gaussian Parsimonious Clustering Model

    double cdet;                // restriction to be applied to the determinants
    double shw;                 // restriction to be applied to the elements of
                                //   the shape matrices inside each group
    double shb;                 // restriction which has to be applied to the elements of
                                //   the shape matrices across each group

    double tolDSR;              // tolerance to use to exit the loop for obtaining 
                                // the requested restricted determinants, shape matrices and rotation             
    double tolR;                // tolerance to use to exit the iterations to obtain the 
                                //   common rotation matrix in presence of varying shape.
    double tolS;                // tolerance to use to exit the iterative procedure for estimating the shape

    double zerotol;             // tolerance value to declare all input values equal 
                                //   to 0 in the eigenvalues restriction routine

    int maxiterDSR;             // maximum number of iterations to obtain the requested restricted
                                //   determinants, shape matrices and rotation.
    int maxiterR;               // maximum number of iterations to obtain the common rotation matrix
                                //   in presence of varying shape
    int maxiterS;               // maximum number of iterations to obtain the restricted shape matrix.

    int sortsh;                 // wheather tp sort the shape matrix when comparing
};

//  Helper for instantiation of the GPCMPars structure
inline GPCMPars parsePars(const Rcpp::List& pa)
{
    GPCMPars x;

    x.p = Rcpp::as<int>(pa["p"]);
    x.k = Rcpp::as<int>(pa["k"]);

    x.pars = Rcpp::as<std::string>(pa["pars"]);

    x.cdet = Rcpp::as<double>(pa["cdet"]);
    x.shw  = Rcpp::as<double>(pa["shw"]);
    x.shb  = Rcpp::as<double>(pa["shb"]);

    x.tolDSR = Rcpp::as<double>(pa["tolDSR"]);
    x.tolR   = Rcpp::as<double>(pa["tolR"]);
    x.tolS   = Rcpp::as<double>(pa["tolS"]);

    x.zerotol = Rcpp::as<double>(pa["zerotol"]);

    x.maxiterDSR = Rcpp::as<int>(pa["maxiterDSR"]);
    x.maxiterR   = Rcpp::as<int>(pa["maxiterR"]);
    x.maxiterS   = Rcpp::as<int>(pa["maxiterS"]);

    x.sortsh = Rcpp::as<int>(pa["sortsh"]);

    return x;
}

namespace Rcpp {

// Support for wrap
template <> SEXP wrap(const iteration& iter);
template <> SEXP wrap(const params& pa);

}
