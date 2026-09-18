/**
 * tkmeans2 - K-means clustering with trimming: C++ Functions
 *
 */

//  Rcpp::compileAttributes("C:/users/valen/onedrive/myrepo/r/tclust")

// #include <RcppCommon.h>
// #include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

#include "tclust_types.h"

#include <Rcpp.h>

using namespace Rcpp;

/**
 * Calculates the initial cluster assignment and initial values for the parameters.
 *
 * @param x: matrix with observations and features.
 * @param iter: a reference to the cluster information. Its values are modified.
 * @param pa: a reference to the procedure parameters.
 */
void tkmeans_initClusters(arma::mat x, iteration &iter, params &pa, 
    Rcpp::Nullable<Rcpp::NumericMatrix> points = R_NilValue)
{
    arma::mat iter_center = arma::mat(pa.p, pa.k);        // Column ki stores the centers of cluster p
    arma::vec size;                                       // Cluster sizes
    arma::vec weights;                                    // Cluster weigths
    
    if((points.isNotNull())) {
        arma::mat points_ = Rcpp::as<arma::mat>(points);
        iter.centers = points_;

        // Rcout << "tkmeans_initClusters_c1: points is NOT NULL" << std::endl;
        // points_.print("Points ...");

    } else {
        // Rcout << "tkmeans_initClusters_c1: points is NULL" << std::endl;

        arma::ivec idx = arma::randi(pa.k * (pa.p + 1), arma::distr_param(0, pa.n - 1));
    
        for(int ki = 0; ki < pa.k; ki++) {
            // Select the p+1 next elements of idx
            arma::uvec rows_to_select = arma::conv_to<arma::uvec>::from(idx.rows(ki * (pa.p + 1), ki * (pa.p + 1) + pa.p));
        
            // Set cluster centers to random observation
            iter_center.col(ki) = x.row(rows_to_select(0)).t();
        }
    
        iter.centers = iter_center.t();
    }
    
    size = arma::vec(pa.k, arma::fill::value(pa.no_trim / pa.k));
    weights = arma::vec(pa.k, arma::fill::value(1 / (double)pa.k));

    iter.size = size;
    iter.weights = weights;
}

/**
 * Estimate the model parameters.
 *
 * @param x: matrix with observations and features.
 * @param iter: a reference to the cluster information. Its values are modified.
 * @param pa: a reference to the procedure parameters.
 */
void tkmeans_estimClustPar(arma::mat x, iteration &iter, params &pa)
{
    for (int ki = 0; ki < pa.k; ki++) {
        if(iter.size(ki) > pa.zero_tol) {
            iter.centers.row(ki) = (iter.posterior.col(ki).t() * x) / iter.size(ki);
        }
        else {
          iter.centers.row(ki) = arma::mat(1, pa.p);
        }
    }
}

arma::vec rowSums(const arma::mat &X){
    int nRows = X.n_rows;
    arma::vec out(nRows);
    for(int i = 0; i < nRows; i++){
        out(i) = sum(X.row(i));
    }
    return(out);
}

/**
 * Apply niter concentration steps to initial solution given by iter and pa.
 *
 * @param niter: number of concentration steps.
 * @param x: matrix with observations and features.
 * @param iter: a reference to the cluster information. Its values are modified.
 * @param pa: a reference to the procedure parameters.
 */
void tkmeans_csteps(int niter, arma::mat x, iteration &iter, params &pa)
{
    int k = pa.k;
    int n = pa.n;
    int no_trim = pa.no_trim;

    for(int i1 = 0; i1 < niter; i1++) {
    
        // estimates the cluster's assigment and trimming
        // findClustAssig(x, iter, pa); 
        
        arma::mat ll(n, k);
        for(int ki = 0; ki < k; ki++)   {
            arma::mat X_c = x;
            X_c.each_row() -= iter.centers.row(ki);
            X_c = arma::pow(X_c, 2);
            
            ll.col(ki) = rowSums(X_c);
        }

        arma::uvec old_assig = iter.cluster;
        arma::vec pre_z = arma::min(ll, 1);

        // Find outliers: determine elements to trim
        arma::uvec tc_set;
        arma::uvec sorted_index = arma::sort_index(pre_z, "ascending");
        arma::uvec last_indexes = arma::linspace<arma::uvec>(no_trim, n - 1, n - no_trim);
        arma::uvec obs_to_trim = sorted_index.elem(last_indexes);
        pre_z.elem(obs_to_trim) = arma::zeros<arma::vec>(n - no_trim);
        tc_set = pre_z > 0;

        // Cluster assignment with trimming
        iter.cluster = (arma::index_min(ll, 1) + 1) % tc_set;
        iter.disttom = pre_z;
        
        // Find posterior matrix
        arma::uvec one_to_n = arma::linspace<arma::uvec>(0, n - 1, n);
        arma::uvec aux_assig = iter.cluster - 1 + (iter.cluster == 0);
        
        // 2xn matrix containing (observation, cluster) pairs in each column
        arma::umat subscripts = (arma::join_rows(one_to_n, aux_assig)).t();
        
        iter.posterior = arma::mat(n, k);
        iter.posterior.elem(arma::sub2ind(arma::size(iter.posterior), subscripts)) = arma::ones<arma::vec>(n);

        iter.posterior.each_col() %= arma::conv_to<arma::vec>::from(tc_set); // Set to 0 all trimmed rows
    
        // calculate cluster size
        iter.size = (arma::sum(iter.posterior, 0)).t();
         
        // stop criterion
        if(arma::all(old_assig == iter.cluster)) {
            iter.code = 2;
            break;
        }
      
        // estimates the cluster's parameters
        tkmeans_estimClustPar(x, iter, pa); 
    }
  
    // calculates the objective function value
    // calcObj(x, iter, pa); 
    iter.obj = arma::sum(iter.disttom)/no_trim;
}

// Internal function for concentration steps (initializations) in tkmeans2
// @name tkmeans_c1
// @param x Rcpp::NumericMatrix, The input data.
// @param k The number of clusters initially searched for.
// @param alpha double, The proportion of observations to be trimmed.
// @param niter1 int, The number of concentration steps to be performed for the 
//     nstart initializations. 
// @param zero_tol The zero tolerance used. By default set to 1e-16.
// @param points Optional initial mean vectors, \code{NULL} or a matrix with \code{k} 
//  vectors used as means to initialize the algorithm. If initial mean vectors are 
//  specified, \code{nstart} should be 1 (otherwise the same initial means are 
//  used for all runs).
// @export
// [[Rcpp::export]]
Rcpp::List tkmeans_c1(arma::mat x, int k, double alpha = 0.05,
                     int niter1 = 3, double zero_tol = 1e-16, 
                     Rcpp::Nullable<Rcpp::NumericMatrix> points = R_NilValue)
{

    int n = x.n_rows;
    int p = x.n_cols;
    int no_trim = std::floor(n * (1 - alpha));
    
    params pa;
    pa.n = n;
    pa.p = p;
    pa.alpha = alpha;
    pa.no_trim = no_trim;
    pa.trimm = n - no_trim;
    pa.k = k;
    pa.equal_weights = true;
    pa.zero_tol = zero_tol;
    
    iteration iter;
    iter.obj = 0.0,
    iter.cluster = arma::uvec(n);
    iter.size = arma::vec(k);
    iter.centers = arma::mat(k, p);
    iter.weights = arma::vec(k);
    iter.code = 0;
    
    tkmeans_initClusters(x, iter, pa, points);     // Cluster random initialization
    tkmeans_csteps(niter1, x, iter, pa);            // Apply niter1 concentration steps
    
    return Rcpp::List::create(
      _["obj"] = iter.obj,
      _["cluster"] = iter.cluster);
}

// Internal function for concentration steps (refinement) in tkmeans2
// @name tkmeans_c2
// @param x Rcpp::NumericMatrix, The input data.
// @param k The number of clusters initially searched for.
// @param cluster arma::uvec A numerical vector of size \code{n} containing the 
//      cluster assignment for each observation. Cluster names are integer numbers 
//      from 1 to k, 0 indicates trimmed observations. 
// @param alpha double, The proportion of observations to be trimmed.
// @param niter2 The maximum number of concentration steps to be performed for the 
//  \code{nkeep} solutions kept for further iteration. The concentration steps are 
//  stopped, whenever two consecutive steps lead to the same data partition.
// @param zero_tol The zero tolerance used. By default set to 1e-16.
// @export
// [[Rcpp::export]]
iteration tkmeans_c2(arma::mat x, int k, arma::uvec cluster, double alpha = 0.05,
                     int niter2 = 20, double zero_tol = 1e-16)
{

    int n = x.n_rows;
    int p = x.n_cols;
    int no_trim = std::floor(n * (1 - alpha));
    
    // Build the posterior matrix
    arma::uvec tc_set;
    arma::uvec one_to_n = arma::linspace<arma::uvec>(0, n - 1, n);
    arma::uvec aux_assig = cluster - 1 + (cluster == 0);
    tc_set = cluster > 0;         // not trimmed observations
    
    // 2xn matrix containing (observation, cluster) pairs in each column
    arma::umat subscripts = (arma::join_rows(one_to_n, aux_assig)).t();
    
    arma::mat posterior(n, k);
    posterior.elem(arma::sub2ind(arma::size(posterior), subscripts)) = arma::ones<arma::vec>(n);

    //VT::04.05.2024 - Fix the issue reported by Domenico:
    //    all trimmed observations were assigned to the first class (0)
    posterior.each_col() %= arma::conv_to<arma::vec>::from(tc_set); // Set to 0 all trimmed rows
    
    arma::vec size = (arma::sum(posterior, 0)).t();

    params pa;
    pa.n = n;
    pa.p = p;
    pa.alpha = alpha;
    pa.no_trim = no_trim;
    pa.trimm = n - no_trim;
    pa.k = k;
    pa.equal_weights = TRUE;
    pa.zero_tol = zero_tol;
    
    iteration iter;
    iter.centers = arma::mat(k, p);
    iter.cluster = cluster;
    iter.obj = 0.0,
    iter.size = size;
    iter.weights = size/no_trim;
    iter.code = 0;
    iter.posterior = posterior;
    
    tkmeans_estimClustPar(x, iter, pa);
    tkmeans_csteps(niter2, x, iter, pa);
    return iter;
}

