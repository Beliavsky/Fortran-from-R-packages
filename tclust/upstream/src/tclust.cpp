/**
 * Tclust - Clustering with trimming: C++ Functions
 *
 */

//  Rcpp::compileAttributes("C:/users/valen/onedrive/myrepo/r/tclust")

// VT::27.10.2024 / suppress the RcppArmadillo warnings, particularly 
//  chol() warning about non-symetric matrix
#define ARMA_WARN_LEVEL 0

// #include <RcppCommon.h>
// #include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

#include "tclust_types.h"

#include <Rcpp.h>

using namespace Rcpp;

Rcpp::List restrSigmaGPCM(arma::cube SigmaB, arma::vec niini, GPCMPars &pa,
                          bool trace=false, 
                          Rcpp::Nullable<arma::vec> lmd_=R_NilValue,
                          Rcpp::Nullable<arma::cube> OMG_=R_NilValue);

namespace Rcpp {
  template <> SEXP wrap(const iteration& iter) {
    return Rcpp::List::create(
      _["obj"] = iter.obj,
      _["NlogL"] = iter.NlogL,
      _["cluster"] = iter.cluster,
      _["disttom"] = iter.disttom,
      _["size"] = iter.size,
      _["weights"] = iter.weights,
      _["centers"] = iter.centers,
      _["cov"] = iter.cov,
      _["code"] = iter.code,
      _["posterior"] = iter.posterior);
  }
  
  template <> SEXP wrap(const params& pa) {
    return Rcpp::List::create(
      _["n"] = pa.n,
      _["p"] = pa.p,
      _["alpha"] = pa.alpha,
      _["trimm"] = pa.trimm,
      _["no_trim"] = pa.no_trim,
      _["k"] = pa.k,
      _["equal_weights"] = pa.equal_weights,
      _["zero_tol"] = pa.zero_tol);
      _["restrC"] = pa.restrC;
      _["deterC"] = pa.deterC;
      _["restr_fact"] = pa.restr_fact;
      _["cshape"] = pa.cshape;
      _["opt"] = pa.opt;
  }
}

/**
 * Convert an iteration object to Rcpp List
 *
 * @param iter: iteration object.
 *
 * @returns Rcpp List with same fields as the iteration object.
 */
Rcpp::List iter_to_list(iteration &iter)
{
  return Rcpp::List::create(
      _["centers"] = iter.centers,
      _["cov"] = iter.cov,
      _["cluster"] = iter.cluster,
      _["disttom"] = iter.disttom,
      _["obj"] = iter.obj,
      _["NlogL"] = iter.NlogL,
      _["size"] = iter.size,
      _["weights"] = iter.weights,
      _["code"] = iter.code,
      _["posterior"] = iter.posterior);
}

// Internal function to perform the E-step for a Gaussian mixture distribution
//  see the documentation of estepRR() below.
//
double estepXX(const arma::mat& ll, arma::mat& postprob, arma::vec& logpdf) {

    // Maximum log-likelihood along rows
    arma::vec maxll = arma::max(ll, 1);

    // Subtract maxll from each row to avoid underflow
    arma::mat post = arma::exp(ll.each_col() - maxll);

    // Compute density: row-wise sum of exponentiated values
    arma::vec density = arma::sum(post, 1);

    // Posterior probabilities: normalize each row
    postprob = post.each_col() / density;

    // Log PDF for each observation
    logpdf = arma::log(density) + maxll;

    // Compute objective (total log-likelihood)
    double obj = arma::sum(logpdf);

    return obj;
}

/**
 * Apply restrictions to eigenvalues.
 *
 * @param autovalues: p x k matrix containin eigenvalues.
 * @param ni_ini: current sample size of the clusters.
 * @param factor_e: the level of the constraints.
 * @param zero_tol: tolerance level.
 *
 * @returns matrix with constrained eigenvalues.
 */
arma::mat restr2Eigenv(arma::mat autovalues, arma::vec ni_ini, double factor_e, double zero_tol)
{

    // Initializations
    double c = factor_e;
    arma::mat d = autovalues.t();
    int p = autovalues.n_rows;
    int k = autovalues.n_cols;
    int n = arma::accu(ni_ini);
    arma::mat nis(k, p); 
    nis.each_col() = ni_ini;
    
    //VT::20.01.2024: set any negative values in d to 0
    arma::uvec idx = find(d < 0);
    d.elem(idx).fill(0);
  
    // d_ is the ordered set of values in which the restriction objective function change the definition
    // points in d_ correspond to the frontiers for the intervals in which this objective function has the same definition
    // ed is a set with the middle points of these intervals
    arma::vec d_ = arma::sort(arma::vectorise(arma::join_cols(d, d / c)));
    int dim = d_.n_elem;
    
    arma::vec d_1 = arma::join_cols(d_, arma::vec({d_.back() * 2}));
    arma::vec d_2 = arma::join_cols(arma::vec({0}), d_);
    arma::vec ed = (d_1 + d_2) / 2;
    
    dim++;
    
    // the only relevant eigenvalues are those belonging to a clusters with sample size greater than 0.
    // eigenvalues corresponding to a clusters with 0 individuals has no influence in the objective function
    // if all the eigenvalues are 0 during the smart initialization we assign to all the eigenvalues the value 1
    if(arma::max(d.elem(arma::find(nis > 0))) <= zero_tol) {
        return arma::mat(p, k, arma::fill::zeros); // solution corresponds to 0 matrix
    }
    
    // we check if the eigenvalues verify the restrictions
    if(std::abs(arma::max(d.elem(arma::find(nis > 0))) / arma::min(d.elem(arma::find(nis > 0)))) <= c) {
        d.elem(arma::find(nis == 0)).fill(arma::mean(d.elem(arma::find(nis > 0))));
        return d.t(); // the solution corresponds to the input because it verifies the constraints
    }
    
    arma::mat t(k, dim);
    arma::mat s(k, dim);
    arma::mat r(k, dim);
    
    arma::vec sol(dim);
    arma::vec sal(dim);

    for (int mp_ = 0; mp_ < dim; mp_++) {
        for (int i = 0; i < k; i++) {
            r(i, mp_) = arma::accu(d.row(i) < ed(mp_)) + arma::accu(d.row(i) > ed(mp_) * c);
            s(i, mp_) = arma::accu(d.row(i) % (d.row(i) < ed(mp_)));
            t(i, mp_) = arma::accu(d.row(i) % (d.row(i) > ed(mp_) * c));
        }
        
        sol(mp_) = arma::accu(ni_ini / n % (s.col(mp_) + t.col(mp_) / c)) / (arma::accu(ni_ini / n % (r.col(mp_))));
        
        arma::mat sol_mp_matrix(k, p, arma::fill::value(sol(mp_)));
        arma::mat c_sol_mp_matrix(k, p, arma::fill::value(c * sol(mp_)));
        
        arma::mat e = sol(mp_) * arma::conv_to<arma::mat>::from(d < sol(mp_)) +
                      d % arma::conv_to<arma::mat>::from((d >= sol(mp_)) % (d <= c * sol(mp_))) +
                      (c * sol(mp_)) * arma::conv_to<arma::mat>::from(d > c * sol(mp_));
        
        arma::mat o = -1.0 / 2.0 * nis / n % (arma::log(e) + d / e);
        
        sal(mp_) = arma::accu(o);
    }

    // m is the optimum value for the eigenvalues procedure
    int eo = sal.index_max();
    double m = sol(eo);
    
    // based on the m value we get the restricted eigenvalues
    return(m * arma::conv_to<arma::mat>::from(d < m) +
           d % arma::conv_to<arma::mat>::from(d >= m) % arma::conv_to<arma::mat>::from(d <= c * m) +
           (c * m) * arma::conv_to<arma::mat>::from(d > c * m)).t();
}

arma::mat HandleSmallEv(arma::mat autovalues, double zero_tol)
{

    // Initializations
    int p = autovalues.n_rows;
    // int k = autovalues.n_cols;       //  not used here

    // Replace all ev less than zero_tol by zero_tol    
    autovalues.elem(find(autovalues <= zero_tol)).fill(zero_tol);   // "<= zero.tol" for checking for zero
    
    arma::mat mi = arma::min(autovalues);                           // the minimum eigenvalue of each cluster
    arma::mat ma = arma::max(autovalues);                           // the maximum eigenvalue of each cluster
    
    //  (mi/ma).print("mi/ma");
    //  (mi/zero_tol).print("mi/zero_tol");
    arma::uvec idx_iter = find(mi/ma <= zero_tol);
    
    //  idx_iter.print("idx_iter");
    //  autovalues.print("autovalues");

    // for each of these clusters. set all "normal" eigenvalues to a high value
    for(arma::uword ki = 0; ki < idx_iter.n_elem; ki++) { 
        arma::uvec idx = arma::find(autovalues.cols(idx_iter(ki), idx_iter(ki)) > mi(idx_iter(ki))/zero_tol);
        //  idx.print("To change index");
        for(arma::uword ji = 0; ji < idx.n_elem; ji++)
            autovalues(idx(ji), idx_iter(ki)) = mi(idx_iter(ki))/zero_tol;
    }
    
    //  autovalues.print("Modified autovalues");
   
	arma::mat det = prod(autovalues);								// the new determinants
    //  det.print("The new determinants");

    //	autovalues_det contains the autovalues corrected by the determinant
    //	the product of the eigenvalues of each column in autovalues_det is equal to 1
    
    //  arma::diagmat(det).print("The diagonal of det");
    //  arma::diagmat(arma::pow(det, -1.0/p)).print("det^");

	arma::mat autovalues_det = autovalues * arma::diagmat(arma::pow(det, -1.0/p));

    //  autovalues_det.print("Returning autovalues corrected by the determinant");
    
    return autovalues_det;  
}

/**
 * Apply restrictions to determinants.
 *
 * @param autovalues: p x k matrix containin eigenvalues.
 * @param ni_ini: current sample size of the clusters.
 * @param factor_e: the level of the constraints.
 * @param zero_tol: tolerance level.
 *
 * @returns matrix with constrained eigenvalues.
 */
arma::mat restr2Deter_old(arma::mat autovalues, arma::vec ni_ini, double factor_e, double zero_tol)
{

    // Initializations
    double c = factor_e;
    int p = autovalues.n_rows;
    int k = autovalues.n_cols;
    //int n = arma::accu(ni_ini);
    //arma::mat nis(k, p);
    //nis.each_col() = ni_ini;
    
//    Rcout << "Entering restr2Deter" << std::endl;

    if(p == 1)
        return restr2Eigenv(autovalues, ni_ini, factor_e, zero_tol);
    
    arma::mat es = prod(autovalues, 0);
    if(arma::max(es(find(ni_ini))) <= zero_tol) {
        arma::mat aret(p, k, arma::fill::zeros);
        // aret.print("!!!! - restr2Deter() returning zeros!");
        return aret;
    }
    
    arma::mat d = es.t();
	arma::mat autovalues_det = HandleSmallEv(autovalues, zero_tol);
	
    double ma = max(d(find(ni_ini > 0)));
    double mi = min(d(find(ni_ini > 0)));
 
//    Rcout << "Max: " << ma << std::endl;
//    Rcout << "Min: " << mi << std::endl;
  
    arma::mat dfin; 		
	if(ma/mi <= c) {
//        Rcout << "ma/mi <= c " << std::endl;
//        Rcout << "Fill with " << mean(d.elem(find(ni_ini  >= 0))) << std::endl;

		d.elem(find(ni_ini <= 0)).fill(mean(d.elem(find(ni_ini  >= 0))));
		dfin = pow(d, 1.0/p);
	} else {
//        arma::pow(d, 1.0/p).print("Passing to restr2Eigenv");
//        Rcout << "pow(c, 1.0/p)" << pow(c, 1.0/p) << std::endl;
		
        dfin = restr2Eigenv(pow(d, 1.0/p).t(), ni_ini, pow(c, 1.0/p), zero_tol);
    }
    
//    dfin.print("This is dfin");
//    diagmat(dfin).print("This is the diag matrix of dfin");
    
//    Rcout << "Exiting restr2Deter" << std::endl;
	return autovalues_det * diagmat(dfin);      							
}

/**
 * Apply restrictions to determinants.
 *
 * @param autovalues: p x k matrix containin eigenvalues.
 * @param ni_ini: current sample size of the clusters.
 * @param restr_factor: constraint level for the determinant constraints
 * @param cshape: constraint level for the eigenvalues constraints, default is 1e10
 * @param zero_tol: tolerance level.
 *
 * @returns matrix with constrained eigenvalues.
 */
arma::mat restr2Deter(arma::mat autovalues, arma::vec ni_ini, double restr_factor, double cshape, double zero_tol)
{

    // Initializations
    int p = autovalues.n_rows;
    int K = autovalues.n_cols;
    arma::vec ni_ini_1 = {1};
    
    //  Rcout << "Entering restr2Deter2" << std::endl;

    arma::uvec idx = find(autovalues < 1e-16);
    autovalues.elem(idx).fill(0);
    arma::mat autovalues_ = autovalues; 
    
    if(p == 1) return restr2Eigenv(autovalues, ni_ini, restr_factor, zero_tol);
    for(int k = 0; k < K; k++)  {
        autovalues_.col(k) = restr2Eigenv(autovalues.col(k), ni_ini_1, cshape, zero_tol);
    }

    arma::mat es = prod(autovalues_, 0);                        // apply(autovalues_, 2, prod)
   
    idx = find(es == 0);                                        // es[es==0] <- 1
    es.elem(idx).fill(1);                                       // ...
    arma::mat esmat(arma::pow(es, 1.0/p));
    esmat.set_size(1, K);
    arma::mat gm = autovalues_/repmat(esmat, p, 1);
    arma::mat d = sum(autovalues/gm, 0)/p;
    
    idx = find_nan(d);                                           // d[is.nan(d)] <- 0
    d.elem(idx).fill(0);                                         // ...
    
    arma::mat dfin = restr2Eigenv(d, ni_ini, pow(restr_factor, 1.0/p), zero_tol);
    
    // gm = gm * (gm > 0) + 1 * (gm == 0)
    idx = find(gm < 0);                                         // gm[gm<0] <- 0
    gm.elem(idx).fill(0);                                       // ...
    idx = find(gm == 0);                                        // gm[gm==0] <- 1
    gm.elem(idx).fill(1);                                       // ...

    dfin.set_size(1, K);
    dfin = repmat(dfin, p, 1);
    
	return dfin % gm;
}

/**
 * Manages constraints of a clustering stage.
 *
 * @param iter: a reference to the cluster information. Its values are modified.
 * @param pa: a reference to the procedure parameters.
 */
void fRestr(iteration &iter, params &pa, GPCMPars &pars) {

    if(pa.restrC == 0) {
        arma::cube u(pa.p, pa.p, pa.k); // Eigenvectors
        arma::mat d(pa.p, pa.k);        // Eigenvalues
    
        for(int ki = 0; ki < pa.k; ki++) {
    
            arma::vec eigval;
            arma::mat eigvec;
            arma::eig_sym(eigval, eigvec, iter.cov.slice(ki));
        
            u.slice(ki) = eigvec;
            d.col(ki) = eigval;
        }
    
        d.elem(find(d < 0)).fill(0);  // all eigenvalue < 0 are assigned to 0, this issue appears for numerical errors
        
        if(!pa.deterC)
            d = restr2Eigenv(d, iter.size, pa.restr_fact, pa.zero_tol);
        else {
            //  iter.size.print("Autovalues size");
            //  d.print("Autovalues to be restricted");
            d = restr2Deter(d, iter.size, pa.restr_fact, pa.cshape, pa.zero_tol);
        }
        
        // checking for singularity in all clusters.
        int code = d.max() > pa.zero_tol;
        iter.code = code;
        
        if(code)  {
            // Reconstructing the cov matrices
            // d.print();
            for (int ki = 0; ki < pa.k; ki++) { 
                iter.cov.slice(ki) = u.slice(ki) * arma::diagmat(d.col(ki)) * u.slice(ki).t();
            }
        }
    }else {
    
        // VT:::23.07.2026: Check if one or more of the covariance matrices are sigular.
        // If so, skip the call to restrSigmaGPCM() and return code=0
        arma::vec cond(pa.k);
        for(int ki = 0; ki < pa.k; ki++) {
            cond(ki) = arma::rcond(iter.cov.slice(ki));
        }
        bool code = cond.min() >= arma::datum::eps;
        
        // DEBUG
        // Rcout << "In fRestr(): cond=" << cond << std::endl;
        // cond.print("cond");
          
        if(code) {
            Rcpp::List tmp = restrSigmaGPCM(iter.cov, iter.size, pars);
            iter.cov = Rcpp::as<arma::cube>(tmp["Sigma"]);
            iter.OMG = Rcpp::as<arma::cube>(tmp["OMG"]);
            iter.lmd = Rcpp::as<arma::vec>(tmp["lmd"]);
            iter.code = code;
        }
    }    
}

/**
 * Calculates the initial cluster assignment and initial values for the parameters.
 *
 * @param x: matrix with observations and features.
 * @param iter: a reference to the cluster information. Its values are modified.
 * @param pa: a reference to the procedure parameters.
 */
void initClusters(arma::mat x, iteration &iter, params &pa)
{
    arma::mat iter_center = arma::mat(pa.p, pa.k);        // Column ki stores the centers of cluster p
    arma::cube iter_sigma = arma::cube(pa.p, pa.p, pa.k); // Covariance matrix of each cluster
    arma::vec size;                                       // Cluster sizes
    arma::vec weights;                                    // Cluster weigths
    
    arma::ivec idx = arma::randi(pa.k * (pa.p + 1), arma::distr_param(0, pa.n - 1));
    
    for (int ki = 0; ki < pa.k; ki++) {
        // Select the p+1 next elements of idx
        arma::uvec rows_to_select = arma::conv_to<arma::uvec>::from(idx.rows(ki * (pa.p + 1), ki * (pa.p + 1) + pa.p));
        
        // Set cluster centers to random observation
        iter_center.col(ki) = x.row(rows_to_select(0)).t();
        
        arma::mat X_ini = x.rows(rows_to_select);
        
        // Calculate cov matrix of cluster k
        iter_sigma.slice(ki) = pa.p / (float)(pa.p + 1) * arma::cov(X_ini);
    }
    
    iter.centers = iter_center.t();
    iter.cov = iter_sigma;
    
    if (pa.equal_weights) {
        size = arma::vec(pa.k, arma::fill::value(pa.no_trim / pa.k));
        weights = arma::vec(pa.k, arma::fill::value(1 / (double)pa.k));
    } else {
        weights = runif(pa.k);
        weights /= arma::accu(weights);
        size = arma::round(pa.n * weights);
    }
    
    iter.size = size;
    iter.weights = weights;
}

/**
 * FUNCTIONS TO CALCULATE THE DENSITY OF A MULTIVARIATE NORMAL
 *
 * Taken from: Nino Hardt, Dicko Ahmadou, Benjamin Christoffersen. 
 *  Faster Multivariate Normal densities with RcppArmadillo and OpenMP
 *  https://gallery.rcpp.org/articles/dmvnorm_arma/
 */

static double const log2pi = std::log(2.0 * M_PI);

/* C++ version of the dtrmv BLAS function */
void inplace_tri_mat_mult(arma::rowvec &x, arma::mat const &trimat) {
  arma::uword const n = trimat.n_cols;

  for(unsigned j = n; j-- > 0;) {
    double tmp(0.);
    for (unsigned i = 0; i <= j; ++i)
      tmp += trimat.at(i, j) * x[i];
    x[j] = tmp;
  }
}

/*
//  This version of dmvnrm_arma_fast() will crash in chol() if one of the covariances 
//  is slightly non-positive definite.
//  If this happens, we will replaace chol() by eigendecomposition.
//  Alternatively, SVD can be used - more stable than eigen but also slower.

arma::vec dmvnrm_arma_fast(arma::mat const &x, arma::rowvec const &mean, arma::mat const &cov, 
        bool const logd=false) {
        
    using arma::uword;
    uword const n = x.n_rows, xdim = x.n_cols;
    arma::vec out(n);
    arma::mat const rooti = arma::inv(trimatu(arma::chol(cov)));
    double const rootisum = arma::sum(log(rooti.diag())),
               constants = -(double)xdim / 2.0 * log2pi,
               other_terms = rootisum + constants;
    
    arma::rowvec z;
    for(uword i = 0; i < n; i++) {
        z = (x.row(i) - mean);
        inplace_tri_mat_mult(z, rooti);
        out(i) = other_terms - 0.5 * arma::dot(z, z);
    }
    
    if(logd)
        return out;
    
    return exp(out);
}
*/

arma::vec dmvnrm_arma_fast(const arma::mat &x,
                           const arma::rowvec &mean,
                           const arma::mat &cov,
                           const bool logd = false)
{
    using arma::uword;
    const uword n = x.n_rows;
    const uword xdim = x.n_cols;

    arma::vec out(n);

    arma::mat U;
    arma::mat rooti;
    double rootisum;
    bool chol_ok = arma::chol(U, cov);
    
    if (chol_ok) {
        rooti = arma::inv(arma::trimatu(U));
        rootisum = arma::sum(log(rooti.diag()));
    } else {
        
        /*
            // Fall back to SVD
            arma::mat U, V;
            arma::vec s;
            arma::svd(U, s, V, cov);
            
            double tol = std::max(s.max(), 1.0) * 1e-12;
            s.transform([tol](double x){ return std::max(x, tol); });
            
            rooti = U * arma::diagmat(1.0 / arma::sqrt(s));
            rootisum = -0.5 * arma::sum(log(s));
        */
        
        // Fall back to eigendecomposition
        arma::vec eigval;
        arma::mat eigvec;
    
        arma::eig_sym(eigval, eigvec, cov);
    
        // Tolerance for tiny negative eigenvalues
        double tol = std::max(eigval.max(), 1.0) * 1e-12;
    
        eigval.for_each([tol](double &v) {
            if (v < tol) v = tol;
        });
    
        arma::vec invsqrt = 1.0 / arma::sqrt(eigval);
    
        rooti = eigvec * arma::diagmat(invsqrt);
    
        // log(det(rooti)) = -1/2 log(det(cov))
        rootisum = -0.5 * arma::sum(log(eigval));
    } 
    
    const double constants =
        -(double)xdim / 2.0 * log2pi;

    const double other_terms = rootisum + constants;

    arma::rowvec z;

    for (uword i = 0; i < n; ++i) {

        z = x.row(i) - mean;

        if (chol_ok) {
            inplace_tri_mat_mult(z, rooti);
        } else {
            z *= rooti;
        }

        out(i) = other_terms - 0.5 * arma::dot(z, z);
    }

    return logd ? out : arma::exp(out);
} 

/**
 * Calculate the objective function value.
 *
 * @param x: matrix with observations and features.
 * @param iter: a reference to the cluster information. Its values are modified.
 * @param pa: a reference to the procedure parameters.
 */
void calcObj(arma::mat x, iteration &iter, params &pa)
{
    int n = pa.n;
    int k = pa.k;
    Rcpp::String opt = pa.opt;

    // VT::25.09.2024 - Compute always the classification log-likelihood and store NlogL
    //  which will be used to compute the BIC CLACLA and MIXCLA.

    // VT::11.02.2025 - use the logs of the likelihood
/*    
    arma::vec ww(n);
    arma::vec w;
    arma::vec ww1(n);
    arma::vec w1;

    //  Rcout << "Entering calcObj()" << std::endl;

    for(int ki = 0; ki < k; ki++) {
        w = iter.weights(ki) * dmvnrm_arma_fast(x, iter.centers.row(ki), iter.cov.slice(ki));
        
        //	calculates each individual contribution for the obj function mixture
        ww = w % (w >= 0) + ww;
        
        w1 = w % arma::conv_to<arma::mat>::from(iter.cluster == ki + 1);

        //	calculates each individual contribution for the obj function hard
        ww1 = w1 % (w1 >= 0) + ww1;
    }
    iter.NlogL = -2 * arma::accu(arma::log(ww1.elem(arma::find(iter.cluster > 0))));

    if(opt == "HARD") 
        iter.obj = -1 * iter.NlogL / 2;
    else
        iter.obj = arma::accu(arma::log(ww.elem(arma::find(iter.cluster > 0))));
*/        
        
    arma::vec wx2(n);
    arma::mat ll(n, k);
    for(int ki = 0; ki < k; ki++)   {  
        ll.col(ki) = log(iter.weights(ki)) + dmvnrm_arma_fast(x, iter.centers.row(ki), iter.cov.slice(ki), true);
            
        if(iter.size(ki) > pa.zero_tol) {

            //	calculates each individual contribution for the obj function hard
            wx2 += ll.col(ki) % arma::conv_to<arma::mat>::from(iter.cluster == ki+1);
       }
    }
    iter.NlogL = -2 * arma::accu(wx2.elem(arma::find(iter.cluster > 0)));
    
    if(opt == "HARD") 
        iter.obj = -1 * iter.NlogL / 2;
    else {
        arma::mat postprob;
        arma::vec logpdf;
        iter.obj = estepXX(ll, postprob, logpdf);
    }
}

/**
 * Estimate the model parameters.
 *
 * @param x: matrix with observations and features.
 * @param iter: a reference to the cluster information. Its values are modified.
 * @param pa: a reference to the procedure parameters.
 */
void estimClustPar(arma::mat x, iteration &iter, params &pa)
{
    for(int ki = 0; ki < pa.k; ki++) {
        if(iter.size(ki) > pa.zero_tol) {
            iter.centers.row(ki) = (iter.posterior.col(ki).t() * x) / iter.size(ki);
            // x centered
            arma::mat X_c = x;
            X_c.each_row() -= iter.centers.row(ki);
            X_c.each_col() %= iter.posterior.col(ki);
            iter.cov.slice(ki) = (X_c.t() * X_c) / iter.size(ki);
        } else {
            iter.centers.row(ki) = arma::mat(1, pa.p);
            iter.cov.slice(ki) = arma::mat(pa.p, pa.p, arma::fill::eye);
        }
    }
}

/**
 * Find cluster assignment and trimming.
 *
 * @param x: matrix with observations and features.
 * @param iter: a reference to the cluster information. Its values are modified.
 * @param pa: a reference to the procedure parameters.
 */
void findClustAssig(arma::mat x, iteration &iter, params &pa)
{
    int k = pa.k;
    int n = pa.n;
    int no_trim = pa.no_trim;
    bool equal_weights = pa.equal_weights;
    Rcpp::String opt = pa.opt;
    
    //Rcout << "Entering findClustAssig()" << std::endl;
    //iter.size.print("iter.size");
    
    arma::mat ll(n, k);
    for(int ki = 0; ki < k; ki++)   {        
        
        // VT::11.02.2025 - use the logs of the likelihood
        // ll.col(ki) = iter.weights(ki) + dmvnrm_arma_fast(x, iter.centers.row(ki), iter.cov.slice(ki), false);
        
        // iter.cov.slice(ki).print();
        
        ll.col(ki) = log(iter.weights(ki)) + dmvnrm_arma_fast(x, iter.centers.row(ki), iter.cov.slice(ki), true);
    }

    arma::uvec old_assig = iter.cluster;
    arma::vec pre_z;
    arma::uvec tc_set;
    arma::mat postprob;
    arma::vec logpdf;
    estepXX(ll, postprob, logpdf);
    
    if (opt == "HARD") {
        pre_z = arma::max(ll, 1);
    } else {
        //  pre_z = arma::sum(ll, 1);
        pre_z = logpdf;
    }
    
    // Rcout << "Before trimming ..." << std::endl;

    // Determine elements to trim
    arma::uvec sorted_index = arma::sort_index(pre_z, "descending");
    arma::uvec last_indexes = arma::linspace<arma::uvec>(no_trim, n - 1, n - no_trim);
    arma::uvec obs_to_trim = sorted_index.elem(last_indexes);
    arma::uvec ind = (arma::index_max(ll, 1) + 1);
    ind.elem(obs_to_trim) = arma::zeros<arma::uvec>(n - no_trim);
    tc_set = ind != 0;

    // Cluster assignment with trimming
    // Rcout << "New cluster assignment ..." << std::endl;
    iter.cluster = (arma::index_max(ll, 1) + 1) % tc_set;

    // Calculate posterior matrix
    if(opt == "MIXT") {
        //VT::11.02.2025 - use the logs of the likelihood
        //  iter.posterior = ll;
        //  iter.posterior.each_col() /= (pre_z + (pre_z == 0));

        iter.posterior = postprob;

    } else {
    
        arma::uvec one_to_n = arma::linspace<arma::uvec>(0, n - 1, n);
        arma::uvec aux_assig = iter.cluster - 1 + (iter.cluster == 0);
        
        // 2xn matrix containing (observation, cluster) pairs in each column
        arma::umat subscripts = (arma::join_rows(one_to_n, aux_assig)).t();
        iter.posterior = arma::mat(n, k);
        iter.posterior.elem(arma::sub2ind(arma::size(iter.posterior), subscripts)) = arma::ones<arma::vec>(n);
    }

    iter.posterior.each_col() %= arma::conv_to<arma::vec>::from(tc_set); // Set to 0 all trimmed rows
    
    if(opt == "HARD" && arma::all(old_assig == iter.cluster)) {
        iter.code = 2;
    }

    // Obtain the clusters size
    iter.size = (arma::sum(iter.posterior, 0)).t();
    
    if(!equal_weights)
        iter.weights = iter.size/no_trim;

}

/**
 * Apply niter concentration steps to initial solution given by iter and pa.
 *
 * @param niter: number of concentration steps.
 * @param x: matrix with observations and features.
 * @param iter: a reference to the cluster information. Its values are modified.
 * @param pa: a reference to the procedure parameters.
 */
void concentration_steps(int niter, arma::mat x, iteration &iter, params &pa, GPCMPars &pars)
{
    // Rcout << "Enter concentration_steps ..." << std::endl;

    for(int i1 = 0; i1 < niter; i1++) {

    // Rcout << i1 << "   " << "Before fRestr() ..." << std::endl;
    // iter.cov.print("cov before");

        fRestr(iter, pa, pars); // restricting the clusters' scatter structure (Changes the iter object)
    
    // Rcout << i1 << "   " << "After fRestr(): iter.code=" << iter.code << std::endl;
    // iter.cov.print("cov after");
                
        if(iter.code == 0)  {
            if(i1 > 0) {
                calcObj(x, iter, pa);
                Rcpp::warning("Data not in general position");
            } else {
                arma::cube cov = arma::cube(pa.p, pa.p, pa.k);
                cov.each_slice() = arma::eye(pa.p, pa.p);
                iter.cov = cov;
            }
        }
        
        // Estimate the cluster's assigment and TRIMMING (mixture models and HARD)
        // Rcout << "Before findClustAssig():" << std::endl; 
        findClustAssig(x, iter, pa); 
        
        // Rcout << "After findClustAssig(): iter.code=" << iter.code << ", i1=" << i1 << std::endl; 
        
        if((int)iter.code == 2 || (i1 == niter - 1))
            break;
        
        // Rcout << "Estimate Cluster Par:" << std::endl; 
        estimClustPar(x, iter, pa); // estimates the cluster's parameters
    }
    
    // Rcout << "Before calcObj():" << std::endl; 
    calcObj(x, iter, pa); // calculates the objcetive function value
    // Rcout << "After calcObj():" << std::endl; 
}

// Internal function for concentration steps (refinement) in tclust()
// @name tclust_c2
// @param x Rcpp::NumericMatrix, The input data.
// @param k The number of clusters initially searched for.
// @param cluster arma::uvec A numerical vector of size \code{n} containing the 
//      cluster assignment for each observation. Cluster names are integer numbers 
//      from 1 to k, 0 indicates trimmed observations. Note that it could be empty 
//      clusters with no observations when \code{equal_weights=FALSE}.
// @param alpha double, The proportion of observations to be trimmed.
// @param restrC Restriction type (0 for restriction on eigenvalues or determinant)
// @param deterC Restriction on determinants (true or false)
// @param restr_fact The constant \code{restr_fact >= 1} constrains the allowed 
//  differences among group scatters in terms of eigenvalues ratio. Larger values 
//  imply larger differences of group scatters, a value of 1 specifies the 
//  strongest restriction.
// @param cshape Shape constraint
// @param niter2 The maximum number of concentration steps to be performed for the 
//  \code{nkeep} solutions kept for further iteration. The concentration steps are 
//  stopped, whenever two consecutive steps lead to the same data partition.
// @param opt Define the target function to be optimized. A classification likelihood 
//  target function is considered if \code{opt="HARD"} and a mixture classification 
//  likelihood if \code{opt="MIXT"}.
// @param equal_weights A logical value, specifying whether equal cluster weights 
//  shall be considered in the concentration and assignment steps.
// @param zero_tol The zero tolerance used. By default set to 1e-16.
// @export
// [[Rcpp::export]]
iteration tclust_c2(arma::mat x, int k, arma::uvec cluster, double alpha = 0.05,
                     int restrC=0, bool deterC=false, double restr_fact = 12, double cshape=1e10, Nullable<List> paList_ = R_NilValue,
                     int niter2 = 20, Rcpp::String opt = "HARD",
                     bool equal_weights = false, double zero_tol = 1e-16)
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
    
    GPCMPars pars;
    if(paList_.isNotNull()){
        Rcpp::List paList(paList_);    // casting to underlying type List
        pars = parsePars(paList);
    }

    params pa;
    pa.n = n;
    pa.p = p;
    pa.alpha = alpha;
    pa.no_trim = no_trim;
    pa.trimm = n - no_trim;
    pa.k = k;
    pa.equal_weights = equal_weights;
    pa.zero_tol = zero_tol;
    pa.restrC = restrC;
    pa.deterC = deterC;
    pa.restr_fact = restr_fact;
    pa.cshape = cshape;
    pa.opt = opt;
    
    iteration iter;
    iter.centers = arma::mat(k, p);
    iter.cov = arma::cube(p, p, k);
    iter.cluster = cluster;
    iter.obj = 0.0,
    iter.NlogL = 0.0,
    iter.size = size;
    iter.lmd = arma::vec(k);
    iter.OMG = arma::cube(p, p, k);

    // DEBUG
    // Rcout << "Enter tclust_c2()=================================================" << std::endl;

    
    // VT::29.10.2024 - the results with or without equal_weights were identical
    //  because the condition 'if(!equal_weights)' was missing
    iter.weights = arma::ones<arma::vec>(size.n_elem) / no_trim;
    if(!equal_weights)
        iter.weights = size / no_trim;
        
    iter.code = 0;
    iter.posterior = posterior;
    
    estimClustPar(x, iter, pa);
    
    // VT::06.05.2024 - If niter2==0 (no iterations in the second step) tclust_c2 could
    //  break with "chol(): decomposition failed" - this happens if nkeep = nstart, i.e.
    // call tclust_c2 on all tried solutions (nstart) and at least one of the solutions
    // has a class with one member which means that its covariance is the zero matrix.
    //  To solve this, we do an exctra call to frestr(), if niter2 == 0.

    if(niter2 == 0) {
        fRestr(iter, pa, pars); 
    }
    
    concentration_steps(niter2, x, iter, pa, pars);

    // DEBUG
    // Rcout << "Returning from tclust_c2():" << std::endl; 
    
    return iter;
}

// Internal function for concentration steps (initializations) in tclust()
// @name tclust_c1
// @param x Rcpp::NumericMatrix, The input data.
// @param k The number of clusters initially searched for.
// @param alpha double, The proportion of observations to be trimmed.
// @param restrC Restriction type (0 for restriction on eigenvalues or determinant)
// @param deterC Restriction on determinants (true or false)
// @param restr_fact The constant \code{restr_fact >= 1} constrains the allowed 
//  differences among group scatters in terms of eigenvalues ratio. Larger values 
//  imply larger differences of group scatters, a value of 1 specifies the 
//  strongest restriction.
// @param cshape Shape constraint
// @param niter1 int, The number of concentration steps to be performed for the 
//     nstart initializations. 
// @param opt Define the target function to be optimized. A classification likelihood 
//  target function is considered if \code{opt="HARD"} and a mixture classification 
//  likelihood if \code{opt="MIXT"}.
// @param equal_weights A logical value, specifying whether equal cluster weights 
//  shall be considered in the concentration and assignment steps.
// @param zero_tol The zero tolerance used. By default set to 1e-16.
// @export
// [[Rcpp::export]]
Rcpp::List tclust_c1(arma::mat x, int k, double alpha = 0.05,
                     int restrC=0, bool deterC=false, double restr_fact = 12, double cshape=1e10, Nullable<List> paList_ = R_NilValue,
                     int niter1 = 3, Rcpp::String opt = "HARD",
                     bool equal_weights = false, double zero_tol = 1e-16)
{

  int n = x.n_rows;
  int p = x.n_cols;
  int no_trim = std::floor(n * (1 - alpha));

    GPCMPars pars;
    if(paList_.isNotNull()){
        Rcpp::List paList(paList_);    // casting to underlying type List
        pars = parsePars(paList);
    }
    
  params pa;
  pa.n = n;
  pa.p = p;
  pa.alpha = alpha;
  pa.no_trim = no_trim;
  pa.trimm = n - no_trim;
  pa.k = k;
  pa.equal_weights = equal_weights;
  pa.zero_tol = zero_tol;
  pa.restrC = restrC;
  pa.deterC = deterC;
  pa.restr_fact = restr_fact;
  pa.cshape = cshape;
  pa.opt = opt;

  iteration iter;
  iter.obj = 0.0,
  iter.NlogL = 0.0,
  iter.cluster = arma::uvec(n);
  iter.size = arma::vec(k);
  iter.weights = arma::vec(k);
  iter.centers = arma::mat(k, p);
  iter.cov = arma::cube(p, p, k);
  iter.code = 0;
  iter.posterior = arma::mat(n, p);
    iter.lmd = arma::vec(k);
    iter.OMG = arma::cube(p, p, k);

    // DEBUG
    // Rcout << "Enter tclust_c1()=================================================" << std::endl;

  initClusters(x, iter, pa);                        // Cluster random initialization
  concentration_steps(niter1, x, iter, pa, pars);   // Apply niter1 concentration steps

    // DEBUG
    // Rcout << "Returning from tclust_c1():" << std::endl; 

  return Rcpp::List::create(
      _["obj"] = iter.obj,
      _["cluster"] = iter.cluster);
}


// [[Rcpp::export]]
arma::mat tclust_restr2_eigenv(arma::mat autovalues, arma::vec ni_ini, double factor_e=12, double zero_tol=1e-16) {
    return restr2Eigenv(autovalues, ni_ini, factor_e, zero_tol);
}

// [[Rcpp::export]]
arma::mat tclust_restr2_deter_old(arma::mat autovalues, arma::vec ni_ini, double factor_e=12, double zero_tol=1e-16) {
    return restr2Deter_old(autovalues, ni_ini, factor_e, zero_tol);
}

// [[Rcpp::export]]
arma::mat tclust_restr2_deter(arma::mat autovalues, arma::vec ni_ini, double restr_factor=12, double cshape=1e10, double zero_tol=1e-16) {
    return restr2Deter(autovalues, ni_ini, restr_factor, cshape, zero_tol);
}

// [[Rcpp::export]]
arma::mat tclust_HandleSmallEv(arma::mat autovalues, double zero_tol=1e-16) {
    return HandleSmallEv(autovalues, zero_tol);
}

// [[Rcpp::export]]
arma::vec dmvnrm(arma::mat x, arma::rowvec mean, arma::mat cov) {
    return dmvnrm_arma_fast(x, mean, cov, false);
}

//' Function to perform the E-step for a Gaussian mixture distribution
//' @name estepRR
//' @description Compute the log PDF for each observation, the posterior probabilities
//'     and the objective function (total log-likelihood) for a Gaussian 
//'     mixture distribution
//' @param ll Rcpp::NumericMatrix, n-by-k where \code{n} is the number of 
//'  observations and \code{k} is the number of clusters.
//' @return The function returns a list with the following elements:
//' \itemize{
//'    \item obj The value of the objective function (total log-likelihood)
//'    \item postprob an \code{n-by-k} matrix with the posterior probablilities
//'    \item logpdf a vector of length \code{n} containing the log PDF for 
//'         each observation
//' }
//' @details Formally a mixture model corresponds to the mixture distribution that
//'   represents the probability distribution of observations in the overall
//'   population. Mixture models are used
//'   to make statistical inferences about the properties of the
//'   sub-populations given only observations on the pooled population, without
//'   sub-population-identity information.
//'   Mixture modeling approaches assume that data at hand $\eqn{y_1, ..., y_n} in
//'   \eqn{R^p} come from a probability distribution with density given by the sum of k components
//'   \deqn{\sum_{j=1}^k  \pi_j \phi( \cdot, \theta_j)}
//'   with \eqn{\phi( \cdot, \theta_j)} being the
//'   \eqn{p}-variate  (generally multivariate normal) densities with parameters
//'   \eqn{\theta_j}, \eqn{j=1, \ldots, k}. Generally \eqn{\theta_j= (\mu_j, \Sigma_j)}
//'   where \eqn{\mu_j} is the population mean  and   \eqn{\Sigma_j} is the covariance
//'   matrix for component \eqn{j}.
//'   \eqn{\pi_j} is the (prior) probability of component \eqn{j}.
//'   The objective function is obj is equal to
//'    \deqn{ obj = \log   \left( \prod_{i=1}^n  \sum_{j=1}^k \pi_j \phi (y_i; \; \theta_j)    \right) 
//'     }
//'
//'    or
//'  
//'     \deqn{ obj =  \sum_{i=1}^n  \log   \left( \sum_{j=1}^k \pi_j \phi (y_i; \; \theta_j)    \right)
//'     }
//'     where \eqn{k} is the number of components of the mixture and \eqn{\pi_j} are the 
//'     component probabilitites and \eqn{\theta_j} are the  parameters of the \eqn{j}-th 
//'     mixture component.
//'
//'
//' @references
//'
//'   McLachlan, G.J.; Peel, D. (2000). Finite Mixture Models. Wiley. ISBN 0-471-00626-2
//'
//' @examples
//'##      Generate two Gaussian normal distributions
//'##      and do not produce plots
//'
//'        mu1 = c(1,2)
//'        sigma1 = matrix(c(2, 0, 0, .5), nrow=2, byrow=TRUE)    #[2 0; 0 .5];
//'        mu2 = c(-3, -5)
//'        sigma2 = matrix(c(1, 0, 0, 1), nrow=2, byrow=TRUE)
//'        n1 = 100
//'        n2 = 200
//'        Y = rbind(MASS::mvrnorm(n1, mu1, sigma1), 
//'                  MASS::mvrnorm(n2, mu2, sigma2))
//'        k = 2
//'        pi = c(1/3, 2/3)
//'        mu = rbind(mu1, mu2)
//'        sigma = array(0, dim=c(2,2,2))
//'        sigma[,,1] = sigma1
//'        sigma[,,2] = sigma2
//'        
//'        ll = matrix(0, nrow=n1+n2, ncol=2)
//'        for(j in 1:k)
//'            ll[,j] = log(pi[j]) +  tclust:::dmvnrm(Y, mu[j,], sigma[,,j])
//'
//'        dd = tclust:::estepRR(ll)
//'        dd$obj
//'        dd$logpdf
//'        dd$postprob
        
// [[Rcpp::export]]
Rcpp::List estepRR(const arma::mat& ll) {

    arma::mat postprob;
    arma::vec logpdf;
    
    double obj = estepXX(ll, postprob, logpdf);

    return Rcpp::List::create(
        Rcpp::Named("obj") = obj,
        Rcpp::Named("postprob") = postprob,
        Rcpp::Named("logpdf") = logpdf);
}

