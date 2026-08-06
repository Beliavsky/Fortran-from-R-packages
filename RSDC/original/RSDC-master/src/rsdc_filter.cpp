// Fast Hamilton-filter log-likelihood for the RSDC models (added in 1.4-0).
// Mirrors the R reference rsdc_hamilton()/rsdc_likelihood() exactly (same ridge,
// same PD check, same log-sum-exp filter) but runs the hot inner loop in C++.
//
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace arma;

// Build one time-t transition matrix from covariates (TVTP), matching the R code:
//  - N == 2: logistic stay-probability on the diagonal, off-diagonal by complement;
//  - N >= 3: softmax over (N-1) free logit vectors per row, N-th logit fixed at 0.
static mat tvtp_P(const rowvec& xt, const mat& beta, int N, int p) {
  mat P(N, N, fill::zeros);
  for (int i = 0; i < N; ++i) {
    if (N == 2) {
      double lin = dot(xt, beta.row(i));
      double pii = 1.0 / (1.0 + std::exp(-lin));   // plogis
      P(i, i) = pii;
      for (int j = 0; j < N; ++j) if (j != i) P(i, j) = (1.0 - pii) / (N - 1);
    } else {
      vec logits(N, fill::zeros);
      for (int j = 0; j < N - 1; ++j)
        logits(j) = dot(xt, beta.row(i).subvec(j * p, (j + 1) * p - 1));
      logits(N - 1) = 0.0;
      double c = logits.max();
      vec el = exp(logits - c);
      P.row(i) = (el / accu(el)).t();
    }
  }
  return P;
}

// Full Hamilton filter + smoother. Returns filtered/smoothed probabilities, the
// log-likelihood and the per-observation contributions. Mirrors the pure-R
// rsdc_hamilton() engine = "r" path exactly (same ridge, PD check, log-sum-exp,
// smoothing recursion). On a non-PD regime or a degenerate step it returns
// ok = FALSE / log_likelihood = -Inf, which the R wrapper maps to NULL outputs.
// [[Rcpp::export]]
Rcpp::List rsdc_filter_full_cpp(const arma::mat& y, const arma::cube& sigma,
                                int mode, const arma::mat& P0,
                                const arma::mat& X, const arma::mat& beta,
                                const arma::vec& xi_init) {
  const int T = y.n_rows, K = y.n_cols, N = sigma.n_slices;
  const double eps = std::numeric_limits<double>::epsilon();   // .Machine$double.eps
  Rcpp::List fail = Rcpp::List::create(Rcpp::Named("ok") = false,
                                       Rcpp::Named("log_likelihood") = -arma::datum::inf);

  arma::mat logdens(N, T);
  const double cst = K * std::log(2.0 * M_PI);
  for (int m = 0; m < N; ++m) {
    arma::mat Rm = sigma.slice(m);
    arma::vec ev;
    if (!eig_sym(ev, Rm) || ev.min() < 1e-8) return fail;
    arma::mat Rr = Rm + 1e-8 * eye(K, K);
    arma::mat invR;
    if (!inv(invR, Rr)) return fail;
    double ld, sgn; log_det(ld, sgn, Rr);
    arma::vec quad = sum((y * invR) % y, 1);
    logdens.row(m) = (-0.5 * (ld + quad + cst)).t();
  }

  const int p = (mode == 1) ? (int) X.n_cols : 0;
  arma::mat filtered(N, T), predicted(N, T), smoothed(N, T);
  arma::vec loglik_t(T);
  arma::vec xi = (xi_init.n_elem == (arma::uword) N) ? xi_init / accu(xi_init)
                                                     : arma::vec(N).fill(1.0 / N);
  double loglik = 0.0;
  for (int t = 0; t < T; ++t) {
    arma::mat Pt = (mode == 0) ? P0 : tvtp_P(X.row(t), beta, N, p);
    arma::vec pred = Pt.t() * xi;
    predicted.col(t) = pred;
    arma::vec ld = logdens.col(t);
    double c = ld.max();
    arma::vec w = pred % exp(ld - c);
    double sw = accu(w);
    if (!std::isfinite(sw) || sw <= 0.0) return fail;
    filtered.col(t) = w / sw;
    loglik_t(t) = std::log(sw) + c;
    loglik += loglik_t(t);
    xi = filtered.col(t);
  }

  smoothed.col(T - 1) = filtered.col(T - 1);
  for (int t = T - 2; t >= 0; --t) {
    arma::mat Pt1 = (mode == 0) ? P0 : tvtp_P(X.row(t + 1), beta, N, p);
    arma::vec denom = predicted.col(t + 1);
    denom.transform([eps](double v) { return v < eps ? eps : v; });
    arma::vec ratio = smoothed.col(t + 1) / denom;
    arma::vec temp = filtered.col(t) % (Pt1 * ratio);
    double s = accu(temp);
    if (!std::isfinite(s) || s < eps) return fail;
    smoothed.col(t) = temp / s;
  }

  return Rcpp::List::create(
    Rcpp::Named("filtered_probs") = filtered,
    Rcpp::Named("smoothed_probs") = smoothed,
    Rcpp::Named("log_likelihood") = loglik,
    Rcpp::Named("loglik_t")       = loglik_t,
    Rcpp::Named("ok")             = true);
}

// [[Rcpp::export]]
double rsdc_loglik_cpp(const arma::mat& y,        // T x K observations
                       const arma::cube& sigma,    // K x K x N regime correlation matrices
                       int mode,                   // 0 = fixed P, 1 = TVTP
                       const arma::mat& P0,         // N x N fixed transition (mode 0)
                       const arma::mat& X,          // T x p covariates (mode 1)
                       const arma::mat& beta,       // N x . TVTP coefficients (mode 1)
                       const arma::vec& xi_init) {  // length N or empty (-> uniform)
  const int T = y.n_rows, K = y.n_cols, N = sigma.n_slices;
  const double neginf = -datum::inf;

  // Per-regime log-densities (N x T). PD check on the un-ridged matrix (as in R),
  // then a 1e-8 ridge for the inverse and log-determinant (consistently).
  mat logdens(N, T);
  const double cst = K * std::log(2.0 * M_PI);
  for (int m = 0; m < N; ++m) {
    mat Rm = sigma.slice(m);
    vec ev;
    if (!eig_sym(ev, Rm) || ev.min() < 1e-8) return neginf;
    mat Rr = Rm + 1e-8 * eye(K, K);
    mat invR;
    if (!inv(invR, Rr)) return neginf;
    double ld, sgn;
    log_det(ld, sgn, Rr);
    vec quad = sum((y * invR) % y, 1);            // rowSums((y %*% invR) * y)
    logdens.row(m) = (-0.5 * (ld + quad + cst)).t();
  }

  vec xi = (xi_init.n_elem == (uword) N) ? xi_init / accu(xi_init)
                                         : vec(N).fill(1.0 / N);
  double loglik = 0.0;
  const int p = (mode == 1) ? (int) X.n_cols : 0;

  for (int t = 0; t < T; ++t) {
    mat Pt = (mode == 0) ? P0 : tvtp_P(X.row(t), beta, N, p);
    vec pred = Pt.t() * xi;                        // predicted state probabilities
    vec ld = logdens.col(t);
    double c = ld.max();
    vec w = pred % exp(ld - c);
    double sw = accu(w);
    if (!std::isfinite(sw) || sw <= 0.0) return neginf;
    xi = w / sw;
    loglik += std::log(sw) + c;
  }
  return loglik;
}
