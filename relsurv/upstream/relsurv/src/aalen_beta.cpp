// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace arma;  // use the Armadillo library for matrix computations
using namespace Rcpp;

// [[Rcpp::export]]
List Yt(DataFrame data, NumericVector times) {
  // Obtain an at-risk list which is of length length(times).
  // Every element is a vector of length nrow(data).

  int n = times.size();
  int m = data.nrows();

  List mat_list(n);

  for (int i = 0; i < n; i++) {
    double t = times[i];

    IntegerVector mat(m);
    NumericVector tstart = data["start"];
    NumericVector tstop = data["Y"];

    for (int j = 0; j < m; j++) {
      if (t <= tstop[j] && tstart[j] < t)
        mat[j] = 1;
      else
        mat[j] = 0;
    }

    mat_list[i] = mat;
  }

  return mat_list;
}


// [[Rcpp::export]]
List dNt(DataFrame data, NumericVector times) {
  // Obtain the dNt list which is of length length(times).
  // Every element is a vector of length nrow(data).

  int nrow = data.nrows();
  int ncol = times.size();

  // Initialize dNt matrix with zeros
  List dNt(ncol);
  NumericVector dY = data["Y"];
  NumericVector dstat = data["stat"];

  for (int i = 0; i < ncol; i++) {
    NumericVector dNt_tmp(nrow);

    for (int j = 0; j < nrow; j++) {
      // Check if data$Y is equal to times
      if (dY[j] == times[i] && dstat[j] == 1) {
        dNt_tmp[j] = 1;
      } else{
        dNt_tmp[j] = 0;
      }
    }
    dNt[i] = dNt_tmp;
  }

  return dNt;
}


// [[Rcpp::export]]
NumericMatrix prepareX(IntegerVector Yt, NumericMatrix xt) {
  // Prepare a (n x p+1) matrix for (Intercept, X). Used for estimating beta.

  int nrow = Yt.size();
  int ncol = xt.ncol()+1;

  NumericMatrix xt_tmp(nrow, ncol);
  for (int i = 0; i < nrow; i++) {

    xt_tmp(i, 0) = Yt(i);

    for (int j = 1; j < ncol; j++) {
      xt_tmp(i, j) = Yt(i) * xt(i, j-1);
    }
  }

  return xt_tmp;
}


// [[Rcpp::export]]
arma::colvec fitOLS(arma::mat mX, arma::vec dNt, IntegerVector Yt) {
  // Run OLS, doesn't save Xminus in the output.

  int no_cov = mX.n_cols;
  int no_at_risk = sum(Yt);

  arma::vec vBeta = arma::vec(no_cov);

  // compute the OLS estimator
  if(no_at_risk >= no_cov){
    arma::mat mXtX = mX.t()*mX;
    double rcf = rcond(mXtX);

    if(rcf != 0){
    vBeta =  solve(mXtX, mX.t()*dNt);
    // arma::colvec vBeta =  solve(mX, dNt);
    }
  }

  return vBeta;
}

// // [[Rcpp::export]]
// arma::mat matrixProduct(arma::mat mX, int no_cov) {
//
//   arma::mat mat_prod = arma::mat(no_cov, no_cov);
//   arma::vec v1=mX.col(0);
//   arma::vec v2=mX.col(0);
//
//   for (int i = 0; i < no_cov; ++i) {
//     v1 = mX.col(i);
//     for (int j = i; j < no_cov; ++j) {
//       v2 = mX.col(j);
//
//       // mat_prod(i, j) = dot(mX.col(i), mX.col(j));
//      mat_prod(i, j) = dot(v1, v2);
//      mat_prod(j, i) = mat_prod(i, j);
//     }
//   }
//   return mat_prod;
// }
//
// // [[Rcpp::export]]
// arma::mat matrixProduct2(arma::mat mX) {
//   arma::mat mat_prod=mX.t()*mX;
//   return mat_prod;
// }



// [[Rcpp::export]]
List fitOLS2(arma::mat mX, arma::vec dNt, IntegerVector Yt) {
  // Run OLS, save Xminus in the output.

  int no_cov = mX.n_cols;
  int sample_size = mX.n_rows;
  int no_at_risk = sum(Yt);

  arma::vec vBeta = arma::vec(no_cov);
  arma::mat Xminus = arma::mat(no_cov, sample_size);

  List out(2);
  out[0] = vBeta;
  out[1] = Xminus;

  // compute the OLS estimator
  if(no_at_risk >= no_cov){
    arma::mat mXtX = mX.t()*mX;

    double rcf = rcond(mXtX);

    if(rcf != 0){
      Xminus = arma::inv(mXtX)*mX.t();
      vBeta = Xminus*dNt;

      out[0] = vBeta;
      out[1] = Xminus;
    }
  }
  return out;
}

// [[Rcpp::export]]
List fitOLSconst(arma::mat mX, arma::mat mZ, arma::vec dNt, IntegerVector Yt) {
  // Run estimation with constant effects.

  int no_cov = mX.n_cols;
  int no_cov_Z = mZ.n_cols;
  int sample_size = mX.n_rows;
  int no_at_risk = sum(Yt);

  // arma::vec vBeta = arma::vec(no_cov);
  arma::mat Xminus = arma::mat(no_cov, sample_size);
  arma::mat H = arma::mat(no_cov, no_cov);
  arma::mat Identity = arma::mat(sample_size, sample_size, fill::eye);

  arma::mat prvaKomponenta = arma::mat(no_cov_Z, no_cov_Z);
  arma::mat drugaKomponenta = arma::mat(no_cov_Z, 1);

  List out(3);
  out[0] = prvaKomponenta;
  out[1] = drugaKomponenta;
  out[2] = Xminus;

  // compute the OLS estimator
  if(no_at_risk >= no_cov){
    arma::mat mXtX = mX.t()*mX;

    double rcf = rcond(mXtX);

    if(rcf != 0){
      Xminus = arma::inv(mXtX)*mX.t();

      H = Identity-mX*Xminus;
      prvaKomponenta = mZ.t()*H*mZ;
      drugaKomponenta = mZ.t()*H*dNt;

      // vBeta = Xminus*dNt;

      out[0] = prvaKomponenta;
      out[1] = drugaKomponenta;
      out[2] = Xminus;
    }
  }
  return out;
}

// Rcpp implementation of unlist-like functionality
// [[Rcpp::export]]
NumericVector rcpp_unlist(List listObject) {
  int total_length = 0;
  int sajz = listObject.size();

  for (int i = 0; i < sajz; ++i) {
    total_length += as<NumericVector>(listObject[i]).size();
  }

  NumericVector result(total_length);
  int pos = 0;

  for (int i = 0; i < sajz; ++i) {
    NumericVector current = as<NumericVector>(listObject[i]);
    for (int j = 0; j < current.size(); ++j) {
      result[pos] = current[j];
      ++pos;
    }
  }

  return result;
}

// // [[Rcpp::export]]
// arma::cube build_array(arma::vec x, arma::vec dimensions) {
//   // TALE VARIANTA DELA SAMO ZA 3D objekte.
//
//   // Initialize empty cube (Order: slices, columns, rows)
//   arma::cube Cube1(dimensions[0], dimensions[1], dimensions[2]);
//
//   // Fill cube by values of vector x
//   std::copy(x.begin(), x.end(), Cube1.begin());
//
//   return Cube1;
// }

// [[Rcpp::export]]
NumericVector build_array3(NumericVector x, IntegerVector dimensions) {

  // int dimenzije = dimensions.n_elem;

  x.attr("dim") = dimensions;

  return x;

}

// // [[Rcpp::export]]
// NumericMatrix matrix_subset(arma::mat x,
//                         arma::uvec wrow,
//                         arma::uvec wcol
//                         ) {
//   // Take subset of matrix.
//
//   // y must be an integer between 0 and columns - 1
//   // Allows for repeated draws from same columns.
//   x = x.cols( wcol );
//   x = x.rows( wrow );
//
//   NumericMatrix x2 = wrap(x);
//   // NumericMatrix x2 = arma::conv_to<NumericMatrix>::from(wrap(x));
//   return x2;
// }

// declare expc:
extern "C" SEXP expc(SEXP   efac2,   SEXP edims2, SEXP   ecut2,     SEXP   expect2, SEXP   x2, 	SEXP   y2);

// // [[Rcpp::export]]
// List dLambdaP(NumericMatrix data, NumericVector all_times, NumericVector event_times,
//                        NumericVector ratetable, List atts) {
//
//   // Parameters for expc:
//   LogicalVector fk = as<IntegerVector>(atts["factor"]) != 1;
//   int nfk = fk.length();
//   List cuts = atts["cutpoints"];
//
//   IntegerVector atts_type = as<IntegerVector>(atts["type"]);
//   int ltype = atts_type.length();
//   IntegerVector rfac(ltype);
//   for (int i = 0; i < ltype; ++i) {
//     if (atts_type[i] == 1){
//       rfac[i] = 1;
//     }
//   }
//
//   IntegerVector adim = as<IntegerVector>(atts["dim"]);
//   NumericVector acuts = as<NumericVector>(rcpp_unlist(cuts));
//
//   // Yt:
//   List Yt_all = Yt(data, all_times);
//
//   int ltimes = all_times.size();
//   int nr = data.nrow();
//
//   // Prepare ratetable:
//   NumericVector ratetable2 = build_array3(ratetable, adim);
//
//   NumericMatrix outcome(nr, ltimes);
//
//   // Prepare objects:
//   double tstart;
//   double tstop;
//   int sY;
//   int lY;
//   int jj;
//   IntegerVector Yt_all_i;
//   IntegerVector choose_cols = seq(3, nfk + 2);
//   arma::mat data_armamat = as<arma::mat>(data);
//
//   // Go through all times:
//   for (int i = 0; i < ltimes; ++i) {
//     if (i == 0) {
//       if (all_times[i] == 0) continue;
//       else tstart = 0;
//     } else {
//       tstart = all_times[i - 1];
//     }
//
//     tstop = all_times[i];
//
//     // Yt at i-th time:
//     Yt_all_i = Yt_all[i];
//     // At-risk size:
//     sY = sum(Yt_all_i);
//     // Sample size:
//     lY = Yt_all_i.length();
//     // Find those that are at-risk (which):
//     IntegerVector at_risk(sY);
//
//     jj = 0;
//     for (int k = 0; k < lY; ++k) {
//       if (Yt_all_i[k]==1) {
//         at_risk[jj] = k;
//         jj = jj + 1;
//       }
//     }
//
//     // Convert objects:
//     arma::uvec roows = arma::conv_to<arma::uvec>::from(as<arma::vec>(wrap(at_risk)));
//     arma::uvec cools = arma::conv_to<arma::uvec>::from(as<arma::vec>(wrap(choose_cols)));
//
//     // Find data subset:
//     NumericMatrix data_tmp = matrix_subset(data_armamat, roows, cools);
//     // NumericMatrix data_tmp(sY, nfk);
//
//     NumericVector tstart_vec = rep(tstart, sY);
//     // Increase age and year by tstart:
//     for (int j = 0; j < nfk; ++j) {
//       if(fk[j]==true){
//         for(int ji = 0; ji < nfk; ++ji){
//           data_tmp(ji, j) += tstart_vec[ji];
//         }
//       }
//     }
//
//     // Prepare times vector:
//     NumericVector times(sY, tstop - tstart);
//
//     // Run expc:
//     List pop_survs0 = as<List>(expc(wrap(rfac), wrap(adim), wrap(acuts), wrap(ratetable2), wrap(data_tmp), wrap(times)));
//     NumericVector pop_survs = pop_survs0["surv"];
//     // NumericVector pop_survs = as<NumericVector>(at_risk);
//
//     // Save hazards:
//     for(int ji = 0; ji < sY; ++ji){
//       int at_risk_ji = at_risk[ji];
//       outcome(at_risk_ji, i) = -log(pop_survs[ji]);
//     }
//   }
//
//   // Cumulative hazards:
//   for (int i = 0; i < nr; ++i) {
//     for (int j = 1; j < ltimes; ++j) {
//       outcome(i, j) += outcome(i, j - 1);
//     }
//   }
//
//   // // Hazards at event times only:
//   // List outcome_l(event_times.size());
//   // int j = 0;
//   // for (int i = 0; i < ltimes; ++i) {
//   //   if (std::find(event_times.begin(), event_times.end(), all_times[i]) != event_times.end()) {
//   //     outcome_l[j] = outcome(_, i);
//   //     ++j;
//   //   }
//   // }
//
//   // Hazards at event times only:
//   List outcome_l(all_times.size());
//   for (int i = 0; i < ltimes; ++i) {
//     outcome_l[i] = outcome(_, i);
//   }
//
//   return outcome_l;
// }

// // [[Rcpp::export]]
// List calculateBetas(DataFrame data, NumericMatrix xt, NumericVector event_times, int var_estimator) {
// // List calculateBetas(DataFrame data, NumericMatrix xt, NumericVector event_times) {
//   // Run OLS at all event times.
//
//   int ncol = event_times.size();
//
//   List Yt_val = Yt(data, event_times);
//   List dNt_val = dNt(data, event_times);
//
//   List betas_list(ncol);
//   // NumericMatrix betas_list(ncol, nrow+1);
//
//   int sample_size = xt.nrow();
//   int number_covs = xt.ncol();
//   arma::mat diag_dNt = arma::mat(sample_size, sample_size);
//   arma::mat beta_var = arma::mat(number_covs, number_covs);
//   List betas_var_list(ncol);
//
//   for (int i = 0; i < ncol; i++) {
//
//     NumericMatrix xx1 = prepareX(Yt_val[i], xt);
//
//     arma::mat xx2 = as<arma::mat>(xx1);
//     arma::vec dNti = as<arma::vec>(dNt_val[i]);
//
//     // arma::vec betas = fitOLS(xx2, dNti, Yt_val[i]);
//     List betas = fitOLS2(xx2, dNti, Yt_val[i]);
//
//     // betas_list[i] = betas;
//     betas_list[i] = betas[0];
//
//     if(var_estimator == 1){
//       diag_dNt.diag() = dNti;
//     }
//     if(var_estimator == 2){
//       arma::vec betas0_vec = betas[0];
//       arma::mat betas0(betas0_vec);
//       betas0.reshape(betas0_vec.size(), 1);
//
//       diag_dNt.diag() = xx2*betas0;
//     }
//
//     arma::mat betas1 = betas[1];
//     betas_var_list[i] = betas1*diag_dNt*betas1.t();
//   }
//
//   List out(2);
//   out[0] = betas_list;
//   out[1] = betas_var_list;
//
//   // return betas_list;
//   return out;
// }

// // [[Rcpp::export]]
// List calculateBetasRelsurv(DataFrame data, NumericMatrix xt, NumericVector event_times, NumericVector all_times,
//                            NumericVector ratetable, List atts, NumericMatrix data_mat) {
//   // Run OLS at all event times.
//
//   int ncol = event_times.size();
//
//   List Yt_val = Yt(data, event_times);
//   List dNt_val = dNt(data, event_times);
//   List dLambdaP_val = dLambdaP(data_mat, all_times, event_times, ratetable, atts);
//
//   List betas_list(ncol);
//
//   arma::vec dL(data.nrows());
//
//   for (int i = 0; i < ncol; i++) {
//
//     NumericMatrix xx1 = prepareX(Yt_val[i], xt);
//
//     arma::mat xx2 = as<arma::mat>(xx1);
//     arma::vec dNti = as<arma::vec>(dNt_val[i]);
//     arma::vec dLambdaPi = as<arma::vec>(dLambdaP_val[i]);
//
//     if(i>0){
//       dL = as<arma::vec>(dLambdaP_val[i-1]);
//     }
//
//     arma::vec betas = fitOLS(xx2, dNti-dLambdaPi+dL, Yt_val[i]);
//
//     betas_list[i] = betas;
//   }
//
//   return betas_list;
// }

// // [[Rcpp::export]]
// List calculateBetasRelsurv2(DataFrame data, NumericMatrix xt, NumericVector event_times, NumericVector all_times,
//                            NumericVector ratetable, List atts, NumericMatrix data_mat) {
//   // Run OLS at all times.
//
//   int ncol = all_times.size();
//
//   List Yt_val = Yt(data, all_times);
//   List dNt_val = dNt(data, all_times);
//   List dLambdaP_val = dLambdaP(data_mat, all_times, event_times, ratetable, atts);
//
//   List betas_list(ncol);
//
//   arma::vec dL(data.nrows());
//
//   for (int i = 0; i < ncol; i++) {
//
//     NumericMatrix xx1 = prepareX(Yt_val[i], xt);
//
//     arma::mat xx2 = as<arma::mat>(xx1);
//     arma::vec dNti = as<arma::vec>(dNt_val[i]);
//     arma::vec dLambdaPi = as<arma::vec>(dLambdaP_val[i]);
//
//     if(i>0){
//       dL = as<arma::vec>(dLambdaP_val[i-1]);
//     }
//
//     arma::vec betas = fitOLS(xx2, dNti-dLambdaPi+dL, Yt_val[i]);
//     // arma::vec betas_null = fitOLS(xx2, dNti, Yt_val[i]);
//     //
//     // double dPi = sum(dLambdaPi-dL);
//     // IntegerVector Yti_vec = Yt_val[i];
//     // double Yti = sum(Yti_vec);
//     //
//     // if(Yti>0){
//     //   if(betas[0] != 0){
//     //     betas[0] = betas_null[0] - dPi*365.241/Yti;
//     //   }
//     // }
//
//     betas_list[i] = betas;
//   }
//
//   return betas_list;
// }


// // [[Rcpp::export]]
// List calculateBetasRelsurv22(DataFrame data, NumericMatrix xt, NumericVector event_times, NumericVector all_times,
//                             NumericVector ratetable, List atts, NumericMatrix data_mat, int var_estimator) {
//   // Run OLS at all times.
//
//   int ncol = all_times.size();
//
//   List Yt_val = Yt(data, all_times);
//   List dNt_val = dNt(data, all_times);
//   List dLambdaP_val = dLambdaP(data_mat, all_times, event_times, ratetable, atts);
//
//   List betas_list(ncol);
//   arma::vec dL(data.nrows());
//
//   int sample_size = xt.nrow();
//   int number_covs = xt.ncol();
//   arma::mat diag_dNt = arma::mat(sample_size, sample_size);
//   arma::mat beta_var = arma::mat(number_covs, number_covs);
//   List betas_var_list(ncol);
//
//   for (int i = 0; i < ncol; i++) {
//
//     NumericMatrix xx1 = prepareX(Yt_val[i], xt);
//
//     arma::mat xx2 = as<arma::mat>(xx1);
//     arma::vec dNti = as<arma::vec>(dNt_val[i]);
//     arma::vec dLambdaPi = as<arma::vec>(dLambdaP_val[i]);
//
//     if(i>0){
//       dL = as<arma::vec>(dLambdaP_val[i-1]);
//     }
//
//     // arma::vec betas = fitOLS(xx2, dNti-dLambdaPi+dL, Yt_val[i]);
//     List betas = fitOLS2(xx2, dNti-dLambdaPi+dL, Yt_val[i]);
//
//     // betas_list[i] = betas;
//     betas_list[i] = betas[0];
//
//     if(var_estimator == 1){
//       diag_dNt.diag() = dNti-dLambdaPi+dL;
//     }
//     if(var_estimator == 2){
//       arma::vec betas0_vec = betas[0];
//       arma::mat betas0(betas0_vec);
//       betas0.reshape(betas0_vec.size(), 1);
//       diag_dNt.diag() = xx2*betas0;
//     }
//     arma::mat betas1 = betas[1];
//     betas_var_list[i] = betas1*diag_dNt*betas1.t();
//
//   }
//
//   List out(2);
//   out[0] = betas_list;
//   out[1] = betas_var_list;
//   // return betas_list;
//   return out;
// }
