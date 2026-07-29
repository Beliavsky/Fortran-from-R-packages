#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

// Functions written by Dominik Schulz, 12/05/2020

// Function to help obtain a sequence from 'from' to 'to' by step 1
// in C++
arma::vec seqCpp(const int from, const int to){
  const int n = to - from + 1;
  arma::vec seqOut(n);
  for (int i = from; i < to + 1; ++i) {
    seqOut(i - from) = i;
  }
  return seqOut;
}

// Function to help obtain a sequence from 'from' to 'to' by step 1
// in C++ (returns a rowvec instead)
arma::rowvec rseqCpp(const int from, const int to) {
  const int n = to - from + 1;
  arma::rowvec seqOut(n);
  for (int i = from; i < to + 1; ++i) {
    seqOut(i - from) = i;
  }
  return seqOut;
}

// Function to calculate the factorial of an 'int' in C++
int factorialCpp(const int k) {
  int fac = 1;
  if (k > 1) {
    for (int i = 2; i < k + 1; ++i) {
      fac *= i;
    }
  }
  return fac;
}

// C++ version of gsmoothCalc
// [[Rcpp::export]]
arma::vec gsmoothCalcCpp(const arma::vec& y, const int v, const int p,
  const int mu, const double b, const int bb) {

  const int n = y.size();
  arma::vec gr(n);
  const int hh = std::trunc(n * b + 0.5);
  const int htm = 2 * hh + 1;
  arma::mat ws(htm, htm);
  arma::rowvec wk = arma::zeros<arma::rowvec>(htm);
  arma::mat xt = arma::zeros<arma::mat>((p + 1), htm);
  arma::mat xw = xt;

  const arma::vec seqhh = seqCpp(0, hh);
  const arma::vec allLowB = hh - seqhh;
  const arma::vec hr = hh + bb * allLowB;
  const arma::vec hr2 = arma::pow(hr + 1, 2.0);
  const arma::vec ht = seqhh + hr;
  arma::mat xa = xt;
  xt.submat(0, 0, 0, ht(0) - 1).fill(1.0);

  const arma::rowvec allsequ = rseqCpp(-hh, hr(0));
  const int m = allsequ.size();
  const arma::rowvec allsequ2 = arma::pow(allsequ, 2.0);
  const arma::rowvec allsequhh = allsequ / hh;
  arma::mat allseqmat(p, m);
  allseqmat.row(0) = allsequhh;
  for (int i = 2; i < p + 1; ++i) {
    allseqmat.row(i - 1) = arma::pow(allsequhh, i);
  }

  const arma::vec allUpB = hh + hr;

  int lowB;
  int upB;

  for (int i = 0; i < hh + 1; ++i) {
    lowB = allLowB(i);
    upB = allUpB(i);
    if (mu == 0) {
      wk.subvec(0, ht(i)).fill(1.0);
    } else if (mu == 1) {
      wk.subvec(0, ht(i)) = 1.0 - allsequ2.subvec(lowB, upB) / hr2(i);
    } else {
      wk.subvec(0, ht(i)) = arma::pow(1.0 - allsequ2.subvec(lowB, upB) / hr2(i), mu);
    }

    xt(0, ht(i)) = 1.0;
    xt.submat(1, 0, p, ht(i)) = allseqmat.cols(lowB, upB);
    xw = xt * diagmat(wk);
    xa = arma::solve(xw * xt.t(), xw);
    ws.row(i) = xa.row(v);
  }
  ws.rows(hh + 1, htm - 1) = std::pow(-1.0, v) * arma::flipud(arma::fliplr(ws.submat(0, 0, hh - 1, htm - 1)));
  ws = factorialCpp(v) * ws * std::pow(n / double(hh), v);
  arma::mat ym(htm, n - htm + 1);
  for (int i = hh; i < n - hh; ++i) {
    ym.col(i - hh) = y.subvec(i - hh, i + hh);
  }
  gr.subvec(0, hh - 1) = ws.rows(0, hh - 1) * y.subvec(0, htm - 1);
  gr.subvec(hh, n - hh - 1) = (ws.row(hh) * ym).t();
  gr.subvec(n - hh, n - 1) = ws.rows(hh + 1, htm - 1) * y.subvec(n - htm, n - 1);
  return gr;
}

// C++ version of gsmoothCalc2
// [[Rcpp::export]]
List gsmoothCalc2Cpp(const arma::vec& y, const int v, const int p,
  const int mu, const double b, const int bb) {

  const int n = y.size();
  arma::vec gr(n);
  const int hh = std::trunc(n * b + 0.5);
  const int htm = 2 * hh + 1;
  arma::mat ws(htm, htm);
  arma::rowvec wk = arma::zeros<arma::rowvec>(htm);
  arma::mat xt = arma::zeros<arma::mat>((p + 1), htm);
  arma::mat xw = xt;

  const arma::vec seqhh = seqCpp(0, hh);
  const arma::vec allLowB = hh - seqhh;
  const arma::vec hr = hh + bb * allLowB;
  const arma::vec hr2 = arma::pow(hr + 1, 2.0);
  const arma::vec ht = seqhh + hr;
  arma::mat xa = xt;
  xt.submat(0, 0, 0, ht(0) - 1).fill(1.0);

  const arma::rowvec allsequ = rseqCpp(-hh, hr(0));
  const int m = allsequ.size();
  const arma::rowvec allsequ2 = arma::pow(allsequ, 2.0);
  const arma::rowvec allsequhh = allsequ / hh;
  arma::mat allseqmat(p, m);
  allseqmat.row(0) = allsequhh;
  for (int i = 2; i < p + 1; ++i) {
    allseqmat.row(i - 1) = arma::pow(allsequhh, i);
  }

  const arma::vec allUpB = hh + hr;

  int lowB;
  int upB;

  for (int i = 0; i < hh + 1; ++i) {

    lowB = allLowB(i);
    upB = allUpB(i);

    if (mu == 0) {
      wk.subvec(0, ht(i)).fill(1.0);
    } else if (mu == 1) {
      wk.subvec(0, ht(i)) = 1.0 - allsequ2.subvec(lowB, upB) / hr2(i);
    } else {
      wk.subvec(0, ht(i)) = arma::pow(1.0 - allsequ2.subvec(lowB, upB) / hr2(i), mu);
    }

    xt(0, ht(i)) = 1.0;
    xt.submat(1, 0, p, ht(i)) = allseqmat.cols(lowB, upB);

    xw = xt * diagmat(wk);
    xa = arma::solve(xw * xt.t(), xw);
    ws.row(i) = xa.row(v);
  }

  ws.rows(hh + 1, htm - 1) = std::pow(-1.0, v) * arma::flipud(arma::fliplr(ws.submat(0, 0, hh - 1, htm - 1)));
  ws = factorialCpp(v) * ws * std::pow(n / double(hh), v);
  arma::mat ym(htm, n - htm + 1);
  for (int i = hh; i < n - hh; ++i) {
    ym.col(i - hh) = y.subvec(i - hh, i + hh);
  }
  gr.subvec(0, hh - 1) = ws.rows(0, hh - 1) * y.subvec(0, htm - 1);
  gr.subvec(hh, n - hh - 1) = (ws.row(hh) * ym).t();
  gr.subvec(n - hh, n - 1) = ws.rows(hh + 1, htm - 1) * y.subvec(n - htm, n - 1);
  List listOut = List::create(_["ye"] = gr,
                              _["ws"] = ws);

  return listOut;

}

// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically
// run after the compilation.
//

/*** R

*/
