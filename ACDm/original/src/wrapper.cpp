// File: src/wrapper.cpp

#include <Rcpp.h>

// forward‐declare your old C function
extern "C" void computeDurationsSubSec(int *y,
                                       int *M,
                                       int *d,
                                       int *h,
                                       int *m,
                                       double *s,
                                       int *yDur,
                                       int *MDur,
                                       int *dDur,
                                       int *hDur,
                                       int *mDur,
                                       double *sDur,
                                       int *vol,
                                       double *price,
                                       int *volDur,
                                       double *priceDur,
                                       int *Ntrans,
                                       double *dur,
                                       int *Ntime,
                                       int *Ndur,
                                       double *open,
                                       double *close,
                                       int *durType,
                                       int *zeroDurHandeling,
                                       double *priceChange,
                                       int *culmVol);

// [[Rcpp::export]]
Rcpp::List computeDurationsSubSec_wrap(
    Rcpp::IntegerVector y,
    Rcpp::IntegerVector M,
    Rcpp::IntegerVector d,
    Rcpp::IntegerVector h,
    Rcpp::IntegerVector m,
    Rcpp::NumericVector s,
    Rcpp::IntegerVector vol,
    Rcpp::NumericVector price,
    double open,
    double close,
    int durType_scalar,
    int zeroDurHandeling_scalar,
    double priceChange_scalar,
    int cumVol_scalar
) {
  int Ntime = y.size();
  
  // allocate maximal buffers
  Rcpp::IntegerVector yDur(Ntime), MDur(Ntime), dDur(Ntime),
  hDur(Ntime), mDur(Ntime), volDur(Ntime),
  Ntrans_vec(Ntime);
  Rcpp::NumericVector sDur(Ntime), priceDur(Ntime), dur(Ntime);
  Rcpp::IntegerVector Ndur(1);
  
  // call the C function
  computeDurationsSubSec(
    y.begin(), M.begin(), d.begin(), h.begin(), m.begin(), s.begin(),
    yDur.begin(), MDur.begin(), dDur.begin(), hDur.begin(), mDur.begin(), sDur.begin(),
    vol.begin(), price.begin(), volDur.begin(), priceDur.begin(),
    Ntrans_vec.begin(), dur.begin(), &Ntime, Ndur.begin(),
    &open, &close,
    &durType_scalar, &zeroDurHandeling_scalar,
    &priceChange_scalar, &cumVol_scalar
  );
  
  int nout = Ndur[0];
  Rcpp::Range idx(0, nout - 1);
  
  // build the list item-by-item
  Rcpp::List out;
  out["y"                   ] = y;
  out["M"                   ] = M;
  out["d"                   ] = d;
  out["h"                   ] = h;
  out["m"                   ] = m;
  out["s"                   ] = s;
  
  out["yDur"                ] = yDur[idx];
  out["MDur"                ] = MDur[idx];
  out["dDur"                ] = dDur[idx];
  out["hDur"                ] = hDur[idx];
  out["mDur"                ] = mDur[idx];
  out["sDur"                ] = sDur[idx];
  
  out["vol"                 ] = vol;
  out["price"               ] = price;
  out["volDur"              ] = volDur[idx];
  out["priceDur"            ] = priceDur[idx];
  
  out["Ntrans"              ] = Ntrans_vec[idx];
  out["dur"                 ] = dur[idx];
  
  out["Ntime"               ] = Ntime;
  out["Ndur"                ] = nout;
  
  out["open"                ] = open;
  out["close"               ] = close;
  out["durType"             ] = durType_scalar;
  out["zeroDurHandeling"    ] = zeroDurHandeling_scalar;
  out["priceChange"         ] = priceChange_scalar;
  out["cumVol"              ] = cumVol_scalar;
  
  return out;
}


#include <Rcpp.h>
// forward‐declare the old C function
extern "C" void computeDurationsShort(int *y,
                                     int *M,
                                     int *d,
                                     int *h,
                                     int *m,
                                     double *s,
                                     int *yDur,
                                     int *MDur,
                                     int *dDur,
                                     int *hDur,
                                     int *mDur,
                                     double *sDur,
                                     double *dur,
                                     int *Ndur,
                                     int *Ntrans,
                                     int *Ntime,
                                     int *open,
                                     int *close,
                                     int *zeroDurHandeling);

// [[Rcpp::export]]
Rcpp::List computeDurationsShort_wrap(
    Rcpp::IntegerVector y,
    Rcpp::IntegerVector M,
    Rcpp::IntegerVector d,
    Rcpp::IntegerVector h,
    Rcpp::IntegerVector m,
    Rcpp::NumericVector s,
    int open,
    int close,
    int zeroDurHandeling
) {
  int Ntime = y.size();
  
  // allocate max‐length buffers for per‐duration outputs
  Rcpp::IntegerVector yDur(Ntime),
  MDur(Ntime),
  dDur(Ntime),
  hDur(Ntime),
  mDur(Ntime),
  Ntrans_vec(Ntime);
  Rcpp::NumericVector sDur(Ntime),
  dur(Ntime);
  Rcpp::IntegerVector Ndur(1);
  
  // call the C routine
  computeDurationsShort(
    y.begin(), M.begin(), d.begin(), h.begin(), m.begin(), s.begin(),
    yDur.begin(), MDur.begin(), dDur.begin(), hDur.begin(), mDur.begin(), sDur.begin(),
    dur.begin(), Ndur.begin(),   // outputs #13–14
    Ntrans_vec.begin(),         // output #15
    &Ntime,                     // output #16
    &open, &close, &zeroDurHandeling
  );
  
  // how many durations were actually written?
  int nout = Ndur[0];
  Rcpp::Range idx(0, nout - 1);
  
  // build the list of 19 elements in exact C‐arg order:
  Rcpp::List out;
  // 1–6: inputs
  out["y"]               = y;
  out["M"]               = M;
  out["d"]               = d;
  out["h"]               = h;
  out["m"]               = m;
  out["s"]               = s;
  
  // 7–12: per-duration date/time
  out["yDur"]            = yDur[idx];
  out["MDur"]            = MDur[idx];
  out["dDur"]            = dDur[idx];
  out["hDur"]            = hDur[idx];
  out["mDur"]            = mDur[idx];
  out["sDur"]            = sDur[idx];
  
  // 13: per-duration durations
  out["dur"]             = dur[idx];
  // 14: scalar Ndur
  out["Ndur"]            = nout;
  // 15: per-duration Ntrans
  out["Ntrans"]          = Ntrans_vec[idx];
  // 16: scalar Ntime
  out["Ntime"]           = Ntime;
  
  // 17–19: the three scalar parameters
  out["open"]            = open;
  out["close"]           = close;
  out["zeroDurHandeling"] = zeroDurHandeling;
  
  return out;
}

