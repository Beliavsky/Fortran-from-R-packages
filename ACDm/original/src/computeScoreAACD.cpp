#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export(name = ".sumOuter")]]
List sumOuter(NumericMatrix d_LL_d_theta,
              NumericMatrix d_psi_d_theta,
              NumericVector psi){
  
  int N = psi.size(), k = d_LL_d_theta.ncol();
  NumericMatrix A(k, k), B(k, k);
  
  for(int i = 0; i < N ; i++){
    for(int col = 0; col < k; col++){  //column
      for(int row = col; row < k ; row++){ //row
        
        B[row + col * k] += d_LL_d_theta[i + row * N] * d_LL_d_theta[i + col * N];
        A[row + col * k] -= d_psi_d_theta[i + row * N] * d_psi_d_theta[i + col * N] / psi[i] / psi[i];

      }
    }
  }
  
  // fills the upper right half:
  for(int col = 1; col < k; col++){  //column
    for(int row = 0; row < col ; row++){ //row
      B[row + col * k] = B[col + row * k];
      A[row + col * k] = A[col + row * k];
    }
  }
  
  List returnList;
  returnList["B"] = B;
  returnList["A"] = A;
  
  return(returnList);
}

double abs_prim(double x){
  
  if(x > 0) return(1.0);
  if(x < 0) return(-1.0);
  else return(0.0);
}

NumericVector get_d_LL_d_psi_distCpp(NumericVector x,
                                     NumericVector psi,
                                     NumericVector logPsi,
                                     NumericVector e,
                                     NumericVector d_psi_d_theta,
                                     NumericVector d_LL_d_theta,
                                     int dist,
                                     NumericVector distPara,
                                     NumericVector LLi,
                                     int NpsiPara,
                                     NumericVector score){
  
  NumericVector LL(1);
  int N = x.size();
  
  switch (dist) {
  case 1: // Exponential distribution
  {
    double multi3;
    for (int i = 0; i < N; i++) {
      
      // calculates LLi and LL:
      LLi[i] = -logPsi[i] - e[i];
      LL[0] += LLi[i];
      
      // calculates d_LL_d_theta:
      multi3 = (e[i] - 1)/ psi[i];
      for(int col = 0; col < NpsiPara; col++) d_LL_d_theta(i, col) = multi3 * d_psi_d_theta(i, col);
    }
    break;
  }
  case 2: //Weibull distribution
  {
    double gammaFncValue = tgamma(1 + 1 / distPara[0]);
    double phi_pow, factor1;
    for (int i = 0; i < N; i++) {
      
      phi_pow = pow(gammaFncValue * e[i], distPara[0]);
      
      // calculates LLi and LL:
      LLi[i] = log(distPara[0] / x[i]) + 
        log(phi_pow) - phi_pow;
      LL[0] += LLi[i];
      
      // calculates d_LL_d_theta:
      factor1 = - distPara[0] / psi[i] * (1 + phi_pow);
      for(int col = 0; col < NpsiPara; col++) d_LL_d_theta(i, col) = factor1 * d_psi_d_theta(i, col);
      d_LL_d_theta(i, NpsiPara) = (1 + (log(phi_pow) - Rf_digamma(1 / distPara[0] + 1)) 
                                     * (1 - phi_pow)) / distPara[0];
      
    }
  }
    break;
  default:
    Rcpp::stop("distribution not yet implemented in 'get_d_LL_d_psi_distCpp'");
  break;
  }
  
  //adds to the score:
  for (int i = 0; i < N; i++) {
    for(int col = 0; col < NpsiPara; col++){
      score[col] += d_LL_d_theta(i, col);
    }
  }
  
  return(LL);
}


NumericVector getLL_distCpp(NumericVector x,
                            NumericVector psi,
                            NumericVector logPsi,
                            NumericVector e,
                            int dist,
                            NumericVector distPara){
  
  NumericVector LL(1);
  int N = x.size();
  
  switch (dist) {
  case 1: //Exponential distribution
  {
    for (int i = 0; i < N; i++) {
      LL[0] += -logPsi[i] - e[i];
    }
    
  }
    break;
  case 2: //Weibull distribution
  {
    double gammaFncValue = tgamma(1 + 1 / distPara[0]);
    for (int i = 0; i < N; i++) {
      LL[0] += log(distPara[0] / x[i]) + 
        distPara[0] * log(gammaFncValue * e[i]) - pow(gammaFncValue * e[i], distPara[0]);
    }
  }
    break;
  default:
    Rcpp::stop("distribution not yet implemented in 'getLL_distCpp'");
  break;
  }
  
  return(LL);
}

// [[Rcpp::export(name = ".computeScoreTACD")]]
List computeScoreTACD(NumericVector x,
                      NumericVector param,
                      NumericVector order,
                      double mean,
                      int dist,
                      NumericVector distPara,
                      IntegerVector newDayR,
                      int forceErrExpec,
                      int returnIndex,
                      int startType,
                      NumericVector bp){
  
  int N = x.size();
  
  NumericVector psi(N); // the conditional durations 
  NumericVector logPsi(N); // log(psi)
  NumericVector e(N); // the errors (x/psi)
  NumericVector LL(1); // log-likelihood
  
  NumericVector regime(N); //the current regime
  
  List returnList; 
  
  // extracts the parameters:
  int p = order(0);
  int q = order(1);
  int M = bp.size();
  int maxPQ = std::max(p, q);
  NumericVector omega = param[Rcpp::seq(0, M)];
  NumericVector alpha; if(p > 0) alpha = param[Rcpp::seq(M + 1, (p+1) * (M+1) - 1)];
  NumericVector beta; if(q > 0) beta = param[Rcpp::seq((p+1) * (M+1), (q+p+1) * (M+1) - 1)];
  
  IntegerVector newDayC(newDayR.size() + 2);   
  if(newDayR.size() > 0) newDayC[Rcpp::Range(1, newDayR.size())] = newDayR - 1; 
  newDayC[newDayC.size() - 1] = N;
  
  int dayIndex = 0;
  
  for(int i = 0; i < N; i++){
    
    if(newDayC[dayIndex + 1] == i) dayIndex++;
    
    // gets the current regime:
    regime[i] = 0;
    double tempThreshVar; 
    tempThreshVar = (i == newDayC[dayIndex]) ? mean : x[i - 1]; // if first dur on new day: previous dur assumed to be mean
    if(tempThreshVar >= bp[0]){
      for(int bpIndex = 1; bpIndex < M; bpIndex++){
        if(tempThreshVar < bp[bpIndex]){
          regime[i] = bpIndex;
          break;
        }
      }
      if(tempThreshVar >= bp[M - 1]) regime[i] = M;
    }
    	
    psi[i] = omega[regime[i]]; //adds the constant
    
    //adds the "p"-part:
    for(int j = 0; j < std::min(p, i - newDayC[dayIndex]); j++) psi[i] += alpha[j * M + regime[i]] * x[i - j - 1];
    for(int j = std::min(p, i - newDayC[dayIndex]); j < p; j++) psi[i] += alpha[j * M + regime[i]] * mean;
    
    //adds the "q"-part:
    for(int j = 0; j < std::min(q, i - newDayC[dayIndex]); j++) psi[i] += beta[j * M + regime[i]] * psi[i - j - 1]; 
    for(int j = std::min(q, i - newDayC[dayIndex]); j < q; j++) psi[i] += beta[j * M + regime[i]] * mean;
    
    // if 'startType' == 2, the first max(p,q) psi[i] are set to the mean:
    if(startType == 2 && i - newDayC[dayIndex] >= 0 && i - newDayC[dayIndex] < maxPQ) psi[i] = mean;
    
    e[i] = x[i]/psi[i];
    logPsi[i] = log(psi[i]);
  }
  
  //computes the derivatives:
  if(returnIndex == 3 || returnIndex == 4){
    
    int Npara = (dist == 1) ? param.size() : param.size() + distPara.size();
    NumericMatrix d_psi_d_theta(N, Npara);
    NumericMatrix d_LL_d_theta(N, Npara);
    NumericVector score(Npara);
    NumericVector LLi(N);
    
    int dayIndex = 0;
    
    for(int i = 0; i < N; i++){
      
      if(newDayC[dayIndex + 1] == i) dayIndex++;
      
      // d_psi_d_omega:
      d_psi_d_theta(i, regime[i]) = 1;
      
      // d_psi_d_alpha:
      for(int j = 0; j < std::min(p, i - newDayC[dayIndex]); j++)
        d_psi_d_theta(i, (j+1) * (M+1) + regime[i]) = x[i - j - 1];
      for(int j = std::min(p, i - newDayC[dayIndex]); j < p; j++)
        d_psi_d_theta(i, (j+1) * (M+1) + regime[i]) = mean;
      
      // d_psi_d_beta:
      for(int j = 0; j < std::min(q, i - newDayC[dayIndex]); j++)
        d_psi_d_theta(i, (p+j+1) * (M+1) + regime[i]) = psi[i - j - 1];
      for(int j = std::min(q, i - newDayC[dayIndex]); j < q; j++)
        d_psi_d_theta(i, (p+j+1) * (M+1) + regime[i]) = mean;
      
      // adds "common" part:
      for(int col = 0; col < param.size(); col++){
        for(int j = 0; j < std::min(q, i - newDayC[dayIndex]); j++)
          d_psi_d_theta(i, col) += beta[j * (M+1) + regime[i]] * d_psi_d_theta(i - j - 1, col);
      }
      
      // if 'startType' == 2, the first max(p,q) psi[i] are set to mean (so that d_psi_d_theta = 0)
      if(startType == 2 && i - newDayC[dayIndex] >= 0 && i - newDayC[dayIndex] < maxPQ){
        for(int col = 0; col < Npara; col++){
          d_psi_d_theta(i, col) = 0;
        }
      }
      
    }
    
    // computes d_LL_d_theta:
    LL = get_d_LL_d_psi_distCpp(x, psi, logPsi, e, d_psi_d_theta,
                                d_LL_d_theta, dist, distPara,
                                LLi, param.size(), score);
    
    if(returnIndex == 3 || returnIndex == 4) returnList["score"] = score; 
    if(returnIndex == 4){
      
      returnList["psi"] = psi;
      returnList["e"] = e;
      returnList["d_LL_d_theta"] = d_LL_d_theta; 
      returnList["d_psi_d_theta"] = d_psi_d_theta;
      returnList["LLi"] = LLi;
    }
    
  }
  
  // computes the log-likelihood:
  if(returnIndex != 3 && returnIndex != 4) 
    LL = getLL_distCpp(x, psi, logPsi, e, dist, distPara);
  
  returnList["LL"] = LL;
  if(returnIndex == 2){
    returnList["psi"] = psi;
    returnList["e"] = e;
  }
  
  return(returnList);
}

// [[Rcpp::export(name = ".computeScoreTAMACD")]]
List computeScoreTAMACD(NumericVector x,
                      NumericVector param,
                      NumericVector order,
                      double mean,
                      int dist,
                      NumericVector distPara,
                      IntegerVector newDayR,
                      int forceErrExpec,
                      int returnIndex,
                      int startType,
                      NumericVector bp){
  
  int N = x.size();
  
  NumericVector psi(N); // the conditional durations 
  NumericVector logPsi(N); // log(psi)
  NumericVector e(N); // the errors (x/psi)
  NumericVector LL(1); // log-likelihood
  
  NumericVector regime(N); //the current regime
  
  List returnList; 
  
  // extracts the parameters:
  int p = order(0);
  int r = order(1);
  int q = order(2);
  int M = bp.size();
  int max_prq = std::max(p, std::max(r, q));
  NumericVector omega = param[Rcpp::seq(0, M)];
  NumericVector alpha; if(p > 0) alpha = param[Rcpp::seq(M + 1, (p + 1) * (M + 1) - 1)];
  NumericVector nu; if(r > 0) nu = param[Rcpp::seq((p + 1) * (M + 1), (r + p + 1) * (M + 1) - 1)];
  NumericVector beta; if(q > 0) beta = param[Rcpp::seq((p + r + 1) * (M + 1), (q + p + r + 1) * (M + 1) - 1)];
  
  IntegerVector newDayC(newDayR.size() + 2);   
  if(newDayR.size() > 0) newDayC[Rcpp::Range(1, newDayR.size())] = newDayR - 1; 
  newDayC[newDayC.size() - 1] = N;
  
  int dayIndex = 0;
  
  for(int i = 0; i < N; i++){
    
    if(newDayC[dayIndex + 1] == i) dayIndex++;
    
    // gets the current regime:
    regime[i] = 0;
    double tempThreshVar; 
    tempThreshVar = (i == newDayC[dayIndex]) ? mean : x[i - 1]; // if first dur on new day: previous dur assumed to be mean
    if(tempThreshVar >= bp[0]){
      for(int bpIndex = 1; bpIndex < M; bpIndex++){
        if(tempThreshVar < bp[bpIndex]){
          regime[i] = bpIndex;
          break;
        }
      }
      if(tempThreshVar >= bp[M - 1]) regime[i] = M;
    }
    
    psi[i] = omega[regime[i]]; //adds the constant
    
    //adds the "p"-part:
    for(int j = 0; j < std::min(p, i - newDayC[dayIndex]); j++) psi[i] += alpha[j * M + regime[i]] * x[i - j - 1];
    for(int j = std::min(p, i - newDayC[dayIndex]); j < p; j++) psi[i] += alpha[j * M + regime[i]] * mean;
    
    //adds the "r"-part:
    for(int j = 0; j < std::min(r, i - newDayC[dayIndex]); j++) psi[i] += nu[j * M + regime[i]] * e[i - j - 1];
    for(int j = std::min(r, i - newDayC[dayIndex]); j < r; j++) psi[i] += nu[j * M + regime[i]] * 1;
    
    //adds the "q"-part:
    for(int j = 0; j < std::min(q, i - newDayC[dayIndex]); j++) psi[i] += beta[j * M + regime[i]] * psi[i - j - 1]; 
    for(int j = std::min(q, i - newDayC[dayIndex]); j < q; j++) psi[i] += beta[j * M + regime[i]] * mean;
    
    // if 'startType' == 2, the first max(p,r,q) psi[i] are set to the mean:
    if(startType == 2 && i - newDayC[dayIndex] >= 0 && i - newDayC[dayIndex] < max_prq) psi[i] = mean;
    
    e[i] = x[i] / psi[i];
    logPsi[i] = log(psi[i]);
  }
  
  //computes the derivatives:
  if(returnIndex == 3 || returnIndex == 4){
    
    int Npara = (dist == 1) ? param.size() : param.size() + distPara.size();
    NumericMatrix d_psi_d_theta(N, Npara);
    NumericMatrix d_LL_d_theta(N, Npara);
    NumericVector score(Npara);
    NumericVector LLi(N);
    
    int dayIndex = 0;
    
    for(int i = 0; i < N; i++){
      
      if(newDayC[dayIndex + 1] == i) dayIndex++;
      
      // d_psi_d_omega:
      d_psi_d_theta(i, regime[i]) = 1;
      
      // d_psi_d_alpha:
      for(int j = 0; j < std::min(p, i - newDayC[dayIndex]); j++)
        d_psi_d_theta(i, (j+1) * (M+1) + regime[i]) = x[i - j - 1];
      for(int j = std::min(p, i - newDayC[dayIndex]); j < p; j++)
        d_psi_d_theta(i, (j+1) * (M+1) + regime[i]) = mean;
      
      // d_psi_d_nu:
      for(int j = 0; j < std::min(r, i - newDayC[dayIndex]); j++)
        d_psi_d_theta(i, (p+j+1) * (M+1) + regime[i]) = e[i - j - 1];
      for(int j = std::min(r, i - newDayC[dayIndex]); j < r; j++)
        d_psi_d_theta(i, (p+j+1) * (M+1) + regime[i]) = 1;
      
      // d_psi_d_beta:
      for(int j = 0; j < std::min(q, i - newDayC[dayIndex]); j++)
        d_psi_d_theta(i, (p+r+j+1) * (M+1) + regime[i]) = psi[i - j - 1];
      for(int j = std::min(q, i - newDayC[dayIndex]); j < q; j++)
        d_psi_d_theta(i, (p+r+j+1) * (M+1) + regime[i]) = mean;
      
      // adds "common" part:
      for(int col = 0; col < param.size(); col++){
        for(int j = 0; j < std::min(q, i - newDayC[dayIndex]); j++)
          d_psi_d_theta(i, col) += d_psi_d_theta(i - j - 1, col) * (beta[j * (M+1) + regime[i]]  
                                   - nu[j * (M+1) + regime[i]] * e[i - j - 1] / psi[i - j - 1]);
      }
      
      // if 'startType' == 2, the first max(p,q) psi[i] are set to mean (so that d_psi_d_theta = 0)
      if(startType == 2 && i - newDayC[dayIndex] >= 0 && i - newDayC[dayIndex] < max_prq){
        for(int col = 0; col < Npara; col++){
          d_psi_d_theta(i, col) = 0;
        }
      }
      
    }
    
    // computes d_LL_d_theta:
    LL = get_d_LL_d_psi_distCpp(x, psi, logPsi, e, d_psi_d_theta,
                                d_LL_d_theta, dist, distPara,
                                LLi, param.size(), score);
    
    if(returnIndex == 3 || returnIndex == 4) returnList["score"] = score; 
    if(returnIndex == 4){
      
      returnList["psi"] = psi;
      returnList["e"] = e;
      returnList["d_LL_d_theta"] = d_LL_d_theta; 
      returnList["d_psi_d_theta"] = d_psi_d_theta;
      returnList["LLi"] = LLi;
    }
    
  }
  
  // computes the log-likelihood:
  if(returnIndex != 3 && returnIndex != 4) 
    LL = getLL_distCpp(x, psi, logPsi, e, dist, distPara);
  
  returnList["LL"] = LL;
  if(returnIndex == 2){
    returnList["psi"] = psi;
    returnList["e"] = e;
  }
  
  return(returnList);
}

// [[Rcpp::export(name = ".computeScoreLSNIACD")]]
List computeScoreLSNIACD(NumericVector x,
                     NumericVector param,
                     NumericVector order,
                     double mean,
                     int dist,
                     NumericVector distPara,
                     IntegerVector newDayR,
                     int forceErrExpec,
                     int returnIndex,
                     int startType,
                     NumericVector bp){
  
  int N = x.size();
  
  NumericVector psi(N); // the conditional durations 
  NumericVector logPsi(N); // log(psi)
  NumericVector e(N); // the errors (x/psi)
  NumericVector LL(1); // log-likelihood
  
  List returnList; 
  
  // extracts the parameters:
  int p = order(0);
  int q = order(1);
  int maxPQ = std::max(p, q);
  int M = bp.size();
  double omega = param(0);
  NumericVector c; if(p > 0) c = param[Rcpp::seq(1, M + 1)];
  NumericVector alpha; if(p > 1) alpha = param[Rcpp::seq(M + 2, M + p)];
  NumericVector beta; if(q > 0) beta = param[Rcpp::seq(1 + M + p, M + p + q)];
  
  IntegerVector newDayC(newDayR.size() + 2);   
  if(newDayR.size() > 0) newDayC[Rcpp::Range(1, newDayR.size())] = newDayR - 1; 
  newDayC[newDayC.size() - 1] = N;
  
  int dayIndex = 0; double holder = 0;
  
  for(int i = 0; i < N; i++){
    
    if(newDayC[dayIndex + 1] == i) dayIndex++;
    
    logPsi[i] = omega; //adds the constant
    
    //adds the "p"-part:
    for(int j = 0; j < std::min(p, i - newDayC[dayIndex]); j++){
      holder = c[0] * e[i - j - 1];
      
      for(int k = 1; k < M + 1; k++){
        if(e[i - j - 1] > bp[k - 1]) holder += c[k] * (e[i - j - 1] - bp[k - 1]);
        else break;
      }
      
      if(j > 0) holder *= alpha[j - 1];
      
      logPsi[i] += holder;
    }
    for(int j = std::min(p, i - newDayC[dayIndex]); j < p; j++){
      holder = c[0];
      
      for(int k = 1; k < M + 1; k++){
        if(1 > bp[k - 1]) holder += c[k] * (1 - bp[k - 1]);
        else break;
      }
      
      if(j > 0) holder *= alpha[j - 1];
      
      logPsi[i] += holder;
    }
    
    //adds the "q"-part:
    for(int j = 0; j < std::min(q, i - newDayC[dayIndex]); j++) logPsi[i] += beta[j] * logPsi[i - j - 1];
    for(int j = std::min(q, i - newDayC[dayIndex]); j < q; j++) logPsi[i] += beta[j] * log(mean);
    
    // if 'startType' == 2, the first max(p,q) psi[i] are set to the mean:
    if(startType == 2 && i - newDayC[dayIndex] >= 0 && i - newDayC[dayIndex] < maxPQ) logPsi[i] = log(mean);
    
    
    psi[i] = exp(logPsi[i]);
    e[i] = x[i] / psi[i];
  }
  
  //computes the derivatives:
  if(returnIndex == 3 || returnIndex == 4){
    
    int Npara = (dist == 1) ? param.size() : param.size() + distPara.size();
    NumericMatrix d_psi_d_theta(N, Npara);
    NumericMatrix d_LL_d_theta(N, Npara);
    NumericVector score(Npara);
    NumericVector LLi(N);
    
    int dayIndex = 0;
    
    for(int i = 0; i < N; i++){
      
      if(newDayC[dayIndex + 1] == i) dayIndex++;
      
      // d_psi_d_omega:
      d_psi_d_theta(i, 0) = 1;
      
      // d_psi_d_c:
      for(int j = 0; j < std::min(p, i - newDayC[dayIndex]); j++){
        if(j == 0) d_psi_d_theta(i, 1) += e[i - j - 1]; //"k = 0"
        //^since alpha[0] is set to 1
        else d_psi_d_theta(i, 1) += alpha[j - 1] * e[i - j - 1]; //"k = 0"
      }
      for(int k = 1; k < M + 1; k++){
        for(int j = 0; j < std::min(p, i - newDayC[dayIndex]); j++){
          if(e[i - j - 1] > bp[k - 1]){
            if(j == 0) d_psi_d_theta(i, 1 + k) += (e[i - j - 1] - bp[k - 1]);
            else d_psi_d_theta(i, 1 + k) += alpha[j - 1] * (e[i - j - 1] - bp[k - 1]);
          }
          else break;
        }
      }
      for(int j = std::min(p, i - newDayC[dayIndex]); j < p; j++){
        if(j == 0) d_psi_d_theta(i, 1) += 1; //"k = 0"
        else d_psi_d_theta(i, 1) += alpha[j - 1]; //"k = 0"
      }
      for(int k = 1; k < M + 1; k++){
        for(int j = std::min(p, i - newDayC[dayIndex]); j < p; j++){
          if(1 > bp[k - 1]){
            if(j == 0) d_psi_d_theta(i, 1 + k) += 1 - bp[k - 1];
            else d_psi_d_theta(i, 1 + k) += alpha[j - 1] * (1 - bp[k - 1]);
          }
          else break;
        }
      }
      
      // d_psi_d_alpha:
      for(int j = 1; j < std::min(p, i - newDayC[dayIndex]); j++){
        d_psi_d_theta(i, M + 1 + j) = c[0] * e[i - j - 1]; //"k = 0"
        for(int k = 1; k < M + 1; k++){
          if(e[i - j - 1] > bp[k - 1])
            d_psi_d_theta(i, M + 1 + j) += c[k] * (e[i - j - 1] - bp[k - 1]);
          else break;
        }
      }
      for(int j = std::min(p, i - newDayC[dayIndex]) + 1; j < p; j++){
        d_psi_d_theta(i, M + 1 + j) = c[0]; //"k = 0"
        for(int k = 1; k < M + 1; k++){
          if(1 > bp[k - 1])
            d_psi_d_theta(i, M + 1 + j) += c[k] * (1 - bp[k - 1]);
          else break;
        }
      }
      
      // d_psi_d_beta:
      for(int j = 0; j < std::min(q, i - newDayC[dayIndex]); j++) d_psi_d_theta(i, M + p + j + 1) = logPsi[i - j - 1];
      for(int j = std::min(q, i - newDayC[dayIndex]); j < q; j++) d_psi_d_theta(i,  M + p + j + 1) = log(mean);
      
      // adds "common" part:
      for(int j = 0; j < std::min(maxPQ, i - newDayC[dayIndex]); j++){
        double factor = 0;
        
        if(j < p && j == 0){
          factor -= c[0] * e[i - j - 1]; //"k = 0"
          for(int k = 1; k < M + 1; k++){
            if(e[i - j - 1] > bp[k - 1])
              factor -= c[k] * e[i - j - 1];
            else break;
          }
        }
        if(j < p && j > 0){
          factor -= alpha[j - 1] * c[0] * e[i - j - 1]; //"k = 0"
          for(int k = 1; k < M + 1; k++){
            if(e[i - j - 1] > bp[k - 1])
              factor -= alpha[j - 1] * c[k] * e[i - j - 1];
            else break;
          }
        }
        
        if(j < q) factor += beta[j];
        
        factor /= psi[i - j - 1];
        
        for(int col = 0; col < param.size(); col++){
          d_psi_d_theta(i, col) += factor * d_psi_d_theta(i - j - 1, col);
        }
      }
      
      // multiplies by 'psi[i]':
      for(int col = 0; col < param.size(); col++) d_psi_d_theta(i, col) *= psi[i];
      
      // if 'startType' == 2, the first max(p,q) psi[i] are set to mean (so that d_psi_d_theta = 0)
      if(startType == 2 && i - newDayC[dayIndex] >= 0 && i - newDayC[dayIndex] < maxPQ){
        for(int col = 0; col < Npara; col++){
          d_psi_d_theta(i, col) = 0;
        }
      }
    }
    
    // computes d_LL_d_theta:
    LL = get_d_LL_d_psi_distCpp(x, psi, logPsi, e, d_psi_d_theta,
                                d_LL_d_theta, dist, distPara,
                                LLi, param.size(), score);
    
    if(returnIndex == 3 || returnIndex == 4) returnList["score"] = score; 
    if(returnIndex == 4){
      
      returnList["psi"] = psi;
      returnList["e"] = e;
      returnList["d_LL_d_theta"] = d_LL_d_theta; 
      returnList["d_psi_d_theta"] = d_psi_d_theta;
      returnList["LLi"] = LLi;
    }
    
  }
  
  // computes the log-likelihood:
  if(returnIndex != 3 && returnIndex != 4) 
    LL = getLL_distCpp(x, psi, logPsi, e, dist, distPara);
  
  returnList["LL"] = LL;
  if(returnIndex == 2){
    returnList["psi"] = psi;
    returnList["e"] = e;
  }
  
  return(returnList);
}

// [[Rcpp::export(name = ".computeScoreSNIACD")]]
List computeScoreSNIACD(NumericVector x,
                        NumericVector param,
                        NumericVector order,
                        double mean,
                        int dist,
                        NumericVector distPara,
                        IntegerVector newDayR,
                        int forceErrExpec,
                        int returnIndex,
                        int startType,
                        NumericVector bp){
  
  int N = x.size();
  
  NumericVector psi(N); // the conditional durations 
  NumericVector logPsi(N); // log(psi)
  NumericVector e(N); // the errors (x/psi)
  NumericVector LL(1); // log-likelihood
  
  List returnList; 
  
  // extracts the parameters:
  int p = order(0);
  int q = order(1);
  int maxPQ = std::max(p, q);
  int M = bp.size();
  double omega = param(0);
  NumericVector c; if(p > 0) c = param[Rcpp::seq(1, M + 1)];
  NumericVector alpha; if(p > 1) alpha = param[Rcpp::seq(M + 2, M + p)];
  NumericVector beta; if(q > 0) beta = param[Rcpp::seq(1 + M + p, M + p + q)];
  
  IntegerVector newDayC(newDayR.size() + 2);   
  if(newDayR.size() > 0) newDayC[Rcpp::Range(1, newDayR.size())] = newDayR - 1; 
  newDayC[newDayC.size() - 1] = N;
  
  int dayIndex = 0; double holder = 0;
  
  for(int i = 0; i < N; i++){
    
    if(newDayC[dayIndex + 1] == i) dayIndex++;
    
    psi[i] = omega; //adds the constant
    
    //adds the "p"-part:
    for(int j = 0; j < std::min(p, i - newDayC[dayIndex]); j++){
      holder = c[0] * e[i - j - 1];
      
      for(int k = 1; k < M + 1; k++){
        if(e[i - j - 1] > bp[k - 1]) holder += c[k] * (e[i - j - 1] - bp[k - 1]);
        else break;
      }
      
      if(j > 0) holder *= alpha[j - 1];
      
      psi[i] += holder;
    }
    for(int j = std::min(p, i - newDayC[dayIndex]); j < p; j++){
      holder = c[0];
      
      for(int k = 1; k < M + 1; k++){
        if(1 > bp[k - 1]) holder += c[k] * (1 - bp[k - 1]);
        else break;
      }
      
      if(j > 0) holder *= alpha[j - 1];
      
      psi[i] += holder;
    }
    
    //adds the "q"-part:
    for(int j = 0; j < std::min(q, i - newDayC[dayIndex]); j++) psi[i] += beta[j] * psi[i - j - 1];
    for(int j = std::min(q, i - newDayC[dayIndex]); j < q; j++) psi[i] += beta[j] * mean;
    
    // if 'startType' == 2, the first max(p,q) psi[i] are set to the mean:
    if(startType == 2 && i - newDayC[dayIndex] >= 0 && i - newDayC[dayIndex] < maxPQ) psi[i] = mean;
    
    
    double psi_min = .000001;
    if(psi[i] < psi_min) psi[i] = psi_min;
    
    
    e[i] = x[i]/psi[i];
    logPsi[i] = log(psi[i]);
  }
  
  //computes the derivatives:
  if(returnIndex == 3 || returnIndex == 4){
    
    int Npara = (dist == 1) ? param.size() : param.size() + distPara.size();
    NumericMatrix d_psi_d_theta(N, Npara);
    NumericMatrix d_LL_d_theta(N, Npara);
    NumericVector score(Npara);
    NumericVector LLi(N);
    
    int dayIndex = 0;
    
    for(int i = 0; i < N; i++){
      
      if(newDayC[dayIndex + 1] == i) dayIndex++;
      
      // d_psi_d_omega:
      d_psi_d_theta(i, 0) = 1;
      
      // d_psi_d_c:
      for(int j = 0; j < std::min(p, i - newDayC[dayIndex]); j++){
        if(j == 0) d_psi_d_theta(i, 1) += e[i - j - 1]; //"k = 0"
        //^since alpha[0] is set to 1
        else d_psi_d_theta(i, 1) += alpha[j - 1] * e[i - j - 1]; //"k = 0"
      }
      for(int k = 1; k < M + 1; k++){
        for(int j = 0; j < std::min(p, i - newDayC[dayIndex]); j++){
          if(e[i - j - 1] > bp[k - 1]){
            if(j == 0) d_psi_d_theta(i, 1 + k) += (e[i - j - 1] - bp[k - 1]);
            else d_psi_d_theta(i, 1 + k) += alpha[j - 1] * (e[i - j - 1] - bp[k - 1]);
          }
          else break;
        }
      }
      for(int j = std::min(p, i - newDayC[dayIndex]); j < p; j++){
        if(j == 0) d_psi_d_theta(i, 1) += 1; //"k = 0"
        else d_psi_d_theta(i, 1) += alpha[j - 1]; //"k = 0"
      }
      for(int k = 1; k < M + 1; k++){
        for(int j = std::min(p, i - newDayC[dayIndex]); j < p; j++){
          if(1 > bp[k - 1]){
            if(j == 0) d_psi_d_theta(i, 1 + k) += 1 - bp[k - 1];
            else d_psi_d_theta(i, 1 + k) += alpha[j - 1] * (1 - bp[k - 1]);
          }
          else break;
        }
      }
      
      // d_psi_d_alpha:
      for(int j = 1; j < std::min(p, i - newDayC[dayIndex]); j++){
        d_psi_d_theta(i, M + 1 + j) = c[0] * e[i - j - 1]; //"k = 0"
        for(int k = 1; k < M + 1; k++){
          if(e[i - j - 1] > bp[k - 1])
            d_psi_d_theta(i, M + 1 + j) += c[k] * (e[i - j - 1] - bp[k - 1]);
          else break;
        }
      }
      for(int j = std::min(p, i - newDayC[dayIndex]) + 1; j < p; j++){
        d_psi_d_theta(i, M + 1 + j) = c[0]; //"k = 0"
        for(int k = 1; k < M + 1; k++){
          if(1 > bp[k - 1])
            d_psi_d_theta(i, M + 1 + j) += c[k] * (1 - bp[k - 1]);
          else break;
        }
      }
      
      // d_psi_d_beta:
      for(int j = 0; j < std::min(q, i - newDayC[dayIndex]); j++) d_psi_d_theta(i, M + p + j + 1) = psi[i - j - 1];
      for(int j = std::min(q, i - newDayC[dayIndex]); j < q; j++) d_psi_d_theta(i,  M + p + j + 1) = mean;
      
      // adds "common" part:
      for(int j = 0; j < std::min(maxPQ, i - newDayC[dayIndex]); j++){
        double factor = 0;
        
        if(j < p && j == 0){
          factor -= c[0] * e[i - j - 1] / psi[i - j - 1]; //"k = 0"
          for(int k = 1; k < M + 1; k++){
            if(e[i - j - 1] > bp[k - 1])
              factor -= c[k] * e[i - j - 1] / psi[i - j - 1];
            else break;
          }
        }
        if(j < p && j > 0){
          factor -= alpha[j - 1] * c[0] * e[i - j - 1] / psi[i - j - 1]; //"k = 0"
          for(int k = 1; k < M + 1; k++){
            if(e[i - j - 1] > bp[k - 1])
              factor -= alpha[j - 1] * c[k] * e[i - j - 1] / psi[i - j - 1];
            else break;
          }
        }
        
        if(j < q) factor += beta[j];
        
        for(int col = 0; col < param.size(); col++){
          d_psi_d_theta(i, col) += factor * d_psi_d_theta(i - j - 1, col);
        }
      }
      
      // if 'startType' == 2, the first max(p,q) psi[i] are set to mean (so that d_psi_d_theta = 0)
      if(startType == 2 && i - newDayC[dayIndex] >= 0 && i - newDayC[dayIndex] < maxPQ){
        for(int col = 0; col < Npara; col++){
          d_psi_d_theta(i, col) = 0;
        }
      }
    }
    
    // computes d_LL_d_theta:
    LL = get_d_LL_d_psi_distCpp(x, psi, logPsi, e, d_psi_d_theta,
                                d_LL_d_theta, dist, distPara,
                                LLi, param.size(), score);
    
    if(returnIndex == 3 || returnIndex == 4) returnList["score"] = score; 
    if(returnIndex == 4){
      
      returnList["psi"] = psi;
      returnList["e"] = e;
      returnList["d_LL_d_theta"] = d_LL_d_theta; 
      returnList["d_psi_d_theta"] = d_psi_d_theta;
      returnList["LLi"] = LLi;
    }
    
  }
  
  // computes the log-likelihood:
  if(returnIndex != 3 && returnIndex != 4) 
    LL = getLL_distCpp(x, psi, logPsi, e, dist, distPara);
  
  returnList["LL"] = LL;
  if(returnIndex == 2){
    returnList["psi"] = psi;
    returnList["e"] = e;
  }
  
  return(returnList);
}


// [[Rcpp::export(name = ".computeScoreABACD")]]
List computeScoreABACD(NumericVector x,
                      NumericVector param,
                      NumericVector order,
                      double mean,
                      int dist,
                      NumericVector distPara,
                      IntegerVector newDayR,
                      int forceErrExpec,
                      int returnIndex,
                      int startType){
  
  int N = x.size();
  
  NumericVector psi(N); // the conditional durations 
  NumericVector psi_d1(N); // psi^d1
  NumericVector logPsi(N); // log(psi)
  NumericVector e(N); // the errors (x/psi)
  NumericVector LL(1); // log-likelihood
  
  NumericVector phi(N);
  NumericVector phi_d2(N);
  
  List returnList; 
  
  // extracts the parameters:
  int p = order(0);
  int q = order(1);
  int maxPQ = std::max(p, q);
  double omega = param(0);
  NumericVector alpha; if(p > 0) alpha = param[Rcpp::seq(1, p)];
  NumericVector beta; if(q > 0) beta = param[Rcpp::seq(1 + p, p + q)];
  double c = param(p + q + 1);
  double v = param(p + q + 2);
  double d1 = param(p + q + 3);
  double d2 = param(p + q + 4);

  IntegerVector newDayC(newDayR.size() + 2);   
  if(newDayR.size() > 0) newDayC[Rcpp::Range(1, newDayR.size())] = newDayR - 1; 
  newDayC[newDayC.size() - 1] = N;
  
  int dayIndex = 0;
  
  for(int i = 0; i < N; i++){
    
    if(newDayC[dayIndex + 1] == i) dayIndex++;
    
    psi_d1[i] = omega; //adds the constant
    
    //adds the "p"-part:
    for(int j = 0; j < std::min(p, i - newDayC[dayIndex]); j++) psi_d1[i] += alpha[j] * phi_d2[i - j - 1];
    for(int j = std::min(p, i - newDayC[dayIndex]); j < p; j++) psi_d1[i] += alpha[j] * pow(std::abs(1 - v) + c * (1 - v), d2);
    
    //adds the "q"-part:
    for(int j = 0; j < std::min(q, i - newDayC[dayIndex]); j++) psi_d1[i] += beta[j] * psi_d1[i - j - 1]; 
    for(int j = std::min(q, i - newDayC[dayIndex]); j < q; j++) psi_d1[i] += beta[j] * pow(mean, d1);
    
    // if 'startType' == 2, the first max(p,q) psi[i] are set to the mean:
    if(startType == 2 && i - newDayC[dayIndex] >= 0 && i - newDayC[dayIndex] < maxPQ) psi_d1[i] = pow(mean, d1);
    
    psi[i] = pow(psi_d1[i], 1/d1);
    e[i] = x[i]/psi[i];
    // if 'startType' == 2, the first max(p,q) residuals are set to 1:
    if(startType == 2 && i - newDayC[dayIndex] >= 0 && i - newDayC[dayIndex] < maxPQ) e[i] = 1;
    logPsi[i] = log(psi[i]);
    
    phi[i] = std::abs(e[i] - v) + c * (e[i] - v);
    phi_d2[i] = pow(phi[i], d2);
  }
  
  //computes the derivatives:
  if(returnIndex == 3 || returnIndex == 4){
    
    NumericVector logPhi(N);
    int Npara = (dist == 1) ? param.size() : param.size() + distPara.size();
    NumericMatrix d_psi_d_theta(N, Npara);
    NumericMatrix d_LL_d_theta(N, Npara);
    NumericVector score(Npara);
    NumericVector LLi(N);
    
    int dayIndex = 0;
    
    for(int i = 0; i < N; i++){
      
      if(newDayC[dayIndex + 1] == i) dayIndex++;
      
      logPhi[i] = log(phi[i]);
      
      // d_psi_d_omega:
      d_psi_d_theta(i, 0) = 1;
      // d_psi_d_alpha:
      for(int k = 1; k <= std::min(p, i - newDayC[dayIndex]); k++) 
        d_psi_d_theta(i, k) = phi_d2[i - k];
      // d_psi_d_beta:
      for(int k = 1; k <= std::min(q, i - newDayC[dayIndex]); k++) 
        d_psi_d_theta(i, p + k) = psi_d1[i - k];
      
      // d_psi_d_c
      for(int j = 0; j < std::min(p, i - newDayC[dayIndex]); j++)
        d_psi_d_theta(i, p + q + 1) += d2 * alpha[j] * phi_d2[i - j - 1] / phi[i - j - 1] *
          (e[i - j - 1] - v);
      //d_psi_d_v:
      for(int j = 0; j < std::min(p, i - newDayC[dayIndex]); j++)
        d_psi_d_theta(i, p + q + 2) -= d2 * alpha[j] * phi_d2[i - j - 1] / phi[i - j - 1] *
          (abs_prim(e[i - j - 1] - v) + c);
      
      // d_psi_d_d1:
      d_psi_d_theta(i, p + q + 3) = - psi_d1[i] * logPsi[i];
      for(int j = 0; j < std::min(q, i - newDayC[dayIndex]); j++) d_psi_d_theta(i, p + q + 3) += 
        beta[j] * psi_d1[i - j - 1] * logPsi[i - j - 1];
      // d_psi_d_d2:
      for(int j = 0; j < std::min(p, i - newDayC[dayIndex]); j++) d_psi_d_theta(i, p + q + 4) += 
        alpha[j] * phi_d2[i - j - 1] * logPhi[i - j - 1];
      
      // adds "common" part:
      for(int j = 0; j < std::min(maxPQ, i - newDayC[dayIndex]); j++){
        
        double factor = 0;
        
        if(j < p) factor -= 
          d2 * alpha[j] * phi_d2[i - j - 1] / phi[i - j - 1] * (abs_prim(e[i - j - 1] - v) + c) *
          e[i - j - 1] / psi[i - j - 1];
        if(j < q) factor += d1 * beta[j] * psi_d1[i - j - 1] / psi[i - j - 1];
        
        for(int col = 0; col < param.size(); col++) 
          d_psi_d_theta(i, col) += factor * d_psi_d_theta(i - j - 1, col);
      }
      
      // multiplies with common factor 'multi2':
      double multi2 = psi[i] / (d1 * psi_d1[i]);
      for(int col = 0; col < param.size(); col++) d_psi_d_theta(i, col) *= multi2;
      
      // if 'startType' == 2, the first max(p,q) psi[i] are set to mean (so that d_psi_d_theta = 0)
      if(startType == 2 && i - newDayC[dayIndex] >= 0 && i - newDayC[dayIndex] < maxPQ){
        for(int col = 0; col < Npara; col++){
          d_psi_d_theta(i, col) = 0;
        }
      }
    }
    
    // computes d_LL_d_theta:
    LL = get_d_LL_d_psi_distCpp(x, psi, logPsi, e, d_psi_d_theta,
                                d_LL_d_theta, dist, distPara,
                                LLi, param.size(), score);
    
    if(returnIndex == 3 || returnIndex == 4) returnList["score"] = score; 
    if(returnIndex == 4){
      
      returnList["psi"] = psi;
      returnList["e"] = e;
      returnList["d_LL_d_theta"] = d_LL_d_theta; 
      returnList["d_psi_d_theta"] = d_psi_d_theta;
      returnList["LLi"] = LLi;
    }
    
  }
  
  // computes the log-likelihood:
  if(returnIndex != 3 && returnIndex != 4) 
    LL = getLL_distCpp(x, psi, logPsi, e, dist, distPara);
  
  returnList["LL"] = LL;
  if(returnIndex == 2){
    returnList["psi"] = psi;
    returnList["e"] = e;
  }
  
  return(returnList);
}

// [[Rcpp::export(name = ".computeScoreAMACD")]]
List computeScoreAMACD(NumericVector x,
                       NumericVector param,
                       NumericVector order,
                       double mean,
                       int dist,
                       NumericVector distPara,
                       IntegerVector newDayR,
                       int forceErrExpec,
                       int returnIndex,
                       int startType){
  
  int N = x.size();
  
  NumericVector psi(N); // the conditional durations 
  NumericVector logPsi(N); // log(psi)
  NumericVector e(N); // the errors (x/psi)
  NumericVector LL(1); // log-likelihood
  
  List returnList; 
  
  // extracts the parameters:
  int p = order(0);
  int r = order(1);
  int q = order(2);
  int maxPQ = std::max(p, q); maxPQ = std::max(maxPQ, r);
  double omega = param(0);
  NumericVector alpha; if(p > 0) alpha = param[Rcpp::seq(1, p)];
  NumericVector nu; if(r > 0) nu = param[Rcpp::seq(1 + p, p + r)];
  NumericVector beta; if(q > 0) beta = param[Rcpp::seq(1 + p + r, p + r + q)];
    
  IntegerVector newDayC(newDayR.size() + 2);   
  if(newDayR.size() > 0) newDayC[Rcpp::Range(1, newDayR.size())] = newDayR - 1; 
  newDayC[newDayC.size() - 1] = N;
  
  int dayIndex = 0;
  
  for(int i = 0; i < N; i++){
    
    if(newDayC[dayIndex + 1] == i) dayIndex++;
    
    psi[i] = omega; //adds the constant
    
    //adds the "p"-part:
    for(int j = 0; j < std::min(p, i - newDayC[dayIndex]); j++) psi[i] += alpha[j] * x[i - j - 1];
    for(int j = std::min(p, i - newDayC[dayIndex]); j < p; j++) psi[i] += alpha[j] * mean;
    
    //adds the "r"-part:
    for(int j = 0; j < std::min(r, i - newDayC[dayIndex]); j++) psi[i] += nu[j] * e[i - j - 1];
    for(int j = std::min(r, i - newDayC[dayIndex]); j < r; j++) psi[i] += nu[j];
    
    //adds the "q"-part:
    for(int j = 0; j < std::min(q, i - newDayC[dayIndex]); j++) psi[i] += beta[j] * psi[i - j - 1]; 
    for(int j = std::min(q, i - newDayC[dayIndex]); j < q; j++) psi[i] += beta[j] * mean;
    
    // if 'startType' == 2, the first max(p,q) psi[i] are set to the mean:
    if(startType == 2 && i - newDayC[dayIndex] >= 0 && i - newDayC[dayIndex] < maxPQ) psi[i] = mean;
    
    logPsi[i] = log(psi[i]);
    e[i] = x[i] / psi[i];
  }
  
  //computes the derivatives:
  if(returnIndex == 3 || returnIndex == 4){
    
    int Npara = (dist == 1) ? param.size() : param.size() + distPara.size();
    NumericMatrix d_psi_d_theta(N, Npara);
    NumericMatrix d_LL_d_theta(N, Npara);
    NumericVector score(Npara);
    NumericVector LLi(N);
    
    int dayIndex = 0;
    
    for(int i = 0; i < N; i++){
      
      if(newDayC[dayIndex + 1] == i) dayIndex++;
      
      // d_psi_d_omega:
      d_psi_d_theta(i, 0) = 1;
      // d_psi_d_alpha:
      for(int j = 0; j < std::min(p, i - newDayC[dayIndex]); j++) d_psi_d_theta(i, j + 1) = x[i - j - 1];
      for(int j = std::min(p, i - newDayC[dayIndex]); j < p; j++) d_psi_d_theta(i, j + 1) = mean;
      // d_psi_d_nu:
      for(int j = 0; j < std::min(r, i - newDayC[dayIndex]); j++) d_psi_d_theta(i, p + j + 1) = e[i - j - 1];
      for(int j = std::min(r, i - newDayC[dayIndex]); j < r; j++) d_psi_d_theta(i, p + j + 1) = 1;
      // d_psi_d_beta:
      for(int j = 0; j < std::min(q, i - newDayC[dayIndex]); j++) d_psi_d_theta(i, p + r + j + 1) = psi[i - j - 1];
      for(int j = std::min(q, i - newDayC[dayIndex]); j < q; j++) d_psi_d_theta(i, p + r + j + 1) = mean;
      
      // adds "common" part:
      for(int j = 0; j < std::min(maxPQ, i - newDayC[dayIndex]); j++){
        
        double factor = 0;
        
        if(j < r) factor -= 
          nu[j] * e[i - j - 1] / psi[i - j - 1];
        if(j < q) factor += beta[j];
        
        for(int col = 0; col < param.size(); col++) 
          d_psi_d_theta(i, col) += factor * d_psi_d_theta(i - j - 1, col);
      }
      
      // if 'startType' == 2, the first max(p,q) psi[i] are set to mean (so that d_psi_d_theta = 0)
      if(startType == 2 && i - newDayC[dayIndex] >= 0 && i - newDayC[dayIndex] < maxPQ){
        for(int col = 0; col < Npara; col++){
          d_psi_d_theta(i, col) = 0;
        }
      }
    }
    
    // computes d_LL_d_theta:
    LL = get_d_LL_d_psi_distCpp(x, psi, logPsi, e, d_psi_d_theta,
                                d_LL_d_theta, dist, distPara,
                                LLi, param.size(), score);
    
    if(returnIndex == 3 || returnIndex == 4) returnList["score"] = score; 
    if(returnIndex == 4){
      
      returnList["psi"] = psi;
      returnList["e"] = e;
      returnList["d_LL_d_theta"] = d_LL_d_theta; 
      returnList["d_psi_d_theta"] = d_psi_d_theta;
      returnList["LLi"] = LLi;
    }
    
  }
  
  // computes the log-likelihood:
  if(returnIndex != 3 && returnIndex != 4) 
    LL = getLL_distCpp(x, psi, logPsi, e, dist, distPara);
  
  returnList["LL"] = LL;
  if(returnIndex == 2){
    returnList["psi"] = psi;
    returnList["e"] = e;
  }
  
  return(returnList);
}

// [[Rcpp::export(name = ".computeScoreBCACD")]]
List computeScoreBCACD(NumericVector x,
                      NumericVector param,
                      NumericVector order,
                      double mean,
                      int dist,
                      NumericVector distPara,
                      IntegerVector newDayR,
                      int forceErrExpec,
                      int returnIndex,
                      int startType){
  
  int N = x.size();
  
  NumericVector psi(N); // the conditional durations 
  NumericVector logPsi(N); // log(psi)
  NumericVector e(N); // the errors (x/psi)
  NumericVector e_d(N); // psi^d
  NumericVector LL(1); // log-likelihood
  
  List returnList; 
  
  // extracts the parameters:
  int p = order(0);
  int q = order(1);
  int maxPQ = std::max(p, q);
  double omega = param(0);
  NumericVector alpha; if(p > 0) alpha = param[Rcpp::seq(1, p)];
  NumericVector beta; if(q > 0) beta = param[Rcpp::seq(1 + p, p + q)];
  double d = param(p + q + 1);
  
  IntegerVector newDayC(newDayR.size() + 2);   
  if(newDayR.size() > 0) newDayC[Rcpp::Range(1, newDayR.size())] = newDayR - 1; 
  newDayC[newDayC.size() - 1] = N;
  
  int dayIndex = 0;
  
  for(int i = 0; i < N; i++){
    
    if(newDayC[dayIndex + 1] == i) dayIndex++;
    
    logPsi[i] = omega; //adds the constant
    
    //adds the "p"-part:
    for(int j = 0; j < std::min(p, i - newDayC[dayIndex]); j++)
      logPsi[i] += alpha[j] * e_d[i - j - 1];
    for(int j = std::min(p, i - newDayC[dayIndex]); j < p; j++)
      logPsi[i] += alpha[j];
    
    //adds the "q"-part:
    for(int j = 0; j < std::min(q, i - newDayC[dayIndex]); j++) logPsi[i] += beta[j] * logPsi[i - j - 1];
    for(int j = std::min(q, i - newDayC[dayIndex]); j < q; j++) logPsi[i] += beta[j] * log(mean);
    
    // if 'startType' == 2, the first max(p,q) psi[i] are set to the mean:
    if(startType == 2 && i - newDayC[dayIndex] >= 0 && i - newDayC[dayIndex] < maxPQ) logPsi[i] = log(mean);
    
    psi[i] = exp(logPsi[i]);
    e[i] = x[i]/psi[i];
    e_d[i] = pow(e[i], d);
  }
  
  //computes the derivatives:
  if(returnIndex == 3 || returnIndex == 4){
    
    int Npara = (dist == 1) ? param.size() : param.size() + distPara.size();
    NumericMatrix d_psi_d_theta(N, Npara);
    NumericMatrix d_LL_d_theta(N, Npara);
    NumericVector score(Npara);
    NumericVector LLi(N);
    
    int dayIndex = 0;
    
    for(int i = 0; i < N; i++){
      
      if(newDayC[dayIndex + 1] == i) dayIndex++;
      
      // d_psi_d_omega:
      d_psi_d_theta(i, 0) = 1;
      // d_psi_d_alpha:
      for(int j = 0; j < std::min(p, i - newDayC[dayIndex]); j++) d_psi_d_theta(i, j + 1) = e_d[i - j - 1];
      for(int j = std::min(p, i - newDayC[dayIndex]); j < p; j++) d_psi_d_theta(i, j + 1) = 1;
      // d_psi_d_beta:
      for(int j = 0; j < std::min(q, i - newDayC[dayIndex]); j++) d_psi_d_theta(i, p + j + 1) = logPsi[i - j - 1];
      for(int j = std::min(q, i - newDayC[dayIndex]); j < q; j++) d_psi_d_theta(i, p + j + 1) = log(mean);
      // d_psi_d_d:
      for(int j = 0; j < std::min(p, i - newDayC[dayIndex]); j++) 
        d_psi_d_theta(i, p + q + 1) += alpha[j] * e_d[i - j - 1] * log(e[i - j - 1]);

      // adds "common" part:
      for(int j = 0; j < std::min(maxPQ, i - newDayC[dayIndex]); j++){
        
        double factor = 0;
        
        if(j < p) factor -= 
          d * alpha[j] * e_d[i - j - 1];
        if(j < q) factor += beta[j];
        factor /= psi[i - j - 1];
        
        for(int col = 0; col < param.size(); col++) 
          d_psi_d_theta(i, col) += factor * d_psi_d_theta(i - j - 1, col);
      }
      
      for(int col = 0; col < param.size(); col++) 
        d_psi_d_theta(i, col) *= psi[i];
      
      // if 'startType' == 2, the first max(p,q) psi[i] are set to mean (so that d_psi_d_theta = 0)
      if(startType == 2 && i - newDayC[dayIndex] >= 0 && i - newDayC[dayIndex] < maxPQ){
        for(int col = 0; col < Npara; col++){
          d_psi_d_theta(i, col) = 0;
        }
      }
    }
    
    // computes d_LL_d_theta:
    LL = get_d_LL_d_psi_distCpp(x, psi, logPsi, e, d_psi_d_theta,
                                d_LL_d_theta, dist, distPara,
                                LLi, param.size(), score);
    
    if(returnIndex == 3 || returnIndex == 4) returnList["score"] = score; 
    if(returnIndex == 4){
      
      returnList["psi"] = psi;
      returnList["e"] = e;
      returnList["d_LL_d_theta"] = d_LL_d_theta; 
      returnList["d_psi_d_theta"] = d_psi_d_theta;
      returnList["LLi"] = LLi;
    }
    
  }
  
  // computes the log-likelihood:
  if(returnIndex != 3 && returnIndex != 4) 
    LL = getLL_distCpp(x, psi, logPsi, e, dist, distPara);
  
  returnList["LL"] = LL;
  if(returnIndex == 2){
    returnList["psi"] = psi;
    returnList["e"] = e;
  }
  
  return(returnList);
}


// [[Rcpp::export(name = ".computeScoreBACD")]]
List computeScoreBACD(NumericVector x,
                       NumericVector param,
                       NumericVector order,
                       double mean,
                       int dist,
                       NumericVector distPara,
                       IntegerVector newDayR,
                       int forceErrExpec,
                       int returnIndex,
                       int startType){
  
  int N = x.size();
  
  NumericVector psi(N); // the conditional durations 
  NumericVector psi_d1(N); // psi^d1
  NumericVector logPsi(N); // log(psi)
  NumericVector e(N); // the errors (x/psi)
  NumericVector e_d2(N); // psi^d2
  NumericVector LL(1); // log-likelihood
  
  List returnList; 
  
  // extracts the parameters:
  int p = order(0);
  int q = order(1);
  int maxPQ = std::max(p, q);
  double omega = param(0);
  NumericVector alpha; if(p > 0) alpha = param[Rcpp::seq(1, p)];
  NumericVector beta; if(q > 0) beta = param[Rcpp::seq(1 + p, p + q)];
  double d1 = param(p + q + 1);
  double d2 = param(p + q + 2);
  
  IntegerVector newDayC(newDayR.size() + 2);   
  if(newDayR.size() > 0) newDayC[Rcpp::Range(1, newDayR.size())] = newDayR - 1; 
  newDayC[newDayC.size() - 1] = N;
  
  int dayIndex = 0;
  
  for(int i = 0; i < N; i++){
    
    if(newDayC[dayIndex + 1] == i) dayIndex++;
    
    psi_d1[i] = omega; //adds the constant
    
    //adds the "p"-part:
    for(int j = 0; j < std::min(p, i - newDayC[dayIndex]); j++)
      psi_d1[i] += alpha[j] * e_d2[i - j - 1];
    for(int j = std::min(p, i - newDayC[dayIndex]); j < p; j++)
      psi_d1[i] += alpha[j];
    
    //adds the "q"-part:
    for(int j = 0; j < std::min(q, i - newDayC[dayIndex]); j++) psi_d1[i] += beta[j] * psi_d1[i - j - 1];
    for(int j = std::min(q, i - newDayC[dayIndex]); j < q; j++) psi_d1[i] += beta[j] * pow(mean, d1);
    
    // if 'startType' == 2, the first max(p,q) psi[i] are set to the mean:
    if(startType == 2 && i - newDayC[dayIndex] >= 0 && i - newDayC[dayIndex] < maxPQ) psi_d1[i] = pow(mean, d1);
    
    
    // double psiMin = .0001, powPsiMin =  pow(psiMin, d1);
    // if(psi_d1[i] < powPsiMin) psi[i] = psiMin;
    // else psi[i] = pow(psi_d1[i], 1/d1);
    psi[i] = pow(psi_d1[i], 1/d1);
    logPsi[i] = log(psi[i]);
    e[i] = x[i] / psi[i];
    e_d2[i] = pow(e[i], d2);
  }
  
  //computes the derivatives:
  if(returnIndex == 3 || returnIndex == 4){
    
    int Npara = (dist == 1) ? param.size() : param.size() + distPara.size();
    NumericMatrix d_psi_d_theta(N, Npara);
    NumericMatrix d_LL_d_theta(N, Npara);
    NumericVector score(Npara);
    NumericVector LLi(N);
    
    int dayIndex = 0;
    
    for(int i = 0; i < N; i++){
      
      if(newDayC[dayIndex + 1] == i) dayIndex++;
      
      // d_psi_d_omega:
      d_psi_d_theta(i, 0) = 1;
      // d_psi_d_alpha:
      for(int j = 0; j < std::min(p, i - newDayC[dayIndex]); j++) d_psi_d_theta(i, j + 1) = e_d2[i - j - 1];
      for(int j = std::min(p, i - newDayC[dayIndex]); j < p; j++) d_psi_d_theta(i, j + 1) = 1;
      // d_psi_d_beta:
      for(int j = 0; j < std::min(q, i - newDayC[dayIndex]); j++) d_psi_d_theta(i, p + j + 1) = psi_d1[i - j - 1];
      for(int j = std::min(q, i - newDayC[dayIndex]); j < q; j++) d_psi_d_theta(i, p + j + 1) = pow(mean, d1);
      // d_psi_d_d1:
      d_psi_d_theta(i, p + q + 1) = - psi_d1[i] * logPsi[i];
      for(int j = 0; j < std::min(q, i - newDayC[dayIndex]); j++) 
        d_psi_d_theta(i, p + q + 1) += beta[j] * psi_d1[i - j - 1] * logPsi[i - j - 1];
      // d_psi_d_d2:
      for(int j = 0; j < std::min(p, i - newDayC[dayIndex]); j++) 
        d_psi_d_theta(i, p + q + 2) += alpha[j] * e_d2[i - j - 1] * log(e[i - j - 1]);
      
      // adds "common" part:
      for(int j = 0; j < std::min(maxPQ, i - newDayC[dayIndex]); j++){
        
        double factor = 0;
        
        if(j < p) factor -= 
          d2 * alpha[j] * e_d2[i - j - 1];
        if(j < q) factor += d1 * beta[j] * psi_d1[i - j - 1];
        factor /= psi[i - j - 1];
        
        for(int col = 0; col < param.size(); col++) 
          d_psi_d_theta(i, col) += factor * d_psi_d_theta(i - j - 1, col);
      }
      
      for(int col = 0; col < param.size(); col++) 
        d_psi_d_theta(i, col) *= psi[i] / (d1 * psi_d1[i]);
      
      // if 'startType' == 2, the first max(p,q) psi[i] are set to mean (so that d_psi_d_theta = 0)
      if(startType == 2 && i - newDayC[dayIndex] >= 0 && i - newDayC[dayIndex] < maxPQ){
        for(int col = 0; col < Npara; col++){
          d_psi_d_theta(i, col) = 0;
        }
      }
    }
    
    // computes d_LL_d_theta:
    LL = get_d_LL_d_psi_distCpp(x, psi, logPsi, e, d_psi_d_theta,
                                d_LL_d_theta, dist, distPara,
                                LLi, param.size(), score);
    
    if(returnIndex == 3 || returnIndex == 4) returnList["score"] = score; 
    if(returnIndex == 4){
      
      returnList["psi"] = psi;
      returnList["e"] = e;
      returnList["d_LL_d_theta"] = d_LL_d_theta; 
      returnList["d_psi_d_theta"] = d_psi_d_theta;
      returnList["LLi"] = LLi;
    }
    
  }
  
  // computes the log-likelihood:
  if(returnIndex != 3 && returnIndex != 4) 
    LL = getLL_distCpp(x, psi, logPsi, e, dist, distPara);
  
  returnList["LL"] = LL;
  if(returnIndex == 2){
    returnList["psi"] = psi;
    returnList["e"] = e;
  }
  
  return(returnList);
}

// [[Rcpp::export(name = ".computeScoreEXACD")]]
List computeScoreEXACD(NumericVector x,
                       NumericVector param,
                       NumericVector order,
                       double mean,
                       int dist,
                       NumericVector distPara,
                       IntegerVector newDayR,
                       int forceErrExpec,
                       int returnIndex,
                       int startType){
  
  int N = x.size();
  
  NumericVector psi(N); // the conditional durations 
  NumericVector logPsi(N); // log(psi)
  NumericVector e(N); // the errors (x/psi)
  NumericVector LL(1); // log-likelihood
  
  List returnList; 
  
  // extracts the parameters:
  int p = order(0);
  int q = order(1);
  int maxPQ = std::max(p, q);
  double omega = param(0);
  NumericVector alpha; if(p > 0) alpha = param[Rcpp::seq(1, p)];
  NumericVector delta; if(p > 0) delta = param[Rcpp::seq(1 + p, 2 * p)];
  NumericVector beta; if(q > 0) beta = param[Rcpp::seq(1 + 2 * p, 2 * p + q)];
  
  IntegerVector newDayC(newDayR.size() + 2);   
  if(newDayR.size() > 0) newDayC[Rcpp::Range(1, newDayR.size())] = newDayR - 1; 
  newDayC[newDayC.size() - 1] = N;
  
  int dayIndex = 0;
  
  for(int i = 0; i < N; i++){
    
    if(newDayC[dayIndex + 1] == i) dayIndex++;
    
    logPsi[i] = omega; //adds the constant

    //adds the "p"-part:
    for(int j = 0; j < std::min(p, i - newDayC[dayIndex]); j++)
      logPsi[i] += alpha[j] * e[i - j - 1] + delta[j] * std::abs(e[i - j - 1] - 1);
    for(int j = std::min(p, i - newDayC[dayIndex]); j < p; j++)
      psi[i] += alpha[j];
    
    //adds the "q"-part:
    for(int j = 0; j < std::min(q, i - newDayC[dayIndex]); j++) logPsi[i] += beta[j] * logPsi[i - j - 1]; 
    for(int j = std::min(q, i - newDayC[dayIndex]); j < q; j++) logPsi[i] += beta[j] * log(mean);
    
    // if 'startType' == 2, the first max(p,q) psi[i] are set to the mean:
    if(startType == 2 && i - newDayC[dayIndex] >= 0 && i - newDayC[dayIndex] < maxPQ) logPsi[i] = log(mean);
    
    psi[i] = exp(logPsi[i]);
    e[i] = x[i]/psi[i];
  }
  
  //computes the derivatives:
  if(returnIndex == 3 || returnIndex == 4){
    
    int Npara = (dist == 1) ? param.size() : param.size() + distPara.size();
    NumericMatrix d_psi_d_theta(N, Npara);
    NumericMatrix d_LL_d_theta(N, Npara);
    NumericVector score(Npara);
    NumericVector LLi(N);
    
    int dayIndex = 0;
    
    for(int i = 0; i < N; i++){
      
      if(newDayC[dayIndex + 1] == i) dayIndex++;
      
      // d_psi_d_omega:
      d_psi_d_theta(i, 0) = 1;
      // d_psi_d_alpha:
      for(int j = 0; j < std::min(p, i - newDayC[dayIndex]); j++) d_psi_d_theta(i, j + 1) = e[i - j - 1];
      for(int j = std::min(p, i - newDayC[dayIndex]); j < p; j++) d_psi_d_theta(i, j + 1) = 1;
      // d_psi_d_delta:
      for(int j = 0; j < std::min(p, i - newDayC[dayIndex]); j++) d_psi_d_theta(i, p + j + 1) = std::abs(e[i - j - 1] - 1);
      // d_psi_d_beta:
      for(int j = 0; j < std::min(q, i - newDayC[dayIndex]); j++) d_psi_d_theta(i, 2 * p + j + 1) = logPsi[i - j - 1];
      for(int j = std::min(q, i - newDayC[dayIndex]); j < q; j++) d_psi_d_theta(i, 2 * p + j + 1) = log(mean);
      
      
      // adds "common" part:
      for(int j = 0; j < std::min(maxPQ, i - newDayC[dayIndex]); j++){
        
        double factor = 0;
        
        if(j < p) factor -= (alpha[j] + delta[j] * abs_prim(e[i - j - 1] - 1)) * e[i - j - 1];
        if(j < q) factor += beta[j];
        factor /= psi[i - j - 1];
          
        for(int col = 0; col < param.size(); col++) 
          d_psi_d_theta(i, col) += factor * d_psi_d_theta(i - j - 1, col);
      }
      
      for(int col = 0; col < param.size(); col++) 
        d_psi_d_theta(i, col) *= psi[i];
      
      // if 'startType' == 2, the first max(p,q) psi[i] are set to mean (so that d_psi_d_theta = 0)
      if(startType == 2 && i - newDayC[dayIndex] >= 0 && i - newDayC[dayIndex] < maxPQ){
        for(int col = 0; col < Npara; col++){
          d_psi_d_theta(i, col) = 0;
        }
      }
    }
    
    // computes d_LL_d_theta:
    LL = get_d_LL_d_psi_distCpp(x, psi, logPsi, e, d_psi_d_theta,
                                d_LL_d_theta, dist, distPara,
                                LLi, param.size(), score);
    
    if(returnIndex == 3 || returnIndex == 4) returnList["score"] = score; 
    if(returnIndex == 4){
      
      returnList["psi"] = psi;
      returnList["e"] = e;
      returnList["d_LL_d_theta"] = d_LL_d_theta; 
      returnList["d_psi_d_theta"] = d_psi_d_theta;
      returnList["LLi"] = LLi;
    }
    
  }
  
  // computes the log-likelihood:
  if(returnIndex != 3 && returnIndex != 4) 
    LL = getLL_distCpp(x, psi, logPsi, e, dist, distPara);
  
  returnList["LL"] = LL;
  if(returnIndex == 2){
    returnList["psi"] = psi;
    returnList["e"] = e;
  }
  
  return(returnList);
}

// [[Rcpp::export(name = ".computeScoreLACD2")]]
List computeScoreLACD2(NumericVector x,
                       NumericVector param,
                       NumericVector order,
                       double mean,
                       int dist,
                       NumericVector distPara,
                       IntegerVector newDayR,
                       int forceErrExpec,
                       int returnIndex,
                       int startType){
  
  int N = x.size();
  
  NumericVector psi(N); // the conditional durations 
  NumericVector logPsi(N); // log(psi)
  NumericVector e(N); // the errors (x/psi)
  NumericVector LL(1); // log-likelihood
  
  List returnList; 
  
  // extracts the parameters:
  int p = order(0);
  int q = order(1);
  int maxPQ = std::max(p, q);
  double omega = param(0);
  NumericVector alpha; if(p > 0) alpha = param[Rcpp::seq(1, p)];
  NumericVector beta; if(q > 0) beta = param[Rcpp::seq(1 + p, p + q)];
  
  IntegerVector newDayC(newDayR.size() + 2);   
  if(newDayR.size() > 0) newDayC[Rcpp::Range(1, newDayR.size())] = newDayR - 1; 
  newDayC[newDayC.size() - 1] = N;
  
  int dayIndex = 0;
  
  for(int i = 0; i < N; i++){
    
    if(newDayC[dayIndex + 1] == i) dayIndex++;
    
    logPsi[i] = omega; //adds the constant
    
    //adds the "p"-part:
    for(int j = 0; j < std::min(p, i - newDayC[dayIndex]); j++) logPsi[i] += alpha[j] * e[i - j - 1];
    for(int j = std::min(p, i - newDayC[dayIndex]); j < p; j++) psi[i] += alpha[j];
    
    //adds the "q"-part:
    for(int j = 0; j < std::min(q, i - newDayC[dayIndex]); j++) logPsi[i] += beta[j] * logPsi[i - j - 1]; 
    for(int j = std::min(q, i - newDayC[dayIndex]); j < q; j++) logPsi[i] += beta[j] * log(mean);
    
    // if 'startType' == 2, the first max(p,q) psi[i] are set to the mean:
    if(startType == 2 && i - newDayC[dayIndex] >= 0 && i - newDayC[dayIndex] < maxPQ) logPsi[i] = log(mean);
    
    psi[i] = exp(logPsi[i]);
    e[i] = x[i]/psi[i];
  }
  
  //computes the derivatives:
  if(returnIndex == 3 || returnIndex == 4){
    
    int Npara = (dist == 1) ? param.size() : param.size() + distPara.size();
    NumericMatrix d_psi_d_theta(N, Npara);
    NumericMatrix d_LL_d_theta(N, Npara);
    NumericVector score(Npara);
    NumericVector LLi(N);
    
    int dayIndex = 0;
    
    for(int i = 0; i < N; i++){
      
      if(newDayC[dayIndex + 1] == i) dayIndex++;
      
      // d_psi_d_omega:
      d_psi_d_theta(i, 0) = 1;
      // d_psi_d_alpha:
      for(int j = 0; j < std::min(p, i - newDayC[dayIndex]); j++) d_psi_d_theta(i, j + 1) = e[i - j - 1];
      for(int j = std::min(p, i - newDayC[dayIndex]); j < p; j++) d_psi_d_theta(i, j + 1) = 1;
      // d_psi_d_beta:
      for(int j = 0; j < std::min(q, i - newDayC[dayIndex]); j++) d_psi_d_theta(i, p + j + 1) = logPsi[i - j - 1];
      for(int j = std::min(q, i - newDayC[dayIndex]); j < q; j++) d_psi_d_theta(i,  p + j + 1) = log(mean);
      
      
      // adds "common" part:
      for(int col = 0; col < param.size(); col++){
        for(int j = 1; j <= std::min(p, i - newDayC[dayIndex]); j++)
          d_psi_d_theta(i, col) -= alpha[j - 1] * e[i - j] / psi[i - j] * d_psi_d_theta(i - j, col);
        for(int j = 1; j <= std::min(q, i - newDayC[dayIndex]); j++)
          d_psi_d_theta(i, col) += beta[j - 1] / psi[i - j] * d_psi_d_theta(i - j, col);
        
        d_psi_d_theta(i, col) *= psi[i];
      }
      
      // if 'startType' == 2, the first max(p,q) psi[i] are set to mean (so that d_psi_d_theta = 0)
      if(startType == 2 && i - newDayC[dayIndex] >= 0 && i - newDayC[dayIndex] < maxPQ){
        for(int col = 0; col < Npara; col++){
          d_psi_d_theta(i, col) = 0;
        }
      }
    }
    
    // computes d_LL_d_theta:
    LL = get_d_LL_d_psi_distCpp(x, psi, logPsi, e, d_psi_d_theta,
                                d_LL_d_theta, dist, distPara,
                                LLi, param.size(), score);
    
    if(returnIndex == 3 || returnIndex == 4) returnList["score"] = score; 
    if(returnIndex == 4){
      
      returnList["psi"] = psi;
      returnList["e"] = e;
      returnList["d_LL_d_theta"] = d_LL_d_theta; 
      returnList["d_psi_d_theta"] = d_psi_d_theta;
      returnList["LLi"] = LLi;
    }
    
  }
  
  // computes the log-likelihood:
  if(returnIndex != 3 && returnIndex != 4) 
    LL = getLL_distCpp(x, psi, logPsi, e, dist, distPara);
  
  returnList["LL"] = LL;
  if(returnIndex == 2){
    returnList["psi"] = psi;
    returnList["e"] = e;
  }
  
  return(returnList);
}

// [[Rcpp::export(name = ".computeScoreLACD1")]]
List computeScoreLACD1(NumericVector x,
                     NumericVector param,
                     NumericVector order,
                     double mean,
                     int dist,
                     NumericVector distPara,
                     IntegerVector newDayR,
                     int forceErrExpec,
                     int returnIndex,
                     int startType){
  
  int N = x.size();
  
  NumericVector psi(N); // the conditional durations 
  NumericVector logPsi(N); // log(psi)
  NumericVector e(N); // the errors (x/psi)
  NumericVector LL(1); // log-likelihood
  
  List returnList; 
  
  // extracts the parameters:
  int p = order(0);
  int q = order(1);
  int maxPQ = std::max(p, q);
  double omega = param(0);
  NumericVector alpha; if(p > 0) alpha = param[Rcpp::seq(1, p)];
  NumericVector beta; if(q > 0) beta = param[Rcpp::seq(1 + p, p + q)];
  
  IntegerVector newDayC(newDayR.size() + 2);   
  if(newDayR.size() > 0) newDayC[Rcpp::Range(1, newDayR.size())] = newDayR - 1; 
  newDayC[newDayC.size() - 1] = N;
  
  int dayIndex = 0;
  
  for(int i = 0; i < N; i++){
    
    if(newDayC[dayIndex + 1] == i) dayIndex++;
    
    logPsi[i] = omega; //adds the constant
    
    //adds the "p"-part:
    for(int j = 0; j < std::min(p, i - newDayC[dayIndex]); j++) logPsi[i] += alpha[j] * log(e[i - j - 1]);
    
    //adds the "q"-part:
    for(int j = 0; j < std::min(q, i - newDayC[dayIndex]); j++) logPsi[i] += beta[j] * logPsi[i - j - 1]; 
    for(int j = std::min(q, i - newDayC[dayIndex]); j < q; j++) logPsi[i] += beta[j] * log(mean);
    
    // if 'startType' == 2, the first max(p,q) psi[i] are set to the mean:
    if(startType == 2 && i - newDayC[dayIndex] >= 0 && i - newDayC[dayIndex] < maxPQ) logPsi[i] = log(mean);
    
    psi[i] = exp(logPsi[i]);
    e[i] = x[i]/psi[i];
  }
  
  //computes the derivatives:
  if(returnIndex == 3 || returnIndex == 4){
    
    int Npara = (dist == 1) ? param.size() : param.size() + distPara.size();
    NumericMatrix d_psi_d_theta(N, Npara);
    NumericMatrix d_LL_d_theta(N, Npara);
    NumericVector score(Npara);
    NumericVector LLi(N);
    
    int dayIndex = 0;
    
    for(int i = 0; i < N; i++){
      
      if(newDayC[dayIndex + 1] == i) dayIndex++;
      
      // d_psi_d_omega:
      d_psi_d_theta(i, 0) = 1;
      // d_psi_d_alpha:
      for(int j = 0; j < std::min(p, i - newDayC[dayIndex]); j++) d_psi_d_theta(i, j + 1) = log(e[i - j - 1]);
      // d_psi_d_beta:
      for(int j = 0; j < std::min(q, i - newDayC[dayIndex]); j++) d_psi_d_theta(i, p + j + 1) = logPsi[i - j - 1];
      for(int j = std::min(q, i - newDayC[dayIndex]); j < q; j++) d_psi_d_theta(i,  p + j + 1) = log(mean);
      
      // adds "common" parts:
      for(int col = 0; col < param.size(); col++){
        for(int j = 1; j <= std::min(p, i - newDayC[dayIndex]); j++)
          d_psi_d_theta(i, col) -= alpha[j - 1] / psi[i - j] * d_psi_d_theta(i - j, col);
        for(int j = 1; j <= std::min(q, i - newDayC[dayIndex]); j++)
          d_psi_d_theta(i, col) += beta[j - 1] / psi[i - j] * d_psi_d_theta(i - j, col);
        
        d_psi_d_theta(i, col) *= psi[i];
      }
      
      // if 'startType' == 2, the first max(p,q) psi[i] are set to mean (so that d_psi_d_theta = 0)
      if(startType == 2 && i - newDayC[dayIndex] >= 0 && i - newDayC[dayIndex] < maxPQ){
        for(int col = 0; col < Npara; col++){
          d_psi_d_theta(i, col) = 0;
        }
      }
    }
    
    // computes d_LL_d_theta:
    LL = get_d_LL_d_psi_distCpp(x, psi, logPsi, e, d_psi_d_theta,
                                d_LL_d_theta, dist, distPara,
                                LLi, param.size(), score);
    
    if(returnIndex == 3 || returnIndex == 4) returnList["score"] = score; 
    if(returnIndex == 4){
      
      returnList["psi"] = psi;
      returnList["e"] = e;
      returnList["d_LL_d_theta"] = d_LL_d_theta; 
      returnList["d_psi_d_theta"] = d_psi_d_theta;
      returnList["LLi"] = LLi;
    }
    
  }
  
  // computes the log-likelihood:
  if(returnIndex != 3 && returnIndex != 4) 
    LL = getLL_distCpp(x, psi, logPsi, e, dist, distPara);
  
  returnList["LL"] = LL;
  if(returnIndex == 2){
    returnList["psi"] = psi;
    returnList["e"] = e;
  }
  
  return(returnList);
}

// [[Rcpp::export(name = ".computeScoreACD")]]
List computeScoreACD(NumericVector x,
                      NumericVector param,
                      NumericVector order,
                      double mean,
                      int dist,
                      NumericVector distPara,
                      IntegerVector newDayR,
                      int forceErrExpec,
                      int returnIndex,
                      int startType){
  
  int N = x.size();
  
  NumericVector psi(N); // the conditional durations 
  NumericVector logPsi(N); // log(psi)
  NumericVector e(N); // the errors (x/psi)
  NumericVector LL(1); // log-likelihood
  
  List returnList; 
  
  // extracts the parameters:
  int p = order(0);
  int q = order(1);
  int maxPQ = std::max(p, q);
  double omega = param(0);
  NumericVector alpha; if(p > 0) alpha = param[Rcpp::seq(1, p)];
  NumericVector beta; if(q > 0) beta = param[Rcpp::seq(1 + p, p + q)];
    
  IntegerVector newDayC(newDayR.size() + 2);   
  if(newDayR.size() > 0) newDayC[Rcpp::Range(1, newDayR.size())] = newDayR - 1; 
  newDayC[newDayC.size() - 1] = N;
    
  int dayIndex = 0;  
  
  for(int i = 0; i < N; i++){
	  
	if(dayIndex + 1 > newDayC.size() - 1) 
	  Rcout << "i: " << i << ", dayIndex: " << dayIndex << std::endl;
    
    if(newDayC[dayIndex + 1] == i) dayIndex++;
    
    psi[i] = omega; //adds the constant
    
    //adds the "p"-part:
    for(int j = 0; j < std::min(p, i - newDayC[dayIndex]); j++) psi[i] += alpha[j] * x[i - j - 1];
    for(int j = std::min(p, i - newDayC[dayIndex]); j < p; j++) psi[i] += alpha[j] * mean;
    
    //adds the "q"-part:
    for(int j = 0; j < std::min(q, i - newDayC[dayIndex]); j++) psi[i] += beta[j] * psi[i - j - 1]; 
    for(int j = std::min(q, i - newDayC[dayIndex]); j < q; j++) psi[i] += beta[j] * mean;
    
    // if 'startType' == 2, the first max(p,q) psi[i] are set to the mean:
    if(startType == 2 && (i - newDayC[dayIndex]) >= 0 && (i - newDayC[dayIndex]) < maxPQ) psi[i] = mean;
    
    e[i] = x[i]/psi[i];
    logPsi[i] = log(psi[i]);
  }
  
  //computes the derivatives:
  if(returnIndex == 3 || returnIndex == 4){
    
    int Npara = (dist == 1) ? param.size() : param.size() + distPara.size();
    NumericMatrix d_psi_d_theta(N, Npara);
    NumericMatrix d_LL_d_theta(N, Npara);
    NumericVector score(Npara);
    NumericVector LLi(N);
    
    int dayIndex = 0;
    
    for(int i = 0; i < N; i++){
      
      if(newDayC[dayIndex + 1] == i) dayIndex++;
      
      // d_psi_d_omega:
      d_psi_d_theta(i, 0) = 1;
      // d_psi_d_alpha:
      for(int j = 0; j < std::min(p, i - newDayC[dayIndex]); j++) d_psi_d_theta(i, j + 1) = x[i - j - 1];
      for(int j = std::min(p, i - newDayC[dayIndex]); j < p; j++) d_psi_d_theta(i, j + 1) = mean;
      // d_psi_d_beta:
      for(int j = 0; j < std::min(q, i - newDayC[dayIndex]); j++) d_psi_d_theta(i, p + j + 1) = psi[i - j - 1];
      for(int j = std::min(q, i - newDayC[dayIndex]); j < q; j++) d_psi_d_theta(i,  p + j + 1) = mean;
      
      // adds "common" part:
      for(int col = 0; col < param.size(); col++){
        for(int j = 1; j <= std::min(q, i - newDayC[dayIndex]); j++)
          d_psi_d_theta(i, col) += beta[j - 1] * d_psi_d_theta(i - j, col);
      }
      
      // if 'startType' == 2, the first max(p,q) psi[i] are set to mean (so that d_psi_d_theta = 0)
      if(startType == 2 && i - newDayC[dayIndex] >= 0 && i - newDayC[dayIndex] < maxPQ){
        for(int col = 0; col < Npara; col++){
          d_psi_d_theta(i, col) = 0;
        }
      }
    }
    
    // computes d_LL_d_theta:
    LL = get_d_LL_d_psi_distCpp(x, psi, logPsi, e, d_psi_d_theta,
                                d_LL_d_theta, dist, distPara,
                                LLi, param.size(), score);
    
    if(returnIndex == 3 || returnIndex == 4) returnList["score"] = score; 
    if(returnIndex == 4){
      
      returnList["psi"] = psi;
      returnList["e"] = e;
      returnList["d_LL_d_theta"] = d_LL_d_theta; 
      returnList["d_psi_d_theta"] = d_psi_d_theta;
      returnList["LLi"] = LLi;
    }
    
  }
  
  // computes the log-likelihood:
  if(returnIndex != 3 && returnIndex != 4) 
    LL = getLL_distCpp(x, psi, logPsi, e, dist, distPara);
  
  returnList["LL"] = LL;
  if(returnIndex == 2){
    returnList["psi"] = psi;
    returnList["e"] = e;
  }
  
  return(returnList);
}


// [[Rcpp::export(name = ".computeScoreAACD")]]
List computeScoreAACD(NumericVector x,
                      NumericVector param,
                      NumericVector order,
                      double mean,
                      int dist,
                      NumericVector distPara,
                      IntegerVector newDayR,
                      int forceErrExpec,
                      int returnIndex,
                      int startType){
  
  int N = x.size();
  
  NumericVector psi(N); // the conditional durations 
  NumericVector psi_d1(N); // psi^d1
  NumericVector logPsi(N); // log(psi)
  NumericVector e(N); // the errors (x/psi)
  NumericVector LL(1); // log-likelihood
  
  NumericVector phi(N);
  NumericVector phi_d2(N);
  
  List returnList; 
  
  // extracts the parameters:
  int p = order(0);
  int q = order(1);
  int maxPQ = std::max(p, q);
  double omega = param(0);
  NumericVector alpha; if(p > 0) alpha = param[Rcpp::seq(1, p)];
  NumericVector beta; if(q > 0) beta = param[Rcpp::seq(1 + p, p + q)];
  double c = param(p + q + 1);
  double v = param(p + q + 2);
  double d1 = param(p + q + 3);
  double d2 = param(p + q + 4);
  
  IntegerVector newDayC(newDayR.size() + 2);   
  if(newDayR.size() > 0) newDayC[Rcpp::Range(1, newDayR.size())] = newDayR - 1; 
  newDayC[newDayC.size() - 1] = N;
  
  int dayIndex = 0;
  
  for(int i = 0; i < N; i++){
    
    if(newDayC[dayIndex + 1] == i) dayIndex++;
    
    psi_d1[i] = omega; //adds the constant
    
    //adds the "p"-part:
    for(int j = 0; j < std::min(p, i - newDayC[dayIndex]); j++) psi_d1[i] += alpha[j] * psi_d1[i - j - 1] * phi_d2[i - j - 1];
    for(int j = std::min(p, i - newDayC[dayIndex]); j < p; j++) psi_d1[i] += alpha[j] * pow(mean, d2) * pow(std::abs(1 - v) + c * (1 - v), d2);
    
    //adds the "q"-part:
    for(int j = 0; j < std::min(q, i - newDayC[dayIndex]); j++) psi_d1[i] += beta[j] * psi_d1[i - j - 1]; 
    for(int j = std::min(q, i - newDayC[dayIndex]); j < q; j++) psi_d1[i] += beta[j] * pow(mean, d1);
    
    // if 'startType' == 2, the first max(p,q) psi[i] are set to the mean:
    if(startType == 2 && i - newDayC[dayIndex] >= 0 && i - newDayC[dayIndex] < maxPQ) psi_d1[i] = pow(mean, d1);
    
    psi[i] = pow(psi_d1[i], 1/d1);
    e[i] = x[i]/psi[i];
    // if 'startType' == 2, the first max(p,q) residuals are set to 1:
    if(startType == 2 && i - newDayC[dayIndex] >= 0 && i - newDayC[dayIndex] < maxPQ) e[i] = 1;
	
    logPsi[i] = log(psi[i]);
    
    phi[i] = std::abs(e[i] - v) + c * (e[i] - v);
    phi_d2[i] = pow(phi[i], d2);
  }
  
  //computes the derivatives:
  if(returnIndex == 3 || returnIndex == 4){
    
    NumericVector logPhi(N);
    int Npara = (dist == 1) ? param.size() : param.size() + distPara.size();
    NumericMatrix d_psi_d_theta(N, Npara);
    NumericMatrix d_LL_d_theta(N, Npara);
    NumericVector score(Npara);
    NumericVector LLi(N);
	std::vector<double> multi1(maxPQ);
	double multi2 = 0;
    
    int dayIndex = 0;
    
    for(int i = 0; i < N; i++){
      
      if(newDayC[dayIndex + 1] == i) dayIndex++;
      
      logPhi[i] = log(phi[i]);
      
      // common multiplier for all derivatives:
      for(int col = 0; col < maxPQ; col++) multi1[col] = 0;
      for(int j = 0; j < std::min(p, i - newDayC[dayIndex]); j++) multi1[j] += d1 * alpha[j] * psi_d1[i - j - 1] / psi[i - j - 1] *
        phi_d2[i - j - 1] -
        d2 * alpha[j] * phi_d2[i - j - 1] / phi[i - j - 1] * (abs_prim(e[i - j - 1] - v) + c) * 
        e[i - j - 1] * psi_d1[i - j - 1] / psi[i - j - 1];
      for(int j = 0; j < std::min(q, i - newDayC[dayIndex]); j++) multi1[j] += d1 * beta[j] * psi_d1[i - j - 1] / psi[i - j - 1];
      
      // d_psi_d_omega:
      d_psi_d_theta(i, 0) = 1;
      for(int j = 0; j < std::min(maxPQ, i - newDayC[dayIndex]); j++) d_psi_d_theta(i, 0) += multi1[j] * d_psi_d_theta(i - j - 1, 0);
      
      // d_psi_d_alpha:
      for(int k = 1; k <= std::min(p, i - newDayC[dayIndex]); k++) 
        d_psi_d_theta(i, k) = psi_d1[i - k] * phi_d2[i - k];
      for(int k = 1; k <= p; k++){
        for(int j = 0; j < std::min(maxPQ, i - newDayC[dayIndex]); j++) 
          d_psi_d_theta(i, k) += multi1[j] * d_psi_d_theta(i - j - 1, k);
      }
      
      // d_psi_d_beta:
      for(int k = 1; k <= std::min(q, i - newDayC[dayIndex]); k++) 
        d_psi_d_theta(i, p + k) = psi_d1[i - k];
      for(int k = 1; k <= q; k++){
        for(int j = 0; j < std::min(maxPQ, i - newDayC[dayIndex]); j++) 
          d_psi_d_theta(i, p + k) += multi1[j] * d_psi_d_theta(i - j - 1, p + k);
      }
      
      // d_psi_d_c and d_psi_d_v:
      for(int j = 0; j < std::min(p, i - newDayC[dayIndex]); j++){
        
        d_psi_d_theta(i, p + q + 1) += d2 * alpha[j] * psi_d1[i - j - 1] * phi_d2[i - j - 1] / phi[i - j - 1] *
          (e[i - j - 1] - v); //c
        d_psi_d_theta(i, p + q + 2) -= d2 * alpha[j] * psi_d1[i - j - 1] * phi_d2[i - j - 1] / phi[i - j - 1] *
          (abs_prim(e[i - j - 1] - v) + c); //v
      }
      for(int j = 0; j < std::min(maxPQ, i - newDayC[dayIndex]); j++){
        d_psi_d_theta(i, p + q + 1) +=
          multi1[j] * d_psi_d_theta(i - j - 1, p + q + 1); //c
        d_psi_d_theta(i, p + q + 2) +=
          multi1[j] * d_psi_d_theta(i - j - 1, p + q + 2); //v
      }
      
      // d_psi_d_d1:
      d_psi_d_theta(i, p + q + 3) = - psi_d1[i] * logPsi[i];
      for(int j = 0; j < std::min(q, i - newDayC[dayIndex]); j++) d_psi_d_theta(i, p + q + 3) += 
        (beta[j] * psi_d1[i - j - 1] + alpha[j] * psi_d1[i - j - 1] * phi_d2[i - j - 1]) *
        logPsi[i - j - 1];
      for(int j = 0; j < std::min(maxPQ, i - newDayC[dayIndex]); j++) d_psi_d_theta(i, p + q + 3) += 
        multi1[j] * d_psi_d_theta(i - j - 1, p + q + 3);
      
      // d_psi_d_d2:
      for(int j = 0; j < std::min(p, i - newDayC[dayIndex]); j++) d_psi_d_theta(i, p + q + 4) += 
        alpha[j] * psi_d1[i - j - 1] * phi_d2[i - j - 1] * logPhi[i - j - 1];
      for(int j = 0; j < std::min(maxPQ, i - newDayC[dayIndex]); j++) d_psi_d_theta(i, p + q + 4) += 
        multi1[j] * d_psi_d_theta(i - j - 1, p + q + 4);
      
      // multiplies with common factor 'multi2':
      multi2 = psi[i] / (d1 * psi_d1[i]);
      for(int col = 0; col < param.size(); col++) d_psi_d_theta(i, col) *= multi2;
      
      // if 'startType' == 2, the first max(p,q) psi[i] are set to mean (so that d_psi_d_theta = 0)
      if(startType == 2 && i - newDayC[dayIndex] >= 0 && i - newDayC[dayIndex] < maxPQ){
        for(int col = 0; col < Npara; col++){
          d_psi_d_theta(i, col) = 0;
        }
      }
    }
    
    // computes d_LL_d_theta:
    LL = get_d_LL_d_psi_distCpp(x, psi, logPsi, e, d_psi_d_theta,
                                d_LL_d_theta, dist, distPara,
                                LLi, param.size(), score);
    
    if(returnIndex == 3 || returnIndex == 4) returnList["score"] = score; 
    if(returnIndex == 4){
      
      returnList["psi"] = psi;
      returnList["e"] = e;
      returnList["d_LL_d_theta"] = d_LL_d_theta; 
      returnList["d_psi_d_theta"] = d_psi_d_theta;
      returnList["LLi"] = LLi;
    }
    
  }
  
  // computes the log-likelihood:
  if(returnIndex != 3 && returnIndex != 4) 
    LL = getLL_distCpp(x, psi, logPsi, e, dist, distPara);
  
  returnList["LL"] = LL;
  if(returnIndex == 2){
    returnList["psi"] = psi;
    returnList["e"] = e;
  }
  
  return(returnList);
}
