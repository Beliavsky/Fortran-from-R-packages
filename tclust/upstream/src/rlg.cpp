//  Rcpp::compileAttributes("C:/users/valen/onedrive/myrepo/r/tclust")

#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]


Rcpp::NumericVector media(Rcpp::NumericMatrix a, Rcpp::NumericVector cluster, int grupo) {
  int nr = a.nrow();
  int nc = a.ncol();
  Rcpp::NumericVector medias(nc);
  double acumulado;
  double contador;
  for(int i = 0; i < nc; i++){
    acumulado = 0.0;
    contador = 0.0;
    for(int j = 0; j < nr; j++){
      if(cluster[j] == grupo){
        acumulado += a(j,i);
        contador += 1;
      }
    }
    medias[i] = acumulado/contador;
  }
  return medias;
}


bool compara(Rcpp::NumericVector a, Rcpp::NumericVector b) {
  for(int i = 0; i < a.length(); i++){
    if(a[i] != b[i]) return false;
  }
  return true;
}


Rcpp::NumericMatrix localpca(Rcpp::NumericMatrix matr) {
  int n = matr.nrow(), k = matr.ncol();
  arma::mat mat(matr.begin(), n, k, false);
  arma::mat pca = arma::princomp(mat);
  return Rcpp::wrap(pca);
}


Rcpp::NumericMatrix selecciona_puntos(Rcpp::NumericMatrix a, Rcpp::NumericVector cluster, int grupo){
  int acumulado = 0;
  for(int i = 0; i < cluster.length(); i++){
    if(cluster[i] == grupo){
      acumulado += 1;
    }
  }
  Rcpp::NumericMatrix puntos2(acumulado,a.ncol());
  int j = 0;
  for(int i = 0; i < a.nrow(); i++){
    if(cluster[i] == grupo){
      puntos2.row(j) = a.row(i);
      j += 1;
    }
  }
  return(puntos2);
}


// Internal function for concentration steps (initializations) in rlg
// @name rlg_c1
// @param x Rcpp::NumericMatrix, The input data.
// @param d Rcpp::NumericVector, A numeric vector of length equal to the number 
//  of clusters to be detected. 
//  Each component of vector \code{d} indicates the intrinsic dimension of the affine subspace 
//  where observations on that cluster are going to be clustered. All the elements 
//  of vector \code{d} should be smaller than the problem dimension minus 1.
// @param alpha double, The proportion of observations to be trimmed.
// @param niter1 int, The number of concentration steps to be performed for the 
// nstart initializations. 
// @export
// [[Rcpp::export]]
Rcpp::List rlg_c1(Rcpp::NumericMatrix x, Rcpp::NumericVector d, double alpha = 0.05, int niter1 = 3) {
  int k = d.length();
  int n = x.nrow();
  int p = x.ncol();
  int n_atip = floor(alpha*n);
  Rcpp::NumericVector cluster(n);
  
  Rcpp::NumericMatrix x_aux(k,p);
  Rcpp::NumericMatrix di(n,k);
  int numRand;
  Rcpp::NumericVector minimos(n);
  Rcpp::NumericVector minimos2(n);
  double lim;
  Rcpp::NumericVector temporal;
  int asignados = 0;
  double ecm = 0;
  
  for(int i = 0; i < k; i++){
    if(d[i] == 0){
      numRand = Rcpp::runif(1, 0, n).at(0);
      x_aux.row(i) = x.row(numRand);
      for(int j = 0; j < n; j++){
        
        // VT::18.03.2024 - fix the objective function value, see mail of Luis Angel from 5.3.2024
        //  di(j,i) = sqrt(Rcpp::sum(pow(x.row(j)-x_aux.row(i),2)));
        di(j,i) = Rcpp::sum(pow(x.row(j)-x_aux.row(i),2));
      }
    } else {
      
      Rcpp::NumericMatrix ci(d[i]+1,p);
      for(int j = 0; j < (d[i]+1); j++){
        numRand = Rcpp::runif(1, 0, n).at(0);
        ci.row(j) = x.row(numRand);
      }
      for(int j = 0; j < p; j++){
        x_aux(i,j) = Rcpp::mean(ci.column(j));
      }
      Rcpp::NumericMatrix U(d[i],p);
      Rcpp::NumericMatrix eigen_vec(p,p);
      eigen_vec = localpca(ci);
      for(int j = 0; j < d[i]; j++){
        U.row(j) = eigen_vec.column(j);
      }
      Rcpp::NumericMatrix productoMat(p,p);
      arma::mat arma_identidad = arma::eye(p,p);
      arma::mat arma_U = Rcpp::as<arma::mat>(U);
      arma::mat arma_productoMat = arma_identidad-arma_U.t()*arma_U;
      arma::mat arma_puntos = Rcpp::as<arma::mat>(x);
      arma::mat arma_x = Rcpp::as<arma::mat>(x_aux);
      
      for(int j = 0; j < n; j++){
        // VT::18.03.2024 - fix the objective function value, see mail of Luis Angel from 5.3.2024
        di(j,i) = arma::norm2est(arma_productoMat*(arma_puntos.row(j)-arma_x.row(i)).t());
        di(j,i) = di(j,i) * di(j,i);
      }
      
    }
  }
  
  for(int j = 0; j < n; j++){
    minimos[j] = Rcpp::min(di.row(j));
    minimos2[j] = Rcpp::min(di.row(j));
    temporal = di.row(j);
    cluster[j] = std::min_element(temporal.begin(),temporal.end()) - temporal.begin();
  }
  std::sort(minimos.begin(), minimos.end());
  if(n_atip>0){
    lim = minimos[n-n_atip]; 
  } else {
    lim = 100000.0;
  }
  for(int j = 0; j < n; j++){
    if(minimos2[j] > lim){
      cluster[j] = 1000;
      asignados += 1;
    }
  }
  for(int j = 0; j < n; j++){
    if(minimos2[j] == lim && asignados < n_atip){
      cluster[j] = 1000;
      asignados += 1;
    }
    if(cluster[j] != 1000){
      ecm += minimos2[j];
    }
  }
  ecm = ecm/(n-n_atip);
  int iter_ext = 0;
  Rcpp::NumericVector gruposAnt(n);
  
  while(iter_ext < niter1 && !compara(cluster,gruposAnt)){
    ecm = 0;
    asignados = 0;
    for(int ind = 0; ind < n; ind++){
      gruposAnt[ind] = cluster[ind];
    }
    
    for(int i = 0; i < k; i++){
      x_aux.row(i) = media(x,cluster,i);
      if(d[i] == 0){
        for(int j = 0; j < n; j++){
            // VT::18.03.2024 - fix the objective function value, see mail of Luis Angel from 5.3.2024
            //  di(j,i) = sqrt(Rcpp::sum(pow(x.row(j)-x_aux.row(i),2)));
            di(j,i) = Rcpp::sum(pow(x.row(j)-x_aux.row(i),2));
        }
      } else {
        
        Rcpp::NumericMatrix U(d[i],p);
        Rcpp::NumericMatrix eigen_vec(p,p);
        Rcpp::NumericMatrix puntos_grupo = selecciona_puntos(x,cluster,i);
        eigen_vec = localpca(puntos_grupo);
        for(int j = 0; j < d[i]; j++){
          U.row(j) = eigen_vec.column(j);
        }
        Rcpp::NumericMatrix productoMat(p,p);
        arma::mat arma_identidad = arma::eye(p,p);
        arma::mat arma_U = Rcpp::as<arma::mat>(U);
        arma::mat arma_productoMat = arma_identidad-arma_U.t()*arma_U;
        arma::mat arma_puntos = Rcpp::as<arma::mat>(x);
        arma::mat arma_x = Rcpp::as<arma::mat>(x_aux);
        
        for(int j = 0; j < n; j++){
        // VT::18.03.2024 - fix the objective function value, see mail of Luis Angel from 5.3.2024
          di(j,i) = arma::norm2est(arma_productoMat*(arma_puntos.row(j)-arma_x.row(i)).t());
          di(j,i) = di(j,i) * di(j,i);
        }
        
      }
    }
    
    
    for(int j = 0; j < n; j++){
      minimos[j] = Rcpp::min(di.row(j));
      minimos2[j] = Rcpp::min(di.row(j));
      temporal = di.row(j);
      cluster[j] = std::min_element(temporal.begin(),temporal.end()) - temporal.begin();
    }
    std::sort(minimos.begin(), minimos.end());
    if(n_atip>0){
      lim = minimos[n-n_atip]; 
    } else {
      lim = 100000.0;
    }
    for(int j = 0; j < n; j++){
      if(minimos2[j] > lim){
        cluster[j] = 1000;
        asignados += 1;
      }
    }
    for(int j = 0; j < n; j++){
      if(minimos2[j] == lim && asignados < n_atip){
        cluster[j] = 1000;
        asignados += 1;
      }
      if(cluster[j] != 1000){
        ecm += minimos2[j];
      }
    }
    //ecm = ecm/(n-n_atip);
    iter_ext += 1;
  }
  
  cluster = cluster + 1;
  cluster[cluster == 1001] = 0;
  
  return Rcpp::List::create(
    Rcpp::_["cluster"] = cluster,
    Rcpp::_["obj"] = ecm
  );
}

// Internal function for concentration steps (refinement) in rlg
// @name rlg_c2
// @param x Rcpp::NumericMatrix, The input data.
// @param d Rcpp::NumericVector, A numeric vector of length equal to the number 
//  of clusters to be detected. 
//  Each component of vector \code{d} indicates the intrinsic dimension of the affine subspace 
//  where observations on that cluster are going to be clustered. All the elements 
//  of vector \code{d} should be smaller than the problem dimension minus 1.
// @param cluster Rcpp::NumericVector  A numeric vector of size \code{n} containing 
//  the cluster assignment for each observation. Cluster names are integer numbers 
//  from 1 to \code{k}, 0 indicates trimmed observations.
// @param alpha double, The proportion of observations to be trimmed.
// @param niter2 int, The maximum number of concentration steps to be performed for 
//  the \code{nkeep} solutions kept for further iteration. The concentration steps 
//  are stopped, whenever two consecutive steps lead to the same data partition. 
// @export
// [[Rcpp::export]]
Rcpp::List rlg_c2(Rcpp::NumericMatrix x, Rcpp::NumericVector d, Rcpp::NumericVector cluster, double alpha = 0.05, int niter2 = 20) {
  int k = d.length();
  int n = x.nrow();
  int p = x.ncol();
  int n_atip = floor(alpha*n);
  
  Rcpp::NumericMatrix x_aux(k,p);
  Rcpp::NumericMatrix di(n,k);
  Rcpp::NumericVector minimos(n);
  Rcpp::NumericVector minimos2(n);
  double lim;
  Rcpp::NumericVector temporal;
  Rcpp::NumericVector gruposAnt(n);
  int asignados = 0;
  double ecm = 0;
  int flag = 0;
  
  cluster = cluster - 1; // Clusters should be 0..k-1 and 1000 for trimmed cluster for compatibility with previous versions
  cluster[cluster == -1] = 1000;
  
  while(flag < niter2 && !compara(cluster,gruposAnt)){
    
    
    
    ecm = 0;
    asignados = 0;
    for(int ind = 0; ind < n; ind++){
      gruposAnt[ind] = cluster[ind];
    }
    
    for(int i = 0; i < k; i++){
      x_aux.row(i) = media(x,cluster,i);
      if(d[i] == 0){
        for(int j = 0; j < n; j++){
            // VT::18.03.2024 - fix the objective function value, see mail of Luis Angel from 5.3.2024
            //  di(j,i) = sqrt(Rcpp::sum(pow(x.row(j)-x_aux.row(i),2)));
            di(j,i) = Rcpp::sum(pow(x.row(j)-x_aux.row(i),2));
        }
      } else {
        Rcpp::NumericMatrix U(d[i],p);
        Rcpp::NumericMatrix eigen_vec(p,p);
        Rcpp::NumericMatrix puntos_grupo = selecciona_puntos(x,cluster,i);
        eigen_vec = localpca(puntos_grupo);
        for(int j = 0; j < d[i]; j++){
          U.row(j) = eigen_vec.column(j);
        }
        Rcpp::NumericMatrix productoMat(p,p);
        arma::mat arma_identidad = arma::eye(p,p);
        arma::mat arma_U = Rcpp::as<arma::mat>(U);
        arma::mat arma_productoMat = arma_identidad-arma_U.t()*arma_U;
        arma::mat arma_puntos = Rcpp::as<arma::mat>(x);
        arma::mat arma_x = Rcpp::as<arma::mat>(x_aux);
        
        for(int j = 0; j < n; j++){
          // VT::18.03.2024 - fix the objective function value, see mail of Luis Angel from 5.3.2024
          di(j,i) = arma::norm2est(arma_productoMat*(arma_puntos.row(j)-arma_x.row(i)).t());
          di(j,i) = di(j,i) * di(j,i);
        }
      }
    }
    
    for(int j = 0; j < n; j++){
      minimos[j] = Rcpp::min(di.row(j));
      minimos2[j] = Rcpp::min(di.row(j));
      temporal = di.row(j);
      cluster[j] = std::min_element(temporal.begin(),temporal.end()) - temporal.begin();
    }
    std::sort(minimos.begin(), minimos.end());
    if(n_atip>0){
      lim = minimos[n-n_atip]; 
    } else {
      lim = 100000.0;
    }
    for(int j = 0; j < n; j++){
      if(minimos2[j] > lim){
        cluster[j] = 1000;
        asignados += 1;
      }
    }
    for(int j = 0; j < n; j++){
      if(minimos2[j] == lim && asignados < n_atip){
        cluster[j] = 1000;
        asignados += 1;
      }
      if(cluster[j] != 1000){
        ecm += minimos2[j];
      }
    }
    //ecm = ecm/(n-n_atip);
    flag += 1;
    
  }
  
  cluster = cluster + 1;
  cluster[cluster == 1001] = 0;
  
  return Rcpp::List::create(
    Rcpp::_["cluster"] = cluster,
    Rcpp::_["obj"] = ecm
    );
}
