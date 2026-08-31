#include <Rcpp.h>
using namespace Rcpp;
using namespace std;

double euclidCpp(double x, double y, double p) {
  return pow(x-y,2);
}

double minkowskiCpp(double x, double y, double p) {
  return pow(abs(x-y),p);
}

List closestEdistXYCpp(int ncol, NumericMatrix x, int xnrow, NumericMatrix y, int ynrow,
                      int part, double p, int fct_var, double eta,
                      std::vector<int> colindices, IntegerVector rowpointers, std::vector<double> entries){

  int jja, i,j,k, jfrom, jto;
  double etap, tmp, pinv;

  etap = pow(eta,p);
  pinv = 1/p;

  jja = 1;

  rowpointers[0] = 1;
  jfrom = 0;
  jto = ynrow-1;

  for (i=0; i<xnrow; i++) {

    if (part < 0) {
      jto = i;
    }
    if (part > 0) {
      jfrom = i;
    }

    for (j = jfrom; j<=jto; j++) {

      tmp = 0.0;
      int flag = 0;
      for (k = 0; k<ncol; k++) {
        if (fct_var == 1) {
          tmp = tmp + euclidCpp(x(i,k),y(j,k),p);
        } else{
          tmp = tmp + minkowskiCpp(x(i,k),y(j,k),p);
        }
        if( tmp>etap){
          flag = 1;
          break;
        }
      }

      if (flag == 1){
        continue;
      }

      colindices.push_back(j+1);
      if (abs(p-2) <= 0.0) {
        entries.push_back(sqrt(tmp));
      } else {
        if (abs(p-1) <= 0.0) {
          entries.push_back(tmp);
        } else {
          entries.push_back(pow(tmp,pinv));
        }
      }
      jja = jja + 1;
    }
    rowpointers[i+1]=jja;
  }

  if (part>0) {
    rowpointers[xnrow]=jja;
  }

  return List::create(Named("entries") = entries,
                      Named("colindices") = colindices,
                      Named("rowpointers")= rowpointers);
}

List closestMAXdistXYCpp(int ncol, NumericMatrix x, int xnrow, NumericMatrix y, int ynrow,
                        int part,
                        double eta,
                        std::vector<int> colindices, IntegerVector rowpointers, std::vector<double> entries) {


  int jja, i,j,k, jfrom, jto;
  double tmp;

  jja=1;

  rowpointers[0]=1;
  jfrom = 0;
  jto = ynrow-1;

  for (i=0; i<xnrow; i++) {

    if (part < 0) {
      jto = i;
    }
    if (part > 0) {
      jfrom = i;
    }

    for (j = jfrom; j<=jto; j++) {

      tmp = 0.0;
      int flag = 0;
      for (k = 0; k<ncol; k++) {
        tmp = std::max(tmp, abs(x(i,k)-y(j,k)));
        if( tmp>eta){
          flag = 1;
          break;
        }
      }

      if (flag == 1){
        continue;
      }

      colindices.push_back(j+1);
      entries.push_back(tmp);
      jja = jja + 1;
    }
    rowpointers[i+1]=jja;
  }

  if (part>0) {
    rowpointers[xnrow]=jja;
  }

  return List::create(Named("entries") = entries,
                      Named("colindices") = colindices,
                      Named("rowpointers")= rowpointers);
}


List closestGCdistXYCpp(NumericMatrix x, int nx, NumericMatrix y, int ny,
                       int part, double p,
                       double eta,
                       std::vector<int> colindices, IntegerVector rowpointers, std::vector<double> entries) {


  bool equi;
  int jja, i,j, jfrom, jto;
  double etap, tmp, tmp1, tmp2;
  double rad, thres;
  NumericVector scy12(ny), ccy12(ny), sy2(ny);
  double scx12, ccx12, sx2;

  rad = 0.01745329251994329;
  thres = 0.99999999999;

  if (p < 0) {
    equi=true;
    p=-p;
  } else {
    equi= false;
  }

  jja=1;

  etap=cos(eta*rad);
  rowpointers[0]=1;
  jfrom = 0;
  jto = ny-1;

  for (j=0; j<ny; j++) {
    tmp1=y(j,0)*rad;
    tmp2=y(j,1)*rad;
    ccy12[j]=cos(tmp1)*cos(tmp2);
    scy12[j]=sin(tmp1)*cos(tmp2);
    sy2[j]=sin(tmp2);
  }

  for (i=0; i<nx; i++) {

    if (equi) {
      ccx12=ccy12[i];
      scx12=scy12[i];
      sx2=sy2[i];
    } else {
      tmp1=x(i,0)*rad;
      tmp2=x(i,1)*rad;
      ccx12=cos(tmp1)*cos(tmp2);
      scx12=sin(tmp1)*cos(tmp2);
      sx2=sin(tmp2);
    }

    if (part < 0) {
      jto = i;
    }
    if (part > 0) {
      jfrom = i;
    }

    for (j = jfrom; j<=jto; j++) {

      tmp = ccx12 * ccy12[j] + scx12 * scy12[j]  + sx2*sy2[j];

      if (tmp < etap) continue;

      if  (tmp >= thres) {
        tmp = 0.0;
      } else {
        tmp = acos( tmp);
      }
      colindices.push_back(j+1);
      entries.push_back(tmp*p);
      jja = jja + 1;
    }
    rowpointers[i+1]=jja;
  }

  if (part>0) {
    rowpointers[nx]=jja;
  }

  return List::create(Named("entries") = entries,
                      Named("colindices") = colindices,
                      Named("rowpointers")= rowpointers);
}


// [[Rcpp::export]]
List closestdistCpp(NumericMatrix x, NumericMatrix y,
                   int part, double p, int method, double eta) {

  int ncol = x.ncol();
  int nrowx = x.nrow();
  int nrowy = y.nrow();

  std::vector<int> colindices;
  std::vector<double> entries;
  IntegerVector rowpointers(nrowx+1);

  if (method==1) {
    p=2.0;
    int fct_var = 1;
    return closestEdistXYCpp(ncol, x, nrowx, y, nrowy,
                            part, p, fct_var, eta,
                            colindices, rowpointers, entries);
  }
  if (method==2) {
    return closestMAXdistXYCpp(ncol, x, nrowx, y, nrowy,
                              part,
                              eta, colindices, rowpointers, entries);
  }
  if (method==3) {
    int fct_var = 2;
    return closestEdistXYCpp(ncol, x, nrowx, y, nrowy,
                            part, p, fct_var, eta,
                            colindices, rowpointers, entries);
  }
  if (method==4) {
    return closestGCdistXYCpp(x, nrowx, y, nrowy,
                             part, p,
                             eta, colindices, rowpointers, entries);
  }
  return List::create(Named("entries") = entries,
                      Named("colindices") = colindices,
                      Named("rowpointers")= rowpointers);
}
