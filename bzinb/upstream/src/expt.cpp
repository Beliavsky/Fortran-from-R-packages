// [[Rcpp::depends(BH)]]
#include <Rcpp.h>
#include <iostream>
#include <math.h>
#include <numeric>
#include <vector>
#include <string>
#include "bzinb.h"
#include <boost/math/special_functions/digamma.hpp>
#include <boost/math/special_functions/trigamma.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
using namespace std;
using namespace Rcpp;
using namespace boost::numeric::ublas;

// 1. Expt
void l1(int& x, int& y, double& a0, double& a1, double& a2, int &k, int& m, 
        double& result, double adjj)
{
  result = lgamma(a1 + k)- lgamma(k+1)- lgamma(a1) + lgamma(x + y + a0 -m -k) 
  - lgamma(x -k +1) - lgamma(a0 + y - m) + lgamma(m + a2) - lgamma(m+1) 
  - lgamma(a2) + lgamma(y +a0 -m) - lgamma(y -m +1) - lgamma(a0) - adjj;
  result = exp(result);
}
void l1_c (double& t1, double& t2, int& k, int& m, double& result, double adjj)
{
  result = exp(k * log(t1) + m * log(t2) - adjj);
}
void l1_AC (double& t1, double& t2, int& x, int& y, 
            double& a0, double& a1, double& a2, 
            int& k, int& m, double& result,  double adjj)
{
  result = lgamma(a1 + k)- lgamma(k+1)- lgamma(a1) + lgamma(x + y + a0 -m -k) 
  - lgamma(x -k +1) - lgamma(a0 + y - m) + lgamma(m + a2) - lgamma(m+1) 
  - lgamma(a2) + lgamma(y +a0 -m) - lgamma(y -m +1) - lgamma(a0) 
  + k *log(t1) + m *log(t2) - adjj;
  result = exp(result);
}
void l2_A (int& x, double& a0, double& a1, double& a2, int& k, 
           double& result, double adjj)
{
  result = lgamma(x +a0 -k) + lgamma(k + a1) - lgamma(a0) - lgamma(x-k+1) - lgamma(a1) - lgamma(k+1) - adjj;
  result = exp(result);
}
void l3_A (int& y, double& a0, double& a1, double& a2, int& m, double& result, double adjj)
{
  result = lgamma(y +a0 -m) + lgamma(m + a2) - lgamma(a0) - lgamma(y-m+1) - lgamma(a2) - lgamma(m+1) - adjj;
  result = exp(result);
}
void R0_E1(int& x, int& y, int& k, int& m, double& a0, double& result)
{
  result = x - k + y - m + a0;
}
double log_R0_E1(int& x, int& y, int& k, int& m, double& a0)
{
  return(boost::math::digamma(x - k + y - m + a0));
}
double	log_R0_E2(int& x, double& a0, int& k)
{
  return(boost::math::digamma(x - k + a0));
}
double	log_R0_E3(int& y, double& a0, int& m)
{
  return(boost::math::digamma(y - m + a0));
}
void R1_E1(int& k, double& a1, double& result)
{
  result = k + a1;
}
double log_R1_E1(int& k, double& a1)
{
  return(boost::math::digamma(k + a1));
}
// double log_R1_E2(int& k, double& a1)  //not used
// {
//   return(boost::math::digamma(k + a1));
// }
void R2_E1(int& m, double& a2, double& result)
{
  result = m + a2;
}
double log_R2_E1(int& m, double& a2)
{
  return(boost::math::digamma(m + a2));
}
// double log_R2_E3(int m, double a2)  //not used
// {
//   return(boost::math::digamma(m + a2));
// }

// [[Rcpp::export]]
void dBvZINB_Expt(int &x, int &y, int &freq, double &a0, double &a1, double &a2,
                  double &b1, double &b2, double &p1, double &p2, 
                  double &p3, double &p4,
                  NumericVector &expt, NumericVector &s_i, NumericVector &info, int se, int bnb)
{
  double t1 = (float)(b1 + b2 + 1) /(b1 + 1);
  double t2 = (float)(b1 + b2 + 1) /(b2 + 1);
  double zeta_xy = 0.0;
  if (x == 0 && y == 0) {zeta_xy = 1.0;}
  
  double adj_A = 0.0;
  double adj_B1 = 0.0;
  double adj_C = 0.0;
  double adj_sum = 0.0;
  double l1_B = - (x + y + a0) * log(1.0 + b1 + b2) + x * log(b1) + y * log(b2) - a1 * log(1.0 + b1) - a2 * log(1.0 + b2);
  double l2_B = (y!=0?0.0:(exp(- (x + a0 + a1) * log(1.0 + b1) + x * log(b1) + adj_B1) * p2));
  double l3_B = (x!=0?0.0:(exp(- (y + a0 + a2) * log(1.0 + b2) + y * log(b2) + adj_B1) * p3));
  double l4_B = ((x+y!=0)?0.0:( p4 * exp(adj_B1)));
  
  double sum_A = 0, sum_C = 0, sum_AC = 0;
  matrix <double> R0_mat(x+1, y+1); 
  matrix <double> R1_mat(x+1, y+1); 
  matrix <double> R2_mat(x+1, y+1); 
  matrix <double> log_R0_mat(x+1, y+1); 
  matrix <double> log_R1_mat(x+1, y+1); 
  matrix <double> log_R2_mat(x+1, y+1);
  matrix <double> l_A_mat(x+1, y+1); 
  matrix <double> l_C_mat(x+1, y+1); 
  matrix <double> l_AC_mat(x+1, y+1); 
  std::vector <double> log_R0_mat2(x+1); 
  std::vector <double> log_R0_mat3(y+1); 
  std::vector <double> l2_A_mat(x+1); 
  std::vector <double> l3_A_mat(y+1);
  std::vector <double> log_R1_mat2(x+1); 
  std::vector <double> log_R2_mat2(x+1); 
  std::vector <double> log_R1_mat3(y+1); 
  std::vector <double> log_R2_mat3(y+1);
  
  // cout << endl << endl << "line 132; x = "<< x << ", y = " << y << endl;
  for(int i = 0;i <= x;i++)
  {
    // cout << "line 139; i,j = " << i <<" "<< 1 <<" l_A_mat(i,1)= " << l_A_mat(i, 1)<<" l_C_mat(i,1)= " << l_C_mat(i, 1) << endl;
    for(int j = 0;j <= y;j++)
    {
      l1(x, y, a0, a1, a2, i, j, l_A_mat(i, j), adj_A);  // updating l_A_mat(i, j)
      l1_c(t1, t2, i, j, l_C_mat(i, j), adj_C);
      sum_A += l_A_mat(i, j);
      sum_C += l_C_mat(i, j);
      //system("pause");
    }
  }
  // cout << "line 139; l1B = " << l1_B << endl;
  
  if (l1_B < -200 && log(l2_B + l3_B + l4_B) < 0)
  {
    // adj_B1 = ((-l1_B - 200)*1.0 / 100) * 100; // prevent exp(l1_B) from being 0
    adj_B1 = -l1_B - 200.0; // prevent exp(l1_B) from being 0
    l1_B += adj_B1;
  }
  l1_B = exp(l1_B) * p1;

  while (log(sum_A) > 250)
  {
    sum_A = 0;
    adj_A += 200;
    for(int i = 0;i <= x;i++)
    {
      for(int j = 0;j <= y;j++)
      {
        l1(x, y, a0, a1, a2, i, j, l_A_mat(i, j), adj_A);
        sum_A += l_A_mat(i, j);
      }
    }
  }
  for(int i = 0;i <= x;i++)
  {
    l2_A(x, a0, a1, a2, i, l2_A_mat[i], adj_A);  //update l2_A_mat[i]
  }
  for(int j = 0;j <= y;j++)
  {
    l3_A(y, a0, a1, a2, j, l3_A_mat[j], adj_A) ; //update l3_A_mat[i]
  }
  while (log(sum_C) > 250)
  {
    sum_C = 0;
    adj_C += 200;
    for(int i = 0;i <= x;i++)
    {
      for(int j = 0;j <= y;j++)
      {
        l1_c(t1, t2, i, j, l_C_mat(i, j), adj_C); //update l_C_mat(i, j) 
        sum_C += l_C_mat(i, j);
      }
    }
  }
  sum_A = 0;
  sum_AC = 0;
  for(int i = 0;i <= x;i++)
  {
    for(int j = 0;j <= y;j++)
    {
      sum_AC += l_A_mat(i, j)*l_C_mat(i, j); //code simplified;
      sum_A += l_A_mat(i, j);                //code missed;
    }
  }

  if(log(sum_AC) > 200)
  {
    adj_A += 100;
    adj_C += 100;
    sum_AC = 0;
    sum_A = 0;
    for(int i = 0;i <= x;i++)
    {
      for(int j = 0;j <= y;j++)
      {
        l1(x, y, a0, a1, a2, i, j, l_A_mat(i, j), adj_A); //update l_A_mat(i, j)
        l1_c(t1, t2, i, j, l_C_mat(i, j), adj_C); //update l_C_mat(i, j) 
        sum_AC += l_A_mat(i, j)*l_C_mat(i, j);
        sum_A +=  l_A_mat(i, j);
      }
    }
  }
  else if(log(sum_AC) < - 100)
  {
    adj_A = adj_A - 200;
    adj_C = adj_C - 200;
    sum_AC =0;
    sum_A = 0;
    for(int i = 0;i <= x;i++)
    {
      for(int j = 0;j <= y;j++)
      {
        l1(x, y, a0, a1, a2, i, j, l_A_mat(i, j), adj_A);  // update l_A_mat(i, j)
        l1_c(t1, t2, i, j, l_C_mat(i, j), adj_C); //update l_C_mat(i, j) 
        l1_AC(t1, t2, x, y, a0, a1, a2, i, j, l_AC_mat(i, j), adj_A+adj_C);  // update l1_AC
        sum_AC += l_AC_mat(i, j);
        sum_A +=  l_A_mat(i, j);
      }
    }
  }
  
  // updating l2_A and l3_A with updated adj_A.
  for(int i = 0;i <= x;i++)
  {
    l2_A(x, a0, a1, a2, i, l2_A_mat[i], adj_A);  //update l2_A_mat[i]
  }
  for(int j = 0;j <= y;j++)
  {
    l3_A(y, a0, a1, a2, j, l3_A_mat[j], adj_A) ; //update l3_A_mat[i]
  }
  
  
  double l_sum =0;
  l_sum = sum_AC * l1_B + sum_A * (l2_B +  l3_B +  l4_B) * exp(-adj_C);
  // cout << l_sum << " " << sum_AC << " " << l1_B << " " << sum_A << " " << l2_B << " " << l3_B << " " << l4_B << " " << adj_C << endl;
  if (l_sum == 0)
  {
    adj_sum = -floor(log(sum_AC)*2*1.0/3 + log(l1_B)*2*1.0/3);
    // cout << adj_sum<<"adj";
    // cout << sum_AC <<"AC"<< adj_sum <<"@"<< l1_B<<"c"<<adj_C<<endl;
    l_sum = sum_AC * exp(adj_sum) * l1_B + sum_A * (exp(adj_sum) * (l2_B +  l3_B +  l4_B)) * exp(-adj_C);
  }
  double R0_E1_B, R0_E2_B, R0_E3_B, R0_E4_B, R1_E1_B, R1_E2_B, R1_E3_B, R1_E4_B;
  double R2_E1_B, R2_E2_B, R2_E3_B, R2_E4_B;
  // expectation components
  R0_E1_B = (double)b1/(1 + b1 + b2);
  R0_E2_B = (double)b1/(1 + b1);
  R0_E3_B = (double)b1/(1 + b2);
  R0_E4_B = (double)b1;
  
  R1_E1_B = (double)b1/(1 + b1);
  R1_E2_B = (double)b1/(1 + b1);
  R1_E3_B = (double)b1;
  R1_E4_B = (double)b1;
  
  R2_E1_B = (double)b1/(1 + b2);
  R2_E2_B = (double)b1;
  R2_E3_B = (double)b1/(1 + b2);
  R2_E4_B = (double)b1;
  
  double log_R0_E = 0, log_R1_E = 0, log_R2_E = 0, R0_E = 0, R1_E = 0, R2_E = 0;
  for(int i = 0;i <= x;i++)
  {
    log_R0_mat2[i] = l2_A_mat[i] * (log_R0_E2(x, a0, i) + log(R0_E2_B));
    log_R1_mat2[i] = l2_A_mat[i] * (log_R1_E1(i,a1) + log (R1_E2_B));
    //cout << "line 238; log_R1_mat2[i] = " << log_R1_mat2[i] << ", l2_A_mat[i] = " << l2_A_mat[i]  << ", logR1E2(x,a0,i) = " << log_R1_E2(i, a1) << ", log(R1E2B)= "<< log(R1_E2_B) << endl;
    log_R2_mat2[i] = l2_A_mat[i] * (boost::math::digamma(a2) + log(R2_E2_B));
    for(int j = 0;j <= y;j++)
    {
      R0_E1(x, y, i, j, a0, R0_mat(i, j));  //update R0_mat(i, j)
      R0_mat(i, j) = R0_mat(i, j)*l_A_mat(i, j);
      R0_E += R0_mat(i, j) * l_C_mat(i, j) * exp(adj_sum) * l1_B * R0_E1_B + R0_mat(i, j) * (l2_B * R0_E2_B + l3_B * R0_E3_B + l4_B * R0_E4_B)*exp(-adj_C + adj_sum);
      R1_E1(i, a1, R1_mat(i, j));  //update R1_mat(i, j)
      R2_E1(j, a2, R2_mat(i, j));  //update R2_mat(i, j)
      R1_mat(i, j) *= l_A_mat(i, j);
      R2_mat(i, j) *= l_A_mat(i, j);
      log_R0_mat(i, j) = (log_R0_E1(x, y, i, j, a0) + log(R0_E1_B)) * l_A_mat(i, j);
      log_R1_mat(i, j) = l_A_mat(i, j) * (log_R1_E1(i,a1) + log (R1_E1_B));
      log_R2_mat(i, j) = (log_R2_E1(j,a2) + log(R2_E1_B)) * l_A_mat(i, j);
      //cout << "line252; log_R1_mat(i, j): " << log_R1_mat(i, j) << ", (log_R1_E1(i,a1) + log (R1_E1_B)): " << (log_R1_E1(i,a1) + log (R1_E1_B)) << ", log_R1_E1(i,a1):" << log_R1_E1(i,a1) << ", log (R1_E1_B): " << log (R1_E1_B) <<endl;
      log_R0_E += log_R0_mat(i, j) * l_C_mat(i, j) * exp(adj_sum) * l1_B;
      //cout << "i = " << i << ", j = " << j << ", logR0E = " << log_R0_E << " " << endl;
      log_R1_E += log_R1_mat(i, j) * l_C_mat(i, j) * exp(adj_sum) * l1_B;
      //cout << "line257; sum = " << log_R1_mat[i][0] * l_C_mat[i][0] * exp(adj_sum - adj_C) * l1_B << ", lCmat: " <<l_C_mat[i][0] << ", l1B: " << l1_B <<endl;
      log_R2_E += log_R2_mat(i, j) * l_C_mat(i, j) * exp(adj_sum) * l1_B;
      // cout << "line286; log_R1_E = " << log_R1_E << ", ";
      // cout << "log_R2_E = " << log_R2_E << endl;
      
      R1_E += R1_mat(i, j) * l_C_mat(i, j) * exp(adj_sum) * l1_B * R1_E1_B + R1_mat(i, j) * (l2_B * R1_E2_B + l3_B * R1_E3_B + l4_B * R1_E4_B)*exp(-adj_C + adj_sum);
      R2_E += R2_mat(i, j) * l_C_mat(i, j) * exp(adj_sum) * l1_B * R2_E1_B + R2_mat(i, j) * (l2_B * R2_E2_B + l3_B * R2_E3_B + l4_B * R2_E4_B)*exp(-adj_C + adj_sum);
    }
    log_R0_E += (log_R0_mat2[i] * l2_B) * exp(adj_sum - adj_C);
    //cout << "line262; log_R1_E(cum1) = " << log_R1_E << ", line sum = " << (log_R1_mat2[i] * l2_B) * exp(adj_sum) << "log_R1_E(cum2) = " << log_R1_E + (log_R1_mat2[i] * l2_B) * exp(adj_sum) << endl;
    log_R1_E += (log_R1_mat2[i] * l2_B) * exp(adj_sum - adj_C);
    log_R2_E += (log_R2_mat2[i] * l2_B) * exp(adj_sum - adj_C);
    // cout << "line296; log_R1_E = " << log_R1_E << ", ";
    // cout << "log_R2_E = " << log_R2_E << endl;
    
    //cout << "line263; cum = " << log_R1_E << " sum = " << (log_R1_mat2[i] * l2_B) * exp(adj_sum) << ", logR1E = " << log_R1_E << ", l2_B = " << l2_B << ", exp(adj_sum) = " << exp(adj_sum) << endl;
    //cout << "logR1mat2.i = " << log_R1_mat2[i] << endl;
  }
  //cout << log_R0_E << " " << endl;
  for(int j = 0; j <= y;j++)
  {
    log_R0_mat3[j] = l3_A_mat[j] * (log_R0_E3(y, a0, j) + log(R0_E3_B));
    log_R1_mat3[j] = l3_A_mat[j] * (boost::math::digamma(a1) + log(R1_E3_B));
    log_R2_mat3[j] = l3_A_mat[j] * (log_R2_E1(j,a2) + log (R2_E3_B));
    log_R0_E += (log_R0_mat3[j] * l3_B) * exp(adj_sum - adj_C);
    log_R1_E += (log_R1_mat3[j] * l3_B) * exp(adj_sum - adj_C);
    log_R2_E += (log_R2_mat3[j] * l3_B) * exp(adj_sum - adj_C);
  // cout << "line311; log_R1_mat3= " << log_R1_mat3[j] << "cum log_R1_E= " << log_R1_E << "/ log_R2_mat3= " << log_R2_mat3[j] << ", cum log_R2_E = " << log_R2_E << endl;
  }
  
  
  R0_E /= l_sum;
  R1_E /= l_sum;
  R2_E /= l_sum;

  log_R0_E += (boost::math::digamma(a0) + log(b1)) * exp(adj_sum - adj_C) * l4_B;
  log_R0_E /= l_sum;
  log_R1_E += (boost::math::digamma(a1) + log(b1)) * exp(adj_sum - adj_C) * l4_B;
  log_R1_E /= l_sum;
  log_R2_E += (boost::math::digamma(a2) + log(b1)) * exp(adj_sum - adj_C) * l4_B;
  log_R2_E /= l_sum;

  double E_E1, E_E2, E_E3, E_E4;
  if (bnb == 0) {
    // cout << sum_AC << l1_B << sum_A << l1_B << l2_B << l3_B << l4_B << adj_C << adj_sum << endl;
    E_E1 = sum_AC * exp(adj_sum) * l1_B;
    E_E2 = sum_A * l2_B * exp(-adj_C + adj_sum);
    E_E3 = sum_A * l3_B * exp(-adj_C + adj_sum);
    E_E4 = sum_A * l4_B * exp(-adj_C + adj_sum);
    
    double su = E_E1 + E_E2 + E_E3 + E_E4;
    E_E1 /= su;
    E_E2 /= su;
    E_E3 /= su;
    E_E4 /= su;
  } else {
    E_E1 = 1;
    E_E2 = 0;
    E_E3 = 0;
    E_E4 = 0;
  }
  
  
  double v_E = (y==0 ? (a0 + a2) * b2 *(E_E2 + E_E4) : y);
  //double v_E = (y==0 ? (a0 + a2) * b2 /(1 + (p1 + p3 * zeta_xy * exp(-(a0+a2) * log(1 + b2)))/(p2 + p4 * zeta_xy)) : y);
  
  /*double xx = sum_AC * exp(adj_sum) * l1_B * y ;
  cout<<xx<<"&"<<endl;
  double yy = sum_A * l2_B * a2 * b2*exp(-adj_C + adj_sum) +
  dnbinom(x, a0 + a1 + 1, b1*1.0/(1+b1)) * exp(-adj_A - adj_C + adj_sum) * a0 * b2 * p2 * y==0?1:0 +
  sum_A * l3_B * y *exp(-adj_C + adj_sum) +
  sum_A * l4_B * (a0 + a2) * b2 *exp(-adj_C + adj_sum);
  cout << yy <<"%"<<endl;
  double v_E = xx + yy;
  v_E= v_E*1.0/l_sum;*/
  //cout << "adj.a " << adj_A << " adj.b " << adj_B1 << " adj.c " << adj_C << " adj.sum " << adj_sum << endl;  
  expt[0] += (log(l_sum) + adj_A -adj_B1 + adj_C - adj_sum) * freq;
  expt[1] += R0_E * freq;
  expt[2] += R1_E * freq;
  expt[3] += R2_E * freq;
  expt[4] += log_R0_E * freq;
  expt[5] += log_R1_E * freq;
  expt[6] += log_R2_E * freq;
  expt[7] += E_E1 * freq;
  expt[8] += E_E2 * freq;
  expt[9] += E_E3 * freq;
  expt[10] += E_E4 * freq;
  expt[11] += v_E * freq; //change lines between E_E and v_E;
  
  if (se == 1) {
    double Dvec[9];
    
    for (int i=0; i<9; i++) {
      Dvec[i] = 0;
    }   
    //Dvec[0] ~ Dvec[4];
    // part I (functional components only)
    for (int i = 0; i <= x; i++) {
      for (int j = 0; j <= y; j++) {
        Dvec[0] += l_A_mat(i, j) * exp(adj_sum) * l_C_mat(i, j) * boost::math::digamma(x - i + y - j + a0);
//cout << "elements["  << i << "]["<< j << "]: "  << l_A_mat(i, j) * exp(adj_sum) * l_C_mat(i, j) * boost::math::digamma(x - i + y - j + a0) << endl;
        Dvec[1] += l_A_mat(i, j) * exp(adj_sum) * l_C_mat(i, j) * boost::math::digamma(a1 + i);
        Dvec[2] += l_A_mat(i, j) * exp(adj_sum) * l_C_mat(i, j) * boost::math::digamma(a2 + j);
        Dvec[3] += l_A_mat(i, j) * exp(adj_sum) * l_C_mat(i, j) * (x/b1 - (x + y - i - j + a0) / (b1 + b2 + 1) - (a1 + i)/(b1 + 1));
        Dvec[4] += l_A_mat(i, j) * exp(adj_sum) * l_C_mat(i, j) * (y/b2 - (x + y - i - j + a0) / (b1 + b2 + 1) - (a2 + j)/(b2 + 1));
      }
    }
//cout << "Dvec (line 343): " << Dvec[0] << endl; 
    // part II (adding & multiplying constants)
    //l1_b has pi_1 multiplied, but for Dvec[0] ~ Dvec[4] pi_1 should not be removed, but later this will be cancelled by not doing lines around 384.
    Dvec[0] = (Dvec[0] * l1_B - sum_AC * exp(adj_sum) * l1_B * (boost::math::digamma(a0) + log(b1 + b2 + 1)));
    Dvec[1] = (Dvec[1] * l1_B - sum_AC * exp(adj_sum) * l1_B * (boost::math::digamma(a1) + log(b1 + 1)));
    Dvec[2] = (Dvec[2] * l1_B - sum_AC * exp(adj_sum) * l1_B * (boost::math::digamma(a2) + log(b2 + 1)));
    Dvec[3] *= l1_B;
    Dvec[4] *= l1_B;
//cout << "Dvec (line 357): " << Dvec[0] << endl;
    // part III (back-adjustments)
    for (int i = 0; i < 5; i++) {
      Dvec[i] *= exp(adj_A + adj_C - adj_B1 - adj_sum);
    }
//cout << "Dvec (line 362): " << Dvec[0] << endl;
    
    //Dvec[5] ~ Dvec[8]; They are only needed when the counterpart (y and x) are zero: so zeta(y) or zeta(x) are already applied.
    double fx = 0, fy = 0;
    if (y == 0) {
      fx = lgamma(x + a0 + a1) - lgamma(x + 1) - lgamma(a0 + a1) + x * log(b1) - (x + a0 + a1) * log(b1 + 1);
      fx = exp(fx);
      Dvec[5] = fx * (boost::math::digamma(x + a0 + a1) - boost::math::digamma(a0 + a1) - log(b1 + 1));
      Dvec[6] = fx * (x/b1 - (x + a0 + a1)/(b1 + 1));
    }
    if (x == 0) {
      fy = lgamma(y + a0 + a2) - lgamma(y + 1) - lgamma(a0 + a2) + y * log(b2) - (y + a0 + a2) * log(b2 + 1);
      fy = exp(fy);
      Dvec[7] = fy * (boost::math::digamma(y + a0 + a2) - boost::math::digamma(a0 + a2) - log(b2 + 1));
      Dvec[8] = fy * (y/b2 - (y + a0 + a2)/(b2 + 1));
    }
    
    // part IV Dvec[?] x pi
    
    // Dvec[0] ~ Dvec[4] has been already multiplied by pi_1. So this should not be done again.
    // for (int i = 0; i < 5; i++) {
    //   Dvec[i] *= p1;
    // }
    for (int i = 5; i < 7; i++) {
      Dvec[i] *= p2;
    }
    for (int i = 7; i < 9; i++) {
      Dvec[i] *= p3;
    }
    
// cout << "Dvec: " ;  
// for (int i = 0; i <9; i++) {
//   cout << "[" << i << "]: " << Dvec[i] << " ";
// }
// cout << endl;
    
    s_i[0] = Dvec[0] + Dvec[5] + Dvec[7];
    s_i[1] = Dvec[1] + Dvec[5];
    s_i[2] = Dvec[2] + Dvec[7];
    s_i[3] = Dvec[3] + Dvec[6];
    s_i[4] = Dvec[4] + Dvec[8];
    s_i[5] = sum_AC * exp(adj_sum) * l1_B * exp(adj_A + adj_C - adj_B1 - adj_sum) - zeta_xy; //f_BNB(x, y)
    s_i[6] = fy - zeta_xy;  // fy is by default 0 if x == 0;
    s_i[7] = fx - zeta_xy;  // fx is by default 0 if y == 0;
    
    // Finally, divide s_i by the density (f_bzinb)
    for (int i = 0; i <8; i++) {
//cout << "si[" << i << "]: " << s_i[i] << " before / after ";      
      s_i[i] /= (l_sum * exp(adj_A - adj_B1 + adj_C - adj_sum));
//cout << s_i[i] << endl;
    }
    
    for (int i = 0; i < 8; i++) {
      for (int j = 0; j < 8; j++) {
        info[i + 8*j] += s_i[i] * s_i[j] * freq;
      }
    }
  }
}

// [[Rcpp::export]]
void dBvZINB_Expt_vec(IntegerVector &xvec, IntegerVector &yvec, IntegerVector &freq, 
                      int &n, double &a0, double &a1, double &a2,
                      double &b1, double &b2, double &p1, double &p2, 
                      double &p3, double &p4,
                      NumericVector &expt, NumericVector &s_i, NumericVector &info, int se, int bnb) {
  int sumFreq = 0, x, y, f;
  
  // initialize result
  for (int i = 0; i < 12; i++) {
    expt[i] = 0.0;
  }
  
  for (int i = 0; i < n; i++) {
    x = xvec[i]; 
    y = yvec[i];
    f = freq[i];
    sumFreq += freq[i];
    
    // add expectations with weights (freq)
    dBvZINB_Expt(x, y, f, a0, a1, a2,
                 b1, b2, p1, p2, p3, p4, expt, s_i, info, se, bnb);
  }
  
  // from sum to mean
  for (int i = 1; i < 12; i++) {  //expt[0] (likelihood) is not averaged but should be left as sum.
    expt[i] /= sumFreq;
  }
}
