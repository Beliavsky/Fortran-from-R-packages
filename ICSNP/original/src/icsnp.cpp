#include <cmath>


extern "C" {

using namespace std;

  double **prepmat(double *X, int n, int k)
  //be careful to free this memory!!!
  {
    int i;
    int j;
    double **Y = new double* [n];
    for (i=0; i<n; i++) Y[i]=new double [k];
    for (i=0;i<n; i++) 
      for (j=0;j<k;j++)
    Y[i][j]=X[j*n+i];
    return Y;
  }




  void norming(double *X, int *nk, double *result)
  {
    int i;
    int j;
    int n=nk[0]; 
    int k=nk[1]; 
    for (i=0; i<n; i++)
      result[i]=0.0;
    for (i=0; i<n; i++) {
      for (j=0; j<k; j++) 
    result[i]+=X[j*n+i]*X[j*n+i];
      result[i]=sqrt(result[i]);
    }
  }

  void pairdiff(double *X, int *nk, double *result)
  {
    int n=nk[0]; 
    int k=nk[1]; 
    int i;
    int j;
    int e;
    int r=0;
    for (i=0; i<(n-1); i++)
      for (j=(i+1); j<n; j++)
    for (e=0;e<k;e++) {
      result[r]=X[e*n+i]-X[e*n+j];
      r++;}
  }
  
  void pairprod(double *X, int *nk, double *result)
  {
    int n=nk[0]; 
    int k=nk[1]; 
    int i;
    int j;
    int e;
    int r=0;
    for (i=0; i<(n-1); i++)
      for (j=(i+1); j<n; j++)
    for (e=0;e<k;e++) {
      result[r]=X[e*n+i]*X[e*n+j];
      r++;}
  }

  void pairsum(double *X, int *nk, double *result)
  {
    int n=nk[0]; 
    int k=nk[1]; 
    int i;
    int j;
    int e;
    int r=0;
    for (i=0; i<(n-1); i++)
      for (j=(i+1); j<n; j++)
    for (e=0;e<k;e++) {
      result[r]=X[e*n+i]+X[e*n+j];
      r++;}
  }

   void sum_of_sign_outers(double *X, int *nk, double *result)
  {
    int n=nk[0]; 
    int k=nk[1]; 
    int i;
    int j;
    int m;  
    int p=0;
    double *r; 
    r = new double[n];

    double **signs=prepmat(X,n,k);
    //compute norms:
    for (i=0; i<n; i++)
      r[i]=0.0;
    for (i=0; i<n; i++) {
      for (j=0; j<k; j++) 
    r[i]+=(signs[i][j]*signs[i][j]);
      r[i]=sqrt(r[i]);
    }
    
    //compute signs:
    for(i=0; i<n; i++) 
      for(j=0; j<k; j++)
    signs[i][j]=signs[i][j]/r[i];


   //prepare the result:
   for (i=0;i<(k*k);i++) result[i]=0.0;

   //compute the sum of outer products:
   for(j=0;j<k;j++)
     for(m=0;m<k;m++)
       { 
     for(i=0;i<n;i++)
       result[p]+= (signs[i][j]*signs[i][m]);
     p++;
       }
   for(i=0;i<n;i++)   
     delete [] signs[i]; 
   delete [] signs;
   delete [] r;
  }



void spatial_ranks(double *X, int *nk, double *result)
{
  int n=nk[0]; 
  int k=nk[1]; 
  int i;
  int j;
  int m;
  double d; 
  double **data=prepmat(X,n,k);

  double **ranks = new double * [n];
    for (i=0; i<n; i++) ranks[i]=new double [k];

  for(i=0;i<n;i++)
    for(m=0;m<k;m++)
      ranks[i][m]=0.0;

  double *temp = new double[k];

  for(i=0;i<n;i++) {
    for(j=0;j<n;j++) {
      if(i!=j) {
    //compute the difference and its norm:
    for(m=0;m<k;m++) 
      temp[m]=data[i][m]-data[j][m];
    d=0.0;
    for(m=0;m<k;m++)
      d+=(temp[m]*temp[m]);
    d=sqrt(d);
    //add to the sum of signs:
    for(m=0;m<k;m++)
      ranks[i][m]+=(temp[m]/d);
      }
    }
  }

  j=0;
  for(i=0;i<n;i++)
    for(m=0;m<k;m++) {
      result[j]=ranks[i][m]/n;
      j++;
    }

  for(i=0;i<n;i++) {
    delete [] data[i];
    delete [] ranks[i];
  } 
  delete [] data;
  delete [] ranks;
  delete [] temp;
}

void signed_ranks(double *X, int *nk, double *result)
{
  int n=nk[0]; 
  int k=nk[1]; 
  int i;
  int j;
  int m;
  double dm; // m as in minus
  double dp; // p as in plus
  double **data=prepmat(X,n,k);

  double **ranks = new double * [n];
    for (i=0; i<n; i++) ranks[i]=new double [k];

  for(i=0;i<n;i++)
    for(m=0;m<k;m++)
      ranks[i][m]=0.0;

  double *tempm = new double[k]; // see double dm above
  double *tempp = new double[k];

  for(i=0;i<n;i++) {
    for(j=0;j<n;j++) {
      if(i!=j) {
    //compute the difference and its norm:
    for(m=0;m<k;m++) 
      tempm[m]=data[i][m]-data[j][m];
    dm=0.0;
    for(m=0;m<k;m++)
      dm+=(tempm[m]*tempm[m]);
    dm=sqrt(dm);
    //compute the sum and its norm:
    for(m=0;m<k;m++) 
      tempp[m]=data[i][m]+data[j][m];
    dp=0.0;
    for(m=0;m<k;m++)
      dp+=(tempp[m]*tempp[m]);
    dp=sqrt(dp);
    //add to the sum of signs:
    for(m=0;m<k;m++)
      ranks[i][m]+=(tempm[m]/dm+tempp[m]/dp)/2;
      }
    }
  }

  j=0;
  for(i=0;i<n;i++)
    for(m=0;m<k;m++) {
      result[j]=ranks[i][m]/n;
      j++;
    }

  for(i=0;i<n;i++) {
    delete [] data[i];
    delete [] ranks[i];
  } 
  delete [] data;
  delete [] ranks;
  delete [] tempm;
  delete [] tempp;
}

  void outer(double *u, double *v, int k, double *uvT)
  {
    int i;
    int j;

    for(i=0;i<k;i++)
      for(j=0;j<k;j++)
    uvT[i*k+j]=u[i]*v[j];
  }

void outer2(double *u, int k, double *uut)
{
  int i;
  int j;
  
  for(i=0;i<k;i++)
    for(j=i;j<k;j++) {
      uut[i*k+j]=u[i]*u[j];
      if(j>i)
    uut[j*k+i]=uut[i*k+j];
    }
}

void touij(double *xi, double *xj, int k, double *uij)
{
  int i;
  double s=0.0;
  for(i=0;i<k;i++) {
    uij[i]= xi[i]-xj[i];
    s+=(uij[i]*uij[i]);
  }
  s=sqrt(s);
  for(i=0;i<k;i++) {
    uij[i]=uij[i]/s;
  }
}

void Q2internals(double *X, int *nk, double *result)
{
  int n=nk[0];
  int k=nk[1];
  int k2=k*k;
  int k4=k2*k2;
  int i;
  int j;
  int m;
  int ii;
  double **data=prepmat(X,n,k);
  double *uijuik = new double [k4];
  double *uij2 = new double [k2];
  double *uik2 = new double [k2];
  double *uij = new double [k];
  double *uik = new double [k];
  double *S2 = new double [k2];
  double *ave = new double [k4];


  for(i=0;i<(k2);i++) S2[i]=0.0;
  for(i=0;i<k4;i++) ave[i]=0.0;
  for (i=0;i<(n-1);i++) {
    for(j=(i+1);j<n;j++) {
      //compute u_ij
      touij(data[i],data[j],k,uij);
      //compute vec(u_ij u_ij^T)
      outer2(uij,k,uij2);
      //sum to S2
      for(ii=0;ii<k2;ii++)
    S2[ii]+=uij2[ii];

      for (m=0;m<n;m++) {
    if(m!=i) {
      touij(data[i],data[m],k,uik);
      outer2(uik,k,uik2);
      outer(uij2,uik2,k2,uijuik);
      for(ii=0;ii<k4;ii++)
        ave[ii]+=uijuik[ii];
    }
    if(m!=j) {
      touij(data[j],data[m],k,uik);
      outer2(uik,k,uik2);
      outer(uij2,uik2,k2,uijuik);
      for(ii=0;ii<k4;ii++)
        ave[ii]+=uijuik[ii];
    }
      }
    }
  }
  //scale S2 correctly, there were n*(n-1)/2 terms to sum
  for(i=0;i<k2;i++) {
    S2[i]=S2[i]/(n*(n-1)/2);
    result[i]=S2[i];
  }
  //scale ave correctly, there were n*(n-1)^2 terms to sum
  for(i=0;i<k4;i++) {
    ave[i]=ave[i]/(n*(n-1)*(n-1));
    result[k2+i]=ave[i];
  }


  for(i=0;i<n;i++) {
    delete [] data[i];
  }
  delete [] data;
  delete [] uij2;
  delete [] uik2;
  delete [] uik;
  delete [] uij;
  delete [] uijuik;
  delete [] S2;
  delete [] ave;
}



void sum_of_diff_sign_outers(double *X, int *nk, double *result)
  {
    int n=nk[0]; 
    int k=nk[1]; 
    int i;
    int j;
    int m;  
    int p=0;
    double r; 
    double *u;
    u = new double[k];
    double **Y = new double* [k];
    for (i=0; i<k; i++) Y[i]=new double [k];

    for(i=0;i<k;i++)
      for(j=0;j<k;j++)
	Y[i][j]=0.0;

    for(i=0;i<(n-1);i++)
      for(j=(i+1);j<n;j++) {
	//go over the pairs
	r=0;
	for(m=0;m<k;m++) {
	  //compute the difference and its squared norm
	  u[m] = X[m*n+i]-X[m*n+j];
	  r+=(u[m]*u[m]);
	}
	r=sqrt(r);
	for(m=0;m<k;m++)
	  //compose the sign
	  u[m]=u[m]/r;
	for(m=0;m<k;m++)
	  for(p=0;p<k;p++)
	    if(p<(m+1)) 
	      //compute the outer product elements
	      //and sum to the sum matrix
	      Y[m][p]+=(u[m]*u[p]);
      } // j, i

    for(m=0;m<(k-1);m++)
      for(p=(m+1);p<k;p++)
	//fill up the matrix
	Y[m][p]=Y[p][m];

    i=0;
    for(m=0;m<k;m++)
      for(p=0;p<k;p++) {
	//put to the result  
	result[i]=Y[m][p];
	i++;
      }
  
   
    delete [] u;
    for(i=0;i<k;i++)
      delete [] Y[i];
    delete [] Y;
  }


   void sum_of_diff_sign_select(double *X, int *nk, int *num, double *result)
  {
    int n=nk[0]; 
    int k=nk[1];
    int nu=num[0]; 
    int i;
    int j;
    int m;  
    int p=0;
    double r; 
    double *u;
    u = new double[k];
    double **Y = new double* [k];
    for (i=0; i<k; i++) Y[i]=new double [k];

    for(i=0;i<k;i++)
      for(j=0;j<k;j++)
	Y[i][j]=0.0;

    double **data=prepmat(X,n,k);

    // n-nu first
    for(i=0;i<(n-nu);i++)
      for(j=1;j<(nu+1);j++) {
	//go over the pairs
	r=0;
	for(m=0;m<k;m++) {
	  //compute the difference and its squared norm
	  u[m] = data[i][m]-data[i+j][m];
	  r+=(u[m]*u[m]);
	}
	r=sqrt(r);
	for(m=0;m<k;m++)
	  //compose the sign
	  u[m]=u[m]/r;
	for(m=0;m<k;m++)
	  for(p=0;p<k;p++)
	    if(p<(m+1)) 
	      //compute the outer product elements
	      //and sum to the sum matrix
	      Y[m][p]+=(u[m]*u[p]);
      } // j, i
    
    // the rest
    for(i=(n-nu);i<(n-1);i++)
      for(j=(i+1);j<n;j++) {
	//go over the pairs
	r=0;
	for(m=0;m<k;m++) {
	  //compute the difference and its squared norm
	  u[m] = data[i][m]-data[j][m];
	  r+=(u[m]*u[m]);
	}
	r=sqrt(r);
	for(m=0;m<k;m++)
	  //compose the sign
	  u[m]=u[m]/r;
	for(m=0;m<k;m++)
	  for(p=0;p<k;p++)
	    if(p<(m+1)) 
	      //compute the outer product elements
	      //and sum to the sum matrix
	      Y[m][p]+=(u[m]*u[p]);
      } // j, i


    for(m=0;m<(k-1);m++)
      for(p=(m+1);p<k;p++)
	//fill up the matrix
	Y[m][p]=Y[p][m];

    i=0;
    for(m=0;m<k;m++)
      for(p=0;p<k;p++) {
	//put to the result  
	result[i]=Y[m][p];
	i++;
      }
  
    for(i=0;i<n;i++)   
      delete [] data[i]; 
    delete [] data;
    delete [] u;
    for(i=0;i<k;i++)
      delete [] Y[i];
    delete [] Y;
  }

  void symm_huber(double *X, double *V, int *nk, double *cs, double *sigs, double *result)
  {
    int n=nk[0]; 
    int k=nk[1]; 
    int i;
    int j;
    int l;
    int m;  
    int p=0;
    double csq=cs[0];
    double sigsq=sigs[0];
    double rd=0.0; 
    double *u;
    u = new double[k];
    double *r;
    r = new double[k];
    double w;
    double **Y = new double* [k];
    for (i=0; i<k; i++) Y[i]=new double [k];

    for(i=0;i<k;i++)
      for(j=0;j<k;j++)
	Y[i][j]=0.0;
 
    for(i=0;i<k;i++)
       r[i]=0.0;

    for(i=0;i<(n-1);i++) 
      for(j=(i+1);j<n;j++) {
	//go over the pairs
	for(m=0;m<k;m++) {
	  //compute the difference and its squared norm
	  u[m] = X[m*n+i]-X[m*n+j];
          //compute the Mahalanobis distance
          for(l=0;l<k;l++) 
           r[l] += u[m]*V[l*k+m];
	}


       for(m=0;m<k;m++)
        rd += r[m]*u[m];

       if(rd<=csq){ w=1/sigsq;
       }else w=(csq/rd)/sigsq;
 
       for(m=0;m<k;m++)
        r[m]=0.0;
       
       rd=0.0;
        
	for(m=0;m<k;m++){
	  for(p=0;p<k;p++){
	    if(p<(m+1)) {
	      //compute the outer product elements
	      //and sum to the sum matrix
	      Y[m][p]+=w*u[m]*u[p];
            }
          }
        }   
      } // j, i

    for(m=0;m<(k-1);m++)
      for(p=(m+1);p<k;p++)
	//fill up the matrix
	Y[m][p]=Y[p][m];

    i=0;
    for(m=0;m<k;m++)
      for(p=0;p<k;p++) {
	//put to the result  
	result[i]=Y[m][p];
	i++;
      }
  
    delete [] r;
    delete [] u;
    for(i=0;i<k;i++)
      delete [] Y[i];
    delete [] Y;
  }


}
