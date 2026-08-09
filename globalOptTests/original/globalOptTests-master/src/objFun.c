#include <math.h>
#define PI 3.14159265359

void Ackleys(double *sum, double *x, int *N)
/*Calculate objective function value of x[].*/
{
  int i,nn;
  double subsum1=0., subsum2=0.;
 
  for (i=0; i<N[0]; i++)
    {
      subsum1 = subsum1 + pow(x[i],2);
      subsum2 = subsum2 + cos(2*PI*x[i]);
    }
  nn = N[0];
  /* note that GAMS uses factor of -0.02 not 0.2  */
  sum[0] = -20.*exp(-.2*sqrt(subsum1/(double)nn))-exp(subsum2/(double)nn)+20.+exp(1);
}
void AluffiPentini(double *sum, double *x, int *N)
/*Calculate objective function value of x[].*/
{
  sum[0] = .25*pow(x[0],4) -.5*x[0]*x[0] + .1*x[0] +.5*x[1]*x[1];

}
void BeckerLago(double *sum, double *x, int *N)
/*Calculate objective function value of x[].*/
{
  sum[0] = pow(fabs(x[0])-5.,2) + pow(fabs(x[1])-5.,2);
}


void Bohachevsky1(double *sum, double *x, int *N)
/*Calculate objective function value of x[].*/
{
  sum[0] = pow(x[0],2) + 2.*x[1]*x[1] - .3*cos(3*PI*x[0]) -.4*cos(4*PI*x[1]) + .7;
}

void Bohachevsky2(double *sum, double *x, int *N)
/*Calculate objective function value of x[].*/
{
  sum[0] = pow(x[0],2) + 2.*x[1]*x[1] - .3*cos(3*PI*x[0])*cos(4*PI*x[1]) + .3;
}

void Branin(double *sum, double *x, int *N)
/*Calculate objective function value of x[].*/
{

  double a=1., b=5.1/(4*PI*PI), c=5/PI;
  double d=6., e=10., f=1/(8*PI);

  sum[0] = a*pow((x[1]-b*x[0]*x[0]+c*x[0]-d),2) + e*(1-f)*cos(x[0]) + e;
}

void Camel3(double *sum, double *x, int *N)
/*Calculate objective function value of x[].*/
{
  sum[0] = (2 - 1.05*x[0]*x[0] + pow(x[0],4)/6)*x[0]*x[0] + x[0]*x[1] + x[1]*x[1];
}

void Camel6(double *sum, double *x, int *N)
/*Calculate objective function value of x[].*/
{
  sum[0] = (4 - 2.1*x[0]*x[0] + pow(x[0],4)/3)*x[0]*x[0] + x[0]*x[1] + (-4 + 4*x[1]*x[1])*x[1]*x[1];
}

void CosMix2(double *sum, double *x, int *N)
/*Calculate objective function value of x[].*/
{
  int i;
  double sum1=0., sum2=0;


  for (i=0; i<N[0]; i++)
    {
      sum1 = sum1 + cos(5*PI*x[i]);
      sum2 = sum2 + pow(x[i],2);
    }
  sum[0] = 0.1*sum1 - sum2;
  sum[0]=-sum[0];
}

void CosMix4(double *sum, double *x, int *N)
/*Calculate objective function value of x[].*/
{
  int i;
  double sum1=0., sum2=0;


  for (i=0; i<N[0]; i++)
    {
      sum1 = sum1 + cos(5*PI*x[i]);
      sum2 = sum2 + pow(x[i],2);
    }
  sum[0] = 0.1*sum1 - sum2;
  sum[0]=-sum[0];

}

void DekkersAarts(double *sum, double *x, int *N)
/*Calculate objective function value of x[].*/
{
  sum[0] = 100000.*x[0]*x[0] + x[1]*x[1] - pow((x[0]*x[0]+x[1]*x[1]),2) +pow((x[0]*x[0]+x[1]*x[1]),4)/100000;



}


void Easom(double *sum, double *x, int *N)
/*Calculate objective function value of x[].*/
{

  sum[0] = -cos(x[0])*cos(x[1])*exp(-pow(x[0]-PI,2)-pow(x[1]-PI,2));



}


void EMichalewicz(double *sum, double *x, int *N)
/*Calculate objective function value of x[].*/
{

  int j, nn;
  double y[10];
  double cost,sint,m=10.;
  nn = N[0];

  cost = cos(PI/6.);
  sint = sin(PI/6.);

  for (j=0; j<nn-1; j+=2 )		/* Corrects errors in original
					   ICEO test bed */
    {
      y[j] = x[j]*cost - x[j+1]*sint;
      y[j+1] = x[j]*sint + x[j+1]*cost;
    }

  if (j==nn-1)  y[j]=x[j];

  for (j=0; j<N[0]; j++)
    {
      sum[0] -= sin(y[j]) * pow(sin((j+1)*y[j]*y[j]/PI),2.0*m);

    }

  //        for (j=0; j<N; j++) printf(" %f", y[j]);
  //for (j=0; j<N; j++) printf(" %f", x[j]);
  //printf(" %f", sum);

}

void Expo(double *sum, double *x, int *N)
/*Calculate objective function value of x[].*/
{
  int j;

  for (j=0; j<N[0]; j++)
    {
      sum[0]+=x[j]*x[j];
    }
  sum[0] = -exp(-.5*sum[0]);



}

void GoldPrice(double *sum, double *x, int *N)
/*Calculate objective function value of x[].*/
{

  sum[0] = (1+(x[0]+x[1]+1)*(x[0]+x[1]+1)*(19-14*x[0]+3*x[0]*x[0]-14*x[1]+6*x[0]*x[1]+3*x[1]*x[1]));
  sum[0] = sum[0]*(30+(2*x[0]-3*x[1])*(2*x[0]-3*x[1])*(18-32*x[0]+12*x[0]*x[0]+48*x[1]-36*x[0]*x[1]+27*x[1]*x[1]));
  
  

}


void Griewank(double *sum, double *x, int *N)
/*Calculate objective function value of x[].*/
{
  int j;
  double prod=1.;

  for (j=0; j<N[0]; j++)
    {
      sum[0]+=x[j]*x[j];
      prod*=cos(x[j]/sqrt((double)(j+1)));
    }
  sum[0]=sum[0]/4000.-prod+1.;



}


void Gulf(double *sum, double *x, int *N)
/*Calculate objective function value of x[].*/
{
  int j;
  double subsum=0.,u;

  for (j=0; j<99; j++)
    {
      u=25.+pow(-50*log(.01*j),.66666);
      subsum=exp(-pow(u-x[1],x[2])/x[0])-0.01*j; 
      sum[0]+= pow(subsum,2);
      //                  printf("j %d %f %f\n", j,u,subsum);
    }



}


void Hartman3(double *sum, double *x, int *N)
/*Calculate objective function value of x[].*/
{
  int i,j;
  double subsum,dist=0.;
  static double a[4][3] = {
    {3, 10, 30},
    {.1, 10, 35},
    {3, 10, 30},
    {.1, 10, 35}};
  
  static double c[5] = {1, 1.2, 3, 3.2};
  static double p[4][3] = {
    {.3689, .117, .2673},
    {.4699, .4387, .747},
    {.1091, .8732, .5547},
    {.03815, .5743, .8828}};

  sum[0]=0.;
  for (i=0; i<5; i++)
    {
      subsum=0.;
      for (j=0; j<N[0]; j++)
	{
	  subsum = subsum - a[i][j]*pow((x[j]-p[i][j]),2);
	}
      sum[0] = sum[0] - c[i]*exp(subsum);
    }

}

void Hartman6(double *sum, double *x, int *N)
/*Calculate objective function value of x[].*/
{
  int i,j;
  double subsum,dist=0.;
  static double a[4][6] = {
    {10, 3, 17, 3.5, 1.7, 8},
    {.05, 10, 17, .1, 8, 14},
    {3, 3.5, 1.7, 10, 17, 8},
    {17, 8, .05, 10, .1, 14}};
  
  static double c[5] = {1, 1.2, 3, 3.2};
  static double p[4][6] = {
    {.1312, .1696, .5569, .0124, .8283, .5886},
    {.2329, .4135, .8307, .3736, .1004, .9991},
    {.2348, .1451, .3522, .2883, .3047, .6650},
    {.4047, .8828, .8732, .5743, .1091, .0381}};

  sum[0]=0.;
  for (i=0; i<5; i++)
    {
      subsum=0.;
      for (j=0; j<N[0]; j++)
	{
	  subsum = subsum - a[i][j]*pow((x[j]-p[i][j]),2);
	}
      sum[0] = sum[0] - c[i]*exp(subsum);
    }
}

void Hosaki(double *sum, double *x, int *N)
/*Calculate objective function value of x[].*/
{

  sum[0] = (1. - 8.*x[0] + 7.*x[0]*x[0] - 7./3*x[0]*x[0]*x[0] + 1./4*x[0]*x[0]*x[0]*x[0]);
  sum[0] = sum[0]*x[1]*x[1]*exp(-x[1]);

}


void Kowalik(double *sum, double *x, int *N)
/*Calculate objective function value of x[].*/
{
  int i;
  static double a[11] = {.1957, .1947, .1735, .16, .0844, .0627, .0456, .0342, .0323, .0235, .0246};

  static double b[11] = {.25, .5, 1, 2, 4, 6, 8, 10, 12, 14, 16};
  sum[0]=0.;
 
  for (i=0; i<11; i++){
    sum[0] = sum[0] + pow((a[i]-x[0]*(1+x[1]*b[i])/(1+x[2]*b[i]+x[3]*b[i]*b[i])),2);
  }



}


void LM1(double *sum, double *x, int *N)
/*Calculate objective function value of x[].*/
{
  int j, nn;
  double zw1, zw2;
  nn = N[0];
  zw1 = 10*pow(sin(PI*(1+.25*(x[0]+1))),2);
  zw2 = pow(.25*(x[(int)nn-1]+1),2);

  for (j=0; j<(nn-1); j++)
    {
      zw1 += pow(.25*(x[j]+1),2)*(1+pow(sin(PI*(.25*x[j+1])),2));
    }
  
  sum[0] = PI/(double)nn*(zw1+zw2);

}

void LM2n10(double *sum, double *x, int *N)
/*Calculate objective function value of x[].*/
{
  int j,nn;
  double zw1, zw2;

  nn = N[0];
  zw1 = pow(sin(3*PI*x[0]),2);
  zw2 = pow(x[nn-1]-1,2)*(1+pow(sin(2*PI*x[nn-1]),2));

  for (j=0; j<(nn-1); j++)
    {
      zw1 += pow(x[j]-1,2)*(1+pow(sin(3*PI*x[j+1]),2));
    }
  sum[0] = .1*(zw1+zw2);



}


void LM2n5(double *sum, double *x, int *N)
/*Calculate objective function value of x[].*/
{
  int j, nn;
  double zw1, zw2;
  nn = N[0];
  zw1 = pow(sin(3*PI*x[0]),2);
  zw2 = pow(x[nn-1]-1,2)*(1+pow(sin(2*PI*x[nn-1]),2));

  for (j=0; j<(nn-1); j++)
    {
      zw1 += pow(x[j]-1,2)*(1+pow(sin(3*PI*x[j+1]),2));
    }
  sum[0] = .1*(zw1+zw2);



}

void McCormic(double *sum, double *x, int *N)
/*Calculate objective function value of x[].*/
{
  sum[0] = sin((x[0]+x[1])) + pow(x[0]-x[1],2) - 1.5*x[0] + 2.5*x[1] + 1.;  
	
  

}


void MeyerRoth(double *sum, double *x, int *N)
/*Calculate objective function value of x[].*/
{
  int i;
  double num=0., den=0.;
  static double t[5] = {1., 2., 1., 2., .1};
  static double v[5] = {1., 1., 2., 2., 0.};
  static double y[5] = {.126, .219, .076, .126, .186};
  sum[0]=0.;
  for (i=0; i<5; i++)
    {
      num = x[0]*x[2]*t[i];
      den = 1.+x[0]*t[i]+x[1]*v[i];
      sum[0]+= pow(num/den-y[i],2);
    }


}

void MieleCantrell(double *sum, double *x, int *N)
/*Calculate objective function value of x[].*/
{
  double sum1=0., sum2=0.;

  sum1 = pow(exp(x[0])-x[1],4) + 100.*pow(x[1]-x[2],6);
  sum2 = pow(tan(x[2]-x[3]),4) + pow(x[0],8);

  sum[0] = sum1+sum2;



}

void Modlangerman(double *sum, double *x, int *N)
/*Calculate objective function value of x[].*/
{
  int i,j;
  double dx,dist;
  static double a[5][10] = {
    {9.681, 0.667, 4.783, 9.095, 3.517, 9.325, 6.544, 0.211, 5.122, 2.020},
    {9.400, 2.041, 3.788, 7.931, 2.882, 2.672, 3.568, 1.284, 7.033, 7.374},
    {8.025, 9.152, 5.114, 7.621, 4.564, 4.711, 2.996, 6.126, 0.734, 4.982},
    {2.196, 0.415, 5.649, 6.979, 9.510, 9.166, 6.304, 6.054, 9.377, 1.426},
    {8.074, 8.777, 3.467, 1.867, 6.708, 6.349, 4.534, 0.276, 7.633, 1.567}};
 
  static double c[5] = {0.806,0.517,0.1,0.908,0.965};

  for (i=0; i<5; i++)
    {
      dist=0;
      for (j=0; j<N[0]; j++)
	{
	  dx=x[j]-a[i][j];
	  dist+=dx*dx;
	}
      sum[0]-=c[i]*(exp(-dist/PI)*cos(PI*dist));
    }



}

void ModRosenbrock(double *sum, double *x, int *N)
/*Calculate objective function value of x[].*/
{
  sum[0]=100.*pow((x[1]-x[0]*x[0]),2) + pow((6.4*pow(x[1]-0.5,2) - x[0] - 0.6),2);



}

void MultiGauss(double *sum, double *x, int *N)
/*Calculate objective function value of x[].*/
{
  int i;
  static double a[5] = {.5, 1.2, 1., 1., 1.2};
  static double b[5] = {0., 1., 0., -.5, 0.};
  static double c[5] = {0., 0., -.5, 0., 1.};
  static double d[5] = {.1, .5, .5, .5, .5};
  sum[0]=0.;


  for (i=0; i<5; i++)
    {
      sum[0] = sum[0] - a[i]*exp(-(pow(x[0]-b[i],2)+pow((x[1]-c[i]),2))/pow(d[i],2));
    }



}

void Neumaier2(double *sum, double *x, int *N)
/*Calculate objective function value of x[].*/
{
  int i,k;
  double sum1=0.0;
  static double b[4] = {8.0, 18.0, 44.0, 114.0};
  sum[0]=0.;
  for (k=0; k<N[0]; k++)
    {
      sum1 = 0.;
      for (i=0; i<N[0]; i++)
	{
	  sum1+= pow(x[i],k+1);
	}
      sum[0]+= pow(b[k]-sum1,2);
    }



}

void Neumaier3(double *sum, double *x, int *N)
/*Calculate objective function value of x[].*/
{
  int i;
  double obj1=0.0, obj2=0.0;

  for (i=0; i<N[0]; i++)
    {
      obj1 += pow(x[i]-1.0,2);
    }
  for (i=1; i<N[0]; i++)
    {
      obj2 += x[i]*x[i-1];
    }

  sum[0] = obj1-obj2;



}

void Paviani(double *sum, double *x, int *N)
/*Calculate objective function value of x[].*/
{
  int j;
  double sum1=0.0, prod1=1.0, prod=1.0;
	
  for (j=0; j<N[0]; j++)
    {
      prod1 = prod1*x[j];
      sum1 = sum1 + pow((double) log(x[j]-2),2) + pow((double) log(10-x[j]),2);
    }

  prod = pow((double) prod1,0.2);
  sum[0] = sum1-prod;

}


void Periodic(double *sum, double *x, int *N)
/*Calculate objective function value of x[].*/
{

  sum[0] = 1. + pow(sin(x[0]),2) + pow(sin(x[1]),2) - 0.1*exp(-x[0]*x[0] - x[1]*x[1]);



}

void PowellQ(double *sum, double *x, int *N)
/*Calculate objective function value of x[].*/
{

  sum[0] = (x[0]+10.*x[0])*(x[0]+10.*x[0]) + 5.*(x[2]-x[3])*(x[2]-x[3]);
  sum[0] = sum[0] + pow((x[1]-2.*x[2]),4) + 10.*pow((x[0]-x[3]),4);

}


void PriceTransistor(double *sum, double *x, int *N)
/*Calculate objective function value of x[].*/
{
  int k;
  double sumsqr=0.0, alpha, beta;
  static double g[5][4] = {
    {0.485, 0.752, 0.869, 0.982},
    {0.369, 1.254, 0.703, 1.455},
    {5.2095, 10.0677, 22.9274, 20.2153},
    {23.3037, 101.779, 111.461, 191.267},
    {28.5132, 111.8467, 134.3884, 211.4823}};

  for (k=0; k<4; k++)
    {
      alpha = (1.0-x[0]*x[1])*x[2]*(exp(x[4]*(g[0][k]-0.001*g[2][k]*x[6]-0.001*x[7]*g[4][k]))-1.0) - g[4][k] + g[3][k]*x[1];
      beta = (1.0-x[0]*x[1])*x[3]*(exp(x[5]*(g[0][k]-g[1][k]-0.001*g[2][k]*x[6]+g[3][k]*0.001*x[8]))-1.0)- g[4][k]*x[0] + g[3][k];
      sumsqr += alpha*alpha + beta*beta;
    }
  sum[0] = pow(x[0]*x[2] - x[1]*x[3],2) + sumsqr;



}


void Rastrigin(double *sum, double *x, int *N)
/*Calculate objective function value of x[].*/
{

  int j;
  sum[0]=0;

  for (j=0; j<N[0]; j++)
    {
      sum[0]+=x[j]*x[j]-10.*cos(2.*PI*x[j])+10.;
    }



}


void Rosenbrock(double *sum, double *x, int *N)
/*Calculate objective function value of x[].*/
{
  int j, nn;
  double a=0.,b=0.;
  sum[0]=0.;
  nn = N[0];
  for (j=0; j<nn-1; j++)
    { 
      a=x[j]*x[j]-x[j+1];
      b=1.-x[j];
      sum[0]+=100.*a*a+b*b;	  
    }
  
  
}

void Salomon(double *sum, double *x, int *N)
/*Calculate objective function value of x[].*/
{
  int j;
  sum[0]=0;

  for (j=0; j<N[0]; j++)
    { 
      sum[0]+=x[j]*x[j];
    }
  sum[0]=sqrt(sum[0]);
  sum[0] = -cos(2.*PI*sum[0])+.1*sum[0]+1.;

  
}


void Schaffer1(double *sum, double *x, int *N)
/*Calculate objective function value of x[].*/
{
  double num=0.,den=0.;

  num = pow((sin(sqrt(x[0]*x[0]+x[1]*x[1]))),2) - 0.5;
  den = pow((1+.001*(x[0]*x[0]+x[1]*x[1])),2);
  sum[0] = 0.5 + num/den;



}


void Schaffer2(double *sum, double *x, int *N)
/*Calculate objective function value of x[].*/
{
  double prod1=0., prod2=0.;

  prod1 = pow(x[0]*x[0]+x[1]*x[1],0.25);
  prod2 = pow(50*(x[0]*x[0]+x[1]*x[1]),0.1);
  sum[0] = prod1*(sin(sin(prod2))+1.0);



}


void Schubert(double *sum, double *x, int *N)
/*Calculate objective function value of x[].*/
{
  int i,j;
  double prod=1.;


  for (i=0; i<N[0]; i++)
    {
      sum[0] = 0.;
      for (j=1; j<=5; j++)
	sum[0] = sum[0] + j*cos((j+1)*x[i]+j);
      prod = prod*sum[0];
    }
  sum[0] = prod;

}


void Schwefel(double *sum, double *x, int *N)
/*Calculate objective function value of x[].*/
{
  int j;
  sum[0]=0.;

  for (j=0; j<N[0]; j++)
    { 
      sum[0]+=x[j]*sin(sqrt(fabs(x[j])));	  
    }
  sum[0] = - sum[0];
  
}

void Shekel10(double *sum, double *x, int *N)
/*Calculate objective function value of x[].*/
{
  int i,j;
  double den;
  static double a[10][4] = {
    {4, 4, 4, 4},{1, 1, 1, 1},
    {8, 8, 8, 8},{6, 6, 6, 6},
    {3, 7, 3, 7},{2, 9, 2, 9},
    {5, 5, 3, 3},{8, 1, 8, 1},
    {6, 2, 6, 2},{7, 3.6, 7, 3.6}};
  static double c[10]={.1,.2,.2,.4,.4,.6,.3,.7,.5,.5};
  sum[0]=0.;
  for (i=0; i<10; i++) {
    den = 0.;
    for (j=0; j<N[0]; j++) {
      den = den + pow((x[j]-a[i][j]),2);
    }
    sum[0] = sum[0] - 1.0/(den + c[i]);
  }


}
	

void Shekel5(double *sum, double *x, int *N)
/*Calculate objective function value of x[].*/
{
  int i,j;
  double den;
  static double a[5][4] = {{4, 4, 4, 4},{1, 1, 1, 1},{8, 8, 8, 8},{6, 6, 6, 6},{3, 7, 3, 7}};
  static double c[5]={0.1,0.2,0.2,0.4,0.4};
  sum[0]=0.;
  for (i=0; i<5; i++) {
    den = 0.;
    for (j=0; j<N[0]; j++) {
      den = den + pow((x[j]-a[i][j]),2);
    }
    sum[0] = sum[0] - 1.0/(den + c[i]);
  }


}

void Shekel7(double *sum, double *x, int *N)
/*Calculate objective function value of x[].*/
{
  int i,j;
  double den;
  static double a[7][4] = {
    {4, 4, 4, 4},{1, 1, 1, 1},
    {8, 8, 8, 8},{6, 6, 6, 6},
    {3, 7, 3, 7},{2, 9, 2, 9},
    {5, 5, 3, 3}};
  static double c[7]={.1,.2,.2,.4,.4,.6,.3};
  sum[0]=0.;
  for (i=0; i<7; i++) {
    den = 0.;
    for (j=0; j<N[0]; j++) {
      den = den + pow((x[j]-a[i][j]),2);
    }
    sum[0] = sum[0] - 1.0/(den + c[i]);
  }


}
	

void Shekelfox5(double *sum, double *x, int *N)
/*Calculate objective function value of x[].*/
{

  int i,j;
  double sp, h;
  static double a[30][10] = {
    {9.681, 0.667, 4.783, 9.095, 3.517, 9.325, 6.544, 0.211, 5.122, 2.020},
    {9.400, 2.041, 3.788, 7.931, 2.882, 2.672, 3.568, 1.284, 7.033, 7.374},
    {8.025, 9.152, 5.114, 7.621, 4.564, 4.711, 2.996, 6.126, 0.734, 4.982},
    {2.196, 0.415, 5.649, 6.979, 9.510, 9.166, 6.304, 6.054, 9.377, 1.426},
    {8.074, 8.777, 3.467, 1.863, 6.708, 6.349, 4.534, 0.276, 7.633, 1.567},
    {7.650, 5.658, 0.720, 2.764, 3.278, 5.283, 7.474, 6.274, 1.409, 8.208},
    {1.256, 3.605, 8.623, 6.905, 4.584, 8.133, 6.071, 6.888, 4.187, 5.448},
    {8.314, 2.261, 4.224, 1.781, 4.124, 0.932, 8.129, 8.658, 1.208, 5.762},
    {0.226, 8.858, 1.420, 0.945, 1.622, 4.698, 6.228, 9.096, 0.972, 7.637},
    {7.305, 2.228, 1.242, 5.928, 9.133, 1.826, 4.060, 5.204, 8.713, 8.247},
    {0.652, 7.027, 0.508, 4.876, 8.807, 4.632, 5.808, 6.937, 3.291, 7.016},
    {2.699, 3.516, 5.874, 4.119, 4.461, 7.496, 8.817, 0.690, 6.593, 9.789},
    {8.327, 3.897, 2.017, 9.570, 9.825, 1.150, 1.395, 3.885, 6.354, 0.109},
    {2.132, 7.006, 7.136, 2.641, 1.882, 5.943, 7.273, 7.691, 2.880, 0.564},
    {4.707, 5.579, 4.080, 0.581, 9.698, 8.542, 8.077, 8.515, 9.231, 4.670},
    {8.304, 7.559, 8.567, 0.322, 7.128, 8.392, 1.472, 8.524, 2.277, 7.826},
    {8.632, 4.409, 4.832, 5.768, 7.050, 6.715, 1.711, 4.323, 4.405, 4.591},
    {4.887, 9.112, 0.170, 8.967, 9.693, 9.867, 7.508, 7.770, 8.382, 6.740},
    {2.440, 6.686, 4.299, 1.007, 7.008, 1.427, 9.398, 8.480, 9.950, 1.675},
    {6.306, 8.583, 6.084, 1.138, 4.350, 3.134, 7.853, 6.061, 7.457, 2.258},
    {0.652, 2.343, 1.370, 0.821, 1.310, 1.063, 0.689, 8.819, 8.833, 9.070},
    {5.558, 1.272, 5.756, 9.857, 2.279, 2.764, 1.284, 1.677, 1.244, 1.234},
    {3.352, 7.549, 9.817, 9.437, 8.687, 4.167, 2.570, 6.540, 0.228, 0.027},
    {8.798, 0.880, 2.370, 0.168, 1.701, 3.680, 1.231, 2.390, 2.499, 0.064},
    {1.460, 8.057, 1.336, 7.217, 7.914, 3.615, 9.981, 9.198, 5.292, 1.224},
    {0.432, 8.645, 8.774, 0.249, 8.081, 7.461, 4.416, 0.652, 4.002, 4.644},
    {0.679, 2.800, 5.523, 3.049, 2.968, 7.225, 6.730, 4.199, 9.614, 9.229},
    {4.263, 1.074, 7.286, 5.599, 8.291, 5.200, 9.214, 8.272, 4.398, 4.506},
    {9.496, 4.830, 3.150, 8.270, 5.079, 1.231, 5.731, 9.494, 1.883, 9.732},
    {4.138, 2.562, 2.532, 9.661, 5.611, 5.500, 6.886, 2.341, 9.699,
     6.500}};

  static double c[30] = {
    0.806,
    0.517,
    0.1,
    0.908,
    0.965,
    0.669,
    0.524,
    0.902,
    0.531,
    0.876,
    0.462,
    0.491,
    0.463,
    0.714,
    0.352,
    0.869,
    0.813,
    0.811,
    0.828,
    0.964,
    0.789,
    0.360,
    0.369,
    0.992,
    0.332,
    0.817,
    0.632,
    0.883,
    0.608,
    0.326};
  sum[0]=0.;
  for (j=0; j<30; j++)
    {
      sp=0.0;
      for (i=0; i<N[0]; i++)
	{
	  h=x[i]-a[j][i];
	  sp+=h*h;
	}
      sum[0]-=1.0/(sp+c[j]);
    }



}

/* void STChebychev9(double *sum, double *x, int *N) */
/* /\*Calculate objective function value of x[].*\/ */
/* { */
/*   int i, j, m=60; */
/*   double d=72.661; */
/*   double u=0.0, v=0.0, w=0.0; */
/*   double p1, p2, p3=0.0; */
/*   double pj; */
 
/*   for (i=0; i<9; i++) u+= pow(1.2,(double)9.0-i-1)*x[i]; */
/*   if (u < d) */
/*     { */
/*       p1 = (u-d)*(u-d); */
/*     } */
/*   else p1 = 0; */

/*   for (i=0; i<9; i++) v+= pow(-1.2,(double)9.0-i-1)*x[i]; */
/*   if (v < d) */
/*     { */
/*       p2 = (v-d)*(v-d); */
/*     } */
/*   else p2 = 0; */

/*   for (j=0; j<=m; j++)  */
/*     { */
/*       w = 0.0; */
/*       for (i=0; i<9; i++) { */
/* 	w+= pow(2.0/m*j,(double)9.0-i-1)*x[i]; */
/* 	//                  printf(">> %f %f\n", 2.0/m*j, pow(2.0/m*j,&N-i-1)); */
/*       } */
/*       if (w>1) pj = (w-1)*(w-1); */
/*       else if (w<-1) pj = (w+1)*(w+1); */
/*       else pj = 0.0; */
/*       //                printf("j %d %f %f\n", j, w, pj); */
/*       p3 += pj; */
/*     } */

/*   sum[0] = p1+p2+p3; */
/*   //        printf("%f %f %f %f %f\n", u, p1, v, p2, p3); */

    

/* } */

/* void STChebychev17(double *sum, double *x, int *N) */
/* /\*Calculate objective function value of x[].*\/ */
/* { */
/*   int i, j,  m=100; */
/*   double d=10558.145; */
/*   double u=0.0, v=0.0, w=0.0; */
/*   double p1, p2, p3=0.0; */
/*   double pj; */
 
/*   for (i=0; i<17; i++) u+= pow(1.2,(double)17.0-i-1)*x[i]; */
/*   if (u < d) */
/*     { */
/*       p1 = (u-d)*(u-d); */
/*     } */
/*   else p1 = 0; */

/*   for (i=0; i<17; i++) v+= pow(-1.2,(double)17.0-i-1)*x[i]; */
/*   if (v < d) */
/*     { */
/*       p2 = (v-d)*(v-d); */
/*     } */
/*   else p2 = 0; */

/*   for (j=0; j<=m; j++)  */
/*     { */
/*       w = 0.0; */
/*       for (i=0; i<17; i++) w+= pow(2.0/m*j,(double)17.0-i-1)*x[i]; */
/*       if (w>1) pj = (w-1)*(w-1); */
/*       else if (w<-1) pj = (w+1)*(w+1); */
/*       else pj = 0.0; */
/*       p3 += pj; */
/*     } */

/*   sum[0] = p1+p2+p3; */

    

/* } */


void Wood(double *sum, double *x, int *N)
/*Calculate objective function value of x[].*/
{
  sum[0] = 100*pow((x[1]-x[0]*x[0]),2) + (1-x[0])*(1-x[0])
    +90*pow((x[3]-x[2]*x[2]),2) + (1-x[2])*(1-x[2])
    +10.1*(pow((x[1]-1),2)+pow((x[3]-1),2))
    +19.8*(x[1]-1)*(x[3]-1);


}


void Zeldasine10(double *sum, double *x, int *N)
/*Calculate objective function value of x[].*/
{
  int j;
  double prod1=1., prod2=1.;
  double A=2.5, B=5., z=PI/6;

  for (j=0; j<N[0]; j++)
    {
      prod1 *= sin(x[j]-z);
      prod2 *= sin(B*(x[j]-z));
    }
  sum[0] = -(A*prod1 + prod2);

}


void Zeldasine20(double *sum, double *x, int *N)
/*Calculate objective function value of x[].*/
{
  int j;
  double prod1=1.0, prod2=1.0;
  double A=2.5, B=5.0, z=PI/6;

  for (j=0; j<N[0]; j++)
    {
      prod1 *= sin(x[j]-z);
      prod2 *= sin(B*(x[j]-z));
    }
  sum[0] = -1*(A*prod1 + prod2);


}
