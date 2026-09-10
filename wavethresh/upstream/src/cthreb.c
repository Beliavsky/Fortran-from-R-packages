#include <stdio.h>
#include <math.h>

#define PI 3.141592653589793

void Ccthrnegloglik(double *parvec, double *SigVec, double *di, double *dr, 
	long *pnd, double *pans)
{
 double p, Sig11, Sig12, Sig22, VpS11, VpS12, VpS22;
 double sum=0.0, detVpS, twopirtdetVpS, detSig, twopirtdetSig;
 double SigInv11, SigInv12, SigInv22, VpSInv11, VpSInv12, VpSInv22;
 double den1, den2, V11, V12, V22;
 int i;

 /* 
  * Evaluate the -ve log likelihood assuming a two-component mixture 
  * prior with mixing parameter p and Normal component variance matrix V. 
  *
  * Data consists of a vector of length nd; di and dr contain the 
  * imaginary and real parts of the data respectively.
  */
 
 p = parvec[0]; 
 Sig11 = SigVec[0];
 Sig12 = SigVec[1];
 Sig22 = SigVec[2];
 
 V11 = parvec[1];
 V22 = parvec[3];
 V12 = parvec[2] * sqrt(V11*V22);
 
 VpS11 = V11 + Sig11;
 VpS12 = V12 + Sig12;
 VpS22 = V22 + Sig22;
 
 detVpS = VpS11 * VpS22 - pow(VpS12, 2.0);
 twopirtdetVpS = 2.0 * PI * sqrt(detVpS);
 detSig = Sig11 * Sig22 - pow(Sig12, 2.0);
 twopirtdetSig = 2.0 * PI * sqrt(detSig);
 
 SigInv11 =  Sig22/detSig;
 SigInv12 = -Sig12/detSig;
 SigInv22 =  Sig11/detSig;
 
 VpSInv11 =  VpS22/detVpS;
 VpSInv12 = -VpS12/detVpS;
 VpSInv22 =  VpS11/detVpS;
 
 for(i = 0; i < (*pnd); i++){
   den1 = VpSInv11 * pow(dr[i], 2.0) + 2.0 *VpSInv12 * dr[i] * di[i] + 
     VpSInv22 * pow(di[i], 2.0);
   den1 = exp(-0.5 * den1) / twopirtdetVpS;
   
   den2 = SigInv11 * pow(dr[i], 2.0) + 2.0 * SigInv12 * dr[i] * di[i] + 
     SigInv22 * pow(di[i], 2.0);
   den2 = exp(-0.5 * den2) / twopirtdetSig;
   
   sum += log(p * den1 + (1.0 - p) * den2);
 } /* End for(i) */
 (*pans) = -sum;
}


void Ccthrcalcodds(long *pnd, double *dr, double *di, double *VVec,
	double *SigVec, double *pp, double *ans, double *odds)
{
  int k;
  double mult, detS, detVpS, tmp, V11, V12, V22;
  double Sig11, Sig12, Sig22, m11, m12, m22;

  /* 
   * Compute posterior weights of non-zero components given:
   *
   * nd coefficients whose real and imaginary parts are in 
   *     dr and di respectively;
   * 
   * prior and noise covariance matrices in VVec and SigVec;
   *
   * and prior weight in pp.
   *
   * Return answers in ans
   */

 Sig11 = SigVec[0];
 Sig12 = SigVec[1];
 Sig22 = SigVec[2];
 
 V11 = VVec[0];
 V12 = VVec[1];
 V22 = VVec[2];

 detS = Sig11 * Sig22 - pow(Sig12, 2.0);
 detVpS = (V11 + Sig11) * (V22 + Sig22) - pow((V12 + Sig12), 2.0);

 m11 = Sig22/detS - (V22 + Sig22)/detVpS;
 m12 = -Sig12/detS + (V12 + Sig12)/detVpS;
 m22 = Sig11/detS - (V11 + Sig11)/detVpS;

 mult = (*pp)/(1.0 - (*pp)) * sqrt(detS/detVpS);

 for(k = 0; k < (*pnd); k++){
   tmp = m11*pow(dr[k], 2.0) + 2.0 * m12 * dr[k] * di[k] + 
           m22 * pow(di[k], 2.0);
   if(tmp > 1400.0)
     tmp = 1400.0;
   odds[k] = mult * exp(tmp/2.0);
   ans[k] = odds[k] / (1 + odds[k]);
 }
}

void Cpostmean(long *pnd, double *dr, double *di, double *VVec,
	double *SigVec, double *w, double *ansr, double *ansi)
{
  int k;
  double detS, detV, tmp, V11, V12, V22, m11, m12, m22;
  double Sig11, Sig12, Sig22, SigI11, SigI12, SigI22, mi11, mi12, mi22;

  /* 
   * Compute posterior means of wavelet coefficients given:
   *
   * nd coefficients whose real and imaginary parts are in 
   *     dr and di respectively;
   * 
   * prior and noise covariance matrices in VVec and SigVec;
   *
   * posterior mixing weights in w.
   *
   * Return answers in ansr and ansi (re and im respectively).
   */

 Sig11 = SigVec[0];
 Sig12 = SigVec[1];
 Sig22 = SigVec[2];
 
 V11 = VVec[0];
 V12 = VVec[1];
 V22 = VVec[2];

 detS = Sig11 * Sig22 - pow(Sig12, 2.0);
 detV = V11 * V22 - pow(V12, 2.0);

 /* Invert Sigma */
 SigI11 = Sig22/detS;
 SigI12 = -Sig12/detS;
 SigI22 = Sig11/detS;

 /* Add Sigma^{-1} to V^{-1} */
 m11 = SigI11 + V22/detV;
 m12 = SigI12 - V12/detV;
 m22 = SigI22 + V11/detV;

 /* Now invert that sum */
 tmp = m11 * m22 - pow(m12, 2.0);
 mi11 = m22 / tmp;
 mi12 = -m12 / tmp;
 mi22 = m11 / tmp;

 for(k = 0; k < (*pnd); k++){
   ansr[k] = w[k] * (dr[k] * (mi11 * SigI11 + mi12 * SigI12) + 
                     di[k] * (mi11 * SigI12 + mi12 * SigI22));
   ansi[k] = w[k] * (dr[k] * (mi12 * SigI11 + mi22 * SigI12) + 
                     di[k] * (mi12 * SigI12 + mi22 * SigI22));
 }
}


