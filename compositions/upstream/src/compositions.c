#include <string.h>
#include <R.h>
#include <Rdefines.h>
#include <Rinternals.h>
#include <limits.h>     
#include <R_ext/Utils.h>
#define USE_FC_LEN_T
#include <Rconfig.h>
#include <R_ext/BLAS.h>
#include <R_ext/Lapack.h>
#include <R_ext/Rdynload.h>
#ifndef FCONE
# define FCONE
#endif

#define isNil(x) ((x)==NULL || (SEXP) (x) == R_NilValue )


/* Geostatistics */
static double sq(double x) {return x*x;}
#define Tvg(s,i,j) vg[(s)+nbins*((i)+D*(j))]
#define Tn(s,i,j) n[(s)+nbins*((i)+D*(j))]
#define Th(s,i,j) h[(s)+nbins*((i)+D*(j))]
#define MX(i,k) X[(i)+N*(k)]
#define MZ(i,k) Z[(i)+N*(k)]


// variogram
static R_NativePrimitiveArgType gsiCGSvariogram_t[] = {
  INTSXP,REALSXP,INTSXP,REALSXP,INTSXP,REALSXP,REALSXP,REALSXP,REALSXP,REALSXP,INTSXP
};

extern void gsiCGSvariogram(
    int *dimZ,  
    double *Z,
    int *dimX,
    double *X,
    int *nbinsP,
    double *bins,
    double *aziVec,
    double *cosMin,
    double *vg,
    double *h,
    int    *n
    ) {
    const int nbins=*nbinsP;
    const int N=dimZ[0];
    const int D=dimZ[1];
    const int S=dimX[1];
    int i,j,s,k,l;
    // Zero arrays
    for(i=0;i<D;i++)
	for(j=0;j<D;j++) 
	    for(s=0;s<nbins;s++) 
		Tvg(s,i,j)=0.0;
    for(i=0;i<D;i++)
	for(j=0;j<D;j++) 
	    for(s=0;s<nbins;s++) 
		Tn(s,i,j)=0.0;
    for(i=0;i<D;i++)
	for(j=0;j<D;j++) 
	    for(s=0;s<nbins;s++) 
		Th(s,i,j)=0.0;
    // Count 
    for(i=0;i<N;i++)
	for(j=0;j<N;j++) {
	    double h2=0.0;
	    for(k=0;k<S;k++)
		h2+= sq(MX(i,k)-MX(j,k));
	    double hh=sqrt(h2);
	    if(*cosMin > -1) {
		h2=0.0;
		for(k=0;k<S;k++)
		    h2+= (MX(i,k)-MX(j,k))*aziVec[k];
		if( h2/hh < *cosMin )
		    break;
	    }
	    for(s=0;s<nbins;s++) {
		if( hh>bins[s] && hh<=bins[s+nbins]) {
		    for(k=0;k<D;k++) {
			double a1=MZ(i,k);
			double a2=MZ(j,k);
			if( R_finite(a1) && R_finite(a2))
			    for(l=0;l<D;l++) {
				double b1=MZ(i,l);
				double b2=MZ(j,l);
				if( R_finite(b1) && R_finite(b2)) {
				    Tn(s,k,l)++;
				    Tvg(s,k,l)+=sq((a1-b1)-(a2-b2));
				    Th(s,k,l)+=hh;
				}
			    }
		    }
		} 
	    }
	}
    // Rescale
    for(i=0;i<D;i++)
	for(j=0;j<D;j++) 
	    for(s=0;s<nbins;s++) 
		Tvg(s,i,j)/=Tn(s,i,j);
    for(i=0;i<D;i++)
	for(j=0;j<D;j++) 
	    for(s=0;s<nbins;s++) 
		Th(s,i,j)/=Tn(s,i,j);
    
    
}
#undef Tvg
#undef TN
#undef MX

#define Tvg(k,i,j) vg[(k)+N*((i)+D*(j))]
#define Tlr(k,i,j) lr[(k)+N*((i)+D*(j))]

// vg2lrvg
static R_NativePrimitiveArgType gsiCGSvg2lrvg_t[] = {
  INTSXP,REALSXP,REALSXP
};



extern void gsiCGSvg2lrvg(
    int *dimVg,
    double *vg,
    double *lr
    ) {
    int i,j,k;
    const int D=dimVg[1];
    const int N=dimVg[0];
    if( dimVg[2]!=D )
	error("gsiCGSvg2lrvg: wrong dimensions");
    for(k=0;k<N;k++) {
	for(i=0;i<D;i++)
	    for(j=0;j<D;j++) {
		Tlr(k,i,j) = Tvg(k,i,i)+Tvg(k,j,j)-Tvg(k,i,j)-Tvg(k,j,i);
	    }
    }
}

void printMat(char *name,char *format,int n,int m,double *A,int lda) {
    Rprintf("%s \n",name);
    int i,j;
    for(i=0;i<n;i++) {
	for(j=0;j<m;j++) {
	    Rprintf(format,A[i+lda*j]);
	    Rprintf("\t");
	}
	Rprintf("\n");
    }
}

void printBlockMat(char *name,char *format,int n,int b,double *A) {
    Rprintf("%s \n",name);
    int i,j,l,k;
    for(i=0;i<n;i++) 
	for( k=0;k<b;k++){
	    for(j=0;j<n;j++)
		for(l=0;l<b;l++) {
		    Rprintf(format,A[i+n*(j+n*(l+b*k))]);
		    Rprintf("\t");
		}
	Rprintf("\n");
	}
}


void checkBlockMatSymmetry(int n,int b,double *x,double tol) {
    int i,j,k,l;
    for(i=0;i<n;i++)
	for(j=0;j<=i;j++)
	    for(k=0;k<b;k++)
		for(l=0;l<b;l++)
		    if( fabs(x[i+n*(j+n*(k+b*l))]-x[j+n*(i+n*(l+b*k))])>tol ) {
			//printBlockMat("What is the error?","%1.4lf",n,b,x);
			Rprintf("\n%lf %lf\n",x[i+n*(j+n*(k+b*l))],x[j+n*(i+n*(l+b*k))]);
			error("checkBlockMatSymmetry: Not symmetric. %d %d %d %d %d %d",i,j,k,l,n,b);
		    }
}

void checkMatSymmetry(int n,double *a,int lda,double tol) {
    int i,j;
    for(i=0;i<n;i++)
	for(j=0;j<i;j++)
	    if( fabs(a[i+lda*j]-a[j+lda*i])>tol)
		error("checkMatSymmetry: Not symmetric. %d %d %d %d",i,j,n,lda);
}




extern void gsiCGSkrigingPrep(
    int *zDim,          // 
    double *z,          // nxD Data
    int *FDim,          // nxp
    double *F,          // n * p: Trend regressors of observed 
    int *GammaDim,      
    double *Gamma,      // (n*n)xDxD  Variogram of observed
    double *W,          // ldN
    int *lmaxP,         // 1: output
    double *Nmat,       // ldN
    int    *ldNp,        // 
    int    *ipiv,
    int    *type,
    int    *ref,
    int    *nmv,
    int    *tki
    ) {
    const int D=zDim[1];// D of comp
    const int D2=D*D;   // size of a error matrix
    const int n=zDim[0];// number of observations
    const int FD=FDim[1];//number of covariables+1
    const int ldN=*ldNp;//MaxFramesize of central matrix 
    int lmax=0;               // Actuall Framesize of central matrix 
    const char transA='N';    // for BLAS
    const char transB='T';    // for BLAS
    const int eins=1;         // for BLAS
    const double zero=0.0;    // for BLAS
    const double one =1.0;    // for BLAS
    int i,j;                  // cases
    int k,l;                  // actual matrix position
    int p,q;                  // variables in cases         
    //const char uplo='U';     // for LAPACK
    int lwork=ldN*60;        // for LAPACK (dsytrf)
    double *work=(double*) R_alloc(lwork,sizeof(double));// for LAPACK
    int info=0;              // for LAPACK
    int    mt;              // temporary for missingType
    int iIdx,jIdx,iRef,jRef;// the variables to actually use
    double v;               // current value in missingness Typing
    double tmp;
    if( (D-1)*(n+FD)>ldN ) error("gsiCGSkrigingPrep: ldN might be to small");
    if( GammaDim[0]!=n*n ) error("gsiCGSkrigingPrep: GammaDim[0] wrong");
    if( GammaDim[1]!=D ) error("gsiCGSkrigingPrep: GammaDim[1] wrong");
    if( GammaDim[2]!=D ) error("gsiCGSkrigingPrep: GammaDim[2] wrong");
    if( FDim[0]!=n ) error("gsiCGSkrigingPrep: FDim[0] wrong");
    checkBlockMatSymmetry(n,D,Gamma,1E-10);
    R_CheckUserInterrupt();
    // Create type
    for(i=0;i<n;i++)
	for(p=0;p<D;p++) {
	    v=z[i+n*p];
	    if( R_FINITE(v) ) {
		if( v > 0 ) { //NMV
		    mt=0;
		} else if( v < 0 ) { // BDL
		    mt=1;
		} else  { // BDL
		    mt=1;
		}
	    } else if( ISNA(v) ) { // MNAR
		mt=5;
	    } else if( v==R_PosInf ) { // ERR oder MAR (later)
		mt=4;
	    } else if( v==R_NegInf ) { // SZ
		mt=3;
	    }else if( ISNAN(v) ) { // MCAR (MAR)
		mt=2;
	    }
	    type[i+n*p]=mt;
	}
    // Initialize 
    for(i=0;i<n;i++) {
	for(p=0;p<D;p++) {
	    tki[i+n*p]=0;
	}
    }
    // Decide on alr
    for(i=0;i<n;i++) {
	ref[i]=0;
	nmv[i]=0;
	for(p=0;p<D;p++) {
	    if( type[i+n*p] == 0 ) {
		ref[i]=p;
		tki[i+n*(nmv[i]++)]=p;
	    } 
	}
    }
    // Framesize of Matrix
    lmax=(D-1)*FD;
    for(i=0;i<n;i++)
	if( nmv[i] > 1 )
	    lmax+=nmv[i]-1;
    // Build central Matrix
    // Gamma
    R_CheckUserInterrupt();
    k=0;
    for(i=0;i<n;i++) {
	iRef = ref[i];
	for(p=0;p<nmv[i]-1;p++) {
	    l=0;
	    iIdx=tki[i+n*p];
	    for(j=0;j<n;j++) {
		jRef = ref[j];
		for(q=0;q<nmv[j]-1;q++) {
		    jIdx = tki[j+n*q];
		    int lg=GammaDim[0]*GammaDim[1]*GammaDim[2];
//		    Rprintf("%d %d %d %d %d %d %d %d\n",i,n,j,D,iIdx,D,jIdx,lg);
		    if( i<0 || n < 0 || j <0 || iIdx < 0 || D<0 || jIdx <0 
			|| iRef<0 || jRef<0 )
			error("gsiCGSkriging: Out ouf bounds 1");
		    if( i+n*(j+D*(iIdx+D*jIdx))>=lg ) {
			Rprintf("%d %d %d %d %d %d %d %d\n",i,n,j,D,iIdx,D,jIdx,lg);
			error("gsiCGSkriging: Out ouf bounds 2");
		    }
		    if(	i+n*(j+D*(iRef+D*jRef))>=lg ) 
			error("gsiCGSkriging: Out ouf bounds 3");
		    if( i+n*(j+D*(iRef+D*jIdx))>=lg )
			error("gsiCGSkriging: Out ouf bounds 4");
		    if( i+n*(j+D*(iIdx+D*jRef))>=lg )
			error("gsiCGSkriging: Out ouf bounds 5");
		    tmp=Gamma[i+n*(j+n*(iIdx+D*jIdx))]+
			Gamma[i+n*(j+n*(iRef+D*jRef))]-
			Gamma[i+n*(j+n*(iRef+D*jIdx))]-
			Gamma[i+n*(j+n*(iIdx+D*jRef))];
		    if( k<0 || ldN<0 || l<0 ||k+ldN*l>=ldN*ldN)
			error("gsiCGSkriging: Out ouf bounds 6");
		    Nmat[k+ldN*l] = tmp;
//			Gamma[i+n*(j+n*(iIdx+D*jIdx))]+
//			Gamma[i+n*(j+n*(iRef+D*jRef))]-
//			Gamma[i+n*(j+n*(iRef+D*jIdx))]-
//			Gamma[i+n*(j+n*(iIdx+D*jRef))];
		    if( ++l > lmax )
			error("gsiCGSkriging: Error in framesize calculation l1.");
		}
	    }
	    // F 
	    //Rprintf("%lf\n",tmp);
	    for(j=0;j<FD;j++) {
		jRef = D-1;
		for(q=0;q<D-1;q++) {
		    jIdx=q;
		    if( k<0 || ldN<0 || l<0 ||k+ldN*l>=ldN*ldN)
			error("gsiCGSkriging: Out ouf bounds 7");
		    Nmat[k+ldN*l] = F[i+n*j]*(
			(iIdx==jIdx?1:0)+
			(iRef==jRef?1:0)-
			(iIdx==jRef?1:0)-
			(iRef==jIdx?1:0));
		    if( ++l > lmax )
			error("gsiCGSkriging: Error in framesize calculation l2.");
		}
	    }
	    // next line
	    if( ++k > lmax )
		error("gsiCGSkriging: Error in framesize calculation k1.");
	}
    }
     // F^t
    R_CheckUserInterrupt();
    for(i=0;i<FD;i++) {
	iRef = D-1;
	for(p=0;p<D-1;p++) {
	    l=0;
	    iIdx=p;
	    for(j=0;j<n;j++) {
		jRef = ref[j];
		for(q=0;q<nmv[j]-1;q++) {
		    jIdx = tki[j+n*q];
		    Nmat[k+ldN*l] = F[j+n*i]*(
			(iIdx==jIdx?1:0)+
			(iRef==jRef?1:0)-
			(iIdx==jRef?1:0)-
			(iRef==jIdx?1:0));
		    if( ++l > lmax )
			error("gsiCGSkriging: Error in framesize calculation l3.");
		}
	    }
	    // 0
	    for(j=0;j<FD;j++) {
		jRef = D-1;
		for(q=0;q<D-1;q++) {
		    jIdx=q;
		    Nmat[k+ldN*l] = 0;
		    if( ++l > lmax )
			error("gsiCGSkriging: Error in framesize calculation l4.");
		}
	    }
	    if( ++k > lmax )
		error("gsiCGSkriging: Error in framesize calculation k2.");
	}
    }
    checkMatSymmetry(lmax,Nmat,ldN,1E-10);
    R_CheckUserInterrupt();
    // Copy a column 
    //for(k=0;k<lmax;k++)
    //	W[k]=Nmat[k+ldN*4];
    //printMat("Main Mat","%2.4lf",lmax,lmax,Nmat,ldN);
    // Factorize (Multiple right hand sides and errors!!!)
    F77_CALL(dsytrf)("U", //const char* uplo, 
		     &lmax, //const int* n,
		     Nmat,  //double* a, 
		     &ldN,  //const int* lda, 
		     ipiv,  //int* ipiv,
		     work,  //double* work, 
		     &lwork,//const int* lwork, 
		     &info //int* info
		     FCONE
	);
    if( info > 0 )
	error("gsiCGSkriging: Kriging Matrix is singular");
    R_CheckUserInterrupt();
    // Covariogram Weights
    //  Prepare RHS of dual kriging
    k=0;
    for(i=0;i<n;i++) {
	double refLog=log(z[i+n*ref[i]]);
	for(p=0;p<nmv[i]-1;p++) {
	    iIdx=tki[i+n*p];
	    W[k++]=log(z[i+n*iIdx])-refLog;
	}
    }
    for(;k<lmax;)
	W[k++]=0;
    //printMat("Wpre","%2.4lf",lmax,1,W,ldN);
    // Compute CgramWeights
    F77_CALL(dsytrs)(
		     "U",	     //&uplo, //const char* uplo,
	&lmax, //const int* n,
	&eins, //const int* nrhs,
	Nmat,  //const double* a, 
	&ldN,  //const int* lda,
	ipiv,  // const int* ipiv,
	W,     // double* b, 
	&ldN,  //const int* ldb, 
	&info  // int* info
	FCONE
	);
    //printMat("Wpost","%2.4lf",lmax,1,W,ldN);
    *lmaxP=lmax;
}


extern void gsiCGSkrigingPredict(
    int *fDim,             
    double *f,          // np x p: Trend of predicted 
    int *gammaDim,      
    double *gamma,      // (n*np)xDxD: Variogram observed predicted
    int *predDim,
    double *pred,       // np x D: output predicted compositions
    double *W,
    int *lmaxP,         // 1: output
    int *nParam,       // ldN
    int    *ldNp,        // 
    int    *ipiv,
    int    *type,
    int    *ref,
    int    *nmv,
    int    *tki
    ) {
    const int D=predDim[1];// D of comp
    const int D2=D*D;   // size of a error matrix
    const int FD=fDim[1];
    const int n=*nParam; // number of observations
    const int ldN=*ldNp; //MaxFramesize of central matrix 
    int lmax=*lmaxP;    // Actuall Framesize of central matrix 
    const int np=predDim[0];//number of cases to predict
    const char transA='N';    // for BLAS
    const char transB='T';    // for BLAS
    const int eins=1;         // for BLAS
    const double zero=0.0;    // for BLAS
    const double one =1.0;    // for BLAS
    int i,j;                  // cases
    int k,l;                  // actual matrix position
    int p,q;                  // variables in cases         
    int iIdx,jIdx,iRef,jRef;// the variables to actually use
    double *RHS=(double*) R_alloc(lmax*D,sizeof(double));    // right vector
    double tmp;
    for(j=0;j<np;j++) {
	jRef = D-1;
	// gamma
	k=0;
	for(i=0;i<n;i++) {
	    iRef = ref[i];
	    for(p=0;p<nmv[i]-1;p++) {
		l=0;
		iIdx=tki[i+n*p];
		for(q=0;q<D-1;q++) {
		    jIdx = q;
		    RHS[k+lmax*q] = 
			gamma[i+n*(j+np*(iIdx+D*jIdx))]+
			gamma[i+n*(j+np*(iRef+D*jRef))]-
			gamma[i+n*(j+np*(iRef+D*jIdx))]-
			gamma[i+n*(j+np*(iIdx+D*jRef))];
		    l++;
		}
		RHS[k+lmax*q] = 0; // alr Reference
		k++;
	    }
	}
	// f
	for(i=0;i<FD;i++) {
	    iRef = D-1;
	    for(p=0;p<D-1;p++) {
		l=0;
		iIdx=p;
		for(q=0;q<D-1;q++) {
		    jIdx = q;
		    RHS[k+lmax*q] = f[j+np*i]*(			
			(iIdx==jIdx?1:0)+
			(iRef==jRef?1:0)-
			(iIdx==jRef?1:0)-
			(iRef==jIdx?1:0));
		    l++;
		}
		RHS[k+lmax*q] = 0; // alr Reference
		k++;
	    }
	}
	// Krige
	//printMat("RHS","%2.4lf",lmax,D,RHS,lmax);
	F77_CALL(dgemm)( // stimmt die Zuordnung?
			"T",		//&transB,  //const char *transa, 
			"N",		//&transA,  //const char *transb, 
	    &eins,    //const int *m,
	    &D,       //const int *n, 
	    &lmax,    //const int *k, 
	    &one,     //const double *alpha,
	    W,        //const double *a, 
	    &lmax,     //const int *lda,
	    RHS,      //const double *b, 
	    &lmax,     //const int *ldb,
	    &zero,    //const double *beta, 
	    pred+j,   //double *c, 
	    &np       //const int *ldc
			FCONE FCONE
	    );
	// Transform into composition
	tmp=0;
	for(q=0;q<D;q++)
	    tmp += (pred[j+np*q]=exp(pred[j+np*q]));
	for(q=0;q<D;q++)
	    pred[j+np*q]/=tmp;
    }    
}


// Kriging
static R_NativePrimitiveArgType gsiCGSkriging_t[] = {
  INTSXP,REALSXP,INTSXP,REALSXP,INTSXP,REALSXP,INTSXP,REALSXP,INTSXP,REALSXP,INTSXP,REALSXP,REALSXP,INTSXP
};



extern void gsiCGSkriging(
    int *zDim,
    double *z,          // n * D: observed compositions
    int *FDim,         
    double *F,          // n * p: Trend regressors of observed 
    int *GammaDim,
    double *Gamma,      // (n*n)xDxD  Variogram of observed
    int *fDim,             
    double *f,          // np x p: Trend of predicted 
    int *gammaDim,      
    double *gamma,      // (n*np)xDxD: Variogram observed predicted
    int *predDim,
    double *pred,       // np x D: output predicted compositions
    double *err,        // D x D x np or 0: output kriging errors 
    int *computeErrors  //boolean -> should err be computed
    ) {
    const int D=zDim[1];// D of comp
    const int D2=D*D;   // size of a error matrix
    const int n=zDim[0];// number of observations
    const int FD=FDim[1];//number of covariables+1
    const int ldN=(D-1)*(n+FD);//MaxFramesize of central matrix 
    int lmax=0;               // Actuall Framesize of central matrix 
    const int np=predDim[0];//number of cases to predict
    const char transA='N';    // for BLAS
    const char transB='T';    // for BLAS
    const int eins=1;         // for BLAS
    const double zero=0.0;    // for BLAS
    const double one =1.0;    // for BLAS
    int i,j;                  // cases
    int k,l;                  // actual matrix position
    int p,q;                  // variables in cases         
    double *Nmat=(double*) R_alloc(ldN*ldN,sizeof(double)); // central mat
    double *RHS=(double*) R_alloc(ldN*D,sizeof(double));    // right vector
    double *RHS2=(double*) R_alloc(ldN*D,sizeof(double));   // left vector err
    double *W  =(double*) R_alloc(ldN,sizeof(double));      // left vector pred
    //double *Ninv=(double*) R_alloc(ldN*ldN,sizeof(double));
    const char uplo='U';     // for LAPACK
    int lwork=ldN*60;        // for LAPACK (dsytrf)
    double *work=(double*) R_alloc(lwork,sizeof(double));// for LAPACK
    int info=0;              // for LAPACK
    int    *ipiv=(int *)   R_alloc(ldN,sizeof(int));// for LAPACK (dsytrf) 
    int    *type=(int *)   R_alloc(n*D,sizeof(int));//Missingnesstype
    int    *ref =(int *)   R_alloc(n,sizeof(int));  // denom in alr
    int    *nmv =(int *)   R_alloc(n,sizeof(int));  // # NonMissingValues
    int    *tki =(int *)   R_alloc(n*D,sizeof(int));// idx NonMissingValues
    double tmp;             // for sums
    int    mt;              // temporary for missingType
    int iIdx,jIdx,iRef,jRef;// the variables to actually use
    double v;               // current value in missingness Typing
    // Check 
    //printBlockMat("What is the error?","%1.4lf",n,D,Gamma);
    if( GammaDim[0]!=n*n ) error("gsiCGSkriging: GammaDim[0] wrong");
    if( GammaDim[1]!=D ) error("gsiCGSkriging: GammaDim[1] wrong");
    if( GammaDim[2]!=D ) error("gsiCGSkriging: GammaDim[2] wrong");
    if( gammaDim[0]!=n*np ) error("gsiCGSkriging: gammaDim[0] wrong");
    if( gammaDim[1]!=D ) error("gsiCGSkriging: gammaDim[1] wrong");
    if( gammaDim[2]!=D ) error("gsiCGSkriging: gammaDim[2] wrong");
    if( predDim[1]!=D ) error("gsiCGSkriging: gammaDim[2] wrong");
    if( FDim[0]!=n ) error("gsiCGSkriging: FDim[0] wrong");
    if( fDim[0]!=np ) error("gsiCGSkriging: fDim[0] wrong");
    if( fDim[1]!=FD ) error("gsiCGSkriging: fDim[1] wrong");
    checkBlockMatSymmetry(n,D,Gamma,1E-10);
    R_CheckUserInterrupt();
    // Create type
    for(i=0;i<n;i++)
	for(p=0;p<D;p++) {
	    v=z[i+n*p];
	    if( R_FINITE(v) ) {
		if( v > 0 ) { //NMV
		    mt=0;
		} else if( v < 0 ) { // BDL
		    mt=1;
		} else  { // BDL
		    mt=1;
		}
	    } else if( ISNA(v) ) { // MNAR
		mt=5;
	    } else if( v==R_PosInf ) { // ERR oder MAR (later)
		mt=4;
	    } else if( v==R_NegInf ) { // SZ
		mt=3;
	    }else if( ISNAN(v) ) { // MCAR (MAR)
		mt=2;
	    } // end if (*5)
	    type[i+n*p]=mt;
	} // end loop (i,p)
    // Initialize 
    for(i=0;i<n;i++) {
	for(p=0;p<D;p++) {
	    tki[i+n*p]=0;
	} // end loop p
    }// end loop i 
    // Decide on alr
    for(i=0;i<n;i++) {
	ref[i]=0;
	nmv[i]=0;
	for(p=0;p<D;p++) {
	    if( type[i+n*p] == 0 ) {
		ref[i]=p;
		tki[i+n*(nmv[i]++)]=p;
	    }  // end if type[i+np]
	} // end loop p: data variables 
    } // end loop i: data locations
    // Framesize of Matrix
    lmax=(D-1)*FD;
    for(i=0;i<n;i++)
	if( nmv[i] > 1 )
	    lmax+=nmv[i]-1;
    // Build central Matrix
    // Gamma
    R_CheckUserInterrupt();
    k=0;
    for(i=0;i<n;i++) {
	iRef = ref[i];
	for(p=0;p<nmv[i]-1;p++) {
	    l=0;
	    iIdx=tki[i+n*p];
	    for(j=0;j<n;j++) {
		jRef = ref[j];
		for(q=0;q<nmv[j]-1;q++) {
		    jIdx = tki[j+n*q];
		    int lg=GammaDim[0]*GammaDim[1]*GammaDim[2];
//		    Rprintf("%d %d %d %d %d %d %d %d\n",i,n,j,D,iIdx,D,jIdx,lg);
		    if( i<0 || n < 0 || j <0 || iIdx < 0 || D<0 || jIdx <0 
			|| iRef<0 || jRef<0 )
			error("gsiCGSkriging: Out ouf bounds 1");
		    if( i+n*(j+D*(iIdx+D*jIdx))>=lg ) {
			Rprintf("%d %d %d %d %d %d %d %d\n",i,n,j,D,iIdx,D,jIdx,lg);
			error("gsiCGSkriging: Out ouf bounds 2");
		    } //endif
		    if(	i+n*(j+D*(iRef+D*jRef))>=lg ) 
			error("gsiCGSkriging: Out ouf bounds 3");
		    if( i+n*(j+D*(iRef+D*jIdx))>=lg )
			error("gsiCGSkriging: Out ouf bounds 4");
		    if( i+n*(j+D*(iIdx+D*jRef))>=lg )
			error("gsiCGSkriging: Out ouf bounds 5");
		    tmp=Gamma[i+n*(j+n*(iIdx+D*jIdx))]+
			Gamma[i+n*(j+n*(iRef+D*jRef))]-
			Gamma[i+n*(j+n*(iRef+D*jIdx))]-
			Gamma[i+n*(j+n*(iIdx+D*jRef))];
		    if( k<0 || ldN<0 || l<0 ||k+ldN*l>=ldN*ldN)
			error("gsiCGSkriging: Out ouf bounds 6");
		    Nmat[k+ldN*l] = tmp;
//			Gamma[i+n*(j+n*(iIdx+D*jIdx))]+
//			Gamma[i+n*(j+n*(iRef+D*jRef))]-
//			Gamma[i+n*(j+n*(iRef+D*jIdx))]-
//			Gamma[i+n*(j+n*(iIdx+D*jRef))];
		    if( ++l > lmax )
			error("gsiCGSkriging: Error in framesize calculation l1.");
		} // end loop q: data2 variables non-missing
	    } // end loop j: data2 locations
	    // F 
	    //Rprintf("%lf\n",tmp);
	    for(j=0;j<FD;j++) {
		jRef = D-1;
		for(q=0;q<D-1;q++) {
		    jIdx=q;
		    if( k<0 || ldN<0 || l<0 ||k+ldN*l>=ldN*ldN)
			error("gsiCGSkriging: Out ouf bounds 7");
		    Nmat[k+ldN*l] = F[i+n*j]*(
			(iIdx==jIdx?1:0)+
			(iRef==jRef?1:0)-
			(iIdx==jRef?1:0)-
			(iRef==jIdx?1:0));
		    if( ++l > lmax )
			error("gsiCGSkriging: Error in framesize calculation l2.");
		} // end loop q (variables in predictions)
	    } // end loop j (covariables) 
	    // next line
	    if( ++k > lmax )
		error("gsiCGSkriging: Error in framesize calculation k1.");
	} // end loop p: data variables (not missing)
    } // end loop i data locations
     // F^t
    R_CheckUserInterrupt();
    for(i=0;i<FD;i++) {
	iRef = D-1;
	for(p=0;p<D-1;p++) {
	    l=0;
	    iIdx=p;
	    for(j=0;j<n;j++) {
		jRef = ref[j];
		for(q=0;q<nmv[j]-1;q++) {
		    jIdx = tki[j+n*q];
		    Nmat[k+ldN*l] = F[j+n*i]*(
			(iIdx==jIdx?1:0)+
			(iRef==jRef?1:0)-
			(iIdx==jRef?1:0)-
			(iRef==jIdx?1:0));
		    if( ++l > lmax )
			error("gsiCGSkriging: Error in framesize calculation l3.");
		} // end loop q (data2 non-missing observations)
	    } // end loop j (data locations)
	    // 0
	    for(j=0;j<FD;j++) {
		jRef = D-1;
		for(q=0;q<D-1;q++) {
		    jIdx=q;
		    Nmat[k+ldN*l] = 0;
		    if( ++l > lmax )
			error("gsiCGSkriging: Error in framesize calculation l4.");
		} //end loop q (data variables)
	    } // end loop j (covariables)
	    if( ++k > lmax )
		error("gsiCGSkriging: Error in framesize calculation k2.");
	} // endloop p (variables in predictions)
    } // end loop i (covariables)
    checkMatSymmetry(lmax,Nmat,ldN,1E-10);
    R_CheckUserInterrupt();
    // Copy a column 
    //for(k=0;k<lmax;k++)
    //	W[k]=Nmat[k+ldN*4];
    //printMat("Main Mat","%2.4lf",lmax,lmax,Nmat,ldN);
    // Factorize (Multiple right hand sides and errors!!!)
    F77_CALL(dsytrf)("U",//&uplo, //const char* uplo, 
		     &lmax, //const int* n,
		     Nmat,  //double* a, 
		     &ldN,  //const int* lda, 
		     ipiv,  //int* ipiv,
		     work,  //double* work, 
		     &lwork,//const int* lwork, 
		     &info //int* info
		     FCONE
	);
    if( info > 0 )
	error("gsiCGSkriging: Kriging Matrix is singular");
    R_CheckUserInterrupt();
    // Covariogram Weights
    //  Prepare RHS of dual kriging
    k=0;
    for(i=0;i<n;i++) {
	double refLog=log(z[i+n*ref[i]]);
	for(p=0;p<nmv[i]-1;p++) {
	    iIdx=tki[i+n*p];
	    W[k++]=log(z[i+n*iIdx])-refLog;
	} // end loop p (non-missing data locations)
    } // end loop i (data locations)
    for(;k<lmax;)
	W[k++]=0;
    //printMat("Wpre","%2.4lf",lmax,1,W,ldN);
    // Compute CgramWeights
    F77_CALL(dsytrs)(
		     "U",	     //&uplo, //const char* uplo,
	&lmax, //const int* n,
	&eins, //const int* nrhs,
	Nmat,  //const double* a, 
	&ldN,  //const int* lda,
	ipiv,  // const int* ipiv,
	W,     // double* b, 
	&ldN,  //const int* ldb, 
	&info  // int* info
	FCONE
	);
    //printMat("Wpost","%2.4lf",lmax,1,W,ldN);
    // Kriging loop
    R_CheckUserInterrupt();
    for(j=0;j<np;j++) {
	jRef = D-1;
	// gamma
	k=0;
	for(i=0;i<n;i++) {
	    iRef = ref[i];
	    for(p=0;p<nmv[i]-1;p++) {
		l=0;
		iIdx=tki[i+n*p];
		for(q=0;q<D-1;q++) {
		    jIdx = q;
		    RHS[k+ldN*q] = 
			gamma[i+n*(j+np*(iIdx+D*jIdx))]+
			gamma[i+n*(j+np*(iRef+D*jRef))]-
			gamma[i+n*(j+np*(iRef+D*jIdx))]-
			gamma[i+n*(j+np*(iIdx+D*jRef))];
		    l++;
		}
		RHS[k+ldN*q] = 0; // alr Reference
		k++;
	    }
	}
	// f
	for(i=0;i<FD;i++) {
	    iRef = D-1;
	    for(p=0;p<D-1;p++) {
		l=0;
		iIdx=p;
		for(q=0;q<D-1;q++) {
		    jIdx = q;
		    RHS[k+ldN*q] = f[j+np*i]*(			
			(iIdx==jIdx?1:0)+
			(iRef==jRef?1:0)-
			(iIdx==jRef?1:0)-
			(iRef==jIdx?1:0));
		    l++;
		}
		RHS[k+ldN*q] = 0; // alr Reference
		k++;
	    }
	}
	// Krige
	//printMat("RHS","%2.4lf",lmax,D,RHS,ldN);
	F77_CALL(dgemm)(
			"T", //&transB,  //const char *transa, 
			"N", //&transA,  //const char *transb, 
	    &eins,    //const int *m,
	    &D,       //const int *n, 
	    &lmax,    //const int *k, 
	    &one,     //const double *alpha,
	    W,        //const double *a, 
	    &ldN,     //const int *lda,
	    RHS,      //const double *b, 
	    &ldN,     //const int *ldb,
	    &zero,    //const double *beta, 
	    pred+j,   //double *c, 
	    &np       //const int *ldc
	    FCONE FCONE
	    );
	// Transform into composition
	tmp=0;
	for(q=0;q<D;q++)
	    tmp += (pred[j+np*q]=exp(pred[j+np*q]));
	for(q=0;q<D;q++)
	    pred[j+np*q]/=tmp;
	// Errors
	if( *computeErrors ) {
	    //Rprintf("Calculating Errors\n");
	    // copy RHS
	    R_CheckUserInterrupt();
	    for(i=0;i<lmax;i++)
		for(q=0;q<D;q++)
		    RHS2[i+ldN*j]=RHS[i+ldN*j];
	    // Solve
	    F77_CALL(dsytrs)(
			     "U",//&uplo, //const char* uplo,
		&lmax, //const int* n,
		&D, //const int* nrhs,
		Nmat,  //const double* a, 
		&ldN,  //const int* lda,
		ipiv,  // const int* ipiv,
		RHS2,  // double* b, 
		&ldN,  //const int* ldb, 
		&info  // int* info
		FCONE
		);
	    // Multiply
	    F77_CALL(dgemm)(
			    "T", //&transB,  //const char *transa, 
			    "N", //&transA,  //const char *transb, 
		&D,       //const int *m,
		&D,       //const int *n, 
		&lmax,    //const int *k, 
		&one,     // const double *alpha,
		RHS,      //const double *a, 
		&ldN,     //const int *lda,
		RHS2,     //const double *b, 
		&ldN,     //const int *ldb,
		&zero,    //const double *beta, 
		err+D2*j, //double *c, 
		&D       //const int *ldc
		FCONE FCONE	    
		);
	    // alr to clr errors
	    tmp=0.0;
	    for(p=0;p<D;p++)
		for(q=0;q<D;q++)
		    tmp+=err[D2*j+p+D*q];
	    tmp/=D2;
	    for(p=0;p<D;p++)
		for(q=0;q<D;q++)
		    err[D2*j+p+D*q]+=tmp;
	}
    } // end j-loop (predicted locations)
    // End Kriging
}
    





/* Linear Model imputation */

#define Mcomp(i,j) comp[s[i]+N*(j)]
#define Mmt(i,j) mt[s[i]+N*(j)]
#define Mmtn(i,j) mt[snew[i]+N*(j)]
#define Mdl(i,j) dl[s[i]+N*(j)]
#define Mclr(i,j) clr[s[i]+N*(j)]
#define Mxlr(i,j) xlr[s[i]+N*(j)]
#define MalrDL(i,j) alrDL[s[i]+N*(j)]
#define MxlrImp(i,j) xlrImp[s[i]+N*(j)]
#define Mpred(i,j) pred[s[i]+N*(j)]
#define MX(i,j) X[s[i]+N*(j)]
#define Mbeta(i,j) beta[(i)+p*(j)]
#define MclrV(i,j) clrV[(i)+D*(j)]
#define MXXinv(i,j) XXinv[(i)+p*(j)]
#define Mnorm(i,j) norm[(i)+dimNorm[0]*(j)]
#define VmtType(i) mtType[s[i]]
#define Vused(i) used[s[i]]
#define Vs(i) s[i]
#define Vsnew(i) snew[i]
#define MimputationIdx(i,j) imputationIdx[(i)+nMissingTypes*(j)]
#define VimputationN(i) imputationN[i]
#define MhpA(i,j) hpA[(i)+D*(j)]
#define TXbig(i,j,c,l) (Xbig[(c)+n*((i)+D*((l)+p*(j)))])
/* gsiCImpAcompGetTypes

This function creates and assigns the missingTypes for the ImpAcomp subsystem.
A missingType is defined as a combination of missing and nonmissing values,
regardless of the type of missingness. 


*/
extern void gsiCImpAcompGetTypes(
    int *Dp,      // 1: Number of dimensions
    int *np,      // 1: Number of cases
    int *Np,      // 1: leading dimension of mt
    int *s,       // n: Cases to be used
    int *mt,      // N x D a missingness signature
    int *snew,    // Out: np or nMissingTypes: the lines of mt for each situation
    int *mtType,  // Out: N: The missingType of each case
    int *nMissingTypes // Out 1: 
    ) {
    const int D=*Dp; // columns of dataset
    const int n=*np; // number cases used
    const int N=*Np; // rows of dataset
    int i,j,k; // cases,columns,missingTypes
    int mtTypes=0; // number of missingTypes
    int fits=0; // has a existing Missingtype been found
    // Loop over cases
    for(i=0;i<n;i++) {
	// Loop over already found missingTypes
	for(k=0;k<mtTypes;k++) {
	    fits = 1;
	    // Check each component
	    for(j=0;j<D;j++) {
		if( (Mmt(i,j)?0:1)+(Mmtn(k,j)?0:1) == 1 ) {
		    fits=0; break;
		}
	    }
	    if( ! fits ) continue;
	    // Assign missingtype and go to next case
	    VmtType(i)=k;
	    break;
	}
	if( ! fits ) // if nothing is found, create new missingType 
	    Vsnew(mtTypes++)=Vs(i);
    }
    // return number of missingTypes
    *nMissingTypes=mtTypes;
}

/* gsiCImpAcompGetIdx 

   computes a reordering of indices for each missingType for the ImpAcomp
   subsystem.  The missing colums come first.

*/

extern void gsiCImpAcompGetIdx(
    int *Dp,      // 1: Number of dimensions
    int *np,      // 1: Number of missingTypes
    int *Np,      // 1: leading dimension of mt
    int *s,       // n: Cases to be used for each missingType
    int *mt,      // NxD a missingness signature
    int    *imputationIdx,   // Out: nxD: new sequence of columns for each missingType 
    int    *imputationN      // Out: n: Number of missings
    ) {
    const int D=*Dp;
    const int n=*np;
    const int N=*Np;
    const int nMissingTypes=n;
    int i,j,k;
    for(i=0;i<n;i++) {
	imputationN[i]=0;
	k=0;
	for(j=0;j<D;j++)
	    if(Mmt(i,j)) {
		VimputationN(i)++;
		MimputationIdx(i,k++)=j;
	    }
	for(j=0;j<D;j++)
	    if(!Mmt(i,j)) 
		MimputationIdx(i,k++)=j;
    }
}


/* gsiCImpAcompNewImputationVariance

   Key function of the ImpAcomp subsystem

   computes a the predictors and the residual variance for each missingType, 
   which allow to give the conditional expectation of the missing part of the
   alr from the observed part of the alr. 
   
   The alr is stored with a zero for the reference component. The reference is
   the last nonmissing component. In case of only missing components it is the
   last component.

   The predictors are stored in imputationCache according to the block matrix
   form
   resVar  lambda 
   Undef   Undef

*/


extern void gsiCImpAcompNewImputationVariance(
    int *Dp,       // Number of dimensions
    int *np,       // Number of different imputation situations
    double *var, // DxD the current estimate of prediction variance
    double *imputationCache, // Output: DxDxnx2 stores predictors and conditional variance
    int    *imputationIdx,   // Dxn indices of missings 
    int    *imputationN  ,   // Number of missings
    double *hpA,   // DxD helper
    double *work,  // large helper
    int *lwork,    // size of the helper
    int *iwork,    // helper for dgelsd
    int *liwork    // size of  iwork
    ) {
    const int n=*np; // number of missingTypes
    const int nMissingTypes=n; // number of missingTypes
    const int D=*Dp; // number of variables 
    const int DD=D*D;// step in imputationCache to get to the next matrix
    const int DDn=D*D*n; // step in imputationCache to get second type of info
    int   rank=0;    // parameter for dgelsd
    char transA='N'; // parameter for dgelsd
    int  info; // parameter for dgelsd
    double rcond=1E-4;// parameter for dgelsd
    double minusOne=-1.0; // parameter for dgemm 
    double one=1.0;       // parameter for dgemm
    double *resVar; // missings x missings residual variance (part of imputationCache)
    double *lambda; // missings x given    prediction weights (part of imputationCache)
    //double *lambdaT; // given x missings lambda^t (part of hpA)
    // Organisation of hpA as blockmatrix:
    // Vmm Vgm
    // Vmg Vgg
    // e.g. Vmm=Variance missing missing
    // The formulae computd in this routine
    // lambda = Vmg * VggI 
    // Vp = Vmm - Vmg * VggI * Vgm
    int i,j;  // index on permuted components, index on permuted components
    int k;    // index of missingType
    // loop over missingTypes
    for(k=0;k<n;k++) {
	int nMissing=imputationN[k];
	int nonMissingValues=D-nMissing;
	lambda = imputationCache+(DD*k+D*nMissing);
	resVar = imputationCache+(DD*k);
	// zero the rest of imputationCache
	for(i=nMissing;i<D;i++)
	    for(j=0;j<D;j++)
		imputationCache[DD*k+i+j*D]=0.0;
	if( nMissing == 0 ) {
	    // Nothing to be done (not prediction, no residuals)
	} else if( nonMissingValues < 2 ) {
	    // No not enough nonmissing values to observe anything
	    // Simply create alr variance in as residual variance
	    for(i=0;i<D;i++)
		for(j=0;j<D;j++)
		    imputationCache[DD*k+i+j*D]=
			var[MimputationIdx(k,i)+D*MimputationIdx(k,j)]+
			var[MimputationIdx(k,D-1)+D*MimputationIdx(k,D-1)]-
			var[MimputationIdx(k,D-1)+D*MimputationIdx(k,j)]-
			var[MimputationIdx(k,i)+D*MimputationIdx(k,D-1)]
			;
	} else {
	    // prepare alr variance for conditioning in alr
	    for(i=0;i<D;i++)
		for(j=0;j<D;j++)
		    MhpA(i,j)=
			var[MimputationIdx(k,i)+D*MimputationIdx(k,j)]+
			var[MimputationIdx(k,D-1)+D*MimputationIdx(k,D-1)]-
			var[MimputationIdx(k,D-1)+D*MimputationIdx(k,j)]-
			var[MimputationIdx(k,i)+D*MimputationIdx(k,D-1)]
			;
	    // Compute lambda = Vmg*VggI as 
	    // lambdaT= VggI * Vgm into Vgm
	    if( *liwork < 1 || *liwork < 3*D*23)
		error("liwork to small in gsiCNewImputationVariance ");
	    int Dnew=D;
	    F77_CALL(dgelsd)(
		&nonMissingValues,//const int* m, 
		&nonMissingValues,//const int* n, 
		&nMissing, //const int* nrhs,
		hpA+nMissing*(D+1),//double* a, 
		&Dnew,//const int* lda, 
		hpA+nMissing,//double* b, 
		&Dnew,//const int* ldb,
		lambda,//double *s, 
		&rcond,//double *rcond, 
		&rank,//int *rank, 
		work,//double *work, 
		lwork,//int *lwork,
		iwork,//int *iwork, 
		&info//int *info
		);
	    // copy result to lambda
	    for(i=0;i<nonMissingValues;i++)
		for(j=0;j<nMissing;j++)
		    lambda[j+i*D]=hpA[nMissing+i+D*j];
	    /* -> This would make it general (not alr specific)
                  However it is not possible to make the residual variance
		  unspecific.
	    // Correction of lambda on the input side
	    for(i=0;i<nonMissingvalues;i++)
		for(j=0;j<nMissing)
		    lambda[j+(nonMissingvalues-1)*D]-=lambda[j+i*D];
	    // Correction of lambda on the output side
	    for(j=0;j<nMissing)
		lambda[j+(nonMissingvalues-1)*D]+=1;
	    */
	    // Compute resVar = Vmm - Vmg *lambdaT into Vmm
	    F77_CALL(dgemm)(
			    "N", //&transA,//const char *transa, 
			    "N", //&transA,//const char *transb, 
		&nMissing,//const int *m,
		&nMissing,//const int *n, 
		&nonMissingValues,//const int *k, 
		&minusOne,//const double *alpha,
		hpA+nMissing*D,//const double *a, 
		&Dnew,//const int *lda,
		hpA+nMissing,//const double *b, 
		&Dnew,//const int *ldb,
		&one,//const double *beta, 
		hpA,//double *c, 
		&D//const int *ldc
		FCONE FCONE
		);
	    // copy result to resVar
	    for(i=0;i<nMissing;i++)
		for(j=0;j<nMissing;j++)
		    resVar[i+j*D]=hpA[i+D*j];
	}
    }
}

/* gsiCImpAcompFillAlr
  
   Computes the alr of the nonmissing part.

*/ 
extern void gsiCImpAcompFillAlr(
    int *Dp,            // 1: Dimension of the compostion
    int *np,            // 1: number of cases
    int *Np,            // 1: size of matrices
    int *s,             // the lines to process
    double *comp,       // N*D: The observed values 
    double *xlr,        // nxD in/out the xlr values 
    int *mtType,        // Out: N: The missingType of each case
    int    *imputationIdx,   // Dxn indices of missings 
    int    *imputationN,      // Number of missings
    int    *nMissingTypesP
    ) {
    const int nMissingTypes=*nMissingTypesP;
    const int D=*Dp;    // number of components
//    const int DD=D*D;   // first two dimensions of imputationCache
    const int n=*np;    // number of cases used
    const int N=*Np;    // number of cases in dataset
    int i,j;          // case,component 
    int type;         // the missingType
    int nMissings;    // number of missings in missingType
    for(i=0;i<n;i++) {
	type=VmtType(i);
	nMissings =VimputationN(type);
	if( nMissings < D ) {
	    double logdenominator=log(Mcomp(i,MimputationIdx(type,D)));
	    for(j=nMissings;j<D;j++) {
		int jj=MimputationIdx(type,j);
		Mxlr(i,jj)=log(Mcomp(i,jj))-logdenominator;
	    }
	}
    }
}

/* gsiCImpAcompCompleteAlr
   
   Computes the a imputationIdx[type,D-1] alr of the nonmissing part.

   and 
   
   And predicts with the missing parts with according to prediction and other
   parts.
*/

extern void gsiCImpAcompCompleteAlr(    
    int *Dp,            // 1: Dimension of the compostion
    int *np,            // 1: number of cases
    int *Np,            // 1: size of matrices
    int *s,             // the lines to process
    double *comp,       // N*D: The observed values 
    double *xlr,        // Out: nxD in/out the xlr values 
    double *pred,       // nxD prediction from model
    int *mtType,        // N: The missingType of each case
    double *imputationCache, // DxDxnx2 stores predictors and conditional variance
    int    *imputationIdx,   // Dxn indices of missings 
    int    *imputationN,      // Number of missings
    int    *nMissingTypesP,
    double *work,        // at least D, used for storing the reordered alr
    int *lwork           // size of work
    ) {
    const int nMissingTypes=*nMissingTypesP;
    const int D=*Dp;    //  number of components
    const int DD=D*D;   //  increment of last dimension in imputationCache
    const int n=*np;    //  number of cases to process
    const int N=*Np;    //  number of cases in dataset
    char transA='N';    //  for dgemv
    double one=1.0;     //  for dgemv
    int eins=1;         //  for dgemv
    int type;           //  missingType
    int nMissings;      //  number of missings
    int nonMissingValues;// D-nMissings
    int i,j;            //  case,component
    double logDenominator,logPredDenominator; // log of reference  in alr and prediction. 
    if( *lwork < D )     // Check lwork 
	error("work to small in gsiCImpAcompCompleteAlr");
    // loop over cases
    for(i=0;i<n;i++) {
	type=VmtType(i);              
	nMissings =VimputationN(type);
	nonMissingValues=D-nMissings;
	if( nMissings == D ) {  // completely missing, use prediction directly
	    logPredDenominator=log(Mcomp(i,MimputationIdx(type,D-1)));
	    for(j=0;j<nMissings;j++) {
		int jj=MimputationIdx(type,j);
		Mxlr(i,jj)=(Mpred(i,jj)-logPredDenominator);
	    }
	} else if( nMissings == 0 ) { // nothing missing, use observation directly
	    logDenominator=log(Mcomp(i,MimputationIdx(type,D-1)));
	    // alr
	    for(j=nMissings;j<D;j++) {
		int jj=MimputationIdx(type,j);
		Mxlr(i,jj)=log(Mcomp(i,jj))-logDenominator;
	    }
	} else if( nMissings < D ) { // something missing, impute type 1
	    logDenominator=log(Mcomp(i,MimputationIdx(type,D-1)));
	    logPredDenominator=log(Mcomp(i,MimputationIdx(type,D-1)));
 	    // the only predicted part of alr
	    for(j=0;j<nMissings;j++) {
		int jj=MimputationIdx(type,j);
		work[j]=(Mpred(i,jj)-logPredDenominator);
	    }
	    // alr of observed part and difference to prediction to work
	    for(j=nMissings;j<D;j++) {
		int jj=MimputationIdx(type,j);
		work[j]=(Mxlr(i,jj)=log(Mcomp(i,jj))-logDenominator)-
		    (Mpred(i,jj)-logPredDenominator);
	    }
 	    // correct the prediction from other observations
	    double *lambda=imputationCache+DD*type+D*nMissings;
	    // workm += lambda work2
	    F77_CALL(dgemv)(
			    "N",//&transA, //const char *trans, 
		&nMissings,//const int *m, 
		&nonMissingValues,//const int *n,
		&one,//const double *alpha, 
		lambda,//const double *a, 
		&D,//const int *lda,
		work+nMissings,//const double *x, 
		&eins,//const int *incx, 
		&one,//const double *beta,
		work,//double *y, 
		&eins//const int *incy
		FCONE
		);
	    // copy result to output
	    for(j=0;j<nMissings;j++) {
		int jj=MimputationIdx(type,j);
		Mxlr(i,jj)=work[j];
	    }
	}
    }
}

/* gsiCImpAcompAlrDetectionlimit 
   
   for every BDL component it computes the alr value the detectionlimit would have.


 */ 

extern void gsiCImpAcompAlrDetectionlimit(
    int *Dp,            // 1: Dimension of the compostion
    int *np,            // 1: number of cases
    int *Np,            // 1: size of matrices
    int *s,             // the lines to process
    double *comp,       // N*D: The observed values 
    int    *mt,         // N*D: 0=ok, 1=BDL, 2=MAR, 3-5=other (treated as MAR) 
    double *dl,         // N*D: Detectionlimit if BDL
    double *alrDL,       // output: N*D: alr of detectionlimit
    int *mtType,        // Out: N: The missingType of each case
    int    *imputationIdx,   // Dxn indices of missings 
    int    *imputationN      // Number of missings
    ) {
    const int D=*Dp;    //  number of components
    const int n=*np;    //  number of cases to process
    const int N=*Np;    //  number of cases in dataset
    int type;           //  missingType
    int nMissings;      //  number of missings
    int nonMissingValues;// D-nMissings
    int i,j;            //  case,component
    double logDenominator; // log of reference in alr and prediction. 
    // zero output
    for(i=0;i<n;i++)
	for(j=0;j<n;j++)
	    MalrDL(i,j)=0;
    // loop over cases
    for(i=0;i<n;i++) {
	type=VmtType(i);              
	nMissings =VimputationN(type);
	nonMissingValues=D-nMissings;
	if( nonMissingValues == 0 ) {
	    // There is no information at all 
	} else if(nMissings==0 ) {
	    // Nothing to be done
	} else {
	    logDenominator = log(Mcomp(i,D-1));
	    // loop over components
	    for(j=0;j<D;j++)
		if(Mmt(i,j)==1) // is a BDL 
		    MalrDL(i,j)=log(Mdl(i,j))-logDenominator;
	}
    }
}

/*
  gsiCImpAcompAddToXtY

  Adds an observation to a XtY matrix




 */
extern void gsiCImpAcompAddToXtY(
    const int p,
    const int incx,
    const double *X,
    const int D,
    const int incy,
    const double *Y,
    double *C,
    int *idx
    ) {
    int i,j;
    for(i=0;i<p;i++) 
	for(j=0;j<D;j++)
	    C[i+p*idx[j]]=X[i*incx]*Y[j*incy];
}
    

/* gsiCImpAcompClrExpectation

   Computes expectation (type 2) for BDL imputed value by approximating the
   integral by evaluation in given Montecarlo points and gives type 1
   imputation values.
   
   


*/
extern void gsiCImpAcompClrExpectation(
    int *Dp,            // 1: Dimension of the compostion
    int *np,            // 1: number of cases
    int *Np,            // 1: size of matrices
    int *s,             // n: the cases to be used
    int    *mt,         // N*D: 0=ok, 1=BDL, 2=MAR, 3-5=other (treated as MAR) 
    double *xlr,        // nxD  the type 1 MAR imputed alr values 
    double *alrDL,      // N*D: Detectionlimit if BDL as alr value
    double *xlrImp,     // Out: NxD  the type 1 MAR or BDL imputed alr values 
    double *used,       // Out: N: the number of MonteCarlo points fitting condition
    int    *mtType,     // Out: N: The missingType of each case
    double *imputationCache, //  DxDxnx2 stores predictors and conditional variance
    int    *imputationIdx,   // Dxn indices of missings 
    int    *imputationN,  // Number of missings
    int    *NmissingTypes, // number of different missing types
    double *clrE,       // out: D: computed expectation of clr
    double *clrVar,      // out: DxD: computed expectation of clr clr^t
    int *dimNorm,        // 2: the dimension of the normalized residuals 
    double *norm,        // normalized residuals    
    double *work,        // at least D*D+4*D
    double *lwork        // size of work
    ) {
    const int nMissingTypes=*NmissingTypes;
    const int n=*np;    // number of cases to be processesed 
    const int D=*Dp;     // columns in dataset
    const int DD=D*D;    // first two dimension of imputationCache
    const int DDn=D*D* *NmissingTypes; // first three dimension of imputationCache
// step in imputationCache to get second type of info
    const int lmax=dimNorm[0]; // number of cases in norm
    const int N=*Np;    // leading dimension of mt,xlr,alrDl,xlrImp
    char transA='N';    //  for dgemv
    char transAt='T';   //  for dgemv
    double one=1.0;     //  for dgemv
    int eins=1;         //  for dgemv
    double *chol=imputationCache+DDn; // the second half of imputationCache
    int type;           //  missingType, the current missingType
    int nMissings;      //  number of missings
    int nonMissingValues;// D-nMissings
    int i,j,k,l;           //  component,component,cases,normcases
    double *smallClrE=work;   // D: For local mean of alr
    double *smallClr2E=smallClrE+D;  // D*D: for local mean of alr alr^t
    double *smallAlr=smallClr2E+DD;  // D: the full reconstructed ALR reordered
    double *smallAlrDL=smallAlr+D;   // D: the detection limit reordered
    double *smallAlrPre=smallAlrDL+D;// D: the prediction reordered
    int fitters;       // the number of cases fitting the description
    int fits;          // true if the case in Norm fits description 
    int truei,truej;   // the true columns (if i and j contain reordered columns)
    double tmp;
    // check prerequirements work
    if( *lwork < D*D+4*D )
	error("gsiCImpAcompClrExpectation: workspace insufficient");
    // zero output
    for(i=0;i<D;i++)
	clrE[i]=0.0;
    for(i=0;i<D;i++)
	for(j=0;j<D;j++)
	    clrVar[i+D*j]=0.0;
    for(k=0;k<n;k++)
	Vused(k)=-1;
    // loop over cases
    for(k=0;k<n;k++) {
	type=VmtType(k);              
	nMissings =VimputationN(type);
	nonMissingValues=D-nMissings;
	if( nMissings==0) { // simply use
	    // no reordering
	    for(j=0;j<D;j++) {
		clrE[j]+=Mxlr(k,j);
		MxlrImp(k,j)=Mxlr(k,j);
	    }
	    for(i=0;i<D;i++)
		for(j=0;j<D;j++)
		    clrVar[i+D*j]+=Mxlr(k,i)*Mxlr(k,j);
	} else if( nonMissingValues == 0 ) { // use predicted distribution
	    // no reordering
	    for(j=0;j<D;j++) {
		clrE[j]+=Mxlr(k,j);
		MxlrImp(k,j)+=Mxlr(k,j);
	    }
	    // it is not reordered we can use the variance directly 
	    for(i=0;i<D;i++)
		for(j=0;j<D;j++)
		    clrVar[i+D*j]+=Mxlr(k,i)*Mxlr(k,j)+imputationCache[DD*type+i+D*j];
	} else {
	    // Are there any BDL?
	    fits=1;
	    for(i=0;fits&&i<nMissings;i++) {
		truei=MimputationIdx(k,i);
		switch(Mmt(k,truei)) {
		    case 0:error("gsiCImpAcompClrExpectation: internal error: missing not missing (1)");
		    case 1:fits=0;
			break;
		    case 2:
		    default:
			break;
		}
	    }
	    if( fits ) { // No, only MARs
		// Use prediction as mean
		for(j=0;j<D;j++)
		    clrE[j]+=Mxlr(k,j);
	    // we can use the variance directly (but reordered) 
		for(i=0;i<D;i++) {
		    truei=MimputationIdx(k,i);
		    for(j=0;j<D;j++) {
			truej=MimputationIdx(k,j);
			clrVar[truei+D*truej]+=Mxlr(k,i)*Mxlr(k,j)+
			    (i<nMissings&&j<nMissings)
			    ?imputationCache[DD*type+i+D*j]
			    :0;
		    }
		}
		continue;
	    }
	    // -> there are BDLs
	    // here comes the core algorithm for BDL
	    //  It computes conditional 1 and second noncentered moment 
	    //  on the subset of Norm that fits that after transformation
	    //  of mean and variance fits the description
	    // a)zero the small arrays
	    fitters=0;
	    for(i=0;i<D;i++)
		smallClrE[i]=0.0;
	    for(i=0;i<D;i++)
		for(j=0;j<D;j++)
		    smallClr2E[i+D*j]=0.0;
	    // b) Gather information
	    double *myL = chol+DD*type;
	    for(i=0;i<nMissings;i++)
		smallAlrDL[i]=MalrDL(k,MimputationIdx(k,i));
	    for(i=0;i<D;i++)
		smallAlrPre[i]=Mxlr(k,MimputationIdx(k,i));
	    for(i=nMissings;i<D;i++) 
		smallAlr[i]=Mxlr(k,MimputationIdx(k,i));
	    // c) Loop over example residuals from Norm
	    if( dimNorm[1] < nMissings )
		error("gsiCImpAcompClrExpectation: Not enough columns in Norm");
	    for(l=0;l<lmax;l++) {
		// complete alr
		//   copy mean
		for(i=0;i<nMissings;i++)
		    smallAlr[i]=smallAlrPre[i];
		// smallAlrm = mean+myL t(norm[l,])
		F77_CALL(dgemv)(
				"N", //&transA, // const char *trans, 
		    &nMissings,//const int *m, 
		    &nMissings,//const int *n,
		    &one, //const double *alpha, 
		    myL,//const double *a, 
		    &D, //const int *lda,
		    norm+l,//const double *x, 
		    &lmax,//const int *incx, 
		    &one,//const double *beta,
		    smallAlr,//double *y, 
		    &eins//const int *incy
		    FCONE
		    );
		// check conditions
		fits=1;
		for(i=0;fits&&i<nMissings;i++) {
		    truei=MimputationIdx(k,i);
		    switch(Mmt(k,truei)) {
			case 0: 
			    error("gsiCImpAcompClrExpectation: internal error: missing not missing");
			    break;
			case 1: // BDL 
			    if( smallAlr[i] >= smallAlrDL[i] )
				fits=0;
			    break;
			case 2: // MAR
			default: // other treated as MAR
			    break; // no restriction
		    }
		}
		if( fits ) {
		    fitters++;
		    for(i=0;i<D;i++)
			smallClrE[i]+=smallAlr[i];
		    for(i=0;i<D;i++)
			for(j=0;j<D;j++)
			    smallClr2E[i+D*j]+=smallAlr[i]*smallAlr[j];
		}
	    }
	    // d) add final results to clrE, clrVar and copy to xlrImp
	    Vused(k)=fitters;
	    if( fitters > 0 ) { // There were MonteCarlo points
		// add conditional expectation into main matrices
		for(i=0;i<D;i++)
		    MxlrImp(k,MimputationIdx(k,i))=smallClrE[i]/fitters;
		for(i=0;i<D;i++)
		    clrE[MimputationIdx(k,i)]+=smallClrE[i]/fitters;
		for(i=0;i<D;i++)
		    for(j=0;j<D;j++)
			clrVar[MimputationIdx(k,i)+D*MimputationIdx(k,j)]+=
			    smallClr2E[i+D*j]/fitters;
	    } else { // There were no MonteCarlo points
		// Create guess in smallAlr and copy and add
		for(i=0;i<D;i++) {
		    truei=MimputationIdx(k,i);
		    switch( Mmt(k,truei) ) {
			case 0: // use smallAlr
			    break;
			case 1:
			    if( smallAlr[i]<smallAlrDL[i] ) {
				// keep smallAlr, since E[alr|.]<smallAlr
			    } else {
				// keep DL, since E[alr|.] <DL
				smallAlr[i]=smallAlrDL[i];
			    }
			    break;
			case 2:
			default:
			    // This is incorrect if dependent bdls exist
			    smallAlr[i]=smallAlrPre[i];
		    }
		    clrE[truei]+=smallAlr[i];
		    MxlrImp(k,truei)=smallAlr[i];
		}
		for(i=0;i<D;i++)
		    for(j=0;j<D;j++)
			clrVar[MimputationIdx(k,i)+D*MimputationIdx(k,j)]+=
			    smallAlr[i+D*j];
	    }
	}
    }
    // Correct output
    // normalize
    for(j=0;j<D;j++)
	clrE[j]/=n;
    // E^2 to var
    for(i=0;i<D;i++)
	for(j=0;j<D;j++)
	    clrVar[i+D*j]=clrVar[i+D*j]/n-clrE[i]*clrE[j];
    // center 
    tmp=0;
    for(j=0;j<D;j++)
	tmp+=clrE[j];
    tmp/=D;
    for(j=0;j<D;j++)
	clrE[j]-=tmp;
}

/* gsiCImpAcompCreateArrays
   A function to compute the structure variables used in the ImAcomp
   subsystem.


*/

extern void gsiCImpAcompCreateArrays(
    int *Dp,     // ncol(comp)
    int *Np,     // nrow(comp)
    double *comp,// comp
    int *mt,     // output: MissingType
    double *dl,  // output: detectionLimit
    double *dlGen// A detectionlimit attribute
    ) {
    const int D=*Dp;
    const int N=*Np;
    int i,j;
    for(i=0;i<N;i++)
	for(j=0;j<D;j++) {
	    int idx=i+N*j;
	    double v=comp[idx];
	    dl[idx]=0.0;
	    if( R_FINITE(v) ) {
		if( v > 0 ) { //NMV
		    mt[idx]=0;
		} else if( v < 0 ) { // BDL
		    mt[idx]=1;
		    dl[idx]=-v;
		} else  { // BDL
		    mt[idx]=1;
		    dl[idx]=*dlGen;
		}
	    } else if( ISNA(v) ) { // MNAR
		mt[idx]=5;
	    } else if( v==R_PosInf ) { // ERR oder MAR (later)
		mt[idx]=4;
	    } else if( v==R_NegInf ) { // SZ
		mt[idx]=3;
	    }else if( ISNAN(v) ) { // MCAR (MAR)
		mt[idx]=2;
	    }
	}
}

static int myIntLog2(int n) {
    int i=0;
    for(i=0;i<64;i++)
	if( ! (n>>=1) )
	    return(i);
    error("Internal Error myIntLog2");
    return(0);
}


static int NLVL(int m,int n) {
    int MINMN= m<n?m:n;
    int log2=myIntLog2( MINMN/2)+1;
    if( log2 < 0 ) 
	return 0; 
    else 
	return log2;
} 

static int iworkForDgelsd(int m,int n) {
    int MINMN= m<n?m:n;
    return(3*MINMN*NLVL(m,n)+11*MINMN);
} 

extern void gsiCImpAcompFitWithProjection(
    int *Dp,            // 1: Dimension of the compostion
    int *np,            // 1: number of cases
    int *Np,            // 1: size of matrices
    int *s,             // np: case selector (max(s[])<N)
    double *comp,       // N*D: The observed values replaced if missing
    int    *mt,         // N*D: 0=ok, 1=BDL, 2=MAR, 3-5=other (treated as MAR) 
    double *clr,        // Out: N*D: predicted clr on output
    int    *mtType,     // Out: N: The missingType of each case
    int    *imputationIdx,   // Dxn indices of missings 
    int    *imputationN,  // Number of missings
    int    *NmissingTypes, // number of different missing types
    int    *dimX,        // dimension of X
    double *X,           // nxp Designmatrix of model
    double *beta         // output pxD parameter matrix
    ) {
    const int nMissingTypes=*NmissingTypes;
    const int D=*Dp;     // columns in dataset
    const int n=*np;    // number of cases to be processesed 
    const int N=*Np;    // leading dimension of mt,xlr,alrDl,xlrImp
    const int p=dimX[1];
    const char transA='N';
    int eins=1;
    const double one=1.0;
    const double zero=0.0;
    int i,j;            // parts
    int c,l;            // cases, parameters
    double *Xbig = (double*) R_alloc(D*D*n*p,sizeof(double));
    int *iwork= (int *) R_alloc(iworkForDgelsd(D*n,D*p)*3,sizeof(int));
    // zero Xbig
    for(c=0;c<n;c++)
	for(l=0;l<p;l++)
	    for(i=0;i<D;i++)
		for(j=0;j<D;j++)
		    TXbig(i,j,c,l)=0.0;
    // Create Equation system
    for(c=0;c<n;c++) {   
	int type=VmtType(c);
	int nMissings=VimputationN(type);
	int nonMissing=D-nMissings;
	// Create big projector an
	double obnm= -(nonMissing>0?1.0/nonMissing:1.0);
	for(l=0;l<p;l++) {
	    double xv = MX(c,l);
	    for(i=0;i<D;i++)
		for(j=0;j<D;j++)
		    if(Mmt(c,i)==0 && Mmt(c,j)==0 )
			TXbig(i,j,c,l)=(i==j)?obnm + 1.0:obnm;
	}
	// create clr
	if( nonMissing== 0) {
	    for(i=0;i<D;i++) 
		Mclr(c,i)=0.0;
	} else {
	    double tmp=0.0;
	    for(i=nMissings;i<D;i++)
		tmp+=log(Mcomp(c,i));
	    for(i=0;i<D;i++) {
		int truei=MimputationIdx(type,i);
		if( Mmt(c,truei) == 0 )
		    Mclr(c,truei)=log(Mcomp(c,truei))-tmp;
		else
		    Mclr(c,truei)=0.0;
	    }
	}
    }
    // Solve equation system
    int rank=0;
    double rcond=1E-8;
    int info=0;
    int XbigM=D*n;
    int XbigN=D*p;
    int ssize=XbigM<XbigN?XbigM:XbigN;
    int lwork=-1;
    double myWork[3]={0.0,0.0,0.0};
    double *work=myWork;
    F77_CALL(dgelsd)(
	&XbigM, // int *m, 
	&XbigN, // int *n, 
	&eins,  // int *nrhs,
	Xbig,   // double *a, 
	&XbigM, //int *lda, 
	clr,    // double *b, 
	&XbigM, // int *ldb, 
	work,   //double *s, 
	&rcond, //const double* rcond, 
	&rank, //int* rank,
	work+1, //double* work, 
	&lwork,//const int* lwork, 
	iwork,// integer working buffer 3 * MINMN * NLVL + 11 * MINMN
	&info //int* info
	);
    if( info != 0)
	error("gsiCImpAcompFitWithProjection: Problem in workspace query");
    lwork=(int)(work[0]);
    work=(double*)R_alloc(lwork+ssize,sizeof(double));
    F77_CALL(dgelsd)(
	&XbigM, // int *m, 
	&XbigN, // int *n, 
	&eins,  // int *nrhs,
	Xbig,   // double *a, 
	&XbigM, //int *lda, 
	clr,    // double *b, 
	&XbigM, // int *ldb, 
	work+lwork,   //double *s, 
	&rcond, //const double* rcond, 
	&rank, //int* rank,
	work, //double* work, 
	&lwork,//const int* lwork, 
	iwork,// integer working buffer 3 * MINMN * NLVL + 11 * MINMN
	&info //int* info
	);
    if( info!= 0 )
	error("gsiCImpAcompFitWithProjection: Problem in finding solution");
    // computed prediction values (for all cases!)
    F77_CALL(dgemm)(
		    "N",//&transA, // const char *transa, 
		    "N",//&transA, // const char *transb, 
	&N,      //const int *m,
	&D,      // const int *n, 
	&p,      // const int *k, 
	&one,    // const double *alpha,
	X,       // const double *a, 
	&N,      // const int *lda,
	beta,    // const double *b, 
	&p,      // const int *ldb,
	&zero,   // const double *beta, 
	clr,     // double *c, 
	&N       // const int *ldc
	FCONE FCONE
	);
    
	
} 



/* gsiCFitWithProjection 
   to be called from R to use the Fitting with projection technique.
   This technique is only really valid for MCAR. 

 */ 
extern void gsiCFitWithProjection(
    int *dimComp, // dim(comp)=c(N,D)
    double *comp, // the data in compositions format
    double *pred, // Output: the predicted data (including type 1 imputations)
    int *lenS,    // typically N
    int *s,       // typicall 1:N the case to use for the fitting
    int *dimX,    // dim(X) =c(N,p) 
    double *X,    // Designmatrix
    double *beta, // Output: coefficients
    double *dlGen // 1: the general detection limit
    ) {
    int D=dimComp[1];     // columns in dataset
    int n=lenS[0];    // number of cases to be processesed 
    int N=dimComp[0];    // leading dimension of mt,xlr,alrDl,xlrImp
    int p=dimX[1];
    int NmissingTypes=0;
    int *mt=(int *) R_alloc(N*D,sizeof(int));
    int *sMiss=(int *) R_alloc(N,sizeof(int));
    int *mtType=(int *) R_alloc(N,sizeof(int));
    int nTypes=0;
    double *dl=(double *) R_alloc(N*D,sizeof(double));
    gsiCImpAcompGetTypes(
	&D, //int *Dp,      // 1: Number of dimensions
	&n, //int *np,      // 1: Number of cases
	&N, //int *Np,      // 1: leading dimension of mt
	s,  //int *s,       // n: Cases to be used
	mt, //int *mt,      // N x D a missingness signature
	sMiss,//int *snew,    // Out: np or nMissingTypes: the lines of mt for each situation
	mtType,//int *mtType,  // Out: N: The missingType of each case
	&nTypes //int *nMissingTypes // Out 1: 
	);
    int *imputationIdx=(int *) R_alloc(nTypes*D,sizeof(int));
    int *imputationN  =(int *) R_alloc(nTypes,sizeof(int));
    gsiCImpAcompGetIdx(
	&D,     //int *Dp,      // 1: Number of dimensions
	&nTypes,//int *np,      // 1: Number of missingTypes
	&N,     //int *Np,      // 1: leading dimension of mt
	sMiss, //int *s,       // n: Cases to be used for each missingType
	mt,    //int *mt,      // NxD a missingness signature
	imputationIdx, // int    *imputationIdx,   // Out: nxD: new sequence of columns for each missingType 
	imputationN  //int    *imputationN      // Out: n: Number of missings
	);
    //double *imputationCache=(double *) R_alloc(D*D*nTypes*2,sizeof(double));
    double *clr=(double *) R_alloc(N*D,sizeof(double));
    gsiCImpAcompFitWithProjection(
	&D, //int *Dp,            // 1: Dimension of the compostion
	&n, //int *np,            // 1: number of cases
	&N, //int *Np,            // 1: size of matrices
	s,  //int *s,             // np: case selector (max(s[])<N)
	comp,//double *comp,       // N*D: The observed values replaced if missing
	mt,  //int    *mt,         // N*D: 0=ok, 1=BDL, 2=MAR, 3-5=other (treated as MAR) 
	clr, //double *clr,        // Out: N*D: predicted clr on output
	mtType, //int    *mtType,     // Out: N: The missingType of each case
	imputationIdx,//int    *imputationIdx,   // Dxn indices of missings 
	imputationN,  //int    *imputationN,  // Number of missings
	&nTypes,  //int    *NmissingTypes, // number of different missing types
	dimX,     //int    *dimX,        // dimension of X
	X,        //double *X,           // nxp Designmatrix of model
	beta      //double *beta         // output pxD parameter matrix
	);
    int i,j;
    double tmp=0.0;
    for(i=0;i<N;i++) {
	int idx=i+N*j;
	tmp=0.0;
	for(j=0;j<D;j++) 
	    tmp += (pred[i+N*j] = exp(clr[i+N*j]));
	for(j=0;j<D;j++) 
	    pred[i+N*j]/=tmp;
    }
}



extern void gsiCFitWithEM(
    int *steps,   // number of EM-steps
    int *dimComp, // dim(comp)=c(N,D)
    double *comp, // the data in compositions format
    double *pred, // Output: the predicted data (including type 1 imputations)
    int *lenS,    // typically N
    int *s,       // typicall 1:N the case to use for the fitting
    int *dimX,    // dim(X) =c(N,p) 
    double *X,    // Designmatrix
    double *beta, // Output: coefficients
    double *dlGen // 1: the general detection limit
    ) {
    int D=dimComp[1];     // columns in dataset
    int n=lenS[0];    // number of cases to be processesed 
    int N=dimComp[0];    // leading dimension of mt,xlr,alrDl,xlrImp
    int p=dimX[1];
    int NmissingTypes=0;
    int *mt=(int *) R_alloc(N*D,sizeof(int));
    int *sMiss=(int *) R_alloc(N,sizeof(int));
    int *mtType=(int *) R_alloc(N,sizeof(int));
    int nTypes=0;
    double *dl=(double *) R_alloc(N*D,sizeof(double));
    gsiCImpAcompGetTypes(
	&D, //int *Dp,      // 1: Number of dimensions
	&n, //int *np,      // 1: Number of cases
	&N, //int *Np,      // 1: leading dimension of mt
	s,  //int *s,       // n: Cases to be used
	mt, //int *mt,      // N x D a missingness signature
	sMiss,//int *snew,    // Out: np or nMissingTypes: the lines of mt for each situation
	mtType,//int *mtType,  // Out: N: The missingType of each case
	&nTypes //int *nMissingTypes // Out 1: 
	);
    int *imputationIdx=(int *) R_alloc(nTypes*D,sizeof(int));
    int *imputationN  =(int *) R_alloc(nTypes,sizeof(int));
    gsiCImpAcompGetIdx(
	&D,     //int *Dp,      // 1: Number of dimensions
	&nTypes,//int *np,      // 1: Number of missingTypes
	&N,     //int *Np,      // 1: leading dimension of mt
	sMiss, //int *s,       // n: Cases to be used for each missingType
	mt,    //int *mt,      // NxD a missingness signature
	imputationIdx, // int    *imputationIdx,   // Out: nxD: new sequence of columns for each missingType 
	imputationN  //int    *imputationN      // Out: n: Number of missings
	);
    double *imputationCache=(double *) R_alloc(D*D*nTypes*2,sizeof(double));
    double *clr=(double *) R_alloc(N*D,sizeof(double));
    gsiCImpAcompFitWithProjection(
	&D, //int *Dp,            // 1: Dimension of the compostion
	&n, //int *np,            // 1: number of cases
	&N, //int *Np,            // 1: size of matrices
	s,  //int *s,             // np: case selector (max(s[])<N)
	comp,//double *comp,       // N*D: The observed values replaced if missing
	mt,  //int    *mt,         // N*D: 0=ok, 1=BDL, 2=MAR, 3-5=other (treated as MAR) 
	clr, //double *clr,        // Out: N*D: predicted clr on output
	mtType, //int    *mtType,     // Out: N: The missingType of each case
	imputationIdx,//int    *imputationIdx,   // Dxn indices of missings 
	imputationN,  //int    *imputationN,  // Number of missings
	&nTypes,  //int    *NmissingTypes, // number of different missing types
	dimX,     //int    *dimX,        // dimension of X
	X,        //double *X,           // nxp Designmatrix of model
	beta      //double *beta         // output pxD parameter matrix
	);
    int iter;
    for(iter=0;iter<*steps;iter++) {
	
    }
    int i,j;
    double tmp=0.0;
    for(i=0;i<N;i++) {
	int idx=i+N*j;
	tmp=0.0;
	for(j=0;j<D;j++) 
	    tmp += (pred[i+N*j] = exp(clr[i+N*j]));
	for(j=0;j<D;j++) 
	    pred[i+N*j]/=tmp;
    }
}

     
/*
extern void gsiCAcompFitAndImpute(
    int *steps,         // 1: number of steps to compute, if negativ continue
    int *Dp,            // 1: Dimension of the compostion
    int *np,            // 1: number of cases
    int *Np,            // 1: size of matrices
    int *s,             // np: case selector (max(s[])<N)
    double *comp,       // N*D: The observed values replaced if missing
    int    *mt,         // N*D: 0=ok, 1=BDL, 2=MAR, 3-5=other (treated as MAR) 
    double *dl,         // N*D: Detectionlimit if BDL
    double *clr,        // output n*D:  imputed clr 
    double *clrV,       // output D*D:  the variance estimate
    int *pp,            // 1: the columns of the design matrix of linear model 
    double *X,          // nxp the design matrix
    double *beta,       // output: p*D the parameter estimate
    double *XXinv,       // pxp (X^tX)^-1 might be a g-inverse
    int *dimNorm,        // 2: the dimension of the normalized residuals 
    double *norm        // normalized residuals    
    ) {
    int n=*np;
    int D=*Dp;
    int N=*Np;
    int p=*pp;
    int i,j,k;
    if( *steps >= 0 ) {
	// Initialize with detection limit.
	for(i=0;i<n;i++)
	    for(j=0;j<D;j++)
		switch(Mmt(i,j)) {
		    case 0: break;
		    case 1: Mcomp(i,j)=Mdl(i,j); break;
		    default: 
	
    } else {
	*steps = -*steps;
    } 
    for(;steps>0;steps--) {
	for(i=0;i<n;i++)
	    for(j=0;j<D;j++)
		;

    }
}
*/

/* Aitchison Distribution integrals */ 

static double gridValue(int gc,int grid) {
   return log((gc+(double)1)/grid);
   // return log((double)(gc+1)/grid);
}

// AitchisonDistributionIntegral
static R_NativePrimitiveArgType gsiAitchisonDistributionIntegral_t[] = {
  INTSXP,INTSXP,INTSXP,REALSXP,REALSXP,REALSXP,REALSXP,REALSXP,REALSXP
};


extern void gsiAitchisonDistributionIntegral(
    int *Dp,               // 1
    int *gridp,            // 1
    int *modep,            // 1
    double *theta,        // D
    double *beta,         // DxD
    double *OneIntegral,  // 1
    double *kappaIntegral, // 1
    double *clrIntegral,  // D
    double *SqIntegral // DxD
    ) {
    int D=*Dp;
    int mode=*modep;
    int grid=*gridp;
    int gridPlus=grid+D;
    int *gc=(int*)R_alloc(D,sizeof(int));
    double *gv=(double*)R_alloc(D,sizeof(double));
    int i,j,k;
    int p=0;
    int steped=0;
    int carry=0;
    double logMean=0.0;
    double tmp=0.0;
    double phi=0.0;
    long gridSize=0;
    double density=0.0;
    /*zero*/
    *OneIntegral   = 0.0;
    *kappaIntegral = 0.0;
    if( mode < 0 )
	  return;
    if( mode >= 1) 
	  for(i=0;i<D;i++)
	    clrIntegral[i] = 0.0;
    if( mode >= 2) 
   	  for(i=0;i<D;i++)
	    for(j=0;j<D;j++)
		  SqIntegral[i+D*j] = 0.0;
    /* Init Grid */
    for(p=0;p<D;p++)
	  gc[p]=0.0;
    gc[0]=grid;
    for(p=0;p<D;p++)
	  gv[p]=gridValue(gc[p],gridPlus);
    /* Check Matrices */
    for(i=0;i<D;i++) {
	  tmp=0.0;
	  for(j=0;j<D;j++) {
	    if( fabs(beta[i+j*D]-beta[j+i*D])>1E-6 )
		  error("gsiAitchisonDistributionIntegral: beta not symmetric");
	    tmp+=beta[i+j*D];
	  } // end for j
	  if(fabs(tmp)>1E-10) {
		error("gsiAitchisonDistributionIntegral: beta not clr matrix");
	  }
    } // end for i
    /* Loop over grid */
    p = 0;
    for(;;) {
	/* Generate Grid */
  	  steped=0;
	  for(p=0;p<D-1;p++) {
	    if( gc[p] > 0 ) {
		  carry=gc[p]-1;
		  gc[p+1]++;
		  gc[p]=0;
		  gc[0]=carry;
		  // Rprintf("%d", gc[0]);
		  gv[p+1]=gridValue(gc[p+1],gridPlus);
		  gv[p]=gridValue(gc[p],gridPlus);
		  gv[0]=gridValue(gc[0],gridPlus);
		  steped=1;
		  break;
	    } // end if
	  } // end for p
	      //Rprintf("%lf %lf %lf\n",gv[0],gv[1],gv[2]);
	  if( ! steped )
	    break;
	  /* Compute CLR const */
	  logMean=0.0;
	  phi = 0.0;
	  for(i=0;i<D;i++) {
	    logMean+=gv[i];
	    phi += gv[i]*(theta[i]-1);
	    for(j=0;j<D;j++) {
		  phi+= gv[j]*gv[i]*beta[i+j*D];
	    } // end for j
	  } // end for i
	  logMean/=D;
	  /* Cummulate*/
	  gridSize++;
	  density=exp(phi);
	  *OneIntegral   += density;
	  *kappaIntegral += density*logMean;
	  if(mode>=1)
	    for(i=0;i<D;i++) 
		  clrIntegral[i]+=density*(gv[i]-logMean);
	  if(mode>=2)
	    for(i=0;i<D;i++) 
		  for(j=0;j<D;j++)
		    SqIntegral[i+D*j]+=density*(gv[i]-logMean)*(gv[j]-logMean);
    } // end for (empty)
    /* Rescale*/
    if(mode>=1)
  	  for(i=0;i<D;i++) 
	    clrIntegral[i]/=*OneIntegral;
    if(mode>=2)
	  for(i=0;i<D;i++) 
	    for(j=0;j<D;j++) {
		  SqIntegral[i+D*j]/=*OneIntegral;
		  if( mode >=3 )
		    SqIntegral[i+D*j]-=clrIntegral[i]*clrIntegral[j];
	    } // end for j
    *kappaIntegral/=*OneIntegral;
    *OneIntegral /= gridSize;
}

/*
extern void gsiMCAitchisonDistributionIntegral(
    int *D,               // 1
    int *mode,            // 1
    double *nMC,          // 1
    double *compMC,       // nMC x D
    double *dMC,          // nMC
    double *theta,        // D
    double *beta,         // DxD
    double *OneIntegral,  // 1
    double *kappaIntegral, // 1
    double *clrIntegral,  // D
    double *SqIntegral // DxD
    ) {
    double sf=1/(double)*grid;
    int *gc=R_alloc(D,sizeof(int));
    int *gv=R_alloc(D,sizeof(double));
    int i,j,k;
    int p;
    // zero/
    *OneIntegral   = 0.0;
    *kappaIntegral = 0.0;
    if( *mode >= 1) 
	for(i=0;i<D;i++)
	    clrIntegral[i] = 0.0;
    if( *mode >= 2) 
	for(i=0;i<D;i++)
	    for(j=0;j<D;j++)
		SqIntegral[i+D*j] = 0.0;
    // Init Grid 
    for(p=0;p<D;p++)
	gc[p]=0;
    gc[0]=*grid;
    for(p=0;p<D;p++)
	gv[p]=gridvalue(gv[p],*grid);
    // Loop over grid 
    p = 0;
    

}

*/

// AitchisonDistributionIntegral
static R_NativePrimitiveArgType gsirandomClr1Aitchison_t[] = {
  INTSXP,INTSXP,REALSXP,REALSXP,REALSXP,REALSXP
};


extern void gsirandomClr1Aitchison(
    int *Dp, // 1
    int *np, // 1
    double *erg, // n x D
    double *theta, // D
    double *mu,  // D
    double *sqrtSigma // DxD
    ) {
    int D=*Dp;
    int n=*np;
    int i,j,k;
    double z=0.0;
    double logD,kappa;
    double alphaD=0.0;
    // Sum over theta
    for(i=0;i<D;i++)
	  alphaD+=theta[i];  // alphaD = sum of thetas = alpha * D
      // Find Maximum Density of Dirichlet
    if( alphaD<0 )
      error("gsirandomClr1Aitchison: theta must have a positive sum");
    if( alphaD<1E-5)
      alphaD=1;
    logD=0.0;
    for(j=0;j<D;j++) 
      logD+=theta[j]*(log(theta[j])-log(alphaD));
    double MaxDir  = exp(logD);  // maximal value of Ait/Normal-equivalence class
    // Do Simulation 
    GetRNGstate();
    // Loop over cases
    for(k=0;k<n;) {
	  // zero erg
	  for(j=0;j<D;j++)
	    erg[n*j]=mu[j];  // for the moment, that mean is zero; possibly generalizable in the future to an arbitrary splitting of Ait = Dirichlet x Normal
	  // Simulate the random clr
	  for(i=0;i<D;i++) {
	    z = norm_rand();
	    for(j=0;j<D;j++){
		  erg[n*j]+=sqrtSigma[j+i*D]*z;  // given that sqrtSigma is a (clr-variance)^(1/2) it adds up to zero by rows and columns, and the result of this matrix product is thus a clr-vector
          // Rprintf("%f\n", erg[n*j]);
	    }
	  } // end for i
  	  // compute Dirichlet density
	  logD = 0.0;  
	  kappa = 0.0;  
	  for(j=0;j<D;j++) {
	    kappa += exp(erg[n*j]);  // accumulate exp(clr(x_i))
	    logD += theta[j]*erg[n*j]; // accumulate theta_i*clr_i(x)
	  } // end for j
	  kappa = log(kappa);  
	  logD -= alphaD*kappa;  // logD += (alphaD-D)*kappa;
	  double Dir  = exp(logD);
	  if( Dir > MaxDir )
	    error("randomClrAitchison: Internal Error Density exceeds Maximum, please report to package authors");
	  // 
      // Rprintf("(%f, %f, %f) \n", erg[n*0],erg[n*1],erg[n*2]);
	  // Rprintf("%f %f \n", Dir, MaxDir);
	  if( Dir/MaxDir >= unif_rand() ) {
	    k++;
	    erg++;
	  } // end if 
    } // end for k
    PutRNGstate();
}



extern void gsiDensityCheck1(
    int *dimX,        //input 
    double *X,        //input
    double *fX,       //output
    double *fXq,       //output
    double *fhatX,    //output
    int *dimY,        //
    double *Y,        //Vergleichsdatensatz
    double *bw,       //Bandweite
    double *DQ        // Ergebnis
    ) {
    int i,j,k;
    int n=dimX[0];
    int N=dimY[0];
    int D=dimX[1];
    double normQ;
    double tmp;
    double sigmaQ=2 * *bw * *bw;
    for(i=0;i<n;i++) {
	fX[i]=0.0;
	fXq[i]=0.0;
    }
    for(i=0;i<n;i++) {
	for(j=0;j<N;j++) {
	    normQ=0;
	    for(k=0;k<D;k++) {
		tmp = X[i+k*n]-Y[j+k*N];
		normQ += tmp*tmp;
	    }
	    tmp = exp(-normQ/sigmaQ); // Kerndichte
	    fX[i]+=tmp;
	    fXq[i]+=tmp*tmp;
	}
    }
    for(i=0;i<n;i++) {
	fX[i]/=N;
	fXq[i]/=N;
	fXq[i]-=fXq[i]*fXq[i]; // Varianz
    }
    for(i=0;i<n;i++)
	fhatX[i]=0.0;
    for(i=0;i<n;i++) {
	for(j=0;j<i;j++) {
	    normQ=0;
	    for(k=0;k<D;k++) {
		tmp = X[i+k*n]-X[j+k*n];
		normQ += tmp*tmp;
	    }
	    tmp = exp(-normQ/sigmaQ); // Kerndichte
	    fhatX[i]+=tmp;
	    fhatX[j]+=tmp;
	}
    }
    for(i=0;i<n;i++)
	fhatX[i]/=(n-1);
    *DQ = 0.0;
    for(i=0;i<n;i++) {
	tmp = (fhatX[i]-fX[i]);
	*DQ += tmp*tmp/(fXq[i]/(n-1)+fXq[i]/N);
    }    
}

extern void gsiDensityCheck2(
    int *dimX,        //input 
    double *X,        //input
    double *fhatX,    //output
    double *fX,        //input
    double *bw,       //Bandweite
    double *lQ        // Ergebnis
    ) {
    int i,j,k;
    int n=dimX[0];
    int D=dimX[1];
    double sigmaQ=3 * *bw * *bw;
    double Pd=sqrt(M_PI*sigmaQ);
    double Konstante=1;
    double normQ;
    double tmp;
    for(i=0;i<n;i++)
	fhatX[i]=0.0;
    for(i=0;i<n;i++) {
	for(j=0;j<i;j++) {
	    normQ=0;
	    for(k=0;k<D;k++) {
		tmp = X[i+k*n]-X[j+k*n];
		normQ += tmp*tmp;
	    }
	    tmp = exp(-normQ/sigmaQ)/Konstante; // Kerndichte
	    fhatX[i]+=tmp;
	    fhatX[j]+=tmp;
	}
    }
    for(i=0;i<n;i++)
	fhatX[i]/=(n-1);
    *lQ = 0.0;
    for(i=0;i<n;i++) {
	tmp = log(fhatX[i])-log(fX[i]);
	*lQ += tmp;
    }    
}

double gsiKkernel(const int *dimX,
		  const double *X,
		  const double *Y,
		  const double *bw 
    ) {
    double tmp=0;
    double normQ=0;
    double k=1.0;
    double TwoSigmaQ=2* *bw * *bw;
    double sQ=TwoSigmaQ*M_PI;
    int i;
    for(i=0;i < dimX[1];i++) {
	tmp=X[i*dimX[0]]-Y[i*dimX[0]];
	normQ+=tmp*tmp;
	k*=sQ;
    }
    return exp(-normQ/ TwoSigmaQ)/sqrt(k);
}

void gsiSelectN(const int *l,const int *n, int *s) {
    const int N2=*l - *n;
    int i=0;
    int maxtries = 3* *n;
    int found=0;
    int j,k;
    if( 2* *n > *l ) {
	gsiSelectN(l,&N2,s);
	for(i=0;i<*l;i++)
	    s[i]=!s[i];
    } else {
	for(i=0;i<*l;i++)
	    s[i]=0;
	for(i=0;i<maxtries&&found< *n;i++) {
	    j= (int)(unif_rand()* *l);
	    if( j >= *l || j<0 || s[j] )
		continue;
	    s[j]=1;found++;
	}
	if( found <*n ) {
	    warning("gsiSelectN: Slow sampling used");
	    for(;found< *n;found++) {
		j= (int)(unif_rand()* (*l - found));
		for(k=0;k<*l;k++) {
		    if(!s[k]) {
			if( j )
			    j--;
			else {
			    s[k]=1;
			    break;
			}
		    }
		}
	    }
	}
    }
    found=0;
    for(i=0;i<*l;i++)
	if(s[i]) found++;
    if( found != *n )
	error("gsiSelectN: failed %d %d",found,*n);
}


#define D(i,j) DD[(ik=(i),jk=(j),((ik)<=(jk)?((ik-1)+((jk)*((jk)+1))/2) : ((jk-1)+(((ik)*((ik)+1))/2))))]

static R_NativePrimitiveArgType gsiDensityCheck_t[] = {
  INTSXP,REALSXP,INTSXP,REALSXP,REALSXP,REALSXP,INTSXP,REALSXP
};


extern void gsiDensityCheck(
    int *dimX,
    double *X,
    int *dimY,
    double *Y,
    double *bw,
    double *stat,
    int *Nreps,
    double *reps
    ) {
    int i,j,k,r;
    int ik,jk;
    int n=dimX[0];
    int m=dimY[0];
    int l=n+m;
    int D=dimX[1];
    int *tmpIntPtr;
    int tmpInt;
    double *tmpDoublePtr;
    double P1=0,P2=0,P3=0;
    double tmp;
    double lwork=(l*(l+1))/2;
    double *DD=(double*) R_alloc(lwork,sizeof(double));
    int    *take = (int*) R_alloc(lwork,sizeof(int));
/*    if( dimX[0] > dimY[0] ) {
	// swap
	tmpIntPtr=dimX;dimX=dimY;dimY=tmpIntPtr;
	tmpDoublePtr=X;X=Y;Y=tmpDoublePtr;
	tmpInt=n;n=m;m=tmpInt;
	}    */
    if( dimX[1]!=dimY[1] )
	error("gsiDensityCheck: Error");
    if( DD == NULL || take==NULL) {
	error("gsiDensityCheck= Out of memory");
    }
    for(i=0;i<n;i++) {
	for(j=0;j<=i;j++) 
	    P1+=(D(i,j)=gsiKkernel(dimX,X+i,X+j,bw));
    }
    for(;i<l;i++) {
	for(j=0;j<n;j++)
	    P2+=(D(i,j)=gsiKkernel(dimX,Y+(i-n),X+j,bw));
	for(j=n;j<=i;j++)
	    P3+=(D(i,j)=gsiKkernel(dimX,Y+(i-n),Y+(j-n),bw));
    }
    P1/=(n*n);P2/=(n*m);P3/=(m*m);
    *stat = P2/sqrt(P1*P3);
    if( *Nreps > 0 ) {
	GetRNGstate();
	for(r=0;r<*Nreps;r++) {
	    P1=P2=P3=0;
	    gsiSelectN(&l,&n,take);
	    for(i=0;i<l;i++) {
		if( take[i] ) {
		    for(j=0;j<=i;j++) {
			if( take[j] ) {
			    P1+=D(i,j);
			} else {
			    P2+=D(i,j);
			}
		    }
		} else {
		    for(j=0;j<=i;j++) {
			if( take[j] ) {
			    P2+=D(i,j);
			} else {
			    P3+=D(i,j);
			}
		    }
		}
	    }
	    P1/=(n*n);P2/=(n*m);P3/=(m*m);
	    reps[r] = P2/sqrt(P1*P3);
	}
	PutRNGstate();
    }
}

void gsiSpeedShuffel(const int *intern,const int *n,const int *nmax,int *s) {
    int i,j,k;
    if( !*intern )
	GetRNGstate();
    for(i=0;i<*n;i++) {
	j=(int)(unif_rand()* (*nmax-i));
	k=s[i];s[i]=s[j];s[j]=k;
    }
    if( !*intern )
	GetRNGstate();
}


extern void gsiFitAitchison(
    int *dimX,
    double *X,
    int *dimY,
    double *Y,
    double *theta,
    double *beta,
    double *E1,
    double *E2,
    int * maxIter,
    int * dimV,
    double *V) {

/*Yet to be implemented
  Fits a Aitchison Distribution Density with respect to the reference 
  measure underlying Y.

*/

}

extern void gsiMLClogLik(
    int *dimX,
    double *X,
    int *ng,
    int *g,
    double *means,
    double *vars
    ) {
    



}

// KS Poisson
static R_NativePrimitiveArgType gsiKSPoisson_t[] = {
  INTSXP,INTSXP,INTSXP,REALSXP,INTSXP,REALSXP
};



extern void gsiKSPoisson(
    int *nd,    // Number Datapoints,number datasets
    int *data,  // the Datapoints
    int   *nps, // maximal used probability value
    double *ps, // probabilities
    int *n,     // [max+1]output: the actual count of observed values
    double *statistic // the statistic
    ) {
    const int N=*nd;
    const int M=*nps;
    const int K=nd[1];
    double D,maxD;
    unsigned int v;
    int i,k;
    for(k=0;k<K;k++) {
	D=maxD=0.0;
	for(i=0;i<M;i++) n[i]=0;
	for(i=0;i<N;i++) {
	    v=(unsigned int)(data[i]);
	    if( v<M )
		n[v]++;
	}
	for(i=0;i<M;i++) {
	    D += N*ps[i] - n[i];
	    if( D > maxD )
		maxD=D;
	    else if( -D > maxD )
		maxD=-D;
	}
	*statistic = maxD/N;
	data+=N;
	statistic++;
    }
}

// KS Sorted Uniforms
static R_NativePrimitiveArgType gsiKSsortedUniforms_t[] = {
  INTSXP,REALSXP,INTSXP
};


extern void gsiKSsortedUniforms(
    const int *n,
    double *data,
    const int *getRng
    ){
    int i;
    int N=*n;
    double tmp=0.0;
    if( *getRng ) {
	GetRNGstate();
    }
    for(i=0;i<N;i++)
	(data[i] = (tmp += -log(unif_rand())));
    tmp += -log(unif_rand());
    for(i=0;i<N;i++)
	data[i] /= tmp;
    if( *getRng ) {
	PutRNGstate();
    }
}


// KS Poisson Sample
static R_NativePrimitiveArgType gsiKSPoissonSample_t[] = {
  INTSXP,REALSXP,INTSXP,REALSXP,INTSXP,REALSXP
};

extern void gsiKSPoissonSample(
    const int *n,          // number of data values 
    double *data,    // space for the data
    const int *Nps,        // The number of probabilities 
    double *ps,      // probabilities
    const int *nsamples,   // Number of samples to be drawn
    double *statistics
    ) {
    int i,j,k;
    double tmp=0;
    const int noDo=0;
    double stat;
    GetRNGstate();
    // Cumulate probabilities
    for(i=0;i<*Nps;i++) 
	tmp = (ps[i]+=tmp);
    for(k=0;k<*nsamples;k++) {
	gsiKSsortedUniforms(
	    n,data,&noDo
	    );
	i=j=0;
	stat=0;
	for(i=0;i<*Nps;i++) {
	    for(;j<*n;j++) {
		if( data[j] > ps[i] )
		    break;
	    }
	    tmp=fabs(j/(double)(*n)-ps[i]);
	    if( tmp > stat )
		stat=tmp;
	}
	statistics[k]=stat;
    }
    PutRNGstate();
}



static R_CMethodDef cMethods[] = {
  {"gsiDensityCheck", (DL_FUNC) &gsiDensityCheck, 8, gsiDensityCheck_t},
  {"gsiCGSvariogram", (DL_FUNC) &gsiCGSvariogram, 11, gsiCGSvariogram_t},
  {"gsiCGSvg2lrvg", (DL_FUNC) &gsiCGSvg2lrvg, 3, gsiCGSvg2lrvg_t},
  {"gsiCGSkriging", (DL_FUNC) &gsiCGSkriging, 14, gsiCGSkriging_t},
  {"gsiAitchisonDistributionIntegral", (DL_FUNC) &gsiAitchisonDistributionIntegral, 9, gsiAitchisonDistributionIntegral_t},
  {"gsirandomClr1Aitchison", (DL_FUNC) &gsirandomClr1Aitchison, 6, gsirandomClr1Aitchison_t},
  {"gsiKSPoisson", (DL_FUNC) &gsiKSPoisson, 6, gsiKSPoisson_t},
  {"gsiKSsortedUniforms", (DL_FUNC) &gsiKSPoissonSample, 3, gsiKSsortedUniforms_t},
  {"gsiKSPoissonSample", (DL_FUNC) &gsiKSPoissonSample, 6, gsiKSPoissonSample_t},
  {NULL, NULL, 0}
};


void R_init_compositions(DllInfo *info)
  {
    R_registerRoutines(info, cMethods, NULL, NULL, NULL);
    R_useDynamicSymbols(info, FALSE);
    R_forceSymbols(info, TRUE);
  }




