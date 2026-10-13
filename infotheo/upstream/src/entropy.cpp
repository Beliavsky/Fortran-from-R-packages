#include "infotheo.h"

SEXP entropyR (SEXP Rdata, SEXP Rnrows, SEXP Rncols, SEXP Rchoice)
{
      const int *data;
      const int *nrows, *ncols, *choice;
         SEXP res;
         PROTECT(Rdata = AS_INTEGER(Rdata));
         PROTECT(Rnrows= AS_INTEGER(Rnrows));
         PROTECT(Rncols= AS_INTEGER(Rncols));
	 PROTECT(Rchoice= AS_INTEGER(Rchoice));
         data = INTEGER_POINTER(Rdata);
         nrows= INTEGER_POINTER(Rnrows);
         ncols= INTEGER_POINTER(Rncols);   
	 choice= INTEGER_POINTER(Rchoice);    
         PROTECT(res = NEW_NUMERIC(1));
		 bool *sel = new bool[*ncols];
		 for( int i=0; i<*ncols; ++i )
			sel[i] = true;
		 REAL(res)[0] = entropy(data, *nrows, *ncols, *choice, sel);
		 UNPROTECT(5);
      return res;
}

SEXP multiinformationR (SEXP Rdata, SEXP Rnrows, SEXP Rncols, SEXP Rchoice)
{
      const int *data;
      const int *nrows, *ncols, *choice;
         SEXP res;
         PROTECT(Rdata = AS_INTEGER(Rdata));
         PROTECT(Rnrows= AS_INTEGER(Rnrows));
         PROTECT(Rncols= AS_INTEGER(Rncols));
		 PROTECT(Rchoice= AS_INTEGER(Rchoice));
         data = INTEGER_POINTER(Rdata);
         nrows= INTEGER_POINTER(Rnrows);
         ncols= INTEGER_POINTER(Rncols);   
		 choice= INTEGER_POINTER(Rchoice);    
         PROTECT(res = NEW_NUMERIC(1));
		 REAL(res)[0] = multiinformation(data, *nrows, *ncols, *choice);
		 UNPROTECT(5);
      return res;
}

SEXP interactionR (SEXP Rdata, SEXP Rnrows, SEXP Rncols, SEXP Rchoice)
{
      const int *data;
      const int *nrows, *ncols, *choice;
         SEXP res;
         PROTECT(Rdata = AS_INTEGER(Rdata));
         PROTECT(Rnrows= AS_INTEGER(Rnrows));
         PROTECT(Rncols= AS_INTEGER(Rncols));
		 PROTECT(Rchoice= AS_INTEGER(Rchoice));
         data = INTEGER_POINTER(Rdata);
         nrows= INTEGER_POINTER(Rnrows);
         ncols= INTEGER_POINTER(Rncols);   
		 choice= INTEGER_POINTER(Rchoice);    
         PROTECT(res = NEW_NUMERIC(1));
		 REAL(res)[0] = interaction(data, *nrows, *ncols, *choice);
		 UNPROTECT(5);
      return res;
}

SEXP buildMIM(SEXP Rdata, SEXP Rnrows, SEXP Rncols, SEXP Rchoice)
{
      const int *data;
      const int *nrows, *ncols, *choice;
      double *ent, *res, mi;
	  bool *sel;
         SEXP Rres;
         PROTECT(Rdata = AS_INTEGER(Rdata));
         PROTECT(Rnrows = AS_INTEGER(Rnrows));
         PROTECT(Rncols = AS_INTEGER(Rncols));
		 PROTECT(Rchoice= AS_INTEGER(Rchoice));
         data = INTEGER_POINTER(Rdata);
         nrows= INTEGER_POINTER(Rnrows);
         ncols= INTEGER_POINTER(Rncols);     
		 choice= INTEGER_POINTER(Rchoice); 
         PROTECT(Rres = NEW_NUMERIC((*ncols)*(*ncols)));
         res = NUMERIC_POINTER(Rres);
		 ent = new double[*ncols];
		 sel = new bool[*ncols];
		 for( int i=0; i<*ncols; ++i ){
			res[i*(*ncols)+i]=0;
			sel[i] = false;
		 }
		 for( int i=0; i<*ncols; ++i ){
			sel[i] = true;
			res[i*(*ncols)+i] = entropy(data, *nrows, *ncols, *choice, sel);
			sel[i] = false;
		 }
	     for( int i=1; i<*ncols; ++i ){
			sel[i] = true;
			for( int j=0; j<i; ++j ) {
				  sel[j] = true;
                  res[j*(*ncols)+i] = res[i*(*ncols)+j] = res[i*(*ncols)+i] + res[j*(*ncols)+j] - entropy(data, *nrows, *ncols, *choice, sel);
				  sel[j] = false;
			}
			sel[i] = false;
         }
         UNPROTECT(5);
      return Rres;
}


double digamma(double z) {
      if(z<=0) return 0;
      double zp, zpr, zprs, digam = 0;
         zp = z;
      while(zp < 30) {
             zpr = 1/zp;
             digam -= zpr;
             zp++;
      }
         zpr = 1/zp;
         zprs = zpr * zpr;
         digam += log(zp)+zpr*(-0.5+zpr*(-1.0/12.0+zprs*(1.0/120.0-zprs/252.0)));
      return digam;
}

double entropy_empirical(std::map< std::vector<int> ,int > frequencies, int nb_samples) {
      double e = 0;
      for (std::map< std::vector<int> ,int>::const_iterator iter = frequencies.begin(); iter != frequencies.end(); ++iter)
            e -= iter->second * log((double)iter->second);
      return log((double)nb_samples) + e/nb_samples;
}

double entropy_miller_madow(std::map< std::vector<int> ,int > frequencies, int nb_samples) {
      return entropy_empirical(frequencies,nb_samples) + (int(frequencies.size())-1)/(2.0*nb_samples);
}

double entropy_dirichlet(std::map< std::vector<int> ,int > frequencies, int nb_samples, double beta) {
      double e = 0;
      for (std::map< std::vector<int> ,int>::const_iterator iter = frequencies.begin(); iter != frequencies.end(); ++iter)
            e+=(iter->second+beta)*(digamma(nb_samples+(frequencies.size()*beta)+1)-digamma(iter->second+beta+1));
      return e/(nb_samples+(frequencies.size()*beta));
}

double entropy_shrink(std::map< std::vector<int> ,int > frequencies, int nb_samples) 
{
      double w = 0;
      int p = frequencies.size(), n2 = nb_samples*nb_samples; 
      double lambda, beta;
      for (std::map< std::vector<int> ,int>::const_iterator iter = frequencies.begin(); iter != frequencies.end(); ++iter) 
            w += iter->second*iter->second;
         lambda = p*(n2 - w)/((nb_samples-1)*(w*p - n2));
      if(lambda >= 1)
        return -log(1.0/p);
      else {
            beta = (lambda/(1-lambda))*nb_samples/frequencies.size();
        return entropy_dirichlet(frequencies, nb_samples, beta);
      }
}

double entropy(const int *d, int nsamples, int nvars, int c, bool *v) { 
// H(d) using estimator c
	std::map< std::vector<int> ,int > freq;
	std::vector<int> sel;
	bool ok = true;
	int nsamples_ok = 0;
	double H = 0;
	for(int s = 0; s < nsamples; ++s) {
		ok = true;
		sel.clear();
		for(int i = 0; i < nvars; ++i) {
			if(v[i]){
				if(d[s+i*nsamples]!=NA_INTEGER)
					sel.push_back(d[s+i*nsamples]); 
				else
					ok = false;
			}
		}
		if(ok) {
			freq[sel]++; 
			nsamples_ok++;
		}
	}
	if( c == 0 ) //empirical
		H = entropy_empirical(freq,nsamples_ok);
	else if( c == 1 ) //miller-madow
		H = entropy_miller_madow(freq,nsamples_ok);
	else if( c == 2 ) //dirichlet Schurmann-Grassberger
		H = entropy_dirichlet(freq,nsamples_ok, 1/freq.size());
	else if( c == 3 ) // shrink
		H = entropy_shrink(freq,nsamples_ok);
	return H;
}

double multiinformation(const int *d, int nsamples, int nvars, int c) {
		 bool *sel = new bool[nvars];
		 double sum = 0;
		 for( int i=0; i<nvars; ++i )
			sel[i] = false;
		 for(int i=0;i<nvars; ++i) {
			sel[i] = true;
			sum += entropy(d, nsamples, nvars, c, sel);
			sel[i] = false;
		 }	
		 for( int i=0; i<nvars; ++i )
			sel[i] = true;
		 sum -= entropy(d, nsamples, nvars, c, sel);
		 return sum;
}



double interaction(const int *d, int nsamples, int nvars, int c) {
	double sum = 0; //
	int sign = 1;//
	bool *sel= new bool[nvars];//
	for(int i = 0; i < nvars; ++i)//
		sel[i] = false; //
	
	std::vector<int> indx;
	int n = nvars;
	int j=1;
	int k=n;
	for(int twk=j;twk<=k;twk++){
		int r=twk;
		bool done=true;
		for(int iwk=0;iwk<r;iwk++)indx.push_back(iwk);
		while(done){
			done=false;
			for(int owk=0;owk<r;owk++) //
				sel[indx[owk]]=true;	//
			sum += sign*entropy(d, nsamples, nvars, c, sel); //
			for(int owk=0;owk<r;owk++) //
				sel[indx[owk]]=false;	//
			for(int iwk=r-1;iwk>=0;iwk--){
				if(indx[iwk]<=(n-1)-(r-iwk)){
					indx[iwk]++;
					for(int swk=iwk+1;swk<r;swk++)
						indx[swk]=indx[swk-1]+1;
					iwk=-1;
					done=true;
				}	
			}	
		}
		sign = -sign; //
		indx.clear();
	}
	return sum;
}
