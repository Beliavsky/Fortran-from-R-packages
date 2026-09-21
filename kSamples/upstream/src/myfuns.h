/* functions prototypes */

double *dmalloc(unsigned long n);
int getCount(double z, double *dat, int n);
int getSum(int *x, int n);
int *imalloc(unsigned long n);
void mcopy(double *x, double *copy, int n, int m);
void imcopy(int *x, int *copy, int n, int m);
double **pdmalloc(unsigned long n);
int exactcomb(int now, int *position, int m, int dimst, 
		void(*testStatFun)(double *teststat, int k, 
           double *x, int *ns));
void initvals(int kk, double *xx, int *nns, double *zzstar, 
           int ll, int *iivec, double *tteststat, double *ppval,
           int nnsum, int nncomb, int ggetSmat, double *ssmat,
	   int ddimst);
void initvals1(int kk, double *xx, int *nns, double *zzstar, 
              int LL);
int runCount(void);
void adkTestStat(double *adk, int k, double *x, int *ns);
void adkPVal(double *pval, int Nsim, int k, double *x, int *ns,
		  double *zstar, int l, int useExact, int getA2mat, 
		  double ncomb, double *a2mat);
void contingency2xtExact(int *Avec, int *Bvec, int tnum, 
                         int ncomb, int getDist, 
                         double *results);
void contingency2xtSim(int *Avec, int *Bvec, int tnum, int nsim, 
		            int getDist, double *results);
void QNraw(double *QN, int k, double *rx, int *ns);
void QNinitvals(int kk, double *rrvec, int *nns, 
			int *iivec, double *tteststat, double *ppval, 
                int nnsum, int ggetQNdist, double *QQNvec, int ddimst); 
void QNinitvals1(int kk, double *rrvec, int *nns);
// int QNexact(int now, int *position, int m);
void QNpvalue(double *pval, int Nsim, int k, double *rx, 
              int *ns, int useExact, int getQNdist, 
              double ncomb, double *QNobs, double *QNvec);
void SteelTestStat(double *Steel, int k, double *rx, int *ns);
              // int alt, double *mu, double *tau);
void Steelinitvals(int kk, double *rrvec, int *nns, 
		        int *iivec, double *tteststat, double *ppval, 
                   int nnsum, int aalt, double *mmu, 
                   double *ttau, int ggetSteeldist, 
                   double *SSteelvec, int ddimst);
void Steelinitvals1(int kk, double *xxvec, int *nns, 
		         int aalt, double *mmu, double *ttau);
int Steelexact(int now, int *position, int m);
void Steelpvalue(double *pval, int Nsim, int k, double *rx, 
                 int *ns, int useExact, int getSteeldist, 
		      double ncomb, int alt, double *mu, 
			 double *tau, double *Steelobs, 
                 double *Steelvec);
void SteelVec(int Nsim, int k, double *rx, int *ns,
		int useExact, double *MannWhitneyStats);
void SteelrawVec(double *SteelVec, int k, double *rx, int *ns);
                 
int SteelexactVec(int now, int *position, int m);
void SteelinitvalsVec(int k, double *rvec, int *ns, int *ivec, 
				 int nsum, double *MMannWhitneyStats);
void insertxp(double xnew, double pnew, int k, int  *Lt, 
              double *values, double *probs);
void convaddtotable(double xnew, double pnew, int *Lt, int M, 
			    double *values, double *probs);
void conv(double *x1, double *p1, int *n1, double *x2, 
          double *p2, int *n2, double *x, double *p, int *n);
void convvec(double *x1, int *n1, double *x2, int *n2, 
             double *x, int *n);

void Harding(int k, int L1, int *nn, int *nvec, double *freq);

void JTinitvals(int kk, double *xxvec, int *nns, 
			int *iivec, double *tteststat, double *ppval, 
                int nnsum, int ggetJTdist, double *JJTvec, int ddimst);
void JTinitvals1(int kk, double *xxvec, int *nns);
void JTTestStat(double *JT, int k, double *rx, int *ns);
void JTpvalue(double *pval, int Nsim, int k, double *rx, 
              int *ns, int useExact, int getJTdist, 
              double ncomb, double *JTobs, double *JTvec);
void randPerm(int nsum, double *rx, double *rc, int *ns);
