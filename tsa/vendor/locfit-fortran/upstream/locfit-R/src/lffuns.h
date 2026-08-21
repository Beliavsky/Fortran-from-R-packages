/*
 *   Copyright (c) 1998-2001 Lucent Technologies.
 *   See README file for details.
 *
 *
 *
 *   Function definitions for Locfit.
 */

/* FILES IN THE src DIRECTORY */

/* lf_adap.c */
double adcri(double lk, double t0, double t2, double pen);
double mmse(lfdata *lfd, smpar *sp, deriv *dv, design *des);
int ainitband(lfdata *lfd, smpar *sp, deriv *dv, design *des);
double aband2(lfdata *lfd, smpar *sp, deriv *dv, design *des, double h0);
double aband3(lfdata *lfd, smpar *sp, deriv *dv, design *des, double h0);
extern int alocfit(lfdata *lfd, smpar *sp, deriv *dv, design *des);

/* band.c */
/* extern void band(), kdeselect(), kdecri(); */
int procvbind(design *des, lfit *lf, int v);
double bcri(double h, int c, int cri);
void bsel2(double h0, double g0, double ifact, int c, int cri);
void bsel3(double h0, double g0, double ifact, int c, int cri);
void bselect(lfit *lf, design *des, int c, int cri, double pn);
double compsda(double *x, double h, int n);
double widthsj(double *x, double lambda, int n);
extern void kdecri(double *x, double h, double *res, double c, int k, int ker, int n);
double esolve(double *x, int j, double h0, double h1, int k, double c, int ker, int n);
extern void kdeselect(double *band, double *x, Sint *ind, double h0, double h1, int *meth, int nm, int ker, int n);

/* density.c */
void prresp(double *coef, double *resp, int p);
int mif(double *u, int d, double *resp, double *M);
int multint(double *t, double *resp1, double *resp2, double *cf, double h);
int mlinint(double *t, double *resp1, double *resp2, double *cf, double h);
extern void prodintresp(double *resp, double prod_wk[MXDIM][2*MXDEG+1], int dim, int deg, int p);
extern int prodint(double *t, double *resp, double *resp2, double *coef, double h);
int gausint(double *t, double *resp, double *C, double *cf, double h, double *sca);
extern int likeden(double *coef, double *lk0, double *f1, double *A);
int inre(double *x, double *bound, int d);
int setintlimits(lfdata *lfd, double *x, double h, int *ang, int *lset);
int selectintmeth(int itype, int lset, int ang);
extern int densinit(lfdata *lfd, design *des, smpar *sp, double *cf);

extern int fact[];
extern int de_mint, de_itype, de_renorm;

/* dens_haz.c */
int haz_sph_int(double *dfx, double *cf, double h, double *r1);
int hazint_sph(double *t, double *resp, double *r1, double *cf, double h);
int hazint_prod(double *t, double *resp, double *x, double *cf, double h);
extern int hazint(double *t, double *resp, double *resp1, double *cf, double h);
extern void haz_init(lfdata *lfd, design *des, smpar *sp, double *il);

/* dens_int.c */
extern void lforder(Sint *ind, double *x, int l, int r);
double estdiv(double x0, double x1, double f0, double f1, double d0, double d1, int lin);
extern double dens_integrate(lfit *lf, design *des, int z);
extern void dens_renorm(lfit *lf, design *des);
extern void dens_lscv(design *des, lfit *lf);

/* ev_atree.c */
extern void atree_guessnv(evstruc *evs, int *nvm, int *ncm, int *vc, int d, double alp);
int atree_split(lfit *lf, Sint *ce, double *le, double *ll, double *ur);
extern void atree_grow(design *des, lfit *lf, Sint *ce, Sint *ct, Sint *term, double *ll, double *ur);
extern void atree_start(design *des, lfit *lf);
extern double atree_int(lfit *lf, double *x, int what);

/* ev_interp.c */
extern double linear_interp(double h, double d, double f0, double f1);
extern void hermite2(double x, double z, double *phi);
extern double cubic_interp(double h, double f0, double f1, double d0, double d1);
extern double cubintd(double h, double f0, double f1, double d0, double d1);
extern double rectcell_interp(double *x, double vv[64][64], double *ll, double *ur, int d, int nc);
extern int exvval(fitpt *fp, double *vv, int nv, int d, int what, int z);
extern void exvvalpv(double *vv, double *vl, double *vr, int d, int k, double dl, int nc);
double grid_int(fitpt *fp, evstruc *evs, double *x, int what);
double fitp_int(fitpt *fp, double *x, int what, int i);
double xbar_int(fitpt *fp, double *x, int what);
extern double dointpoint(lfit *lf, double *x, int what, int ev, int j);

/* ev_kdtre.c */
extern void kdtre_guessnv(evstruc *evs, int *nvm, int *ncm, int *vc, int n, int d, double alp);
int ksmall(int l, int r, int m, double *x, Sint *pi);
int terminal(lfit *lf, int p, Sint *pi, int fc, int d, int *m, double *split_val);
extern void kdtre_start(design *des, lfit *lf);
void newcell(int *nv, int vc, double *xev, int d, int k, double split_val, Sint *cpar, Sint *clef, Sint *crig);
double blend(fitpt *fp, evstruc *evs, double s, double *x, double *ll, double *ur, int j, int nt, int *t, int what);
extern double kdtre_int(fitpt *fp, evstruc *evs, double *x, int what);

/* ev_sphere.c */
extern void sphere_guessnv(int *nvm, int *ncm, int *vc, int *mg);
extern void sphere_start(design *des, lfit *lf);
extern double sphere_int(lfit *lf, double *x, int what);

/* ev_main.c */
extern void lfit_alloc(lfit *lf);
extern int lfit_reqd(int d, int nvm, int ncm, int geth);
extern int lfit_reqi(int nvm, int ncm, int vc);
extern void trchck(lfit *lf, int nvm, int ncm, int vc);
extern void data_guessnv(int *nvm, int *ncm, int *vc, int n);
extern void dataf(design *des, lfit *lf);
extern void xbar_guessnv(int *nvm, int *ncm, int *vc);
extern void xbarf(design *des, lfit *lf);
extern void preset(design *des, lfit *lf);
extern void crossf(design *des, lfit *lf);
extern void gridf(design *des, lfit *lf);
extern int findpt(fitpt *fp, evstruc *evs, int i0, int i1);
extern int newsplit(design *des, lfit *lf, int i0, int i1, int pv);

/* ev_trian.c */
void solve(double *A, double *b, int d);
extern void triang_guessnv(int *nvm, int *ncm, int *vc, int d, int mk);
int triang_split(lfit *lf, Sint *ce, double *le);
void resort(int *pv, double *xev, int *dig);
extern void triang_grow(design *des, lfit *lf, Sint *ce, Sint *ct, Sint *term);
void triang_descend(lfit *tr, double *xa, Sint *ce);
void covrofdata(lfdata *lfd, double *V, double *mn);
int intri(double *x, Sint *w, double *xev, double *xa, int d);
extern void triang_start(design *des, lfit *lf);
double triang_cubicint(double *v, double *vv, Sint *w, int d, int nc, double *xxa);
double triang_clotoch(double *xev, double *vv, Sint *ce, int p, double *xxa);
int triang_getvertexvals(fitpt *fp, evstruc *evs, double *vv, int i, int what);
extern double triang_int(lfit *lf, double *x, int what);

/* family.c */
extern int defaultlink(int link, int family);
extern int validlinks(int link, int family);
int famdens(double mean, double th, int link, double *res, int cens, double w);
int famgaus(double y, double mean, double th, int link, double *res, int cens, double w);
int famrobu(double y, double mean, double th, int link, double *res, int cens, double w, double rs);
int famcauc(double y, double p, double th, int link, double *res, int cens, double w, double rs);
int famrbin(double y, double p, double th, int link, double *res, int cens, double w);
int fambino(double y, double p, double th, int link, double *res, int cens, double w);
int fampois(double y, double mean, double th, int link, double *res, int cens, double w);
int famgamm(double y, double mean, double th, int link, double *res, int cens, double w);
int famgeom(double y, double mean, double th, int link, double *res, int cens, double w);
int famweib(double y, double mean, double th, int link, double *res, int cens, double w);
int famcirc(double y, double mean, double th, int link, double *res, int cens, double w);
extern int links(double th, double y, int fam, int link, double *res, int c, double w, double rs);
extern int stdlinks(double *res, lfdata *lfd, smpar *sp, int i, double th, double rs);
extern double b2(double th, int tg, double w);
extern double b3(double th, int tg, double w);
extern double b4(double th, int tg, double w);
extern double lf_link(double y, int lin);
extern double invlink(double th, int lin);

/* fitted.c */
double resid(double y, double w, double th, int fam, int ty, double *res);
double studentize(double res, double inl, double var, int ty, double *link);
extern void fitted(lfit *lf, double *fit, int what, int cv, int st, int ty);

/* frend.c */
extern void ressummd(lfit *lf);
void ressumm(lfit *lf, design *des);
extern double rss(lfit *lf, design *des, double *df);

/* lf_dercor.c */
void dercor(lfdata *lfd, smpar *sp, design *des, double *coef);

/* lf_fitfun.c */
extern int calcp(smpar *sp, int d);
extern int coefnumber(deriv *dv, int kt, int d, int deg);
extern void makecfn(smpar *sp, design *des, deriv *dv, int d);
void fitfunangl(double dx, double *ff, double sca, int cd, int deg);
extern void fitfun(lfdata *lfd, smpar *sp, double *x, double *t, double *f, deriv *dv);
extern void designmatrix(lfdata *lfd, smpar *sp, design *des);

/* lf_nbhd.c */

extern double rho(double *x, double *sc, int d, int kt, int *sty);
extern double kordstat(double *x, int k, int n, Sint *ind);
int inlim(lfdata *lfd, int i);
double compbandwid(double *di, Sint *ind, double *x, int n, int d, int nn, double fxh);
extern void nbhd1(lfdata *lfd, smpar *sp, design *des, int k);
void nbhd_zeon(lfdata *lfd, design *des);
void nbhd(lfdata *lfd, design *des, int nn, int redo, smpar *sp);

/* lf_robust.c */
extern double median(double *x, int n);
double nrobustscale(lfdata *lfd, smpar *sp, design *des, double rs);
double robustscale(lfdata *lfd, smpar *sp, design *des);
double update_rs(double x);
extern void lf_robust(lfdata *lfd, smpar *sp, design *des, int mxit);

/* lfstr.c */
int ct_match(char *z1, char *z2);
int pmatch(char *z, char **strings, int *vals, int n, int def);
extern int lffamily(char *z);
extern int lfkernel(char *z);
extern int lfketype(char *z);
extern int lflink(char *z);
extern int deitype(char *z);
extern int lfacri(char *z);
extern int lfevstr(char *z);
extern int restyp(char *z);
extern int ppwhat(char *z);

/* lf_vari.c */
void vmat(lfdata *lfd, smpar *sp, design *des, double *M12, double *M2);
extern void lf_vcov(lfdata *lfd, smpar *sp, design *des);
extern void comp_vari(lfdata *lfd, smpar *sp, design *des, double *tr, double *t0);
extern void local_df(lfdata *lfd, smpar *sp, design *des, double *tr);

/* locfit.c */
extern void lfdata_init(lfdata *lfd);
extern void smpar_init(smpar *sp, lfdata *lfd);
extern void deriv_init(deriv *dv);
extern void des_init(design *des, int n, int p);
void deschk(design *des, int n, int p);
int likereg(double *coef, double *lk0, double *f1, double *Z);
int robustinit(lfdata *lfd, design *des);
int circinit(lfdata *lfd, design *des);
int reginit(lfdata *lfd, design *des);
int lfinit(lfdata *lfd, smpar *sp, design *des);
extern void lfiter(design *des, int maxit);
int use_robust_scale(int tg);
extern int locfit(lfdata *lfd, design *des, smpar *sp, int noit, int nb, int cv);
extern int des_reqd(int n, int p);
extern int des_reqi(int n, int p);

extern int lf_maxit, lf_debug;

/* math.c */
extern double lf_exp(double x);
extern double lferfc(double x);
extern double lferf(double x);
extern double lflgamma(double x);
extern double lfdaws(double x);
extern double ptail(double x);
extern double logit(double x);
extern double expit(double x);
extern int factorial(int n);

/* minmax.c */
extern double ipower(double x, int n);
double setmmwt(design *des, double *a, double gam);
int mmsums(double *coef, double *f, double *z, jacobian *J);
double updatesd(design *des, double *z, int p, double *a, double *a0, double sw0, double gam);
int mm_initial(design *des, double *z, int p, double *coef);
void mmax(double *coef, double *old_coef, double *f1, double *delta, jacobian *J, int p, int maxit, double tol, int *err);
double findab(double gam);
double weightmm(double *coef, double di, double *ff, double gam);
extern double minmax(lfdata *lfd, design *des, smpar *sp);

/* dens_odi.c */
int exbctay(double b, double c, int n, double *z);
double explinjtay(double l0, double l1, int j, double *cf);
void explint1(double l0, double l1, double *cf, double *I, int p);
void explintyl(double l0, double l1, double *cf, double *I, int p);
void solvetrid(double *X, double *y, int m);
void initi0i1(double *I, double *cf, double y0, double y1, double l0, double l1);
void explinsid(double l0, double l1, double *cf, double *I, int p);
void explinbkr(double l0, double l1, double *cf, double *I, int p);
void explinfbk0(double l0, double l1, double *cf, double *I, int p);
void explinfbk(double l0, double l1, double *cf, double *I, int p);
void recent(double *I, double *resp, double *wt, int p, int s, double x);
extern void recurint(double l0, double l2, double *cf, double *resp, int p, int ker);
int onedexpl(double *cf, int deg, double *resp);
int onedgaus(double *cf, int deg, double *resp);
extern int onedint(smpar *sp, double *cf, double l0, double l1, double *resp);

/* pcomp.c */
extern int noparcomp(smpar *sp, int geth);
extern int pc_reqd(int d, int p);
extern void pcchk(paramcomp *pc, int d, int p, int lc);
extern void compparcomp(design *des, lfdata *lfd, smpar *sp, paramcomp *pc, int geth, int nopc);
extern void subparcomp(design *des, lfit *lf, double *coef);
extern void subparcomp2(design *des, lfit *lf, double *vr, double *il);
extern double addparcomp(lfit *lf, double *x, int c);

/* preplot.c */
void predptall(lfit *lf, double *x, int what, int ev, int i);
void prepvector(lfit *lf, double **x, int n, int what);
void prepfitp(lfit *lf, int what);
void prepgrid(lfit *lf, double **x, Sint *mg, int n, int what);
extern void preplot(lfit *lf, double **x, double *f, double *se, char band, Sint *mg, 
             int where, int what);
/* extern void cpreplot();*/
/* extern int setpppoints(); */

/* procv.c */
double vocri(double lk, double t0, double t2, double pen);
extern int procvraw(design *des, lfit *lf, int v);
void set_default_like(fitpt *fp, int v);
extern int procv(design *des, lfit *lf, int v);
double intvo(design *des, lfit *lf, double *c0, double *c1, double a, int p, double t0, double t20, double t21);
extern int procvvord(design *des, lfit *lf, int v);
extern int procvhatm(design *des, lfit *lf, int v);

/* resid.c */
/* extern double resid(); */

/* scb.c */
double covar_par(lfit *lf, design *des, double x1, double x2);
void cumulant(lfit *lf, design *des, double sd);
double q2(double u);
double p2(double u);
double gldn_like(double a);
void get_gldn(fitpt *fp, design *des, double *lo, double *hi, int v);
int procvscb2(design *des, lfit *lf, int v);
extern void scb(design *des, lfit *lf);
/* extern void scb(), cscbsim(); */

/* scb_iface.c */
int scbfitter(double *x, double *l, int reqd);
extern int constants(design *des, lfit *lf);

/* simul.c */
void goldensec(double (*f)(), design *des, lfit *tr, double eps, double *xm, double *ym, int meth);
double dnk(double x, int k);
double locai(double h, design *des, lfit *lf);
double loccp(double h, design *des, lfit *lf, int m);
double cp(design *des, lfit *lf, int meth);
double gkk(design *des, lfit *lf);
double rsw(design *des, lfit *lf);
extern void rband(design *des, lfit *lf, double *hhat, int *meth, int nmeth);
/* extern void liksim(), scbsim(), scbmax(), regband(), rband(); */

/* startlf.c */
void evstruc_init(evstruc *evs);
void fitpt_init(fitpt *fp);
extern void lfit_init(lfit *lf);
void fitdefault(lfit *lf);
extern void set_flim(lfdata *lfd, evstruc *evs);
double vecsum(double *v, int n);
double vvari(double *v, int n);
extern void set_scales(lfdata *lfd);
extern void startlf(design *des, lfit *lf, int (*vfun)(), int nopc);
/*
extern void set_flim(), set_scales(), startlf(), lfit_init();
extern void fitoptions(), clocfit(), endfit();
extern int nofit(); */

/* strings.c */
/*
extern int stm(), pmatch(), matchlf(), matchrt(), checkltor(), checkrtol();
extern void strip();
*/

/* lf_wdiag.c */
void nnresproj(lfdata *lfd, smpar *sp, design *des, double *u, int m, int p);
void wdexpand(double *l, int n, Sint *ind, int m);
extern int wdiagp(lfdata *lfd, smpar *sp, design *des, double *lx, paramcomp *pc, deriv *dv, int deg, int ty, int exp);
extern int wdiag(lfdata *lfd, smpar *sp, design *des, double *lx, deriv *dv, int deg, int ty, int exp);

/* weight.c */
extern double W(double u, int ker);
extern int iscompact(int ker);
double weightprod(lfdata *lfd, double *u, double h, int ker);
double weightsph(lfdata *lfd, double *u, double h, int ker, int hasdi, double di);
extern double weight(lfdata *lfd, smpar *sp, double *x, double *t, double h, int hasdi, double di);
double sgn(double x);
double WdW(double u, int ker);
extern double weightd(double u, double sc, int d, int ker, int kt, double h, int sty, double di);
double weightdd(double *u, double *sc, int d, int ker, int kt, double h, int *sty, double di, int i0, int i1);
extern double Wd(double u, int ker);
extern double Wdd(double u, int ker);
extern double wint(int d, int *j, int nj, int ker);
extern int wtaylor(double *f, double x, int ker);
extern double Wconv(double v, int ker);
extern double Wconv1(double v, int ker);
extern double Wconv4(double v, int ker);
extern double Wconv5(double v, int ker);
extern double Wconv6(double v, int ker);
extern double Wikk(int ker, int deg);
