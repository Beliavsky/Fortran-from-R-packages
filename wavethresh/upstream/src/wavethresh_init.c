#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* .C calls */

/* 1-10 */

extern void accessDwpst(double *coefvec, int *lansvec, int *nlev, int *avixstart, int *primaryindex, int *nwppkt, int *pklength, int *level, double *weave, int *lweave, int *error);

extern void av_basisWRAP(double *wst, double *wstC, int *LengthData, int *level, double *H, int *LengthH, double *answer, int *error);

extern void c2to4(int *l, int *a);

extern void Ccthrcalcodds(long *pnd, double *dr, double *di, double *VVec, double *SigVec, double *pp, double *ans, double *odds);

extern void Ccthrnegloglik(double *parvec, double *SigVec, double *di, double *dr, long *pnd, double *pans);

extern void Cmnv(double *wst, double *wstC, int *LengthData, int *nlevels, int *upperctrl, double *upperl, int *firstl, int *verbose, int *error);

extern void comAB_WRAP(double *wstR, double *wstI, double *wstCR, double *wstCI, int *LengthData, int *level, double *HR, double *HI, double *GR, double *GI, int *LengthH, double *answerR, double *answerI, int *error);

extern void computec(int *n,double *c,int *gridn,double *Gmatrix,int *Gindex,
              double *H, int *LengthH, int *bc, int *error);

extern void comwd(double *CR, double *CI, int *LengthC, double *DR, double *DI, int *LengthD, double *HR, double *HI, double *GR, double *GI, int *LengthH, int *levels, int *firstC, int *lastC, int *offsetC, int *firstD, int *lastD, int *offsetD, int *type, int *bc, int *error);

extern void comwr(double *CR, double *CI, int *LengthC, double *DR, double *DI, int *LengthD, double *HR, double *HI, double *GR, double *GI, int *LengthH, int *levels, int *firstC, int *lastC, int *offsetC, int *firstD, int *lastD, int *offsetD, int *type, int *bc, int *error);

/* 11-20 */

extern void comwst(double *CaR, double *CaI, double *DataR, double *DataI, int *LengthData, int *levels, double *HR, double *HI, double *GR, double *GI, int *LengthH, int *error);

extern void conbarL(double *c_in, int *LengthCin, int *firstCin, double *d_in, int *LengthDin, int *firstDin, double *H, int *LengthH, double *c_out, int *LengthCout, int *firstCout, int *lastCout, int *type, int *bc);

extern void Cpostmean(long *pnd, double *dr, double *pi, double *VVec, double *SigVec, double *w, double *ansr, double *ansi);

extern void Crsswav(double *noisy, int *nnoisy, double *value, double *C, double *D, int *LengthD, double *H, int *LengthH, int *levels, int *firstC, int *lastC, int *offsetC, int *firstD, int *lastD, int *offsetD, int *ntt, int *ll, int *bc, double *ssq, int *error);

extern void CScalFn(double *v, double *ans, int *res, double *H, int *lengthH);

extern void Cthreshold(double *D, int *LengthD, int *firstD, int *lastD, int *offsetD, int *Dlevels, int *ntt, double *value, int *levels, int *qlevels, int *bc, int *error);

extern void CWavDE(double *x, int *n, double *minx, double *maxx, int *Jmax, double *threshold, double *xout, double *fout, int *nout, double *PrimRes, double *SFx, double *SFy, int *lengthSF, double *WVx, double *WVy, int *lengthWV, int *kmin, int *kmax, int *kminW, int *kmaxW, double *xminW, double *xmaxW, double *phiLH, double *phiRH, double *psiLH, double *psiRH, int *verbose, int *error);

extern void CWaveletCV(double *noisy, int *nnoisy, double *UniversalThresh, double *C, double *D, int *LengthD, double *H, int *LengthH, int *levels, int *firstC, int *lastC, int *offsetC, int *firstD, int *lastD, int *offsetD, int *ntt, int *ll, int *bc, double *tol, int *maxits, double *xvthresh, int *interptype, int *error);

extern void dec(double* data, int* size, int* filternumber, int* minscale, int* precond, int* filterhistory);

extern void FullWaveletCV(double *noisy, int *nnoisy, double *UniversalThresh, double *H, int *LengthH, int *ntt, int *ll, double *tol, double *xvthresh, int *error);

/* 21-30 */

extern void getARRel(double *Carray, int *size, int *level, double *GHH, double *HGH, double *GGH, double *HHG, double *GHG, double *HGG, double *GGG);

extern void getpacketwst2D(double *am, int *D1, int *D12, int *maxlevel, int *level, int *index, int *type, double *out, int *sl);

extern void ixtoco(int *level, int *maxlevel, int *index, int *x, int *y);

extern void makegrid(double *x,double *y,int *n,double *gridx,double *gridy,
               int *gridn, double *G, int *Gindex);

extern void multiwd(double *C, int *lengthc, double *D, int *lengthd, int *nlevels, int *nphi, int *npsi, int *ndecim, double *H, double *G, int *NH, int *lowerc, int *upperc, int *offsetc, int *upperd, int *lowerd, int *offsetd, int *nbc);

extern void multiwr(double *C, int *lengthc, double *D, int *lengthd, int *nlevels, int *nphi, int *npsi, int *ndecim, double *H, double *G, int *NH, int *lowerc, int *upperc, int *offsetc, int *upperd, int *lowerd, int *offsetd, int *nbc, int *startlevel);

extern void PLDE2(double *C, double *p, double *filter, int *nf, int *prec, int *kmin, int *kmax, double *gx, double *gy, int *ng, double *philh, double *phirh, int *error);

extern void PsiJ(int *J, double *H, int *LengthH, double *tol, double *wout, int *lwout, int *rlvec, int *error);

extern void putarr(double *Carray, int *truesize, int *level, int *Iarrayix, double *Iarray);

extern void putpacketwst2D(double *am, int *D1, int *D12, int *maxlevel, int *level, int *index, int *type, double *in, int *sl);

/* 31-40 */

extern void rainmatPARENT(int *J, double *H, int *LengthH, double *fmat, double *tol, int *error);


extern void rainmatPARTIAL(int *J, int *donej, double *H, int *LengthH, double *fmat, double *tol, int *error);

extern void rec(double* data, int* size, int* filterhistory, int* currentscale, int* precond);

extern void SAvBasis(double *am, int *D1, int *D12, double *TheSmooth, int *levj, double *H, int *LengthH, int *error);

extern void SFDE5(double *x, int *nx, double *p, double *filter, int *nf, int *prec, double *chat, int *kmin, int *kmax, double *philh, double *phirh, int *error);

extern void SFDE6(double *x, int *nx, double *p, double *filter, int *nf, int *prec, double *chat, double *covar, int *kmin, int *kmax, double *philh, double *phirh, int *error);

extern void StoDCDS(double *C, int *Csize, int *firstCin, double *H, int *LengthH, int *LengthCout, int *firstCout, int *lastCout, int *LengthDout, int *firstDout, int *lastDout, double *ImCC, double *ImDD, int *bc, int *type, int *error);

extern void StoIDS(double *C, int *Csize, int *firstCin, double *H, int *LengthH, int *LengthCout, int *firstCout, int *lastCout, int *LengthDout, int *firstDout, int *lastDout, double *ImCC, double *ImCD, double *ImDC, double *ImDD, int *bc, int *type, int *error);

extern void StoIRS(double *ImCC, double *ImCD, double *ImDC, double *ImDD, int *LengthCin, int *firstCin, int *LengthDin, int *firstDin, double *H, int *LengthH, int *LengthCout, int *firstCout, int *lastCout, double *ImOut, int *bc, int *error);

extern void SWT2Dall(double *m, int *nm, double *am, int *J, double *H, int *LengthH, int *error);

/* 41-50 */

extern void tpwd(double *image, int *nrow, int *ncol, int *levr, int *levc, int *firstCr, int *lastCr, int *offsetCr, int *firstDr, int *lastDr, int *offsetDr, int *firstCc, int *lastCc, int *offsetCc, int *firstDc, int *lastDc, int *offsetDc, int *type, int *bc, double *H, int *LengthH, int *error); 

extern void tpwr(double *image, int *nrow, int *ncol, int *levr, int *levc, int *firstCr, int *lastCr, int *offsetCr, int *firstDr, int *lastDr, int *offsetDr, int *firstCc, int *lastCc, int *offsetCc, int *firstDc, int *lastDc, int *offsetDc, int *type, int *bc, double *H, int *LengthH, int *error); 

extern void wavedecomp(double *C, double *D, double *H, int *LengthH, int *levels, int *firstC, int *lastC, int *offsetC, int *firstD, int *lastD, int *offsetD, int *type, int *bc, int *error);

extern void wavedecomp_dh(double *C, double *D, double *H, int *LengthH, int *levels, int *firstC, int *lastC, int *offsetC, int *firstD, int *lastD, int *offsetD, int *type, int *bc, int *error);

extern void wavepackde(double *Data, int *LengthData, int *levels, double *H, int *LengthH);

extern void wavepackrecon(double *rdata, int *ldata, int *nrsteps, int *rvector, double *H, int *LengthH, int *error);

extern void wavepackst(double *Carray, double *Data, int *LengthData, int *levels, double *H, int *LengthH, int *error);

extern void waverecons(double *C, double *D, double *H, int *LengthH, int *levels, int *firstC, int *lastC, int *offsetC, int *firstD, int *lastD, int *offsetD, int *type, int *bc, int *error);

extern void waverecons_dh(double *C, double *D, double *H, int *LengthH, int *levels, int *firstC, int *lastC, int *offsetC, int *firstD, int *lastD, int *offsetD, int *type, int *bc, int *error);

extern void wd3D(double *Carray, int *size, double *H, int *LengthH, int *error);

/* 51-53 */

extern void wpCmnv(double *wp, int *LengthData, int *nlevels, int *upperctrl, double *upperl, int *firstl, int *verbose, int *error);

extern void wpst(double *ansvec, int *lansvec, int *nlev, int *finish_level, int *avixstart, double *H, int *LengthH, int *error);

extern void wr3D(double *Carray, int *truesize, double *H, int *LengthH, int *error);

static const R_CMethodDef CEntries[] = {
    {"accessDwpst",    (DL_FUNC) &accessDwpst,    11},
    {"av_basisWRAP",   (DL_FUNC) &av_basisWRAP,    8},
    {"c2to4",          (DL_FUNC) &c2to4,           2},
    {"Ccthrcalcodds",  (DL_FUNC) &Ccthrcalcodds,   8},
    {"Ccthrnegloglik", (DL_FUNC) &Ccthrnegloglik,  6},
    {"Cmnv",           (DL_FUNC) &Cmnv,            9},
    {"comAB_WRAP",     (DL_FUNC) &comAB_WRAP,     14},
    {"computec",       (DL_FUNC) &computec,        8},
    {"comwd",          (DL_FUNC) &comwd,          21},
    {"comwr",          (DL_FUNC) &comwr,          21},
    {"comwst",         (DL_FUNC) &comwst,         12},
    {"conbarL",        (DL_FUNC) &conbarL,        14},
    {"Cpostmean",      (DL_FUNC) &Cpostmean,       8},
    {"Crsswav",        (DL_FUNC) &Crsswav,        20},
    {"CScalFn",        (DL_FUNC) &CScalFn,         5},
    {"Cthreshold",     (DL_FUNC) &Cthreshold,     12},
    {"CWavDE",         (DL_FUNC) &CWavDE,         28},
    {"CWaveletCV",     (DL_FUNC) &CWaveletCV,     23},
    {"dec",            (DL_FUNC) &dec,             6},
    {"FullWaveletCV",  (DL_FUNC) &FullWaveletCV,  10},
    {"getARRel",       (DL_FUNC) &getARRel,       10},
    {"getpacketwst2D", (DL_FUNC) &getpacketwst2D,  9},
    {"ixtoco",         (DL_FUNC) &ixtoco,          5},
    {"makegrid",       (DL_FUNC) &makegrid,        8},
    {"multiwd",        (DL_FUNC) &multiwd,        18},
    {"multiwr",        (DL_FUNC) &multiwr,        19},
    {"PLDE2",          (DL_FUNC) &PLDE2,          13},
    {"PsiJ",           (DL_FUNC) &PsiJ,            8},
    {"putarr",         (DL_FUNC) &putarr,          5},
    {"putpacketwst2D", (DL_FUNC) &putpacketwst2D,  9},
    {"rainmatPARENT",  (DL_FUNC) &rainmatPARENT,   6},
    {"rainmatPARTIAL", (DL_FUNC) &rainmatPARTIAL,  7},
    {"rec",            (DL_FUNC) &rec,             5},
    {"SAvBasis",       (DL_FUNC) &SAvBasis,        8},
    {"SFDE5",          (DL_FUNC) &SFDE5,          12},
    {"SFDE6",          (DL_FUNC) &SFDE6,          13},
    {"StoDCDS",        (DL_FUNC) &StoDCDS,        16},
    {"StoIDS",         (DL_FUNC) &StoIDS,         18},
    {"StoIRS",         (DL_FUNC) &StoIRS,         16},
    {"SWT2Dall",       (DL_FUNC) &SWT2Dall,        7},
    {"tpwd",           (DL_FUNC) &tpwd,           22},
    {"tpwr",           (DL_FUNC) &tpwr,           22},
    {"wavedecomp",     (DL_FUNC) &wavedecomp,     14},
    {"wavedecomp_dh",  (DL_FUNC) &wavedecomp_dh,  14},
    {"wavepackde",     (DL_FUNC) &wavepackde,      5},
    {"wavepackrecon",  (DL_FUNC) &wavepackrecon,   7},
    {"wavepackst",     (DL_FUNC) &wavepackst,      7},
    {"waverecons",     (DL_FUNC) &waverecons,     14},
    {"waverecons_dh",  (DL_FUNC) &waverecons_dh,  14},
    {"wd3D",           (DL_FUNC) &wd3D,            5},
    {"wpCmnv",         (DL_FUNC) &wpCmnv,          8},
    {"wpst",           (DL_FUNC) &wpst,            8},
    {"wr3D",           (DL_FUNC) &wr3D,            5},
    {NULL, NULL, 0}
};

void R_init_wavethresh(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
