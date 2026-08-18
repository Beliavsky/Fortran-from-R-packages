#include <R.h>
#include <Rinternals.h>


#define xx(i,j,l) (X[(i)+dimx[0]*((j)+dimx[1]*(l))])
#define yy(i,j,l) (Y[(i)+dimy[0]*((j)+dimy[1]*(l))])
#define ee(i,j,l) (E[(i)+dime[0]*((j)+dime[1]*(l))])

static R_NativePrimitiveArgType tensoramulhelper_t[] = {
  /* int *dimx, int *dimy, int *dime,  double *X, double *Y, double *E*/
  INTSXP,    INTSXP,    INTSXP,     REALSXP,   REALSXP,   REALSXP
};

/**
The function performs many matrix multiplications in parallel.
*/          
extern void tensoramulhelper(int *dimx,int *dimy,int *dime,
		      double *X,double *Y,double *E) {
    int i,j,k,l;
    double tmp;
    if(dimx[1]!=dimy[0]||dimx[2]!=dimy[2]||dimx[0]!=dime[0]||
       dimy[1]!=dime[1]||dimx[2]!=dime[2])
	error("multensorhelper: Dimension missmatch");
    for(l=0;l<dime[2];l++)
	for(i=0;i<dimx[0];i++)
	    for(j=0;j<dimy[1];j++) {
		tmp=0.0;
		for(k=0;k<dimx[1];k++) 
		    tmp += xx(i,k,l)*yy(k,j,l);
		ee(i,j,l)=tmp;
	    }
}
#define realpart r
#define imaginaripart i

#define xxX(i,j,l) (X[((i)+dimx[0]*((j)+dimx[1]*(l)))].realpart)
#define yyX(i,j,l) (Y[((i)+dimy[0]*((j)+dimy[1]*(l)))].realpart)
#define eeX(i,j,l) (E[((i)+dime[0]*((j)+dime[1]*(l)))].realpart)
#define xxY(i,j,l) (X[((i)+dimx[0]*((j)+dimx[1]*(l)))].imaginaripart)
#define yyY(i,j,l) (Y[((i)+dimy[0]*((j)+dimy[1]*(l)))].imaginaripart)
#define eeY(i,j,l) (E[((i)+dime[0]*((j)+dime[1]*(l)))].imaginaripart)


static R_NativePrimitiveArgType tensoraCmulhelper_t[] = {
  /* int *dimx, int *dimy, int *dime, Rcomplex *X, Rcomplex *Y, Rcomplex *E*/
  INTSXP,    INTSXP,    INTSXP,    CPLXSXP,    CPLXSXP,   CPLXSXP
};

extern void tensoraCmulhelper(int *dimx,int *dimy,int *dime,
			      Rcomplex *X,
			      Rcomplex *Y,
			      Rcomplex *E) {
    int i,j,k,l;
    double tmpX,tmpY;
    if(dimx[1]!=dimy[0]||dimx[2]!=dimy[2]||dimx[0]!=dime[0]||
       dimy[1]!=dime[1]||dimx[2]!=dime[2])
	error("multensorhelper: Dimension missmatch");
    for(l=0;l<dime[2];l++)
	for(i=0;i<dimx[0];i++)
	    for(j=0;j<dimy[1];j++) {
		tmpX=0.0;
		tmpY=0.0;
		for(k=0;k<dimx[1];k++) { 
		    tmpX += xxX(i,k,l)*yyX(k,j,l)-xxY(i,k,l)*yyY(k,j,l);
		    tmpY += xxX(i,k,l)*yyY(k,j,l)+xxY(i,k,l)*yyX(k,j,l);
		}
		eeX(i,j,l)=tmpX;
		eeY(i,j,l)=tmpY;
	    }
}


static R_CMethodDef cMethods[] = {
  {"tensoramulhelper", (DL_FUNC) &tensoramulhelper, 6, tensoramulhelper_t},
  {"tensoraCmulhelper", (DL_FUNC) &tensoraCmulhelper, 6, tensoraCmulhelper_t},
  {NULL, NULL, 0}
};


void R_init_tensora(DllInfo *info)
{
  R_registerRoutines(info, cMethods, NULL, NULL, NULL);
  R_useDynamicSymbols(info, FALSE);
  R_forceSymbols(info, TRUE);
}
