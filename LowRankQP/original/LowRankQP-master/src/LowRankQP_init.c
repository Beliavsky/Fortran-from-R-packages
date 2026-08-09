#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* .C calls */
extern void PrintMatrix(char *,double *, int *, int *);
extern double VectorAbsSum(double *, int *);
extern void VectorVectorCopy(double *, double *, int *);
extern void VectorVectorDivide(double *, double *, double *, int *);
extern void VectorVectorMult(double *, double *, double *, int *);
extern void VectorVectorMinus(double *, double *, double *, int *);
extern double VectorVectorDot(double *, double *, int *);
extern void MatrixMatrixDiagSolve(double *, double *, double *, int *, int *);
extern void MatrixVectorMult(double *, double *, int *, double *, double *, double *, int *, int *);
extern void MatrixConstantPlusDiag(double *, double *, int *);
extern void MatrixCholFactorize(double *, int *, int *);
extern void MatrixCholSolve(double *, int *, double *, int *, int *);
extern void MatrixLUFactorize(double *, int *, int *, int *);
extern void MatrixLUSolve(double *, int *, int *, double *, int *);
extern void MatrixMatrixPlusDiag(double *, double *, int *);
extern void MatrixMatrixCopy(double *, double *, int *, int *);
extern void MatrixMatrixMinus(double *, double *, double *, int *, int *);
extern void MatrixMatrixPlus(double *, double *, double *, int *, int *);
extern void MatrixMatrixPlusEquals(double *, double *, int *, int *);
extern void MatrixConstantSet(double *, double *, int *, int *);

extern void MatrixMatrixMult(double *, double *, int *, double *, int *, double *, double *, int *, 
    int *, int *, int *, int *, int *);
    
//extern void LRQPHeader(void);

extern void LRQPInitPoint(int *, int *, int *, double *, double *, double *, double *, double *, 
    double *, double *, double *, double *, double *, double *);

extern void LRQPCalcStats(int *, int *, int *, double *, double *, double *, double *, double *, 
    double *, double *, double *, double *, double *, double *, double *, double *, double *, 
    double *, double *, double *, double *, double *, double *, double *, double *, double *, 
    double *, double *, double *, double *);

extern void LRQPDisplay(int *, double *, double *, double *, double *, double *);
extern void LRQPSummary(int *, int *, int *, int *, int *, double *, double *, double *, double *, double *);
extern void PfcfSolve(int *, double *, double *, double *, int *);
extern void PfcfFactorize(int *, int *, double *, double *, double *, double *, double *, double *, double *);
extern void LRQPFactorize(int *, int *, int *, double *, double *, double *, int *, double *, 
    double *, double *, double *, double *, double *);
    
extern void LRQPSolve(int *, int *, int *, int *, double *, double *, double *, double *, double *, 
    int *, double *, double *, double *, double *);

extern void LRQPCalcDx(int *, int *, int *, int *, double *, double *, double *, double *, double *, 
    double *, double *, double *, double *, double *, double *, double *, double *, double *, 
    double *, double *, double *, double *, double *, double *, int *, double *, double *, double *, 
    double *, double *, double *, double *, double *, double *, double *, double *, double *, 
    double *, double *, double *, int *);
	
extern void LRQPStep(int *, int *, double *, double *, double *, double *, double *, double *, 
    double *, double *, double *, double *);

extern void LowRankQP(int *, int *, int *, int *, int *, int *, double *, double *, double *, 
    double *, double *, double *, double *, double *, double *, double *);

static const R_CMethodDef CEntries[] = {
    {"LowRankQP",              (DL_FUNC) &LowRankQP,	          16},
    {"LRQPCalcDx",             (DL_FUNC) &LRQPCalcDx,	       	  41},
    {"LRQPCalcStats",          (DL_FUNC) &LRQPCalcStats,	  30},     
    {"LRQPDisplay",            (DL_FUNC) &LRQPDisplay,	           6},
    {"LRQPFactorize",          (DL_FUNC) &LRQPFactorize,	  13},     
    {"LRQPInitPoint",          (DL_FUNC) &LRQPInitPoint,	  13},     
    {"LRQPSolve",              (DL_FUNC) &LRQPSolve,	          14},
    {"LRQPStep",               (DL_FUNC) &LRQPStep,	          12},
    {"LRQPSummary",            (DL_FUNC) &LRQPSummary,	          10},      
    {"MatrixConstantPlusDiag", (DL_FUNC) &MatrixConstantPlusDiag,  3},
    {"MatrixConstantSet",      (DL_FUNC) &MatrixConstantSet,       4},  
    {"MatrixCholFactorize",    (DL_FUNC) &MatrixCholFactorize,     3}, 
    {"MatrixCholSolve",        (DL_FUNC) &MatrixCholSolve,	   5},    
    {"MatrixLUFactorize",      (DL_FUNC) &MatrixLUFactorize,       4},
    {"MatrixLUSolve",	       (DL_FUNC) &MatrixLUSolve,	   5},  
    {"MatrixMatrixCopy",       (DL_FUNC) &MatrixMatrixCopy,        4},
    {"MatrixMatrixDiagSolve",  (DL_FUNC) &MatrixMatrixDiagSolve,   5},
    {"MatrixMatrixMinus",      (DL_FUNC) &MatrixMatrixMinus,       5},
    {"MatrixMatrixPlus",       (DL_FUNC) &MatrixMatrixPlus,        5},
    {"MatrixMatrixPlusDiag",   (DL_FUNC) &MatrixMatrixPlusDiag,    3},
    {"MatrixMatrixPlusEquals", (DL_FUNC) &MatrixMatrixPlusEquals,  4},
    {"MatrixMatrixMult",       (DL_FUNC) &MatrixMatrixMult,       13},
    {"MatrixVectorMult",       (DL_FUNC) &MatrixVectorMult,        8},
    {"PfcfFactorize",          (DL_FUNC) &PfcfFactorize,	   9},  
    {"PfcfSolve",	       (DL_FUNC) &PfcfSolve,	       	   5},
    {"PrintMatrix",	       (DL_FUNC) &PrintMatrix,	       	   4},
    {"VectorAbsSum",	       (DL_FUNC) &VectorAbsSum,	       	   2},
    {"VectorVectorCopy",       (DL_FUNC) &VectorVectorCopy,        3},
    {"VectorVectorDivide",     (DL_FUNC) &VectorVectorDivide,      4},
    {"VectorVectorDot",	       (DL_FUNC) &VectorVectorDot,	   3},  
    {"VectorVectorMinus",      (DL_FUNC) &VectorVectorMinus,       4},
    {"VectorVectorMult",       (DL_FUNC) &VectorVectorMult,        4},
    {NULL, NULL, 0}
};

void R_init_LowRankQP(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}



 	    
     
     
     
     
     
     
 	    
 	    
     
     
 	    
     
     
 	    
     
 	    
     
 	    
 	    
     
 	    
 	    
 	    
     
 	    
     
     
 	    
     
     
 	    
       

