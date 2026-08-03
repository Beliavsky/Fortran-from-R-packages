#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

// Declare your C functions here
SEXP calcUpfrontTest
(SEXP baseDate_input,  /* (I) Value date  for zero curve       */
 SEXP types, /* "MMMMMSSSSSSSSS"*/
 SEXP rates, /* rates[14] = {1e-9, 1e-9, 1e-9, 1e-9, 1e-9, 1e-9, 1e-9,
 1e-9, 1e-9, 1e-9, 1e-9, 1e-9, 1e-9, 1e-9};/\* (I)
 Array of swap rates *\/ */
 SEXP expiries,
 SEXP mmDCC,          /* (I) DCC of MM instruments            */
 
 SEXP fixedSwapFreq,   /* (I) Fixed leg freqency/interval               */
 SEXP floatSwapFreq,   /* (I) Floating leg freqency/interval            */
 SEXP fixedSwapDCC,    /* (I) DCC of fixed leg                 */
 SEXP floatSwapDCC,    /* (I) DCC of floating leg              */
 SEXP badDayConvZC, //'M'  badDayConv for zero curve
 SEXP holidays,//'None'
 
 // input for upfront charge calculation
 SEXP todayDate_input, /*today: T (Where T = trade date)*/
   SEXP valueDate_input, /* value date: T+3 Business Days*/
   SEXP benchmarkDate_input,/* start date of benchmark CDS for internal
 ** clean spread bootstrapping;
 ** accrual Begin Date  */
 SEXP startDate_input,/* Accrual Begin Date */
 SEXP endDate_input,/*  Maturity (Fixed) */
 SEXP stepinDate_input,  /* T + 1*/
 
 SEXP dccCDS, 			/* accruedDcc */
 SEXP ivlCDS,
 SEXP stubCDS,
 SEXP badDayConvCDS,
 SEXP calendar,
 
 SEXP parSpread,
 SEXP couponRate,
 SEXP recoveryRate,
 SEXP isPriceClean_input,
 SEXP payAccruedOnDefault_input,
 SEXP notional);

SEXP calcCdsoneSpread
(// variables for the zero curve
    SEXP baseDate_input,  /* (I) Value date  for zero curve       */
    SEXP types,           /*"MMMMMSSSSSSSSS" */
    
    SEXP rates, /* rates[14] = {1e-9, 1e-9, 1e-9, 1e-9, 1e-9, 1e-9, 1e-9,
 1e-9, 1e-9, 1e-9, 1e-9, 1e-9, 1e-9, 1e-9};/\* (I)
 Array of swap rates *\/ */
 
 SEXP expiries, 
 SEXP mmDCC, /* (I) DCC of MM instruments            */
 
 SEXP fixedSwapFreq,   /* (I) Fixed leg freqency               */
 SEXP floatSwapFreq,   /* (I) Floating leg freqency            */
 SEXP fixedSwapDCC,    /* (I) DCC of fixed leg                 */
 SEXP floatSwapDCC,    /* (I) DCC of floating leg              */
 SEXP badDayConvZC, //'M'  badDayConv for zero curve
 SEXP holidays,//'None'
 
 SEXP todayDate_input, /*today: T (Where T = trade date)*/
   SEXP valueDate_input, /* value date: T+3 Business Days*/
   SEXP benchmarkDate_input,  /* start date of benchmark CDS for
 ** internal clean spread bootstrapping;
 ** accrual Begin Date  */
 SEXP startDate_input,  /* Accrual Begin Date */
 SEXP endDate_input,  /*  Maturity (Fixed) */
 SEXP stepinDate_input,  /* T + 1*/
 
 SEXP couponRate_input, 	/* Fixed Coupon Amt */
 SEXP payAccruedOnDefault_input, /* TRUE in new contract */
 
 SEXP dccCDS, //accrueDCC_input,	/* ACT/360 */
 SEXP dateInterval,		  /* Q - 3 months in new contract */
 SEXP stubType, 		/* F/S */
 SEXP badDayConv_input, 	/* (F) Following */ 
 SEXP calendar_input,	/*  None (no holiday calendar) in new contract */
 
 SEXP upfrontCharge_input,
 SEXP recoveryRate_input,
 SEXP payAccruedAtStart_input	/* (True/False), True: Clean Upfront supplied */
);

static const R_CallMethodDef CallEntries[] = {
  {"calcUpfrontTest", (DL_FUNC) &calcUpfrontTest, 28},
  {"calcCdsoneSpread", (DL_FUNC) &calcCdsoneSpread, 27},
  {NULL, NULL, 0}
};

void R_init_creditr(DllInfo *dll) {
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
