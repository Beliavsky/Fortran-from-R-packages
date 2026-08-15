#define R_NO_REMAP
#include <R.h>
#include <Rinternals.h>
#include <math.h>
#include <stddef.h>
#include "f2c.h"
#include "R_ext/RS.h"

/* Table of constant values */

static doublereal c_b8 = 4.;
static doublereal c_b12 = 0.;

/*     ALGORITHM 614 COLLECTED ALGORITHMS FROM ACM. */
/*     ALGORITHM APPEARED IN ACM-TRANS. MATH. SOFTWARE, VOL.10, NO. 2, */
/*     JUN., 1984, P. 152-160. */
/*<       SUBROUTINE INTHP(A, B, D, F, M, P, EPS, INF, QUADR)                >*/
/* Subroutine */ SEXP inthp_sexp(SEXP a, SEXP b, SEXP d__, 
	SEXP f, SEXP m, SEXP p, SEXP eps, SEXP inf, 
	SEXP envir)
{
  double *A = REAL(a);
  double *B = REAL(b);
  double *P = REAL(p);
  double *EPS = REAL(eps);
  double *D__ = REAL(d__);
  int *M   = INTEGER(m);
  int *INF = INTEGER(inf);
  
    /* System generated locals */
    doublereal d__1;

    /* Builtin functions */
    double atan(doublereal), sqrt(doublereal), exp(doublereal), R_pow(
	    doublereal *, doublereal *), log(doublereal);

    /* Local variables */
    static doublereal alfa, exph, exph0, c__, h__;
    static integer i__, k, l, n;
    static doublereal s, t, u, v, w, c0, e1, h0, h1;
    static integer i1, l1, m1, m2, n1;
    static doublereal s1, v0, v1, v2, w1, w2, w3, w4, ba, pi, sr, sq2, cor, 
	    sum;
    static logical inf1, inf2;
    static doublereal eps3, sum1, sum2;

  //static double *tmp;
  // static char *mode[1], *ss[1];
  // static long length[1];
  // static void *args[1];
  static double zz[1];
  double quadr[1];
  SEXP ans;
  Rf_protect(ans = Rf_allocVector(REALSXP, 1));

/*        THIS SUBROUTINE COMPUTES INTEGRAL OF FUNCTIONS WHICH */
/*     MAY HAVE SINGULARITIES AT ONE OR BOTH END-POINTS OF AN */
/*     INTERVAL (A,B), SEE [1, 2]. IT CONTAINS FOUR DIFFERENT */
/*     QUADRATURE ROUTINES: ONE OVER A FINITE INTERVAL (A,B), */
/*     TWO OVER (A,+INFINITY), AND ONE OVER (-INFINITY,+INFINITY). */
/*     OF THE TWO FORMULAS OVER (A,+INFINITY), THE FIRST (INF=2 */
/*     BELOW) IS MORE SUITED TO NON-OSCILLATORY INTEGRANDS, WHILE */
/*     THE SECOND (INF=3) IS MORE SUITED TO OSCILLATORY INTEGRANDS. */
/*        THE USER SUPPLIES THE INTEGRAND FUNCTION, HE SPECIFIES THE */
/*     INTERVAL, AS WELL AS THE RELATIVE ERROR TO WHICH THE INTE- */
/*     GRAL IS TO BE EVALUATED. */
/*        THE FORMULAS ARE OPTIMAL IN CERTAIN HARDY SPACES H(P,DD), */
/*     SEE [1, 2]. HERE DD IS AN OPEN DOMAIN IN THE COMPLEX PLANE, */
/*     A AND B BELONG TO THE BOUNDARY OF DD AND H(P,DD), P.GT.1, IS */
/*     THE SET OF ALL ANALYTIC FUNCTONS IN DD WHOSE P-TH NORM DEFI- */
/*     NED AS IN [2] IS FINITE. */
/*        IF THE USER IS UNABLE TO SPECIFY THE PARAMETERS P AND D */
/*     OF THE SPACE H(P,DD) TO WHICH HIS INTEGRAND BELONGS, THE */
/*     ALGORITHM TERMINATES ACCORDING TO A HEURISTIC CRITERION, SEE */
/*     [2] AND COMMENTS TO EPS. */
/*        IF THE USER CAN SPECIFY THE PARAMETERS P AND D OF THE */
/*     SPACE H(P,DD) TO WHICH HIS INTEGRAND BELONGS, THE ALGORITHM */
/*     TERMINATES WITH AN ANSWER HAVING A GUARANTEED ACCURACY ( DE- */
/*     TEMINISTIC CRITERION, SEE [1, 2] AND COMMENTS TO EPS). */




/*     INPUT PARAMETERS */


/*     A = LOWER LIMIT OF INTEGRATION (SEE COMMENTS TO INF). */

/*     B = UPPER LIMIT OF INTEGRATION (SEE COMMENTS TO INF). */

/*     D = A PARAMETER OF THE CLASS H(P,DD) (SEE COMMENTS TO */
/*         INF). */

/*         USER SETS D: */

/*         HEURISTIC TERMINATION */
/*       = ANY REAL NUMBER */

/*         DETERMINISTIC TERMINATION */
/*       = A NUMBER IN THE RANGE 0.LT.D.LE.PI/2. */

/*     F = A NAME OF AN EXTERNAL INTEGRAND FUNCTION TO BE */
/*         SUPPLIED BY THE USER. F(X) COMPUTES THE VALUE OF */
/*         A FUNCTION F AT A POINT X. THE STATEMENT */
/*         ...EXTERNAL F... MUST APPEAR IN THE MAIN PROGRAM. */

/*     M = MAXIMAL NUMBER OF FUNCTION EVALUATIONS ALLOWED IN */
/*         THE COMPUTATIONS, M.GE.3.( ALTERED ON EXIT ). */

/*     P = 0, 1, .GT.1  A PARAMETER OF THE CLASS H(P,DD). */

/*         USER SETS P: */
/*       = 0 - HEURISTIC TERMINATION. */
/*       = 1 - DETERMINISTIC TERMINATION WITH THE INFINITY */
/*             NORM. */
/*      .GT.1 -DETERMINISTIC TERMINATION WITH THE P-TH NORM. */

/*   EPS = A REAL NUMBER - THE RELATIVE ERROR BOUND - SEE */
/*         REMARKS BELOW. ( ALTERED ON EXIT ). */


/*   INF = 1, 2, 3, 4 - INFORMATION PARAMETER. ( ALTERED ON EXIT ). */

/*       = 1 SIGNIFIES AN INFINITE INTERVAL (A,B)=REAL LINE, */
/*           A AND B ANY NUMBERS. */
/*           DETERMINISTIC TERMINATION - */
/*           DD=STRIP(Z:ABS(IM(Z)).LT.D). */

/*       = 2 SIGNIFIES A SEMI-INFINITE INTERVAL (A, +INFINITY) */
/*           USER SUPPLIES A, B ANY NUMBER. */
/*           QUADRATURE SUITED TO NON-OSCILLATORY INTEGRANDS. */
/*           DETERMINISTIC TERMINATION - */
/*           DD=SECTOR(Z:ABS(ARG(Z-A)).LT.D). */

/*       = 3 SIGNIFIES A SEMI INFINITE INTERVAL (A,+INFINITY) */
/*           USER SUPPLIES A, B ANY NUMBER. */
/*           QUADRATURE SUITED TO OSCILLATORY INTEGRANDS. */
/*           DETERMINISTIC TERMINATION - */
/*           DD=REGION(Z:ABS(ARG(SINH(Z-A))).LT.D). */

/*       = 4 SIGNIFIES A FINITE INTERVAL (A,B). */
/*           USER SUPPLIES A AND B. */
/*           DETERMINISTIC TERMINATION - */
/*           DD=LENS REGION(Z:ABS(ARG((Z-A)/(B-Z))).LT.D). */




/*     OUTPUT PARAMETERS */


/*     M = THE NUMBER OF FUNCTION EVALUATIONS USED IN THE */
/*         QUADRATURE. */

/*   EPS = THE RELATIVE ERROR BOUND (SEE REMARKS BELOW). */

/*         DETERMINISTIC TERMINATION */

/*       = THE RELATIVE ERROR REXA BOUND, I.E., */
/*                 REXA(F,M(OUTPUT)) .LE. EPS. */

/*         HEURISTIC TERMINATION */

/*       = MAX(EPS(INPUT),MACHEP). */

/*   INF = 0, 1 - DETERMINISTIC TERMINATION */

/*       = 0 COMPUTED QUADRATURE QCOM(F,M(EPS)), SEE REMARKS */
/*           BELOW. */

/*       = 1 COMPUTED QUADRATURE QCOM(F,M1), SEE REMARKS */
/*           BELOW. */

/*   INF = 2, 3, 4 - HEURISTIC TERMINATION. */

/*       = 2 INTEGRATION COMPLETED WITH EPS=MAX(EPS(INPUT), */
/*           MACHEP). WE CAN EXPECT THE RELATIVE ERROR */
/*           REXA TO BE OF THE ORDER OF EPS (FOR SOME P.GE.1). */

/*       = 3 INTEGRATION NOT COMPLETED. ATTEMPT TO EXCEED THE */
/*           MAXIMAL ALLOWED NUMBER OF FUNCTION EVALUATIONS M. */
/*           TRUNCATION CONDITIONS (SEE [2]) SATISFIED. QUADR */
/*           SET TO BE EQUAL TO THE LAST TRAPEZOIDAL APPRO- */
/*           XIMATION. IT IS LIKELY THAT QUADR APPROXIMATES THE */
/*           INTEGRAL IF M IS LARGE. */

/*       = 4 INTEGRATION NOT COMPLETED. ATTEMPT TO EXCEED THE */
/*           MAXIMAL ALLOWED NUMBER OF FUNCTION EVALUATIONS M. */
/*           TRUNCATION CONDITIONS (SEE [2]) NOT SATISFIED. */
/*           QUADR SET TO BE EQUAL TO THE COMPUTED TRAPEZOIDAL */
/*           APPROXIMATION. IT IS UNLIKELY THAT QUADR APPROXIMATES */
/*           THE INTEGRAL. */

/*   INF = 10, 11, 12, 13 - INCORRECT INPUT */

/*       = 10  M.LT.3. */

/*       = 11  P DOES NOT SATISFY P=0, P=1 OR P.GT.1 OR IN THE */
/*             CASE OF DETERMINISTIC TERMINATION D DOES NOT */
/*             SATISFY 0.LT.D.LE.PI/2. */

/*       = 12  A.GE.B IN CASE OF A FINITE INTERVAL. */

/*       = 13  INF NOT EQUAL TO 1, 2, 3, OR 4. */


/*   QUADR = THE COMPUTED VALUE OF QUADRATURE. */




/*     REMARKS: */

/*         LET  QEXA(F,M)  ( QCOM(F,M) ) BE THE EXACT (COMPUTED) */
/*         VALUE OF THE QUADRATURE WITH M FUNCTION EVALUATIONS, */
/*         AND LET  REXA(F,M) ( RCOM(F,M) ) BE THE RELATIVE ERROR */
/*         OF QEXA (QCOM) ,I.E., */
/*            REXA(F,M)=ABS(INTEGRAL(F)-QEXA(F,M))/NORM(F), */
/*            RCOM(F,M)=ABS(INTEGRAL(F)-QCOM(F,M))/NORM(F), */
/*         WITH THE NOTATION 0/0=0. */
/*             DUE TO THE ROUNDOFF ONE CANNOT EXPECT THE ERROR */
/*         RCOM TO BE LESS THAN THE RELATIVE MACHINE PRECISION */
/*         MACHEP. THEREFORE THE INPUT VALUE OF EPS IS CHANGED */
/*         ACCORDING TO THE FORMULA */
/*                   EPS=MAX(EPS,MACHEP). */

/*         DETERMINISTIC TERMINATON CASE */

/*             THE NUMBER OF FUNCTON EVALUATIONS M(EPS) IS COMPUTED */
/*         SO THAT THE ERROR REXA IS NO GREATER THAN EPS,I.E., */

/*         (*)     REXA(F,M(EPS)) .LE. EPS . */

/*         IF M(EPS).LE.M THEN THE QUADRATURE QCOM(F,M(EPS)) IS COM- */
/*         PUTED. OTHERWISE, WHICH MEANS THAT EPS IS TOO SMALL WITH */
/*         RESPECT TO M, THE QUADRATURE QCOM(F,M1) IS COMPUTED, WHERE */
/*         M1=2*INT((M-1)/2)+1. IN THIS CASE EPS IS CHANGED TO THE */
/*         SMALLEST NUMBER FOR WHICH THE ESTIMATE (*) HOLDS WITH */
/*         M(EPS)=M1 FUNCTION EVALUATIONS. */

/*         HEURISTIC TERMINATION CASE */

/*             WE CAN EXPECT THE RELATIVE ERROR REXA TO BE OF THE */
/*         ORDER OF EPS, SEE [2]. IF EPS IS TOO SMALL WITH RESPECT */
/*         TO M THEN THE QUADRATURE QCOM(F,M) IS COMPUTED. */

/*         ROUNDOFF ERRORS */

/*             IN BOTH DETERMINISTIC AND HEURISTIC CASES THE ROUND- */
/*         OFF ERROR */
/*                    ROFF=ABS(QEXA(F,M)-QCOM(F,M)) */
/*         CAN BE ESTIMATED BY */

/*         (**)       ROFF .LE. 3*C1*R*MACHEP, */

/*         WHERE  R=QCOM(ABS(F),M)+(1+2*C2)/3*SUM(W(I),I=1,2,...M) */
/*         AND C1 IS OF THE ORDER OF UNITY, C1=1/(1-3*MACHEP), W(I) */
/*         ARE THE WEIGHTS OF THE QUADRATURE, SEE [2], AND C2 IS */
/*         A CONSTANT ESTIMATING THE ACCURACY OF COMPUTING FUNCTION */
/*         VALUES, I.E., */
/*               ABS(EXACT(F(X))-COMPUTED(F(X))).LE.C2*MACHEP. */
/*         IF THE INTEGRAND VALUES ARE COMPUTED INACCURATELY, I.E., */
/*         C2 IS LARGE, THEN THE ESTIMATE (**) IS LARGE AND ONE CAN */
/*         EXPECT THE ACTUAL ERROR ROFF TO BE LARGE. NUMERICAL TESTS */
/*         INDICATE THAT THIS HAPPENS ESPECIALLY WHEN THE INTEGRAND */
/*         IS EVALUATED INACCURATELY NEAR A SINGULARITY. THE WAYS OF */
/*         CIRCUMVENTING SUCH PITFALLS ARE EXPLAINED IN [2]. */

/*     REFERENCES: */

/*     [1] SIKORSKI,K., OPTIMAL QUADRATURE ALGORITHMS IN HP */
/*            SPACES, NUM. MATH., 39, 405-410 (1982). */
/*     [2] SIKORSKI,K., STENGER,F., OPTIMAL QUADRATURES IN */
/*            HP SPACES, ACM TOMS. */





/*     modified by J.K. Lindsey for R September 1998 */
/*<       implicit none >*/
/*<       INTEGER I, I1, INF, K, L, L1, M, M1, M2, N, N1 >*/
/*<    >*/
/*<    >*/
/*<       double precision V2, W, W1, W2, W3, W4 >*/
/*<       LOGICAL INF1, INF2 >*/
SEXP call2, result2, state2, tmp;
  // mode[0] = "double";
  // length[0] = 1;
  // args[0] = (void *)(zz);

/*<       PI = 4.*ATAN(1.0) >*/
    pi = atan(1.f) * 4.f;

/*     CHECK THE INPUT DATA */

/*<    >*/
    if (*INF != 1 && *INF != 2 && *INF != 3 && *INF != 4) {
	goto L300;
    }
/*<       IF (M.LT.3) GO TO 270 >*/
    if (*M < 3) {
	goto L270;
    }
/*<       IF (P.LT.1. .AND. P.NE.0.) GO TO 280 >*/
    if (*P < 1.f && *P != 0.f) {
	goto L280;
    }
/*<       IF (P.GE.1. .AND. (D.LE.0. .OR. D.GT.PI/2.)) GO TO 280 >*/
    if (*P >= 1.f && (*D__ <= 0.f || *D__ > pi / 2.f)) {
	goto L280;
    }
/*<       IF (INF.EQ.4 .AND. A.GE.B) GO TO 290 >*/
    if (*INF == 4 && *A >= *B) {
	goto L290;
    }

/*<       SQ2 = SQRT(2.0) >*/
    sq2 = sqrt(2.f);
/*<       I1 = INF - 2 >*/
    i1 = *INF - 2;
/*<       BA = B - A >*/
    ba = *B - *A;
/*<       N1 = 0 >*/
    n1 = 0;

/*     COMPUTE THE RELATIVE MACHINE PRECISION AND CHECK */
/*     THE VALUE OF EPS.  CAUTION...THIS LOOP MAY NOT WORK ON A */
/*     MACHINE THAT HAS AN ACCURATED ARITHMETIC PROCESS COMPARED */
/*     TO THE STORAGE PRECISION.  THE VALUE OF U MAY NEED TO BE */
/*     SIMPLY DEFINED AS THE RELATIVE ACCURACY OF STORAGE PRECISION. */

/*<       U = 1. >*/
    u = 1.f;
/*<    10 U = U/10. >*/
L10:
    u /= 10.f;
/*<       T = 1. + U >*/
    t = u + 1.f;
/*<       IF (1..NE.T) GO TO 10 >*/
    if (1.f != t) {
	goto L10;
    }
/*<       U = U*10. >*/
    u *= 10.f;
/*<       IF (EPS.LT.U) EPS = U >*/
    if (*EPS < u) {
	*EPS = u;
    }

/*<       IF (P.EQ.0.) GO TO 40 >*/
    if (*P == 0.f) {
	goto L40;
    }

/*     SET UP DATA FOR THE DETERMINISTIC TERMINATION */

/*<       IF (P.EQ.1.) ALFA = 1. >*/
    if (*P == 1.f) {
	alfa = 1.f;
    }
/*<       IF (P.GT.1.) ALFA = (P-1.)/P >*/
    if (*P > 1.f) {
	alfa = (*P - 1.f) / *P;
    }
/*<       C = 2.*PI/(1.-1./EXP(PI*SQRT(ALFA))) + 4.**ALFA/ALFA >*/
    c__ = pi * 2.f / (1.f - 1.f / exp(pi * sqrt(alfa))) + R_pow(&c_b8, &alfa)
	     / alfa;
/*<       W = dLOG(C/EPS) >*/
    w = log(c__ / *EPS);
/*<       W1 = 1./(PI*PI*ALFA)*W*W >*/
    w1 = 1.f / (pi * pi * alfa) * w * w;
/*<       N = INT(W1) >*/
    n = (integer) w1;
/*<       IF (W1.GT.FLOAT(N)) N = N + 1 >*/
    if (w1 > (real) n) {
	++n;
    }
/*<       IF (W1.EQ.0.) N = 1 >*/
    if (w1 == 0.f) {
	n = 1;
    }
/*<       N1 = 2*N + 1 >*/
    n1 = (n << 1) + 1;
/*<       SR = SQRT(ALFA*FLOAT(N)) >*/
    sr = sqrt(alfa * (real) n);
/*<       IF (N1.LE.M) GO TO 20 >*/
    if (n1 <= *M) {
	goto L20;
    }

/*     EPS TOO SMALL WITH RESPECT TO M. COMPUTE THE NEW EPS */
/*     GUARANTEED BY THE VALUE OF M. */

/*<       N1 = 1 >*/
    n1 = 1;
/*<       N = INT(FLOAT((M-1)/2)) >*/
    n = (integer) ((real) ((*M - 1) / 2));
/*<       SR = SQRT(ALFA*FLOAT(N)) >*/
    sr = sqrt(alfa * (real) n);
/*<       M = 2*N + 1 >*/
    *M = (n << 1) + 1;
/*<       EPS = C/EXP(PI*SR) >*/
    *EPS = c__ / exp(pi * sr);
/*<       GO TO 30 >*/
    goto L30;

/*<    20 M = N1 >*/
L20:
    *M = n1;
/*<       N1 = 0 >*/
    n1 = 0;
/*<    30 H = 2.*D/SR >*/
L30:
    h__ = *D__ * 2.f / sr;
/*<       SUM2 = 0. >*/
    sum2 = 0.f;
/*<       L1 = N >*/
    l1 = n;
/*<       K = N >*/
    k = n;
/*<       INF1 = .FALSE. >*/
    inf1 = FALSE_;
/*<       INF2 = .FALSE. >*/
    inf2 = FALSE_;
/*<       H0 = H >*/
    h0 = h__;
/*<       GO TO 50 >*/
    goto L50;

/*     SET UP DATA FOR THE HEURISTIC TERMINATION */

/*<    40 H = 1. >*/
L40:
    h__ = 1.f;
/*<       H0 = 1. >*/
    h0 = 1.f;
/*<       EPS3 = EPS/3. >*/
    eps3 = *EPS / 3.f;
/*<       SR = SQRT(EPS) >*/
    sr = sqrt(*EPS);
/*<       V1 = EPS*10. >*/
    v1 = *EPS * 10.f;
/*<       V2 = V1 >*/
    v2 = v1;
/*<       M1 = M - 1 >*/
    m1 = *M - 1;
/*<       N = INT(FLOAT(M1/2)) >*/
    n = (integer) ((real) (m1 / 2));
/*<       M2 = N >*/
    m2 = n;
/*<       L1 = 0 >*/
    l1 = 0;
/*<       INF1 = .TRUE. >*/
    inf1 = TRUE_;
/*<       INF2 = .FALSE. >*/
    inf2 = FALSE_;

/*     INITIALIZE THE QUADRATURE */

/*<    50 I = 0 >*/
L50:
    i__ = 0;
/*<       IF (INF.EQ.1) SUM = F(0.d0) >*/
    if (*INF == 1) {
      /*	sum = (*f)(&c_b12);*/
    zz[0]=c_b12;
      
      state2 = Rf_protect(Rf_allocVector(REALSXP, 1)); /*.C -> .Call requires this*/
      REAL(state2)[0] = zz[0];                                   /*.C -> .Call requires this*/
      Rf_protect(call2=Rf_lang2(f, state2));        /*.C -> .Call requires this*/
      Rf_protect(result2=Rf_eval(call2, envir));     /*.C -> .Call requires this*/
      Rf_protect(tmp=Rf_coerceVector(result2, REALSXP));/*.C -> .Call requires this*/    
      
    sum=REAL(tmp)[0];
    Rf_unprotect(4);
    }
/*<       IF (INF.EQ.2) SUM = F(A+1.) >*/
    if (*INF == 2) {
	d__1 = *A + 1.f;
	/*	sum = (*f)(&d__1);*/
    zz[0]=d__1;
    
    state2 = Rf_protect(Rf_allocVector(REALSXP, 1)); /*.C -> .Call requires this*/
	REAL(state2)[0] = zz[0];                                   /*.C -> .Call requires this*/
	Rf_protect(call2=Rf_lang2(f, state2));        /*.C -> .Call requires this*/
	Rf_protect(result2=Rf_eval(call2, envir));     /*.C -> .Call requires this*/
	Rf_protect(tmp=Rf_coerceVector(result2, REALSXP));/*.C -> .Call requires this*/    
	
    sum=REAL(tmp)[0];
    Rf_unprotect(4);
    }
/*<       IF (INF.EQ.3) SUM = F(A+dLOG(1.+SQ2))/SQ2 >*/
    if (*INF == 3) {
	d__1 = *A + log(sq2 + 1.f);
	/*	sum = (*f)(&d__1) / sq2;*/
    zz[0]=d__1;
    
    state2 = Rf_protect(Rf_allocVector(REALSXP, 1)); /*.C -> .Call requires this*/
	REAL(state2)[0] = zz[0];                                   /*.C -> .Call requires this*/
	Rf_protect(call2=Rf_lang2(f, state2));        /*.C -> .Call requires this*/
	Rf_protect(result2=Rf_eval(call2, envir));     /*.C -> .Call requires this*/
	Rf_protect(tmp=Rf_coerceVector(result2, REALSXP));/*.C -> .Call requires this*/    
	
    sum=REAL(tmp)[0]/sq2;
    Rf_unprotect(4);
    }
/*<       IF (INF.EQ.4) SUM = F((A+B)/2.)/4.*BA >*/
    if (*INF == 4) {
	d__1 = (*A + *B) / 2.f;
	/*	sum = (*f)(&d__1) / 4.f * ba;*/
    zz[0]=d__1;
    
    state2 = Rf_protect(Rf_allocVector(REALSXP, 1)); /*.C -> .Call requires this*/
	REAL(state2)[0] = zz[0];                                   /*.C -> .Call requires this*/
	Rf_protect(call2=Rf_lang2(f, state2));        /*.C -> .Call requires this*/
	Rf_protect(result2=Rf_eval(call2, envir));     /*.C -> .Call requires this*/
	Rf_protect(tmp=Rf_coerceVector(result2, REALSXP));/*.C -> .Call requires this*/    
	
    sum=REAL(tmp)[0] / 4.f * ba;
    Rf_unprotect(4);
    }

/*     COMPUTE WEIGHTS, NODES AND FUNCTION VALUES */

/*<    60 EXPH = EXP(H) >*/
L60:
    exph = exp(h__);
/*<       EXPH0 = EXP(H0) >*/
    exph0 = exp(h0);
/*<       H1 = H0 >*/
    h1 = h0;
/*<       E1 = EXPH0 >*/
    e1 = exph0;
/*<       U = 0. >*/
    u = 0.f;
/*<       COR = 0. >*/
    cor = 0.f;

/*<    70 IF (I1) 80, 90, 100 >*/
L70:
    if (i1 < 0) {
	goto L80;
    } else if (i1 == 0) {
	goto L90;
    } else {
	goto L100;
    }

/*<    80 V = F(H1) >*/
L80:
    /*    v = (*f)(&h1);*/
    zz[0]=h1;
    
    state2 = Rf_protect(Rf_allocVector(REALSXP, 1)); /*.C -> .Call requires this*/
    REAL(state2)[0] = zz[0];                                   /*.C -> .Call requires this*/
    Rf_protect(call2=Rf_lang2(f, state2));        /*.C -> .Call requires this*/
    Rf_protect(result2=Rf_eval(call2, envir));     /*.C -> .Call requires this*/
    Rf_protect(tmp=Rf_coerceVector(result2, REALSXP));/*.C -> .Call requires this*/    
    
    v=REAL(tmp)[0];
    Rf_unprotect(4);
/*<       H1 = H1 + H >*/
    h1 += h__;
/*<       GO TO 150 >*/
    goto L150;

/*<    90 V = E1*F(A+E1) >*/
L90:
    d__1 = *A + e1;
    /*    v = e1 * (*f)(&d__1);*/
    zz[0]=d__1;
    
    state2 = Rf_protect(Rf_allocVector(REALSXP, 1)); /*.C -> .Call requires this*/
    REAL(state2)[0] = zz[0];                                   /*.C -> .Call requires this*/
    Rf_protect(call2=Rf_lang2(f, state2));        /*.C -> .Call requires this*/
    Rf_protect(result2=Rf_eval(call2, envir));     /*.C -> .Call requires this*/
    Rf_protect(tmp=Rf_coerceVector(result2, REALSXP));/*.C -> .Call requires this*/    
    
    v=e1 * REAL(tmp)[0];
    Rf_unprotect(4);
/*<       E1 = E1*EXPH >*/
    e1 *= exph;
/*<       GO TO 150 >*/
    goto L150;

/*<   100 IF (INF.EQ.4) GO TO 140 >*/
L100:
    if (*INF == 4) {
	goto L140;
    }
/*<       W1 = SQRT(E1+1./E1) >*/
    w1 = sqrt(e1 + 1.f / e1);
/*<       W2 = SQRT(E1) >*/
    w2 = sqrt(e1);
/*<       IF (E1.LT.0.1) GO TO 110 >*/
    if (e1 < .1f) {
	goto L110;
    }
/*<       S = dLOG(E1+W1*W2) >*/
    s = log(e1 + w1 * w2);
/*<       GO TO 130 >*/
    goto L130;
/*<   110 W3 = E1 >*/
L110:
    w3 = e1;
/*<       W4 = E1*E1 >*/
    w4 = e1 * e1;
/*<       C0 = 1. >*/
    c0 = 1.f;
/*<       S = E1 >*/
    s = e1;
/*<       S1 = E1 >*/
    s1 = e1;
/*<       T = 0. >*/
    t = 0.f;
/*<   120 C0 = -C0*(0.5+T)*(2.*T+1.)/(2.*T+3.)/(T+1.) >*/
L120:
    c0 = -c0 * (t + .5f) * (t * 2.f + 1.f) / (t * 2.f + 3.f) / (t + 1.f);
/*<       T = T + 1. >*/
    t += 1.f;
/*<       W3 = W3*W4 >*/
    w3 *= w4;
/*<       S = S + C0*W3 >*/
    s += c0 * w3;
/*<       IF (S.EQ.S1) GO TO 130 >*/
    if (s == s1) {
	goto L130;
    }
/*<       S1 = S >*/
    s1 = s;
/*<       GO TO 120 >*/
    goto L120;
/*<   130 V = W2/W1*F(A+S) >*/
L130:
    d__1 = *A + s;
    /*    v = w2 / w1 * (*f)(&d__1);*/
    zz[0]=d__1;
    
    state2 = Rf_protect(Rf_allocVector(REALSXP, 1)); /*.C -> .Call requires this*/
    REAL(state2)[0] = zz[0];                                   /*.C -> .Call requires this*/
    Rf_protect(call2=Rf_lang2(f, state2));        /*.C -> .Call requires this*/
    Rf_protect(result2=Rf_eval(call2, envir));     /*.C -> .Call requires this*/
    Rf_protect(tmp=Rf_coerceVector(result2, REALSXP));/*.C -> .Call requires this*/    
    
    v = w2 / w1 * REAL(tmp)[0];
    Rf_unprotect(4);
/*<       E1 = E1*EXPH >*/
    e1 *= exph;
/*<       GO TO 150 >*/
    goto L150;

/*<   140 W1 = E1 + 1. >*/
L140:
    w1 = e1 + 1.f;
/*<       V = E1/W1/W1*F((A+B*E1)/W1)*BA >*/
    d__1 = (*A + *B * e1) / w1;
    /*    v = e1 / w1 / w1 * (*f)(&d__1) * ba;*/
    zz[0]=d__1;
    
    state2 = Rf_protect(Rf_allocVector(REALSXP, 1)); /*.C -> .Call requires this*/
    REAL(state2)[0] = zz[0];                                   /*.C -> .Call requires this*/
    Rf_protect(call2=Rf_lang2(f, state2));        /*.C -> .Call requires this*/
    Rf_protect(result2=Rf_eval(call2, envir));     /*.C -> .Call requires this*/
    Rf_protect(tmp=Rf_coerceVector(result2, REALSXP));/*.C -> .Call requires this*/    
    
    v = e1 / w1 / w1 * REAL(tmp)[0] * ba;
    Rf_unprotect(4);
/*<       E1 = E1*EXPH >*/
    e1 *= exph;

/*     SUMMATION ALGORITHM */

/*<   150 I = I + 1 >*/
L150:
    ++i__;
/*<       SUM1 = U + V >*/
    sum1 = u + v;
/*<       IF (ABS(U).LT.ABS(V)) GO TO 160 >*/
    if (abs(u) < abs(v)) {
	goto L160;
    }
/*<       COR = V - (SUM1-U) + COR >*/
    cor = v - (sum1 - u) + cor;
/*<       GO TO 170 >*/
    goto L170;
/*<   160 COR = U - (SUM1-V) + COR >*/
L160:
    cor = u - (sum1 - v) + cor;
/*<   170 U = SUM1 >*/
L170:
    u = sum1;
/*<       IF (I.LT.L1) GO TO 70 >*/
    if (i__ < l1) {
	goto L70;
    }

/*     SWITCH TO CHECK TRUNCATION CONDITION ( HEURISTIC */
/*     TERMINATION) */

/*<       IF (INF1) GO TO 190 >*/
    if (inf1) {
	goto L190;
    }

/*     SWITCH TO COMPUTE THE MIDORDINATE APPROXIMATION */
/*     ( HEURISTIC TERMINATION ) OR TO STOP ( DETERMINIS- */
/*     TIC TERMINATION) */

/*<       IF (INF2) GO TO 210 >*/
    if (inf2) {
	goto L210;
    }

/*     SET UP PARAMETERS TO CONTINUE SUMMATION */

/*<       L1 = K >*/
    l1 = k;
/*<   180 INF2 = .TRUE. >*/
L180:
    inf2 = TRUE_;
/*<       I = 0. >*/
    i__ = 0.f;
/*<       EXPH = 1./EXPH >*/
    exph = 1.f / exph;
/*<       H0 = -H0 >*/
    h0 = -h0;
/*<       E1 = 1./EXPH0 >*/
    e1 = 1.f / exph0;
/*<       H1 = H0 >*/
    h1 = h0;
/*<       H = -H >*/
    h__ = -h__;
/*<       GO TO 70 >*/
    goto L70;

/*     TRUNCATION CONDITION */

/*<   190 V0 = V1 >*/
L190:
    v0 = v1;
/*<       V1 = V2 >*/
    v1 = v2;
/*<       V2 = ABS(V) >*/
    v2 = abs(v);
/*<       IF (V0+V1+V2.LE.EPS3) GO TO 200 >*/
    if (v0 + v1 + v2 <= eps3) {
	goto L200;
    }
/*<       IF (I.LT.M2) GO TO 70 >*/
    if (i__ < m2) {
	goto L70;
    }
/*<       N1 = 5 >*/
    n1 = 5;
/*<   200 IF (INF2) K = I >*/
L200:
    if (inf2) {
	k = i__;
    }
/*<       IF (.NOT.INF2) L = I >*/
    if (! inf2) {
	l = i__;
    }
/*<       V1 = 10.*EPS >*/
    v1 = *EPS * 10.f;
/*<       V2 = V1 >*/
    v2 = v1;
/*<       M2 = M1 - L >*/
    m2 = m1 - l;
/*<       IF (.NOT.INF2) GO TO 180 >*/
    if (! inf2) {
	goto L180;
    }

/*     N1=5 - TRUNCATION CONDITION NOT SATISFIED */

/*<       IF (N1.EQ.5) GO TO 260 >*/
    if (n1 == 5) {
	goto L260;
    }

/*     TRUNCATION CONDITION SATISFIED, SUM2=TRAPEZOIDAL */
/*     APPROXIMATION */

/*<       SUM2 = SUM1 + COR + SUM >*/
    sum2 = sum1 + cor + sum;
/*<       M2 = 2*(K+L) >*/
    m2 = (k + l) << 1;

/*     CHECK THE NUMBER OF FUNCTION EVALUATIONS */

/*<       IF (M2.GT.M1) GO TO 240 >*/
    if (m2 > m1) {
	goto L240;
    }

/*     INITIALIZE ITERATION */

/*<       INF1 = .FALSE. >*/
    inf1 = FALSE_;
/*<       INF2 = .FALSE. >*/
    inf2 = FALSE_;
/*<       L1 = L >*/
    l1 = l;
/*<       I = 0 >*/
    i__ = 0;
/*<       H = -H >*/
    h__ = -h__;
/*<       H0 = H/2. >*/
    h0 = h__ / 2.f;
/*<       GO TO 60 >*/
    goto L60;

/*     P.GE.1 = DETERMINISTIC TERMINATION */

/*<   210 IF (P.GE.1.) GO TO 220 >*/
L210:
    if (*P >= 1.f) {
	goto L220;
    }

/*     COMPUTE THE MIDORDINATE APPROXIMATION SUM1 */

/*<       H = -H >*/
    h__ = -h__;
/*<       SUM1 = (SUM1+COR)*H >*/
    sum1 = (sum1 + cor) * h__;
/*<       W1 = (SUM1+SUM2)/2. >*/
    w1 = (sum1 + sum2) / 2.f;

/*     TERMINATION CONDITION */

/*<       IF (ABS(SUM1-SUM2).LE.SR) GO TO 230 >*/
    if ((d__1 = sum1 - sum2, abs(d__1)) <= sr) {
	goto L230;
    }

/*     SET UP DATA FOR THE NEXT ITERATION */

/*<       M2 = 2*M2 >*/
    m2 <<= 1;
/*<       IF (M2.GT.M1) GO TO 250 >*/
    if (m2 > m1) {
	goto L250;
    }
/*<       I = 0 >*/
    i__ = 0;
/*<       K = 2*K >*/
    k <<= 1;
/*<       L = 2*L >*/
    l <<= 1;
/*<       L1 = L >*/
    l1 = l;
/*<       H = H/2. >*/
    h__ /= 2.f;
/*<       H0 = H/2. >*/
    h0 = h__ / 2.f;
/*<       SUM2 = W1 >*/
    sum2 = w1;
/*<       INF2 = .FALSE. >*/
    inf2 = FALSE_;
/*<       GO TO 60 >*/
    goto L60;

/*     FINAL RESULTS */
Rf_protect(ans = Rf_allocVector(REALSXP, 1));
/*<   220 QUADR = -H*(SUM1+COR+SUM) >*/
L220:
    *quadr = -h__ * (sum1 + cor + sum);
/*<       INF = N1 >*/
    *INF = n1;
/*<       RETURN >*/
    REAL(ans)[0] = quadr[0]; 
    Rf_unprotect(1); 
    return ans;// return 0;

/*<   230 QUADR = W1 >*/
L230:
    *quadr = w1;
/*<       INF = 2 >*/
    *INF = 2;
/*<       M = M2 + 1 >*/
    *M = m2 + 1;
/*<       RETURN >*/
    REAL(ans)[0] = quadr[0]; 
    Rf_unprotect(1); 
    return ans;// return 0;

/*<   240 QUADR = SUM2 >*/
L240:
    *quadr = sum2;
/*<       INF = 3 >*/
    *INF = 3;
/*<       M = K + L + 1 >*/
    *M = k + l + 1;
/*<       RETURN >*/
    REAL(ans)[0] = quadr[0]; 
    Rf_unprotect(1); 
    return ans;// return 0;

/*<   250 QUADR = W1 >*/
L250:
    *quadr = w1;
/*<       INF = 3 >*/
    *INF = 3;
/*<       M = M2/2 + 1 >*/
    *M = m2 / 2 + 1;
/*<       RETURN >*/
    REAL(ans)[0] = quadr[0]; 
    Rf_unprotect(1); 
    Rf_error("INF in toms614_sexp.c is 3 -- integration incomplete - try larger max");
    return ans;// return 0;

/*<   260 QUADR = U + COR + SUM >*/
L260:
    *quadr = u + cor + sum;
/*<       INF = 4 >*/
    *INF = 4;
/*<       M = K + L + 1 >*/
    *M = k + l + 1;
/*<       RETURN >*/
    REAL(ans)[0] = quadr[0]; 
    Rf_unprotect(1); 
    Rf_error("INF in toms614_sexp.c is 4 -- integration incomplete - try larger max");
    return ans;// return 0;

/*<   270 INF = 10 >*/
L270:
    *INF = 10;
/*<       RETURN >*/
    //REAL(ans)[0] = quadr[0]; 
    Rf_unprotect(1); 
    Rf_error("INF in toms614_sexp.c is 10 -- incorrect arguments");
    return ans;// return 0;

/*<   280 INF = 11 >*/
L280:
    *INF = 11;
/*<       RETURN >*/
    //REAL(ans)[0] = quadr[0]; 
    Rf_unprotect(1); 
    Rf_error("INF in toms614_sexp.c is 11 -- incorrect arguments");
    return ans;// return 0;

/*<   290 INF = 12 >*/
L290:
    *INF = 12;
/*<       RETURN >*/
    //REAL(ans)[0] = quadr[0]; 
    Rf_unprotect(1); 
    Rf_error("INF in toms614_sexp.c is 12 -- incorrect arguments");
    return ans;// return 0;

/*<   300 INF = 13 >*/
L300:
    *INF = 13;
/*<       RETURN >*/
    //REAL(ans)[0] = quadr[0]; 
    Rf_unprotect(1); 
    Rf_error("INF in toms614_sexp.c is 13 -- incorrect arguments");
    return ans;// return 0;

/*<       END >*/
} /* inthp_ */

