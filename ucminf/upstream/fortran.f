      SUBROUTINE UCMINF(N,X,DX,EPS,MAXFUN,W,IW,ICONTR,GRAD,GRSTEP,RHO)
************************************************************************
* Unconstrained minimization of a scalar function.
*
* For a User's guide see H.B. Nielsen  "UCMINF -- AN ALGORITHM FOR
* UNCONSTRAINED, NONLINEAR OPTIMIZATION", Report IMM-REP-2000-18,
* Department of Mathematical Modelling, Technical University of Denmark,
* December 2000.
*
* Hans Bruun Nielsen, IMM, DTU.  00.12.19
*
*   Changes for implementation in R:
*     - FDF removed from argument to UCMINF and thus not declared EXTERNAL.
*     - GRAD, GRADSTEP are passed to FDF to choose gradient type.
*     - OUT is changed from INTEGER to CHARACTER to store output from
*       write which is then displayed in R using DBLEPR.
*     - DBLEPR is declared EXTERNAL.
*   Stig Mortensen, IMM, DTU. December 2008.
*
************************************************************************
      IMPLICIT          NONE
C  Parameters
C     EXTERNAL          FDF
      INTEGER           N, MAXFUN, IW, ICONTR, GRAD, RHO(*)
      DOUBLE PRECISION  X(N), DX, EPS(2), W(IW), GRSTEP(2)
C  Local variables
      LOGICAL           DGIVN, OPTIM, REDU, TRACE, USEDEL
      INTEGER           DI,DN,FAIL,GN,GP,HN,I,II,INDX(3),
     &                  MEVAL,NEVAL,NN,NN1,WN
      DOUBLE PRECISION  A,FX,FXN,NMG,NMH,NMX,SL(2),THRX,YH,YV
      INTRINSIC         ABS,DBLE,MAX,MIN
C
C  BLAS functions
      INTEGER           IDAMAX
      DOUBLE PRECISION  DDOT, DNRM2
C If BLAS is available, then remove the * in col. 1 of the next lines
C and delete from line ??? to the end of the file
      EXTERNAL          DCOPY, DDOT, DNRM2, DSCAL, DSPMV, DSPR2, IDAMAX

C       ... What job ?
        OPTIM = (ICONTR .GT. 0)
        DGIVN = (ICONTR .GT. 2)
        TRACE = ((ICONTR .EQ. 2) .OR. (ICONTR .GT. 3))
C       ... Simple checks
        ICONTR = 0
        NN = (N * (N+1)) / 2
        IF  (N .LE. 0)  THEN
          ICONTR = -2
        ELSEIF (OPTIM)  THEN
          IF  (DX .LE. 0D0)  THEN
            ICONTR = -4
          ELSEIF  ((EPS(1) .LE. 0D0) .OR. (EPS(2) .LE. 0D0))  THEN
            ICONTR = -5
          ELSEIF  (MAXFUN .LE. 0)  THEN
            ICONTR = -6
          ENDIF
        ELSEIF  (DX .EQ. 0D0)  THEN
          ICONTR = -4
        ELSEIF  ((IW .LT. MAX((N*(N+11))/2, 7)) .OR.
     +           (DGIVN .AND. (IW .LT. MAX(2*NN, NN+5*N))))  THEN
          ICONTR = -8
        ENDIF
C       ... Exit if error in a parameter
        IF (ICONTR .LT. 0)  RETURN
C
        IF  (.NOT. OPTIM)  THEN
C         ... Check gradient
          GN = 5
          HN = GN + N
          CALL CHKDFN(N,X,DX,W,INDX,W(GN),W(HN),FAIL,GRAD,GRSTEP,RHO)
          IF  (FAIL .GT. 0)  THEN
C           ... DX is too small
            ICONTR = -4
          ELSE
            DO  I = 1, 3
              W(I+4) = DBLE(INDX(I))
            END DO
          ENDIF
          RETURN
        ENDIF
C
C       ... Optimize.  Split workspace
        GP = N + 1
        GN = GP + N
        HN = GN + N
        DN = HN + N
        WN = DN + NN
C
        IF  (DGIVN)  THEN
C         ... Check given  D0
          NN1 = NN + 1
          CALL DCOPY(NN, W(DN),1, W,1)
          CALL DCOPY(NN, W,1, W(NN1),1)
          CALL SPCHOL(N,W(NN1),FAIL)
          IF  (FAIL .NE. 0)  THEN
C           ... Not positive definite
            ICONTR = -7
            RETURN
          ENDIF
C         ... Restore given  D
          CALL DCOPY(NN, W,-1, W(DN),-1)
          USEDEL = .FALSE.
        ELSE
C         ... Initialize inverse Hessian to unit matrix
          DO  I = HN, WN
            W(I) = 0D0
          END DO
          II = DN
          DI = N
          DO  I = 1, N
            W(II) = 1D0
            II = II + DI
            DI = DI - 1
         END DO
          USEDEL = .TRUE.
        ENDIF
C       ... First call of FDF
        CALL FDF(N,X,W(GN),FX,GRAD,GRSTEP,RHO)
        NEVAL = 1
        NMH = 0D0
        NMX = DNRM2(N, X,1)
        NMG = ABS(W(GN-1 + IDAMAX(N, W(GN),1)))
        IF  (NMG .LE. EPS(1))  THEN
          ICONTR = 1
          GOTO 200
        ENDIF
C
C       ... Repeat from here
 100    CONTINUE
        IF  (TRACE)  THEN
           CALL PRTRAC(NEVAL, FX, NMG, N, X)
c$$$          WRITE(OUT,'(A,I3,2(2X,A,1P1D11.3))') 'neval =',NEVAL,
c$$$     +      'F(x) =',FX, 'max|g(x)| =',NMG
c$$$          CALL DBLEPR (OUT, -1, X, 0)
c$$$          CALL PRVCTR('  x',X,1,N,OUT)
        ENDIF
C
C       ... Copy current x and gradient and get new step
        CALL DCOPY(N, X,1, W,1)
        CALL DCOPY(N, W(GN),1, W(GP),1)
        CALL DSPMV('L', N, -1D0, W(DN), W(GN),1, 0D0,W(HN),1)
C       ... Adjust step length to trust region
        REDU = .FALSE.
        NMH = DNRM2(N, W(HN),1)
        IF  (NMH .LE. EPS(2)*(EPS(2) + NMX))  THEN
          ICONTR = 2
          GOTO 200
        ENDIF
        IF  ((NMH .GT. DX) .OR. USEDEL)  THEN
          REDU = .TRUE.
          CALL DSCAL(N, DX/NMH, W(HN),1)
          NMH = DX
          USEDEL = .FALSE.
        ENDIF
C       ... Line search (MEVAL is max iterations in line search)
        MEVAL = 5
        CALL SLINE(N,X,FX,W(GN),W(HN),W(WN),A,FXN,SL,MEVAL,
     &             GRAD,GRSTEP,RHO)
        IF  (TRACE)  THEN
           CALL PRLINE(A, SL)
c$$$          WRITE(OUT,'(A,1P1D11.3,2(2X,A,1P1D11.3))')
c$$$     +      'Line search: alpha =',
c$$$     +      A, 'dphi(0) =',SL(1),'dphi(alpha) =',SL(2)
c$$$          CALL DBLEPR (OUT, -1, X, 0)
        ENDIF
        IF  (A .EQ. 0D0)  THEN
          ICONTR = 4
          NMH = 0D0
          GOTO 200
        ENDIF
C       ... Update  neval, x, f(x) and ||g||
        NEVAL = NEVAL + MEVAL
        NMG = ABS(W(GN-1 + IDAMAX(N, W(GN),1)))
        FX = FXN
        CALL DAXPY(N, A, W(HN),1, X,1)
        NMX = DNRM2(N, X,1)
        CALL DAXPY(N, -1D0, X,1, W,1)
        NMH = DNRM2(N, W,1)
C       ... Update trust region
        IF  (A .LT. 1D0)  THEN
C         ... Reduce Delta
          DX = .35D0 * DX
        ELSEIF  (REDU .AND. (SL(2) .LT. .7D0*SL(1)))  THEN
C         ... Increase Delta
          DX = 3D0 * DX
        ENDIF
C       ... Update of inverse Hessian by BFGS
        CALL DSCAL(N, -1D0, W,1)
        CALL DAXPY(N, -1D0,W(GN),1, W(GP),1)
        YH = -DDOT(N, W(GP),1, W,1)
        IF  (YH .GT. 1D-8 * NMH * DNRM2(N, W(GP),1))  THEN
          CALL DSPMV('L', N, -1D0, W(DN), W(GP),1, 0D0,W(HN),1)
          YV = -DDOT(N, W(GP),1, W(HN),1)
          A = (1D0 + YV/YH)/YH
          CALL DSCAL(N, -1D0/YH, W(HN),1)
          CALL DAXPY(N, .5D0*A, W,1, W(HN),1)
          CALL DSPR2('L', N, 1D0, W,1, W(HN),1, W(DN))
        ENDIF
C       ...  Check stopping criteria
        THRX = EPS(2)*(EPS(2) + NMX)
        DX = MAX(DX, THRX)
        IF  (NEVAL .GE. MAXFUN)  ICONTR = 3
        IF  (NMH .LE. THRX)      ICONTR = 2
        IF  (NMG .LE. EPS(1))    ICONTR = 1
        IF  (ICONTR .EQ. 0)  GOTO 100
C
 200    CONTINUE
C       ... Set return values
        MAXFUN = NEVAL
        W(1) = FX
        W(2) = NMG
        W(3) = NMH
        IF  (TRACE)  THEN
           IF(ICONTR.EQ.1.OR.ICONTR.EQ.2.OR.ICONTR.EQ.4) THEN
c$$$              WRITE(OUT,'(A)') 'Optimization has converged.'
              CALL PRCONV()
           ELSE
c$$$              WRITE(OUT,'(A27,I3,A22)') 'Optimization stopped after '
c$$$     +             ,NEVAL,' function evaluations.'
              CALL PRFAIL(NEVAL)
           ENDIF
c$$$           CALL DBLEPR (OUT, -1, X, 0)
        ENDIF
        RETURN
***************************  End of  UCMINF  ***************************
      END
C
      SUBROUTINE SLINE(N,X,F,G,H,W,ALPHA,FN,SLPS,NEV,GRAD,GRSTEP,RHO)
************************************************************************
* Soft line search
* Hans Bruun Nielsen, IMM, DTU.  00.12.18
************************************************************************
      IMPLICIT          NONE
C  Parameters
C      EXTERNAL          FDF
      INTEGER           N, NEV, GRAD, RHO(*)
      DOUBLE PRECISION  X(N),F,G(N),H(N),W(*),ALPHA,FN,SLPS(2)
      DOUBLE PRECISION  GRSTEP(2)
C  Local variables
      LOGICAL           OK,STOP
      INTEGER           MEVAL,NG,NX
      DOUBLE PRECISION  A,B,C,D,FI0,SL0,SLTHR,XFD(3,3)
      INTRINSIC         ABS,DBLE,MAX,MIN
C  BLAS functions
      DOUBLE PRECISION  DDOT
C If BLAS is available, then remove the * in col. 1 of the next lines
C and delete from line ??? to the end of the file
      EXTERNAL          DAXPY, DCOPY, DDOT
C
C       ... Default return values
        ALPHA = 0D0
        FN = F
        MEVAL = NEV
        NEV = 0
C       ... Get initial slope and check descent direction
        SLPS(1) = DDOT(N, G,1, H,1)
        SLPS(2) = SLPS(1)
        IF  (SLPS(1) .GE. 0D0)  RETURN
C       ... Split work space and finish initialization
        NX = 1
        NG = N + 1
        FI0 = F
        SL0 = 5D-2 * SLPS(1)
        SLTHR = .995D0 * SLPS(1)
        OK = .FALSE.
        STOP = .FALSE.
        XFD(1,1) = 0D0
        XFD(2,1) = F
        XFD(3,1) = SLPS(1)
        NEV = 0
        B = 1D0
  10    CONTINUE
C       ... Evaluate at  x + b*h
        XFD(1,2) = B
        CALL DCOPY(N, X,1, W,1)
        CALL DAXPY(N, B, H,1, W,1)

c        XTMP=X
c        WRITE(6,'(A,2F8.3)') 'X1a =',X
        CALL FDF(N,W,W(NG),XFD(2,2),GRAD,GRSTEP,RHO)
c        WRITE(6,'(A,2F8.3)') 'X1b =',X
c        X=XTMP
c        WRITE(6,'(A,2F8.3)') 'X1c =',X
        NEV = NEV + 1
        XFD(3,2) = DDOT(N, W(NG),1, H,1)
        IF  (B .EQ. 1D0)  SLPS(2) = XFD(3,2)
        IF  (XFD(2,2) .LE. FI0 + SL0*XFD(1,2))  THEN
C         ... New lower bound
          IF  (XFD(3,2) .LE. ABS(SLTHR))  THEN
            OK = .TRUE.
            ALPHA = XFD(1,2)
            FN = XFD(2,2)
            SLPS(2) = XFD(3,2)
            CALL DCOPY(N, W(NG),1, G,1)
            IF  ((B .LT. 2D0) .AND. (XFD(3,2) .LT. SLTHR))  THEN
C             ... Expand
              CALL DCOPY(3, XFD(1,2),1, XFD(1,1),1)
              B = 2D0
              GOTO 10
            ENDIF
          ENDIF
        ENDIF
C
        D = XFD(1,2) - XFD(1,1)
C
  20    IF  (OK .OR. (NEV .EQ. MEVAL))  RETURN
C
C       ... Refine interval.  Min of quadratic interpolator
        C = XFD(2,2) - XFD(2,1) - D*XFD(3,1)
        IF  (C .GT. 1D-15*DBLE(n)*XFD(1,2))  THEN
C         ... Minimizer in interval
          A = XFD(1,1) - .5D0 * XFD(3,1) * (D**2 / C)
          D = .1D0 * D
          XFD(1,3) = MIN(MAX(XFD(1,1)+D,A), XFD(1,2)-D)
        ELSE
          XFD(1,3) = .5D0 * (XFD(1,1) + XFD(1,2))
        ENDIF
        CALL DCOPY(N, X,1, W,1)
        CALL DAXPY(N, XFD(1,3), H,1, W,1)
c        XTMP=X
c        WRITE(6,'(A,2F8.3)') 'X2a =',X
        CALL FDF(N,W,W(NG),XFD(2,3),GRAD,GRSTEP,RHO)
c        WRITE(6,'(A,2F8.3)') 'X2b =',X
c        X=XTMP
        NEV = NEV + 1
        XFD(3,3) = DDOT(N, W(NG),1, H,1)
        IF  (XFD(2,3) .LT. FI0 + SL0*XFD(1,3))  THEN
C         ... New lower bound
          OK = .TRUE.
          ALPHA = XFD(1,3)
          FN = XFD(2,3)
          SLPS(2) = XFD(3,3)
          CALL DCOPY(N, W(NG),1, G,1)
          CALL DCOPY(3, XFD(1,3),1, XFD(1,1),1)
        ELSE
          CALL DCOPY(3, XFD(1,3),1, XFD(1,2),1)
        ENDIF
C       ... Check convergence
        D = XFD(1,2) - XFD(1,1)
        OK = OK .AND. (ABS(XFD(3,3)) .LE. ABS(SLTHR))
        OK = OK .OR. (D .LE. 0D0)
        GOTO 20
C
***************************  End of  SLINE  ****************************
      END
C
      SUBROUTINE SPCHOL(N,A,FAIL)
************************************************************************
* Cholesky factorization of symmetric matrix given in lower triangle,
* packed form.
* FAIL = 0:  The matrix is positive definite and the Cholesky factor
*            has overwritten  A
* FAIL > 0:  The leading minor of order FAIL is not positive definite
*
* Hans Bruun Nielsen, IMM, DTU.  00.12.18
************************************************************************
      IMPLICIT          NONE
C  Parameters
      INTEGER           N,FAIL
      DOUBLE PRECISION  A(*)
C  Local variables
      INTEGER           K,KK,KN,NK
      INTRINSIC         SQRT
C
C If BLAS is available, then remove the * in col. 1 of the next lines
C and delete from line ??? to the end of the file
      EXTERNAL         DSCAL, DSPR
C
        FAIL = 0
        KK = 1
        DO  K = 1, N
C         ... Test for pos def
          IF  (A(KK) .LE. 0D0)  THEN
            FAIL = K
            RETURN
          ENDIF
          A(KK) = SQRT(A(KK))
          IF  (K .LT. N)  THEN
C           ... Compute k'th column and update trailing submatrix
            NK = N - K
            CALL DSCAL(NK, 1D0/A(KK), A(KK+1),1)
            KN = KK + NK + 1
            CALL DSPR('L',NK, -1D0,A(KK+1),1, A(KN))
            KK = KN
          ENDIF
        END DO
        RETURN
***************************  End of  SPCHOL  ***************************
      END
C
      SUBROUTINE CHKDFN(N,X,STEPL,DIFF,INDX,G,G1,FAIL,GRAD,GRSTEP,RHO)
************************************************************************
* Check implementation of gradient of function of N variables
* Hans Bruun Nielsen, IMM, DTU.  00.09.29
************************************************************************
      IMPLICIT          NONE
      INTEGER           N,INDX(3),FAIL, I,GRAD, RHO(*)
      DOUBLE PRECISION  X(N),STEPL,DIFF(4),G(N),G1(N),
     &                  F,F1, XI,H,AF,AB,AE,ER, GRSTEP(2)
C      EXTERNAL          FDF
      INTRINSIC         ABS, MAX
C       ... Initialize
        FAIL = 1
        DO  I = 1, 4
          DIFF(I) = 0D0
        END DO
        DO  I = 1, 3
          INDX(I) = 0
        END DO
        CALL FDF(N,X,G,F,GRAD,GRSTEP,RHO)
C       ... Run through components of X
        DO  I = 1, N
          DIFF(1) = MAX(DIFF(1), ABS(G(I)))
          XI = X(I)
C         ... Forward
          X(I) = XI + STEPL
          H = X(I) - XI
          IF  (H .EQ. 0D0)  RETURN
          CALL FDF(N,X,G1,F1,GRAD,GRSTEP,RHO)
          AF = (F1 - F)/H
          ER = AF - G(I)
          IF  (ABS(ER) .GT. ABS(DIFF(2)))  THEN
            DIFF(2) = ER
            INDX(1) = I
          ENDIF
C         ... Back
          X(I) = XI - .5D0 * STEPL
          H = X(I) - XI
          IF  (H .EQ. 0D0)  RETURN
          CALL FDF(N,X,G1,F1,GRAD,GRSTEP,RHO)
          AB = (F1 - F)/H
          ER = AB - G(I)
          IF  (ABS(ER) .GT. ABS(DIFF(3)))  THEN
            DIFF(3) = ER
            INDX(2) = I
          ENDIF
C         ... Extrapolated
          AE = (2D0*AB + AF)/3D0
          ER = AE - G(I)
          IF  (ABS(ER) .GT. ABS(DIFF(4)))  THEN
            DIFF(4) = ER
            INDX(3) = I
          ENDIF
C         ... Restore x(i)
          X(I) = XI
       END DO
        FAIL = 0
        RETURN
**************************  end of CHKDFN   ****************************
       END
c$$$
c$$$      SUBROUTINE PRVCTR(NAME,X,I1,I2,UNT)
c$$$************************************************************************
c$$$*  Print on UNT elements  I1 to I2  of  X  with name  NAME
c$$$*  Hans Bruun Nielsen, Numerisk Institut, DTH.  89.09.28.
c$$$************************************************************************
c$$$*  Modified for printing in R. Stig B. Mortensen, Dec. 2008.
c$$$************************************************************************
c$$$
c$$$      IMPLICIT NONE
c$$$      CHARACTER*3      NAME
c$$$      INTEGER          I1,I2,J,J1,J2
c$$$      DOUBLE PRECISION X(*)
c$$$      CHARACTER        UNT*80
c$$$      EXTERNAL         DBLEPR
c$$$C
c$$$        J2 = I1 - 1
c$$$  10    IF  (J2 .GE. I2)  RETURN
c$$$        J1 = J2 + 1
c$$$        J2 = J2 + 5
c$$$        IF  (J2 .GT. I2)  J2 = I2
c$$$        WRITE(UNT,'(1X,A,A,I3,A,I3,A,T18,1P1D10.3,1P4D12.3)')
c$$$     /        NAME,'(',J1,'..',J2,') =',(X(J), J=J1,J2)
c$$$        CALL DBLEPR (UNT, -1, X(1), 0)
c$$$        GOTO 10
c$$$**************************  end of  PRVCTR  ****************************
c$$$      END
c$$$


c ------------------------------------------------------------------------------
      SUBROUTINE FDF(N, X, G, F, GRAD, GRSTEP, RHO)
      IMPLICIT NONE
      INTEGER N, GRAD, RHO(*)
      DOUBLE PRECISION X(N), G(N), F, GRSTEP(2)
      CALL FUNC(N, X, F, RHO)
      IF(GRAD==0) THEN
         CALL USRGR(N, X, G, RHO)
      ELSE
         CALL GR(N, X, F, G, GRAD, GRSTEP, RHO)
      ENDIF
      END


c ------------------------------------------------------------------------------
c$$$      SUBROUTINE FUNC(N, X, VALUE)
c$$$      IMPLICIT NONE
c$$$      INTEGER N
c$$$      DOUBLE PRECISION X(N), VALUE
c$$$      CALL CFUNC(N, X, VALUE)
c$$$      END
c$$$
c$$$      SUBROUTINE USRGR(N, X, G)
c$$$      IMPLICIT NONE
c$$$      INTEGER N
c$$$      DOUBLE PRECISION X(N), G(N)
c$$$      CALL CGRAD(N, X, G)
c$$$      END

      SUBROUTINE GR(N, X, F, G, GRAD, GRSTEP, RHO)
      IMPLICIT NONE
      INTEGER N,I,J, GRAD, RHO(*)
      LOGICAL FWDIFF
      DOUBLE PRECISION X(N), G(N), F, DX, X2(N), F2, F3, GRSTEP(2)
      INTRINSIC ABS
!---- FLAG FOR FORWARD OR CENTRAL DIFF
      FWDIFF = (GRAD==1)
      DO I=1,N
         DO J=1,N
            X2(J) = X(J)
         END DO
         DX = ABS(X2(I)) * GRSTEP(1) + GRSTEP(2)
         X2(I) = X2(I) + DX
         CALL FUNC(N, X2, F2, RHO)
         IF(FWDIFF) THEN
            G(I) = (F2-F)/DX
         ELSE
            X2(I) = X2(I) - 2*DX
            CALL FUNC(N, X2, F3, RHO)
            G(I) = (F2-F3)/(2*DX)
         ENDIF
      END DO
      END

