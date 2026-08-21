      SUBROUTINE BOBYQA (N,NPT,X,XL,XU,RHOBEG,RHOEND,IPRINT,
     1 MAXFUN,W,IERR)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION X(*),XL(*),XU(*),W(*)
C
C     This subroutine seeks the least value of a function of many variables,
C     by applying a trust region method that forms quadratic models by
C     interpolation. There is usually some freedom in the interpolation
C     conditions, which is taken up by minimizing the Frobenius norm of
C     the change to the second derivative of the model, beginning with the
C     zero matrix. The values of the variables are constrained by upper and
C     lower bounds. The arguments of the subroutine are as follows.
C
C     N must be set to the number of variables and must be at least two.
C     NPT is the number of interpolation conditions. Its value must be in
C       the interval [N+2,(N+1)(N+2)/2]. Choices that exceed 2*N+1 are not
C       recommended.
C     Initial values of the variables must be set in X(1),X(2),...,X(N). They
C       will be changed to the values that give the least calculated F.
C     For I=1,2,...,N, XL(I) and XU(I) must provide the lower and upper
C       bounds, respectively, on X(I). The construction of quadratic models
C       requires XL(I) to be strictly less than XU(I) for each I. Further,
C       the contribution to a model from changes to the I-th variable is
C       damaged severely by rounding errors if XU(I)-XL(I) is too small.
C     RHOBEG and RHOEND must be set to the initial and final values of a trust
C       region radius, so both must be positive with RHOEND no greater than
C       RHOBEG. Typically, RHOBEG should be about one tenth of the greatest
C       expected change to a variable, while RHOEND should indicate the
C       accuracy that is required in the final values of the variables. An
C       error return occurs if any of the differences XU(I)-XL(I), I=1,...,N,
C       is less than 2*RHOBEG.
C     The value of IPRINT should be set to 0, 1, 2 or 3, which controls the
C       amount of printing. Specifically, there is no output if IPRINT=0 and
C       there is output only at the return if IPRINT=1. Otherwise, each new
C       value of RHO is printed, with the best vector of variables so far and
C       the corresponding value of the objective function. Further, each new
C       value of F with its variables are output if IPRINT=3.
C     MAXFUN must be set to an upper bound on the number of calls of CALFUN.
C     The array W will be used for working space. Its length must be at least
C       (NPT+5)*(NPT+N)+3*N*(N+5)/2.
CJN   Add IERR to tell what the error is to minqer.
C     IERR gives an error code that is interpreted by minqer interface routine.
C       and passed back to minqa.R
C
C     DOUBLE PRECISION FUNCTION CALFUN (N,X,IP) has to be provided by
C     the user. It returns the value of the objective function for
C     the current values of the variables X(1),X(2),...,X(N), which are
C     generated automatically in a way that satisfies the bounds given
C     in XL and XU.
C
C     Return if the value of NPT is unacceptable.
C
C  Modified by John Nash to put the setup controls into the R code.
      NP=N+1
CJ  Comment out to END IF as in R code
      IF (NPT .LT. N+2 .OR. NPT .GT. ((N+2)*NP)/2) THEN
CJN         CALL minqer(10)
c$$$          PRINT 10
c$$$   10     FORMAT (/4X,'Return from BOBYQA because NPT is not in',
c$$$     1      ' the required interval')
c$$$          GO TO 40
          IERR = 10
CJN    Fail out NPT has unacceptable value
          GO TO 40
      END IF
C
C     Partition the working space array, so that different parts of it can
C     be treated separately during the calculation of BOBYQB. The partition
C     requires the first (NPT+2)*(NPT+N)+3*N*(N+5)/2 elements of W plus the
C     space that is taken by the last array in the argument list of BOBYQB.
C
      NDIM=NPT+N
      IXB=1
      IXP=IXB+N
      IFV=IXP+N*NPT
      IXO=IFV+NPT
      IGO=IXO+N
      IHQ=IGO+N
      IPQ=IHQ+(N*NP)/2
      IBMAT=IPQ+NPT
      IZMAT=IBMAT+NDIM*N
      ISL=IZMAT+NPT*(NPT-NP)
      ISU=ISL+N
      IXN=ISU+N
      IXA=IXN+N
      ID=IXA+N
      IVL=ID+N
      IW=IVL+NDIM
CJN 100807 to ensure a defined value for the return code
      IERR = 0
C
C     Return if there is insufficient space between the bounds. Modify the
C     initial X if necessary in order to avoid conflicts between the bounds
C     and the construction of the first quadratic model. The lower and upper
C     bounds on moves from the updated X are set now, in the ISL and ISU
C     partitions of W, in order to provide useful and exact information about
C     components of X that become within distance RHOBEG from their bounds.
C
      ZERO=0.0D0
      DO J=1,N
         TEMP=XU(J)-XL(J)
         IF (TEMP .LT. RHOBEG+RHOBEG) THEN
C     JN         CALL minqer(20)
c$$$  PRINT 20
c$$$  20     FORMAT (/4X,'Return from BOBYQA because one of the',
c$$$  1      ' differences XU(I)-XL(I)'/6X,' is less than 2*RHOBEG.')
c$$$  GO TO 40
            IERR = 20
            GOTO 40
         END IF
         JSL=ISL+J-1
         JSU=JSL+N
         W(JSL)=XL(J)-X(J)
         W(JSU)=XU(J)-X(J)
         IF (W(JSL) .GE. -RHOBEG) THEN
            IF (W(JSL) .GE. ZERO) THEN
               X(J)=XL(J)
               W(JSL)=ZERO
               W(JSU)=TEMP
            ELSE
               X(J)=XL(J)+RHOBEG
               W(JSL)=-RHOBEG
               W(JSU)=DMAX1(XU(J)-X(J),RHOBEG)
            END IF
         ELSE IF (W(JSU) .LE. RHOBEG) THEN
            IF (W(JSU) .LE. ZERO) THEN
               X(J)=XU(J)
               W(JSL)=-TEMP
               W(JSU)=ZERO
            ELSE
               X(J)=XU(J)-RHOBEG
               W(JSL)=DMIN1(XL(J)-X(J),-RHOBEG)
               W(JSU)=RHOBEG
         END IF
      END IF
      END DO
C     
C     Make the call of BOBYQB.
C
      CALL BOBYQB (N,NPT,X,XL,XU,RHOBEG,RHOEND,IPRINT,MAXFUN,W(IXB),
     1  W(IXP),W(IFV),W(IXO),W(IGO),W(IHQ),W(IPQ),W(IBMAT),W(IZMAT),
     2  NDIM,W(ISL),W(ISU),W(IXN),W(IXA),W(ID),W(IVL),W(IW), IERR)
 40   RETURN
CJN Added 100807 -- put label 40     
      END

