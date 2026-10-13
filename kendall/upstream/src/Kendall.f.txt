C TAUK2.
C CALCULATE KENDALL'S TAU AND ITS P-VALUE FOR THE GENERAL
C CASE OF 2 VARIATES.
C KENDALL'S METHOD IS USED FOR TIES.
C
C REFERENCE: M.G. KENDALL, "RANK CORRELATION METHODS" PUBLISHED BY
C            GRIFFIN & CO.
C
      Subroutine TAUK2(X,Y,N,TAU,PROB,SLTAU,S,VARS,DENOM,IW,IER)
C
      Dimension X(N), Y(N), IW(N)
      Logical SW, SWX, SWY, TIES
      Integer SCOREK
      IER = 0
      TIES = .FALSE.
      PROB = 1.0
      SLTAU = 1.0
      TAU = 1.0
      If (N.LT.2) Then
       IER = 10
       return
      endif
       SWX = .TRUE.
       SWY = .TRUE.
       Do I = 2, N
        If (X(I).NE.X(I-1)) SWX = .FALSE.
        If (Y(I).NE.Y(I-1)) SWY = .FALSE.
       Enddo
c sw is true if at least one of X or Y has no ties
       SW = SWX .OR. SWY
       If (.NOT.SW) Then
        NM1 = N - 1
        IS = 0
        Do I = 1, NM1
         IP1 = I + 1
         Do J = IP1, N
          IS = IS + SCOREK(X(I),Y(I),X(J),Y(J))
         Enddo
        Enddo
        S = FLOAT(IS)
        CN = FLOAT(N) * FLOAT(N-1)
        DN = 0.5 * CN
        Do I = 1, N
         IW(I) = 0
        Enddo
        SUMT = 0.0
        SUMA1 = 0.0
        SUMA2 = 0.0
        Do I = 1, NM1
         IP1 = I + 1
         TMP = 1.0
         Do J = IP1, N
          If (X(I).EQ.X(J)) Then
           If (IW(J).NE.1) Then
            IW(J) = 1
            TMP = TMP + 1.0
            TIES = .TRUE.
           Endif
          Endif
         Enddo
         SUMT = SUMT + TMP * (TMP-1.0)
         SUMA1 = SUMA1 + TMP * (TMP-1.0) * (2.0*TMP+5.0)
         SUMA2 = SUMA2 + TMP * (TMP-1.0) * (TMP-2.0)
        Enddo
        D1 = DN - 0.5 * SUMT
        D1 = SQRT(D1)
        Do I = 1, N
         IW(I) = 0
        Enddo
        SUMU = 0.0
        SUMB1 = 0.0
        SUMB2 = 0.0
        Do I = 1, NM1
         IP1 = I + 1
         TMP = 1.0
         Do J = IP1, N
          If (Y(I).EQ.Y(J)) Then
           If (IW(J).NE.1) Then
            IW(J) = 1
            TMP = TMP + 1.0
            TIES = .TRUE.
           Endif
          Endif
         Enddo
         SUMU = SUMU + TMP * (TMP-1.0)
         SUMB1 = SUMB1 + TMP * (TMP-1.0) * (2.0*TMP+5.0)
         SUMB2 = SUMB2 + TMP * (TMP-1.0) * (TMP-2.0)
        Enddo
        D2 = DN - 0.5 * SUMU
        D2 = SQRT(D2)
        If ((D1.LE.0.0).OR.(D2.LE.0.0)) Then
		ier = 30
		return
	endif
         DENOM = D1*D2
         TAU = S / DENOM
         VARS1 = (CN*FLOAT(2*N+5)-SUMA1-SUMB1) / 18.0
         VARS2 = SUMA2 * SUMB2 / (9.0*CN*FLOAT(N-2))
         VARS3 = SUMT * SUMU / (CN+CN)
         VARS = VARS1 + VARS2 + VARS3
         SDS = SQRT(VARS)
C
C USE EXACT METHOD IF THERE ARE NO TIES.
C
         If (sw) Then
C
C CALCULATE P-VALUE USING EXACT METHOD
C
          IS = IFIX(S)
          PROB = 1.0 - PRTAUS(IS,N,IER1)
          If (IER1.EQ.0) Then
           SLTAU = 2.0 * PROB
           If (PROB.GE.0.5) SLTAU = 2.0 * (1.0-PROB)
           Return
          Endif
C
C IF THERE ARE TIES, NEED ATLEAST SAMPLE SIZE OF ATLEAST 3
C
         Elseif (N.GT.3) Then
C
C USE CONTINUITY CORRECTION FOR S
          SCOR = 0.0
          If (S.GT.0) SCOR = S - 1.0
          If (S.LT.0) SCOR = S + 1.0
C CALCULATE P-VALUE USING NORMAL APPROXIMATION
          ZSCORE = SCOR / SDS
          PROB = ALNORM(ZSCORE,.FALSE.)
          SLTAU = 2.0 * PROB
          If (PROB.GE.0.5) SLTAU = 2.0 * (1.0-PROB)
          Return
       Endif
      Endif
C THIS ONLY HAPPENS WHEN N<4. aim. 06/07/2009
      IER = 12
      Return
      End
      Integer Function SCOREK(X1,Y1,X2,Y2)
      SCOREK = 1
      If ((X1.GT.X2.AND.Y1.LT.Y2).OR.(X1.LT.X2.AND.Y1.GT.Y2)) Then
       SCOREK = -1
      Elseif ((X1.EQ.X2).OR.(Y1.EQ.Y2)) Then
       SCOREK = 0
      Endif
      Return
      End
      Function ALNORM(X,UPPER)
C
C   ALGORITHM AS 66 APPL. STATIST. (1973) VOL.22, NO.3
C
C   EVALUATES THE TAIL AREA OF THE STANDARDISED NORMAL
C   CURVE FROM X TO INFINITY IF UPPER IS .TRUE. OR
C   FROM MINUS INFINITY TO X IF UPPER IS .FALSE.
C
      Real LTONE, UTZERO, ZERO, HALF, ONE, CON, Z, Y, X
      Logical UPPER, UP
C
C   LTONE AND UTZERO MUST BE SET TO SUIT THE PARTICULAR
C   COMPUTER (SEE INTRODUCTORY TEXT)
C
      Data LTONE, UTZERO /7.0, 18.66/
      Data ZERO, HALF, ONE, CON /0.0, 0.5, 1.0, 1.28/
C
      ALNORM = ONE
      If (X.GT.ZERO) ALNORM = ZERO
      If (ABS(X).LT.10.0) Then
C
       UP = UPPER
       Z = X
       If (Z.LT.ZERO) Then
        UP = .NOT. UP
        Z = -Z
       Endif
       If (Z.LE.LTONE.OR.UP.AND.Z.LE.UTZERO) Then
        Y = HALF * Z * Z
        If (Z.GT.CON) Then
C
         ALNORM = 0.398942280385 * EXP(-Y) / (Z-3.8052E-8+1.00000615302/
     *       (Z+3.98064794E-4+1.98615381364/(Z-0.151679116635+
     *       5.29330324926/(Z+4.8385912808-15.1508972451/(Z+
     *       0.742380924027+30.789933034/(Z+3.99019417011))))))
        Else
C
         ALNORM = HALF - Z * (0.398942280444-0.399903438504*Y/(Y+
     *       5.75885480458-29.8213557808/(Y+2.62433121679+48.6959930692/
     *       (Y+5.92885724438))))
        Endif
       Else
        ALNORM = ZERO
       Endif
C
       If (.NOT.UP) ALNORM = ONE - ALNORM
       Return
      Else
       If (UPPER) Return
       ALNORM = 1.0 - ALNORM
       Return
      Endif
      End
      Function PRTAUS(IS,N,IFAULT)
C
C        ALGORITHM AS 71  APPL. STATIST. (1974) VOL.23, NO.1
C
C        GIVEN A VALUE OF IS CALCULATED FROM TWO RANKINGS (WITHOUT TIES)
C        OF N OBJECTS, THE FUNCTION COMPUTES THE PROBABILITY OF
C        OBTAINING A VALUE GREATER THAN, OR EQUAL TO, IS.
C
      Dimension H(15), L(2,15)
      Equivalence (H(1),L(1,1))
C
C        CHECK ON THE VALIDITY OF IS AND N VALUES
C
      PRTAUS = 1.0
      IFAULT = 2
      If (N.LT.1) Return
      IFAULT = 1
      M = N * (N-1) / 2 - IABS(IS)
      If (M.LT.0.OR.M.GT.(M/2)*2) Return
      IFAULT = 0
      If (M.EQ.0.AND.IS.LE.0) Return
      If (N.GT.8) Then
C
C        CALCULATION OF TCHEBYCHEFF-HERMITE POLYNOMIALS
C
       X = FLOAT(IS-1) / SQRT((FLOAT(6+N*(5-N*(3+2*N))))/(-18.0))
       H(1) = X
       H(2) = X * X - 1.0
       Do I = 3, 15
        H(I) = X * H(I-1) - FLOAT(I-1) * H(I-2)
       Enddo
C
C        PROBABILITIES CALCULATED BY MODIFIED EDGEWORTH SERIES FOR
C        N GREATER THAN 8
C
       R = 1.0 / FLOAT(N)
       SC = R * (H(3)*(-9.0000E-2+R*(4.5000E-2+R*(-5.325E-1+R*5.06E-1)))
     *     +R*(H(5)*(3.6735E-2+R*(-3.6735E-2+R*3.214E-1))+H(7)*(
     *     4.0500E-3+R*(-2.3336E-2+R*7.787E-2))+R*(H(9)*(-3.3061E-3-R*
     *     6.5166E-3)+H(11)*(-1.2150E-4+R*2.5927E-3)+R*(H(13)*1.4878E-4+
     *     H(15)*2.7338E-6))))
C
C        CALL TO ALGORITHM AS 66
C
       PRTAUS = ALNORM(X,.TRUE.) + SC * 0.398942 * EXP(-0.5*X*X)
       If (PRTAUS.LT.0.0) PRTAUS = 0.0
       If (PRTAUS.GT.1.0) PRTAUS = 1.0
       Return
      Else
C
C        PROBABILITIES CALCULATED BY RECURRENCE RELATION FOR
C     N LESS THAN 9
C
       If (IS.LT.0) M = M - 2
       IM = M / 2 + 1
       L(1,1) = 1
       L(2,1) = 1
       If (IM.GE.2) Then
        Do I = 2, IM
         L(1,I) = 0
         L(2,I) = 0
        Enddo
       Endif
       IL = 1
       I = 1
       M = 1
       J = 1
       JJ = 2
       Do While (I.NE.N)
        IL = IL + I
        I = I + 1
        M = M * I
        J = 3 - J
        JJ = 3 - JJ
        IN = 1
        IO = 0
        K = MIN0(IM,IL)
        Do While (.TRUE.)
         IN = IN + 1
         If (IN.GT.K) Goto 10
         L(JJ,IN) = L(JJ,IN-1) + L(J,IN)
         If (IN.GT.I) Then
          IO = IO + 1
          L(JJ,IN) = L(JJ,IN) - L(J,IO)
         Endif
        Enddo
   10  Continue
       Enddo
       K = 0
       Do I = 1, IM
        K = K + L(JJ,I)
       Enddo
       PRTAUS = FLOAT(K) / FLOAT(M)
       If (IS.LT.0) PRTAUS = 1.0 - PRTAUS
       Return
      Endif
      End

