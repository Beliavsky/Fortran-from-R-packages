

       SUBROUTINE BLOCK(N,TOPBLK,NRWTOP,NOVRLP,ARRAY,NRWBLK,
     *         NCLBLK,NBLOKS,BOTBLK,NRWBOT,PIVOT,rB,cB,B,X,IFLAG,
     *         TMP1, TMP2)

        DOUBLE PRECISION TOPBLK,ARRAY,BOTBLK,B,X,TMP1,TMP2
        INTEGER PIVOT(*),nrwtop,NRWBLK,NCLBLK,NRWBOT,N,NOVRLP,NBLOKS
        INTEGER I, J, IFLAG, rB, cB
        DIMENSION TOPBLK(NRWTOP,*),ARRAY(NRWBLK,NCLBLK,*),
     *          BOTBLK(NRWBOT,*),B(rB,cB),X(N,cB),TMP1(*),TMP2(*)


         DO J = 1,N
           TMP1(j) = B(J,1)
         ENDDO
         CALL COLROW(N,TOPBLK,NRWTOP,NOVRLP,ARRAY,NRWBLK,NCLBLK,NBLOKS,
     *          BOTBLK,NRWBOT,PIVOT,TMP1,TMP2,IFLAG)
         DO J = 1,N
           X(J,1) = TMP2(j)
         ENDDO

        IF (IFLAG.NE.0)RETURN
        DO I = 2,cB

         DO J = 1,N
           TMP1(j) = B(J,I)
         ENDDO

         CALL CRSLVE(TOPBLK,NRWTOP,NOVRLP,ARRAY,NRWBLK,NCLBLK,NBLOKS,
     *          BOTBLK,NRWBOT,PIVOT,TMP1,TMP2)
         DO J = 1,N
           X(J,I) = TMP2(j)
         ENDDO

        ENDDO
        
        RETURN
        END

       SUBROUTINE COLROW(N,TOPBLK,NRWTOP,NOVRLP,ARRAY,NRWBLK,
     *             NCLBLK,NBLOKS,BOTBLK,NRWBOT,PIVOT,B,X,IFLAG)

C
C***************************************************************
C
C  THIS PROGRAM SOLVES THE LINEAR SYSTEM  A*X = B  WHERE  A IS
C  AN ALMOST BLOCK DIAGONAL MATRIX OF THE FORM
C
C               TOPBLK
C               ARRAY(1)
C                     ARRAY(2)
C                          .
C                             .
C                                .
C                                   .
C                                    ARRAY(NBLOKS)
C                                           BOTBLK
C
C  WHERE
C           TOPBLK IS  NRWTOP  BY NOVRLP
C           ARRAY(K), K=1,NBLOKS, ARE NRWBLK BY NRWBLK+NOVRLP
C           BOTBLK IS NRWBOT BY NOVRLP,
C  AND
C           NOVRLP = NRWTOP + NRWBOT
C  WITH
C           NOVRLP.LE.NRWBLK .
C
C  THE LINEAR SYSTEM IS OF ORDER  N = NBLOKS*NRWBLK + NOVRLP.
C
C  THE METHOD IMPLEMENTED IS BASED ON GAUSS ELIMINATION WITH
C  ALTERNATE ROW AND COLUMN ELIMINATION WITH PARTIAL PIVOTING,
C  WHICH PRODUCES A STABLE DECOMPOSITION OF THE MATRIX  A
C  WITHOUT INTRODUCING FILL-IN.
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C               *****  PARAMETERS  *****
C
C       *** ON ENTRY ...
C
C               N      - INTEGER
C                         THE ORDER OF THE LINEAR SYSTEM,
C                         GIVEN BY NBLOKS*NRWBLK + NOVRLP
C
C               TOPBLK - DOUBLE PRECISION(NRWTOP,NOVRLP)
C                         THE FIRST BLOCK OF THE ALMOST BLOCK
C                         DIAGONAL MATRIX A
C
C               NRWTOP - INTEGER
C                         NUMBER OF ROWS IN THE BLOCK TOPBLK
C
C               NOVRLP - INTEGER
C                         THE NUMBER OF COLUMNS IN WHICH SUCC-
C                         ESSIVE BLOCKS OVERLAP, WHERE
C                                NOVRLP = NRWTOP + NRWBOT
C
C               ARRAY  - DOUBLE PRECISION(NRWBLK,NCLBLK,NBLOKS)
C                         ARRAY(,,K) CONTAINS THE K-TH NRWBLK
C                         BY NCLBLK BLOCK OF THE MATRIX A
C
C               NRWBLK - INTEGER
C                         NUMBER OF ROWS IN K-TH BLOCK
C
C               NCLBLK - INTEGER
C                         NUMBER OF COLUMNS IN K-TH BLOCK
C
C               NBLOKS - INTEGER
C                         NUMBER OF NRWBLK BY NCLBLK BLOCKS IN
C                         THE MATRIX A
C
C               BOTBLK - DOUBLE PRECISION(NRWBOT,NOVRLP)
C                         THE LAST BLOCK OF THE MATRIX A
C
C               NRWBOT - INTEGER
C                         NUMBER OF ROWS IN THE BLOCK BOTBLK
C
C                PIVOT - INTEGER(N)
C                         WORK SPACE
C
C                    B - DOUBLE PRECISION(N)
C                         THE RIGHT HAND SIDE VECTOR
C
C                    X - DOUBLE PRECISION(N)
C                         WORK SPACE
C
C       *** ON RETURN  ...
C
C               TOPBLK,ARRAY,BOTBLK - ARRAYS CONTAINING THE
C                        DESIRED DECOMPOSITION OF THE MATRIX A
C                        (IF IFLAG = 0)
C
C                PIVOT - INTEGER(N)
C                         RECORDS THE PIVOTING INDICES DETER-
C                         MINED IN THE DECOMPOSITION
C
C                    X - DOUBLE PRECISION(N)
C                         THE SOLUTION VECTOR (IF IFLAG = 0)
C
C               IFLAG  - INTEGER
C                         =  1, IF INPUT PARAMETERS ARE INVALID
C                         = -1, IF MATRIX IS SINGULAR
C                         =  0, OTHERWISE
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C               *****  AUXILIARY PROGRAMS  *****
C
C       CRDCMP(N,TOPBLK,NRWTOP,NOVRLP,ARRAY,NRWBLK,NCLBLK,NBLOKS,
C    *     BOTBLK,NRWBOT,PIVOT,IFLAG)
C            - DECOMPOSES THE MATRIX  A  USING MODIFIED
C              ALTERNATE ROW AND COLUMN ELIMINATON WITH
C              PARTIAL PIVOTING, AND IS USED FOR THIS
C              PURPOSE IN C O L R O W.
C              THE ARGUMENTS ARE AS IN C O L R O W.
C
C       CRSLVE(N,TOPBLK,NRWTOP,NOVRLP,ARRAY,NRWBLK,NCLBLK,NBLOKS,
C    *     BOTBLK,NRWBOT,PIVOT,B,X)
C            - SOLVES THE SYSTEM A*X = B ONCE A IS DECOMPOSED.
C              THE ARGUMENTS ARE ALLAS IN C O L R O W.
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C       THE SUBROUTINE  C O L R O W  AUTOMATICALLY SOLVES THE
C  INPUT SYSTEM WHEN IFLAG=0.  C O L R O W  IS CALLED ONLY ONCE
C  FOR A GIVEN SYSTEM. THE SOLUTION FOR A SEQUENCE OF P RIGHT
C  HAND SIDES CAN BE OBTAINED BY ONE CALL TO  C O L R O W  AND
C  P-1 CALLS TO CRSLVE ONLY. SINCE THE ARRAYS TOPBLK,ARRAY,
C  BOTBLK AND PIVOT CONTAIN THE DECOMPOSITION OF THE GIVEN
C  COEFFICIENT MATRIX AND PIVOTING INFORMATION ON RETURN FROM
C  C O L R O W , THEY MUST NOT BE ALTERED BETWEEN SUCCESSIVE
C  CALLS TO CRSLVE WITH THE SAME LEFT HAND SIDES. FOR THE
C  SAME REASON, IF THE USER WISHES TO SAVE THE COEFFICIENT
C  MATRIX, THE ARRAYS TOPBLK,ARRAY,BOTBLK MUST BE COPIED
C  BEFORE A CALL TO  C O L R O W .
C
C*************************************************************************
C
        DOUBLE PRECISION TOPBLK,ARRAY,BOTBLK,B,X
        INTEGER PIVOT(*),nrwtop,NRWBLK,NCLBLK,NRWBOT,N,NOVRLP,NBLOKS
        INTEGER IFLAG
        DIMENSION TOPBLK(NRWTOP,*),ARRAY(NRWBLK,NCLBLK,*),
     *          BOTBLK(NRWBOT,*),B(*),X(*)
        CALL CRDCMP(N,TOPBLK,NRWTOP,NOVRLP,ARRAY,NRWBLK,NCLBLK,NBLOKS,
     *          BOTBLK,NRWBOT,PIVOT,IFLAG)
        IF(IFLAG.NE.0)RETURN
        CALL CRSLVE(TOPBLK,NRWTOP,NOVRLP,ARRAY,NRWBLK,NCLBLK,NBLOKS,
     *          BOTBLK,NRWBOT,PIVOT,B,X)
        RETURN
        END



        SUBROUTINE CRDCMP(N,TOPBLK,NRWTOP,NOVRLP,ARRAY,NRWBLK,
     *             NCLBLK,NBLOKS,BOTBLK,NRWBOT,PIVOT,IFLAG)
C
C***************************************************************
C
C  C R D C M P DECOMPOSES THE ALMOST BLOCK DIAGONAL MATRIX A
C  USING MODIFIED ALTERNATE ROW AND COLUMN ELIMINATION WITH
C  PARTIAL PIVOTING.  THE MATRIX  A  IS STORED IN THE ARRAYS
C  TOPBLK, ARRAY, AND BOTBLK.
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C               *****  PARAMETERS  *****
C
C       *** ON ENTRY ...
C
C               N      - INTEGER
C                         THE ORDER OF THE LINEAR SYSTEM,
C                         GIVEN BY NBLOKS*NRWBLK + NOVRLP
C
C               TOPBLK - DOUBLE PRECISION(NRWTOP,NOVRLP)
C                         THE FIRST BLOCK OF THE ALMOST BLOCK
C                         DIAGONAL MATRIX A TO BE DECOMPOSED
C
C               NRWTOP - INTEGER
C                         NUMBER OF ROWS IN THE BLOCK TOPBLK
C
C               NOVRLP - INTEGER
C                         THE NUMBER OF COLUMNS IN WHICH SUCC-
C                         ESSIVE BLOCKS OVERLAP, WHERE
C                                NOVRLP = NRWTOP + NRWBOT
C
C               ARRAY  - DOUBLE PRECISION(NRWBLK,NCLBLK,NBLOKS)
C                         ARRAY(,,K) CONTAINS THE K-TH NRWBLK
C                         BY NCLBLK BLOCK OF THE MATRIX A
C
C               NRWBLK - INTEGER
C                         NUMBER OF ROWS IN K-TH BLOCK
C
C               NCLBLK - INTEGER
C                         NUMBER OF COLUMNS IN K-TH BLOCK
C
C               NBLOKS - INTEGER
C                         NUMBER OF NRWBLK BY NCLBLK BLOCKS IN
C                         THE MATRIX A
C
C               BOTBLK - DOUBLE PRECISION(NRWBOT,NOVRLP)
C                         THE LAST BLOCK OF THE MATRIX A
C
C               NRWBOT - INTEGER
C                         NUMBER OF ROWS IN THE BLOCK BOTBLK
C
C                PIVOT - INTEGER(N)
C                         WORK SPACE
C
C       *** ON RETURN  ...
C
C               TOPBLK,ARRAY,BOTBLK - ARRAYS CONTAINING THE
C                        DESIRED DECOMPOSITION OF THE MATRIX A
C                        (IF IFLAG = 0)
C
C                PIVOT - INTEGER(N)
C                         RECORDS THE PIVOTING INDICES DETER-
C                         MINED IN THE DECOMPOSITION
C
C               IFLAG  - INTEGER
C                         =  1, IF INPUT PARAMETERS ARE INVALID
C                         = -1, IF MATRIX IS SINGULAR
C                         =  0, OTHERWISE
C
C***************************************************************
C
        DOUBLE PRECISION TOPBLK,ARRAY,BOTBLK
        DOUBLE PRECISION ROWMAX,ROWPIV,ROWMLT,COLMAX,COLPIV
        DOUBLE PRECISION SWAP,COLMLT,PIVMAX,ZERO,TEMPIV
        INTEGER PIVOT(*),NRWTOP,NRWBLK,NCLBLK,NRWBOT,IFLAG
        INTEGER NRWTP1,NROWEL,NRWEL1,NVRLP0,NOVRLP,N,NBLOKS
        INTEGER I,IPLUS1,IPVT,J,L,INCR,K,KPLUS1,JPLUS1,JMINN,LOOP
        INTEGER INCRJ,IPLUSN,INCRN,IRWBLK,IPVBLK,JRWBLK
        DIMENSION TOPBLK(NRWTOP,*),ARRAY(NRWBLK,NCLBLK,*),
     *          BOTBLK(NRWBOT,*)
C        DATA ZERO/0.0D0/
C
C***************************************************************
C
C          ****  DEFINE THE CONSTANTS USED THROUGHOUT  ****
C
C***************************************************************
C
        ZERO  = 0.D0
        IFLAG = 0
        PIVMAX = ZERO
        NRWTP1 = NRWTOP+1
        NROWEL = NRWBLK-NRWTOP
        NRWEL1 = NROWEL+1
        NVRLP0 = NOVRLP-1
C
C***************************************************************
C
C          ****  CHECK VALIDITY OF THE INPUT PARAMETERS....
C
C               IF PARAMETERS ARE INVALID THEN TERMINATE AT 10;
C                                         ELSE CONTINUE AT 100.
C
C***************************************************************
C
        IF(N.NE.NBLOKS*NRWBLK+NOVRLP)GO TO 10
        IF(NOVRLP.NE.NRWTOP+NRWBOT)GO TO 10
        IF(NCLBLK.NE.NOVRLP+NRWBLK)GO TO 10
        IF(NOVRLP.GT.NRWBLK)GO TO 10
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C          PARAMETERS ARE ACCEPTABLE - CONTINUE AT 100.
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
        GO TO 100
10      CONTINUE
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C          PARAMETERS ARE INVALID.  SET IFLAG = 1, AND TERMINATE
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
        IFLAG = 1
        RETURN
100     CONTINUE
C
C***************************************************************
C
C               ****  FIRST, IN TOPBLK....
C
C***************************************************************
C
C          ***  APPLY NRWTOP COLUMN ELIMINATIONS WITH COLUMN
C                 PIVOTING ....
C
C***************************************************************
C
        DO 190 I = 1,NRWTOP
           IPLUS1 = I+1
C
C          CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C               DETERMINE COLUMN PIVOT AND PIVOT INDEX
C
C          CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
           IPVT = I
           COLMAX = DABS(TOPBLK(I,I))
           DO 110 J = IPLUS1,NOVRLP
              TEMPIV = DABS(TOPBLK(I,J))
              IF(TEMPIV.LE.COLMAX)GO TO 110
                 IPVT = J
                 COLMAX = TEMPIV
110        CONTINUE
C
C          CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C               TEST FOR SINGULARITY:
C
C                       IF SINGULAR THEN TERMINATE AT 1000;
C                                   ELSE CONTINUE.
C
C          CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
           IF(PIVMAX+COLMAX.EQ.PIVMAX)GO TO 1000
           PIVMAX = DMAX1(COLMAX,PIVMAX)
C
C          CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C               IF NECESSARY INTERCHANGE COLUMNS
C
C          CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
           PIVOT(I) = IPVT
           IF(IPVT.EQ.I)GO TO 140
              DO 120 L = I,NRWTOP
                 SWAP = TOPBLK(L,IPVT)
                 TOPBLK(L,IPVT) = TOPBLK(L,I)
                 TOPBLK(L,I) = SWAP
120           CONTINUE
              DO 130 L = 1,NRWBLK
                 SWAP = ARRAY(L,IPVT,1)
                 ARRAY(L,IPVT,1) = ARRAY(L,I,1)
                 ARRAY(L,I,1) = SWAP
130           CONTINUE
140        CONTINUE
C
C          CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C               COMPUTE MULTIPLIERS AND PERFORM COLUMN
C                       ELIMINATION
C
C          CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
           COLPIV = TOPBLK(I,I)
           DO 180 J = IPLUS1,NOVRLP
              COLMLT = TOPBLK(I,J)/COLPIV
              TOPBLK(I,J) = COLMLT
              IF(IPLUS1.GT.NRWTOP)GO TO 160
                 DO 150 L = IPLUS1,NRWTOP
                    TOPBLK(L,J) = TOPBLK(L,J)-COLMLT*TOPBLK(L,I)
150              CONTINUE
160           CONTINUE
              DO 170 L = 1,NRWBLK
                 ARRAY(L,J,1) = ARRAY(L,J,1)-COLMLT*ARRAY(L,I,1)
170           CONTINUE
180        CONTINUE
190     CONTINUE
C
C***************************************************************
C
C          ****  IN EACH BLOCK ARRAY(,,K)....
C
C***************************************************************
C
        INCR = 0
        DO 395 K = 1,NBLOKS
           KPLUS1 = K+1
C
C          *****************************************************
C
C          ***  FIRST APPLY NRWBLK-NRWTOP ROW ELIMINATIONS WITH
C                       ROW PIVOTING....
C
C          *****************************************************
C
           DO 270 J = NRWTP1,NRWBLK
              JPLUS1 = J+1
              JMINN = J-NRWTOP
C
C             CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C               DETERMINE ROW PIVOT AND PIVOT INDEX
C
C             CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
              IPVT = JMINN
              ROWMAX = DABS(ARRAY(JMINN,J,K))
              LOOP = JMINN+1
              DO 210 I = LOOP,NRWBLK
                 TEMPIV = DABS(ARRAY(I,J,K))
                 IF(TEMPIV.LE.ROWMAX)GO TO 210
                 IPVT = I
                 ROWMAX = TEMPIV
210           CONTINUE
C
C             CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C               TEST FOR SINGULARITY:
C
C                       IF SINGULAR THEN TERMINATE AT 1000;
C                                   ELSE CONTINUE.
C
C             CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
              IF(PIVMAX+ROWMAX.EQ.PIVMAX)GO TO  1000
              PIVMAX = DMAX1(ROWMAX,PIVMAX)
C
C             CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C               IF NECESSARY INTERCHANGE ROWS
C
C             CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
              INCRJ = INCR+J
              PIVOT(INCRJ) = INCR+IPVT+NRWTOP
              IF(IPVT.EQ.JMINN)GO TO 230
                 DO 220 L = J,NCLBLK
                    SWAP = ARRAY(IPVT,L,K)
                    ARRAY(IPVT,L,K) = ARRAY(JMINN,L,K)
                    ARRAY(JMINN,L,K) = SWAP
220              CONTINUE
230           CONTINUE
C
C             CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C               COMPUTE MULTIPLERS
C
C             CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
              ROWPIV = ARRAY(JMINN,J,K)
              DO 240 I = LOOP,NRWBLK
                 ARRAY(I,J,K) = ARRAY(I,J,K)/ROWPIV
240           CONTINUE
C
C             CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C               PERFORM ROW ELIMINATION WITH COLUMN INDEXING
C
C             CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
              DO 260 L = JPLUS1,NCLBLK
                 ROWMLT = ARRAY(JMINN,L,K)
                 DO 250 I = LOOP,NRWBLK
                    ARRAY(I,L,K) = ARRAY(I,L,K)
     *                                -ROWMLT*ARRAY(I,J,K)
250              CONTINUE
260           CONTINUE
270        CONTINUE
C
C          *****************************************************
C
C          ***  NOW APPLY NRWTOP COLUMN ELIMINATIONS WITH
C                      COLUMN PIVOTING....
C
C          *****************************************************
C
           DO 390 I = NRWEL1,NRWBLK
              IPLUSN = I+NRWTOP
              IPLUS1 = I+1
C
C             CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C               DETERMINE COLUMN PIVOT AND PIVOT INDEX
C
C             CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
              IPVT = IPLUSN
              COLMAX = DABS(ARRAY(I,IPVT,K))
              LOOP = IPLUSN+1
              DO 310 J = LOOP,NCLBLK
                 TEMPIV = DABS(ARRAY(I,J,K))
                 IF(TEMPIV.LE.COLMAX)GO TO 310
                 IPVT = J
                 COLMAX = TEMPIV
310           CONTINUE
C
C             CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C               TEST FOR SINGULARITY:
C
C                       IF SINGULAR THEN TERMINATE AT 1000;
C                                   ELSE CONTINUE.
C
C             CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
              IF(PIVMAX+COLMAX.EQ.PIVMAX)GO TO 1000
              PIVMAX = DMAX1(COLMAX,PIVMAX)
C
C             CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C               IF NECESSARY INTERCHANGE COLUMNS
C
C             CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
              INCRN = INCR+IPLUSN
              PIVOT(INCRN) = INCR+IPVT
              IRWBLK = IPLUSN-NRWBLK
              IF(IPVT.EQ.IPLUSN)GO TO 340
                 DO 315 L = I,NRWBLK
                    SWAP = ARRAY(L,IPVT,K)
                    ARRAY(L,IPVT,K) = ARRAY(L,IPLUSN,K)
                    ARRAY(L,IPLUSN,K) = SWAP
315              CONTINUE
                 IPVBLK = IPVT-NRWBLK
                 IF(K.EQ.NBLOKS)GO TO 330
                    DO 320 L = 1,NRWBLK
                       SWAP = ARRAY(L,IPVBLK,KPLUS1)
                       ARRAY(L,IPVBLK,KPLUS1)
     *                                 = ARRAY(L,IRWBLK,KPLUS1)
                       ARRAY(L,IRWBLK,KPLUS1) = SWAP
320                 CONTINUE
                    GO TO 340
330              CONTINUE
                 DO 335 L = 1,NRWBOT
                    SWAP = BOTBLK(L,IPVBLK)
                    BOTBLK(L,IPVBLK) = BOTBLK(L,IRWBLK)
                    BOTBLK(L,IRWBLK) = SWAP
335              CONTINUE
340           CONTINUE
C
C             CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C               COMPUTE MULTIPLIERS AND PERFORM COLUMN
C                       ELIMINATION
C
C             CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
              COLPIV = ARRAY(I,IPLUSN,K)
              DO 380 J = LOOP,NCLBLK
                 COLMLT = ARRAY(I,J,K)/COLPIV
                 ARRAY(I,J,K) = COLMLT
                 IF(I.EQ.NRWBLK)GO TO 350
                    DO 345 L = IPLUS1,NRWBLK
                       ARRAY(L,J,K) = ARRAY(L,J,K)
     *                                -COLMLT*ARRAY(L,IPLUSN,K)
345                 CONTINUE
350              CONTINUE
                 JRWBLK = J-NRWBLK
                 IF(K.EQ.NBLOKS)GO TO 370
                    DO 360 L = 1,NRWBLK
                       ARRAY(L,JRWBLK,KPLUS1) =
     *                                  ARRAY(L,JRWBLK,KPLUS1)
     *                         -COLMLT*ARRAY(L,IRWBLK,KPLUS1)
360                 CONTINUE
                    GO TO 380
370              CONTINUE
                 DO 375 L = 1,NRWBOT
                    BOTBLK(L,JRWBLK) = BOTBLK(L,JRWBLK)
     *                              -COLMLT*BOTBLK(L,IRWBLK)
375              CONTINUE
380           CONTINUE
390        CONTINUE
           INCR = INCR + NRWBLK
395     CONTINUE
C
C***************************************************************
C
C          ****  FINALLY, IN BOTBLK....
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C          ***  APPLY NRWBOT ROW ELIMINATIONS WITH ROW
C                  PIVOTING....
C
C               IF BOT HAS JUST ONE ROW GO TO 500
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
        IF(NRWBOT.EQ.1)GO TO 500
           DO 470 J = NRWTP1,NVRLP0
              JPLUS1 = J+1
              JMINN = J-NRWTOP
C
C             CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C               DETERMINE ROW PIVOT AND PIVOT INDEX
C
C             CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
              IPVT = JMINN
              ROWMAX = DABS(BOTBLK(JMINN,J))
              LOOP = JMINN+1
              DO 410 I = LOOP,NRWBOT
                 TEMPIV = DABS(BOTBLK(I,J))
                 IF(TEMPIV.LE.ROWMAX) GO TO 410
                 IPVT = I
                 ROWMAX = TEMPIV
410           CONTINUE
C
C             CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C               TEST FOR SINGULARITY:
C
C                       IF SINGULAR THEN TERMINATE AT 1000;
C                                   ELSE CONTINUE.
C
C             CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
              IF(PIVMAX+ROWMAX.EQ.PIVMAX)GO TO 1000
              PIVMAX = DMAX1(ROWMAX,PIVMAX)
C
C             CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C               IF NECESSARY INTERCHANGE ROWS
C
C             CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
              INCRJ = INCR+J
              PIVOT(INCRJ) = INCR+IPVT+NRWTOP
              IF(IPVT.EQ.JMINN)GO TO 430
                 DO 420 L = J,NOVRLP
                    SWAP = BOTBLK(IPVT,L)
                    BOTBLK(IPVT,L) = BOTBLK(JMINN,L)
                    BOTBLK(JMINN,L) = SWAP
420              CONTINUE
430           CONTINUE
C
C             CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C               COMPUTE MULTIPLIERS
C
C             CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
              ROWPIV = BOTBLK(JMINN,J)
              DO 440 I = LOOP,NRWBOT
                 BOTBLK(I,J) = BOTBLK(I,J)/ROWPIV
440           CONTINUE
C
C             CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C               PERFORM ROW ELIMINATION WITH COLUMN INDEXING
C
C             CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
              DO 460 L = JPLUS1,NOVRLP
                 ROWMLT = BOTBLK(JMINN,L)
                 DO 450 I = LOOP,NRWBOT
                    BOTBLK(I,L) = BOTBLK(I,L)-ROWMLT*BOTBLK(I,J)
450              CONTINUE
460           CONTINUE
470        CONTINUE
500     CONTINUE
C
C***************************************************************
C
C          DONE PROVIDED THE LAST ELEMENT IS NOT ZERO
C
C***************************************************************
C
        IF(PIVMAX+DABS(BOTBLK(NRWBOT,NOVRLP)).NE.PIVMAX) RETURN
C
C***************************************************************
C
C       ****  MATRIX IS SINGULAR - SET IFLAG = - 1.
C                                  TERMINATE AT 1000.
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
1000    CONTINUE
        IFLAG = -1
        RETURN
        END
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc



        SUBROUTINE CRSLVE(TOPBLK,NRWTOP,NOVRLP,ARRAY,NRWBLK,
     *             NCLBLK,NBLOKS,BOTBLK,NRWBOT,PIVOT,B,X)
C
C***************************************************************
C
C  C R S L V E  SOLVES THE LINEAR SYSTEM
C                       A*X = B
C  USING THE DECOMPOSITION ALREADY GENERATED IN  C R D C M P.
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C               *****  PARAMETERS  *****
C
C       *** ON ENTRY  ...
C
C               TOPBLK - DOUBLE PRECISION(NRWTOP,NOVRLP)
C                         OUTPUT FROM  C R D C M P
C
C               NOVRLP - INTEGER
C                         THE NUMBER OF COLUMNS IN WHICH SUCC-
C                         ESSIVE BLOCKS OVERLAP, WHERE
C                                NOVRLP = NRWTOP + NRWBOT
C
C               NRWTOP - INTEGER
C                         NUMBER OF ROWS IN THE BLOCK TOPBLK
C
C               ARRAY  - DOUBLE PRECISION(NRWBLK,NCLBLK,NBLOKS)
C                         OUTPUT FROM  C R D C M P
C
C               NRWBLK - INTEGER
C                         NUMBER OF ROWS IN K-TH BLOCK
C
C               NCLBLK - INTEGER
C                         NUMBER OF COLUMNS IN K-TH BLOCK
C
C               NBLOKS - INTEGER
C                         NUMBER OF NRWBLK BY NCLBLK BLOCKS IN
C                         THE MATRIX A
C
C               BOTBLK - DOUBLE PRECISION(NRWBOT,NOVRLP)
C                         OUTPUT FROM  C R D C M P
C
C               NRWBOT - INTEGER
C                         NUMBER OF ROWS IN THE BLOCK BOTBLK
C
C                PIVOT - INTEGER(N)
C                         THE PIVOT VECTOR FROM  C R D C M P
C
C                    B - DOUBLE PRECISION(N)
C                         THE RIGHT HAND SIDE VECTOR
C
C                    X - DOUBLE PRECISION(N)
C                         WORK SPACE
C
C       *** ON RETURN  ...
C
C          
C                X - DOUBLE PRECISION(N)
C                         THE SOLUTION VECTOR
C
C***************************************************************
C
        DOUBLE PRECISION TOPBLK,ARRAY,BOTBLK,X,B
        DOUBLE PRECISION DOTPRD,XJ,XINCRJ,BINCRJ,SWAP
        INTEGER PIVOT(*),NRWTOP,NRWBLK,NCLBLK,NRWBOT
        INTEGER NRWTP1,NRWBK1,NVRLP1,NOVRLP,NRWTP0,NRWBT1
        INTEGER NROWEL,NRWEL1,NVRLP0,NBLKS1,NBLOKS,NBKTOP
        INTEGER J,LOOP,I,INCR,K,INCRTP,INCRJ,INCRI,JPIVOT
        INTEGER JRWTOP,LL,NRWBTL,L,L1,IPLUSN,INCRN,IPVTN
        INTEGER NRWELL,IPVTI
        DIMENSION TOPBLK(NRWTOP,*),ARRAY(NRWBLK,NCLBLK,*),
     *          BOTBLK(NRWBOT,*),B(*),X(*)
C
C***************************************************************
C
C          ****  DEFINE THE CONSTANTS USED THROUGHOUT  ****
C
C***************************************************************
C
        NRWTP1 = NRWTOP+1
        NRWBK1 = NRWBLK+1
        NVRLP1 = NOVRLP+1
        NRWTP0 = NRWTOP-1
        NRWBT1 = NRWBOT+1
        NROWEL = NRWBLK-NRWTOP
        NRWEL1 = NROWEL+1
        NVRLP0 = NOVRLP-1
        NBLKS1 = NBLOKS+1
        NBKTOP = NRWBLK+NRWTOP
C
C***************************************************************
C
C               ****  FORWARD RECURSION  ****
C
C***************************************************************
C
C          ***  FIRST, IN TOPBLK....
C
C       CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C               FORWARD SOLUTION
C
C       CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
        DO 130 J = 1,NRWTOP
           X(J) = B(J)/TOPBLK(J,J)
           IF(J.EQ.NRWTOP)GO TO 120
              XJ = -X(J)
              LOOP = J+1
              DO 110 I = LOOP,NRWTOP
                 B(I) = B(I)+TOPBLK(I,J)*XJ
110           CONTINUE
120        CONTINUE
130     CONTINUE
C
C       ********************************************************
C
C          ***  IN EACH BLOCK ARRAY(,,K)....
C
C       CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
        INCR = 0
        DO 280 K = 1,NBLOKS
           INCRTP = INCR+NRWTOP
C
C          CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C               FORWARD MODIFICATION
C
C          CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
           DO 220 J = 1,NRWTOP
              INCRJ = INCR+J
              XINCRJ = -X(INCRJ)
              DO 210 I = 1,NRWBLK
                 INCRI = INCRTP+I
                 B(INCRI) = B(INCRI)+ARRAY(I,J,K)*XINCRJ
210           CONTINUE
220        CONTINUE
C
C          CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C               FORWARD ELIMINATION
C
C          CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
           DO 240 J = NRWTP1,NRWBLK
              INCRJ = INCR+J
              JPIVOT = PIVOT(INCRJ)
              IF(JPIVOT.EQ.INCRJ)GO TO 225
                 SWAP = B(INCRJ)
                 B(INCRJ) = B(JPIVOT)
                 B(JPIVOT) = SWAP
225           CONTINUE
              BINCRJ = -B(INCRJ)
              LOOP = J-NRWTP0
              DO 230 I = LOOP,NRWBLK
                 INCRI = INCRTP+I
                 B(INCRI) = B(INCRI)+ARRAY(I,J,K)*BINCRJ
230           CONTINUE
240        CONTINUE
C
C          CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C               FORWARD SOLUTION
C
C          CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
           DO 270 J = NRWBK1,NBKTOP
              INCRJ = INCR+J
              JRWTOP = J -NRWTOP
              X(INCRJ) = B(INCRJ)/ARRAY(JRWTOP,J,K)
              IF(J.EQ.NBKTOP)GO TO 260
                 XINCRJ = -X(INCRJ)
                 LOOP = J-NRWTP0
                 DO 250 I = LOOP,NRWBLK
                    INCRI = INCRTP+I
                    B(INCRI) = B(INCRI)+ARRAY(I,J,K)*XINCRJ
250              CONTINUE
260           CONTINUE
270        CONTINUE
           INCR = INCR+NRWBLK
280     CONTINUE
C
C       ********************************************************
C
C          ***  FINALLY, IN BOTBLK....
C
C       CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C               FORWARD MODIFICATION
C
C       CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
        INCRTP = INCR+NRWTOP
        DO 320 J = 1,NRWTOP
           INCRJ = INCR+J
           XINCRJ = -X(INCRJ)
           DO 310 I = 1,NRWBOT
              INCRI = INCRTP+I
              B(INCRI) = B(INCRI)+BOTBLK(I,J)*XINCRJ
310        CONTINUE
320     CONTINUE
C
C       CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C               FORWARD ELIMINATION
C
C       CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
        IF(NRWBOT.EQ.1)GO TO 350
           DO 340 J = NRWTP1,NVRLP0
              INCRJ = INCR+J
              JPIVOT = PIVOT(INCRJ)
              IF(JPIVOT.EQ.INCRJ)GO TO 325
                 SWAP = B(INCRJ)
                 B(INCRJ) = B(JPIVOT)
                 B(JPIVOT) = SWAP
325           CONTINUE
              BINCRJ = -B(INCRJ)
              LOOP = J-NRWTP0
              DO 330 I = LOOP,NRWBOT
                 INCRI = INCRTP+I
                 B(INCRI) = B(INCRI)+BOTBLK(I,J)*BINCRJ
330           CONTINUE
340        CONTINUE
350     CONTINUE
C
C***************************************************************
C
C               ****  BACKWARD RECURSION  ****
C
C***************************************************************
C
C          ***  FIRST IN BOTBLK....
C
C       CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C               BACKWARD SOLUTION
C
C       CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
        DO 430 LL = 1,NRWBOT
           J = NVRLP1-LL
           INCRJ = INCR+J
           NRWBTL = NRWBT1-LL
           X(INCRJ) = B(INCRJ)/BOTBLK(NRWBTL,J)
           IF(LL.EQ.NRWBOT)GO TO 420
              XINCRJ = -X(INCRJ)
              LOOP = NRWBOT-LL
              DO 410 I = 1,LOOP
                 INCRI = INCRTP+I
                 B(INCRI) = B(INCRI)+BOTBLK(I,J)*XINCRJ
410           CONTINUE
420        CONTINUE
430     CONTINUE
C
C       ********************************************************
C
C          ***  THEN IN EACH BLOCK ARRAY(,,K)....
C
C       CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
        DO 490 L = 1,NBLOKS
C
C          CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C               BACKWARD ELIMINATION
C
C          CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
           K = NBLKS1-L
           INCR = INCR-NRWBLK
           DO 450 L1 = NRWEL1,NRWBLK
              I = NRWBLK+NRWEL1-L1
              IPLUSN = I+NRWTOP
              LOOP = IPLUSN+1
              INCRN = INCR+IPLUSN
              DOTPRD = X(INCRN)
              DO 440 J = LOOP,NCLBLK
                 INCRJ = INCR+J
                 DOTPRD = DOTPRD-ARRAY(I,J,K)*X(INCRJ)
440           CONTINUE
              X(INCRN) = DOTPRD
              IPVTN = PIVOT(INCRN)
              IF(INCRN.EQ.IPVTN)GO TO 445
                 SWAP = X(INCRN)
                 X(INCRN) = X(IPVTN)
                 X(IPVTN) = SWAP
445           CONTINUE
450        CONTINUE
C
C          CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C               BACKWARD MODIFICATION
C
C          CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
           INCRTP = INCR+NRWTOP
           DO 460 J = NRWBK1,NCLBLK
              INCRJ = INCR+J
              XINCRJ = -X(INCRJ)
              DO 455 I = 1,NROWEL
                 INCRI = INCRTP+I
                 B(INCRI) = B(INCRI)+ARRAY(I,J,K)*XINCRJ
455           CONTINUE
460        CONTINUE
C
C          CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C               BACKWARD SOLUTION
C
C          CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
           DO 480 LL = 1,NROWEL
              J = NRWBK1-LL
              INCRJ = INCR+J
              NRWELL = NRWEL1-LL
              X(INCRJ) = B(INCRJ)/ARRAY(NRWELL,J,K)
              IF(LL.EQ.NROWEL)GO TO 470
                 XINCRJ = -X(INCRJ)
                 LOOP = NROWEL-LL
                 DO 465 I = 1,LOOP
                    INCRI = INCRTP+I
                    B(INCRI) = B(INCRI)+ARRAY(I,J,K)*XINCRJ
465              CONTINUE
470           CONTINUE
480        CONTINUE
490     CONTINUE
C
C       ********************************************************
C
C          ***  IN TOPBLK FINISH WITH....
C
C       CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C               BACKWARD ELIMINATION
C
C       CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
        DO 520 L = 1,NRWTOP
           I = NRWTP1-L
           LOOP = I+1
           DOTPRD = X(I)
           DO 510 J = LOOP,NOVRLP
              DOTPRD = DOTPRD-TOPBLK(I,J)*X(J)
510        CONTINUE
           X(I) = DOTPRD
           IPVTI = PIVOT(I)
           IF(I.EQ.IPVTI)GO TO 515
                 SWAP = X(I)
                 X(I) = X(IPVTI)
                 X(IPVTI) = SWAP
515        CONTINUE
520     CONTINUE
        RETURN
        END
