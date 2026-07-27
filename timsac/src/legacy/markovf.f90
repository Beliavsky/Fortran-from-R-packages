! Part of the experimental modern Fortran port of timsac 1.3.8-6.
! This file was created or modified for the Fortran project on 2026-07-23.
! Original TIMSAC credits are retained; see NOTICE and ORIGIN.md.
! SPDX-License-Identifier: GPL-2.0-or-later

subroutine markovf(n,lagh3,id,cyy0,k,nh,jaw,aw1,b1,icont,idd,&
&ir,ij,ik,ipq,g,a1,a,b,vd,iqm,bm,au,zz,c0,aicd,mj3,mj4,mj6,mj7)
  use timsac_kinds, only: dp
  implicit none
!
!c       PROGRAM MARKOV
!     PROGRAM 74.2.2. MAXIMUM LIKELIHOOD COMPUTATION OF MARKOVIAN MODEL
!-----------------------------------------------------------------------
!     ** DESIGNED BY H. AKAIKE, THE INSTITUTE OF STATISTICAL MATHEMATICS
!     ** PROGRAMMED BY E. ARAHATA, THE INSTITUTE OF STATISTICAL MATHEMAT
!         TOKYO
!     ** DATE OF THE LATEST REVISION: DECEMBER 6, 1978
!     ** THIS PROGRAM WAS ORIGINALLY PUBLISHED IN
!        "TIMSAC-74 A TIME SERIES ANALYSIS AND CONTROL PROGRAM PACKAGE(1
!        BY H. AKAIKE, E. ARAHATA AND T. OZAKI, COMPUTER SCIENCE MONOGRA
!        NO.5, MARCH 1975, THE INSTITUTE OF STATISTICAL MATHEMATICS
!     ** FOR THE BASIC THEORY SEE "CANONICAL CORRELATION ANALYSIS OF TIM
!        AND THE USE OF AN INFORMATION CRITERION" BY H. AKAIKE, IN
!        "SYSTEM IDENTIFICATION: ADVANCES AND CASE STUDIES" R. K. MEHRA
!        D. G. LAINIOTIS EDS. ACADEMIC PRESS, NEW YORK, 1976
!-----------------------------------------------------------------------
!     IN THIS PROGRAM THE MATRICES A,B,C, STAND FOR TRANSITION MATRIX (F
!     INPUT MATRIX (G), OUTPUT MATRIX (H), RESPECTIVELY.
!
!     THE INPUTS REQUIRED ARE AS FOLLOWS:
!     (N,LAGH0,ID0):
!         N, LENGTH OF ORIGINAL DATA
!         LAGH0, MAXIMUM LAG OF COVARIANCE
!         ID0, DIMENSION OF Y(I)
!     (CYY(I)(I=0,LAGH0): COVARIANCE MATRIX SEQUENCE OF Y(I). CYY(I) ARE
!                           THE OUTPUTS OF THE PROGRAM MULCOR OF TIMSAC AND
!                           ARE USED AS THE INPUT TO THE PROGRAM CANOCA OF
!                           TIMSAC-74. FOR THE CONSTRUCTION OF CYY(I),
!                           SEE PROGRAM CANOCA.
!     THE OUTPUTS OF PROGRAM CANACA:
!          (ID,K):
!              ID, DIMENSION OF THE TIME SERIES Y(I) (NOT GREATER THAN 5)
!              K, DIMENSION OF THE STATE VECTOR (NOT GREATER THAN 10)
!          (NH(I))(I=1,K): STRUCTURAL CHARACTERISTIC VECTOR
!          (AW(I))(I=1,IAW): INITIAL ESTIMATE OF THE VECTOR OF FREE
!                               PARAMETERS IN F (=A)
!          B(I,J)(I=1,ID+1;J=1,ID): INITIAL ESTIMATES OF THE FREE
!                                       PARAMETERS IN G (=B)
!     ICONT: OUTPUT CONTROL
!          = 0, FOR AR-MA COEFFICIENTS
!          = 1, FOR SIMCON INPUT
!          = 2, FOR BOTH
!
!     THE OUTPUTS OF THIS PROGRAM ARE;
!     THE MAXIMUM LIKELIHOOD ESTIMATES
!          (ID,K):
!          (NH(I))(I=1,K):
!          (AW(I))(I=1,IAW):
!          (B(I,J))(I=1,ID+1;J=1,ID):
!     AND THE OUTPUTS TO BE USED BY THE PROGRAMS PRDCTR AND SIMCON
!          (ID,Q,Q-1):
!          (B(I,J,L))(I,J=1,ID,L=1,Q): AR-COEFFICIENT MATRICES
!          AND, WHEN ICONT=1 OR 2,
!          (W(I,J,L))(I,J=1,ID,L=1,Q-1): IMPULSE RESPONSE MATRICES
!          AND, WHEN ICONT=0 OR 1,
!          (A(I,J,L))(I,J=1,ID;L=1,Q-1): MA-COEFFICIENT MATRICES
!     AND
!          (C0(I,J))(I,J=1,ID): INNOVATION COVARIANCE.
!
!     COMMENTS:
!     ***** IT IS ASSUMED THAT ID0=ID. OTHERWISE USE THE LOADING
!     PROGRAM OF CYY(I) USED BY THE PROGRAM CANOCA.
!     ***** THE INPUT TO THIS PROGRAM IS NOT LIMITED TO THE OUTPUT OF
!     PROGRAM CANOCA. THE FINAL DECISION ON THE STRUCTURAL
!     CHARACTERISTIC VECTOR NH(I) (I=1,K) SHOULD BE MADE BY COMPARING
!     THE AIC VALUES OF VARIOUS POSSIBLE CHOICES OF NH WITH
!     THE MINIMUM AIC OF THE AR-MODEL.
!     ***** TO KEEP THE ACCURACY OF COMPUTATION THE COMPONENTS OF Y(I)
!     SHOULD BE SCALED SO THAT THE VARIANCES ARE NEARLY EQUAL TO EACH
!     OTHER.
!     *****
!
!xx      IMPLICIT REAL*8(A-H,O-Z)
!c      COMMON /COM10/CYY1
!c      COMMON /COM50/VD
!c      COMMON /COM60/IDD,IR
!c      COMMON /COM61/IJ,IK
!c      COMMON /COM70/NH
!c      COMMON /COM80/A
!c      COMMON /COM81/B
!c      COMMON /COM82/AW
! common /COM99/
integer iswr
common /com99/iswr
!c      COMMON /COM87/ICONT
! common /COM101/
real(dp) aico
common /com101/aico
! common /COM102/
integer kswr
common /com102/kswr
!c      COMMON /COM199/AICD
!c      DIMENSION CYY1(160,5,5)
!c      DIMENSION A(10,10),B(10,5),OB(10,5)
!c      DIMENSION VD(50,50)
!c      DIMENSION X(100),C0(5,5),G(100),R(100)
!c      DIMENSION AW(50),OAW(50)
!c      DIMENSION NH(10),IDD(10),IR(10),IJ(5),IK(5)
!xx      DIMENSION CYY0(LAGH3,ID,ID),CYY1(MJ3,ID,ID)
!xx      DIMENSION A(K,K),B(K,ID),OB(K,ID)
!xx      DIMENSION VD(MJ4,MJ4)
!xx      DIMENSION X(MJ4),C0(ID,ID),G(MJ4),R(MJ4)
!xx      DIMENSION AW(JAW),OAW(JAW)
!xx      DIMENSION A1(K,K),B1(K,ID),AW1(JAW)
!xx      DIMENSION NH(K),IDD(K),IR(K),IJ(ID),IK(ID)
!xx      DIMENSION BM(ID,ID,MJ6),AU(ID,ID,MJ7),ZZ(ID,ID,MJ7)
!
!xx      DIMENSION O(50,ID,K),Q(50,K,ID),X1(50,K,K),X2(50,ID,ID)
!xx      DIMENSION CXY1(100,ID,ID),CXX1(51,ID,ID),CXV1(51,ID,K)
integer n, lagh3, id, k, nh(k), jaw, icont, idd(k), ir(k), ij(id),&
&ik(id), ipq, iqm, mj3, mj4, mj6, mj7
real(dp) cyy0(lagh3,id,id), aw1(jaw), b1(k,id), g(mj4),&
&a1(k,k), a(k,k), b(k,id), vd(mj4,mj4),&
&bm(id,id,mj6), au(id,id,mj7), zz(id,id,mj7),&
&c0(id,id), aicd
! local
integer i, ii, iaw, idp1, ig, isfin, iswro, j, jj, jl, js, kk,&
&l, m
real(dp) cyy1(mj3,id,id), ob(k,id), x(mj4), r(mj4),&
&aw(jaw), oaw(jaw), o(50,id,k), q(50,k,id),&
&x1(50,k,k), x2(50,id,id), cxy1(100,id,id),&
&cxx1(51,id,id), cxv1(51,id,k), cst0, cst1, cst05,&
&aipq, sum
!
!x      INTEGER*1  TMP(1)
!x      CHARACTER  CNAME*80
!
!     INPUT / OUTPUT DATA FILE OPEN
!c      CHARACTER(100) DFNAM
!c      CALL SETWND
!c      DFNAM='markov.out'
!c      CALL FLOPN3(DFNAM,NFL)
!c      IF (NFL.EQ.0) GO TO 999
!
!x      IER=0
!x      LU=3
!x      DO 10 I = 1,80
!x   10 CNAME(I:I) = ' '
!x      I = 1
!x      IFG = 1
!x      DO WHILE( (IFG.EQ.1) .AND. (I.LE.80) )
!x         IF ( TMP(I).NE.ICHAR(' ') ) THEN
!x            CNAME(I:I) = CHAR(TMP(I))
!x            I = I+1
!x         ELSE
!x            IFG = 0
!x         END IF
!x      END DO
!x      IF ( I.GT.1 ) THEN
!x         IFG = 1
!x         OPEN (LU,FILE=CNAME,IOSTAT=IVAR)
!x         IF (IVAR .NE. 0) THEN
!xcx            WRITE(*,*) ' ***  markov temp FILE OPEN ERROR :',CNAME,IVAR
!x            IER=IVAR
!x            IFG=0
!x         END IF
!x      END IF
!
!     ABSOLUTE DIMENSIONS USED FOR SUBROUTINE CALL
!c      MJ1=10
!c      MJ2=5
cst0=0.0d-00
cst1=1.0d-00
cst05=0.00005d-00
!     INITIAL INPUT
!c	READ(5,1) N,LAGH0,ID
!c	LAGH3=LAGH0+1
!     COVARIANCE INPUT
!xx      DO 1010 II=1,ID
do 1011 ii=1,id
!c      READ(5,1) IR0,IC0
!c      READ(5,2) (CYY1(I,IR0,IC0),I=1,LAGH3)
!c      IF(II.EQ.1) GO TO 1010
!c      IIM1=II-1
!c      DO 1011 JJ=1,IIM1
!c      READ(5,1) IR1,IC1
!c      READ(5,2) (CYY1(I,IR1,IC1),I=1,LAGH3)
!c      READ(5,1) IR2,IC2
!c      READ(5,2) (CYY1(I,IR2,IC2),I=1,LAGH3)
!c 1011 CONTINUE
do 1010 jj=1,id
do 1008 kk=1,mj3
cyy1(kk,jj,ii)=0
1008 continue
do 1009 kk=1,lagh3
cyy1(kk,jj,ii)=cyy0(kk,jj,ii)
1009 continue
1010 continue
1011 continue
!
!c      READ(5,1) ID,K
!c      READ(5,1) (NH(I),I=1,K)
!c      CALL SUBIDR(K,ID,IAW)
call subidr(nh,idd,ir,ij,ik,k,id,iaw)
!     INITIAL PRINT OUT
!c	WRITE(6,900)
!c	WRITE(6,901)
!c	WRITE(6,903) ID,N,LAGH0
!c	WRITE(6,904) K,IAW
!
!     A MATRIX ARRANGEMENT
!c	READ(5,2) (AW(I),I=1,IAW)
do 4009 i=1,jaw
aw(i)=aw1(i)
4009 continue
js=0
!xx      DO 4010 I=1,K
!xx      DO 4010 J=1,K
!c 4010 A(I,J)=CST0
!xx      A(I,J)=CST0
!xx      A1(I,J)=CST0
!xx 4010 CONTINUE
a(1:k,1:k)=cst0
a1(1:k,1:k)=cst0
do 4011 i=1,k
if(idd(i).ne.0) go to 4012
jj=ir(i)
a(i,jj)=cst1
a1(i,jj)=cst1
go to 4011
4012 jl=ir(i)
do 4013 j=1,jl
js=js+1
!c 4013 A(I,J)=AW(JS)
a(i,j)=aw(js)
a1(i,j)=aw(js)
4013 continue
4011 continue
!c	DO 4020 I=1,ID
!c	DO 4019 J=1,ID
!c 4019 B(I,J)=CST0
!c 4020 B(I,I)=CST1
idp1=id+1
!c      DO 4021 I=IDP1,K
!c 4021 READ(5,2) (B(I,J),J=1,ID)
do 4020 i=1,k
do 4019 j=1,id
b(i,j)=b1(i,j)
4019 continue
4020 continue
!c      READ(5,1) ICONT
!     NH(I) PRINT OUT
!c      WRITE(6,4300) (I,I=1,K)
!c      WRITE(6,4310) (NH(I),I=1,K)
!     IDD(I),IR(I) PRINT OUT
!c      WRITE(6,3610) (IDD(I),I=1,K)
!c      WRITE(6,3620) (IR(I),I=1,K)
!     IJ(I),IK(I) PRINT OUT
!c      WRITE(6,3650) (IJ(I),I=1,ID)
!c      WRITE(6,3660) (IK(I),I=1,ID)
!c      WRITE(6,6000)
!c      WRITE(6,6001)
!c      WRITE(6,6002)
!c      WRITE(6,6004)
!c      WRITE(6,6006)
!c      WRITE(6,6008)
!c      WRITE(6,6100)
!     A,B PRINT OUT
!c      WRITE(6,52)
!c      DO 180 I=1,K
!c  180 WRITE(6,2200) I,(A(I,J),J=1,K)
!c      WRITE(6,53)
!c      DO 181 I=1,K
!c  181 WRITE(6,2200) I,(B(I,J),J=1,ID)
!c      WRITE(6,340)
!c      DO 360 I=1,IAW
!c  360 WRITE(6,351) I,AW(I)
iswr=0
iswro=0
!
!     IPQ=TOTAL NUMBER OF FREE PARAMETERS
400 ipq=iaw+(k-id)*id
aipq=ipq
idp1=id+1
!c	WRITE(6,411) ISWRO
!x      IF (IFG.NE.0) WRITE(LU,411) ISWRO
!xx  411 FORMAT(/1H ,'ISWRO=',I5)
do 500 i=1,iaw
!xx  500 X(I)=AW(I)
x(i)=aw(i)
500 continue
jj=iaw
!xx      DO 510 I=IDP1,K
do 511 i=idp1,k
do 510 j=1,id
jj=jj+1
!xx  510 X(JJ)=B(I,J)
x(jj)=b(i,j)
510 continue
511 continue
!     INNOVATION COVARIANCE MATRIX (CXX0(0)) AND GRADIENT G(X)
!     OF LOG(DET.CXX0(0)) COMPUTATION
kswr=0
ig=0
!c      CALL C0GR(X,C0,G,M,L,K,ID,IAW,IPQ,MJ1,MJ2,IG)
!xx      CALL C0GR(CYY1,MJ3,B,AW,IAW,NH,IDD,IR,IJ,IK,X,C0,G,M,L,K,ID,IPQ,
call c0gr(cyy1,mj3,b,aw,iaw,nh,idd,ir,ij,ik,x,c0,g,m,l,k,id,&
!x     *                MJ4,IG,O,Q,X1,X2,CXY1,CXX1,CXV1,IFG,LU)
&mj4,ig,o,q,x1,x2,cxy1,cxx1,cxv1)
do 1500 i=1,iaw
!xx 1500 X(I)=AW(I)
x(i)=aw(i)
1500 continue
jj=iaw
!xx      DO 1510 I=IDP1,K
do 1511 i=idp1,k
do 1510 j=1,id
jj=jj+1
!xx 1510 X(JJ)=B(I,J)
x(jj)=b(i,j)
1510 continue
1511 continue
!c	WRITE(6,340)
!c	DO 1360 I=1,IAW
!c 1360 WRITE(6,351) I,AW(I)
!c	WRITE(6,53)
!c	DO 2181 I=1,K
!c 2181 WRITE(6,2200) I,(B(I,J),J=1,ID)
!x      IF (IFG.NE.0) THEN
!x	 WRITE(LU,340)
!x	 DO 1360 I=1,IAW
!x 1360	 WRITE(LU,351) I,AW(I)
!x	 WRITE(LU,53)
!x	 DO 2181 I=1,K
!x 2181	 WRITE(LU,2200) I,(B(I,J),J=1,ID)
!x      END IF
!
!     HESSIAN COMPUTATION
!c	CALL SUBHES(M,L,K,ID,IAW,MJ1,MJ2)
call subhes(cyy1,ij,ik,vd,m,l,k,id,iaw,mj3,mj4,o,q,x1,x2,&
&cxy1,cxx1,cxv1)
if(iswro.gt.0) go to 2651
!c	WRITE(6,2000)
!c	DO 2100 I=1,IPQ
!c 2100 WRITE(6,2200) I,(VD(I,J),J=1,IPQ)
!x      IF (IFG.NE.0) THEN
!x	   WRITE(LU,2000)
!x	 DO 2100 I=1,IPQ
!x 2100	 WRITE(LU,2200) I,(VD(I,J),J=1,IPQ)
!x	 WRITE(LU,'(/)')
!x      END IF
!     (INVERSE OF HESSIAN) COMPUTATION
!c 2651 CALL MATINV(IPQ)
!xx 2651 CALL MATINV(IPQ,VD,MJ4,IFG,LU)
2651 call matinv(ipq,vd,mj4)
if(iswro.gt.0) go to 2650
!c	WRITE(6,2500)
!c	DO 2600 I=1,IPQ
!c 2600 WRITE(6,2200) I,(VD(I,J),J=1,IPQ)
!x      IF (IFG.NE.0) THEN
!x	 WRITE(LU,2500)
!x	 DO 2600 I=1,IPQ
!x 2600	 WRITE(LU,2200) I,(VD(I,J),J=1,IPQ)
!x      END IF
!
!     CORRECTION R(X)=V*G(X) COMPUTATION
2650 do 100 i=1,ipq
sum=cst0
do 110 j=1,ipq
!xx  110 SUM=SUM+VD(I,J)*G(J)
sum=sum+vd(i,j)*g(j)
110 continue
!xx  100 R(I)=SUM
r(i)=sum
100 continue
!     ORIGINAL PARAMETER VALUES
do 2655 i=1,iaw
!xx 2655 OAW(I)=AW(I)
oaw(i)=aw(i)
2655 continue
!xx      DO 2656 I=IDP1,K
do 2657 i=idp1,k
do 2656 j=1,id
!xx 2656 OB(I,J)=B(I,J)
ob(i,j)=b(i,j)
2656 continue
2657 continue
!     DAVIDON'S PROCEDURE
!c	CALL SUBDAV(X,C0,G,R,N,K,ID,IAW,IPQ,MJ1,MJ2,ISWRO)
call subdav(cyy1,mj3,b,aw,iaw,nh,idd,ir,ij,ik,vd,x,c0,g,r,n,k,&
!x     * ID,IPQ,AICD,MJ1,MJ4,ISWRO,O,Q,X1,X2,CXY1,CXX1,CXV1,IFG,LU)
!xx     * ID,IPQ,AICD,MJ1,MJ4,ISWRO,O,Q,X1,X2,CXY1,CXX1,CXV1)
&id,ipq,aicd,mj4,iswro,o,q,x1,x2,cxy1,cxx1,cxv1)
isfin=0
if(iswro.ge.20) go to 599
if (iswro.eq.0) go to 599
do 150 i=1,iaw
if(dabs(oaw(i)-x(i)).ge.cst05) go to 153
150 continue
jj=iaw
!xx      DO 151 I=IDP1,K
do 152 i=idp1,k
do 151 j=1,id
jj=jj+1
if(dabs(ob(i,j)-x(jj)).ge.cst05) go to 153
151 continue
152 continue
go to 599
153 isfin=1
!
!
599 do 600 i=1,iaw
!xx  600 AW(I)=X(I)
aw(i)=x(i)
600 continue
idp1=id+1
jj=iaw
!xx      DO 610 I=IDP1,K
do 611 i=idp1,k
do 610 j=1,id
jj=jj+1
!xx  610 B(I,J)=X(JJ)
b(i,j)=x(jj)
610 continue
611 continue
if(isfin.eq.0) go to 1450
go to 400
!
!     FINAL PRINT OUT
!c 1450 WRITE(6,1460)
!c 1460 FORMAT(//1H ,'FINAL RESULTS')
1450 continue
!
!     GRADIENT PRINT OUT
!c	WRITE(6,1400)
!c	CALL SUBVCP(G,IPQ)
!
!     A ARRANGEMENT
!c	CALL SUBA(K,ID,MJ1,MJ2,IAW)
call suba(a,aw,ij,ik,k,id,iaw)
!     A, B PRINT OUT
!c	WRITE(6,52)
!c	DO 1180 I=1,K
!c 1180 WRITE(6,2200) I,(A(I,J),J=1,K)
!c	WRITE(6,53)
!c	DO 1181 I=1,K
!c 1181 WRITE(6,2200) I,(B(I,J),J=1,ID)
!     VD(I,J) PRINT OUT
!c	WRITE(6,2550)
!c	DO 2800 I=1,IPQ
!c 2800 WRITE(6,2200) I,(VD(I,J),J=1,IPQ)
!     A, B PUNCH OUT
!c	WRITE(7,1) ID,K
!c	WRITE(7,1) (NH(I),I=1,K)
!c	WRITE(7,2) (AW(I),I=1,IAW)
!c	DO 410 I=IDP1,K
!c  410 WRITE(7,2) (B(I,J),J=1,ID)
!c	CALL ARMACO(ID,K,MJ1,MJ2)
!xx      CALL ARMACO(A,B,AW,IAW,NH,IDD,IR,IJ,IK,ID,ICONT,K,BM,AU,ZZ,
!xx     *            IQM,MJ6,MJ7)
call armaco(b,aw,iaw,nh,idd,ir,ik,id,icont,k,bm,au,zz,&
&iqm,mj6,mj7)
!     INNOVATION VARIANCE PUNCH OUT
!c	DO 460 I=1,ID
!c  460 WRITE(7,2) (C0(I,J),J=1,ID)
!c	WRITE(6,6610)
!c	DO 6620 I=1,ID
!c 6620 WRITE(6,2200) I,(C0(I,J),J=1,ID)
!c	WRITE(6,3609) AICD
!
!c	CALL FLCLS3(NFL)
!xx  999 CONTINUE
!x      IF (IFG.NE.0) CLOSE(LU)
return
!xx    1 FORMAT(16I5)
!xx    2 FORMAT(4D20.10)
!xx   52 FORMAT(/1H ,'A(I,J)')
!xx   53 FORMAT(/1H ,'B(I,J)')
!xx  340 FORMAT(/1H ,'AW(I)')
!xx  351 FORMAT(1H ,I5,D17.5)
!xx  900 FORMAT(/1H ,'PROGRAM 74.2.2. MARKOV')
!xx  901 FORMAT(/1H ,'MAXIMUM LIKELIHOOD COMPUTATION OF ',
!xx     A'MULTIVARIATE MARKOVIAN MODEL.')
!xx  903 FORMAT(/1H ,'DIMENSION OF Y(I) (ID)=',I5,5X,'DATA LENGTH (N)=',
!xx     AI5,5X,'MAXIMUM LAG OF COVARIANCE (LAGH0)=',I5)
!xx  904 FORMAT(1H ,'DIMENSION OF THE STATE VECTOR (K)=',I5,5X,
!xx     A'TOTAL NUMBER OF FREE PARAMETERS (IAW)=',I5)
!xx 1400 FORMAT(/1H ,'GRADIENT')
!xx 4300 FORMAT(/1H ,5X,'I',16I5)
!xx 4310 FORMAT(1H ,1X,'NH(I)',16I5)
!xx 3610 FORMAT(1H ,'IDD(I)',16I5)
!xx 3620 FORMAT(1H ,1X,'IR(I)',16I5)
!xx 3650 FORMAT(1H ,1X,'IJ(I)',16I5)
!xx 3660 FORMAT(1H ,1X,'IK(I)',16I5)
!xx 2000 FORMAT(///1H ,'HESSIAN')
!c 2200 FORMAT(/1H ,I5,10D12.5,/(1H ,5X,10D12.5))
!xx 2200 FORMAT(/I5,5(1X,D12.5)/5X,5(1X,D12.5))
!xx 2500 FORMAT(///1H ,'INVERSE OF HESSIAN')
!xx 2550 FORMAT(///1H ,'DAVIDON VARIANCE')
!xx 6000 FORMAT(/1H ,'COMMENTS:')
!xx 6001 FORMAT(1H ,'NH(I) IS THE STRUCTURAL CHARACTERISTIC VECTOR.')
!xx 6002 FORMAT(1H ,'ID(I)=1 MEANS THAT THE I-TH ROW OF A CONTAINS FREE',
!xx     A' PARAMETERS.')
!xx 6004 FORMAT(1H ,'IR(I) DENOTES THE POSITION OF THE LAST NON-ZERO',
!xx     A' ELEMENT WITHIN THE I-TH ROW OF A.')
!xx 6006 FORMAT(1H ,'IJ(I) DENOTES THE POSITION OF THE I-TH NON-TRIVIAL',
!xx     A' ROW WITHIN A.')
!xx 6008 FORMAT(1H ,'IK(I) DENOTES THE NUMBER OF FREE PARAMETERS WITHIN',
!xx     A' THE I-TH NON-TRIVIAL ROW OF A.')
!xx 6100 FORMAT(/1H ,'INITIAL ESTIMATES')
!xx 6610 FORMAT(/1H ,'INNOVATION VARIANCE')
!xx 3609 FORMAT(/1H ,'AIC=N*LOG(DET(CXX0))+2.0*IPQ=',D12.5)
end
!
!c	SUBROUTINE C0GR(X,C0,GR,M,L,K,ID,IAW,IPQ,MJ1,MJ2,IG)
subroutine c0gr(cyy1,mj3,b,aw,iaw,nh,idd,ir,ij,ik,x,c0,gr,m,l,k,&
!x     *                   ID,IPQ,MJ4,IG,O,Q,X1,X2,CXY1,CXX1,CXV1,IFG,LU)
!xx     *                   ID,IPQ,MJ4,IG,O,Q,X1,X2,CXY1,CXX1,CXV1)
&id,mj4,ig,o,q,x1,x2,cxy1,cxx1,cxv1)
  use timsac_kinds, only: dp
  implicit none
!     THIS SUBROUTINE COMPUTES THE INNOVATION MATRIX CXX(0) AND
!     THE GRADIENT OF LOG(DET CXX(0))
!xx      IMPLICIT REAL*8(A-H,O-Z)
!c	COMMON /COM81/B
!c	COMMON /COM82/AW
! common /COM102/
integer kswr
common /com102/kswr
!c	DIMENSION CYY1(160,5,5),B(10,5),AW(50)
!c	DIMENSION IDD(10),IR(10),IJ(5),IK(5)
!c	DIMENSION X(100),C0(10,10),GR(100)
!xx      DIMENSION CYY1(MJ3,ID,ID),B(K,ID),AW(IAW)
!xx      DIMENSION NH(K),IDD(K),IR(K),IJ(ID),IK(ID)
!xx      DIMENSION X(MJ4),C0(ID,ID),GR(MJ4)
!xx      DIMENSION O(50,ID,K),Q(50,K,ID),X1(50,K,K),X2(50,ID,ID)
!xx      DIMENSION CXY1(100,ID,ID),CXX1(51,ID,ID),CXV1(51,ID,K)
integer mj3, iaw, m, l, k, id, mj4, ig, nh(k), idd(k),&
&ir(k), ij(id), ik(id)
real(dp) cyy1(mj3,id,id), b(k,id), aw(iaw), x(mj4),&
&c0(id,id), gr(mj4), o(50,id,k), q(50,k,id),&
&x1(50,k,k), x2(50,id,id), cxy1(100,id,id),&
&cxx1(51,id,id), cxv1(51,id,k)
! local
integer i, idp1, ig1, ig2, j, jj, kk
!
do 600 i=1,iaw
!xx  600 AW(I)=X(I)
aw(i)=x(i)
600 continue
idp1=id+1
jj=iaw
!xx      DO 610 I=IDP1,K
do 611 i=idp1,k
do 610 j=1,id
jj=jj+1
!xx  610 B(I,J)=X(JJ)
b(i,j)=x(jj)
610 continue
611 continue
!     IMPULSE RESPONSE COMPUTATION
ig1=ig
ig2=ig
!c	CALL NSUBX1(K,ID,M,MJ1,MJ2,IG1)
!c	CALL NSUBX2(K,ID,KK,MJ1,MJ2,IG2)
!x      CALL NSUBX1(B,AW,IAW,NH,IDD,IR,IJ,K,ID,M,IG1,X1,IFG,LU)
!x      CALL NSUBX2(B,AW,IAW,IDD,IR,K,ID,KK,IG2,Q,X2,IFG,LU)
!xx      CALL NSUBX1(B,AW,IAW,NH,IDD,IR,IJ,K,ID,M,IG1,X1)
call nsubx1(b,aw,iaw,nh,idd,ir,k,id,m,ig1,x1)
call nsubx2(b,aw,iaw,idd,ir,k,id,kk,ig2,q,x2)
ig=ig1+ig2
if(kswr.eq.0) go to 620
if(ig.ne.0) go to 1000
620 l=kk
!     CXY1(LS) LS=0,M+L+1 COMPUTATION
!c	CALL SBCXY1(M,L,ID,MJ2)
call sbcxy1(cyy1,mj3,m,l,id,x2,cxy1)
!     CXX1(NS) NS=0,M+1.
!c	CALL SUBCXX(C0,M,KK,ID,MJ2)
call subcxx(c0,m,kk,id,x2,cxy1,cxx1)
!     CXV1(MS) MS=0,M+1 COMPUTATION
!c	CALL GCXV1(M,L,K,ID,MJ1,MJ2)
call gcxv1(m,l,k,id,q,cxy1,cxv1)
!     O(NS) COMPUTATION
!c	CALL NSUBO(C0,K,ID,M,MJ1,MJ2)
call nsubo(c0,k,id,m,o,x1)
!     GRADIENT GA AND GB COMPUTATION
!c	CALL GRAD(GR,M,K,ID,IAW,IPQ,MJ1,MJ2)
!x      CALL GRAD(IJ,IK,GR,M,K,ID,IAW,IPQ,MJ4,O,CXX1,CXV1,IFG,LU)
!xx      CALL GRAD(IJ,IK,GR,M,K,ID,IAW,IPQ,MJ4,O,CXX1,CXV1)
call grad(ij,ik,gr,m,k,id,iaw,mj4,o,cxx1,cxv1)
1000 kswr=1
return
end
!
!c	SUBROUTINE NSUBX1(K,ID,M,MJ1,MJ2,IG1)
!x      SUBROUTINE NSUBX1(B,AW,IAW,NH,IDD,IR,IJ,K,ID,M,IG1,X1,IFG,LU)
!xx      SUBROUTINE NSUBX1(B,AW,IAW,NH,IDD,IR,IJ,K,ID,M,IG1,X1)
subroutine nsubx1(b,aw,iaw,nh,idd,ir,k,id,m,ig1,x1)
  use timsac_kinds, only: dp
  implicit none
!     THIS SUBROUTINE COMPUTES THE IMPULSE RESPONSE OF
!     (INVERSE OF (C*(INVERSE OF PHAI)*B))*C*
!     (INVERSE OF PHAI). PHAI = INVERSE OF (I-A*EXP(-I*2*PI*F)).
!     C*(INVERSE OF PHAI): FROM V TO Y
!     INVERSE OF (C*(INVERSE OF PHAI)*B): FROM Y TO X
!     K: DIMENSION OF THE STATE VECTOR
!     ID: DIMENSION OF Y(I)
!     M: UPPER LIMIT OF THE NUMBER OF IMPULSE RESPONSES.
!
!xx      IMPLICIT REAL*8(A-H,O-Z)
!c      COMMON /COM15/X1
!c	COMMON /COM60/IDD,IR
!c	COMMON /COM81/B
!c	COMMON /COM82/AW
! common /COM99/
integer iswr
common /com99/iswr
!c	DIMENSION IDD(10),IR(10)
!c	DIMENSION AW(50),B(10,5)
!xx      DIMENSION NH(K),IDD(K),IR(K),IJ(ID)
!xx      DIMENSION AW(IAW),B(K,ID)
!c      DIMENSION X1(50,10,10)
!c	DIMENSION V(10,10),AV(10,10)
!c	DIMENSION X(5,10),Y(5,10),U(10,10)
!c	DIMENSION AU(10,10)
!xx      DIMENSION X1(50,K,K)
!xx      DIMENSION V(K,K),AV(K,K)
!xx      DIMENSION X(ID,K),Y(ID,K),U(K,K)
!xx      DIMENSION AU(K,K)
integer iaw, k, id, m, ig1, nh(k), idd(k), ir(k)
real(dp) b(k,id), aw(iaw), x1(50,k,k)
! local
integer i, ih, itinp, j, jj, lh, ns, nsm1
real(dp) v(k,k), av(k,k), x(id,k), y(id,k), u(k,k),&
&au(k,k), cst0, cst1, constm, constd, sum, xmax
!
!
!
cst0=0.0d-00
cst1=1.0d-00
itinp=0
!     V(0)=I
!xx  199 DO 200 I=1,K
!xx      DO 210 J=1,K
!xx  210 V(I,J)=CST0
!xx  200 V(I,I)=CST1
199 continue
v(1:k,1:k)=cst0
do 200 i=1,k
v(i,i)=cst1
200 continue
!
!     C=(I,0)
!     I=(ID,ID), 0:(ID,(K-ID))
!     Y(0)=C
y(1:id,1:k)=cst0
do 300 i=1,id
!xx      DO 310 J=1,K
!xx  310 Y(I,J)=CST0
!xx  300 Y(I,I)=CST1
y(i,i)=cst1
300 continue
!
!     X(0)=Y(0)
!     X:(ID,K)
!xx      DO  400  I=1,ID
do 401 i=1,id
do 400 j=1,k
x(i,j)=y(i,j)
!xx  400 X1(1,I,J)=X(I,J)
x1(1,i,j)=x(i,j)
400 continue
401 continue
!
!     U(0)=(B,0)
!     B:(K,ID), 0:(K,K-ID)
u(1:k,1:k)=cst0
do 500 i=1,k
!xx      DO 590 J=1,K
!xx  590 U(I,J)=CST0
do 510 j=1,id
!xx  510 U(I,J)=B(I,J)
u(i,j)=b(i,j)
510 continue
500 continue
!
constm=0.1d-00
constd=0.0001d-00
lh=6
ih=0
ns=1
nsm1=ns-1
if(ig1.eq.0) go to 270
constd=0.01d-00
ig1=0
270 if(iswr.ge.1) go to 280
constd=0.01d-00
!
!     V(NS+1)=A*V(NS)+Z(NS+1)
!c  280 CALL SUBAWZ(V,AV,K,K,MJ1,MJ1)
280 call subawz(aw,iaw,idd,ir,v,av,k,k)
!xx      DO 230 I=1,K
do 231 i=1,k
do 230 j=1,k
!xx  230 V(I,J)=AV(I,J)
v(i,j)=av(i,j)
230 continue
231 continue
!
!     Y(NS)=C*V(NS) COMPUTATION
!     Y:(ID,K), C:(ID,K), V:(K,K)
!xx      DO 320 I=1,ID
do 321 i=1,id
do 320 j=1,k
!xx  320 Y(I,J)=V(I,J)
y(i,j)=v(i,j)
320 continue
321 continue
!
!     X(NS)=Y(NS)-C*A*U(NS-1) COMPUTATION
!c	CALL SUBAWZ(U,AU,K,K,MJ1,MJ1)
call subawz(aw,iaw,idd,ir,u,au,k,k)
!xx      DO 450 I=1,ID
do 451 i=1,id
do 450 j=1,k
!xx  450 X(I,J)=Y(I,J)-AU(I,J)
x(i,j)=y(i,j)-au(i,j)
450 continue
451 continue
!xx      DO 452 I=1,ID
do 453 i=1,id
do 452 j=1,k
!xx  452 X1(NS+1,I,J)=X(I,J)
x1(ns+1,i,j)=x(i,j)
452 continue
453 continue
!
!     U(NS)=A*U(NS-1)+B*X(NS) COMPUTATION
!xx      DO 180 I=1,K
do 182 i=1,k
do 180 j=1,k
sum=cst0
do 181 jj=1,id
!xx  181 SUM=SUM+B(I,JJ)*X(JJ,J)
sum=sum+b(i,jj)*x(jj,j)
181 continue
!xx  180 U(I,J)=AU(I,J)+SUM
u(i,j)=au(i,j)+sum
180 continue
182 continue
!
!     MAXIMUM ABSOLUTE VALUE OF X(I,J) SEARCH
!c	CALL SUBMAX(X,XMAX,ID,K,MJ2,MJ1)
call submax(x,xmax,id,k)
!
if(xmax.ge.constd) go to 740
ih=ih+1
if(ih.lt.lh) go to 720
go to 760
740 ih=0
720 if(ns.ge.49) go to 770
ns=ns+1
go to 280
770 m=49
!c	WRITE(6,2000) M,XMAX
!c	WRITE( 6,2001 )
!x      IF (IFG.NE.0) THEN
!x	 WRITE( LU,2000) M,XMAX
!x	 WRITE( LU,2001 )
!x 2001	 FORMAT( 1H ,'WARNING : INCOMPLETE CONVERGENCE' )
!x      END IF
ig1=1
if(iswr.ge.1) go to 771
ig1=0
!     RESCALING FOR A FEASIBLE INITIAL
if(itinp.ge.50) go to 771
itinp=itinp+1
!c	CALL RESCAL(K,ID)
!x      CALL RESCAL(B,AW,IAW,NH,IDD,IR,IJ,IK,K,ID,IFG,LU)
!xx      CALL RESCAL(B,AW,IAW,NH,IDD,IR,IJ,IK,K,ID)
call rescal(b,aw,iaw,nh,idd,ir,k,id)
!c	WRITE(6,2330) ITINP
!x      IF (IFG.NE.0) WRITE(LU,2330) ITINP
!xx 2330 FORMAT(//1H ,'ITINP=',I5)
go to 199
771 iswr=iswr+1
return
760 m=ns-lh
!c	WRITE(6,2000) M,XMAX
!x      IF (IFG.NE.0) WRITE(LU,2000) M,XMAX
!xx 2000 FORMAT(1H ,'M=',I5,5X,'XMAX1=',D12.5)
iswr=iswr+1
return
end
!
!c	SUBROUTINE NSUBX2(K,ID,KK,MJ1,MJ2,IG2)
!x      SUBROUTINE NSUBX2(B,AW,IAW,IDD,IR,K,ID,KK,IG2,Q,X2,IFG,LU)
subroutine nsubx2(b,aw,iaw,idd,ir,k,id,kk,ig2,q,x2)
  use timsac_kinds, only: dp
  implicit none
!     THIS SUBROUTINE COMPUTES X2(NS,I,J) AND U(NS),
!     NS=0,1,2,...,L(=KK).
!
!xx      IMPLICIT REAL*8(A-H,O-Z)
!c      COMMON /COM3/X2
!c	COMMON /COM60/IDD,IR
!c	COMMON /COM81/B
!c	COMMON /COM82/AW
!c      COMMON /COM2/Q
! common /COM99/
integer iswr
common /com99/iswr
!c      DIMENSION Q(50,10,5)
!c	DIMENSION IDD(10),IR(10)
!c      DIMENSION X2(50,5,5)
!c	DIMENSION AW(50),B(10,5)
!c	DIMENSION U(10,5)
!c	DIMENSION XX(5,5),AU(10,5)
!xx      DIMENSION Q(50,K,ID)
!xx      DIMENSION IDD(K),IR(K)
!xx      DIMENSION X2(50,ID,ID)
!xx      DIMENSION AW(IAW),B(K,ID)
!xx      DIMENSION U(K,ID)
!xx      DIMENSION XX(ID,ID),AU(K,ID)
integer iaw, k, id, kk, ig2, idd(k), ir(k)
real(dp) b(k,id), aw(iaw), q(50,k,id), x2(50,id,id)
! local
integer i, idp1, ih, ipd, j, jj, kmd, lh, ns, nsm1
real(dp) u(k,id), xx(id,id), xmax, au(k,id), cst0, cst1,&
&constd, sum
!
!     X(0)=I
!     X:(ID,ID)
cst0=0.0d-00
cst1=1.0d-00
x2(1,1:id,1:id)=cst0
do 200 i=1,id
!xx      DO 210 J=1,ID
!xx  210 X2(1,I,J)=CST0
!xx  200 X2(1,I,I)=CST1
x2(1,i,i)=cst1
200 continue
!
!     U(0)=B
!     U:(K,ID), B(K,ID)
!xx      DO 300 I=1,K
do 301 i=1,k
do 300 j=1,id
!xx  300 U(I,J)=B(I,J)
u(i,j)=b(i,j)
300 continue
301 continue
!
constd=0.0001d-00
lh=6
ih=0
ns=1
nsm1=ns-1
kmd = k - id
idp1 = id + 1
!xx      DO 191 I=1,KMD
do 192 i=1,kmd
ipd = i + id
do 191 j=1,id
!xx  191 Q(1,I,J) = U(IPD,J)
q(1,i,j) = u(ipd,j)
191 continue
192 continue
if(ig2.eq.0) go to 270
constd=0.01d-00
ig2=0
270 if(iswr.ge.2) go to 280
constd=0.01d-00
!
!     X(NS)=-C*A*U(NS-1) COMPUTATION
!c  280 CALL SUBAWZ(U,AU,K,ID,MJ1,MJ2)
280 call subawz(aw,iaw,idd,ir,u,au,k,id)
!xx      DO 450 I=1,ID
do 451 i=1,id
do 450 j=1,id
!xx  450 XX(I,J)=-AU(I,J)
xx(i,j)=-au(i,j)
450 continue
451 continue
!
!     U(NS)=A*U(NS-1)+B*X(NS) COMPUTATION
!xx      DO 180 I=1,K
do 182 i=1,k
do 180 j=1,id
sum=cst0
do 181 jj=1,id
!xx  181 SUM=SUM+B(I,JJ)*XX(JJ,J)
sum=sum+b(i,jj)*xx(jj,j)
181 continue
!xx  180 U(I,J)=AU(I,J)+SUM
u(i,j)=au(i,j)+sum
180 continue
182 continue
!xx      DO 190 I=1,KMD
do 193 i=1,kmd
ipd = i + id
do 190 j=1,id
!xx  190 Q(NS+1,I,J) = U(IPD,J)
q(ns+1,i,j) = u(ipd,j)
190 continue
193 continue
!xx      DO 530 I=1,ID
do 531 i=1,id
do 530 j=1,id
!xx  530 X2(NS+1,I,J)=XX(I,J)
x2(ns+1,i,j)=xx(i,j)
530 continue
531 continue
!
!     MAXIMUM ABSOLUTE VALUE OF XX(I,J) SEARCH
!c	CALL SUBMAX(XX,XMAX,ID,ID,MJ2,MJ2)
call submax(xx,xmax,id,id)
!
if(xmax.ge.constd) go to 740
ih=ih+1
if(ih.lt.lh) go to 720
go to 760
740 ih=0
720 if(ns.ge.49) go to 770
ns=ns+1
go to 280
770 kk=49
!c	WRITE(6,2200) KK,XMAX
!c	WRITE( 6,2300 )
!x      IF (IFG.NE.0) THEN
!x	 WRITE( LU,2200) KK,XMAX
!x	 WRITE( LU,2300 )
!x 2300	 FORMAT( 1H ,'WARNING : X2 INCOMPLETE CONVERGENCE' )
!x      END IF
ig2=1
return
760 kk=ns-lh
!c	WRITE(6,2200) KK,XMAX
!x      IF (IFG.NE.0) WRITE(LU,2200) KK,XMAX
!xx 2200 FORMAT(1H ,'KK=',I5,5X,'XMAX2=',D12.5)
return
end
!
!c	SUBROUTINE NSUBO(C0,K,ID,M,MJ1,MJ2)
subroutine nsubo(c0,k,id,m,o,x1)
  use timsac_kinds, only: dp
  implicit none
!     THIS SUBROUTINE COMPUTES O(NS).
!     O(NS)=(INVERSE OF C0)*X1(NS)
!xx      IMPLICIT REAL*8(A-H,O-Z)
!c      COMMON /COM15/X1
!c      COMMON /COM1/O
!c      DIMENSION X1(50,10,10),O(50,5,10)
!c	DIMENSION C0(5,5),C0I(5,5)
!c	DIMENSION X(5,10),OO(5,10)
!xx      DIMENSION X1(50,K,K),O(50,ID,K)
!xx      DIMENSION C0(ID,ID),C0I(ID,ID)
!xx      DIMENSION X(ID,K),OO(ID,K)
integer k, id, m
real(dp) c0(id,id), o(50,id,k), x1(50,k,k)
! local
integer i, j, mp1, ns, nsm1
real(dp) c0i(id,id), x(id,k), oo(id,k), c0det
!
!     INVERSE OF C0 COMPUTATION
!xx      DO 600 I=1,ID
do 601 i=1,id
do 600 j=1,id
!xx  600 C0I(I,J)=C0(I,J)
c0i(i,j)=c0(i,j)
600 continue
601 continue
!c	CALL INVDET(C0I,C0DET,ID,MJ2)
call invdet(c0i,c0det,id,id)
!
mp1=m+1
do 1000 ns=1,mp1
!xx      DO 610 J=1,K
do 611 j=1,k
do 610 i=1,id
!xx  610 X(I,J)=X1(NS,I,J)
x(i,j)=x1(ns,i,j)
610 continue
611 continue
!     O(NS)=(INVERSE OF C0)*X(NS) COMPUTATION
!c	CALL MULPLY(C0I,X,OO,ID,ID,K,MJ2,MJ2,MJ1)
call mulply(c0i,x,oo,id,id,k)
!xx      DO 620 I=1,ID
do 621 i=1,id
do 620 j=1,k
!xx  620 O(NS,I,J)=OO(I,J)
o(ns,i,j)=oo(i,j)
620 continue
621 continue
nsm1=ns-1
1000 continue
return
end
!
!c	SUBROUTINE SUBMAX(XX,XMAX,IA,IB,MJ6,MJ7)
subroutine submax(xx,xmax,ia,ib)
  use timsac_kinds, only: dp
  implicit none
!     MAXIMUM ABSOLUTE VALUE OF XX(I,J) SEARCH
!xx      IMPLICIT REAL*8(A-H,O-Z)
!c	DIMENSION XX(MJ6,MJ7)
!xx      DIMENSION XX(IA,IB)
integer ia, ib
real(dp) xx(ia,ib), xmax
! local
integer i, j
real(dp) t
!
xmax=0.0d-00
do 709 i=1,ia
do 710 j=1,ib
t=dabs(xx(i,j))
if(xmax.ge.t) go to 710
xmax=t
710 continue
709 continue
return
end
!
!c	SUBROUTINE GRAD(GR,M,K,ID,IAW,IPQ,MJ1,MJ2)
!x      SUBROUTINE GRAD(IJ,IK,GR,M,K,ID,IAW,IPQ,MJ4,O,CXX1,CXV1,IFG,LU)
!xx      SUBROUTINE GRAD(IJ,IK,GR,M,K,ID,IAW,IPQ,MJ4,O,CXX1,CXV1)
subroutine grad(ij,ik,gr,m,k,id,iaw,mj4,o,cxx1,cxv1)
  use timsac_kinds, only: dp
  implicit none
!     THIS SUBROUTINE COMPUTES GRADIENT (GA AND GB).
!     GA=-2.0*(O(0)'*CXV1(1)+O(1)'*CXV1(2)+O(2)'*CXV1(3)+...+O(NSL)'*CXV
!     1(1+NSL)).
!     GB=-2.0*(O(0)'*CXX1(0)+O(1)'*CXX1(1)+O(2)'*CXX1(2)+...+O(NSL)'*CXX
!     1(NSL)).
!
!xx      IMPLICIT REAL*8(A-H,O-Z)
!c      COMMON /COM1/O
!c      COMMON /COM12/CXX1
!c      COMMON /COM14/CXV1
!c      COMMON/COM61/IJ,IK
!      DIMENSION O(50,5,10),CXX1(51,5,5),CXV1(51,5,10)
!c	DIMENSION OO(10,5),CXX(5,5),CXV(5,10)
!c	DIMENSION ZA(10,10),ZZA(10,10),ZB(10,5),ZZB(10,5)
!c	DIMENSION GR(100)
!c	DIMENSION ZAW(50)
!c	DIMENSION IJ(5),IK(5)
!xx      DIMENSION O(50,ID,K),CXX1(51,ID,ID),CXV1(51,ID,K)
!xx      DIMENSION OO(K,ID),CXX(ID,ID),CXV(ID,K)
!xx      DIMENSION ZA(K,K),ZZA(K,K),ZB(K,ID),ZZB(K,ID)
!xx      DIMENSION GR(MJ4),ZAW(MJ4)
!xx      DIMENSION IJ(ID),IK(ID)
integer m, k, id, iaw, mj4, ij(id), ik(id)
real(dp) gr(mj4), o(50,id,k), cxx1(51,id,id),&
&cxv1(51,id,k)
! local
integer i, idp1, ijk, iksj, j, jj, jl, ks, mp1,  nsp1, nsp2
real(dp) oo(k,id), cxx(id,id), cxv(id,k), za(k,k),&
&zza(k,k), zb(k,id), zzb(k,id), zaw(mj4),&
&cst0, cst1, cst2
!
!
!     GA, GB COMPUTATION
cst0=0.0d-00
cst1=1.0d-00
cst2=2.0d-00
!xx      DO 210 I=1,K
!xx      DO 208 J=1,K
!xx  208 ZA(I,J)=CST0
!xx      DO 209 J=1,ID
!xx  209 ZB(I,J)=CST0
!xx  210 CONTINUE
za(1:k,1:k)=cst0
zb(1:k,1:id)=cst0
mp1=m+1
do 300 nsp1=1,mp1
nsp2=nsp1+1
!xx      DO 310 I=1,K
do 311 i=1,k
do 310 j=1,id
oo(i,j)=o(nsp1,j,i)
!xx  310 CXV(J,I)=CXV1(NSP2,J,I)
cxv(j,i)=cxv1(nsp2,j,i)
310 continue
311 continue
!xx      DO 320 I=1,ID
do 321 i=1,id
do 320 j=1,id
!xx  320 CXX(I,J)=CXX1(NSP1,I,J)
cxx(i,j)=cxx1(nsp1,i,j)
320 continue
321 continue
!c	CALL MULPLY(OO,CXV,ZZA,K,ID,K,MJ1,MJ2,MJ1)
!c	CALL MATADL(ZA,ZZA,K,K,MJ1,MJ1)
!c	CALL MULPLY(OO,CXX,ZZB,K,ID,ID,MJ1,MJ2,MJ2)
!c	CALL MATADL(ZB,ZZB,K,ID,MJ1,MJ2)
call mulply(oo,cxv,zza,k,id,k)
call matadl(za,zza,k,k)
call mulply(oo,cxx,zzb,k,id,id)
call matadl(zb,zzb,k,id)
300 continue
!     GRADIENT ARRANGEMENT
!c	IKSUM=0
iksj=0
do 1300 ks=1,id
jl=ik(ks)
do 1310 j=1,jl
!c	IKSJ=IKSUM+J
iksj=iksj+1
ijk=ij(ks)
!xx 1310 ZAW(IKSJ)=ZA(IJK,J)
zaw(iksj)=za(ijk,j)
1310 continue
!c 1300 IKSUM=IKSUM+JL
1300 continue
iaw=iksj
do 500 i=1,iaw
!xx  500 GR(I)=-CST2*ZAW(I)
gr(i)=-cst2*zaw(i)
500 continue
!     K>ID
idp1=id+1
jj=iaw
!xx      DO 510 I=IDP1,K
do 511 i=idp1,k
do 510 j=1,id
jj=jj+1
!xx  510 GR(JJ)=-CST2*ZB(I,J)
gr(jj)=-cst2*zb(i,j)
510 continue
511 continue
!c	WRITE(6,1400)
!c	CALL SUBVCP(GR,IPQ)
!x      IF (IFG.NE.0) WRITE(LU,1400)
!x      IF (IFG.NE.0) CALL SUBVCP(GR,IPQ,LU)
return
!xx 1400 FORMAT(/1H ,'GRADIENT')
end
!
!c	SUBROUTINE SBCXY1(M,L,ID,MJ2)
subroutine sbcxy1(cyy1,mj3,m,l,id,x2,cxy1)
  use timsac_kinds, only: dp
  implicit none
!     THIS SUBROUTINE COMPUTES CXY1(LS) LS=0,M+1+L.
!     CXY1(LS)=X2(0)*CYY1(LS)+X2(1)*CYY1(LS-1)+X2(2)*CYY1(LS-2)+...
!     +X2(L)*CYY1(LS-L).
!xx      IMPLICIT REAL*8(A-H,O-Z)
!c      COMMON /COM3/X2
!c	COMMON /COM10/CYY1
!c      COMMON /COM11/CXY1
!c	DIMENSION X2(50,5,5),CYY1(160,5,5),CXY1(100,5,5)
!c	DIMENSION XX(5,5),CYY(5,5),Z(5,5),ZZ(5,5)
!xx      DIMENSION X2(50,ID,ID),CYY1(MJ3,ID,ID),CXY1(100,ID,ID)
!xx      DIMENSION XX(ID,ID),CYY(ID,ID),Z(ID,ID),ZZ(ID,ID)
integer mj3, m, l, id
real(dp) cyy1(mj3,id,id), x2(50,id,id), cxy1(100,id,id)
! local
integer i, j, lms, lmsp1, lp1, ls, lsl, lsp1, ms, msp1
real(dp) xx(id,id), cyy(id,id), z(id,id), zz(id,id), cst0
!
!
cst0=0.0d-00
lsl=m+2+l
lp1=l+1
!
do 200 lsp1=1,lsl
ls=lsp1-1
!xx      DO 209 I=1,ID
!xx      DO 209 J=1,ID
!xx  209 Z(I,J)=CST0
z(1:id,1:id)=cst0
do 300 msp1=1,lp1
ms=msp1-1
!xx      DO 310 I=1,ID
do 311 i=1,id
do 310 j=1,id
!xx  310 XX(I,J)=X2(MSP1,I,J)
xx(i,j)=x2(msp1,i,j)
310 continue
311 continue
lms=ls-ms
if(lms.lt.0) go to 319
lmsp1=lms+1
!xx      DO 320 I=1,ID
do 321 i=1,id
do 320 j=1,id
!xx  320 CYY(I,J)=CYY1(LMSP1,I,J)
cyy(i,j)=cyy1(lmsp1,i,j)
320 continue
321 continue
go to 340
319 lmsp1=-lms+1
!xx      DO 330 I=1,ID
do 331 i=1,id
do 330 j=1,id
!xx  330 CYY(I,J)=CYY1(LMSP1,J,I)
cyy(i,j)=cyy1(lmsp1,j,i)
330 continue
331 continue
!c  340 CALL MULPLY(XX,CYY,ZZ,ID,ID,ID,MJ2,MJ2,MJ2)
!c	CALL MATADL(Z,ZZ,ID,ID,MJ2,MJ2)
340 call mulply(xx,cyy,zz,id,id,id)
call matadl(z,zz,id,id)
300 continue
!xx  280 DO 350 I=1,ID
!xx      DO 350 I=1,ID
do 351 i=1,id
do 350 j=1,id
!xx  350 CXY1(LSP1,I,J)=Z(I,J)
cxy1(lsp1,i,j)=z(i,j)
350 continue
351 continue
200 continue
return
end
!
!c	SUBROUTINE SUBCXX(C0,M,KK,ID,MJ2)
subroutine subcxx(c0,m,kk,id,x2,cxy1,cxx1)
  use timsac_kinds, only: dp
  implicit none
!     THIS SUBROUTINE COMPUTES CXX1(NS) NS=0,M+1.
!     CXX1(NS)=CXY1(NS)*X2(0)'+CXY1(NS+1)*X2(1)'+CXY1(NS+2)*X2(2)'+...
!     +CXY1(NS+MSL)*X2(MSL)'.
!xx      IMPLICIT REAL*8(A-H,O-Z)
!c      COMMON /COM3/X2
!c      COMMON /COM11/CXY1
!c      COMMON /COM12/CXX1
!c      DIMENSION X2(50,5,5),CXY1(100,5,5),CXX1(51,5,5)
!c	DIMENSION XX(5,5),CXY(5,5),Z(5,5),ZZ(5,5)
!c	DIMENSION C0(5,5)
!xx      DIMENSION X2(50,ID,ID),CXY1(100,ID,ID),CXX1(51,ID,ID)
!xx      DIMENSION XX(ID,ID),CXY(ID,ID),Z(ID,ID),ZZ(ID,ID)
!xx      DIMENSION C0(ID,ID)
integer m, kk, id
real(dp) c0(id,id), x2(50,id,id), cxy1(100,id,id),&
&cxx1(51,id,id)
! local
integer i, j, kkp1, msp1, ns, nsl, nsp1, nspms1
real(dp) xx(id,id), cxy(id,id), z(id,id), zz(id,id), cst0
!
!
cst0=0.0d-00
nsl=m+2
kkp1=kk+1
!
do 200 nsp1=1,nsl
ns=nsp1-1
!xx      DO 209 I=1,ID
!xx      DO 209 J=1,ID
!xx  209 Z(I,J)=CST0
z(1:id,1:id)=cst0
do 300 msp1=1,kkp1
nspms1=ns+msp1
!xx      DO 310 I=1,ID
do 311 i=1,id
do 310 j=1,id
xx(i,j)=x2(msp1,j,i)
!xx  310 CXY(I,J)=CXY1(NSPMS1,I,J)
cxy(i,j)=cxy1(nspms1,i,j)
310 continue
311 continue
!c	CALL MULPLY(CXY,XX,ZZ,ID,ID,ID,MJ2,MJ2,MJ2)
!c	CALL MATADL(Z,ZZ,ID,ID,MJ2,MJ2)
call mulply(cxy,xx,zz,id,id,id)
call matadl(z,zz,id,id)
300 continue
!xx  280 DO 350 I=1,ID
!xx      DO 350 I=1,ID
do 351 i=1,id
do 350 j=1,id
!xx  350 CXX1(NSP1,I,J)=Z(I,J)
cxx1(nsp1,i,j)=z(i,j)
350 continue
351 continue
200 continue
!xx      DO 400 I=1,ID
do 401 i=1,id
do 400 j=1,id
!xx  400 C0(I,J)=CXX1(1,I,J)
c0(i,j)=cxx1(1,i,j)
400 continue
401 continue
return
end
!
!c	SUBROUTINE GCXV1(M,L,K,ID,MJ1,MJ2)
subroutine gcxv1(m,l,k,id,q,cxy1,cxv1)
  use timsac_kinds, only: dp
  implicit none
!     THIS SUBROUTINE COMPUTES CXV1(MS) MS=0,M+1.
!     CXV1(MS)=CXY1(MS)*Q(0)'+CXY1(MS+1)*Q(1)'+...+CXY1(MS+L)*Q(L)'
!xx      IMPLICIT REAL*8(A-H,O-Z)
!c      COMMON /COM2/Q
!c      COMMON /COM11/CXY1
!c      COMMON /COM14/CXV1
!c      DIMENSION Q(50,10,5),CXY1(100,5,5),CXV1(51,5,10)
!c	DIMENSION QQ(5,10),CXY(5,5),Z(5,10),ZZ(5,10)
!xx      DIMENSION Q(50,K,ID),CXY1(100,ID,ID),CXV1(51,ID,K)
!xx      DIMENSION QQ(ID,K-ID),CXY(ID,ID),Z(ID,K-ID),ZZ(ID,K-ID)
integer m, l, k, id
real(dp) q(50,k,id), cxy1(100,id,id), cxv1(51,id,k)
! local
integer i, idp1, j, jmd, kmd, lp1, mnsp1, mp2, ms, msp1, nsp1
real(dp) qq(id,k-id), cxy(id,id), z(id,k-id), zz(id,k-id),&
&cst0
!
!
cst0=0.0d-00
mp2=m+2
lp1=l+1
kmd = k - id
idp1 = id + 1
!
do 200 msp1=1,mp2
ms=msp1-1
!xx      DO 209 I=1,ID
!c	DO 209 J=1,K
!xx      DO 209 J=1,KMD
!xx  209 Z(I,J)=CST0
z(1:id,1:kmd)=cst0
do 300 nsp1=1,lp1
mnsp1=ms+nsp1
do 210 i=1,id
do 220 j=1,id
!xx  220 CXY(I,J)=CXY1(MNSP1,I,J)
cxy(i,j)=cxy1(mnsp1,i,j)
220 continue
do  230  j=1,kmd
!xx  230 QQ(I,J)=Q(NSP1,J,I)
qq(i,j)=q(nsp1,j,i)
230 continue
210 continue
!c	CALL  MULPLY ( CXY,QQ,ZZ,ID,ID,KMD,MJ2,MJ2,MJ1 )
!c	CALL  MATADL ( Z,ZZ,ID,KMD,MJ2,MJ1 )
call  mulply ( cxy,qq,zz,id,id,kmd )
call  matadl ( z,zz,id,kmd )
300 continue
!xx  280 DO 350 I=1,ID
!xx      DO 350 I=1,ID
do 351 i=1,id
do 349   j=1,id
!xx  349 CXV1(MSP1,I,J)=CXY1(MSP1,I,J)
cxv1(msp1,i,j)=cxy1(msp1,i,j)
349 continue
do 350   j=idp1,k
jmd = j - id
!xx  350 CXV1(MSP1,I,J) = Z(I,JMD)
cxv1(msp1,i,j) = z(i,jmd)
350 continue
351 continue
200 continue
return
end
!
!
!c	SUBROUTINE  SUBIDR ( K,ID,IAW )
subroutine  subidr ( nh,idd,ir,ij,ik,k,id,iaw )
  use timsac_kinds, only: dp
  implicit none
!     CONSTRUCTION OF CHARACTERISTIC VECTORS
!c	COMMON /COM60/IDD,IR
!c	COMMON /COM61/IJ,IK
!c	COMMON /COM70/NH
!c	DIMENSION NH(10),IDD(10),IR(10)
!c	DIMENSION IJ(5),IK(5)
!xx      DIMENSION NH(K),IDD(K),IR(K)
!xx      DIMENSION IJ(ID),IK(ID)
integer k, id, iaw, nh(k), idd(k), ir(k), ij(id), ik(id)
! local
integer i, idh, ifp, j, lr
!
iaw = 0
ifp = 0
do 100 i=1,k
idh = nh(i) + id
if( nh(k) .lt. idh )   go to 111
lr=0
j=i+1
130 lr=lr+1
if(nh(j).lt.idh) go to 110
if(nh(j).gt.idh) go to 120
idd(i)=0
ir(i)=i+lr
go to 100
120 ir(i) = i+ lr - 1
go to 113
110 j=j+1
go to 130
111 ir(i) = k
113 idd(i) = 1
!xx  112 IFP = IFP + 1
ifp = ifp + 1
ij(ifp) = i
ik(ifp) = ir(i)
iaw = iaw + ir(i)
100 continue
return
end
!
!c	SUBROUTINE SUBA(K,ID,MJ1,MJ2,IAW)
subroutine suba(a,aw,ij,ik,k,id,iaw)
  use timsac_kinds, only: dp
  implicit none
!     CONSTRUCTION OF A FROM AW
!xx      IMPLICIT REAL*8(A)
!c	COMMON /COM61/IJ,IK
!c	COMMON /COM80/A
!c	COMMON /COM82/AW
!c	DIMENSION A(10,10),AW(50)
!c	DIMENSION IJ(5),IK(5)
!xx      DIMENSION A(K,K),AW(IAW)
!xx      DIMENSION IJ(ID),IK(ID)
integer k, id, iaw, ij(id), ik(id)
real(dp) a(k,k), aw(iaw)
! local
integer ijk, iksj, j, jl, ks
!
!c	IKSUM=0
iksj=0
do 300 ks=1,id
jl=ik(ks)
do 310 j=1,jl
!c	IKSJ=IKSUM+J
iksj=iksj+1
ijk=ij(ks)
!xx  310 A(IJK,J)=AW(IKSJ)
a(ijk,j)=aw(iksj)
310 continue
!c  300 IKSUM=IKSUM+JL
300 continue
iaw=iksj
return
end
!
!c	SUBROUTINE SUBAWZ(Z,ZP,K,IC,MJ3,MJ4)
subroutine subawz(aw,iaw,idd,ir,z,zp,k,ic)
  use timsac_kinds, only: dp
  implicit none
!     ZP=A*Z
!xx      IMPLICIT REAL*8(A-H,O-Z)
!c	COMMON /COM60/IDD,IR
!c	COMMON /COM82/AW
!c	DIMENSION AW(50),IDD(10),IR(10)
!c	DIMENSION Z(MJ3,MJ4),ZP(MJ3,MJ4)
!xx      DIMENSION AW(IAW),IDD(K),IR(K)
!xx      DIMENSION Z(K,IC),ZP(K,IC)
integer k, ic, iaw, idd(k), ir(k)
real(dp) aw(iaw), z(k,ic), zp(k,ic)
! local
integer i, iri, j, jj, ni, ni0, ni0j
real(dp) cst0, sum
cst0=0.0d-00
ni=0
do 400 i=1,k
if(idd(i).ne.0) go to 410
iri=ir(i)
do 420 jj=1,ic
!xx  420 ZP(I,JJ)=Z(IRI,JJ)
zp(i,jj)=z(iri,jj)
420 continue
go to 400
410 ni0=ni
ni=ni+ir(i)
iri=ir(i)
do 430 jj=1,ic
sum=cst0
do 440 j=1,iri
ni0j=ni0+j
!xx  440 SUM=SUM+AW(NI0J)*Z(J,JJ)
sum=sum+aw(ni0j)*z(j,jj)
440 continue
!xx  430 ZP(I,JJ)=SUM
zp(i,jj)=sum
430 continue
400 continue
return
end
!
!c	SUBROUTINE SUBCM(M,L,K,ID,MJ1,MJ2)
!xx      SUBROUTINE SUBCM(M,L,K,ID,O,CM1)
subroutine subcm(m,k,id,o,cm1)
  use timsac_kinds, only: dp
  implicit none
!     THIS SUBROUTINE COMPUTES CM1(MS).
!xx      IMPLICIT REAL*8(A-H,O-Z)
!c      COMMON /COM1/O
!c      COMMON /COM15/CM1
!c      DIMENSION O(50,5,10) , CM1(50,10,10)
!c	DIMENSION OO(10,5),XX(5,10),Z(10,10),ZZ(10,10)
!xx      DIMENSION O(50,ID,K),CM1(50,K,K)
!xx      DIMENSION OO(K,ID),XX(ID,K),Z(K,K),ZZ(K,K)
integer m, k, id
real(dp) o(50,id,k), cm1(50,k,k)
! local
integer i, j, mnsp1, mp1, ms, msp1, nsl, nsp1
real(dp) oo(k,id), xx(id,k), z(k,k), zz(k,k), cst0
!
!
!     CM1(MS) COMPUTATION
!
cst0=0.0d-00
mp1=m+1
do 200 msp1=1,mp1
ms=msp1-1
!xx      DO 209 I=1,K
!xx      DO 209 J=1,K
!xx  209 Z(I,J)=CST0
z(1:k,1:k)=cst0
nsl=m-ms+1
do 300 nsp1=1,nsl
mnsp1=ms+nsp1
!xx      DO 220 I=1,K
do 221 i=1,k
do 220 j=1,id
oo(i,j)=o(nsp1,j,i)
!xx  220 XX(J,I) = CM1(MNSP1,J,I)
xx(j,i) = cm1(mnsp1,j,i)
220 continue
221 continue
!c	CALL MULPLY(OO,XX,ZZ,K,ID,K,MJ1,MJ2,MJ1)
!c	CALL MATADL(Z,ZZ,K,K,MJ1,MJ1)
call mulply(oo,xx,zz,k,id,k)
call matadl(z,zz,k,k)
300 continue
!xx      DO 350 I=1,K
do 351 i=1,k
do 350 j=1,k
!xx  350 CM1(MSP1,I,J)=Z(I,J)
cm1(msp1,i,j)=z(i,j)
350 continue
351 continue
200 continue
return
end
!
!c	SUBROUTINE SUBHES(M,L,K,ID,IAW,MJ1,MJ2)
subroutine subhes(cyy1,ij,ik,hs,m,l,k,id,iaw,mj3,mj4,o,q,cm1,x2,&
&cxy1,cxx1,cxv1)
  use timsac_kinds, only: dp
  implicit none
!     THIS SUBROUTINE COMPUTES HESSIAN.
!xx      IMPLICIT REAL*8(A-H,O-Z)
!c      COMMON /COM61/IJ,IK
!c      COMMON /COM50/HS
!c      COMMON /COM15/CM1
!c      COMMON /COM12/CXX1
!c      COMMON /COM16/CVV1
!c      COMMON /COM19/CXV2
!c      COMMON /COM14/CXV1
!c      DIMENSION CYY1(160,5,5),IJ(5),IK(5)
!c      DIMENSION CXV1(51,5,10)
!c      DIMENSION CM1(50,10,10),CXX1(51,5,5),CVV1(50,10,10)
!c      DIMENSION HS(50,50)
!xx      DIMENSION CYY1(MJ3,ID,ID),IJ(ID),IK(ID)
!xx      DIMENSION CXV1(51,ID,K)
!xx      DIMENSION CM1(50,K,K),CXX1(51,ID,ID),CVV1(50,K,K)
!xx      DIMENSION CXV2(51,ID,K)
!xx      DIMENSION HS(MJ4,MJ4)
!xx      DIMENSION O(50,ID,K),Q(50,K,ID),X2(50,ID,ID)
!xx      DIMENSION CXY1(100,ID,ID),CXY2(50,ID,ID)
integer m, l, k, id, iaw, mj3, mj4, ij(id), ik(id)
real(dp) cyy1(mj3,id,id), hs(mj4,mj4), o(50,id,k),&
&q(50,k,id), cm1(50,k,k), x2(50,id,id),&
&cxy1(100,id,id), cxx1(51,id,id), cxv1(51,id,k)
! local
integer i, i2, idp1, ijk, ijk2, iksj, iksj2, j, j2, jj, jj2, jl,&
&jl2, ks, ks2, ms, mp1, msp1, msp2
real(dp) cvv1(50,k,k), cxv2(51,id,k), cxy2(50,id,id),&
&cst2, sumbb, sumab, sumaa
!
cst2=2.0d-00
!     CXY2(MS) MS=0,M-1 COMPUTATION
!c	CALL SBCXY2(M,L,ID,MJ2)
call sbcxy2(cyy1,mj3,m,l,id,x2,cxy2)
!     CXV2(MS) MS=0,M-1 COMPUTATION
!c	CALL HCXV2(M,L,K,ID,MJ1,MJ2)
call hcxv2(m,l,k,id,q,cxy1,cxy2,cxv2)
!     CM1(MS) MS=0,M COMPUTATION
!c	CALL SUBCM(M,L,K,ID,MJ1,MJ2)
!xx      CALL SUBCM(M,L,K,ID,O,CM1)
call subcm(m,k,id,o,cm1)
!     HESSIAN COMPUTATION
!     HBB COMPUTATION
mp1=m+1
idp1=id+1
jj=iaw
!xx      DO 460 I=IDP1,K
do 461 i=idp1,k
do 460 j=1,id
jj=jj+1
jj2=iaw
!xx      DO 1460 I2=IDP1,I
do 1461 i2=idp1,i
do 1460 j2=1,id
jj2=jj2+1
sumbb=cm1(1,i,i2)*cxx1(1,j,j2)
do 490 msp1=2,mp1
!xx  490 SUMBB=SUMBB+CM1(MSP1,I,I2)*CXX1(MSP1,J,J2)+CM1(MSP1,I2,I)*CXX1(MSP
!xx     A1,J2,J)
sumbb=sumbb+cm1(msp1,i,i2)*cxx1(msp1,j,j2)+cm1(msp1,i2,i)*cxx1(msp&
&1,j2,j)
490 continue
hs(jj,jj2)=cst2*sumbb
hs(jj2,jj)=hs(jj,jj2)
1460 continue
1461 continue
460 continue
461 continue
!
!     HAB COMPUTATION
iksj = 0
do 800 ks=1,id
jl=ik(ks)
do 810 j=1,jl
iksj = iksj + 1
ijk=ij(ks)
jj2=iaw
!xx      DO 860 I2=IDP1,K
do 861 i2=idp1,k
do 860 j2=1,id
jj2=jj2+1
sumab=cm1(1,ijk,i2)*cxv1(2,j2,j)
do 870 ms=1,m
msp1=ms+1
msp2=msp1+1
sumab=sumab+cm1(msp1,ijk,i2)*cxv2(ms,j2,j)+cm1(msp1,i2,ijk)*cxv1(m&
&sp2,j2,j)
870 continue
hs(iksj,jj2)=cst2*sumab
hs(jj2,iksj)=hs(iksj,jj2)
860 continue
861 continue
810 continue
800 continue
!
!
!     CYV1,CYV2 COMPUTATION
!c	CALL SBCYV1(M,L,K,ID,MJ1,MJ2)
!c	CALL SBCYV2(L,K,ID,MJ1,MJ2)
call sbcyv1(cyy1,mj3,m,l,k,id,q,cxv1)
call sbcyv2(cyy1,mj3,l,k,id,q,cxv2)
!     CVV1 COMPUTATION
!c	CALL SUBCVV(M,L,K,ID,MJ1,MJ2)
call subcvv(m,l,k,id,q,cxv1,cvv1,cxv2)
!
!     HAA COMPUTATION
iksj = 0
!xx      DO 310 KS=1,ID
do 311 ks=1,id
jl=ik(ks)
do 310 j=1,jl
iksj = iksj + 1
ijk=ij(ks)
iksj2 = 0
!xx      DO 1310 KS2=1,KS
do 1311 ks2=1,ks
jl2=ik(ks2)
do 1310 j2=1,jl2
iksj2=iksj2+1
ijk2=ij(ks2)
sumaa=cm1(1,ijk,ijk2)*cvv1(1,j,j2)
do 440 msp1=2,mp1
!xx  440 SUMAA=SUMAA+CM1(MSP1,IJK,IJK2)*CVV1(MSP1,J,J2)+CM1(MSP1,IJK2,IJK)*
!xx     ACVV1(MSP1,J2,J)
sumaa=sumaa+cm1(msp1,ijk,ijk2)*cvv1(msp1,j,j2)+cm1(msp1,ijk2,ijk)*&
&cvv1(msp1,j2,j)
440 continue
hs(iksj,iksj2)=cst2*sumaa
hs(iksj2,iksj)=hs(iksj,iksj2)
1310 continue
1311 continue
310 continue
311 continue
return
end
!
!
!c	SUBROUTINE SUBCVV(M,L,K,ID,MJ1,MJ2)
subroutine subcvv(m,l,k,id,q,cyv1,cvv1,cyv2)
  use timsac_kinds, only: dp
  implicit none
!     THIS SUBROUTINE COMPUTES CVV1(MS) MS=0,M.
!     CVV1(MS)=Q(0)*CYV(MS)+Q(1)*CYV(MS-1)+...+Q(L)*CYV(MS-L)
!xx      IMPLICIT REAL*8(A-H,O-Z)
!c      COMMON /COM2/Q
!c      COMMON /COM14/ CYV1
!c      COMMON /COM16/ CVV1
!c      COMMON /COM19/ CYV2
!c      DIMENSION Q(50,10,5),CYV1(51,5,10),CVV1(50,10,10)
!xx      DIMENSION Q(50,K,ID),CYV1(51,ID,K),CVV1(50,K,K)
!xx      DIMENSION CYV2(51,ID,K)
!c	DIMENSION QQ(10,5),CXV(5,10),Z(10,10),ZZ(10,10)
!xx      DIMENSION QQ(K-ID,ID),CXV(ID,K),Z(K-ID,K),ZZ(K-ID,K)
integer m, l, k, id
real(dp) q(50,k,id), cyv1(51,id,k), cvv1(50,k,k),&
&cyv2(51,id,k)
! local
integer i, idp1, imd, j, kmd, lp1, lsp1, ml, mp1, msp1
real(dp) qq(k-id,id), cxv(id,k), z(k-id,k), zz(k-id,k),&
&cst0, cst1
!
!
!
cst0=0.0d-00
cst1=1.0d-00
!     CVV1(MS) COMPUTATION
mp1=m+1
lp1=l+1
idp1 = id + 1
kmd = k - id
do 200 msp1=1,mp1
!xx      DO 209 I=1,KMD
!xx      DO 209 J=1,K
!xx  209 Z(I,J)=CST0
z(1:kmd,1:k)=cst0
do 300 lsp1=1,lp1
!xx      DO 210 I=1,KMD
do 211 i=1,kmd
do 210 j=1,id
!xx  210 QQ(I,J)=Q(LSP1,I,J)
qq(i,j)=q(lsp1,i,j)
210 continue
211 continue
ml=msp1-lsp1
if(ml.lt.0) go to 319
!xx      DO 320 I=1,ID
do 321 i=1,id
do 320 j=1,k
!xx  320 CXV(I,J)=CYV1(ML+1,I,J)
cxv(i,j)=cyv1(ml+1,i,j)
320 continue
321 continue
go to 340
319 ml=-ml
!xx      DO 330 I=1,ID
do 331 i=1,id
do 330 j=1,k
!xx  330 CXV(I,J)=CYV2(ML+1,I,J)
cxv(i,j)=cyv2(ml+1,i,j)
330 continue
331 continue
!c  340 CALL  MULPLY ( QQ,CXV,ZZ,KMD,ID,K,MJ1,MJ2,MJ1 )
!c	CALL  MATADL ( Z,ZZ,KMD,K,MJ1,MJ1 )
340 call  mulply ( qq,cxv,zz,kmd,id,k )
call  matadl ( z,zz,kmd,k )
300 continue
do 349   i=1,id
!xx      DO 349   J=1,K
do 348   j=1,k
!xx  349 CVV1(MSP1,I,J)=CYV1(MSP1,I,J)
cvv1(msp1,i,j)=cyv1(msp1,i,j)
348 continue
349 continue
!xx      DO 350   I=IDP1,K
do 351   i=idp1,k
do 350   j=1,k
imd = i - id
!xx  350 CVV1(MSP1,I,J) = Z(IMD,J)
cvv1(msp1,i,j) = z(imd,j)
350 continue
351 continue
200 continue
return
end
!
!c	SUBROUTINE SBCYV1(M,L,K,ID,MJ1,MJ2)
subroutine sbcyv1(cyy1,mj3,m,l,k,id,q,cyv1)
  use timsac_kinds, only: dp
  implicit none
!     THIS SUBROUTINE COMPUTES CYV1(MS) (MS=0,M).
!     CYV1(MS)=CYY1(MS)*Q(0)'+CYY1(MS+1)*Q(1)'+CYY1(MS+2)*Q(2)'+...
!     -CYY1(MS+L)*Q(L)'.
!xx      IMPLICIT REAL*8(A-H,O-Z)
!c      COMMON /COM2/Q
!c      COMMON /COM10/CYY1
!c      COMMON /COM14/CYV1
!c	DIMENSION Q(50,10,5),CYY1(160,5,5),CYV1(51,5,10)
!c	DIMENSION QQ(5,10),CYY(5,5),Z(5,10),ZZ(5,10)
!xx      DIMENSION Q(50,K,ID),CYY1(MJ3,ID,ID),CYV1(51,ID,K)
!xx      DIMENSION QQ(ID,K-ID),CYY(ID,ID),Z(ID,K-ID),ZZ(ID,K-ID)
integer mj3, m, l, k, id
real(dp) cyy1(mj3,id,id), q(50,k,id), cyv1(51,id,k)
! local
integer i, idp1, j, jmd, kmd, lp1, lsp1, ml, mp1, msp1
real(dp) qq(id,k-id), cyy(id,id), z(id,k-id), zz(id,k-id),&
&cst0
!
cst0=0.0d-00
mp1=m+1
lp1=l+1
idp1 = id + 1
kmd = k - id
do 200 msp1=1,mp1
!xx      DO 209 I=1,ID
!c	DO 209 J=1,K
!xx      DO 209 J=1,KMD
!xx  209 Z(I,J)=CST0
z(1:id,1:kmd)=cst0
do 300 lsp1=1,lp1
ml=msp1+lsp1-1
!xx      DO 210 I=1,ID
do 211 i=1,id
do 210 j=1,id
!xx  210 CYY(I,J)=CYY1(ML,I,J)
cyy(i,j)=cyy1(ml,i,j)
210 continue
211 continue
!xx      DO 220 I=1,ID
do 221 i=1,id
do 220 j=1,kmd
!xx  220 QQ(I,J)=Q(LSP1,J,I)
qq(i,j)=q(lsp1,j,i)
220 continue
221 continue
!c	CALL  MULPLY ( CYY,QQ,ZZ,ID,ID,KMD,MJ2,MJ2,MJ1 )
!c	CALL  MATADL ( Z,ZZ,ID,KMD,MJ2,MJ1 )
call  mulply ( cyy,qq,zz,id,id,kmd )
call  matadl ( z,zz,id,kmd )
300 continue
!xx      DO 350 I=1,ID
do 351 i=1,id
do 349 j=1,id
!xx  349 CYV1(MSP1,I,J)=CYY1(MSP1,I,J)
cyv1(msp1,i,j)=cyy1(msp1,i,j)
349 continue
do 350 j=idp1,k
jmd = j - id
!xx  350 CYV1(MSP1,I,J)=Z(I,JMD)
cyv1(msp1,i,j)=z(i,jmd)
350 continue
351 continue
200 continue
return
end
!
!c	SUBROUTINE SBCYV2(L,K,ID,MJ1,MJ2)
subroutine sbcyv2(cyy1,mj3,l,k,id,q,cyv2)
  use timsac_kinds, only: dp
  implicit none
!     THIS SUBROUTINE COMPUTES CYV2(MS) (MS=-L,0).
!     CYV2(MS)=CYV1(-MS)
!     =CYY1(-MS)*Q(0)'+CYY1(-MS+1)*Q(1)'+CYY1(-MS+2)*Q(2)'+...
!     +CYY1(-MS+L)*Q(L)'.
!xx      IMPLICIT REAL*8(A-H,O-Z)
!c      COMMON /COM2/Q
!c      COMMON /COM10/CYY1
!c      COMMON /COM19/ CYV2
!c	DIMENSION Q(50,10,5),CYY1(160,5,5),CYV2(51,5,10)
!c	DIMENSION QQ(5,10),CYY(5,5),Z(5,10),ZZ(5,10)
!xx      DIMENSION Q(50,K,ID),CYY1(MJ3,ID,ID),CYV2(51,ID,K)
!xx      DIMENSION QQ(ID,K-ID),CYY(ID,ID),Z(ID,K-ID),ZZ(ID,K-ID)
integer mj3, l, k, id
real(dp) cyy1(mj3,id,id), q(50,k,id), cyv2(51,id,k)
! local
integer i, idp1, j, jmd, kmd, lp1, lsp1, ml, msp1
real(dp) qq(id,k-id), cyy(id,id), z(id,k-id), zz(id,k-id),&
&cst0
!
cst0=0.0d-00
lp1=l+1
idp1 = id + 1
kmd = k - id
do 200 msp1=1,lp1
!xx      DO 209 I=1,ID
!c	DO 209 J=1,K
!xx      DO 209 J=1,KMD
!xx  209 Z(I,J)=CST0
z(1:id,1:kmd)=cst0
do 300 lsp1=1,lp1
!xx      DO 210 I=1,ID
do 211 i=1,id
do 210 j=1,kmd
!xx  210 QQ(I,J)=Q(LSP1,J,I)
qq(i,j)=q(lsp1,j,i)
210 continue
211 continue
ml=lsp1-msp1
if(ml.lt.0) go to 319
!xx      DO 320 I=1,ID
do 321 i=1,id
do 320 j=1,id
!xx  320 CYY(I,J)=CYY1(ML+1,I,J)
cyy(i,j)=cyy1(ml+1,i,j)
320 continue
321 continue
go to 340
319 ml=-ml
!xx      DO 330 I=1,ID
do 331 i=1,id
do 330 j=1,id
!xx  330 CYY(I,J)=CYY1(ML+1,J,I)
cyy(i,j)=cyy1(ml+1,j,i)
330 continue
331 continue
!c  340 CALL  MULPLY ( CYY,QQ,ZZ,ID,ID,KMD,MJ2,MJ2,MJ1 )
!c	CALL  MATADL ( Z,ZZ,ID,KMD,MJ2,MJ1 )
340 call  mulply ( cyy,qq,zz,id,id,kmd )
call  matadl ( z,zz,id,kmd )
300 continue
!xx      DO 350 I=1,ID
do 351 i=1,id
do 349 j=1,id
!xx  349 CYV2(MSP1,I,J)=CYY1(MSP1,J,I)
cyv2(msp1,i,j)=cyy1(msp1,j,i)
349 continue
do 350 j=idp1,k
jmd = j - id
!xx  350 CYV2(MSP1,I,J)=Z(I,JMD)
cyv2(msp1,i,j)=z(i,jmd)
350 continue
351 continue
200 continue
return
end
!
!c	SUBROUTINE HCXV2(M,L,K,ID,MJ1,MJ2)
subroutine hcxv2(m,l,k,id,q,cxy1,cxy2,cxv2)
  use timsac_kinds, only: dp
  implicit none
!     THIS SUBROUTINE COMPUTES CXV2(MS) MS=0,M-1.
!     CXV2(MS)=CXV1(-MS)
!     =CXY1(-MS)*Q(0)'+CXY1(-MS+1)*Q(1)'+...+CXY1(-MS+L)*Q(L)'
!xx      IMPLICIT REAL*8(A-H,O-Z)
!c      COMMON /COM2/Q
!c      COMMON /COM11/CXY1
!c      COMMON /COM110/CXY2
!c      COMMON /COM19/CXV2
!c      DIMENSION Q(50,10,5),CXY1(100,5,5),CXY2(50,5,5),CXV2(51,5,10)
!c	DIMENSION QQ(5,10),CXY(5,5),Z(5,10),ZZ(5,10)
!xx      DIMENSION Q(50,K,ID),CXY1(100,ID,ID),CXY2(50,ID,ID),CXV2(51,ID,K)
!xx      DIMENSION QQ(ID,K-ID),CXY(ID,ID),Z(ID,K-ID),ZZ(ID,K-ID)
integer m, l, k, id
real(dp) q(50,k,id), cxy1(100,id,id), cxy2(50,id,id),&
&cxv2(51,id,k)
! local
integer i, idp1, j, jmd, kmd, lp1, lsp1, ml, msp1
real(dp) qq(id,k-id), cxy(id,id), z(id,k-id), zz(id,k-id),&
&cst0
!
cst0=0.0d-00
lp1=l+1
idp1 = id + 1
kmd = k - id
do 200 msp1=1,m
!xx      DO 209 I=1,ID
!c	DO 209 J=1,K
!xx      DO 209 J=1,KMD
!xx  209 Z(I,J)=CST0
z(1:id,1:kmd)=cst0
do 300 lsp1=1,lp1
!xx      DO 210 I=1,ID
do 211 i=1,id
do 210 j=1,kmd
!xx  210 QQ(I,J)=Q(LSP1,J,I)
qq(i,j)=q(lsp1,j,i)
210 continue
211 continue
ml=lsp1-msp1
if(ml.lt.0) go to 319
!xx      DO 320 I=1,ID
do 321 i=1,id
do 320 j=1,id
!xx  320 CXY(I,J)=CXY1(ML+1,I,J)
cxy(i,j)=cxy1(ml+1,i,j)
320 continue
321 continue
go to 340
319 ml=-ml
!xx      DO 330 I=1,ID
do 331 i=1,id
do 330 j=1,id
!xx  330 CXY(I,J)=CXY2(ML+1,I,J)
cxy(i,j)=cxy2(ml+1,i,j)
330 continue
331 continue
!c  340 CALL MULPLY ( CXY,QQ,ZZ,ID,ID,KMD,MJ2,MJ2,MJ1 )
!c	CALL  MATADL ( Z,ZZ,ID,KMD,MJ2,MJ1 )
340 call mulply ( cxy,qq,zz,id,id,kmd )
call  matadl ( z,zz,id,kmd )
300 continue
!xx      DO 350 I=1,ID
do 351 i=1,id
do 349 j=1,id
!xx  349 CXV2(MSP1,I,J)=CXY2(MSP1,I,J)
cxv2(msp1,i,j)=cxy2(msp1,i,j)
349 continue
do 350 j=idp1,k
jmd = j - id
!xx  350 CXV2(MSP1,I,J) = Z(I,JMD)
cxv2(msp1,i,j) = z(i,jmd)
350 continue
351 continue
200 continue
return
end
!
!c	SUBROUTINE SBCXY2(M,L,ID,MJ2)
subroutine sbcxy2(cyy1,mj3,m,l,id,x2,cxy2)
  use timsac_kinds, only: dp
  implicit none
!     THIS SUBROUTINE COMPUTES CXY2(LS) LS=0,M-1.
!     CXY2(LS)=CXY1(-LS)
!     =X2(0)*CYY1(-LS)+X2(1)*CYY1(-LS-1)+X2(2)*CYY1(-LS-2)+...
!     +X2(L)*CYY1(-LS-L).
!xx      IMPLICIT REAL*8(A-H,O-Z)
!c      COMMON /COM3/X2
!c      COMMON /COM10/CYY1
!c      COMMON /COM110/CXY2
!
!c	DIMENSION X2(50,5,5),CYY1(160,5,5),CXY2(50,5,5)
!c	DIMENSION XX(5,5),CYY(5,5),Z(5,5),ZZ(5,5)
!xx      DIMENSION X2(50,ID,ID),CYY1(MJ3,ID,ID),CXY2(50,ID,ID)
!xx      DIMENSION XX(ID,ID),CYY(ID,ID),Z(ID,ID),ZZ(ID,ID)
integer mj3, m, l, id
real(dp) cyy1(mj3,id,id), x2(50,id,id), cxy2(50,id,id)
! local
integer i, j, lp1, lsp1, ml, msp1
real(dp) xx(id,id), cyy(id,id), z(id,id), zz(id,id), cst0
!
cst0=0.0d-00
lp1=l+1
do 200 lsp1=1,m
!xx      DO 209 I=1,ID
!xx      DO 209 J=1,ID
!xx  209 Z(I,J)=CST0
z(1:id,1:id)=cst0
do 300 msp1=1,lp1
!xx      DO 310 I=1,ID
do 311 i=1,id
do 310 j=1,id
!xx  310 XX(I,J)=X2(MSP1,I,J)
xx(i,j)=x2(msp1,i,j)
310 continue
311 continue
ml=lsp1+msp1-1
!xx      DO 320 I=1,ID
do 321 i=1,id
do 320 j=1,id
!xx  320 CYY(I,J)=CYY1(ML,J,I)
cyy(i,j)=cyy1(ml,j,i)
320 continue
321 continue
!c	CALL MULPLY(XX,CYY,ZZ,ID,ID,ID,MJ2,MJ2,MJ2)
!c	CALL MATADL(Z,ZZ,ID,ID,MJ2,MJ2)
call mulply(xx,cyy,zz,id,id,id)
call matadl(z,zz,id,id)
300 continue
!xx      DO 350 I=1,ID
do 351 i=1,id
do 350 j=1,id
!xx  350 CXY2(LSP1,I,J)=Z(I,J)
cxy2(lsp1,i,j)=z(i,j)
350 continue
351 continue
200 continue
return
end
!
!c	SUBROUTINE SUBDAV(X,C0,G,R,N,K,ID,IAW,IPQ,MJ1,MJ2,ISWRO)
subroutine subdav(cyy1,mj3,b,aw,iaw,nh,idd,ir,ij,ik,vd,x,c0,g,r,&
!x     *  N,K,ID,IPQ,AICD,MJ1,MJ4,ISWRO,O,Q,X1,X2,CXY1,CXX1,CXV1,IFG,LU)
!xx     *  N,K,ID,IPQ,AICD,MJ1,MJ4,ISWRO,O,Q,X1,X2,CXY1,CXX1,CXV1)
&n,k,id,ipq,aicd,mj4,iswro,o,q,x1,x2,cxy1,cxx1,cxv1)
  use timsac_kinds, only: dp
  implicit none
!     DAVIDON'S (MINIMIZATION) PROCEDURE
!xx      IMPLICIT REAL*8(A-H,O-Z)
!c	COMMON /COM50/VD
real(dp) aico
common /com101/aico
!c	COMMON /COM199/AICD
!c	DIMENSION CYY1(160,5,5),B(10,5)
!c	DIMENSION IJ(5),IK(5)
!c	DIMENSION VD(50,50)
!c	DIMENSION X(100),C0(5,5),G(100),R(100),C0D(5,5)
!c     DIMENSION SX(100),SC0(5,5),SC0D(5,5),SG(100),SR(100)
!c	DIMENSION AW(50),IDD(10),IR(10)
!xx      DIMENSION CYY1(MJ3,ID,ID),B(K,ID)
!xx      DIMENSION IJ(ID),IK(ID)
!xx      DIMENSION VD(MJ4,MJ4)
!xx      DIMENSION X(MJ4),C0(ID,ID),G(MJ4),R(MJ4),C0D(ID,ID)
!xx      DIMENSION SX(MJ4),SC0(ID,ID),SC0D(ID,ID),SG(MJ4),SR(MJ4)
!c      DIMENSION AW(IAW),IDD(K),IR(K)
!xx      DIMENSION AW(IAW),NH(K),IDD(K),IR(K)
!xx      DIMENSION O(50,ID,K),Q(50,K,ID),X1(50,K,K),X2(50,ID,ID)
!xx      DIMENSION CXY1(100,ID,ID),CXX1(51,ID,ID),CXV1(51,ID,K)
integer mj3, iaw, n, k, id, ipq, mj4, iswro, nh(k), idd(k),&
&ir(k), ij(id), ik(id)
real(dp) cyy1(mj3,id,id), b(k,id), aw(iaw), vd(mj4,mj4),&
&x(mj4), c0(id,id), g(mj4), r(mj4), aicd,&
&o(50,id,k), q(50,k,id), x1(50,k,k), x2(50,id,id),&
&cxy1(100,id,id), cxx1(51,id,id), cxv1(51,id,k)
! local
integer i, ig, ipq2, iram, itn, itns, isphai, iphai, j, l, m
real(dp) c0d(id,id), sx(mj4), sc0(id,id), sc0d(id,id),&
&sg(mj4), sr(mj4), cst0, cst1, cst2, cst01, cst05,&
&consta, constb, eps1, eps3, eps4, eps5, aipq,&
&an, phaid, ro, phai, ephai1, t1, ram, ramro,&
&ramrot, sum, sro, srod, dgam, gsr, dgam1, ramsro,&
&ramt, sphai, ram1, consdr, ophai, aics, daic
!
!
!     DAVIDON'S PROCEDURE
!
!     CONSTANT
cst1=1.0d-00
cst0=0.0d-00
cst2=2.0d-00
cst01=0.1d-00
cst05=0.5d-00
consta=0.5d-00
constb=2.0d-00
eps1=0.01d-00
eps3=0.000001d-00
eps4=0.1d-10
eps5=1.0d-15
itn=1
isphai=0
iphai=1
aipq=ipq
an=n
phaid=cst1
!
!c  150 WRITE(6,160) ITN
!x  150 IF (IFG.NE.0) WRITE(LU,160) ITN
!xx  160 FORMAT(//1H ,'ITN=',I5)
150 continue
!
!     RO=G'*R COMPUTATION
itns=0
40 call innerp(g,r,ro,ipq)
!c	WRITE(6,2210) RO
!x      IF (IFG.NE.0) WRITE(LU,2210) RO
if(iphai.eq.0) go to 101
!     DETERMINANT OF CXX(0) COMPUTATION
!xx      DO 100 I=1,ID
do 102 i=1,id
do 100 j=1,id
!xx  100 C0D(I,J)=C0(I,J)
c0d(i,j)=c0(i,j)
100 continue
102 continue
!c	CALL SUBDET(C0D,PHAI,ID,MJ2)
call subdetm(c0d,phai,id)
phaid=phai
101 ophai=phai
ephai1=eps1*phai
t1=ro-cst2*phai
if(t1.le.ephai1) go to 140
!
ram=cst2*phai/ro
!
!     V=V+((RAM-1.0)/RO)*(R*R')
ramro=(ram-cst1)/ro
!xx      DO 110 I=1,IPQ
do 111 i=1,ipq
ramrot=ramro*r(i)
do 110 j=1,ipq
!xx  110 VD(I,J)=VD(I,J)+RAMROT*R(J)
vd(i,j)=vd(i,j)+ramrot*r(j)
110 continue
111 continue
!
!     R=RAM*R
do 120 i=1,ipq
!xx  120 R(I)=RAM*R(I)
r(i)=ram*r(i)
120 continue
if(itns.ge.10) go to 140
itns=itns+1
go to 40
!
!     SX=X-R
!c  140 WRITE(6,600) ITNS,T1,EPHAI1
!x  140 IF (IFG.NE.0) WRITE(LU,600) ITNS,T1,EPHAI1
!xx  600 FORMAT(1H ,'ITNS=',I5,5X,'RO-2.0*PHAI=',D12.5,5X,'EPHAI1=',D12.5)
140 continue
ig=0
205 do 210 i=1,ipq
!xx  210 SX(I)=X(I)-R(I)
sx(i)=x(i)-r(i)
210 continue
!
!     SPHAI, SG COMPUTATION
!c	CALL C0GR(SX,SC0,SG,M,L,K,ID,IAW,IPQ,MJ1,MJ2,IG)
call c0gr(cyy1,mj3,b,aw,iaw,nh,idd,ir,ij,ik,sx,sc0,sg,m,l,k,id,&
!x     *        IPQ,MJ4,IG,O,Q,X1,X2,CXY1,CXX1,CXV1,IFG,LU)
!xx     *        IPQ,MJ4,IG,O,Q,X1,X2,CXY1,CXX1,CXV1)
&mj4,ig,o,q,x1,x2,cxy1,cxx1,cxv1)
if(ig.eq.0) go to 215
!c	WRITE(6,211)
!x      IF (IFG.NE.0) WRITE(LU,211)
!xx  211 FORMAT(1H ,'HESSIAN MODIFIED FOR FEASIBILITY')
!xx      DO 213 I=1,IPQ
do 214 i=1,ipq
r(i)=cst05*r(i)
do 213 j=1,ipq
!xx  213 VD(I,J)=CST05*VD(I,J)
vd(i,j)=cst05*vd(i,j)
213 continue
214 continue
go to 205
!
!     DETERMINANT OF (NEW CXX(0)) COMPUTATION
!xx  215 DO 220 I=1,ID
215 do 221 i=1,id
do 220 j=1,id
!xx  220 SC0D(I,J)=SC0(I,J)
sc0d(i,j)=sc0(i,j)
220 continue
221 continue
!c	CALL SUBDET(SC0D,SPHAI,ID,MJ2)
call subdetm(sc0d,sphai,id)
!
!
!     SR=V*SG
do 310 i=1,ipq
sum=cst0
do 311 j=1,ipq
!xx  311 SUM=SUM+VD(I,J)*SG(J)
sum=sum+vd(i,j)*sg(j)
311 continue
!xx  310 SR(I)=SUM
sr(i)=sum
310 continue
!
!     SRO=(SG)'*(SR)
call innerp(sg,sr,sro,ipq)
srod=sro/phai
!c	WRITE(6,340) SRO
!x      IF (IFG.NE.0) WRITE(LU,340) SRO
!c  340 FORMAT(1H ,'SRO=',D17.5)
!xx  340 FORMAT(/1H ,'SRO=',D17.5)
!
!
!     DGAM=-G'*(SR)/SRO
call innerp(g,sr,gsr,ipq)
dgam=-gsr/sro
dgam1=dgam+cst1
dgam1=dabs(dgam1)+0.1d-70
ram=dabs(dgam)/dgam1
!     IF(RAM.LE.CONSTA) RAM=CONSTA
if(ram.gt.consta) go to 430
ram=consta
iram=1
!c	WRITE(6,420) RAM
!x      IF (IFG.NE.0) WRITE(LU,420) RAM
!xx  420 FORMAT(1H ,'RAM=CONSTA=',D17.5)
go to 470
!     IF(RAM.GE.CONSTB) RAM=CONSTB
430 if(ram.lt.constb) go to 450
ram=constb
iram=-1
!c	WRITE(6,440) RAM
!x      IF (IFG.NE.0) WRITE(LU,440) RAM
!xx  440 FORMAT(1H ,'RAM=CONSTB=',D17.5)
go to 470
!     RAM=RAM
450 continue
iram=0
!c	WRITE(6,460) RAM
!x      IF (IFG.NE.0) WRITE(LU,460) RAM
!xx  460 FORMAT(1H ,'RAM=DGAM/(DGAM+1.0)=',D17.5)
!
!     V=V+((RAM-1.0)/SRO)*(SR)*(SR)'
470 ramsro=(ram-cst1)/sro
!xx      DO 480 I=1,IPQ
do 481 i=1,ipq
ramt=ramsro*sr(i)
do 480 j=1,ipq
!xx  480 VD(I,J)=VD(I,J)+RAMT*SR(J)
vd(i,j)=vd(i,j)+ramt*sr(j)
480 continue
481 continue
!
!
!
if(phai.gt.sphai) go to 540
!     SPHAI.GE.PHAI: TEST OF CORRECTION
ram1=ram-cst1
if(dabs(ram1).lt.eps3) go to 555
consdr=dgam*ram1
do 550 i=1,ipq
!xx  550 R(I)=R(I)-CONSDR*SR(I)
r(i)=r(i)-consdr*sr(i)
550 continue
iphai=0
if(srod.gt.eps4) go to 900
!     END OF ITERATION
555 iswro=iswro+1
go to 1000
!
!     SPHAI.LT.PHAI: SUCCESSFUL REDUCTION
540 do 560 i=1,ipq
x(i)=sx(i)
g(i)=sg(i)
!xx  560 R(I)=RAM*SR(I)
r(i)=ram*sr(i)
560 continue
!xx      DO 570 I=1,ID
do 571 i=1,id
do 570 j=1,id
!xx  570 C0(I,J)=SC0(I,J)
c0(i,j)=sc0(i,j)
570 continue
571 continue
phai=sphai
iphai=1
!xx  800 CONTINUE
aico=an*dlog(ophai)+cst2*aipq
aics=an*dlog(sphai)+cst2*aipq
phaid=sphai
aicd=aics
daic=aico-aics
!c	WRITE(6,810) OPHAI,SPHAI,AICO,AICS,DAIC
!c  810 FORMAT(1H ,'OPHAI=',D17.5,5X,'SPHAI=',D17.5,5X,'OAIC=',D17.5,5X,'A
!c     AIC=',D17.5,5X,'DAIC=',D17.5)
!x      IF (IFG.NE.0) WRITE(LU,810) OPHAI,SPHAI,AICO,AICS,DAIC
!xx  810 FORMAT(1H ,'OPHAI=',D17.5,5X,'SPHAI=',D17.5,/,' OAIC=',D17.5,5X,'A
!xx     AIC=',D17.5,5X,'DAIC=',D17.5)
if(iram.ne.0) go to 901
if(srod.lt.eps4) go to 555
!     ITERATION CHECK
900 ipq2=ipq+ipq
if(itn.ge.ipq2) go to 555
isphai=(isphai+(1-iphai))*(1-iphai)
if(isphai.ge.5) go to 555
itn=itn+1
go to 150
901 if(daic.ge.cst01) go to 900
if(srod.lt.eps4) go to 555
go to 900
!     END OF MINIMIZATION
!xx  999 ISWRO=0
iswro=0
!
1000 continue
aicd=an*dlog(phaid)+cst2*aipq
!c	WRITE(6,3609) AICD
!x      IF (IFG.NE.0) WRITE(LU,3609) AICD
!xx 3609 FORMAT(/1H ,'AIC=N*LOG(DET(CXX0))+2.0*IPQ=',D12.5)
!     INNOVATION VARIANCE PRINT OUT
!c	WRITE(6,3610)
!c	DO 3620 I=1,ID
!c 3620 WRITE(6,3630) I,(C0(I,J),J=1,ID)
!x      IF (IFG.NE.0) THEN
!x	 WRITE(LU,3610)
!x	 DO 3620 I=1,ID
!x 3620	 WRITE(LU,3630) I,(C0(I,J),J=1,ID)
!x      END IF
return
!xx 3610 FORMAT(/1H ,'INNOVATION VARIANCE')
!xx 3630 FORMAT(/1H ,I5,10D12.5,/(1H ,5X,10D12.5))
!xx 2210 FORMAT(1H ,'RO=',D12.5)
end
!
!c	SUBROUTINE ARMACO(ID,K,MJ1,MJ2)
!xx      SUBROUTINE ARMACO(A,B,AW,IAW,IH,IDD,IR,IJ,IK,ID,ICONT,K,BBM,
!xx     *                  AAU,ZZ,IQM,MJ6,MJ7)
subroutine armaco(b,aw,iaw,ih,idd,ir,ik,id,icont,k,bbm,&
&aau,zz,iqm,mj6,mj7)
  use timsac_kinds, only: dp
  implicit none
!     AR-MA COEFFICIENTS COMPUTATION
!xx      IMPLICIT REAL*8 (A-H,O-Z)
!c	COMMON /COM60/ IDD,IR
!c	COMMON /COM61/ IJ,IK
!c	COMMON /COM70/ IH
!c	COMMON /COM80/ A
!c	COMMON /COM81/ B
!c	COMMON /COM82/ AW
!c	COMMON /COM87/ICONT
!c	DIMENSION A(10,10),B(10,5)
!c	DIMENSION AW(50)
!c	DIMENSION IH(10),IDD(10),IR(10),IJ(5),IK(5)
!xx      DIMENSION A(K,K),B(K,ID)
!xx      DIMENSION AW(IAW)
!xx      DIMENSION IH(K),IDD(K),IR(K),IJ(ID),IK(ID)
!xx      DIMENSION IO(K)
!c	DIMENSION BA(5,55),B0(5,5),BM(5,55),U(10,5),AU(10,5)
!c	DIMENSION W0(5,5),Z(5,5),ZZ(5,5)
!xx      DIMENSION BA(ID,MJ6),B0(ID,ID),BM(ID,MJ6)
!xx      DIMENSION U(K,ID),AU(K,ID,MJ7),AAU(ID,ID,MJ7)
!xx      DIMENSION BBM(ID,ID,MJ6),W0(ID,ID),Z(ID,ID),ZZ(ID,ID,MJ7)
integer iaw, id, icont, k, iqm, mj6, mj7, ih(k), idd(k), ir(k),&
&ik(id)
real(dp) b(k,id), aw(iaw), bbm(id,id,mj6),&
&aau(id,id,mj7), zz(id,id,mj7)
! local
integer i, idl, idl1, iki, imaj, iq, iqd, iqm1, iqmax, ird, irs,&
&isr, isrd, isum, it, itm1, j, ja, jb, jb1, jbj, jdq,&
&jma, jo, jt, jtb, jtbj, jtw, jtwj, jw, jwj, ns, io(k)
real(dp) ba(id,mj6), b0(id,id), bdet, bm(id,mj6), u(k,id),&
&au(k,id,mj7), w0(id,id), z(id,id), cst0, cst1
!     BA COMPUTATION
!     ORDER DETERMINATION
cst0=0.0d-00
cst1=1.0d-00
iqmax=0
do 100 i=1,k
iq=(ih(i)-1)/id+1
io(i)=iq
if(io(i).le.iqmax)      go to 100
iqmax=iq
100 continue
!
!     TRANSFORMATION OF AW(F MATRIX IN VECTOR FORM) INTO THE MATRICES OF
!     THE RAW AR-COEFFICIENTS.
!
!c	DO 150 I=1,5
!c	DO 150 J=1,55
!xx      DO 150 I=1,ID
!xx      DO 150 J=1,MJ6
!xx  150 BA(I,J)=CST0
ba(1:id,1:mj6)=cst0
isum=0
irs=0
do 201 i=1,k
if(idd(i).eq.0) go to 201
irs=irs+1
iki=ik(irs)
jo=ih(i)-(io(i)-1)*id
idl=(iqmax-io(i))*id
idl1=idl+id
ja=ih(i)+idl1
ba(jo,ja)=-cst1
!
do 200 j=1,iki
isum=isum+1
jb=ih(j)+idl
ba(jo,jb)=aw(isum)
200 continue
201 continue
!     AR-COEFFICIENTS OUTPUT
!c	WRITE(6,7000)
!c 7000 FORMAT(//1H ,'AR-COEFFICIENTS B(MS) (MS=1,Q)')
iqmax=iqmax+1
iqm=iqmax-1
iqm1=iqm-1
!c	WRITE(6,7200) IQM
!c7200	FORMAT(/1H ,'K=',I5)
!c	WRITE(7,1000) ID,IQM,IQM1
iqd=iqm*id
!xx      DO 300 I=1,ID
do 301 i=1,id
do 300 j=1,id
jdq=iqd+j
!xx  300 B0(I,J)=-BA(I,JDQ)
b0(i,j)=-ba(i,jdq)
300 continue
301 continue
!c	MJ5=5
!c	MJ55=55
!c	CALL INVDET(B0,BDET,ID,MJ5)
!c	CALL MULPLY(B0,BA,BM,ID,ID,IQD,MJ5,MJ5,MJ55)
call invdet(b0,bdet,id,id)
call mulply(b0,ba,bm,id,id,iqd)
!c	DO 321 I=1,5
!c	DO 321 J=1,55
!xx      DO 321 I=1,ID
do 322 i=1,id
do 321 j=1,iqd
!xx  321 BM(I,J)=-BM(I,J)
bm(i,j)=-bm(i,j)
321 continue
322 continue
do 320 irs=1,iqm
!c	WRITE(6,7250) IRS
isr=iqm-irs
isrd=isr*id+1
ird=isrd+id-1
!c	DO 310 I=1,ID
!c	WRITE(6,1001) I,(BM(I,ISJ),ISJ=ISRD,IRD)
!c	WRITE(7,1002) (BM(I,ISJ),ISJ=ISRD,IRD)
!c  310 CONTINUE
!xx      DO 310 I=1,ID
do 311 i=1,id
do 310 j=1,id
bbm(i,j,irs) = bm(i,isrd+j-1)
310 continue
311 continue
320 continue
!     IMPULSE RESPONSE MATRICES OUTPUT
!c	WRITE(6,8000)
!c8000	FORMAT(/1H ,'IMPULSE RESPONSE MATRICES W(LS) (LS=1,Q-1)')
!xx      DO 400 I=1,K
do 401 i=1,k
do 400 j=1,id
!xx400   U(I,J)=B(I,J)
u(i,j)=b(i,j)
400 continue
401 continue
iqm1=iqm-1
imaj=-id
do 490 ns=1,iqm1
!c	WRITE(6,8100) NS
imaj=imaj+id
!c	CALL SUBAWZ(U,AU,K,ID,MJ1,MJ2)
call subawz(aw,iaw,idd,ir,u,au(1,1,ns),k,id)
do 470 i=1,id
if(icont.eq.0) go to 440
!c	WRITE(6,1001) I,(AU(I,J),J=1,ID)
!c  450 WRITE(7,1002) (AU(I,J),J=1,ID)
do 450 j=1,id
aau(i,j,ns)=au(i,j,ns)
450 continue
440 do 460 j=1,id
jma=imaj+j
!c  460 BA(I,JMA)=AU(I,J)
!xx  460 BA(I,JMA)=AU(I,J,NS)
ba(i,jma)=au(i,j,ns)
460 continue
470 continue
!xx      DO 480 I=1,K
do 481 i=1,k
do 480 j=1,id
!c480	U(I,J)=AU(I,J)
!xx480   U(I,J)=AU(I,J,NS)
u(i,j)=au(i,j,ns)
480 continue
481 continue
490 continue
if(icont.eq.1) go to 506
!     MA-COEFFICIENT COMPUTATION
!c	WRITE(6,8200)
!c 8200 FORMAT(/1H ,'MA-COEFFICIENTS A(LS) (LS=1,Q-1)')
jw=-id
jb=iqd
jb1=iqd
do 505 it=1,iqm1
jw=jw+id
jb=jb-id
!xx      DO 500 I=1,ID
do 510 i=1,id
do 500 j=1,id
jwj=jw+j
jbj=jb+j
!c  500 ZZ(I,J)=BA(I,JWJ)+BM(I,JBJ)
!xx  500 ZZ(I,J,IT)=BA(I,JWJ)+BM(I,JBJ)
zz(i,j,it)=ba(i,jwj)+bm(i,jbj)
500 continue
510 continue
if(it.eq.1) go to 503
itm1=it-1
jtw=jw
jtb=jb1
do 502 jt=1,itm1
jtw=jtw-id
jtb=jtb-id
!xx      DO 501 I=1,ID
do 511 i=1,id
do 501 j=1,id
jtbj=jtb+j
jtwj=jtw+j
b0(i,j)=bm(i,jtbj)
!xx  501 W0(I,J)=BA(I,JTWJ)
w0(i,j)=ba(i,jtwj)
501 continue
511 continue
!c	CALL MULPLY(B0,W0,Z,ID,ID,ID,MJ2,MJ2,MJ2)
!c	CALL MATADL(ZZ,Z,ID,ID,MJ2,MJ2)
call mulply(b0,w0,z,id,id,id)
call matadl(zz(1,1,it),z,id,id)
502 continue
!c  503 WRITE(6,8300) IT
503 continue
!c	DO 504 I=1,ID
!c	WRITE(6,1001) I,(ZZ(I,J),J=1,ID)
!c	WRITE(7,1002) (ZZ(I,J),J=1,ID)
!c  504 CONTINUE
505 continue
506 continue
return
!xx1000  FORMAT(16I5)
!xx 1001 FORMAT(1H ,I5,4X,10D12.5)
!xx1002  FORMAT(4D20.10)
!xx1003  FORMAT(6I5)
!xx7250  FORMAT(/1H ,'B(MS) MS=',I5)
!xx2003  FORMAT(/1H ,'BA')
!xx8100  FORMAT(/1H ,'W(LS) LS=',I5)
!xx 8300 FORMAT(/1H ,'A(LS) LS=',I5)
end
!
!c	SUBROUTINE  RESCAL ( K,ID )
!x      SUBROUTINE  RESCAL ( B,AW,IAW,NH,IDD,IR,IJ,IDUMMY,K,ID,IFG,LU )
!xx      SUBROUTINE  RESCAL ( B,AW,IAW,NH,IDD,IR,IJ,IDUMMY,K,ID )
subroutine  rescal ( b,aw,iaw,nh,idd,ir,k,id )
  use timsac_kinds, only: dp
  implicit none
!     RESCALING OF F AND G MATRICES FOR A FEASIBLE INITIAL
!xx      IMPLICIT REAL*8(A-H,O-Z)
!c	COMMON /COM60/IDD , IR
!c	COMMON /COM61/IJ, IDUMMY
!c	COMMON /COM70/NH
!c	COMMON /COM82/AW
!c	COMMON /COM81/B
!c	DIMENSION AW(50),B(10,5)
!c	DIMENSION IJ(5),NH(10),IR(10)
!c	DIMENSION IDUMMY(5)
!c	DIMENSION IDD(10)
!c	DIMENSION IO(10)
!xx      DIMENSION AW(IAW),B(K,ID)
!xx      DIMENSION IJ(ID),NH(K),IR(K)
!c      DIMENSION IDUMMY(ID)
!xx      DIMENSION IDD(K)
!xx      DIMENSION IO(K)
integer iaw, k, id, nh(k), idd(k), ir(k)
real(dp) b(k,id), aw(iaw)
! local
integer i, idp1, iksj, im, ioi, ior, j, jl, ks, io(k)
real(dp) ro9, scale
!
ro9=0.95d-00
do 10 i=1,k
!xx   10 IO(I)=(NH(I)-1)/ID
io(i)=(nh(i)-1)/id
10 continue
!
!     F MATRIX RESCALING
iksj = 0
do 20 i=1,k
if(idd(i).eq.0) go to 20
im = io(i) + 1
ks=ks+1
jl = ir(i)
do 30 j=1,jl
ior = im - io(j)
scale=ro9**ior
iksj = iksj + 1
!xx   30 AW(IKSJ)=AW(IKSJ)*SCALE
aw(iksj)=aw(iksj)*scale
30 continue
20 continue
!
!     G MATRIX RESCALING
idp1=id+1
do 40 i=idp1,k
ioi=io(i)
scale=ro9**ioi
do 50 j=1,id
!xx   50 B(I,J) = B(I,J) * SCALE
b(i,j) = b(i,j) * scale
50 continue
40 continue
!c	WRITE(6,2000)
!x      IF (IFG.NE.0) WRITE(LU,2000)
!xx 2000 FORMAT(/1H ,'INITIAL MODIFIED FOR FEASIBILITY')
return
end
!
!c	SUBROUTINE SUBDET(X,XDETMI,MM,MJ)
subroutine subdetm(x,xdetmi,mm)
  use timsac_kinds, only: dp
  implicit none
! *** COMMON SUBROUTINE
!     THIS SUBROUTINE COMPUTES THE DETERMINANT OF UPPER LEFT MM X MM
!     OF X.  FOR GENERAL USE STATEMENTS 20-21 SHOULD BE RESTORED.
!     X: ORIGINAL MATRIX
!     XDETMI: DETERMINANT OF UPPER LEFT MM X MM OF X
!     MJ: ABSOLUTE DIMENSION OF X IN THE MAIN ROUTINE
!xx      IMPLICIT REAL*8(X)
!c	DIMENSION X(MJ,MJ)
!xx      DIMENSION X(MM,MM)
integer mm
real(dp) x(mm,mm), xdetmi
! local
integer i, i1, j, jj, k, mm1
real(dp) cst0, cst1, xxc, xc
cst0=0.0d-00
cst1=1.0d-00
xdetmi=cst1
if(mm.eq.1) go to 18
mm1=mm-1
do 10 i=1,mm1
!xx   20 IF(X(I,I).NE.CST0) GO TO 11
if(x(i,i).ne.cst0) go to 11
do 12 j=i,mm
if(x(i,j).eq.cst0) go to 12
jj=j
go to 13
12 continue
xdetmi=cst0
go to 17
13 do 14 k=i,mm
xxc=x(k,jj)
x(k,jj)=x(k,i)
!xx   14 X(K,I)=XXC
x(k,i)=xxc
14 continue
!xx   21 XDETMI=-XDETMI
xdetmi=-xdetmi
11 xdetmi=xdetmi*x(i,i)
xc=cst1/x(i,i)
i1=i+1
do 15 j=i1,mm
xxc=x(j,i)*xc
do 16 k=i1,mm
!xx   16 X(J,K)=X(J,K)-X(I,K)*XXC
x(j,k)=x(j,k)-x(i,k)*xxc
16 continue
15 continue
10 continue
18 xdetmi=xdetmi*x(mm,mm)
17 return
end
