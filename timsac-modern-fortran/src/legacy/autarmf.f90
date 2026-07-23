! Part of the experimental modern Fortran port of timsac 1.3.8-6.
! This file was created or modified for the Fortran project on 2026-07-23.
! Original TIMSAC credits are retained; see NOTICE and ORIGIN.md.
! SPDX-License-Identifier: GPL-2.0-or-later

subroutine autarm(n,lagh01,cyy1,newl1,iqi1,b1,ipi1,a1,&
!x     *	          NEWN,IQ,B2,IP,A2,STD,CXX2,G,SAIC,AICM,KQ,KP,TMP,IER,
&newn,iq,b2,ip,a2,std,cxx2,g,saic,aicm,kq,kp,&
&lmax,mmax,nmax)
  use timsac_kinds, only: dp
  implicit none
!
!c      PROGRAM AUTARM
!     PROGRAM 74.1.2. AUTOMATIC AR-MA MODEL FITTING; SCALAR CASE.
!-----------------------------------------------------------------------
!     ** DESIGNED BY H. AKAIKE, THE INSTITUTE OF STATISTICAL MATHEMATICS
!     ** PROGRAMMED BY E. ARAHATA, THE INSTITUTE OF STATISTICAL MATHEMAT
!        TOKYO
!     ** DATE OF THE LATEST REVISION: MARCH 25, 1977
!     ** THIS PROGRAM WAS ORIGINALLY PUBLISHED IN
!        "TIMSAC-74 A TIME SERIES ANALYSIS AND CONTROL PROGRAM PACKAGE(1
!        BY H. AKAIKE, E. ARAHATA AND T. OZAKI, COMPUTER SCIENCE MONOGRA
!        NO.5, MARCH 1975, THE INSTITUTE OF STATISTICAL MATHEMATICS
!     ** FOR THE BASIC THEORY SEE "CANONICAL CORRELATION ANALYSIS OF TIM
!        AND THE USE OF AN INFORMATION CRITERION" BY H. AKAIKE, IN
!        "SYSTEM IDENTIFICATION: ADVANCES AND CASE STUDIES" R. K. MEHRA
!        D. G. LAINIOTIS EDS. ACADEMIC PRESS, NEW YORK, 1976
!-----------------------------------------------------------------------
!     THIS PROGRAM PROVIDES AN AUTOMATIC AR-MA MODEL FITTING PROCEDURE.
!     MODELS WITH VARIOUS ORDERS ARE FITTED AND THE BEST CHOICE IS DETER
!     WITH THE AID OF THE STATISTICS AIC.
!     THE MAXIMUM LIKELIHOOD ESTIMATES OF THE COEFFICIENTS OF A SCALAR
!     AUTOREGRESSIVE MOVING AVERAGE MODEL Y(I)+B(1)Y(I-1)+...+B(IQ)Y(I-I
!     =X(I)+A(1)X(I-1)+...+A(IP)X(I-IP) OF A TIME SERIES Y(I)
!     ARE OBTAINED BY USING DAVIDON'S VARIANCE ALGORITHM.
!     PURE AUTOREGRESSION IS NOT ALLOWED.
!     FOR AR-MODELS USE THE INTERMEDIATE OUTPUTS OF CANARM.
!
!     THIS PROGRAM REQUIRES THE FOLLWING INPUTS:
!     (N,LAGH): N, LENGTH OF THE ORIGINAL DATA
!                  LAGH, MAXIMUM LAG OF COVARIANCE, NOT GREATER THAN 500
!     CYY(I) (I=0,LAGH): COVARIANCE SEQUENCE ... LAGH SHOULD BE LARGE EN
!                             TO KEEP THE INNOVATION VARIANCE AND
!                             GRADIENT COMPUTATION MEANINGFUL.
!     ***** WHEN N IS NOT GREATER THAN 501, PUT LAGH EQUAL TO N-1.
!     NEWL: TOTAL NUMBER OF CASES, NOT GREATER THAN 25
!     IQ: INITIAL AR ORDER
!     B(I)(I=1,IQ): INITIAL ESTIMATES OF AR-COEFFICIENTS
!     IP: INITIAL MA ORDER
!     A(I)(I=1,IP): INITIAL ESTIMATES OF MA-COEFFICIENTS
!     (IQI(I),IPI(I))(I=1,NEWL-1): (AR,MA) ORDERS TO BE FITTED SUCCESSIV
!                                       UNNECESSARY WHEN NEWL=1
!     ***** WHEN THE BEST CHOICE IS ON THE BORDER, SOME (AR,MA) ORDERS A
!              AUTOMATICALLY WITHIN THE LIMIT OF THE TOTAL NUMBER 25.
!
!     OUTPUTS: FOR EACH PAIR OF AR-MA ORDERS
!     ONE CARD WITH THE STATEMENTS OF THE PROBLEM
!     IQ: AR ORDER
!     B(I)(I=1,IQ): MAXIMUM LIKELIHOOD ESTIMATES OF AR COEFFICIENTS
!     IP: MA ORDER
!     A(I)(I=1,IP): MAXIMUM LIKELIHOOD ESTIMATES OF MA COEFFICIENTS
!     CXX0: INNOVATION VARIANCE
!
!c      PARAMETER (LMAX=500)
!c      PARAMETER (MMAX=50)
!c      PARAMETER (NMAX=25)
!xxx      PARAMETER (ICST=190)
integer, parameter :: icst=190
!
!xx      IMPLICIT REAL*8(A-H,O-Z)
!c      COMMON /COM50/VD
integer iswro, idos, isik
common /com70/iswro
common /com71/idos
common /com72/isik
!c      DIMENSION CYY(1001)
!c      DIMENSION A(50),B(50)
!c      DIMENSION OA(50),OB(50)
!c      DIMENSION X(50)
!c      DIMENSION G(50),STD(50)
!c      DIMENSION C(50)
!c      DIMENSION CC(190)
!c      DIMENSION VD(50,50)
!c      DIMENSION CN(51)
!c      DIMENSION IQI(25),IPI(25),SAIC(25)
!xx      DIMENSION CYY(LMAX*2+1)
!xx      DIMENSION A(MMAX),B(MMAX)
!xx      DIMENSION OA(MMAX),OB(MMAX)
!xx      DIMENSION X(MMAX)
!xx      DIMENSION G(MMAX,NMAX),STD(MMAX,NMAX)
!xx      DIMENSION C(MMAX)
!xx      DIMENSION CC(ICST)
!xx      DIMENSION VD(MMAX,MMAX)
!xx      DIMENSION CN(MMAX+1)
!xx      DIMENSION IQI(NMAX),IPI(NMAX),SAIC(NMAX)
!
!xx      DIMENSION CYY1(LAGH01)
!xx      DIMENSION IQI1(NEWL1),IPI1(NEWL1)
!xx      DIMENSION A1(IPI1(1)),B1(IQI1(1))
!xx      DIMENSION IP(NMAX),IQ(NMAX),IPO(NMAX),IQO(NMAX)
!xx      DIMENSION A2(MMAX,NMAX),B2(MMAX,NMAX)
!xx      DIMENSION SMAIC2(NMAX),CXX2(NMAX)
integer n, lagh01, newl1, newn, kq, kp, lmax, mmax, nmax,&
&iqi1(newl1), ipi1(newl1), iq(nmax), ip(nmax)
real(dp) cyy1(lagh01), a1(ipi1(1)), b1(iqi1(1)),&
&a2(mmax,nmax), b2(mmax,nmax), std(mmax,nmax),&
&cxx2(nmax), g(mmax,nmax), saic(nmax), aicm
! local
integer i, ii, iim, iip, ib, ig, ide, ido, ipq, ipm1, iqm1, isfin,&
&j, jp, jp1, jpo, jq, jq1, jqo, lagh, lagh1, lagh2, lagh4,&
&newl, newlm1,&
&iqi(nmax), ipi(nmax), ipo(nmax), iqo(nmax)
real(dp) cyy(lmax*2+1), a(mmax), b(mmax), oa(mmax),&
&ob(mmax), x(mmax), c(mmax), cc(icst),&
&vd(mmax,mmax), cn(mmax+1), smaic2(nmax), cst0,&
&cst1, cst2, cst05, smaic, an, cxx0, sum, aipq,&
&dmaic, const1, san
!
!x      INTEGER*1  TMP(1)
!x      CHARACTER  CNAME*80
!
!     INPUT / OUTPUT DATA FILE OPEN
!c	CHARACTER(100) DFNAM
!c	CALL SETWND
!c	DFNAM='autarm.out'
!c	CALL FLOPN3(DFNAM,NFL)
!c	IF (NFL.EQ.0) GO TO 999
!
!x      IER=0
!x      LU=3
!x      DO 7 I = 1,80
!x    7 CNAME(I:I) = ' '
!x      I = 1
!x      IFG = 1
!x      DO WHILE( (IFG.EQ.1) .AND. (I.LE.80) )
!x	   IF ( TMP(I).NE.ICHAR(' ') ) THEN
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
!xcx            WRITE(*,*) ' ***  autarm temp FILE OPEN ERROR :',CNAME,IVAR
!x            IER=IVAR
!x            IFG=0
!x         END IF
!x      ELSE
!x         IFG = 0
!x      END IF
!
cst0=0.0d-00
cst1=1.0d-00
cst2=2.0d-00
cst05=0.00005d-00
!c      DO 8 I=1,1001
do 8 i=1,2*lmax+1
!xx    8 CYY(I)=CST0
cyy(i)=cst0
8 continue
!c	LAGH4=501
!c	LAGH1=500
lagh4=lmax+1
lagh1=lmax
!     AUTOCOVARIANCE LOADING
!c	READ(5,1) N,LAGH
lagh=lagh01-1
lagh2=lagh4+lagh
!c	READ(5,2) (CYY(I),I=LAGH4,LAGH2)
do 88 i=1,lagh01
cyy(lmax+i)=cyy1(i)
88 continue
do 9 i=1,lagh1
iim=lagh4-i
iip=lagh4+i
!xx    9 CYY(IIM)=CYY(IIP)
cyy(iim)=cyy(iip)
9 continue
an=n
smaic=an*dlog(cyy(lagh4))
!c      IQO=0
!c      IPO=0
jqo=0
jpo=0
!c 3200 READ(5,1) NEWL
newl=newl1
newn=1
!     INITIAL CONDITION LOADING FOR AR-MA (IQ,IP)
!c 2000 READ(5,1) IQ
jq=iqi1(1)
jp=ipi1(1)
!c      IF(IQ.LE.0) GO TO 4205
if(jq.le.0) go to 4205
!c      READ(5,2) (B(I),I=1,IQ)
do 10 i=1,jq
b(i)=b1(i)
10 continue
!c 4205 READ(5,1) IP
4205 continue
!c      IF(IP.LE.0) GO TO 4204
if(jp.le.0) go to 4204
!c      READ(5,2) (A(I),I=1,IP)
do 20 i = 1,jp
a(i)=a1(i)
20 continue
!
4204 newlm1=newl-1
!c      IQI(1)=IQ
!c      IPI(1)=IP
iqi(1)=jq
ipi(1)=jp
if(newlm1.eq.0) go to 4208
do 4206 i=2,newl
!c 4206 READ(5,1) IQI(I),IPI(I)
iqi(i)=iqi1(i)
!xx 4206 IPI(I)=IPI1(I)
ipi(i)=ipi1(i)
4206 continue
4208 idos=3
isik=0
4207 continue
iswro=0
!     INITIAL PRINT OUT
!c      WRITE (6,11111)
!x      IF (IFG.NE.0) WRITE(LU,11111)
!c	WRITE(6,1600)
!c	WRITE(6,1601)
!c	WRITE(6,1610) N,LAGH,IQ,IP
!c 4210 IPQ=IP+IQ
4210 ipq=jp+jq
aipq=ipq
!c	WRITE(6,3112) ISWRO
!c 3112 FORMAT(/1H ,'ISWRO=',I5)
!c	IF(IQ.LE.0) GO TO 4215
!c	WRITE(6,1622)
!c	DO 1632 I=1,IQ
!c 1632 WRITE(6,1611) I,B(I)
!c 4215 IF(IP.LE.0) GO TO 4216
!c	WRITE(6,1621)
!c	DO 1631 I=1,IP
!c 1631 WRITE(6,1611) I,A(I)
!c 4216 CONTINUE
!x      IF (IFG.NE.0) THEN
!x	 WRITE(LU,3112) ISWRO
!x 3112	 FORMAT(/' ISWRO=',I5)
!x	 IF (JQ.GT.0) THEN
!x	      WRITE(LU,1622)
!x	    DO 1632 I=1,JQ
!x 1632	    WRITE(LU,1611) I,B(I)
!x      END IF
!x 4215	 IF (JP.GT.0) THEN
!x         WRITE(LU,1621)
!x         DO 1631 I=1,JP
!x 1631	    WRITE(LU,1611) I,A(I)
!x      END IF
!x      END IF
!c      DO 100 I=1,IP
do 100 i=1,jp
!xx  100 X(I)=A(I)
x(i)=a(i)
100 continue
!c      IF(IQ.LE.0) GO TO 420
!c	DO 110 I=1,IQ
!c	II=IP+I
if(jq.le.0) go to 420
do 110 i=1,jq
ii=jp+i
!xx  110 X(II)=B(I)
x(ii)=b(i)
110 continue
420 continue
!     INNOVATION VARIANCE, GRADIENT AND HESSIAN COMPUTATION
!c	CALL SC0GRH(X,CYY,G,CN,CXX0,IP,IQ)
!xx      CALL SC0GRH(X,CYY,G(1,NEWN),CN,CXX0,JP,JQ,VD,MMAX,LMAX,ICST,
!xx     *                   IFG,LU)
call sc0grh(x,cyy,g(1,newn),cn,cxx0,jp,jq,vd,mmax,lmax,icst)
!     INVERSE OF HESSIAN COMPUTATION
!     COMMON SUBROUTINE CALL
!c	CALL MATINV(HESDET,IPQ)
!xx      CALL MATINV(IPQ,VD,MMAX,0,LU)
call matinv(ipq,vd,mmax)
!     CORRECTION TERM C(X)=V*G(X) COMPUTATION
do 900 i=1,ipq
sum=cst0
do 910 j=1,ipq
!c  910 SUM=SUM+VD(I,J)*G(J)
!xx  910 SUM=SUM+VD(I,J)*G(J,NEWN)
sum=sum+vd(i,j)*g(j,newn)
910 continue
!xx  900 C(I)=SUM
c(i)=sum
900 continue
!     DAVODON'S PROCEDURE; MINIMIZATION OF INNOVATION VARIANCE
!c	CALL SDAV1(X,CYY,CXX0,G,C,IP,IQ,N)
!xx      CALL SDAV1(X,CYY,CXX0,G(1,NEWN),C,JP,JQ,N,VD,MMAX,LMAX,ICST,
!xx     *                IFG,LU)
call sdav1(x,cyy,cxx0,g(1,newn),c,jp,jq,n,vd,mmax,lmax,icst)
if(iswro.le.0) go to 940
isfin=0
if(iswro.ge.10) go to 940
!c	DO 902 I=1,IP
!c	IF(DABS(A(I)-X(I)).GE.CST05) GO TO 904
do 902 i=1,jp
if(dabs(a(i)-x(i)).ge.cst05) go to 904
902 continue
!c	IF(IQ.LE.0) GO TO 919
!c	DO 903 I=1,IQ
!c	II=IP+I
!c	IF(DABS(B(I)-X(II)).GE.CST05) GO TO 904
if(jq.le.0) go to 919
do 903 i=1,jq
ii=jp+i
if(dabs(b(i)-x(ii)).ge.cst05) go to 904
903 continue
go to 919
904 isfin=1
!c  919 DO 920 I=1,IP
!c  920 A(I)=X(I)
!c	IF(IQ.LE.0) GO TO 925
!c	DO 930 I=1,IQ
!c	II=IP+I
!c  930 B(I)=X(II)
919 do 920 i=1,jp
!xx  920 A(I)=X(I)
a(i)=x(i)
920 continue
if(jq.le.0) go to 925
do 930 i=1,jq
ii=jp+i
!xx  930 B(I)=X(II)
b(i)=x(ii)
930 continue
925 continue
if(isfin.eq.0) go to 940
!c	WRITE(6,926)
!x      IF (IFG.NE.0) WRITE(LU,926)
!xx  926 FORMAT(/1H ,'HESSIAN RESET')
go to 4210
940 continue
!c	WRITE(6,1008) CXX0
!x      IF (IFG.NE.0) WRITE(LU,1008) CXX0
!     HESSIAN COMPUTATION
iswro=1
!c	CALL SC0GRH(X,CYY,G,CN,CXX0,IP,IQ)
!xx      CALL SC0GRH(X,CYY,G(1,NEWN),CN,CXX0,JP,JQ,VD,MMAX,LMAX,ICST,
!xx     *                   IFG,LU)
call sc0grh(x,cyy,g(1,newn),cn,cxx0,jp,jq,vd,mmax,lmax,icst)
!c	WRITE(6,1008) CXX0
!x      IF (IFG.NE.0) WRITE(LU,1008) CXX0
!     INVERSE OF HESSIAN COMPUTATION
!     COMMON SUBROUTINE CALL
!c	CALL MATINV(HESD2,IPQ)
!xx      CALL MATINV(IPQ,VD,MMAX,0,LU)
call matinv(ipq,vd,mmax)
!     INVERSE OF HESSIAN PRINT OUT
!c	WRITE(6,3000)
!c	DO 3100 I=1,IPQ
!c 3100 WRITE(6,3110) I,(VD(I,J),J=1,IPQ)
!x      IF (IFG.NE.0) THEN
!x	 WRITE(LU,3000)
!x	 DO 3100 I=1,IPQ
!x 3100	 WRITE(LU,3110) I,(VD(I,J),J=1,IPQ)
!x      END IF
!     PARAMETER VARIANCE MATRIX COMPUTATION
an=n
const1=cxx0/an
!xx      DO 6000 I=1,IPQ
do 6001 i=1,ipq
do 6000 j=1,ipq
!xx 6000 VD(I,J)=CONST1*VD(I,J)
vd(i,j)=const1*vd(i,j)
6000 continue
6001 continue
!c	WRITE(6,6100)
!c	DO 6200 I=1,IPQ
!c 6200 WRITE(6,3110) I,(VD(I,J),J=1,IPQ)
!x      IF (IFG.NE.0) THEN
!x	 WRITE(LU,6100)
!x	 DO 6200 I=1,IPQ
!x 6200	 WRITE(LU,3110) I,(VD(I,J),J=1,IPQ)
!x      END IF
do 6400 i=1,ipq
if(vd(i,i).lt.cst0) vd(i,i)=cst0
!c 6400 STD(I)=DSQRT(VD(I,I))
!xx 6400 STD(I,NEWN)=DSQRT(VD(I,I))
std(i,newn)=dsqrt(vd(i,i))
6400 continue
!     CN(I)=CXX(I)/CXX(0) I=1,50
!c	WRITE(6,8000)
san=cst2/dsqrt(an)
!c	WRITE(6,7999) SAN
!c	WRITE(6,7998) (CN(I),I=2,51)
!x      IF (IFG.NE.0) THEN
!x	 WRITE(LU,8000) SAN
!x	 WRITE(LU,7998) (CN(I),I=2,MMAX+1)
!x      END IF
!c	DO 800 I=1,IP
!c  800 A(I)=X(I)
!c	IF(IQ.LE.0) GO TO 820
!c	DO 810 I=1,IQ
!c	II=IP+I
!c  810 B(I)=X(II)
do 800 i=1,jp
!xx  800 A(I)=X(I)
a(i)=x(i)
800 continue
if(jq.le.0) go to 820
do 810 i=1,jq
ii=jp+i
!xx  810 B(I)=X(II)
b(i)=x(ii)
810 continue
820 continue
!c	WRITE(6,7910)
!c	WRITE(6,1014) NEWN
!x      IF (IFG.NE.0) WRITE(LU,1014) NEWN
!x 1014 FORMAT(//1H ,'CASE NO.',I5)
!c	IF(IQ.LE.0) GO TO 4290
!c	IF(IQ(NEWN).LE.0) GO TO 4290
if(jq.le.0) go to 4291
!c	WRITE(6,862)
!c	DO 863 I=1,IQ
!c	II=IP+I
!c  863 WRITE(6,864) I,B(I),STD(II)
do 863 i=1,jq
ii=jp+i
863 continue
!     INVERSE OF AR(B) COMPUTATION
ig=0
!c	CALL INVERS(B,IQ,A,0,CC,IB,IG)
!x      CALL INVERS(B,JQ,A,0,CC,IB,IG,IFG,LU)
!xx      CALL INVERS(B,JQ,A,0,CC,IB,ICST,IG,IFG,LU)
call invers(b,jq,a,0,cc,ib,icst,ig)
!c	WRITE(6,4289) IB
!x      IF (IFG.NE.0) WRITE(LU,4289) IB
!c 4290 IF(IP.LE.0) GO TO 4291
!c	WRITE(6,865)
!c	DO 866 I=1,IP
!c  866 WRITE(6,864) I,A(I),STD(I)
4291 continue
!c	WRITE(6,1008) CXX0
saic(newn)=an*dlog(cxx0)+cst2*aipq
!c	WRITE(6,1001) SAIC(NEWN)
!c 1001 FORMAT(/1H ,'AIC=N*LOG(CXX0)+2.0*(IQ+IP)=',D12.5)
!c	WRITE(7,4) NEWN,IQ,IP,SAIC(NEWN)
!c    4 FORMAT(/'CASE NO.',I2,1X,'AR',I2,1X,'MA',I2,2X,'AIC=',D12.5)
!c	WRITE(7,1) IQ
!c	IF(IQ.LE.0) GO TO 4292
!c	WRITE(7,2) (B(I),I=1,IQ)
if(jq.le.0) go to 4292
iq(newn)=jq
do 4200 i=1,jq
b2(i,newn)=b(i)
4200 continue
!c 4292 WRITE(7,1) IP
4292 continue
!c	IF(IP.LE.0) GO TO 1000
!c	WRITE(7,2) (A(I),I=1,IP)
if(jp.le.0) go to 1000
ip (newn)=jp
do 4201 i=1,jp
a2(i,newn)=a(i)
4201 continue
!c 1000 WRITE(7,2) CXX0
1000 cxx2(newn)=cxx0
!     FINAL GRADIENT PRINT OUT
!c	WRITE(6,1010)
!c 1010 FORMAT(/1H ,'FINAL GRADIENT')
!c	WRITE(6,7998) (G(I),I=1,IPQ)
dmaic=smaic-saic(newn)
if(dmaic.lt.cst0) go to 1013
smaic=saic(newn)
!c	IQO=IQ
!c	IPO=IP
!c      DO 1011 I=1,IQ
!c 1011 OB(I)=B(I)
!c      DO 1012 I=1,IP
!c 1012 OA(I)=A(I)
jqo=jq
jpo=jp
do 1011 i=1,jq
!xx 1011 OB(I)=B(I)
ob(i)=b(i)
1011 continue
do 1012 i=1,jp
!xx 1012 OA(I)=A(I)
oa(i)=a(i)
1012 continue
1013 continue
if(newn.ge.newl)  go to 2100
!c 2085 IQ1=IQ+1
!c      IP1=IP+1
2085 jq1=jq+1
jp1=jp+1
!c      DO 2090 I=IQ1,50
!xx      DO 2090 I=JQ1,MMAX
!xx 2090 B(I)=CST0
b(jq1:mmax)=cst0
!c      DO 2095 I=IP1,50
!xx      DO 2095 I=JP1,MMAX
!xx 2095 A(I)=CST0
a(jp1:mmax)=cst0
newn=newn+1
idos=0
!c 2096 IQ=IQI(NEWN)
!c      IF(IP.LE.IPI(NEWN)) GO TO 2976
jq=iqi(newn)
if(jp.le.ipi(newn)) go to 2976
idos=3
!c 2976 IP=IPI(NEWN)
2976 jp=ipi(newn)
go to 4207
!c 2100 WRITE(6,3) SMAIC,IQO,IPO
!c    3 FORMAT(/1H ,'MINUMUM AIC =',D12.5,' ATTAINED AT THE BEST CHOICE
!c     AAR=',I5,'  MA=',I5)
2100 continue
smaic2(newn)=smaic
iqo(newn)=jqo
ipo(newn)=jpo
!     BORDER CHECK
!c 2111 IF(NEWL.GE.25) GO TO 2120
!c 2112 IQM1=IQO-1
!c      IPM1=IPO-1
if(newl.ge.nmax) go to 2120
iqm1=jqo-1
ipm1=jpo-1
idos=0
ido=-1
!c      IP=IPO+1
!c      IQ=IQO+1
jp=jpo+1
jq=jqo+1
go to 2150
2109 ido=0
!c      IP=MAX0(IPM1,1)
!c      IQ=MAX0(IQM1,0)
!c      IF(IPO.LE.IP) GO TO 2150
jp=max0(ipm1,1)
jq=max0(iqm1,0)
if(jpo.le.jp) go to 2150
idos=3
go to 2150
2110 ido=1
!c      IP=IPO
!c      IQ=MAX0(IQM1,0)
jp=jpo
jq=max0(iqm1,0)
go to 2150
2113 ido=2
!c      IQ=IQO+1
jq=jqo+1
go to 2150
2114 ido=3
!c      IQ=IQO
!c      IP=MAX0(IPM1,1)
!c      IF(IPO.LE.IP) GO TO 2150
jq=jqo
jp=max0(ipm1,1)
if(jpo.le.jp) go to 2150
idos=3
go to 2150
2115 ido=4
!c      IP=IPO+1
jp=jpo+1
go to 2150
2116 ido=5
!c      IQ=MAX0(IQM1,0)
jq=max0(iqm1,0)
go to 2150
2117 ido=6
!c      IQ=IQO+1
!c      IP=MAX0(IPM1,1)
!c      IF(IPO.LE.IP) GO TO 2150
jq=jqo+1
jp=max0(ipm1,1)
if(jpo.le.jp) go to 2150
idos=3
2150 do 2151 i=1,newl
!c      IDE=IABS(IQI(I)-IQ)+IABS(IPI(I)-IP)
ide=iabs(iqi(i)-jq)+iabs(ipi(i)-jp)
if(ide.eq.0) go to 2152
2151 continue
go to 2154
2152 idos=0
if(ido.eq.-1) go to 2109
if(ido.eq.0) go to 2110
if(ido.eq.1) go to 2113
if(ido.eq.2) go to 2114
if(ido.eq.3) go to 2115
if(ido.eq.4) go to 2116
if(ido.eq.5) go to 2117
!c      WRITE(6,2153)
!x      IF (IFG.NE.0) WRITE(LU,2153)
!xx 2153 FORMAT(//1H ,'BORDER CHECK COMPLETED')
go to 2120
2154 newl=newl+1
!c      IQI(NEWL)=IQ
!c      IPI(NEWL)=IP
!c      IQ=IQO
!c      IP=IPO
iqi(newl)=jq
ipi(newl)=jp
jq=jqo
jp=jpo
!c      DO 2155 I=1,IQO
do 2155 i=1,jqo
!xx 2155 B(I)=OB(I)
b(i)=ob(i)
2155 continue
!c      DO 2156 I=1,IPO
do 2156 i=1,jpo
!xx 2156 A(I)=OA(I)
a(i)=oa(i)
2156 continue
go to 2085
!c 2120 CALL FLCLS3(NFL)
2120 continue
!xx  999 CONTINUE
continue
aicm=smaic2(newn)
kq=iqo(newn)
kp=ipo(newn)
!x      IF (IFG.NE.0) CLOSE(LU)
return
!xx    1 FORMAT(16I5)
!xx    2 FORMAT(4D20.10)
!xx  862 FORMAT(/1H ,4X,1HI,13X,5HAR(I),1X,'STANDARD DEVIATION')
!xx  864 FORMAT(1H ,I5,2D17.5)
!xx  865 FORMAT(/1H ,4X,1HI,13X,5HMA(I),1X,'STANDARD DEVIATION')
!xx 1600 FORMAT(/1H ,'AUTOMATIC AR-MA MODEL FITTING; SCALAR CASE')
!xx 1601 FORMAT(/1H ,'DAVIDON''S (MINIMIZATION) PROCEDURE')
!xx 1610 FORMAT(/1H ,'INITIAL CONDITION / N=',I5,',LAGH=',I5,',AR-ORDER=',
!xx     A I5,',MA-ORDER=',I5)
!xx 1622 FORMAT(1H ,4X,1HI,12X,5HAR(I))
!xx 1611 FORMAT(1H ,I5,D17.5)
!xx 1621 FORMAT(/1H ,4X,1HI,12X,5HMA(I))
!xx 1008 FORMAT(/1H ,'CXX0=',D12.5)
!xx 3000 FORMAT(/1H ,'INVERSE OF HESSIAN')
!xx 3110 FORMAT(/1H ,I5,4X,10D12.5,/(1H ,9X,10D12.5))
!xx 6100 FORMAT(/1H ,'PARAMETER VARIANCE MATRIX ESTIMATE')
!c 8000 FORMAT(/1H ,'NORMALIZED AUTOCOVARIANCE OF INNOVATION')
!c 7999 FORMAT(1H ,43X,'(2*(INVERSE OF SQUARE ROOT OF N)=',D12.5,')')
!xx 8000 FORMAT(/1H ,'NORMALIZED AUTOCOVARIANCE OF INNOVATION',
!xx     *'  (2*(INVERSE OF SQUARE ROOT OF N)=',D12.5,')')
!xx 7998 FORMAT(1H ,9X,10D12.5/(1H ,9X,10D12.5))
!xx 7910 FORMAT(/1H )
!xx 4289 FORMAT(1H ,7X,'ORDER OF THE INVERSE OF AR=',I5)
!xx11111 FORMAT(//1H ,'PROGRAM 74.1.2. AUTARM')
end
!
!
!c	SUBROUTINE SC0GRH(X,CYY,G,CN,CXX0,IP,IQ)
!xx      SUBROUTINE SC0GRH(X,CYY,G,CN,CXX0,IP,IQ,AL,MM,LL,ICST,IFG,LU)
subroutine sc0grh(x,cyy,g,cn,cxx0,ip,iq,al,mm,ll,icst)
  use timsac_kinds, only: dp
  implicit none
!     THIS SUBROUTINE COMPUTES CXX0,GRADIENT AND HESSIAN.
!xx      IMPLICIT REAL*8(A-H,O-Z)
!c      COMMON /COM50/AL
integer iswro, idos, isik
common /com70/iswro
common /com71/idos
common /com72/isik
!c      DIMENSION X(50),A(50),B(50),AI(190),G(50)
!c      DIMENSION CYY(1001)
!c      DIMENSION CN(51)
!xx      DIMENSION X(IP+IQ),A(IP),B(IQ),AI(ICST),G(IP+IQ)
!xx      DIMENSION CYY(LL*2+1),CN(MM+1)
!xx      DIMENSION AL(MM,MM)
!
!c      DIMENSION A2(100),A2B(190),AIB(190),Y(1001)
!c      DIMENSION CXX(1001),CXY(1001),CUU(1001)
!c      DIMENSION CYX(1001),CUX(1001)
!c      DIMENSION CUZ(1001),CUY(1001)
!c      DIMENSION CYU(1001),CYZ(1001),CZX(1001)
!c      DIMENSION CZY(1001),CZZ(1001)
!c      DIMENSION AL(50,50)
!xx      DIMENSION A2(IP*2),A2B(ICST),AIB(ICST),Y(LL*2+1)
!xx      DIMENSION CXX(LL*2+1)
!xx      DIMENSION CXY(LL*2+1)
!xx      DIMENSION CYX(LL*2+1)
!xx      DIMENSION CUZ(LL*2+1)
!xx      DIMENSION CYU(LL*2+1)
!xx      DIMENSION CZY(LL*2+1)
!c      EQUIVALENCE (CXX(1),CXY(1),CUU(1))
!c      EQUIVALENCE (CYX(1),CUX(1))
!c      EQUIVALENCE (CUZ(1),CUY(1))
!c      EQUIVALENCE (CZY(1),CZZ(1))
!c      EQUIVALENCE (CYU(1),CYZ(1),CZX(1))
!
integer ip, iq, mm, ll, icst
real(dp) x(ip+iq), cyy(ll*2+1), g(ip+iq), cn(mm+1),&
&cxx0, al(mm,mm)
! local
integer i, ii, ia, ij, ij1, ien, ik, ig, i1, i2p, iaib, ia2b,&
&iorig, ipm1, iqm1, ist, j, jj, l2, luu, lux, luy, luyu,&
&luyz, luz, lxx, lxy, lyu, lyx, lyz, lzx, lzy, lzyz, lzz,&
&m2, muu, mux, muy, muyu, muyz, muz, mxx, mxy, myu, myx,&
&myz, mzx, mzy, mzyz, mzz
real(dp) a(ip), b(iq), ai(icst), a2(ip*2), a2b(icst),&
&aib(icst), y(ll*2+1), cxx(ll*2+1), cxy(ll*2+1),&
&cyx(ll*2+1), cuz(ll*2+1), cyu(ll*2+1),&
&czy(ll*2+1), cst0, cst1, dsr2, cai1
!
cst0=0.0d-00
cst1=1.0d-00
!c      IORIG=501
iorig=ll+1
dsr2=0.95d-00
do 100 i=1,ip
!xx  100 A(I)=X(I)
a(i)=x(i)
100 continue
if(iq.le.0) go to 420
do 110 i=1,iq
ii=ip+i
!xx  110 B(I)=X(II)
b(i)=x(ii)
110 continue
420 continue
ig=1
if(isik.ne.0) go to 25
!xx   24 IG=0
ig=0
!     ADJUSTMENT FOR FEASIBLE INITIAL
!     INVERSE OF A(I) COMPUTATION
!c   25 CALL INVERS(A,IP,B,0,AI,IA,IG)
!x   25 CALL INVERS(A,IP,B,0,AI,IA,IG,IFG,LU)
!xx   25 CALL INVERS(A,IP,B,0,AI,IA,ICST,IG,IFG,LU)
25 call invers(a,ip,b,0,ai,ia,icst,ig)
if(isik.eq.0) go to 26
if(iswro.ne.0) go to 1900
if(idos.ne.3) go to 1900
26 if(ig.eq.0) go to 1900
cai1=cst1
do 1941 i=1,ip
cai1=cai1*dsr2
!xx 1941 A(I)=A(I)*CAI1
a(i)=a(i)*cai1
1941 continue
!c      WRITE(6,1940)
!x      IF (IFG.NE.0) WRITE(LU,1940)
!xx 1940 FORMAT(1H ,'NON-INVERTIBLE MA PART')
if(isik.ne.0) go to 1899
ig=0
1899 go to 25
1900 if(ia.ne.0) go to 1901
ia=1
ai(1)=cst0
1901 continue
!
do 2100 i=1,ip
!xx 2100 X(I)=A(I)
x(i)=a(i)
2100 continue
if(iq.le.0) go to 2420
do 2110 i=1,iq
ii=ip+i
!xx 2110 X(II)=B(I)
x(ii)=b(i)
2110 continue
2420 continue
isik=1
ipm1=ip-1
iqm1=iq-1
!     AIB=(INVERSE OF A)*B
ig=0
!c      CALL INVERS(A,IP,B,IQ,AIB,IAIB,IG)
!x      CALL INVERS(A,IP,B,IQ,AIB,IAIB,IG,IFG,LU)
!xx      CALL INVERS(A,IP,B,IQ,AIB,IAIB,ICST,IG,IFG,LU)
call invers(a,ip,b,iq,aib,iaib,icst,ig)
!     A2B=(INVERSE OF A*A)*B
!     A2=A*A
y(iorig)=cst1
do 502 i=1,ip
ij=iorig+i
ik=iorig-i
a2b(i)=a(i)
y(ik)=a(i)
!xx  502 Y(IJ)=CST0
y(ij)=cst0
502 continue
ik=iorig-ip
do 503 i=1,ip
ik=ik-1
!xx  503 Y(IK)=CST0
y(ik)=cst0
503 continue
l2=-ip-ip
m2=-1
!c      CALL SCONVL(Y,A2B,Y,IP,L2,M2)
call sconvl(y,a2b,y,ip,l2,m2,ll)
i2p=ip+ip
do 504 i=1,i2p
ij=iorig-i
!xx  504 A2(I)=Y(IJ)
a2(i)=y(ij)
504 continue
ig=1
!c      CALL INVERS(A2,I2P,B,IQ,A2B,IA2B,IG)
!x      CALL INVERS(A2,I2P,B,IQ,A2B,IA2B,IG,IFG,LU)
!xx      CALL INVERS(A2,I2P,B,IQ,A2B,IA2B,ICST,IG,IFG,LU)
call invers(a2,i2p,b,iq,a2b,ia2b,icst,ig)
lxx=0
!c      MXX=50
mxx=mm
lxy=0
mxy=mxx+iaib
lyx=-mxy
myx=-lxy
lzx=-iq
mzx=-1
lzy=lzx
mzy=mzx+iaib
lzz=0
mzz=iqm1
lzyz=lzz
mzyz=mzz+ia
lzy=min0(lzy,lzyz)
mzy=max0(mzy,mzyz)
lyz=-mzy
myz=-lzy
lux=-ip
mux=-1
luy=-ip
muy=mux+iaib
luu=0
muu=ipm1
luyu=0
muyu=muu+ia2b
luy=min0(luy,luyu)
muy=max0(muy,muyu)
luz=-ipm1
muz=iqm1
luyz=luz
muyz=muz+ia
luy=min0(luy,luyz)
muy=max0(muy,muyz)
lyu=-muy
myu=-luy
!c      WRITE(6,2502)
!c 2500 WRITE(6,3000) LXX,LXY,LZX,LZY,LUX,LUU,LUZ,LUY,LUYU,LZYZ,IA,IAIB,
!c     AIA2B
!c 2501 WRITE(6,3000) MXX,MXY,MZX,MZY,MUX,MUU,MUZ,MUY,MUYU,MZYZ
!x      IF (IFG.NE.0) THEN
!x      WRITE(LU,2502)
!x 2500	 WRITE(LU,3000) LXX,LXY,LZX,LZY,LUX,LUU,LUZ,LUY,LUYU,LZYZ,
!x     *	 IA,IAIB,IA2B
!x 2501	 WRITE(LU,3000) MXX,MXY,MZX,MZY,MUX,MUU,MUZ,MUY,MUYU,MZYZ
!x      END IF
!     CXX0 COMPUTATION
!     CYX=CYY*AIB'
!c	CALL SCONVL(CYY,AIB,CYX,IAIB,LYX,MYX)
call sconvl(cyy,aib,cyx,iaib,lyx,myx,ll)
!
!c      CALL TURN(CYX,CXY,LYX,MYX)
do 505 i=1,ll*2+1
cxy(i)=cxx(i)
505 continue
call turn(cyx,cxy,lyx,myx,ll)
!     CXX=CXY*AIB'
!c      CALL SCONVL(CXY,AIB,CXX,IAIB,LXX,MXX)
call sconvl(cxy,aib,cxx,iaib,lxx,mxx,ll)
do 506 i=1,ll*2+1
cxx(i)=cxy(i)
506 continue
cxx0=cxx(iorig)
ist=iorig+lxx
ien=iorig+mxx
ij=0
do 510 i=ist,ien
ij=ij+1
!xx  510 CN(IJ)=CXX(I)/CXX0
cn(ij)=cxx(i)/cxx0
510 continue
!     GA COMPUTATION
!     CYU=CYY*A2B'
!c      CALL SCONVL(CYY,A2B,CYU,IA2B,LYU,MYU)
!c      CALL TURN(CYU,CUY,LYU,MYU)
!c      CALL SCONVL(CUY,AIB,CUX,IAIB,LUX,MUX)
call sconvl(cyy,a2b,cyu,ia2b,lyu,myu,ll)
call turn(cyu,cuz,lyu,myu,ll)
call sconvl(cuz,aib,cyx,iaib,lux,mux,ll)
!     HAA
!c      CALL SCONVL(CUY,A2B,CUU,IA2B,LUU,MUU)
call sconvl(cuz,a2b,cxx,ia2b,luu,muu,ll)
if(iq.eq.0) go to 550
!     HAB
!c      CALL SCONVL(CUY,AI,CUZ,IA,LUZ,MUZ)
call sconvl(cuz,ai,cuz,ia,luz,muz,ll)
!     GB COMPUTATION
!     CYZ=CYY*AI'
!     IF(IQ.EQ.0) GO TO 550
!c      CALL SCONVL(CYY,AI,CYZ,IA,LYZ,MYZ)
!c      CALL TURN(CYZ,CZY,LYZ,MYZ)
!c      CALL SCONVL(CZY,AIB,CZX,IAIB,LZX,MZX)
call sconvl(cyy,ai,cyu,ia,lyz,myz,ll)
call turn(cyu,czy,lyz,myz,ll)
call sconvl(czy,aib,cyu,iaib,lzx,mzx,ll)
!     HBB
!c      CALL SCONVL(CZY,AI,CZZ,IA,LZZ,MZZ)
call sconvl(czy,ai,czy,ia,lzz,mzz,ll)
!     HESSIAN ARRANGEMENT FOR U
550 continue
do 211 i=1,ip
do 212 j=1,i
ij1=iorig+i-j
!c      AL(I,J)=CUU(IJ1)
al(i,j)=cxx(ij1)
!xx  212 AL(J,I)=AL(I,J)
al(j,i)=al(i,j)
212 continue
211 continue
if(iq.le.0) go to 4220
!     HESSIAN ARRANGEMENT FOR V
do 231 i=1,iq
ii=ip+i
do 232 j=1,i
jj=ip+j
ij1=iorig+i-j
!c      AL(II,JJ)=CZZ(IJ1)
al(ii,jj)=czy(ij1)
!xx  232 AL(JJ,II)=AL(II,JJ)
al(jj,ii)=al(ii,jj)
232 continue
231 continue
!     HESSIAN ARRANGEMENT FOR -W AND -W'
do 251 i=1,iq
ii=ip+i
do 252 j=1,ip
ij1=iorig+i-j
al(ii,j)=-cuz(ij1)
!xx  252 AL(J,II)=AL(II,J)
al(j,ii)=al(ii,j)
252 continue
251 continue
4220 continue
!     GRADIENT ARRANGEMENT
do 280 i=1,ip
i1=iorig-i
!c  280 G(I)=-CUX(I1)
!xx  280 G(I)=-CYX(I1)
g(i)=-cyx(i1)
280 continue
if(iq.le.0) go to 4230
do 281 i=1,iq
ii=ip+i
i1=iorig-i
!c  281 G(II)=CZX(I1)
!xx  281 G(II)=CYU(I1)
g(ii)=cyu(i1)
281 continue
idos=0
4230 return
!xx 2502 FORMAT(/1H ,'PARAMETER PRINT OUT AT THE STATEMENT NUMBER',
!xx     A' 2500-2501 OF SC0GRH')
!xx 3000 FORMAT(1H ,16I5)
end
!
!xx      SUBROUTINE SC0GR1(X,CYY,G,CXX0,IP,IQ,IG,LL,ICST,IFG,LU)
subroutine sc0gr1(x,cyy,g,cxx0,ip,iq,ig,ll,icst)
  use timsac_kinds, only: dp
  implicit none
!     THIS SUBROUTINE COMPUTES CXX0 AND GRADIENT.
!xx      IMPLICIT REAL*8(A-H,O-Z)
!c      DIMENSION X(50),A(50),B(50),AI(190),G(50)
!c      DIMENSION A2(100),A2B(190),AIB(190),Y(1001)
!c      DIMENSION CYY(1001)
!c      DIMENSION CXX(1001),CXY(1001)
!c      DIMENSION CYX(1001),CUX(1001)
!c      DIMENSION CUY(1001)
!c      DIMENSION CYU(1001),CYZ(1001),CZX(1001)
!c      DIMENSION CZY(1001)
!xx      DIMENSION X(IP+IQ),A(IP),B(IQ),AI(ICST),G(IP+IQ)
!xx      DIMENSION A2(IP*2),A2B(ICST),AIB(ICST),Y(LL*2+1)
!xx      DIMENSION CYY(LL*2+1),CXX(LL*2+1),CYX(LL*2+1)
!xx      DIMENSION CUY(LL*2+1),CYU(LL*2+1),CZY(LL*2+1)
!c      EQUIVALENCE (CXX(1),CXY(1))
!c      EQUIVALENCE (CYX(1),CUX(1))
!c      EQUIVALENCE (CYU(1),CYZ(1),CZX(1))
integer ip, iq, ig, ll, icst
real(dp) x(ip+iq), cyy(ll*2+1), g(ip+iq), cxx0
! local
integer i, ii, i1, ia, ib, ij, ik, iaib, igaib, iorig, iga2b,&
&ia2b, ipq, ipm1, iqm1, i2p, l2, lxx, mxx, lxy, mxy, lyx,&
&myx, lzx, mzx, lzy, mzy, lyz, m2, myz, lux, mux, luy, muy,&
&lyu, myu
real(dp) a(ip), b(iq), ai(icst), a2(ip*2), a2b(icst),&
&aib(icst), y(ll*2+1), cxx(ll*2+1), cyx(ll*2+1),&
&cuy(ll*2+1), cyu(ll*2+1), czy(ll*2+1), cst0, cst1
!
cst0=0.0d-00
cst1=1.0d-00
!c      IORIG=501
iorig=ll+1
iga2b=ig
do 100 i=1,ip
!xx  100 A(I)=X(I)
a(i)=x(i)
100 continue
if(iq.le.0) go to 420
do 110 i=1,iq
ii=ip+i
!xx  110 B(I)=X(II)
b(i)=x(ii)
110 continue
420 continue
ib=iq
!     INVERSE OF A(I) COMPUTATION
!c   24 CALL INVERS(A,IP,B,0,AI,IA,IG)
!x   24 CALL INVERS(A,IP,B,0,AI,IA,IG,IFG,LU)
!xx   24 CALL INVERS(A,IP,B,0,AI,IA,ICST,IG,IFG,LU)
call invers(a,ip,b,0,ai,ia,icst,ig)
if(ig.ne.1) go to 1900
go to 1000
1900 if(ia.ne.0) go to 1901
ia=1
ai(1)=cst0
1901 continue
do 2100 i=1,ip
!xx 2100 X(I)=A(I)
x(i)=a(i)
2100 continue
if(iq.le.0) go to 2420
do 2110 i=1,iq
ii=ip+i
!xx 2110 X(II)=B(I)
x(ii)=b(i)
2110 continue
2420 continue
ipm1=ip-1
iqm1=iq-1
!     AIB=(INVERSE OF A)*B
igaib=0
!c   25 CALL INVERS(A,IP,B,IQ,AIB,IAIB,IGAIB)
!x   25 CALL INVERS(A,IP,B,IQ,AIB,IAIB,IGAIB,IFG,LU)
!xx   25 CALL INVERS(A,IP,B,IQ,AIB,IAIB,ICST,IGAIB,IFG,LU)
call invers(a,ip,b,iq,aib,iaib,icst,igaib)
!     A2B=(INVERSE OF A*A)*B
!     A2=A*A
y(iorig)=cst1
do 502 i=1,ip
ij=iorig+i
ik=iorig-i
a2b(i)=a(i)
y(ik)=a(i)
!xx  502 Y(IJ)=CST0
y(ij)=cst0
502 continue
ik=iorig-ip
do 503 i=1,ip
ik=ik-1
!xx  503 Y(IK)=CST0
y(ik)=cst0
503 continue
l2=-ip-ip
m2=-1
!c      CALL SCONVL(Y,A2B,Y,IP,L2,M2)
call sconvl(y,a2b,y,ip,l2,m2,ll)
i2p=ip+ip
do 504 i=1,i2p
ij=iorig-i
!xx  504 A2(I)=Y(IJ)
a2(i)=y(ij)
504 continue
!c      CALL INVERS(A2,I2P,B,IQ,A2B,IA2B,IGA2B)
!x      CALL INVERS(A2,I2P,B,IQ,A2B,IA2B,IGA2B,IFG,LU)
!xx      CALL INVERS(A2,I2P,B,IQ,A2B,IA2B,ICST,IGA2B,IFG,LU)
call invers(a2,i2p,b,iq,a2b,ia2b,icst,iga2b)
lxx=0
mxx=0
lxy=0
mxy=mxx+iaib
lyx=-mxy
myx=-lxy
lzx=-iq
mzx=-1
lzy=lzx
mzy=mzx+iaib
lyz=-mzy
myz=-lzy
lux=-ip
mux=-1
luy=-ip
muy=mux+iaib
lyu=-muy
myu=-luy
!     CXX0 COMPUTATION
!     CYX=CYY*AIB'
!c      CALL SCONVL(CYY,AIB,CYX,IAIB,LYX,MYX)
!c      CALL TURN(CYX,CXY,LYX,MYX)
call sconvl(cyy,aib,cyx,iaib,lyx,myx,ll)
call turn(cyx,cxx,lyx,myx,ll)
!     CXX=CXY*AIB'
!c      CALL SCONVL(CXY,AIB,CXX,IAIB,LXX,MXX)
call sconvl(cxx,aib,cxx,iaib,lxx,mxx,ll)
cxx0=cxx(iorig)
!     GA COMPUTATION
!     CYU=CYY*A2B'
!c      CALL SCONVL(CYY,A2B,CYU,IA2B,LYU,MYU)
!c      CALL TURN(CYU,CUY,LYU,MYU)
!c      CALL SCONVL(CUY,AIB,CUX,IAIB,LUX,MUX)
call sconvl(cyy,a2b,cyu,ia2b,lyu,myu,ll)
call turn(cyu,cuy,lyu,myu,ll)
call sconvl(cuy,aib,cyx,iaib,lux,mux,ll)
!     GB COMPUTATION
!     CYZ=CYY*AI'
if(iq.eq.0) go to 5279
!c      CALL SCONVL(CYY,AI,CYZ,IA,LYZ,MYZ)
!c      CALL TURN(CYZ,CZY,LYZ,MYZ)
!c      CALL SCONVL(CZY,AIB,CZX,IAIB,LZX,MZX)
call sconvl(cyy,ai,cyu,ia,lyz,myz,ll)
call turn(cyu,czy,lyz,myz,ll)
call sconvl(czy,aib,cyu,iaib,lzx,mzx,ll)
5279 continue
!     GRADIENT ARRANGEMENT
do 5280 i=1,ip
i1=iorig-i
!c 5280 G(I)=-CUX(I1)
!xx 5280 G(I)=-CYX(I1)
g(i)=-cyx(i1)
5280 continue
if(iq.le.0) go to 5290
do 5281 i=1,iq
ii=ip+i
i1=iorig-i
!c 5281 G(II)=CZX(I1)
!xx 5281 G(II)=CYU(I1)
g(ii)=cyu(i1)
5281 continue
5290 continue
!     CXX0, GRADIENT PRINT OUT
ipq=ip+iq
1000 return
end
!
!c      SUBROUTINE SDAV1(X,CYY,CXX0,G,C,IP,IQ,N)
!xx      SUBROUTINE SDAV1(X,CYY,CXX0,G,C,IP,IQ,N,VD,NN,LL,ICST,IFG,LU)
subroutine sdav1(x,cyy,cxx0,g,c,ip,iq,n,vd,nn,ll,icst)
  use timsac_kinds, only: dp
  implicit none
!      DADIDON'S (MINIMIZATION) PROCEDURE
!xx      IMPLICIT REAL*8(A-H,O-Z)
!c      COMMON /COM50/VD
integer iswro
common /com70/iswro
!c      DIMENSION VD(50,50)
!c      DIMENSION X(50),G(50),SX(50),SG(50),SR(50)
!c      DIMENSION C(50)
!xx      DIMENSION X(IP+IQ),CYY(LL*2+1)
!xx      DIMENSION G(IP+IQ),C(IP+IQ),VD(NN,NN)
!xx      DIMENSION SX(IP+IQ),SG(IP+IQ),SR(IP+IQ)
integer ip, iq, n, nn, ll, icst
real(dp) x(ip+iq), cyy(ll*2+1), cxx0, g(ip+iq), c(ip+iq),&
&vd(nn,nn)
! local
integer i, ig, isphai, itn, itns, iphai, ipq, ipq2, iram, j
real(dp) sx(ip+iq), sg(ip+iq), sr(ip+iq), cst0, cst1,&
&cst2, cst05, consta, constb, eps1, eps3, eps4,&
&aipq, an, phai, ephai1, t1, ro, ram, ramro,&
&ramrot, sum, sro, srod, dgam, dgam1, gsr, ramt,&
&ramsro, ram1, consdr, sphai, oaic, ophai, aic,&
&daic
!
!     CONSTANT
cst0=0.0d-00
cst1=1.0d-00
cst2=2.0d-00
cst05=0.5d-00
consta=0.5d-00
constb=2.0d-00
eps1=0.01d-00
eps3=0.000001d-00
eps4=0.1d-10
isphai=0
itn=1
iphai=1
ipq=ip+iq
aipq=ipq
an=n
phai=cxx0
150 continue
!     RO=G'*C COMPUTATION
itns=0
!     COMMON SUBROUTINE CALL
40 call innerp(g,c,ro,ipq)
if(iphai.eq.0) go to 101
phai=cxx0
101 ophai=phai
ephai1=eps1*phai
t1=ro-cst2*phai
if(t1.le.ephai1) go to 140
ram=cst2*phai/ro
!     V=V+((RAM-1.0)/RO)*(C*C')
ramro=(ram-cst1)/ro
!xx      DO 110 I=1,IPQ
do 111 i=1,ipq
ramrot=ramro*c(i)
do 110 j=1,ipq
!xx  110 VD(I,J)=VD(I,J)+RAMROT*C(J)
vd(i,j)=vd(i,j)+ramrot*c(j)
110 continue
111 continue
!     C=RAM*C
do 120 i=1,ipq
!xx  120 C(I)=RAM*C(I)
c(i)=ram*c(i)
120 continue
if(itns.ge.10) go to 140
itns=itns+1
go to 40
!     SX=X-R
140 continue
ig=0
1210 continue
do 210 i=1,ipq
!xx  210 SX(I)=X(I)-C(I)
sx(i)=x(i)-c(i)
210 continue
!     SPHAI=CXX0, SG=GRADIENT COMPUTATION
!c      CALL SC0GR1(SX,CYY,SG,SPHAI,IP,IQ,IG,LL,ICST)
!xx      CALL SC0GR1(SX,CYY,SG,SPHAI,IP,IQ,IG,LL,ICST,IFG,LU)
call sc0gr1(sx,cyy,sg,sphai,ip,iq,ig,ll,icst)
if(ig.ne.1) go to 309
!xx      DO 303 I=1,IPQ
do 304 i=1,ipq
c(i)=cst05*c(i)
do 303 j=1,ipq
!xx  303 VD(I,J)=CST05*VD(I,J)
vd(i,j)=cst05*vd(i,j)
303 continue
304 continue
go to 1210
309 continue
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
!     SRO=(SG)'*(SR)
!     COMMON SUBROUTINE CALL
call innerp(sg,sr,sro,ipq)
srod=sro/phai
!     DGAM=-G'*(SR)/SRO
!     COMMON SUBROUTINE CALL
call innerp(g,sr,gsr,ipq)
dgam=-gsr/sro
dgam1=dgam+cst1
dgam1=dabs(dgam1)+0.1d-70
ram=dabs(dgam)/dgam1
!     IF RAM . LE. CONSTA THEN RAM=CONSTA
if(ram.gt.consta) go to 430
ram=consta
iram=1
go to 470
!     IF RAM . GE. CONSTB THEN RAM=CONSTB
430 if(ram.lt.constb) go to 450
ram=constb
iram=-1
go to 470
!     RAM=RAM
450 continue
iram=0
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
if(phai.ge.sphai) go to 540
!     SPHAI.GT.PHAI: TEST OF CORRECTION
ram1=ram-cst1
if(dabs(ram1).lt.eps3) go to 555
consdr=dgam*ram1
do 550 i=1,ipq
!xx  550 C(I)=C(I)-CONSDR*SR(I)
c(i)=c(i)-consdr*sr(i)
550 continue
iphai=0
if(srod.gt.eps4) go to 900
!     END OF ITERATION
555 iswro=iswro+1
go to 1000
!     SPHAI LE. PHAI: SUCCESSFUL REDUCTION
540 do 560 i=1,ipq
x(i)=sx(i)
g(i)=sg(i)
!xx  560 C(I)=RAM*SR(I)
c(i)=ram*sr(i)
560 continue
cxx0=sphai
phai=sphai
iphai=1
!xx  800 CONTINUE
oaic=an*dlog(ophai)+cst2*aipq
aic=an*dlog(phai)+cst2*aipq
daic=oaic-aic
if(iram.ne.0) go to 901
if(srod.lt.eps4) go to 555
!     ITERATION CHECK
900 ipq2=ipq+ipq
if(itn.ge.ipq2) go to 555
isphai=(isphai+(1-iphai))*(1-iphai)
if(isphai.gt.10) go to 555
itn=itn+1
go to 150
901 if(srod.lt.eps4) go to 555
go to 900
!     END OF MINIMIZATION
!xx  999 ISWRO=0
iswro=0
1000 continue
!xx 1001 RETURN
return
end
!
!c      SUBROUTINE SCONVL(Y,A,Z,K,L,M)
subroutine sconvl(y,a,z,k,l,m,ll)
  use timsac_kinds, only: dp
  implicit none
!     Y(I), Z(I) CENTERED AT I=IORIG
!     Z(I)=Y(I)+Y(I+1)A(1)+...+Y(I+K)A(K) (I=L,M)
!xx      IMPLICIT REAL*8(A-H,O-Z)
!c      DIMENSION Y(1001),A(190),Z(1001)
!xx      DIMENSION
integer k, l, m, ll
real(dp) y(ll*2+1), a(k), z(ll*2+1)
! local
integer i, ij, iorig, ist, ien, j
real(dp) sum
!c      IORIG=501
iorig=ll+1
ist=iorig+l
ien=iorig+m
do 3 i=ist,ien
sum=y(i)
do 2 j=1,k
ij=i+j
!xx    2 SUM=SUM+Y(IJ)*A(J)
sum=sum+y(ij)*a(j)
2 continue
!xx    3 Z(I)=SUM
z(i)=sum
3 continue
return
end
!
!c      SUBROUTINE TURN(Y,Z,L,M)
subroutine turn(y,z,l,m,ll)
  use timsac_kinds, only: dp
  implicit none
!     Z(IORIG+I)=Y(IORIG-I) (I=1,M)
!xx      IMPLICIT REAL*8(A-H,O-Z)
!c      DIMENSION Z(1001),Y(1001)
!xx      DIMENSION Z(LL-L+1),Y(LL+M+1)
integer l, m, ll
real(dp) y(ll+m+1), z(ll-l+1)
! local
integer i, iorig, ist, ien, ij
!c      IORIG=501
iorig=ll+1
ist=iorig+l
ien=iorig+m
do 1 i=ist,ien
ij=iorig-(i-iorig)
!xx    1 Z(IJ)=Y(I)
z(ij)=y(i)
1 continue
return
end
!
!c      SUBROUTINE INVERS(A,IP,B,IQ,X,IX,IG)
!x      SUBROUTINE INVERS(A,IP,B,IQ,X,IX,IG,IFG,LU)
!xx      SUBROUTINE INVERS(A,IP,B,IQ,X,IX,ICST,IG,IFG,LU)
subroutine invers(a,ip,b,iq,x,ix,icst,ig)
  use timsac_kinds, only: dp
  implicit none
!     X=(INVERSE OF B )*A
!     W(I)+B(1)W(I-1)+...B(IQ)W(I-IQ)=X(I)+A(1)X(I-1)+...+A(IP)X(I-IP)
!     INPUT W(0)=1, W(I)=0 FOR I. NE. 0
!     OUTPUT X(I) (I=1,IX)
!xx      IMPLICIT REAL*8(A-H,O-Z)
!x      DIMENSION A(1),B(1),X(1)
!xx      DIMENSION A(IP),B(IQ),X(ICST)
integer ip, iq, ix, icst, ig
real(dp) a(ip), b(iq), x(icst)
! local
integer i, ih, im1, imj, ipq, j, jm, k, lh
real(dp) cst0, gconst, gammax, sum, gam2
cst0=0.0d-00
ipq=ip+iq
if(ipq.le.0) go to 999
gconst=0.0005d-00
gammax=1.0d+10
k=0
lh=6
ih=0
if(ig.eq.0) go to 13
gconst=0.01d-00
ig=0
!x   13 DO 10 I=1,190
13 do 10 i=1,icst
ix=i
sum=cst0
if(i.gt.iq) go to 2
sum=b(i)
2 if(i.gt.ip) go to 3
sum=sum-a(i)
3 im1=i-1
jm=min0(im1,ip)
if(jm.le.0) go to 5
do 4 j=1,jm
imj=i-j
!xx    4 SUM=SUM-X(IMJ)*A(J)
sum=sum-x(imj)*a(j)
4 continue
5 x(i)=sum
gam2=dabs(sum)
if(gam2.ge.gconst) go to 24
ih=ih+1
if(ih.lt.lh) go to 10
go to 1000
24 if(gam2.le.gammax) go to 26
ig=1
!c	WRITE(6,60)
!x      IF (IFG.NE.0) WRITE(LU,60)
go to 1000
26 ih=0
10 continue
if(ih.ge.lh) go to 1000
ig=1
!c      WRITE(6,59)
!x      IF (IFG.NE.0) WRITE(LU,59)
go to 1000
999 ix=0
1000 return
!xx   59 FORMAT(1H ,'INCOMPLETE CONVERGENCE OF INVERSE')
!xx   60 FORMAT(1H ,'DIVERGENT INVERSE')
end
