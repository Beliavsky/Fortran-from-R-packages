! Part of the experimental modern Fortran port of timsac 1.3.8-6.
! This file was created or modified for the Fortran project on 2026-07-23.
! Original TIMSAC credits are retained; see NOTICE and ORIGIN.md.
! SPDX-License-Identifier: GPL-2.0-or-later

subroutine fpec7(n,l,ir,ip,ip0,inw,r1,r2,fpec,rfpec,aic,&
&ifpec,ofpec,orfpec,oaic,osd,ao)
  use timsac_kinds, only: dp
  implicit none
!
!     PROGRAM 5.3.2   FPEC(AR-MODEL FITTING FOR CONTROL)
!-----------------------------------------------------------------------
!      SUBROUTINE FPEC7F(N,L,IR,IP,IP0,LAGH1,INW,MJ,MJ0,
!      SUBROUTINE RECOVA(X,LAGH1,L1,IP0)
!      SUBROUTINE SFPEC(SD,N,K,IR,MS,Z,RZ,OOZ,AIC)
!-----------------------------------------------------------------------
!     ** DESIGNED BY H. AKAIKE, THE INSTITUTE OF STATISTICAL MATHEMATICS
!     ** PROGRAMMED BY E. ARAHATA, THE INSTITUTE OF STATISTICAL MATHEMAT
!         TOKYO
!     ** DATE OF THE LATEST REVISION: MARCH 25, 1977
!     ** THIS PROGRAM WAS ORIGINALLY PUBLISHED IN
!         "DAINAMIKKU SISTEMU NO TOKEI-TEKI KAISEKI TO SEIGYO (STATISTICA
!         ANALYSIS AND CONTROL OF DYNAMIC SYSTEMS)" BY H. AKAIKE AND
!         T. NAKAGAWA, SAIENSU-SHA, TOKYO, 1972 (IN JAPANESE)
!-----------------------------------------------------------------------
!     THIS PROGRAM PERFORMS FPEC(AR-MODEL FITTING FOR CONTROL)
!     COMPUTATION.
!     BESIDES THE OUTPUTS OF PROGRAM 5.1.2   MULCOR, THE FOLLOWING
!     INPUTS ARE REQUIRED:
!     L: UPPER LIMIT OF MODEL ORDER M (LESS THAN 30)
!     IR: NUMBER OF CONTROLLED VARIABLES
!     IL: NUMBER OF MANINPULATED VARIABLES, IL=0 FOR MFPE COMPUTATION
!     INW(I): INDICATOR; FIRST IR INDICATE THE CONTROLLED VARIABLES
!     AND THE REST THE MANIPULATE VARIABLES WITHIN THE IP0 VARIABLES
!     IN THE OUTPUT OF PROGRAM 5.1.2   MULCOR.
!     THE OUTPUTS ARE THE PREDICTION ERROR COVARIANCE MATRIX OSD AND
!     THE SET OF COEFFICIENT MATRICES A AND B TO BE USED IN
!     PROGRAM 5.5.1   OPTIMAL CONTROLLER DESIGN.
!
!xx      IMPLICIT REAL*8(A-H,O-Z)
!      DIMENSION R1(30,10,10),A1(30,10,10),B1(30,10,10),AO(30,10,10)
!      DIMENSION SD(10,10),SE(10,10),SF(10,10),OSD(10,10)
!      DIMENSION XSD(10,10),XSF(10,10),D(10,10),E(10,10),Z1(10,10)
!      DIMENSION INW(10),C1(10,10)
!xx      DIMENSION INW(IP),R1(L+1,IP0,IP0),R2(L+1,IP,IP)
!xx      DIMENSION FPEC(0:L),RFPEC(0:L),AIC(0:L)
!xx      DIMENSION OSD(IR,IR),AO(L,IR,IP)
!xx      DIMENSION A1(L,IP,IP),B1(L,IP,IP),C1(IP0,IP0)
!xx      DIMENSION SD(IP,IP),SE(IP,IP),SF(IP,IP)
!xx      DIMENSION XSD(IP,IP),XSF(IP,IP),D(IP,IP),E(IP,IP),Z1(IP,IP)
integer n, l, ir, ip, ip0, inw(ip), ifpec
real(dp) r1(l+1,ip0,ip0), r2(l+1,ip,ip), fpec(0:l),&
&rfpec(0:l), aic(0:l), ofpec, orfpec, oaic,&
&osd(ir,ir), ao(l,ir,ip)
! local
integer i, ii, j, jj, l1, m, ms
real(dp) a1(l,ip,ip), b1(l,ip,ip), c1(ip0,ip0), sd(ip,ip),&
&se(ip,ip), sf(ip,ip), xsd(ip,ip), xsf(ip,ip),&
&d(ip,ip), e(ip,ip), z1(ip,ip), oofpec, sddet,&
&sfdet
!
!     INPUT / OUTPUT DATA FILE OPEN
!      CHARACTER(100) DFNAM
!      DFNAM='fpec.out'
!      CALL SETWND
!      CALL FLOPN3(DFNAM,NFL)
!      IF (NFL.EQ.0) GO TO 999
!     INITIAL CONDITION INPUT
!      READ(5,1) L,IR,IL
!      IP=IR+IL
!      READ(5,1) (INW(I),I=1,IP)
!     READING THE OUTPUTS OF PROGRAM 5.1.2 MULCOR
!      READ(5,1) N,LAGH,IP0
l1=l+1
!      LAGH1=LAGH+1
!      CALL RECOVA(R1,LAGH1,L1,IP0,MJ0,MJ)
do 10 ii=1,l1
!xx      DO 20 I=1,IP0
do 22 i=1,ip0
do 20 j=1,ip0
!xx   20 C1(I,J)=R1(II,I,J)
c1(i,j)=r1(ii,i,j)
20 continue
22 continue
!     MATRIX REARRANGEMENT BY INW
!        CALL REARRA(C1,INW,IP0,IP,MJ)
call rearra(c1,inw,ip0,ip)
!xx      DO 21 I=1,IP
do 23 i=1,ip
do 21 j=1,ip
!xx   21 R2(II,I,J)=C1(I,J)
r2(ii,i,j)=c1(i,j)
21 continue
23 continue
10 continue
!     INITIAL CONDITION AND COVARIANCE PRINT OUT
!      WRITE(6,39)
!      WRITE(6,40)
!      WRITE(6,41) N,L,IR,IL
!      WRITE(6,259) (INW(I),I=1,1,IP)
!      WRITE(6,42)
!      CALL PRMAT3(R1,L1,IP,IP,1,MJ0,MJ,MJ)
!     INITIAL SD, SF, SE COMPUTATION
!xx      DO 330 II=1,IP
do 331 ii=1,ip
do 330 jj=1,ip
!      SD(II,JJ)=R1(1,II,JJ)
sd(ii,jj)=r2(1,ii,jj)
sf(ii,jj)=sd(ii,jj)
!      SE(II,JJ)=R1(2,II,JJ)
se(ii,jj)=r2(2,ii,jj)
xsd(ii,jj)=sd(ii,jj)
!xx  330 XSF(II,JJ)=SF(II,JJ)
xsf(ii,jj)=sf(ii,jj)
330 continue
331 continue
!     0-TH STEP COMPUTATION
ifpec=0
ms=0
!     OFPEC, ORFPEC COMPUTATION
!      CALL SFPEC(SD,N,IP,IR,MS,OFPEC,ORFPEC,OOFPEC,MJ)
call sfpec(sd,n,ip,ir,ms,fpec(0),rfpec(0),oofpec,aic(0))
!     OFPEC, ORFPEC PRINT OUT
!      WRITE(6,600)
!      WRITE(6,264) MS,OFPEC,ORFPEC,AIC
!      OAIC=AIC
oaic=aic(0)
ofpec=fpec(0)
orfpec=rfpec(0)
!     ITERATION M=1 TO L
do 400 m=1,l
!     INVERSE OF SD, SF COMPUTATION
!      CALL INVDET(XSD,SDDET,IP,MJ)
!      CALL INVDET(XSF,SFDET,IP,MJ)
call invdet(xsd,sddet,ip,ip)
call invdet(xsf,sfdet,ip,ip)
!     D, E, SD, SF COMPUTATION
!      CALL MULPLY(SE,XSF,D,IP,IP,IP,MJ,MJ,MJ)
!      CALL TRAMDL(SE,XSD,E,IP,IP,IP,MJ,MJ,MJ)
!      CALL TRAMDR(D,SE,Z1,IP,IP,IP,MJ,MJ,MJ)
!      CALL SUBTAL(SD,Z1,IP,IP,MJ,MJ)
!      CALL MULPLY(E,SE,Z1,IP,IP,IP,MJ,MJ,MJ)
!      CALL SUBTAL(SF,Z1,IP,IP,MJ,MJ)
call mulply(se,xsf,d,ip,ip,ip)
call tramdl(se,xsd,e,ip,ip,ip)
call tramdr(d,se,z1,ip,ip,ip)
call subtal(sd,z1,ip,ip)
call mulply(e,se,z1,ip,ip,ip)
call subtal(sf,z1,ip,ip)
ms=m
!xx      DO 410 II=1,IP
do 411 ii=1,ip
do 410 jj=1,ip
xsd(ii,jj)=sd(ii,jj)
!xx  410 XSF(II,JJ)=SF(II,JJ)
xsf(ii,jj)=sf(ii,jj)
410 continue
411 continue
!     FPEC,RFPEC COMPUTATION
!      CALL SFPEC(SD,N,IP,IR,MS,FPEC,RFPEC,OOFPEC,MJ)
call sfpec(sd,n,ip,ir,ms,fpec(m),rfpec(m),oofpec,aic(m))
!     FPEC,RFPEC PRINT OUT
!      WRITE(6,264) MS,FPEC,RFPEC,AIC
!     FORWARD AND BACKWARD PREDICTOR COMPUTATION
!      CALL COEFAB(A1,B1,D,E,MS(M),IP,MJ0,MJ)
call coefab(a1,b1,d,e,ms,l,ip)
!     MIN.FPEC, MIN.RFPEC COMPUTATION
!      IF(OFPEC.LE.FPEC) GO TO 440
!      OAIC=AIC
!      OFPEC=FPEC
!      ORFPEC=RFPEC
if(ofpec.le.fpec(m)) go to 440
oaic=aic(m)
ofpec=fpec(m)
orfpec=rfpec(m)
ifpec=m
!xx      DO 560 II=1,IR
do 566 ii=1,ir
do 560 jj=1,ir
!xx  560 OSD(II,JJ)=SD(II,JJ)
osd(ii,jj)=sd(ii,jj)
560 continue
566 continue
do 561 i=1,m
do 562 ii=1,ir
!xx      DO 562 JJ=1,IP
do 563 jj=1,ip
!xx  562 AO(I,II,JJ)=A1(I,II,JJ)
ao(i,ii,jj)=a1(i,ii,jj)
563 continue
562 continue
561 continue
440 if(m.eq.l) go to 400
!     SE COMPUTATION
!      CALL NEWSE(A1,R1,SE,MS(M),IP,MJ0,MJ)
call newse(a1,r2,se,ms,l,ip,l+1)
400 continue
!     MIN.FPEC, MIN.RFPEC PRINT OUT
!      WRITE(6,607) OFPEC,ORFPEC,IFPEC
!      WRITE(6,1607) OAIC
! 1607 FORMAT(1H ,'MINIMUM AIC=',D12.5)
!     OSD, AO PRINT AND PUNCH OUT
!      WRITE(6,608)
!      CALL SUBMPR(OSD,IR,IR,MJ,MJ)
!  690 WRITE(7,1) N,IFPEC,IR,IL
!      DO 680 II=1,IR
!  680 WRITE(7,2) (OSD(II,JJ),JJ=1,IR)
!      IF(IFPEC.LE.0) GO TO 699
!      WRITE(6,609)
!      CALL PRMAT3(AO,IFPEC,IR,IP,0,MJ0,MJ,MJ)
!      DO 581 I=1,IFPEC
!      DO 582 II=1,IR
!  582 WRITE(7,2) (AO(I,II,JJ),JJ=1,IP)
!  581 CONTINUE
!  699 CONTINUE
!      CALL FLCLS3(NFL)
!  999 CONTINUE
!    1 FORMAT(10I5)
!    2 FORMAT(4D20.10)
!   39 FORMAT(1H ,50HPROGRAM 5.3.2   FPEC(AR-MODEL FITTING FOR CONTROL))
!   40 FORMAT(1H ,17HINITIAL CONDITION)
!   41 FORMAT(1H ,2HN=,I5,5X,2HL=,I5,5X,3HIR=,I5,5X,3HIL=,I5)
!   42 FORMAT(//1H ,17HCOVARIANCE MATRIX)
!  264 FORMAT(1H ,I5,2X,3D14.5)
!  259 FORMAT(/1H ,6HINW(I),5X,10I5)
!  600 FORMAT(///1H ,4X,1HI,12X,4HFPEC,9X,5HRFPEC,11X,3HAIC)
!  607 FORMAT(1H ,13HMINIMUM FPEC=,D12.5,2X,14HMINIMUM RFPEC=,D12.5,2X,14
!     AHATTAINED AT M=,I5)
!  608 FORMAT(//1H ,10X,10HOSD(II,JJ))
!  609 FORMAT(//1H ,10X,10H(A(I)B(I)))
return
end subroutine
!
!      SUBROUTINE SFPEC(SD,N,K,IR,MS,Z,RZ,OOZ,MJ,AIC)
subroutine sfpec(sd,n,k,ir,ms,z,rz,ooz,aic)
  use timsac_kinds, only: dp
  implicit none
!     FPEC COMPUTATION
!xx      IMPLICIT REAL*8(A-H,O-Z)
!      COMMON /COMA/AIC
!xx      DIMENSION SD(K,K)
!xx      DIMENSION SD1(IR,IR)
integer n, k, ir, ms
real(dp) sd(k,k), z, rz, ooz, aic
! local
integer i, j, km
real(dp) sd1(ir,ir), an, anp, anm, ap, apr, cst1, sdrm,&
&arm2
an=n
km=k*ms
anp=n+1+km
anm=n-1-km
ap=anp/anm
apr=ap**ir
cst1=1.0d-00
do 9 i=1,ir
!xx      DO 9 J=1,IR
do 8 j=1,ir
!xx    9 SD1(I,J)=SD(I,J)
sd1(i,j)=sd(i,j)
8 continue
9 continue
!      CALL SUBDET(SD1,SDRM,IR,MJ)
call subdet(sd1,sdrm,ir,ir)
z=apr*sdrm
arm2=2*ms*k*ir
aic=an*dlog(sdrm)+arm2
if(ms.ne.0) go to 10
ooz=cst1/z
10 rz=z*ooz
return
end subroutine
