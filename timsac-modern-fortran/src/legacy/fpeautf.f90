! Part of the experimental modern Fortran port of timsac 1.3.8-6.
! This file was created or modified for the Fortran project on 2026-07-23.
! Original TIMSAC credits are retained; see NOTICE and ORIGIN.md.
! SPDX-License-Identifier: GPL-2.0-or-later

subroutine fpeautf(l,n,sd,cxx,ssd,fpe,rfpe,d,chi2,&
&ofpe1,ofpe2,orfpe,mo,osd,a,ao)
  use timsac_kinds, only: dp
  implicit none
!
!     PROGRAM 5.3.1   FPE AUTO
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
!     THIS PROGRAM PERFORMS FPE(FINAL PREDICTION ERROR) COMPUTATION FOR
!     ONE-DIMENSIONAL AR-MODEL. A CARD CONTAINING THE FOLLOWING
!     INFORMATION OF L, UPPER LIMIT OF MODEL ORDER, SHOULD BE ADDED ON
!     TOP OF THE OUTPUT OF PROGRAM 5.1.1 AUTCOR TO FORM THE INPUT TO
!     THIS PROGRAM.
!     CXX(0) IS READ AS INITIAL SD.
!     THE OUTPUTS ARE THE COEFFICIENTS A(I) OF AR-PROCESS
!     X(N)=A(1)X(N-1)+...+A(M)X(N-M)+E(N)
!     AND THE VARIANCE SIGMA**2 OF E(N).
!     CHI**2 SHOWS THE SIGNIFICANCE OF PARCOR=A(M) AS A CHI-SQUARED
!     VARIABLE WITH D.F.=1.
!
!xx      IMPLICIT REAL*8(A-H,O-Z)
!      DIMENSION CXX(501),A(501),B(501),AO(501)
!xx      DIMENSION CXX(L),A(L,L),B(L),AO(L)
!xx      DIMENSION SSD(L),FPE(L),RFPE(L),D(L),CHI2(L)
integer l, n, mo
real(dp) sd, cxx(l), ssd(l), fpe(l), rfpe(l), d(l),&
&chi2(l),ofpe1, ofpe2, orfpe, osd, a(l,l), ao(l)
! local
integer i, im, lm, m, mp1, nm1, np1
real(dp) b(l), an, anp1, anm1, cst1, oofpe, se, sd0, d2
!
!     INPUT / OUTPUT DATA FILE OPEN
!      CHARACTER(100) DFNAM
!      CALL SETWND
!      DFNAM='fpeaut.out'
!      CALL FLOPN3(DFNAM,NFL)
!      IF (NFL.EQ.0) GO TO 999
!     L SPECIFICATION
!      READ(5,1) L
!     READING THE OUTPUT OF PROGRAM 5.1.1 AUTCOR
!      READ(5,1) N,LAGH
!      READ(5,2) SD,(CXX(I),I=1,LAGH)
!
!     COMPUTATION START
an=n
np1=n+1
nm1=n-1
anp1=np1
anm1=nm1
!      OFPE=(ANP1/ANM1)*SD
ofpe1=(anp1/anm1)*sd
cst1=1.0d-00
!      OOFPE=CST1/OFPE
oofpe=cst1/ofpe1
orfpe=cst1
osd=sd
mo=0
ofpe2=ofpe1
!      WRITE(6,155)
!      WRITE(6,156)
!      WRITE(6,57) N,L
!      WRITE(6,140)
!      WRITE(6,141) SD
!      CALL PRCOL1(CXX,1,L,0)
!      WRITE(6,157)
!      WRITE(6,58) OFPE
se=cxx(1)
sd0=sd
!
do 400 m=1,l
mp1=m+1
!      D=SE/SD
!      A(M)=D
!      D2=D*D
!      SD=(CST1-D2)*SD
d(m)=se/sd0
a(m,m)=d(m)
d2=d(m)*d(m)
ssd(m)=(cst1-d2)*sd0
sd0=ssd(m)
anp1=np1+m
anm1=nm1-m
!      FPE=(ANP1/ANM1)*SD
!      RFPE=FPE*OOFPE
!      CHI2=D2*ANM1
fpe(m)=(anp1/anm1)*ssd(m)
rfpe(m)=fpe(m)*oofpe
chi2(m)=d2*anm1
if(m.eq.1) go to 410
!     A(I) COMPUTATION
lm=m-1
do 420 i=1,lm
!  420 A(I)=A(I)-D*B(I)
a(i,m)=a(i,m-1)-d(m)*b(i)
420 continue
!  410 DO 421 I=1,M
410 continue
do 421 i=1,m
im=mp1-i
!  421 B(I)=A(IM)
b(i)=a(im,m)
421 continue
!      WRITE(6,60) M
!      WRITE(6,61) SD,FPE,RFPE
!      WRITE(6,62) D,CHI2
!      WRITE(6,160)
!      CALL PRCOL1(A,1,M,0)
!
!      IF(OFPE.LT.FPE) GO TO 440
!      OFPE=FPE
!      ORFPE=RFPE
!      OSD=SD
if(ofpe2.lt.fpe(m)) go to 440
ofpe2=fpe(m)
orfpe=rfpe(m)
osd=ssd(m)
mo=m
do 430 i=1,m
!  430 AO(I)=A(I)
!xx  430 AO(I)=A(I,M)
ao(i)=a(i,m)
430 continue
440 if(m.eq.l) go to 400
se=cxx(mp1)
do 441 i=1,m
!xx  441 SE=SE-B(I)*CXX(I)
se=se-b(i)*cxx(i)
441 continue
400 continue
!
!      WRITE(6,63) OFPE,MO
!      WRITE(6,64) ORFPE
!      WRITE(7,1) N,MO
!      WRITE(7,2) OSD
!      IF(MO.LE.0) GO TO 699
!      WRITE(7,2) (AO(I),I=1,MO)
!  699 CONTINUE
!      CALL FLCLS3(NFL)
!    1 FORMAT(10I5)
!    2 FORMAT(4D20.10)
!   57 FORMAT(1H ,2HN=,I5,5X,2HL=,I5)
!   58 FORMAT(1H ,5HOFPE=,D12.5)
!   60 FORMAT(1H ,2HM=,I5)
!   61 FORMAT(1H ,9HSIGMA**2=,D12.5,2X,4HFPE=,D12.5,2X,5HRFPE=,D12.5)
!   62 FORMAT(1H ,7HPARCOR=,D14.5,2X,15HCIH**2(D.F.=1)=,D12.5)
!  160 FORMAT(1H ,4X,1HI,12X,4HA(I))
!   63 FORMAT(1H ,13HMINIMUM FPE =,D12.5,2X,14HATTAINED AT M=,I5)
!   64 FORMAT(1H ,13HMINIMUM RFPE=,D12.5)
!  140 FORMAT(1H ,4X,1HI,5X,15HAUTO COVARIANCE)
!  141 FORMAT(1H ,4X,1H0,D16.5)
!  155 FORMAT(1H ,24HPROGRAM 5.3.1   FPE AUTO)
!  156 FORMAT(1H ,17HINITIAL CONDITION)
!  157 FORMAT(1H ,2HM=,4X,1H0)
return
end subroutine
