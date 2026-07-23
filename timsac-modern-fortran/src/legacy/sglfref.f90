! Part of the experimental modern Fortran port of timsac 1.3.8-6.
! This file was created or modified for the Fortran project on 2026-07-23.
! Original TIMSAC credits are retained; see NOTICE and ORIGIN.md.
! SPDX-License-Identifier: GPL-2.0-or-later

subroutine sglfref(inp,iout,n,lagh1,ip,p,p11,p22,c,s,cp11,cp22,&
&cc,cs,r,ph)
  use timsac_kinds, only: dp
  implicit none
!
!     PROGRAM 5.2.3   FREQUENCY RESPONSE FUNCTION (SINGLE CHANNEL)
!-----------------------------------------------------------------------
!      SUBROUTINE SGLARC(C,S,ARC,LAGH1)
!      SUBROUTINE SGLERR(CH,R,N,LAGH1)
!      SUBROUTINE SGLPAC(ARC,PH,LAGH1)
!      SUBROUTINE SPHASE(C,S,PH,LAGH1)
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
!     THIS PROGRAM COMPUTES 1-INPUT,1-OUTPUT FREQUECNY RESPONSE FUNCTION
!     ,GAIN,PHASE,COHERENCY AND RELATIVE ERROR STATISTICS.
!     ONE CARD WITH SPECIFICATION OF INPUT(INP) AND OUTPUT(IOUT)
!     VARIABLE SHOULD BE ADDED ON TOP OF THE OUTPUT OF PROGRAM 5.2.2
!     MULSPE TO FORM THE INPUT TO THIS PROGRAM.
!     WITHIN IP VARIABLES OF MULSPE OUTPUT, INP-TH AND IOUT-TH VARIABLE
!     ARE TAKEN AS INPUT AND OUTPUT VARIABLE.
!
!xx      IMPLICIT REAL*8(A-H,O-Z)
!      DIMENSION P11(501),P22(501),C(501),S(501)
!      DIMENSION R(501),PH(501),P(10,10)
!xx      DIMENSION P11(LAGH1),P22(LAGH1),C(LAGH1),S(LAGH1)
!xx      DIMENSION CP11(LAGH1),CP22(LAGH1),CC(LAGH1),CS(LAGH1)
!xx      DIMENSION R(LAGH1),PH(LAGH1),P(LAGH1,IP,IP)
integer inp, iout, n, lagh1, ip
real(dp) p(lagh1,ip,ip), p11(lagh1), p22(lagh1), c(lagh1),&
&s(lagh1), cp11(lagh1), cp22(lagh1), cc(lagh1),&
&cs(lagh1), r(lagh1), ph(lagh1)
! local
integer i
!     INPUT / OUTPUT DATA FILE OPEN
!      CALL SETWND
!      CALL FLOPN2(NFL)
!      IF (NFL.EQ.0) GO TO 999
!     ABSOLUTE DIMENSION USED FOR SUBROUTINE CALL
!      MJ=10
!     INPUT OUTPUT VARIABLE SPECIFICATION
!      READ(5,1) INP,IOUT
!     READING THE OUTPUT OF PROGRAM 5.2.2 MULSPE
!      READ(5,1) N,LAGH,IP
!      LAGH1=LAGH+1
do 5 i=1,lagh1
!      CALL REMATX(P,IP,IP,1,MJ,MJ)
!     MATRIX REARRANGEMENT
!      P11(I)=P(INP,INP)
!      P22(I)=P(IOUT,IOUT)
p11(i)=p(i,inp,inp)
p22(i)=p(i,iout,iout)
if(inp.lt.iout) go to 7
!      C(I)=P(INP,IOUT)
!      S(I)=-P(IOUT,INP)
c(i)=p(i,inp,iout)
s(i)=-p(i,iout,inp)
go to 5
!    7 C(I)=P(IOUT,INP)
!      S(I)=P(INP,IOUT)
7 c(i)=p(i,iout,inp)
s(i)=p(i,inp,iout)
5 continue
!     INITIAL CONDITION PRINT OUT
!      WRITE(6,55)
!      WRITE(6,56)
!      WRITE(6,57) N,LAGH
!      WRITE(6,57) N,LAGH
!      WRITE(6,259) INP,IOUT
!      WRITE(6,58)
!      WRITE(6,59)
!      WRITE(6,159)
!      CALL PRCOL4(P11,P22,C,S,1,LAGH1,1)
!     FREQUENCY RESPONSE FUNCTION COMPUTATION
do 10 i=1,lagh1
!      C(I)=C(I)/P11(I)
!      S(I)=S(I)/P11(I)
!   10 P22(I)=P22(I)/P11(I)
cc(i)=c(i)/p11(i)
cs(i)=s(i)/p11(i)
!xx   10 CP22(I)=P22(I)/P11(I)
cp22(i)=p22(i)/p11(i)
10 continue
!     GAIN COMPUTATION
do 11 i=1,lagh1
!      R(I)=C(I)**2+S(I)**2
!   11 P11(I)=DSQRT(R(I))
r(i)=cc(i)**2+cs(i)**2
!xx   11 CP11(I)=DSQRT(R(I))
cp11(i)=dsqrt(r(i))
11 continue
!     PHASE COMPUTATION
!      CALL SPHASE(C,S,PH,LAGH1)
call sphase(cc,cs,ph,lagh1)
!    COHERENCY COMPUTATION
do 12 i=1,lagh1
!   12 P22(I)=R(I)/P22(I)
!xx   12 CP22(I)=R(I)/CP22(I)
cp22(i)=r(i)/cp22(i)
12 continue
!     RELATIVE ERROR STATISTICS COMPUTATION
!      CALL SGLERR(P22,R,N,LAGH1)
call sglerr(cp22,r,n,lagh1)
!     FREQUENCY RESPONSE FUNCTION, GAIN, PHASE, COHERENCY AND RELATIVE
!     ERROR STATISTICS PRINT OUT
!      WRITE(6,60)
!      WRITE(6,61)
!      CALL PRCOL6(C,S,P11,PH,P22,R,1,LAGH1,1)
!	CALL FLCLS2(NFL)
!  999 CONTINUE
!    1 FORMAT(10I5)
!   55 FORMAT(1H ,60HPROGRAM 5.2.3   FREQUENCY RESPONSE FUNCTION (SINGLE
!     ACHANNEL))
!   56 FORMAT(1H ,17HINITIAL CONDITION)
!   57 FORMAT(1H ,2HN=,I5,5X,5HLAGH=,I5)
!   58 FORMAT(1H ,22HINITIAL DATA(SPECTRUM))
!   59 FORMAT(1H ,6X,14HPOWER SPECTRUM,1X,14HPOWER SPECTRUM,3X,11HCO-SPEC
!     ATRUM,1X,13HQUAD-SPECTRUM)
!  159 FORMAT(1H ,4X,1HI,10X,6HP(1,1),8X,6HP(2,2),10X,4HC(I),10X,4HS(I))
!   60 FORMAT(//1H ,4X,1HI,3X,27HFREQUENCY RESPONSE FUNCTION,10X,4HGAIN,9
!     AX,5HPHASE,5X,9HCOHERENCY,6X,8HRELATIVE)
!   61 FORMAT(1H ,12X,9HREAL PART,4X,10HIMAG. PART,51X,5HERROR)
!  259 FORMAT(1H ,6HINPUT=,I5,5X,7HOUTPUT=,I5)
return
end subroutine
!
subroutine sglarc(c,s,arc,lagh1)
  use timsac_kinds, only: dp
  implicit none
!     THIS SUBROUTINE COMPUTES RAW PHASES.
!     (SINGLE CHANNEL)
!xx      IMPLICIT REAL*8(A-H,O-Z)
!xx      DIMENSION C(LAGH1),S(LAGH1),ARC(LAGH1)
integer lagh1
real(dp) c(lagh1), s(lagh1), arc(lagh1)
! local
integer i
real(dp) pi, cst5
pi=3.1415926536
cst5=0.5d-00
do 10 i=1,lagh1
!c      IF(C(I)) 11,12,13
!c   11 IF(S(I)) 14,15,16
!c   12 IF(S(I)) 17,18,19
if(c(i).eq.0) go to 12
if(c(i).gt.0) go to 13
!xx   11 IF(S(I).LT.0) GO TO 14
if(s(i).lt.0) go to 14
if(s(i).eq.0) go to 15
if(s(i).gt.0) go to 16
12 if(s(i).lt.0) go to 17
if(s(i).eq.0) go to 18
if(s(i).gt.0) go to 19
13 arc(i)=datan(s(i)/c(i))
go to 10
14 arc(i)=datan(s(i)/c(i))-pi
go to 10
15 arc(i)=-pi
go to 10
16 arc(i)=datan(s(i)/c(i))+pi
go to 10
17 arc(i)=-pi*cst5
go to 10
18 arc(i)=0.0d-00
go to 10
19 arc(i)=pi*cst5
10 continue
return
end
!
subroutine sglerr(ch,r,n,lagh1)
  use timsac_kinds, only: dp
  implicit none
!     THIS SUBROUTINE COMPUTES RELATIVE ERROR STATISTICS.
!     (SINGLE CHANNEL)
!xx      IMPLICIT REAL*8(A-H,O-Z)
!xx      DIMENSION CH(LAGH1),R(LAGH1)
integer n, lagh1
real(dp) ch(lagh1), r(lagh1)
! local
integer i, lagh
real(dp) d1, d2, cst0, cst1, cst100, e1, er
!     CONSTANTS D1,D2 COMPUTATION
lagh=lagh1-1
call subd12(n,lagh,1,d1,d2)
!     RELATIVE ERROR STATISTICS COMPUTATION
cst0=0.0d-00
cst1=1.0d-00
cst100=100.0d-00
do 20 i=1,lagh1
if(ch(i).le.cst0) go to 22
if(ch(i).gt.cst1) go to 22
e1=cst1/ch(i)-cst1
er=dsqrt(e1)
if(i.eq.1) go to 23
if(i.eq.lagh1) go to 23
r(i)=d2*er
go to 20
23 r(i)=d1*er
go to 20
22 r(i)=cst100
20 continue
return
end
!
subroutine sglpac(arc,ph,lagh1)
  use timsac_kinds, only: dp
  implicit none
!     THIS SUBROUTINE MAKES PHASE CURVE CONTINUOUS.
!     (SINGLE CHANNEL)
!xx      IMPLICIT REAL*8(A-H,O-Z)
!xx      DIMENSION ARC(LAGH1),PH(LAGH1)
integer lagh1
real(dp) arc(lagh1), ph(lagh1)
! local
integer i
real(dp) pi, pi2, dk
pi=3.1415926536
pi2=pi+pi
ph(1)=arc(1)
do 10 i=2,lagh1
dk=arc(i)-arc(i-1)
if(dk.gt.pi) go to 11
if(dk.lt.-pi) go to 12
ph(i)=ph(i-1)+dk
go to 10
11 ph(i)=ph(i-1)+dk-pi2
go to 10
12 ph(i)=ph(i-1)+dk+pi2
10 continue
return
end
!
subroutine sphase(c,s,ph,lagh1)
  use timsac_kinds, only: dp
  implicit none
!     THIS SUBROUTINE COMPUTES PHASE.
!     (SINGLE CHANNEL)
!xx      IMPLICIT REAL*8(A-H,O-Z)
!xx      DIMENSION C(LAGH1),S(LAGH1),PH(LAGH1)
!      DIMENSION ARC(501)
!xx      DIMENSION ARC(LAGH1)
integer lagh1
real(dp) c(lagh1), s(lagh1), ph(lagh1)
! local
real(dp) arc(lagh1)
!     ARCTANGENT COMPUTATION
call sglarc(c,s,arc,lagh1)
!     PHASE COMPUTATION
call sglpac(arc,ph,lagh1)
return
end
