! Part of the experimental modern Fortran port of timsac 1.3.8-6.
! This file was created or modified for the Fortran project on 2026-07-23.
! Original TIMSAC credits are retained; see NOTICE and ORIGIN.md.
! SPDX-License-Identifier: GPL-2.0-or-later

subroutine covgenf(l,k,f,g,c,cn)
  use timsac_kinds, only: dp
  implicit none
!
!c	PROGRAM COVGEN
!     PROGRAM 74.1.3. COVARIANCE GENERATION FROM GAIN FUNCTION.
!-----------------------------------------------------------------------
!     ** DESIGNED BY H. AKAIKE, THE INSTITUTE OF STATISTICAL MATHEMATICS
!     ** PROGRAMMED BY E. ARAHATA, THE INSTITUTE OF STATISTICAL MATHEMAT
!         TOKYO
!     ** DATE OF THE LATEST REVISION: MARCH 25, 1977
!     ** THIS PROGRAM WAS ORIGINALLY PUBLISHED IN
!         "TIMSAC-74 A TIME SERIES ANALYSIS AND CONTROL PROGRAM PACKAGE(1
!         BY H. AKAIKE, E. ARAHATA AND T. OZAKI, COMPUTER SCIENCE MONOGRA
!         NO.5, MARCH 1975, THE INSTITUTE OF STATISTICAL MATHEMATICS
!     ** FOR THE BASIC THEORY SEE "CANONICAL CORRELATION ANALYSIS OF TIM
!         AND THE USE OF AN INFORMATION CRITERION" BY H. AKAIKE, IN
!         "SYSTEM IDENTIFICATION: ADVANCES AND CASE STUDIES" R. K. MEHRA
!         D. G. LAINIOTIS EDS. ACADEMIC PRESS, NEW YORK, 1976
!-----------------------------------------------------------------------
!     THIS PROGRAM PRODUCES THE FOURIER TRANSFORM OF A POWER
!     GAIN FUNCTION IN THE FORM OF AN AUTOCOVARIANCE SEQUENCE.
!     THE GAIN FUNCTION IS DEFINED AS A RECTILINEAR FUNCTION WITH
!     THE VALUES G(I) SPECIFIED AT THE FREQUENCIES F(I),I=1,K.
!     THE OUTPUTS OF THIS PROGRAM ARE USED AS THE INPUTS TO THE CANONICA
!     CORRELATION ANALYSIS PROGRAM CANARM, TO REALIZE A FILTER WITH
!     THE DESIRED GAIN FUNCTION.
!
!     THE FOLLOWING INPUTS ARE REQUIRED:
!     (L,K): L, DESIRED MAXIMUM LAG OF COVARIANCE (AT MOST 1024)
!            K, NUMBER OF DATA POINTS (LESS THAN OR EQUAL TO 500)
!     (F(I),G(I))(I=1,K): F(I), FREQUENCY. BY DEFINITION F(1)=0.0 AND F(
!                           F(I)'S ARE ARRANGED IN INCREASING ORDER.
!                           G(I), POWER GAIN OF THE FILTER AT THE FREQUEN
!
!     OUTPUTS:
!     (N,LAGH): N=2048
!     C(I)(I=0,LAGH):
!
!xx      PARAMETER (LMAX=1024)
integer, parameter :: lmax=1024
!xx      IMPLICIT REAL*8 (O-Y)
!c      REAL*8 C(1025),CN(1025),CX0
!c      COMPLEX*16 A(2048)
!c      DIMENSION F(500),G(500)
!c      DATA F / 500*0.0 /, G / 500*0.0 /
!xx      REAL*8  C(L+1),CN(L+1),CX0
!xx      COMPLEX*16 A(LMAX*2)
!xx      REAL*8  F(K),G(K)
integer l, k
real(dp) f(k), g(k), c(l+1), cn(l+1)
! local
integer i, i1, i2, j, j1, lagh1, n, n2p, nn
real(dp) cx0, q00, an, bn, cm, q1, q2, q3, q4, q25
complex(kind(0.0d0)) a(lmax*2)
!
!     INPUT / OUTPUT DATA FILE OPEN
!c      CHARACTER(100) DFNAM
!c      CALL SETWND
!c      DFNAM='covgen.out'
!c      CALL FLOPN3(DFNAM,NFL)
!c      IF (NFL.EQ.0) GO TO 999
!
q00=0.0d-00
!     INITIAL CONDITION INPUT
!c      READ(5,800) L,K
!c      WRITE(6,900) L,K
!c      DO  100 I=1,K
!c  100 READ(5,801) F(I),G(I)
!c      DO  101 I=1,K
!c  101 WRITE(6,901) F(I),G(I)
!c  105 N=2048
n=lmax*2
n2p=11
!c      AN=2048.0
!c      NN=1023
an=n
nn=lmax-1
do  200 i=1,nn
bn=i
cm=bn/an
j=k
!c  106 IF(CM-F(J)) 107,108,108
106 continue
if(cm-f(j) .ge. 0) go to 108
!xx  107 J=J-1
j=j-1
go to 106
108 j1=j+1
q1=(f(j1)-cm)*g(j)
q2=(cm-f(j))*g(j1)
q3=f(j1)-f(j)
q4=(q1+q2)/q3
i1=i+1
q1=q4
!xx      A(I1)=DCMPLX(Q1,Q00)
a(i1)=cmplx(q1,q00,kind(0.0d0))
i2=n-i+1
!xx  200 A(I2)=A(I1)
a(i2)=a(i1)
200 continue
q25=g(1)
!xx      A(1)=DCMPLX(Q25,Q00)
a(1)=cmplx(q25,q00,kind(0.0d0))
q25=g(k)
!xx      A(NN+2)=DCMPLX(Q25,Q00)
a(nn+2)=cmplx(q25,q00,kind(0.0d0))
!     FAST FOURIER TRANSFORM OF A
!     COMMON SUBROUTINE CALL
call mixrad(a,n,n2p,+1)
lagh1=l+1
do  201 i=1,lagh1
!xx  201 C(I)=DREAL(A(I))
!xx      C(I)=DREAL(A(I)) 
c(i)=real(a(i))
201 continue
!     C(I+1),I=0,L     ARE THE DESIRED AUTOCOVARIANCES.
!     NORMARIZED COVARIANCES ARE GIVEN AS CN(I+1),I=0,L.
cx0=c(1)
!     COMMON SUBROUTINE CALL
call cornom(c,cn,lagh1,cx0,cx0)
!c      WRITE(6,902)
!c      DO  210 I=1,LAGH1
!c      J=I-1
!c  210 WRITE(6,903) J,C(I),CN(I)
!c      LAGH=LAGH1-1
!c      WRITE(7,800) N,LAGH
!c      WRITE(7,904) (C(I),I=1,LAGH1)
!c      CALL FLCLS3(NFL)
return
!xx  800 FORMAT(2I5)
!xx  801 FORMAT(2F10.5)
!xx  900 FORMAT(1H ,' PROGRAM 74.1.3. COVGEN: L=',I5,5X,'K=',I5/1H )
!xx  901 FORMAT(1H ,' FREQUENCY(',F10.5,')',20X,'DESIRED GAIN = ',F10.5)
!xx  902 FORMAT(//1H ,20X,'AUTO COVARIANCE',25X,
!xx     *'AUTO COVARIANCE NORMALIZED')
!xx  903 FORMAT(1H ,'LAG =',I5,5X,D20.10,31X,D20.10)
!xx  904 FORMAT(4D20.10)
end
