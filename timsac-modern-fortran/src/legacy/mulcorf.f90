! Part of the experimental modern Fortran port of timsac 1.3.8-6.
! This file was created or modified for the Fortran project on 2026-07-23.
! Original TIMSAC credits are retained; see NOTICE and ORIGIN.md.
! SPDX-License-Identifier: GPL-2.0-or-later

subroutine mulcorf(x1,n,k,lagh1,sm,c,cn)
  use timsac_kinds, only: dp
  implicit none
!
!     PROGRAM 5.1.2   MULTIPLE CORRELATION
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
!     THIS PROGRAM REQUIRES FOLLOWING INPUTS:
!     N: LENGTH OF DATA
!     K: DIMENSION OF THE OBSERVATION VECTOR
!     LAGH: MAXIMUM LAG
!     ISW: ISW=1...ROWWISE DATA INPUT
!           ISW=2...COLUMNWISE DATA INPUT
!     DFORM: INPUT FORMAT SPECIFICATION STATEMENT IN ONE CARD,
!     FOR EXAMPLE
!     (8F10.4)
!     (X1(S,I); S=1,...,N, I=1,...,K): ORIGINAL DATA MATRIX.
!     THE OUTPUTS ARE (CIJ(L): L=0,1,...,LAGH) (I=1,...,K; J=1,...,K),
!     WHERE CIJ(L)=COVARIANCE(XI(S+L),XJ(S)),
!     AND THEIR NORMALIZED (CORRELATION) VALUES.
!
!xx      IMPLICIT REAL*8 (A-H,O-Z)
!      DIMENSION X1(2000,10)
!      DIMENSION X(2000),Y(2000)
!      DIMENSION C1(501),C2(501),CN1(501),CN2(501)
!      DIMENSION SM(10),C0(10)
!      REAL*4 DFORM
!      DIMENSION DFORM(20)
!xx      DIMENSION X1(N,K),X2(N,K)
!xx      DIMENSION X(N),Y(N)
!xx      DIMENSION C(LAGH1,K,K),CN(LAGH1,K,K)
!xx      DIMENSION C1(LAGH1),C2(LAGH1),CN1(LAGH1),CN2(LAGH1)
!xx      DIMENSION SM(K),C0(K)
integer n, k, lagh1
real(dp) x1(n,k), sm(k), c(lagh1,k,k), cn(lagh1,k,k)
! local
integer i, ii, im1, jj
real(dp) x2(n,k), x(n), y(n), c1(lagh1), c2(lagh1),&
&cn1(lagh1), cn2(lagh1), c0(k), cx0, cy0, xmean
!
!     INPUT / OUTPUT DATA FILE OPEN
!      CHARACTER(100) DFNAM
!      DFNAM='mulcor.out'
!      CALL SETWND
!      CALL FLOPN3(DFNAM,NFL)
!      IF (NFL.EQ.0) GO TO 999
!     INITIAL CONDITION INPUT AND OUTPUT
!      READ(5,1) N,LAGH,K,ISW
!      LAGH1=LAGH+1
!      WRITE(6,50)
!      WRITE(6,51)
!      WRITE(6,52) N,LAGH,K
!     INITIAL CONDITION PUNCH OUT
!      WRITE(7,1) N,LAGH,K
!     INPUT FORMAT SPECIFICATION
!      READ(5,4) (DFORM(I),I=1,20)
!    4 FORMAT(20A4)
!     ORIGINAL DATA INPUT AND OUTPUT
!      GO TO(8,9),ISW
!    8 DO 208 I=1,N
!  208 READ(5,DFORM) (X(I,II),II=1,K)
!      GO TO 400
!    9 DO 209 II=1,K
!  209 READ(5,DFORM) (X(I,II),I=1,N)
!  400 WRITE(6,53)
!      WRITE(6,54)
!      WRITE(6,154)
!      DO 220 I=1,N
!  220 WRITE(6,55) I,(X(I,II),II=1,K)
!
!     MEAN DELETION
do 300 ii=1,k
do 310 i=1,n
!xx  310  X(I)=X1(I,II)
x(i)=x1(i,ii)
310 continue
call dmeadl(x,n,xmean)
sm(ii)=xmean
do 320 i=1,n
!xx  320   X2(I,II)=X(I)
x2(i,ii)=x(i)
320 continue
300 continue
!
!     COVARIANCE COMPUTATION
do 10 ii=1,k
do 110 i=1,n
!  110 X(I)=X1(I,II)
!xx  110 X(I)=X2(I,II)
x(i)=x2(i,ii)
110 continue
!     AUTO COVARIANCE COMPUTATION
call crosco(x,x,n,c1,lagh1)
!     NORMALIZATION
c0(ii)=c1(1)
cx0=c0(ii)
call cornom(c1,cn1,lagh1,cx0,cx0)
!     AUTO COVARIANCE PRINT OUT
!      WRITE(6,162) II,II,SM(II)
!      WRITE(6,163)
!      CALL PRCOL2(C1,CN1,1,LAGH1,1)
!     AUTO COVARIANCE PUNCH OUT
!      WRITE(7,1) II,II
!      WRITE(7,2) (C1(I),I=1,LAGH1)
do 115 i=1,lagh1
c(i,ii,ii)=c1(i)
cn(i,ii,ii)=cn1(i)
115 continue
if(ii.eq.1) go to 10
im1=ii-1
do 11 jj=1,im1
do 120 i=1,n
!  120 Y(I)=X1(I,JJ)
!xx  120 Y(I)=X2(I,JJ)
y(i)=x2(i,jj)
120 continue
!     CROSS COVARIANCE COMPUTATION
call crosco(x,y,n,c1,lagh1)
call crosco(y,x,n,c2,lagh1)
! 　　NORMALIZATION
cx0=c0(ii)
cy0=c0(jj)
call cornom(c1,cn1,lagh1,cx0,cy0)
call cornom(c2,cn2,lagh1,cx0,cy0)
!     CROSS COVARIANCE PRINT OUT
!      WRITE(6,165) II,JJ
!      WRITE(6,166)
!      CALL PRCOL4(C1,CN1,C2,CN2,1,LAGH1,1)
!     CROSS COVARIANCE PUNCH OUT
!      WRITE(7,1) II,JJ
!      WRITE(7,2) (C1(I),I=1,LAGH1)
!      WRITE(7,1) JJ,II
!      WRITE(7,2) (C2(I),I=1,LAGH1)
do 125 i=1,lagh1
c(i,ii,jj)=c1(i)
c(i,jj,ii)=c2(i)
cn(i,ii,jj)=cn1(i)
cn(i,jj,ii)=cn2(i)
125 continue
11 continue
10 continue
!      CALL FLCLS3(NFL)
!  999 CONTINUE
!    1 FORMAT(10I5)
!    2 FORMAT(4D20.10)
!   50 FORMAT(1H ,36HPROGRAM 5.1.2   MULTIPLE CORRELATION)
!   51 FORMAT(1H ,17HINITIAL CONDITION)
!   52 FORMAT(1H ,2HN=,I5,5X,5HLAGH=,I5,5X,2HK=,I5)
!   53 FORMAT(1H ,13HORIGINAL DATA)
!   54 FORMAT(1H ,4X,1HI,2X,8HX1(I,II))
!  154 FORMAT(1H ,16X,91H1	   2	     3	       4	 5
!     A	 6	   7	     8	       9	10)
!   55 FORMAT(1H ,I5,2X,10F10.4)
!  162 FORMAT(//1H ,14HAUTOCOVARIANCE,5X,6HCIJ(L),5X,2HI=,I5,5X,2HJ=,I5,5
!     AX,5HMEAN=,D15.5)
!  163 FORMAT(1H ,4X,1HL,5X,6HCIJ(L),8X,10HNORMALIZED)
!  165 FORMAT(//1H ,16HCROSS COVARIANCE,5X,6HCIJ(L),5X,2HI=,I5,5X,2HJ=,I5
!     A)
!  166 FORMAT(1H ,4X,1HL,5X,6HCIJ(L),8X,10HNORMALIZED,4X,6HCJI(L),8X,10HN
!     AORMALIZED)
return
end subroutine
