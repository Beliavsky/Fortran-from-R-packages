! Part of the experimental modern Fortran port of timsac 1.3.8-6.
! This file was created or modified for the Fortran project on 2026-07-23.
! Original TIMSAC credits are retained; see NOTICE and ORIGIN.md.
! SPDX-License-Identifier: GPL-2.0-or-later

subroutine autcorf(x,n,cxx,cn,lagh1,xmean)
  use timsac_kinds, only: dp
  implicit none
!
!     PROGRAM 5.1.1   AUTO CORRELATION
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
!     LAGH: MAXIMUM LAG
!     DFORM: INPUT FORMAT SPECIFICATION STATEMENT IN ONE CARD,
!     FOR EXAMPLE
!     (8F10.4)
!     (X(I),I=1,N): ORIGINAL DATA.
!     THE OUTPUTS ARE AUTOCOVARIANCES (CXX(I); I=0,LAGH) AND
!     AUTO CORRELATIONS (NORMALIZED COVARIANCES).
!
!c      !DEC$ ATTRIBUTES DLLEXPORT :: AUTCORF
!
!xx      IMPLICIT REAL*8 (A-H,O-Z)
!      DIMENSION X(5000),CXX(1001),CN(1001)
!	 REAL*4 DFORM
!      DIMENSION  DFORM(20)
!xx      DIMENSION X(N),X1(N),CXX(LAGH1),CN(LAGH1)
integer n, lagh1
real(dp) x(n), cxx(lagh1), cn(lagh1), xmean
! local
integer i
real(dp) x1(n), cx0
!
!     INPUT / OUTPUT DATA FILE OPEN
!	CHARACTER(100) DFNAM
!	CALL SETWND
!	DFNAM='autcor.out'
!      CALL FLOPN3(DFNAM,NFL)
!      IF (NFL.EQ.0) GO TO 999
!     INITIAL CONDITION INPUT AND PRINT OUT
!      READ(5,1) N,LAGH
!      WRITE(6,50)
!      WRITE(6,51)
!      WRITE(6,52) N,LAGH
!     INPUT FORMAT SPECIFICATION
!      READ(5,4) (DFORM(I),I=1,20)
!    4 FORMAT(20A4)
!     ORIGINAL DATA INPUT AND PRINT OUT
!      READ(5,DFORM) (X(I),I=1,N)
!      WRITE(6,53)
!      WRITE(6,54)
!      DO 220 I=1,N
!  220 WRITE(6,55) I,X(I)
do 220 i=1,n
!xx  220 X1(I)=X(I)
x1(i)=x(i)
220 continue
!     MEAN DELETION
call dmeadl(x1,n,xmean)
!     AUTO COVARIANCE COMPUTATION
!      LAGH1=LAGH+1
call crosco(x1,x1,n,cxx,lagh1)
!     NORMALIZATION
cx0=cxx(1)
call cornom(cxx,cn,lagh1,cx0,cx0)
!     AUTO COVARIANCE PRINT OUT
!      WRITE(6,162) N,LAGH,XMEAN
!      WRITE(6,163)
!      CALL PRCOL2(CXX,CN,1,LAGH1,1)
!     AUTO COVARIANCE PUNCH OUT
!      WRITE(7,1) N,LAGH
!      WRITE(7,2) (CXX(I),I=1,LAGH1)
!	CALL FLCLS3(NFL)
!  999 CONTINUE
!    1 FORMAT(10I5)
!    2 FORMAT(4D20.10)
!   50 FORMAT(1H ,13HPROGRAM 5.1.1,3X,16HAUTO CORRELATION)
!   51 FORMAT(1H ,17HINITIAL CONDITION)
!   52 FORMAT(1H ,2HN=,I5,5X,5HLAGH=,I5)
!   53 FORMAT(1H ,13HORIGINAL DATA)
!   54 FORMAT(1H ,4X,1HI,6X,4HX(I))
!   55 FORMAT(1H ,I5,2X,F10.4)
!CC  162 FORMAT(//1H ,14HAUTOCOVARIANCE,5X,6HCXX(I),5X,2HN=,I5,5X,5HLAGH=,I
!  162 FORMAT(1H ,14HAUTOCOVARIANCE,5X,6HCXX(I),5X,2HN=,I5,5X,5HLAGH=,I
!     A5,5X,5HMEAN=,D15.5)
!  163 FORMAT(1H ,4X,1HI,5X,6HCXX(I),8X,10HNORMALIZED)
return
end
