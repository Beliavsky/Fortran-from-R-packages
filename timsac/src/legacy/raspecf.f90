! Part of the experimental modern Fortran port of timsac 1.3.8-6.
! This file was created or modified for the Fortran project on 2026-07-23.
! Original TIMSAC credits are retained; see NOTICE and ORIGIN.md.
! SPDX-License-Identifier: GPL-2.0-or-later

subroutine raspecf(h,l,k,sgme2,a,b,pxx)
  use timsac_kinds, only: dp
  implicit none
!
!     PROGRAM 5.4.1   RATIONAL SPECTRUM
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
!     THIS PROGRAM COMPUTES POWER SPECTRUM OF AR-MA PROCESS
!     X(N)=A(1)X(N-1)+...+A(L)X(N-L)+E(N)+B(1)E(N-1)+...+B(K)E(N-K)
!     WHERE E(N) IS A WHITE NOISE WITH ZERO MEAN AND VARIANCE EQUAL TO
!     SGME2.  OUTPUTS PXX(I) ARE GIVEN AT FREQUENCIES I/(2*H)
!     I=0,1,...,H.
!     REQUIRED INPUTS ARE:
!     L,K,H,SGME2,(A(I),I=1,L), AND (B(I),I=1,K).
!     0 IS ALLOWABLE AS L AND/OR K.
!
!xx      IMPLICIT REAL*8(A-H,O-Z)
!xx      INTEGER H,H1
!C      DIMENSION A(501),B(501)
!C      DIMENSION G(501),GR1(501),GI1(501),GR2(501),GI2(501)
!C      DIMENSION PXX(510)
!xx      DIMENSION A(L),B(K)
!C      DIMENSION G(MAX(L,K)+1),GR1(H+1),GI1(H+1),GR2(H+1),GI2(H+1)
!xx      DIMENSION G(L+K+1),GR1(H+1),GI1(H+1),GR2(H+1),GI2(H+1)
!xx      DIMENSION PXX(H+1)
integer h, l, k
real(dp) sgme2, a(l), b(k), pxx(h+1)
! local
integer i, i1, k1, l1, h1
real(dp) g(l+k+1), gr1(h+1), gi1(h+1), gr2(h+1), gi2(h+1)
!
!     INPUT / OUTPUT DATA FILE OPEN
!C      CALL SETWND
!C      CALL FLOPN2(NFL)
!C      IF (NFL.EQ.0) GO TO 999
!     H SPECIFICATION
!C      READ(5,1) H
!     SGME2 AND A INPUT
!     THE OUTPUTS OF PROGRAM 5.3.1 FPE AUTO CAN BE USED AS THE FOLLOWING
!     INPUTS WITH K=0.
!C      READ(5,1) N,L
!C      READ(5,2) SGME2
!C      IF(L.LE.0) GO TO 300
!C      READ(5,2) (A(I),I=1,L)
!     K INPUT
!C  300 READ(5,1) K
!C      IF(K.LE.0) GO TO 310
!C      READ(5,2) (B(I),I=1,K)
!C  310 H1=H+1
h1=h+1
l1=l+1
k1=k+1
g(1)=1.0d-00
if(l.le.0) go to 400
do 10 i=1,l
i1=i+1
!xx   10 G(I1)=-A(I)
g(i1)=-a(i)
10 continue
400 call fouger(g,l1,gr1,gi1,h1)
g(1)=1.0d-00
if(k.le.0) go to 410
do 20 i=1,k
i1=i+1
!xx   20 G(I1)=B(I)
g(i1)=b(i)
20 continue
410 call fouger(g,k1,gr2,gi2,h1)
do 30 i=1,h1
!xx   30 PXX(I)=(GR2(I)**2+GI2(I)**2)/(GR1(I)**2+GI1(I)**2)*SGME2
pxx(i)=(gr2(i)**2+gi2(i)**2)/(gr1(i)**2+gi1(i)**2)*sgme2
30 continue
!C      WRITE(6,60)
!C      WRITE(6,160)
!C      WRITE(6,61) L,K,H
!C      WRITE(6,164) SGME2
!C      IF(L.LE.0) GO TO 500
!C      WRITE(6,62)
!C      CALL PRCOL1(A,1,L,0)
!C  500 IF(K.LE.0) GO TO 510
!C      WRITE(6,63)
!C      CALL PRCOL1(B,1,K,0)
!C  510 WRITE(6,64)
!C      CALL PRCOL1(PXX,1,H1,1)
!C      CALL FLCLS2(NFL)
!C  999 CONTINUE
!C    1 FORMAT(10I5)
!C    2 FORMAT(4D20.10)
!C   60 FORMAT(1H ,33HPROGRAM 5.4.1   RATIONAL SPECTRUM)
!C   61 FORMAT(1H ,2HL=,I5,2X,2HK=,I5,2X,2HH=,I5)
!C   62 FORMAT(1H ,4X,1HI,12X,4HA(I))
!C   63 FORMAT(1H ,4X,1HI,12X,4HB(I))
!C   64 FORMAT(1H ,4X,1HI,10X,6HPXX(I))
!C  160 FORMAT(1H ,17HINITIAL CONDITION)
!C  164 FORMAT(1H ,6HSGME2=,D12.5)
return
end
